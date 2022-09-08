mod behave;
mod curve;
mod isogeny;
mod utils;

use crate::polynomial::{Coefficients, PointValue, Polynomial};
pub(crate) use curve::Ep;
use utils::{integrate_evaluation, EcFftCache};

use pairing::bn256::Fq as Fp;
use rayon::{current_num_threads, join};
use std::marker::PhantomData;

// precomputed params for ecfft
#[derive(Clone, Debug)]
pub struct EcFft {
    // polynomial degree 2^k
    max_k: usize,
    // precomputed ecfft params
    caches: Vec<EcFftCache>,
}

impl EcFft {
    pub fn new() -> Self {
        let max_k = 14;
        let n = 1 << max_k;

        let acc = Ep::generator();
        let presentative = Ep::representative();
        let mut caches = Vec::new();
        let mut coset = (0..n)
            .map(|i| {
                let coset_point = presentative + acc * Fp::from_raw([i, 0, 0, 0]);
                coset_point.to_affine().point_projective()
            })
            .collect::<Vec<_>>();

        for i in 0..max_k {
            let cache = EcFftCache::new(max_k - i, coset.clone());
            caches.push(cache);
            coset = coset.into_iter().step_by(2).collect();
        }

        EcFft { max_k, caches }
    }

    pub fn evaluate(
        &self,
        k: usize,
        mut coeffs: Polynomial<Fp, Coefficients>,
    ) -> Polynomial<Fp, PointValue> {
        assert!(k <= self.max_k);

        let thread_log = 1 << current_num_threads();
        self.enter(&mut coeffs.values, k, thread_log);
        Polynomial {
            values: coeffs.values,
            _marker: PhantomData,
        }
    }

    fn enter(&self, coeffs: &mut [Fp], k: usize, thread_log: usize) {
        if k == 0 {
            return;
        }

        let next_k = k - 1;
        let half_n = 1 << next_k;
        let (low, high) = coeffs.split_at_mut(half_n);

        join(
            || self.enter(low, next_k, thread_log),
            || self.enter(high, next_k, thread_log),
        );

        let cache = &self.caches[self.max_k - k];
        let (low_prime, high_prime) = (low.to_vec(), high.to_vec());

        cache.extend(low, high, k);

        let (low, high) = (low.to_vec(), high.to_vec());

        integrate_evaluation(
            coeffs,
            &cache.powered_coset,
            low_prime,
            high_prime,
            low,
            high,
        )
    }

    #[cfg(test)]
    pub(crate) fn get_cache(&self, k: usize) -> EcFftCache {
        self.caches[self.max_k - k].clone()
    }
}

#[cfg(test)]
mod tests {
    use super::EcFft;

    #[test]
    fn test_precomputed_params() {
        let k = 14;
        let ecfft = EcFft::new();

        for i in 1..k {
            let mut n = 1 << (k - i);
            let cache = ecfft.get_cache(k - i);
            let coset = cache.get_coset();

            assert_eq!(coset.len(), n);

            for j in 0..(k - i) {
                let tree = cache.get_tree(j);
                let (s, s_prime) = tree.get_domain();
                let factor = tree.get_factor();
                let inv_factor = tree.get_inv_factor();

                n /= 2;

                assert_eq!(s.len(), n);
                assert_eq!(s_prime.len(), n);
                assert_eq!(factor.len(), n / 2);
                assert_eq!(inv_factor.len(), n / 2);
            }
        }
    }
}
