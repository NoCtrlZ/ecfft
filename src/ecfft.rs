mod arithmetic;
mod behave;
mod curve;
mod fftree;
mod isogeny;
mod utils;

use crate::polynomial::{Coefficients, PointValue, Polynomial};
use arithmetic::{parallel_integrate_evaluation, serial_integrate_evaluation};
pub(crate) use curve::Ep;
use utils::EcFftCache;

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
        let thread_num = current_num_threads();
        let thread_log = thread_num;
        let is_parallel = k > thread_log;

        if is_parallel {
            self.par_enter(&mut coeffs.values, k, thread_log)
        } else {
            self.enter(&mut coeffs.values, k)
        }

        Polynomial {
            values: coeffs.values,
            _marker: PhantomData,
        }
    }

    fn enter(&self, coeffs: &mut [Fp], k: usize) {
        if k == 1 {
            let cache = &self.caches[self.max_k - k];
            let (a, b) = (
                cache.powered_coset[0] * coeffs[1],
                cache.powered_coset[1] * coeffs[1],
            );
            coeffs[1] = coeffs[0];
            coeffs[0] += a;
            coeffs[1] += b;
            return;
        }

        let next_k = k - 1;
        let (low, high) = coeffs.split_at_mut(1 << next_k);

        join(|| self.enter(low, next_k), || self.enter(high, next_k));

        let (low_prime, high_prime) = (low.to_vec(), high.to_vec());
        let cache = &self.caches[self.max_k - k];

        cache.extend(low, high, k);

        let (low, high) = (low.to_vec(), high.to_vec());
        serial_integrate_evaluation(
            coeffs,
            low_prime,
            high_prime,
            low,
            high,
            &cache.powered_coset,
        );
    }

    fn par_enter(&self, coeffs: &mut [Fp], k: usize, thread_log: usize) {
        if k == 1 {
            let powered_coset = &self.caches[self.max_k - k].powered_coset;
            let (a, b) = (powered_coset[0] * coeffs[1], powered_coset[1] * coeffs[1]);
            coeffs[1] = coeffs[0];
            coeffs[0] += a;
            coeffs[1] += b;
            return;
        }

        let next_k = k - 1;
        let (low, high) = coeffs.split_at_mut(1 << next_k);
        let cache = &self.caches[self.max_k - k];

        join(
            || self.par_enter(low, next_k, thread_log),
            || self.par_enter(high, next_k, thread_log),
        );
        let (low_prime, high_prime) = (low.to_vec(), high.to_vec());
        cache.par_extend(low, high, k);
        let (low, high) = (low.to_vec(), high.to_vec());
        serial_integrate_evaluation(
            coeffs,
            low_prime,
            high_prime,
            low,
            high,
            &cache.powered_coset,
        );
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
