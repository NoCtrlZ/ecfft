mod behave;
mod curve;
mod isogeny;
mod utils;

use crate::polynomial::{Coefficients, Polynomial};
use curve::Ep;
use utils::EcFftCache;

use pairing::arithmetic::BaseExt;
use pairing::bn256::Fq as Fp;
use rayon::join;

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
        let mut coset: Vec<Fp> = (0..n)
            .map(|i| {
                let coset_point = presentative + acc * Fp::from_raw([i, 0, 0, 0]);
                coset_point.to_affine().point_projective()
            })
            .collect();

        for i in 0..max_k {
            let cache = EcFftCache::new(max_k - i, coset.clone());
            caches.push(cache);
            coset = coset.into_iter().step_by(2).collect();
        }

        EcFft { max_k, caches }
    }

    pub fn evaluate(&self, coeffs: Polynomial<Fp, Coefficients>) -> Vec<Fp> {
        let n = 1 << (self.max_k - 1);

        let mut coeffs = coeffs.get_values();
        let mut coeffs_prime = coeffs.clone();

        assert_eq!(coeffs.len(), n);

        self.enter(&mut coeffs, &mut coeffs_prime, self.max_k - 1);
        coeffs_prime
    }

    fn enter(&self, coeffs: &mut [Fp], coeffs_prime: &mut [Fp], k: usize) {
        let n = 1 << k;
        assert_eq!(coeffs.len(), n);

        if n == 1 {
            return;
        }

        let (low, high) = coeffs.split_at_mut(n / 2);
        let (low_prime, high_prime) = coeffs_prime.split_at_mut(n / 2);

        join(
            || self.enter(low, low_prime, k - 1),
            || self.enter(high, high_prime, k - 1),
        );

        low_prime.copy_from_slice(low);
        high_prime.copy_from_slice(high);

        let cache = &self.caches[self.max_k - k];

        assert_eq!(n, cache.coset.len());

        (0..(n / 2)).for_each(|i| {
            coeffs[2 * i] =
                low_prime[i] + cache.coset[2 * i].pow(&[n as u64 / 2, 0, 0, 0]) * high_prime[i];
        });

        join(|| cache.extend(low_prime), || cache.extend(high_prime));

        (0..(n / 2)).for_each(|i| {
            coeffs[2 * i + 1] =
                low_prime[i] + cache.coset[2 * i + 1].pow(&[n as u64 / 2, 0, 0, 0]) * high_prime[i];
        });
    }

    #[cfg(test)]
    pub fn get_coset(&self, k: usize) -> Vec<Fp> {
        self.caches[self.max_k - k].coset.clone()
    }
}
