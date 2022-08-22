mod behave;
mod curve;
mod isogeny;
mod utils;

use curve::Ep;
use utils::{matrix_arithmetic, EcFftCache};

use pairing::arithmetic::BaseExt;
use pairing::bn256::Fq as Fp;
use rayon::join;

// precomputed params for ecfft
#[derive(Clone, Debug)]
pub struct EcFft {
    // polynomial degree 2^k
    k: usize,
    // precomputed ecfft params
    caches: Vec<EcFftCache>,
}

impl EcFft {
    pub fn new(k: usize) -> Self {
        assert!(k == 14);
        let n = 1 << k;

        let acc = Ep::generator();
        let presentative = Ep::representative();
        let mut caches = Vec::new();
        let mut coset: Vec<Fp> = (0..n)
            .map(|i| {
                let coset_point = presentative + acc * Fp::from_raw([i, 0, 0, 0]);
                coset_point.to_affine().point_projective()
            })
            .collect();

        for i in 0..k {
            let cache = EcFftCache::new(k - i, coset.clone());
            caches.push(cache);
            coset = coset.into_iter().step_by(2).collect();
        }

        EcFft { k, caches }
    }

    pub fn evaluate(&self, coeffs: &mut [Fp]) -> Vec<Fp> {
        let n = 1 << (self.k - 1);
        assert_eq!(coeffs.len(), n);

        let mut coeffs_prime = coeffs.to_vec().clone();
        self.ecfft(coeffs, &mut coeffs_prime, self.k);
        coeffs_prime
    }

    fn ecfft(&self, coeffs: &mut [Fp], coeffs_prime: &mut [Fp], k: usize) {
        let n = 1 << k;
        assert_eq!(coeffs.len(), n);

        if n == 1 {}
        let (low, high) = coeffs.split_at_mut(n / 2);
        let (low_prime, high_prime) = coeffs_prime.split_at_mut(n / 2);

        join(
            || self.ecfft(low, low_prime, k - 1),
            || self.ecfft(high, high_prime, k - 1),
        );

        low_prime.copy_from_slice(low);
        high_prime.copy_from_slice(high);

        let coset = &self.caches[self.k - k].coset;

        assert_eq!(n, coset.len());

        (0..(n / 2)).for_each(|i| {
            coeffs[2 * i] =
                low_prime[i] + coset[2 * i].pow(&[n as u64 / 2, 0, 0, 0]) * high_prime[i];
        });

        join(|| self.extend(low_prime), || self.extend(high_prime));

        (0..(n / 2)).for_each(|i| {
            coeffs[2 * i + 1] =
                low_prime[i] + coset[2 * i + 1].pow(&[n as u64 / 2, 0, 0, 0]) * high_prime[i];
        });
    }

    // evaluate n/2 size of polynomial on n size coset
    fn extend(&self, coeffs: &mut [Fp]) {
        let n = 1 << (self.k - 1);
        assert_eq!(coeffs.len(), n);
        let log_n = n.trailing_zeros() as usize;

        low_degree_extention(coeffs, n, 0, &self.caches[self.k - log_n])
    }

    // transform n size of polynomial to normal form
    pub fn enter(&self, coeffs: &mut [Fp]) {
        let n = 1 << self.k;
        assert_eq!(coeffs.len(), n);
    }
}

// low degree extention using divide and conquer algorithm
fn low_degree_extention(coeffs: &mut [Fp], n: usize, depth: usize, caches: &EcFftCache) {
    if n == 1 {}
    let cache = caches.get_cache(depth);
    let (left, right) = coeffs.split_at_mut(n / 2);
    matrix_arithmetic(left, right, cache.get_inv_factor());
    join(
        || low_degree_extention(left, n / 2, depth + 1, caches),
        || low_degree_extention(right, n / 2, depth + 1, caches),
    );
    matrix_arithmetic(left, right, cache.get_factor());
}
