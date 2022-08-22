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
    k: u32,
    // precomputed ecfft params
    caches: EcFftCache,
    // precomputed coset
    cosets: Vec<Vec<Fp>>,
}

impl EcFft {
    pub fn new(k: u32) -> Self {
        assert!(k == 14);
        let n = 1 << k;

        let acc = Ep::generator();
        let presentative = Ep::representative();
        let mut cosets = Vec::new();
        let mut coset: Vec<Fp> = (0..n)
            .map(|i| {
                let coset_point = presentative + acc * Fp::from_raw([i, 0, 0, 0]);
                coset_point.to_affine().point_projective()
            })
            .collect();
        let caches = EcFftCache::new(k, coset.clone());

        for _ in 0..k {
            cosets.push(coset.clone());
            coset = coset.into_iter().step_by(2).collect();
        }

        EcFft { k, caches, cosets }
    }

    pub fn evaluate(&self, coeffs: &mut [Fp]) -> Vec<Fp> {
        let n = 1 << (self.k - 1);
        assert_eq!(coeffs.len(), n);

        let mut coeffs_prime = coeffs.to_vec().clone();
        ecfft_arithmetic(coeffs, &mut coeffs_prime, n, 0, &self.cosets);
        coeffs_prime
    }

    // evaluate n/2 size of polynomial on n size coset
    pub fn extend(&self, coeffs: &mut [Fp]) {
        let n = 1 << (self.k - 1);
        assert_eq!(coeffs.len(), n);

        low_degree_extention(coeffs, n, 0, &self.caches)
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

// ecfft using divide and conquer algorithm
fn ecfft_arithmetic(
    coeffs: &mut [Fp],
    coeffs_prime: &mut [Fp],
    n: usize,
    depth: usize,
    cosets: &Vec<Vec<Fp>>,
) {
    if n == 1 {}
    let (low, high) = coeffs.split_at_mut(n / 2);
    let (low_prime, high_prime) = coeffs_prime.split_at_mut(n / 2);
    low_prime.copy_from_slice(low);
    high_prime.copy_from_slice(high);
    join(
        || ecfft_arithmetic(low, low_prime, n / 2, depth + 1, cosets),
        || ecfft_arithmetic(high, high_prime, n / 2, depth + 1, cosets),
    );
    let coset = &cosets[depth];
    assert_eq!(n, coset.len());
    (0..(n / 2)).for_each(|i| {
        coeffs[2 * i] = low_prime[i] + coset[2 * i].pow(&[n as u64 / 2, 0, 0, 0]) * high_prime[i];
    });
}
