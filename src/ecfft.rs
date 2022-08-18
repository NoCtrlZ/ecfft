mod behave;
mod curve;
mod isogeny;
mod utils;

use utils::{swap_bit_reverse, EcFftCache};

use pairing::bn256::Fq as Fp;
use rayon::{join, prelude::*};

// precomputed params for ecfft
#[derive(Clone, Debug)]
pub struct EcFft {
    // polynomial degree 2^k
    k: u32,
    // precomputed ecfft params
    cache: EcFftCache,
}

impl EcFft {
    pub fn new(k: u32) -> Self {
        assert!(k == 14);
        let cache = EcFftCache::new(k);

        EcFft { k, cache }
    }

    // perform extend operation
    pub fn extend(&self, coeffs: &mut [Fp]) {
        let n = 1 << self.k;
        assert_eq!(coeffs.len(), n);

        swap_bit_reverse(coeffs, n, self.k);

        ecfft_arithmetic(coeffs, n, 1, &self.cache)
    }

    // perform enter operation
    pub fn enter(&self, coeffs: &mut [Fp]) {}
}

// ecfft using divide and conquer algorithm
fn ecfft_arithmetic(coeffs: &mut [Fp], n: usize, depth: u32, cache: &EcFftCache) {
    if n == 1 {
    } else {
        let (left, right) = coeffs.split_at_mut(n / 2);
        join(
            || ecfft_arithmetic(left, n / 2, depth + 1, cache),
            || ecfft_arithmetic(right, n / 2, depth + 1, cache),
        );
    }
}
