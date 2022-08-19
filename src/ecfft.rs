mod behave;
mod curve;
mod isogeny;
mod utils;

use utils::{butterfly_arithmetic, swap_bit_reverse, EcFftCache};

use pairing::bn256::Fq as Fp;
use rayon::{join, prelude::*};

// precomputed params for ecfft
#[derive(Clone, Debug)]
pub struct EcFft {
    // polynomial degree 2^k
    k: u32,
    // precomputed ecfft params
    caches: EcFftCache,
}

impl EcFft {
    pub fn new(k: u32) -> Self {
        assert!(k == 14);
        let caches = EcFftCache::new(k);

        EcFft { k, caches }
    }

    // evaluate n/2 size of polynomial on n size coset
    pub fn extend(&self, coeffs: &mut [Fp]) {
        let n = 1 << (self.k - 1);
        assert_eq!(coeffs.len(), n);

        swap_bit_reverse(coeffs, n, self.k - 1);

        ecfft_arithmetic(coeffs, n, 0, &self.caches)
    }

    // transform n size of polynomial to normal form
    pub fn enter(&self, coeffs: &mut [Fp]) {
        let n = 1 << self.k;
        assert_eq!(coeffs.len(), n);
    }
}

// ecfft using divide and conquer algorithm
fn ecfft_arithmetic(coeffs: &mut [Fp], n: usize, depth: usize, caches: &EcFftCache) {
    if n == 1 {
    } else {
        let cache = caches.get_cache(depth);
        let (left, right) = coeffs.split_at_mut(n / 2);
        butterfly_arithmetic(left, right, cache.get_inv_factor());
        join(
            || ecfft_arithmetic(left, n / 2, depth + 1, caches),
            || ecfft_arithmetic(right, n / 2, depth + 1, caches),
        );
        butterfly_arithmetic(left, right, cache.get_factor());
    }
}
