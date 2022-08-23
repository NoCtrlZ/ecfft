mod utils;

use utils::{butterfly_arithmetic, swap_bit_reverse};

use pairing::arithmetic::*;
use pairing::group::ff::{Field, PrimeField};
use rayon::{join, prelude::*};

// classic fft structure
#[derive(Clone, Debug)]
pub struct ClassicFft<G: Group> {
    // polynomial degree 2^k
    k: u32,
    // generator of order 2^{k - 1} multiplicative group used as twiddle factors
    twiddle_factors: Vec<G::Scalar>,
    // multiplicative group generator inverse
    inv_twiddle_factors: Vec<G::Scalar>,
    // n inverse for inverse discrete fourier transform
    n_inv: G::Scalar,
}

impl<G: Group> ClassicFft<G> {
    pub fn new(k: u32) -> Self {
        let n = 1 << k;
        let half_n = n / 2;
        let mut multiplicative_generator = G::Scalar::root_of_unity();
        for _ in 0..G::Scalar::S - k {
            multiplicative_generator = multiplicative_generator.square()
        }
        let inv_multiplicative_generator = multiplicative_generator.invert().unwrap();
        let n_inv = G::Scalar::from(n).invert().unwrap();

        // precompute twiddle factors
        let twiddle_factors = (0..half_n as usize)
            .scan(G::Scalar::one(), |w, _| {
                let tw = *w;
                w.group_scale(&multiplicative_generator);
                Some(tw)
            })
            .collect::<Vec<_>>();

        // precompute inverse twiddle factors
        let inv_twiddle_factors = (0..half_n as usize)
            .scan(G::Scalar::one(), |w, _| {
                let tw = *w;
                w.group_scale(&inv_multiplicative_generator);
                Some(tw)
            })
            .collect::<Vec<_>>();

        ClassicFft {
            k,
            twiddle_factors,
            inv_twiddle_factors,
            n_inv,
        }
    }

    // perform classic discrete fourier transform
    pub fn dft(&self, coeffs: &mut [G]) {
        let n = 1 << self.k;
        assert_eq!(coeffs.len(), n);

        swap_bit_reverse(coeffs, n, self.k);

        classic_fft_arithmetic(coeffs, n, 1, &self.twiddle_factors)
    }

    // perform classic inverse discrete fourier transform
    pub fn idft(&self, coeffs: &mut [G]) {
        let n = 1 << self.k;
        assert_eq!(coeffs.len(), n);

        swap_bit_reverse(coeffs, n, self.k);

        classic_fft_arithmetic(coeffs, n, 1, &self.inv_twiddle_factors);
        coeffs
            .par_iter_mut()
            .for_each(|coeff| coeff.group_scale(&self.n_inv))
    }
}

// classic fft using divide and conquer algorithm
fn classic_fft_arithmetic<G: Group>(
    coeffs: &mut [G],
    n: usize,
    twiddle_chunk: usize,
    twiddles: &[G::Scalar],
) {
    if n == 2 {
        let t = coeffs[1];
        coeffs[1] = coeffs[0];
        coeffs[0].group_add(&t);
        coeffs[1].group_sub(&t);
    } else {
        let (left, right) = coeffs.split_at_mut(n / 2);
        join(
            || classic_fft_arithmetic(left, n / 2, twiddle_chunk * 2, twiddles),
            || classic_fft_arithmetic(right, n / 2, twiddle_chunk * 2, twiddles),
        );
        butterfly_arithmetic(left, right, twiddle_chunk, twiddles)
    }
}

#[cfg(test)]
mod tests {
    use super::ClassicFft;
    use pairing::bn256::Fr as Fp;
    use pairing::group::ff::Field;
    use proptest::prelude::*;
    use rand_core::OsRng;

    fn arb_poly(k: u32) -> Vec<Fp> {
        (0..(1 << k)).map(|_| Fp::random(OsRng)).collect::<Vec<_>>()
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(100))]
        #[test]
        fn classic_fft_test(k in 3_u32..10) {
            let mut poly_a = arb_poly(k);
            let poly_b = poly_a.clone();
            let classic_fft = ClassicFft::new(k);

            classic_fft.dft(&mut poly_a);
            classic_fft.idft(&mut poly_a);

            assert_eq!(poly_a, poly_b)
        }
    }
}
