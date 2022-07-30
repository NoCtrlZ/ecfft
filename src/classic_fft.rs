mod utils;

use ff::Field;
use group::ff::PrimeField;
use pasta_curves::arithmetic::*;
use rayon::{join, prelude::*};
use utils::{butterfly_arithmetic, swap_bit_reverse};

// classic fft structure
#[derive(Clone, Debug)]
pub struct ClassicFft<G: Group> {
    // polynomial degree 2^k
    k: u32,
    // generator of order 2^{k - 1} multiplicative group used as twiddle factors
    multiplicative_generator: G::Scalar,
    // multiplicative group generator inverse
    inv_multiplicative_generator: G::Scalar,
    // n inverse for inverse discrete fourier transform
    n_inv: G::Scalar,
}

impl<G: Group> ClassicFft<G> {
    pub fn new(k: u32) -> Self {
        let mut multiplicative_generator = G::Scalar::root_of_unity();
        let k_diff = G::Scalar::S - k;
        for _ in 0..k_diff {
            multiplicative_generator = multiplicative_generator.square()
        }
        let inv_multiplicative_generator = multiplicative_generator.invert().unwrap();
        let n_inv = G::Scalar::from(1 << k).invert().unwrap();

        ClassicFft {
            k,
            multiplicative_generator,
            inv_multiplicative_generator,
            n_inv,
        }
    }

    // perform classic fft
    pub fn fft(&self, coeffs: &mut [G]) {
        let n = 1 << self.k;
        assert_eq!(coeffs.len(), n);

        swap_bit_reverse(coeffs, n, self.k);

        // precompute twiddle factors
        let twiddle_factors: Vec<_> = (0..(n / 2) as usize)
            .scan(G::Scalar::one(), |w, _| {
                let tw = *w;
                w.group_scale(&self.multiplicative_generator);
                Some(tw)
            })
            .collect();
        println!("{:?}", twiddle_factors);

        classic_fft(coeffs, n, 1, &twiddle_factors)
    }

    // perform classic ifft
    pub fn ifft(&self, coeffs: &mut [G]) {
        let n = 1 << self.k;
        assert_eq!(coeffs.len(), n);

        swap_bit_reverse(coeffs, n, self.k);

        // precompute inverse twiddle factors
        let twiddle_factors: Vec<_> = (0..(n / 2) as usize)
            .scan(G::Scalar::one(), |w, _| {
                let tw = *w;
                w.group_scale(&self.inv_multiplicative_generator);
                Some(tw)
            })
            .collect();

        classic_fft(coeffs, n, 1, &twiddle_factors);
        coeffs
            .par_iter_mut()
            .for_each(|coeff| coeff.group_scale(&self.n_inv))
    }
}

// classic fft using divide and conquer algorithm
fn classic_fft<G: Group>(coeffs: &mut [G], n: usize, twiddle_chunk: usize, twiddles: &[G::Scalar]) {
    if n == 2 {
        let t = coeffs[1];
        coeffs[1] = coeffs[0];
        coeffs[0].group_add(&t);
        coeffs[1].group_sub(&t);
    } else {
        let (left, right) = coeffs.split_at_mut(n / 2);
        join(
            || classic_fft(left, n / 2, twiddle_chunk * 2, twiddles),
            || classic_fft(right, n / 2, twiddle_chunk * 2, twiddles),
        );
        butterfly_arithmetic(left, right, twiddle_chunk, twiddles)
    }
}

#[cfg(test)]
mod tests {
    use super::ClassicFft;
    use ff::Field;
    use pasta_curves::Fp;
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

            classic_fft.fft(&mut poly_a);
            classic_fft.ifft(&mut poly_a);

            assert_eq!(poly_a, poly_b)
        }
    }
}
