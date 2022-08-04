mod classic_fft;
mod common;
mod ecfft;
mod naive;

pub use classic_fft::ClassicFft;
pub use common::point_multiply;
pub use naive::{evaluate, naive_multiply};

#[cfg(test)]
mod tests {
    use super::{naive_multiply, point_multiply, ClassicFft};
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
        fn poly_multiplication_test(k in 3_u32..10) {
            let mut poly_a = arb_poly(k-1);
            let mut poly_b = arb_poly(k-1);
            let classic_fft = ClassicFft::new(k);

            // order(n^2) normal multiplication
            let poly_c = naive_multiply(poly_a.clone(), poly_b.clone());

            // order(nlogn) classic fft multiplication
            poly_a.resize(1<<k, Fp::zero());
            poly_b.resize(1<<k, Fp::zero());
            classic_fft.dft(&mut poly_a);
            classic_fft.dft(&mut poly_b);
            let mut poly_d = point_multiply(poly_a, poly_b);
            classic_fft.idft(&mut poly_d);

            assert_eq!(poly_c, poly_d)
        }
    }
}
