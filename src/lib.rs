#[cfg(test)]
mod test;

mod classic_fft;
mod ecfft;
mod polynomial;

pub use crate::ecfft::EcFft;
pub use classic_fft::ClassicFft;
pub use polynomial::{Coefficients, Polynomial};

#[cfg(test)]
mod tests {
    use super::{ClassicFft, EcFft};
    use crate::test::{arb_poly_fq, arb_poly_fr, point_multiply_fr};
    use pairing::bn256::Fr;
    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(100))]
        #[test]
        fn classic_fft_poly_multiplication_test(k in 3_u32..10) {
            let poly_a = arb_poly_fr(k-1);
            let poly_b = arb_poly_fr(k-1);
            let mut poly_c = poly_a.clone().get_values();
            let mut poly_d = poly_b.clone().get_values();

            let classic_fft = ClassicFft::new(k);
            // order(n^2) normal multiplication
            let poly_e = poly_a.clone().naive_multiply(poly_b.clone());

            // order(nlogn) classic fft multiplication
            poly_c.resize(1<<k, Fr::zero());
            poly_d.resize(1<<k, Fr::zero());
            classic_fft.dft(&mut poly_c);
            classic_fft.dft(&mut poly_d);
            let mut poly_f = point_multiply_fr(poly_c, poly_d);
            classic_fft.idft(&mut poly_f);

            assert_eq!(poly_e.get_values(), poly_f)
        }
    }

    #[test]
    fn ecfft_poly_evaluation_test() {
        let ecfft = EcFft::new();
        let max_k = 14;
        for k in 1..max_k {
            let poly_a = arb_poly_fq(k);
            let cache = ecfft.get_cache(k);
            assert_eq!(cache.coset.len(), poly_a.values.len());

            // order(n^2) normal evaluation
            let poly_b = poly_a.clone().to_point_value(&cache.coset);

            // order(nlog^2n) ecfft evaluation
            let poly_c = ecfft.evaluate(k, poly_a.clone());

            assert_eq!(poly_b, poly_c)
        }
    }
}
