mod classic_fft;
mod common;
mod ecfft;
mod naive;

pub use classic_fft::ClassicFft;
pub use common::{point_multiply_fq, point_multiply_fr};
pub use ecfft::EcFft;
pub use naive::{evaluate, naive_multiply_fq, naive_multiply_fr};

#[cfg(test)]
mod tests {
    use super::{
        naive_multiply_fq, naive_multiply_fr, point_multiply_fq, point_multiply_fr, ClassicFft,
        EcFft,
    };
    use pairing::bn256::{Fq, Fr};
    use pairing::group::ff::Field;
    use proptest::prelude::*;
    use rand_core::OsRng;

    fn arb_poly_fr(k: u32) -> Vec<Fr> {
        (0..(1 << k)).map(|_| Fr::random(OsRng)).collect::<Vec<_>>()
    }

    fn arb_poly_fq(k: u32) -> Vec<Fq> {
        (0..(1 << k)).map(|_| Fq::random(OsRng)).collect::<Vec<_>>()
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(100))]
        #[test]
        fn classic_fft_poly_multiplication_test(k in 3_u32..10) {
            let mut poly_a = arb_poly_fr(k-1);
            let mut poly_b = arb_poly_fr(k-1);

            // order(n^2) normal multiplication
            let poly_e = naive_multiply_fr(poly_a.clone(), poly_b.clone());

            // order(nlogn) classic fft multiplication
            let classic_fft = ClassicFft::new(k);
            poly_a.resize(1<<k, Fr::zero());
            poly_b.resize(1<<k, Fr::zero());
            classic_fft.dft(&mut poly_a);
            classic_fft.dft(&mut poly_b);
            let mut poly_f = point_multiply_fr(poly_a, poly_b);
            classic_fft.idft(&mut poly_f);

            assert_eq!(poly_e, poly_f)
        }
    }

    #[test]
    fn ecfft_poly_multiplication_test() {
        let k = 14;
        let mut poly_a = arb_poly_fq(k - 1);
        let mut poly_b = arb_poly_fq(k - 1);

        // order(n^2) normal multiplication
        let poly_c = naive_multiply_fq(poly_a.clone(), poly_b.clone());

        // order(nlog^2n) ecfft multplication
        let ecfft = EcFft::new(k);
        println!("precomputed");
        ecfft.extend(&mut poly_a);
        ecfft.extend(&mut poly_b);
    }
}
