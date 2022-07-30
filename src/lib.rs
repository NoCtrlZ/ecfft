mod naive;
mod classic_fft;

pub use naive::evaluate;
pub use classic_fft::ClassicFft;

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

    // order(n^2) coefficients multiplication
    fn naive_multiply(a: Vec<Fp>, b: Vec<Fp>) -> Vec<Fp> {
        let mut c = vec![Fp::zero(); a.len() + b.len()];
        a.iter().enumerate().for_each(|(i_a, coeff_a)| {
            b.iter().enumerate().for_each(|(i_b, coeff_b)| {
                c[i_a + i_b] += coeff_a * coeff_b;
            })
        });
        c
    }

    // order(n) point multiplication
    fn point_multiply(a: Vec<Fp>, b: Vec<Fp>) -> Vec<Fp> {
        a.iter().zip(b.iter()).map(|(a, b)| a * b).collect()
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
            classic_fft.fft(&mut poly_a);
            classic_fft.fft(&mut poly_b);
            let mut poly_d = point_multiply(poly_a, poly_b);
            classic_fft.ifft(&mut poly_d);

            assert_eq!(poly_c, poly_d)
        }
    }
}
