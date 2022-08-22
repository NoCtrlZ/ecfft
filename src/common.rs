use pairing::bn256::{Fq, Fr};

// order(n) polynomials points multiplication
pub fn point_multiply_fr(a: Vec<Fr>, b: Vec<Fr>) -> Vec<Fr> {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).collect()
}

// order(n) polynomials points multiplication
pub fn point_multiply_fq(a: Vec<Fq>, b: Vec<Fq>) -> Vec<Fq> {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).collect()
}

// order(n) polynomials points multiplication
pub fn polynomial_evaluation(coeffs: Vec<Fq>, x: Fq) -> Fq {
    coeffs
        .iter()
        .rev()
        .fold(Fq::zero(), |acc, coeff| acc * x + coeff)
}

#[cfg(test)]
mod tests {
    use super::{polynomial_evaluation, Fq};
    use pairing::arithmetic::BaseExt;
    use pairing::group::ff::Field;
    use proptest::{collection::vec, prelude::*};
    use rand_core::OsRng;

    fn arb_poly_fq(k: u32) -> Vec<Fq> {
        (0..(1 << k)).map(|_| Fq::random(OsRng)).collect::<Vec<_>>()
    }

    prop_compose! {
        fn arb_point()(
            bytes in vec(any::<u8>(), 64)
        ) -> Fq {
            Fq::from_bytes_wide(&<[u8; 64]>::try_from(bytes).unwrap())
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(100))]
        #[test]
        fn test_polynomial_evaluation(point in arb_point(),k in 1u32..20) {
            let mut eval = Fq::zero();
            let mut exp = Fq::one();
            let poly_a = arb_poly_fq(k);

            poly_a.iter().for_each(|coeff| {
                eval += coeff * exp;
                exp *= point;
            });

            assert_eq!(polynomial_evaluation(poly_a, point), eval)
        }
    }
}
