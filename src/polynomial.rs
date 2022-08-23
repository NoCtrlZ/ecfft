use pairing::bn256::{Fq, Fr};
use pairing::group::ff::Field;

use std::fmt::Debug;
use std::marker::PhantomData;

/// The basis over which a polynomial is described.
pub trait Basis: Copy + Debug + Send + Sync {}

/// The polynomial coefficient representation
#[derive(Clone, Copy, Debug)]
pub struct Coefficients;
impl Basis for Coefficients {}

/// The polynomial point-value representation
#[derive(Clone, Copy, Debug)]
pub struct PointValue;
impl Basis for PointValue {}

#[derive(Clone, Debug)]
pub struct Polynomial<F, B> {
    values: Vec<F>,
    _marker: PhantomData<B>,
}

impl<F: Field, B: Basis> Polynomial<F, B> {
    pub fn new(coeffs: Vec<F>) -> Polynomial<F, Coefficients> {
        Polynomial {
            values: coeffs,
            _marker: PhantomData,
        }
    }

    pub fn get_values(self) -> Vec<F> {
        self.values
    }

    // order(n) polynomials points multiplication
    pub fn point_multiply(self, b: Polynomial<F, PointValue>) -> Polynomial<F, PointValue> {
        let values = self
            .values
            .iter()
            .zip(b.values.iter())
            .map(|(a, b)| *a * *b)
            .collect::<Vec<_>>();
        Polynomial {
            values,
            _marker: PhantomData,
        }
    }

    // order(n^2) polynomials coefficients multiplication
    pub fn naive_multiply(self, b: Polynomial<F, Coefficients>) -> Polynomial<F, PointValue> {
        let mut c = vec![F::zero(); self.values.len() + b.values.len()];
        self.values.iter().enumerate().for_each(|(i_a, coeff_a)| {
            b.values.iter().enumerate().for_each(|(i_b, coeff_b)| {
                c[i_a + i_b] += *coeff_a * *coeff_b;
            })
        });
        Polynomial {
            values: c,
            _marker: PhantomData,
        }
    }

    // order(n) polynomials points multiplication
    pub fn polynomial_evaluation(self, x: F) -> F {
        self.values
            .iter()
            .rev()
            .fold(F::zero(), |acc, coeff| acc * x + coeff)
    }

    // order(n^2) transform coeffitient to point value representation
    pub fn to_point_value(self, domain: Vec<F>) -> Polynomial<F, PointValue> {
        let values = domain
            .iter()
            .map(|x| self.clone().polynomial_evaluation(*x))
            .collect::<Vec<_>>();
        Polynomial {
            values,
            _marker: PhantomData,
        }
    }
}

// order(n) polynomials points multiplication
pub fn point_multiply_fr(a: Vec<Fr>, b: Vec<Fr>) -> Vec<Fr> {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).collect()
}

#[cfg(test)]
mod tests {
    use super::{Coefficients, Fq, Polynomial};
    use pairing::arithmetic::BaseExt;
    use pairing::group::ff::Field;
    use proptest::{collection::vec, prelude::*};
    use rand_core::OsRng;

    fn arb_poly_fq(k: u32) -> Polynomial<Fq, Coefficients> {
        Polynomial::<Fq, Coefficients>::new(
            (0..(1 << k)).map(|_| Fq::random(OsRng)).collect::<Vec<_>>(),
        )
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

            poly_a.clone().get_values().iter().for_each(|coeff| {
                eval += coeff * exp;
                exp *= point;
            });

            assert_eq!(poly_a.polynomial_evaluation(point), eval)
        }
    }
}
