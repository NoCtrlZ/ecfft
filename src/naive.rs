use pairing::bn256::{Fq, Fr};
use pairing::group::ff::Field;

// evaluate coeffs form polynomial with given at
pub fn evaluate<F: Field>(coeffs: &[F], at: F) -> F {
    coeffs
        .iter()
        .rev()
        .fold(F::zero(), |acc, coeff| acc * at + coeff)
}

// order(n^2) polynomials coefficients multiplication
pub fn naive_multiply_fr(a: Vec<Fr>, b: Vec<Fr>) -> Vec<Fr> {
    let mut c = vec![Fr::zero(); a.len() + b.len()];
    a.iter().enumerate().for_each(|(i_a, coeff_a)| {
        b.iter().enumerate().for_each(|(i_b, coeff_b)| {
            c[i_a + i_b] += coeff_a * coeff_b;
        })
    });
    c
}

// order(n^2) polynomials coefficients multiplication
pub fn naive_multiply_fq(a: Vec<Fq>, b: Vec<Fq>) -> Vec<Fq> {
    let mut c = vec![Fq::zero(); a.len() + b.len()];
    a.iter().enumerate().for_each(|(i_a, coeff_a)| {
        b.iter().enumerate().for_each(|(i_b, coeff_b)| {
            c[i_a + i_b] += coeff_a * coeff_b;
        })
    });
    c
}
