use ff::Field;
use pasta_curves::Fp;

// evaluate coeffs form polynomial with given at
pub fn evaluate<F: Field>(coeffs: &[F], at: F) -> F {
    coeffs
        .iter()
        .rev()
        .fold(F::zero(), |acc, coeff| acc * at + coeff)
}

// order(n^2) polynomials coefficients multiplication
pub fn naive_multiply(a: Vec<Fp>, b: Vec<Fp>) -> Vec<Fp> {
    let mut c = vec![Fp::zero(); a.len() + b.len()];
    a.iter().enumerate().for_each(|(i_a, coeff_a)| {
        b.iter().enumerate().for_each(|(i_b, coeff_b)| {
            c[i_a + i_b] += coeff_a * coeff_b;
        })
    });
    c
}
