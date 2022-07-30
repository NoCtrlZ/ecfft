use ff::Field;

// evaluate coeffs form polynomial with given at
pub fn evaluate<F: Field>(coeffs: &[F], at: F) -> F {
    coeffs
        .iter()
        .rev()
        .fold(F::zero(), |acc, coeff| acc * at + coeff)
}
