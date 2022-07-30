use pasta_curves::Fp;

// order(n) polynomials points multiplication
pub fn point_multiply(a: Vec<Fp>, b: Vec<Fp>) -> Vec<Fp> {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).collect()
}
