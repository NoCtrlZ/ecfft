use pairing::bn256::{Fq, Fr};

// order(n) polynomials points multiplication
pub fn point_multiply_fr(a: Vec<Fr>, b: Vec<Fr>) -> Vec<Fr> {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).collect()
}

// order(n) polynomials points multiplication
pub fn point_multiply_fq(a: Vec<Fq>, b: Vec<Fq>) -> Vec<Fq> {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).collect()
}
