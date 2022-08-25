use crate::ecfft::Ep;
use crate::polynomial::{Coefficients, Polynomial};

use pairing::bn256::{Fq, Fr};
use pairing::group::ff::Field;
use rand_core::OsRng;

pub(crate) fn arb_poly_fr(k: u32) -> Polynomial<Fr, Coefficients> {
    Polynomial::<Fr, Coefficients>::new(
        (0..(1 << k)).map(|_| Fr::random(OsRng)).collect::<Vec<_>>(),
    )
}

pub(crate) fn arb_poly_fq(k: usize) -> Polynomial<Fq, Coefficients> {
    Polynomial::<Fq, Coefficients>::new(
        (0..(1 << k)).map(|_| Fq::random(OsRng)).collect::<Vec<_>>(),
    )
}

pub(crate) fn layer_coset(depth: usize) -> Vec<Fq> {
    let k = 14;
    let n = 1 << k;
    let step = 1 << depth;
    let acc = Ep::generator();
    let presentative = Ep::representative();

    (0..n)
        .step_by(step)
        .map(|i| {
            let coset_point = presentative + acc * Fq::from_raw([i, 0, 0, 0]);
            coset_point.to_affine().point_projective()
        })
        .collect::<Vec<_>>()
}
