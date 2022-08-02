use ff::Field;
use pasta_curves::arithmetic::Group;

#[derive(Clone, Debug)]
pub struct Isogenies<G: Group> {
    numerator: [G::Scalar; 3],
    denominator: [G::Scalar; 2],
}

// constant params for ecfft
#[derive(Clone, Debug)]
pub struct EcFftParams<G: Group> {
    coset_generator: G::Scalar,
}

impl<G: Group> EcFftParams<G> {
    pub fn coset_generator() -> G::Scalar {
        G::Scalar::one()
    }

    pub fn isogenies() -> Isogenies<G> {
        Isogenies {
            numerator: [G::Scalar::one(), G::Scalar::one(), G::Scalar::one()],
            denominator: [G::Scalar::one(), G::Scalar::one()],
        }
    }
}
