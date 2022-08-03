mod behave;
mod curve;
mod params;

use params::EcFftParams;

use ff::Field;
use group::ff::PrimeField;
use pasta_curves::arithmetic::*;

use self::params::Isogenies;

// ecfft structure
#[derive(Clone, Debug)]
pub struct ECFft<G: Group> {
    // polynomial degree 2^k
    k: u32,
    // coset used as evaluation domain
    coset: Vec<G::Scalar>,
    // isogenies used for halving domain
    isogenies: Isogenies<G>,
}

impl<G: Group> ECFft<G> {
    pub fn new(k: u32) -> Self {
        let n = 1 << k;
        let half_n = n / 2;
        let mut coset_generator: G::Scalar = EcFftParams::<G>::coset_generator();
        for _ in 0..G::Scalar::S - k {
            coset_generator = coset_generator.square();
        }

        let coset: Vec<_> = (0..half_n as usize)
            .scan(G::Scalar::one(), |w, _| {
                let tw = *w;
                w.group_scale(&coset_generator);
                Some(tw)
            })
            .collect();

        let isogenies = EcFftParams::isogenies();
        ECFft {
            k,
            coset,
            isogenies,
        }
    }
}
