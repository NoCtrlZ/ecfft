use super::curve::Ep;
use pairing::bn256::Fr as Fp;

// constant params for ecfft
#[derive(Clone, Debug)]
pub struct EcFftParams {
    pub(crate) domain: (Vec<Fp>, Vec<Fp>),
}

impl EcFftParams {
    pub fn new(k: usize) -> Self {
        assert!(k == 12);
        let n = 1 << k;
        // let mut acc = Ep::subgroup_generator();
        let presentative = Ep::representative();

        let mut even = Vec::new();
        let mut odd = Vec::new();

        (0..n).for_each(|i| {
            // let coset = presentative + acc;
            // acc = acc + Ep::subgroup_generator();
            // if i % 2 == 1 {
            //     odd.push(coset.to_affine().point_projective())
            // } else {
            //     even.push(coset.to_affine().point_projective())
            // }
        });

        EcFftParams {
            domain: (odd, even),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::EcFftParams;

    #[test]
    fn test_projective_domain() {
        let k = 12;
        let ecfft_params = EcFftParams::new(k);

        assert_eq!(1 << (k - 1), ecfft_params.domain.0.len());
        assert_eq!(1 << (k - 1), ecfft_params.domain.1.len());
    }
}
