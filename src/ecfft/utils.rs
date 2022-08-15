use super::curve::Ep;
use super::isogeny::Isogeny;
use pairing::bn256::Fq as Fp;

// constant params for ecfft
#[derive(Clone, Debug)]
pub struct EcFftParams {
    pub(crate) domain: (Vec<Fp>, Vec<Fp>),
}

impl EcFftParams {
    pub fn new(k: usize) -> Self {
        assert!(k == 14);
        let n = 1 << k;
        let acc = Ep::generator();
        let presentative = Ep::representative();

        let mut s = Vec::new();
        let mut s_prime = Vec::new();

        (0..n).for_each(|i| {
            let coset = presentative + acc * Fp::from_raw([i, 0, 0, 0]);
            if i % 2 == 0 {
                s.push(coset.to_affine().point_projective())
            } else {
                s_prime.push(coset.to_affine().point_projective())
            }
        });

        EcFftParams {
            domain: (s, s_prime),
        }
    }

    pub fn get_domain(&self) -> (Vec<Fp>, Vec<Fp>) {
        (self.domain.0.clone(), self.domain.1.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::{EcFftParams, Isogeny};

    #[test]
    fn test_projective_domain() {
        let k = 14;
        let ecfft_params = EcFftParams::new(k);
        let (s, s_prime) = ecfft_params.get_domain();

        assert_eq!(1 << (k - 1), s.len());
        assert_eq!(1 << (k - 1), s_prime.len());
    }

    #[test]
    fn test_isogeny_half_sizing() {
        let k = 14;
        let ecfft_params = EcFftParams::new(k);
        let (mut s, mut s_prime) = ecfft_params.get_domain();
        let isogeny = Isogeny::new(0);
        s = s.iter().map(|coeff| isogeny.evaluate(*coeff)).collect();
        s_prime = s_prime
            .iter()
            .map(|coeff| isogeny.evaluate(*coeff))
            .collect();

        for i in 1..k {
            let isogeny = Isogeny::new(i);
            s = s.iter().map(|coeff| isogeny.evaluate(*coeff)).collect();
            s_prime = s_prime
                .iter()
                .map(|coeff| isogeny.evaluate(*coeff))
                .collect();
            s.sort();
            s.dedup();
            s_prime.sort();
            s_prime.dedup();
            assert_eq!(1 << (k - (i + 1)), s.len());
            assert_eq!(1 << (k - (i + 1)), s_prime.len());
        }
    }
}
