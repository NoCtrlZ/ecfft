use super::curve::Ep;
use super::isogeny::Isogeny;
use pairing::bn256::Fq as Fp;

// precomputed params for ecfft
#[derive(Clone, Debug)]
pub(crate) struct EcFftParams {
    params: Vec<EcFftCache>,
}

#[derive(Clone, Debug)]
struct EcFftCache {
    domain: (Vec<Fp>, Vec<Fp>),
}

impl EcFftParams {
    pub fn new(k: usize) -> Self {
        assert!(k == 14);
        let n = 1 << k;
        let acc = Ep::generator();
        let presentative = Ep::representative();

        let mut params = Vec::new();
        let mut s = Vec::new();
        let mut s_prime = Vec::new();

        let isogeny = Isogeny::new(0);
        (0..n).for_each(|i| {
            let coset = presentative + acc * Fp::from_raw([i, 0, 0, 0]);
            if i % 2 == 0 {
                s.push(isogeny.evaluate(coset.to_affine().point_projective()))
            } else {
                s_prime.push(isogeny.evaluate(coset.to_affine().point_projective()))
            }
        });
        params.push(EcFftCache {
            domain: (s.clone(), s_prime.clone()),
        });

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
            params.push(EcFftCache {
                domain: (s.clone(), s_prime.clone()),
            });
        }

        EcFftParams { params }
    }

    pub fn get_domain(&self, depth: usize) -> (Vec<Fp>, Vec<Fp>) {
        (
            self.params[depth].domain.0.clone(),
            self.params[depth].domain.1.clone(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{EcFftParams, Isogeny};

    #[test]
    fn test_projective_domain() {
        let k = 14;
        let ecfft_params = EcFftParams::new(k);

        for i in 0..k {
            let (s, s_prime) = ecfft_params.get_domain(i);

            assert_eq!(1 << (k - (1 + i)), s.len());
            assert_eq!(1 << (k - (1 + i)), s_prime.len());
        }
    }

    #[test]
    fn test_isogeny_half_sizing() {
        let k = 14;
        let ecfft_params = EcFftParams::new(k);
        let (mut s, mut s_prime) = ecfft_params.get_domain(0);

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
