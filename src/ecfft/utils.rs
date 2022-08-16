use super::curve::Ep;
use super::isogeny::Isogeny;
use pairing::bn256::Fq as Fp;

#[derive(Clone, Debug)]
pub(crate) struct EcFftCache {
    params: Vec<EcFftree>,
}

#[derive(Clone, Debug)]
struct EcFftree {
    pub(crate) domain: (Vec<Fp>, Vec<Fp>),
}

impl EcFftCache {
    pub fn new(k: u32) -> Self {
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
        params.push(EcFftree {
            domain: (s.clone(), s_prime.clone()),
        });

        for i in 1..k {
            let isogeny = Isogeny::new(i);
            s = s[..1 << (k - (1 + i))]
                .iter()
                .map(|coeff| isogeny.evaluate(*coeff))
                .collect();
            s_prime = s_prime[..1 << (k - (1 + i))]
                .iter()
                .map(|coeff| isogeny.evaluate(*coeff))
                .collect();
            params.push(EcFftree {
                domain: (s.clone(), s_prime.clone()),
            });
        }

        EcFftCache { params }
    }

    pub fn get_domain(&self, depth: usize) -> (Vec<Fp>, Vec<Fp>) {
        (
            self.params[depth].domain.0.clone(),
            self.params[depth].domain.1.clone(),
        )
    }
}

pub(crate) fn swap_bit_reverse(a: &mut [Fp], n: usize, k: u32) {
    assert!(k <= 64);
    let diff = 64 - k;
    for i in 0..n as u64 {
        let ri = i.reverse_bits() >> diff;
        if i < ri {
            a.swap(ri as usize, i as usize);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{EcFftCache, Isogeny};

    #[test]
    fn test_isogeny_and_domain() {
        let k = 14;
        let ecfft_params = EcFftCache::new(k);
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

            let (mut l, mut l_prime) = ecfft_params.get_domain(i as usize);
            l.sort();
            l_prime.sort();

            assert_eq!(s, l);
            assert_eq!(s_prime, l_prime);
        }
    }
}
