use super::curve::Ep;
use super::isogeny::Isogeny;
use pairing::arithmetic::BaseExt;
use pairing::bn256::Fq as Fp;

#[derive(Clone, Debug)]
pub(crate) struct EcFftCache {
    pub(crate) cache: Vec<FfTree>,
}

#[derive(Clone, Debug)]
pub(crate) struct FfTree {
    // evaluation domain same size with polynomial
    domain: (Vec<Fp>, Vec<Fp>),
    // factor for performing multiplication
    factor: Vec<((Fp, Fp), (Fp, Fp))>,
    // inverse factor for performing multiplication
    inv_factor: Vec<((Fp, Fp), (Fp, Fp))>,
}

impl EcFftCache {
    pub fn new(k: u32) -> Self {
        assert!(k == 14);
        let n = 1 << k;
        let acc = Ep::generator();
        let presentative = Ep::representative();

        let mut cache = Vec::new();
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

        for i in 0..k {
            let isogeny = Isogeny::new(i);
            let degree_n = 1 << (k - (1 + i));
            let half_defree_n = degree_n / 2;
            let exp = &[(degree_n - 1) as u64, 0, 0, 0];
            s = s[..degree_n]
                .iter()
                .map(|coeff| isogeny.evaluate(*coeff))
                .collect();
            s_prime = s_prime[..degree_n]
                .iter()
                .map(|coeff| isogeny.evaluate(*coeff))
                .collect();
            let isogeny = Isogeny::new(i + 1);
            let inv_factor = s[..half_defree_n]
                .iter()
                .zip(&s[half_defree_n..])
                .map(|(a, b)| {
                    assert_eq!(isogeny.evaluate(*a), isogeny.evaluate(*b));
                    let f1 = isogeny.evaluate_with_denominator(*a).pow(exp);
                    let f2 = a * f1;
                    let f3 = isogeny.evaluate_with_denominator(*b).pow(exp);
                    let f4 = b * f3;
                    ((f1, f2), (f3, f4))
                })
                .collect();
            let factor = s_prime[..half_defree_n]
                .iter()
                .zip(&s_prime[half_defree_n..])
                .map(|(a, b)| {
                    assert_eq!(isogeny.evaluate(*a), isogeny.evaluate(*b));
                    let f1 = isogeny.evaluate_with_denominator(*a).pow(exp);
                    let f2 = a * f1;
                    let f3 = isogeny.evaluate_with_denominator(*b).pow(exp);
                    let f4 = b * f3;
                    ((f1, f2), (f3, f4))
                })
                .collect();
            cache.push(FfTree {
                domain: (s.clone(), s_prime.clone()),
                factor,
                inv_factor,
            });
        }

        EcFftCache { cache }
    }

    pub(crate) fn get_cache(&self, depth: usize) -> &FfTree {
        &self.cache[depth]
    }
}

impl FfTree {
    pub(crate) fn get_domain(&self) -> &(Vec<Fp>, Vec<Fp>) {
        &self.domain
    }

    pub(crate) fn get_factor(&self) -> &Vec<((Fp, Fp), (Fp, Fp))> {
        &self.factor
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

pub(crate) fn butterfly_arithmetic(coeffs: &mut [Fp], cache: &FfTree) {
    let (s, s_prime) = cache.get_domain();
    let factor = cache.get_factor();
    assert_eq!(coeffs.len(), s.len());
    assert_eq!(coeffs.len(), s_prime.len());
    assert_eq!(coeffs.len(), factor.len());
}

#[cfg(test)]
mod tests {
    use super::{BaseExt, EcFftCache, Isogeny};

    #[test]
    fn test_isogeny_and_domain() {
        let k = 14;
        let ecfft_params = EcFftCache::new(k);
        let cache = ecfft_params.get_cache(0);
        let (mut s, mut s_prime) = cache.domain.clone();

        for i in 1..k {
            let degree_n = 1 << (k - (1 + i));
            let half_defree_n = degree_n / 2;
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

            let cache = ecfft_params.get_cache(i as usize);
            let (mut l, mut l_prime) = cache.domain.clone();
            l.sort();
            l_prime.sort();

            assert_eq!(s, l);
            assert_eq!(s_prime, l_prime);

            let isogeny = Isogeny::new(i + 1);
            let (s, s_prime) = cache.domain.clone();
            s[..half_defree_n]
                .iter()
                .zip(&s[half_defree_n..])
                .zip(&s_prime[..half_defree_n])
                .zip(&s_prime[half_defree_n..])
                .for_each(|(((a, b), c), d)| {
                    assert_eq!(isogeny.evaluate(*a), isogeny.evaluate(*b));
                    assert_eq!(isogeny.evaluate(*c), isogeny.evaluate(*d));
                });
        }
    }
}
