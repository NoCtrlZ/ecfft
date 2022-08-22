use super::curve::Ep;
use super::isogeny::Isogeny;
use pairing::bn256::Fq as Fp;

#[derive(Clone, Debug)]
pub(crate) struct EcFftCache {
    pub(crate) k: usize,
    pub(crate) cache: Vec<FfTree>,
    pub(crate) coset: Vec<Fp>,
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
    pub fn new(k: usize, coset: Vec<Fp>) -> Self {
        let mut cache = Vec::new();
        let mut s = Vec::new();
        let mut s_prime = Vec::new();

        coset.iter().enumerate().for_each(|(i, elm)| {
            if i % 2 == 0 {
                s.push(*elm)
            } else {
                s_prime.push(*elm)
            }
        });

        for i in 0..k {
            let isogeny = Isogeny::new(i);
            let degree_n = 1 << (k - (1 + i));

            s = isogeny.domain_half_sizing(s, degree_n);
            s_prime = isogeny.domain_half_sizing(s_prime, degree_n);

            let isogeny = Isogeny::new(i + 1);
            let half_defree_n = degree_n / 2;
            let exp = &[(degree_n - 1) as u64, 0, 0, 0];

            let inv_factor = isogeny.get_factor(&s, half_defree_n, exp);
            let factor = isogeny.get_factor(&s_prime, half_defree_n, exp);

            cache.push(FfTree {
                domain: (s.clone(), s_prime.clone()),
                factor,
                inv_factor,
            });
        }

        EcFftCache { k, cache, coset }
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

    pub(crate) fn get_inv_factor(&self) -> &Vec<((Fp, Fp), (Fp, Fp))> {
        &self.inv_factor
    }
}

pub(crate) fn matrix_arithmetic(
    left: &mut [Fp],
    right: &mut [Fp],
    factor: &Vec<((Fp, Fp), (Fp, Fp))>,
) {
    left.iter_mut()
        .zip(right.iter_mut())
        .zip(factor.iter())
        .for_each(|((a, b), c)| {
            let ((f0, f1), (f2, f3)) = c;
            let (x, y) = (f0 * *a + f1 * *b, f2 * *a + f3 * *b);
            *a = x;
            *b = y;
        })
}

#[cfg(test)]
mod tests {
    use super::{EcFftCache, Ep, Fp, Isogeny};

    #[test]
    fn test_isogeny_and_domain() {
        let k = 14;
        let n = 1 << k;

        let acc = Ep::generator();
        let presentative = Ep::representative();
        let coset: Vec<Fp> = (0..n)
            .map(|i| {
                let coset_point = presentative + acc * Fp::from_raw([i, 0, 0, 0]);
                coset_point.to_affine().point_projective()
            })
            .collect();
        let ecfft_params = EcFftCache::new(k, coset);
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
        assert_eq!(s.len(), 1);
        assert_eq!(s_prime.len(), 1);
    }
}
