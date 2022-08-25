use super::isogeny::Isogeny;
use crate::polynomial::{PointValue, Polynomial};

use pairing::bn256::Fq as Fp;
use rayon::join;

#[derive(Clone, Debug)]
pub(crate) struct EcFftCache {
    pub(crate) k: usize,
    pub(crate) trees: Vec<FfTree>,
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
        let mut trees = Vec::new();
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

            let (new_s, new_s_prime) = join(
                || isogeny.domain_half_sizing(s, degree_n),
                || isogeny.domain_half_sizing(s_prime, degree_n),
            );
            s = new_s;
            s_prime = new_s_prime;

            let isogeny = Isogeny::new(i + 1);
            let half_defree_n = degree_n / 2;
            let exp = &[(degree_n - 1) as u64, 0, 0, 0];

            let (inv_factor, factor) = join(
                || isogeny.get_factor(&s, half_defree_n, exp),
                || isogeny.get_factor(&s_prime, half_defree_n, exp),
            );

            trees.push(FfTree {
                domain: (s.clone(), s_prime.clone()),
                factor,
                inv_factor,
            });
        }

        EcFftCache { k, trees, coset }
    }

    pub(crate) fn get_tree(&self, depth: usize) -> &FfTree {
        &self.trees[depth]
    }

    #[cfg(test)]
    pub(crate) fn get_coset(&self) -> &Vec<Fp> {
        &self.coset
    }

    // evaluate n/2 size of polynomial on n size coset
    pub(crate) fn extend(&self, poly: &mut Polynomial<Fp, PointValue>) {
        let n = 1 << (self.k - 1);
        assert_eq!(poly.values.len(), n);

        low_degree_extention(&mut poly.values, n, 0, &self)
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

// low degree extention using divide and conquer algorithm
fn low_degree_extention(coeffs: &mut [Fp], n: usize, depth: usize, caches: &EcFftCache) {
    if n == 1 {
        return;
    }

    let cache = caches.get_tree(depth);
    let (left, right) = coeffs.split_at_mut(n / 2);
    matrix_arithmetic(left, right, cache.get_inv_factor());
    join(
        || low_degree_extention(left, n / 2, depth + 1, caches),
        || low_degree_extention(right, n / 2, depth + 1, caches),
    );
    matrix_arithmetic(left, right, cache.get_factor());
}

pub(crate) fn matrix_arithmetic(
    left: &mut [Fp],
    right: &mut [Fp],
    factor: &Vec<((Fp, Fp), (Fp, Fp))>,
) {
    assert_eq!(left.len(), factor.len());
    assert_eq!(right.len(), factor.len());
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
    use super::{EcFftCache, Isogeny};
    use crate::test::{arb_poly_fq, layer_coset};

    #[test]
    fn test_isogeny_and_domain() {
        let k = 14;

        let coset = layer_coset(0);
        let ecfft_params = EcFftCache::new(k, coset);
        let cache = ecfft_params.get_tree(0);
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

            let cache = ecfft_params.get_tree(i as usize);
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

    #[test]
    fn test_extend_operation() {
        let k = 1;
        let n = 1 << k;
        let depth = 14 - k;
        let coset = layer_coset(depth);
        let cache = EcFftCache::new(k, coset.clone());
        let tree = cache.get_tree(0);
        let (s, s_prime) = tree.get_domain();
        let factor = tree.get_factor();
        let inv_factor = tree.get_inv_factor();

        let coeff_a = arb_poly_fq(k - 1);
        let mut coeff_a_on_s = coeff_a.to_point_value(s);
        cache.extend(&mut coeff_a_on_s);
        let point_value_a_on_s_prime = coeff_a.to_point_value(s_prime);

        assert_eq!(coeff_a.values.len(), n / 2);
        assert_eq!(coset.len(), n);
        assert_eq!(cache.trees.len(), k);
        assert_eq!(s.len(), n / 2);
        assert_eq!(s_prime.len(), n / 2);
        assert_eq!(factor.len(), n / 4);
        assert_eq!(inv_factor.len(), n / 4);
        assert_eq!(coeff_a_on_s, point_value_a_on_s_prime);
    }
}
