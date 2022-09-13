use super::EcFftCache;

use pairing::bn256::Fq as Fp;
use rayon::{join, prelude::*};

// low degree extention using divide and conquer algorithm
pub(crate) fn low_degree_extention(
    coeffs: &mut [Fp],
    coeffs_prime: &mut [Fp],
    n: usize,
    k: usize,
    depth: usize,
    caches: &EcFftCache,
) {
    if k == 2 {
        let cache = caches.get_last_tree();
        let (a, b) = coeffs.split_at_mut(1);
        let (c, d) = coeffs_prime.split_at_mut(1);
        matrix_arithmetic(
            &mut a[0],
            &mut b[0],
            &mut c[0],
            &mut d[0],
            &cache.get_inv_factor()[0],
        );
        matrix_arithmetic(
            &mut a[0],
            &mut b[0],
            &mut c[0],
            &mut d[0],
            &cache.get_factor()[0],
        );
        return;
    }

    let half_n = n / 2;
    let cache = caches.get_tree(depth);
    let (left, right) = coeffs.split_at_mut(half_n);
    let (left_prime, right_prime) = coeffs_prime.split_at_mut(half_n);

    serial_matrix_arithmetic(left, right, left_prime, right_prime, cache.get_inv_factor());
    join(
        || low_degree_extention(left, left_prime, half_n, k - 1, depth + 1, caches),
        || low_degree_extention(right, right_prime, half_n, k - 1, depth + 1, caches),
    );
    serial_matrix_arithmetic(left, right, left_prime, right_prime, cache.get_factor());
}

// matrix arithmetic with factor
pub fn serial_matrix_arithmetic(
    left: &mut [Fp],
    right: &mut [Fp],
    left_prime: &mut [Fp],
    right_prime: &mut [Fp],
    factor: &Vec<((Fp, Fp), (Fp, Fp))>,
) {
    left.iter_mut()
        .zip(right.iter_mut())
        .zip(left_prime.iter_mut())
        .zip(right_prime.iter_mut())
        .zip(factor.iter())
        .for_each(|((((a, b), c), d), e)| matrix_arithmetic(a, b, c, d, e))
}

fn matrix_arithmetic(
    a: &mut Fp,
    b: &mut Fp,
    c: &mut Fp,
    d: &mut Fp,
    factor: &((Fp, Fp), (Fp, Fp)),
) {
    let ((f0, f1), (f2, f3)) = factor;
    let tmp = f2 * *a + f3 * *b;
    *a = f0 * *a + f1 * *b;
    *b = tmp;
    let tmp = f2 * *c + f3 * *d;
    *c = f0 * *c + f1 * *d;
    *d = tmp;
}

pub(crate) fn serial_integrate_evaluation(
    coeffs: &mut [Fp],
    low_prime: Vec<Fp>,
    high_prime: Vec<Fp>,
    low: Vec<Fp>,
    high: Vec<Fp>,
    powered_coset: &Vec<Fp>,
) {
    coeffs
        .chunks_mut(2)
        .zip(low_prime.iter())
        .zip(high_prime.iter())
        .zip(low.iter())
        .zip(high.iter())
        .zip(powered_coset.chunks(2))
        .for_each(|(((((coeffs, a), b), c), d), e)| {
            coeffs[0] = a + e[0] * b;
            coeffs[1] = c + e[1] * d;
        });
}

pub(crate) fn parallel_integrate_evaluation(
    coeffs: &mut [Fp],
    low_prime: Vec<Fp>,
    high_prime: Vec<Fp>,
    low: Vec<Fp>,
    high: Vec<Fp>,
    powered_coset: &Vec<Fp>,
) {
    coeffs
        .par_chunks_mut(2)
        .zip(low_prime.par_iter())
        .zip(high_prime.par_iter())
        .zip(low.par_iter())
        .zip(high.par_iter())
        .zip(powered_coset.par_chunks(2))
        .for_each(|(((((coeffs, a), b), c), d), e)| {
            coeffs[0] = a + e[0] * b;
            coeffs[1] = c + e[1] * d;
        });
}
