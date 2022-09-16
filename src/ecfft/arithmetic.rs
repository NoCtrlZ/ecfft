use super::EcFftCache;

use pairing::bn256::Fq as Fp;
use rayon::{join, prelude::*};

// low degree extention using divide and conquer algorithm
pub(crate) fn serial_low_degree_extention(
    coeffs: &mut [Fp],
    coeffs_prime: &mut [Fp],
    mut n: usize,
    k: usize,
    _: usize,
    caches: &EcFftCache,
) {
    for depth in 0..(k - 2) {
        let inv_factor = &caches.trees[depth].get_inv_factor();
        poly_conversion(coeffs, coeffs_prime, n, n >> 1, inv_factor);
        n >>= 1;
    }
    coeffs
        .chunks_mut(2)
        .zip(coeffs_prime.chunks_mut(2))
        .for_each(|(coeffs, coeffs_prime)| {
            bottom_poly_conversion(coeffs, coeffs_prime, caches.get_last_tree())
        });
    n = 4;
    for depth in 0..(k - 2) {
        let factor = &caches.trees[k - depth - 3].get_factor();
        poly_conversion(coeffs, coeffs_prime, n, n >> 1, factor);
        n <<= 1;
    }
}

// low degree extention using divide and conquer algorithm
pub(crate) fn parallel_low_degree_extention(
    coeffs: &mut [Fp],
    coeffs_prime: &mut [Fp],
    n: usize,
    k: usize,
    depth: usize,
    caches: &EcFftCache,
) {
    if k < 4 {
        serial_low_degree_extention(coeffs, coeffs_prime, n, k, depth, caches);
        return;
    }

    let half_n = n >> 1;
    let cache = caches.get_tree(depth);
    let (left, right) = coeffs.split_at_mut(half_n);
    let (left_prime, right_prime) = coeffs_prime.split_at_mut(half_n);
    let (factor, inv_factor) = cache.get_factors();

    serial_matrix_arithmetic(left, right, left_prime, right_prime, inv_factor);
    join(
        || parallel_low_degree_extention(left, left_prime, half_n, k - 1, depth + 1, caches),
        || parallel_low_degree_extention(right, right_prime, half_n, k - 1, depth + 1, caches),
    );
    serial_matrix_arithmetic(left, right, left_prime, right_prime, factor);
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

pub(crate) fn poly_conversion(
    coeffs: &mut [Fp],
    coeffs_prime: &mut [Fp],
    n: usize,
    half_n: usize,
    factor: &Vec<((Fp, Fp), (Fp, Fp))>,
) {
    coeffs
        .chunks_mut(n)
        .zip(coeffs_prime.chunks_mut(n))
        .for_each(|(coeffs, coeffs_prime)| {
            let (left, right) = coeffs.split_at_mut(half_n);
            let (left_prime, right_prime) = coeffs_prime.split_at_mut(half_n);
            serial_matrix_arithmetic(left, right, left_prime, right_prime, factor);
        });
}

pub(crate) fn bottom_poly_conversion(
    coeffs: &mut [Fp],
    coeffs_prime: &mut [Fp],
    (((f4, f5), (f6, f7)), ((f0, f1), (f2, f3))): (&((Fp, Fp), (Fp, Fp)), &((Fp, Fp), (Fp, Fp))),
) {
    let tmp_a = f2 * coeffs[0] + f3 * coeffs[1];
    let tmp_b = f0 * coeffs[0] + f1 * coeffs[1];
    let tmp_prime_a = f2 * coeffs_prime[0] + f3 * coeffs_prime[1];
    let tmp_prime_b = f0 * coeffs_prime[0] + f1 * coeffs_prime[1];
    coeffs[0] = f4 * tmp_b + f5 * tmp_a;
    coeffs[1] = f6 * tmp_b + f7 * tmp_a;
    coeffs_prime[0] = f4 * tmp_prime_b + f5 * tmp_prime_a;
    coeffs_prime[1] = f6 * tmp_prime_b + f7 * tmp_prime_a;
}
