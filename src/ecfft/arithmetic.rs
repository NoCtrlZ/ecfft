use super::EcFftCache;

use pairing::bn256::Fq as Fp;
use rayon::{join, prelude::*};

// low degree extention using divide and conquer algorithm
pub(crate) fn serial_low_degree_extention(
    coeffs: &mut [Fp],
    coeffs_prime: &mut [Fp],
    n: usize,
    k: usize,
    depth: usize,
    caches: &EcFftCache,
) {
    if k == 2 {
        let cache = caches.get_tree(depth);
        let ((f0, f1), (f2, f3)) = cache.get_inv_factor()[0];
        let tmp = f2 * coeffs[0] + f3 * coeffs[1];
        coeffs[0] = f0 * coeffs[0] + f1 * coeffs[1];
        coeffs[1] = tmp;
        let tmp = f2 * coeffs_prime[0] + f3 * coeffs_prime[1];
        coeffs_prime[0] = f0 * coeffs_prime[0] + f1 * coeffs_prime[1];
        coeffs_prime[1] = tmp;
        let ((f0, f1), (f2, f3)) = cache.get_factor()[0];
        let tmp = f2 * coeffs[0] + f3 * coeffs[1];
        coeffs[0] = f0 * coeffs[0] + f1 * coeffs[1];
        coeffs[1] = tmp;
        let tmp = f2 * coeffs_prime[0] + f3 * coeffs_prime[1];
        coeffs_prime[0] = f0 * coeffs_prime[0] + f1 * coeffs_prime[1];
        coeffs_prime[1] = tmp;
        return;
    }

    let half_n = n / 2;
    let cache = caches.get_tree(depth);
    let (left, right) = coeffs.split_at_mut(half_n);
    let (left_prime, right_prime) = coeffs_prime.split_at_mut(half_n);

    serial_matrix_arithmetic(left, right, left_prime, right_prime, cache.get_inv_factor());
    join(
        || serial_low_degree_extention(left, left_prime, half_n, k - 1, depth + 1, caches),
        || serial_low_degree_extention(right, right_prime, half_n, k - 1, depth + 1, caches),
    );
    serial_matrix_arithmetic(left, right, left_prime, right_prime, cache.get_factor());
}

// low degree extention using divide and conquer algorithm
pub(crate) fn parallel_low_degree_extention(
    coeffs: &mut [Fp],
    coeffs_prime: &mut [Fp],
    n: usize,
    k: usize,
    depth: usize,
    caches: &EcFftCache,
    thread_log: usize,
) {
    if k == 2 {
        let cache = caches.get_tree(depth);
        let ((f0, f1), (f2, f3)) = cache.get_inv_factor()[0];
        let tmp = f2 * coeffs[0] + f3 * coeffs[1];
        coeffs[0] = f0 * coeffs[0] + f1 * coeffs[1];
        coeffs[1] = tmp;
        let tmp = f2 * coeffs_prime[0] + f3 * coeffs_prime[1];
        coeffs_prime[0] = f0 * coeffs_prime[0] + f1 * coeffs_prime[1];
        coeffs_prime[1] = tmp;
        let ((f0, f1), (f2, f3)) = cache.get_factor()[0];
        let tmp = f2 * coeffs[0] + f3 * coeffs[1];
        coeffs[0] = f0 * coeffs[0] + f1 * coeffs[1];
        coeffs[1] = tmp;
        let tmp = f2 * coeffs_prime[0] + f3 * coeffs_prime[1];
        coeffs_prime[0] = f0 * coeffs_prime[0] + f1 * coeffs_prime[1];
        coeffs_prime[1] = tmp;
        return;
    }

    let half_n = n / 2;
    let cache = caches.get_tree(depth);
    let (left, right) = coeffs.split_at_mut(half_n);
    let (left_prime, right_prime) = coeffs_prime.split_at_mut(half_n);

    if k > thread_log {
        parallel_matrix_arithmetic(left, right, left_prime, right_prime, cache.get_inv_factor());
        join(
            || {
                parallel_low_degree_extention(
                    left,
                    left_prime,
                    half_n,
                    k - 1,
                    depth + 1,
                    caches,
                    thread_log,
                )
            },
            || {
                parallel_low_degree_extention(
                    right,
                    right_prime,
                    half_n,
                    k - 1,
                    depth + 1,
                    caches,
                    thread_log,
                )
            },
        );
        parallel_matrix_arithmetic(left, right, left_prime, right_prime, cache.get_factor());
    } else {
        serial_matrix_arithmetic(left, right, left_prime, right_prime, cache.get_inv_factor());
        join(
            || serial_low_degree_extention(left, left_prime, half_n, k - 1, depth + 1, caches),
            || serial_low_degree_extention(right, right_prime, half_n, k - 1, depth + 1, caches),
        );
        serial_matrix_arithmetic(left, right, left_prime, right_prime, cache.get_factor());
    }
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
        .for_each(|((((a, b), c), d), e)| {
            let ((f0, f1), (f2, f3)) = e;
            let tmp = f2 * *a + f3 * *b;
            *a = f0 * *a + f1 * *b;
            *b = tmp;
            let tmp = f2 * *c + f3 * *d;
            *c = f0 * *c + f1 * *d;
            *d = tmp;
        })
}

// matrix arithmetic with factor
pub fn parallel_matrix_arithmetic(
    left: &mut [Fp],
    right: &mut [Fp],
    left_prime: &mut [Fp],
    right_prime: &mut [Fp],
    factor: &Vec<((Fp, Fp), (Fp, Fp))>,
) {
    left.par_iter_mut()
        .zip(right.par_iter_mut())
        .zip(left_prime.par_iter_mut())
        .zip(right_prime.par_iter_mut())
        .zip(factor.par_iter())
        .for_each(|((((a, b), c), d), e)| {
            let ((f0, f1), (f2, f3)) = e;
            let tmp = f2 * *a + f3 * *b;
            *a = f0 * *a + f1 * *b;
            *b = tmp;
            let tmp = f2 * *c + f3 * *d;
            *c = f0 * *c + f1 * *d;
            *d = tmp;
        })
}
