#[macro_use]
extern crate criterion;

use ecfft::{Coefficients, EcFft, Polynomial};

use pairing::bn256::Fq;
use pairing::group::ff::Field;

use criterion::{BenchmarkId, Criterion};
use rand_core::OsRng;

fn criterion_benchmark(c: &mut Criterion) {
    let mut enter_group = c.benchmark_group("ecfft_enter");
    let ecfft = EcFft::new();
    for k in 10..15 {
        enter_group.bench_function(BenchmarkId::new("k", k), |b| {
            let poly_a = Polynomial::<Fq, Coefficients>::new(
                (0..(1 << k)).map(|_| Fq::random(OsRng)).collect::<Vec<_>>(),
            );
            b.iter(|| ecfft.evaluate(k, poly_a.clone()));
        });
    }
    enter_group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
