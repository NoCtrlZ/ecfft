#[macro_use]
extern crate criterion;

use ecfft::{Coefficients, EcFft, Polynomial};

use pairing::bn256::Fq;
use pairing::group::ff::Field;

use criterion::{BenchmarkId, Criterion};
use rand_core::OsRng;

fn criterion_benchmark(c: &mut Criterion) {
    let mut ecfft_group = c.benchmark_group("ecfft");
    let ecfft = EcFft::new();
    for k in 1..15 {
        ecfft_group.bench_function(BenchmarkId::new("k", k), |b| {
            let poly_a = Polynomial::<Fq, Coefficients>::new(
                (0..(1 << k)).map(|_| Fq::random(OsRng)).collect::<Vec<_>>(),
            );
            b.iter(|| ecfft.evaluate(k, poly_a.clone()));
        });
    }
    ecfft_group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
