#[macro_use]
extern crate criterion;

use criterion::{black_box, Criterion};

use bessel_rust;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("bessel_j_smallz_double_double", move |b| {
        b.iter(|| bessel_rust::bessel_j_smallz(black_box(0.0), black_box(1.0)))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
