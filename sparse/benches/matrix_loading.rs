//! Benchmarks for test matrix loading and validation.

use criterion::{Criterion, criterion_group, criterion_main};
use rivrs_sparse::io::registry;
use rivrs_sparse::validate;

fn bench_load_mtx(c: &mut Criterion) {
    c.bench_function("load_mtx arrow-10-indef", |b| {
        b.iter(|| {
            let test = registry::load_test_matrix("arrow-10-indef")
                .expect("registry error")
                .expect("matrix should exist");
            criterion::black_box(test.matrix);
        });
    });
}

fn bench_load_and_reconstruct(c: &mut Criterion) {
    c.bench_function("load_and_reconstruct arrow-10-indef", |b| {
        b.iter(|| {
            let test = registry::load_test_matrix("arrow-10-indef")
                .expect("registry error")
                .expect("matrix should exist");
            let refdata = test.reference.expect("reference should exist");
            let err = validate::reconstruction_error(&test.matrix, &refdata);
            criterion::black_box(err);
        });
    });
}

criterion_group!(benches, bench_load_mtx, bench_load_and_reconstruct);
criterion_main!(benches);
