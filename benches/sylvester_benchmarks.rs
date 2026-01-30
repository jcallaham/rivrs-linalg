//! Benchmarks for Sylvester equation solvers.
//!
//! Measures performance of continuous and discrete solvers at various matrix sizes.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use faer::prelude::*;
use rand::prelude::*;
use rand_distr::StandardNormal;

/// Generate a random dense n×n matrix with eigenvalues having real parts
/// offset by `shift` to ensure well-conditioning.
fn random_matrix(n: usize, shift: f64, rng: &mut impl Rng) -> Mat<f64> {
    let mut m = Mat::from_fn(n, n, |_, _| rng.sample::<f64, _>(StandardNormal));
    // Add shift to diagonal to control eigenvalue locations
    for i in 0..n {
        m[(i, i)] += shift;
    }
    m
}

/// Generate a random n×m matrix.
fn random_rhs(n: usize, m: usize, rng: &mut impl Rng) -> Mat<f64> {
    Mat::from_fn(n, m, |_, _| rng.sample::<f64, _>(StandardNormal))
}

fn bench_continuous(c: &mut Criterion) {
    let mut group = c.benchmark_group("continuous");
    let sizes: Vec<usize> = vec![10, 20, 50, 100, 200];

    for &size in &sizes {
        group.bench_with_input(
            BenchmarkId::new("solve_continuous", size),
            &size,
            |b, &size| {
                let mut rng = StdRng::seed_from_u64(42);
                let a = random_matrix(size, 1.0, &mut rng);
                let bb = random_matrix(size, 5.0, &mut rng);
                let cc = random_rhs(size, size, &mut rng);
                b.iter(|| {
                    csrrs::sylvester::solve_continuous(a.as_ref(), bb.as_ref(), cc.as_ref())
                        .unwrap()
                });
            },
        );
    }
    group.finish();
}

fn bench_discrete(c: &mut Criterion) {
    let mut group = c.benchmark_group("discrete");
    let sizes: Vec<usize> = vec![10, 20, 50, 100, 200];

    for &size in &sizes {
        group.bench_with_input(
            BenchmarkId::new("solve_discrete", size),
            &size,
            |b, &size| {
                let mut rng = StdRng::seed_from_u64(42);
                // Use smaller eigenvalues for discrete case (stable systems)
                let a = random_matrix(size, 0.0, &mut rng) * 0.3;
                let bb = random_matrix(size, 0.0, &mut rng) * 0.3;
                let cc = random_rhs(size, size, &mut rng);
                b.iter(|| {
                    csrrs::sylvester::solve_discrete(a.as_ref(), bb.as_ref(), cc.as_ref()).unwrap()
                });
            },
        );
    }
    group.finish();
}

fn bench_triangular(c: &mut Criterion) {
    let mut group = c.benchmark_group("triangular");
    let sizes: Vec<usize> = vec![10, 20, 50, 100, 200];

    for &size in &sizes {
        group.bench_with_input(
            BenchmarkId::new("solve_triangular", size),
            &size,
            |b, &size| {
                // Create upper triangular matrices directly
                let mut rng = StdRng::seed_from_u64(42);
                let mut a = Mat::zeros(size, size);
                let mut bb = Mat::zeros(size, size);
                for i in 0..size {
                    for j in i..size {
                        a[(i, j)] = rng.sample::<f64, _>(StandardNormal);
                        bb[(i, j)] = rng.sample::<f64, _>(StandardNormal);
                    }
                    a[(i, i)] += 1.0; // Shift for conditioning
                    bb[(i, i)] += 5.0;
                }
                let cc_orig = random_rhs(size, size, &mut rng);

                b.iter(|| {
                    let mut cc = cc_orig.clone();
                    csrrs::sylvester::triangular::solve_triangular_sylvester(
                        a.as_ref(),
                        bb.as_ref(),
                        cc.as_mut(),
                        1.0,
                    )
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_continuous, bench_discrete, bench_triangular);
criterion_main!(benches);
