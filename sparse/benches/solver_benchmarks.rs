//! Criterion benchmarks for the SSIDS sparse symmetric indefinite solver.
//!
//! Groups:
//! - `ssids/factor` — analyze + factor on CI SuiteSparse matrices
//! - `ssids/solve` — solve alone (pre-factored) on CI SuiteSparse matrices
//! - `ssids/symbolic_analysis` — symbolic analysis (AMD ordering)
//! - `ssids/mc64_matching` — MC64 weighted bipartite matching
//! - `ssids/match_order_metis` — combined MC64 + METIS ordering pipeline
//! - `kernel/two_level`, `kernel/single_level` — dense APTP kernel comparison

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

use faer::Col;
use faer::Par;
use faer::dyn_stack::{MemBuffer, MemStack};

use rivrs_sparse::aptp::{
    AnalyzeOptions, FactorOptions, Mc64Job, SparseLDLT, match_order_metis, mc64_matching,
};
use rivrs_sparse::benchmarking::read_peak_rss_kb;
use rivrs_sparse::testing::{TestCaseFilter, load_test_cases};

/// Benchmark analyze+factor and solve on CI SuiteSparse matrices using SparseLDLT.
fn bench_factor_solve(c: &mut Criterion) {
    let cases = match load_test_cases(&TestCaseFilter::ci_subset()) {
        Ok(cases) => cases,
        Err(e) => {
            eprintln!("WARNING: Failed to load test cases: {}", e);
            return;
        }
    };

    let suitesparse: Vec<_> = cases
        .into_iter()
        .filter(|c| c.properties.source == "suitesparse")
        .collect();

    if suitesparse.is_empty() {
        eprintln!("WARNING: No SuiteSparse test cases matched the filter");
        return;
    }

    let opts = AnalyzeOptions::default();
    let factor_opts = FactorOptions::default();

    // Factor benchmark: analyze + factor (analysis is amortized setup)
    {
        let mut group = c.benchmark_group("ssids/factor");
        group.sample_size(10);

        for case in &suitesparse {
            let nnz = case.properties.nnz;
            group.throughput(Throughput::Elements(nnz as u64));

            group.bench_with_input(BenchmarkId::from_parameter(&case.name), &case, |b, case| {
                b.iter(|| {
                    let mut solver = SparseLDLT::analyze_with_matrix(&case.matrix, &opts).unwrap();
                    solver.factor(&case.matrix, &factor_opts).unwrap();
                });
            });
        }
        group.finish();
    }

    // Solve benchmark: solve only (pre-factored)
    {
        let mut group = c.benchmark_group("ssids/solve");
        group.sample_size(10);

        for case in &suitesparse {
            let n = case.matrix.nrows();
            let nnz = case.properties.nnz;
            group.throughput(Throughput::Elements(nnz as u64));

            let mut solver = SparseLDLT::analyze_with_matrix(&case.matrix, &opts)
                .expect("analyze should succeed");
            solver
                .factor(&case.matrix, &factor_opts)
                .expect("factor should succeed");

            let rhs = Col::from_fn(n, |i| (i % 7) as f64 - 3.0);
            let scratch_req = solver.solve_scratch(1);

            group.bench_with_input(
                BenchmarkId::from_parameter(&case.name),
                &case,
                |b, _case| {
                    b.iter(|| {
                        let mut mem = MemBuffer::new(scratch_req);
                        let stack = MemStack::new(&mut mem);
                        solver.solve(&rhs, stack, Par::Seq).unwrap()
                    });
                },
            );
        }
        group.finish();
    }
}

fn bench_symbolic_analysis(c: &mut Criterion) {
    use faer::sparse::linalg::cholesky::SymmetricOrdering;
    use rivrs_sparse::aptp::AptpSymbolic;

    let cases = match load_test_cases(&TestCaseFilter::ci_subset()) {
        Ok(cases) => cases,
        Err(e) => {
            eprintln!("WARNING: Failed to load test cases: {}", e);
            return;
        }
    };

    if cases.is_empty() {
        eprintln!("WARNING: No test cases matched the filter");
        return;
    }

    let mut group = c.benchmark_group("ssids/symbolic_analysis");
    group.sample_size(20);

    for case in &cases {
        let nnz = case.properties.nnz;
        group.throughput(Throughput::Elements(nnz as u64));

        group.bench_with_input(BenchmarkId::from_parameter(&case.name), &case, |b, case| {
            b.iter(|| {
                AptpSymbolic::analyze(case.matrix.symbolic(), SymmetricOrdering::Amd)
                    .expect("symbolic analysis should succeed")
            });
        });
    }
    group.finish();
}

fn bench_mc64_matching(c: &mut Criterion) {
    let rss_before = read_peak_rss_kb();

    let cases = match load_test_cases(&TestCaseFilter::ci_subset()) {
        Ok(cases) => cases,
        Err(e) => {
            eprintln!("WARNING: Failed to load test cases: {}", e);
            return;
        }
    };

    // Filter to SuiteSparse-source only (hand-constructed matrices are too small)
    let suitesparse: Vec<_> = cases
        .into_iter()
        .filter(|c| c.properties.source == "suitesparse")
        .collect();

    if suitesparse.is_empty() {
        eprintln!("WARNING: No SuiteSparse test cases matched the filter");
        return;
    }

    let mut group = c.benchmark_group("ssids/mc64_matching");
    group.sample_size(20);

    for case in &suitesparse {
        let nnz = case.properties.nnz;
        group.throughput(Throughput::Elements(nnz as u64));

        group.bench_with_input(BenchmarkId::from_parameter(&case.name), &case, |b, case| {
            b.iter(|| {
                mc64_matching(&case.matrix, Mc64Job::MaximumProduct)
                    .expect("MC64 matching should succeed")
            });
        });
    }
    group.finish();

    let rss_after = read_peak_rss_kb();
    if let (Some(before), Some(after)) = (rss_before, rss_after) {
        eprintln!(
            "\nMC64 Peak RSS: {} KB -> {} KB (delta: {} KB)",
            before,
            after,
            after.saturating_sub(before)
        );
    } else if let Some(after) = rss_after {
        eprintln!("\nMC64 Peak RSS: {} KB", after);
    }
}

fn bench_match_order(c: &mut Criterion) {
    let cases = match load_test_cases(&TestCaseFilter::ci_subset()) {
        Ok(cases) => cases,
        Err(e) => {
            eprintln!("WARNING: Failed to load test cases: {}", e);
            return;
        }
    };

    // Filter to SuiteSparse-source only (hand-constructed matrices are too small)
    let suitesparse: Vec<_> = cases
        .into_iter()
        .filter(|c| c.properties.source == "suitesparse")
        .collect();

    if suitesparse.is_empty() {
        eprintln!("WARNING: No SuiteSparse test cases matched the filter");
        return;
    }

    let mut group = c.benchmark_group("ssids/match_order_metis");
    group.sample_size(10);

    for case in &suitesparse {
        let nnz = case.properties.nnz;
        group.throughput(Throughput::Elements(nnz as u64));

        group.bench_with_input(BenchmarkId::from_parameter(&case.name), &case, |b, case| {
            b.iter(|| match_order_metis(&case.matrix).expect("match_order_metis should succeed"));
        });
    }
    group.finish();
}

/// Dense APTP kernel benchmark: two-level vs single-level comparison.
///
/// Benchmarks `aptp_factor` on random dense symmetric indefinite matrices,
/// comparing two-level blocking (default nb=256, ib=32) against single-level
/// (forced via outer_block_size=usize::MAX).
fn bench_kernel_two_level(c: &mut Criterion) {
    use faer::Mat;
    use rivrs_sparse::aptp::{AptpOptions, aptp_factor};

    let sizes = [128, 256, 512, 1024];

    // Generate deterministic test matrices
    let matrices: Vec<(usize, Mat<f64>)> = sizes
        .iter()
        .map(|&n| {
            let a = Mat::from_fn(n, n, |i, j| {
                let (i, j) = if i >= j { (i, j) } else { (j, i) };
                let seed = (i * 1000 + j * 7 + 13) as f64;
                let val = (seed * 0.618033988749).fract() * 2.0 - 1.0;
                if i == j { val * 10.0 } else { val }
            });
            (n, a)
        })
        .collect();

    // Two-level kernel (default block sizes)
    {
        let mut group = c.benchmark_group("kernel/two_level");
        group.sample_size(10);

        for (n, a) in &matrices {
            let opts = AptpOptions {
                outer_block_size: 256,
                inner_block_size: 32,
                ..AptpOptions::default()
            };

            group.bench_with_input(BenchmarkId::from_parameter(n), a, |b, a| {
                b.iter(|| aptp_factor(a.as_ref(), &opts).expect("factor should succeed"));
            });
        }
        group.finish();
    }

    // Single-level kernel (force via large outer_block_size)
    {
        let mut group = c.benchmark_group("kernel/single_level");
        group.sample_size(10);

        for (n, a) in &matrices {
            let opts = AptpOptions {
                outer_block_size: usize::MAX,
                inner_block_size: 32,
                ..AptpOptions::default()
            };

            group.bench_with_input(BenchmarkId::from_parameter(n), a, |b, a| {
                b.iter(|| aptp_factor(a.as_ref(), &opts).expect("factor should succeed"));
            });
        }
        group.finish();
    }
}

criterion_group!(factor_solve_benches, bench_factor_solve);
criterion_group!(symbolic_benches, bench_symbolic_analysis);
criterion_group!(mc64_benches, bench_mc64_matching);
criterion_group!(match_order_benches, bench_match_order);
criterion_group!(kernel_benches, bench_kernel_two_level);
criterion_main!(
    factor_solve_benches,
    symbolic_benches,
    mc64_benches,
    match_order_benches,
    kernel_benches
);
