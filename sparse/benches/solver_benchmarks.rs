//! Criterion benchmark binary for SSIDS solver phases.
//!
//! Benchmarks individual solver phases (analyze, factor, solve) and the
//! end-to-end pipeline (roundtrip) across test matrices using the
//! `Benchmarkable` trait.

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

use rivrs_sparse::benchmarking::traits::Benchmarkable;
use rivrs_sparse::benchmarking::{
    BenchmarkConfig, BenchmarkPhase, MockBenchmarkable, SkippedBenchmark, read_peak_rss_kb,
};
use rivrs_sparse::testing::{TestCaseFilter, load_test_cases};

fn run_component_benchmarks(
    c: &mut Criterion,
    config: &BenchmarkConfig,
    solver: &dyn Benchmarkable,
) {
    let cases = match load_test_cases(&config.filter) {
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

    let mut skipped: Vec<SkippedBenchmark> = Vec::new();

    for phase in &config.phases {
        if *phase == BenchmarkPhase::Roundtrip {
            continue; // Roundtrip handled by run_e2e_benchmarks
        }

        let group_name = format!("ssids/{}", phase);
        let mut group = c.benchmark_group(&group_name);

        if let Some(sample_size) = config.sample_size {
            group.sample_size(sample_size);
        }
        if let Some(measurement_time) = config.measurement_time {
            group.measurement_time(measurement_time);
        }
        if let Some(warm_up_time) = config.warm_up_time {
            group.warm_up_time(warm_up_time);
        }

        for case in &cases {
            let nnz = case.properties.nnz;
            group.throughput(Throughput::Elements(nnz as u64));

            // Check if the phase is implemented by the solver
            let phase_available = match phase {
                BenchmarkPhase::Analyze => solver.bench_analyze(&case.matrix).is_some(),
                BenchmarkPhase::Factor => solver.bench_factor(&case.matrix, None).is_some(),
                BenchmarkPhase::Solve => {
                    let rhs = vec![1.0; case.matrix.nrows()];
                    solver.bench_solve(&case.matrix, None, &rhs).is_some()
                }
                BenchmarkPhase::Roundtrip => unreachable!(),
            };

            if !phase_available {
                skipped.push(SkippedBenchmark {
                    matrix_name: case.name.clone(),
                    phase: *phase,
                    reason: "phase not implemented".to_string(),
                });
                continue;
            }

            group.bench_with_input(BenchmarkId::from_parameter(&case.name), &case, |b, case| {
                match phase {
                    BenchmarkPhase::Analyze => {
                        b.iter(|| solver.bench_analyze(&case.matrix));
                    }
                    BenchmarkPhase::Factor => {
                        let analysis = solver.bench_analyze(&case.matrix);
                        b.iter(|| solver.bench_factor(&case.matrix, analysis.as_deref()));
                    }
                    BenchmarkPhase::Solve => {
                        let analysis = solver.bench_analyze(&case.matrix);
                        let factorization = solver.bench_factor(&case.matrix, analysis.as_deref());
                        let rhs: Vec<f64> = vec![1.0; case.matrix.nrows()];
                        b.iter(|| solver.bench_solve(&case.matrix, factorization.as_deref(), &rhs));
                    }
                    BenchmarkPhase::Roundtrip => unreachable!(),
                }
            });
        }
        group.finish();
    }

    if !skipped.is_empty() {
        eprintln!("\nSkipped benchmarks:");
        for s in &skipped {
            eprintln!("  {}", s);
        }
    }
}

fn run_e2e_benchmarks(c: &mut Criterion, config: &BenchmarkConfig, solver: &dyn Benchmarkable) {
    let cases = match load_test_cases(&config.filter) {
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

    let mut group = c.benchmark_group("ssids/roundtrip");

    if let Some(sample_size) = config.sample_size {
        group.sample_size(sample_size);
    }
    if let Some(measurement_time) = config.measurement_time {
        group.measurement_time(measurement_time);
    }
    if let Some(warm_up_time) = config.warm_up_time {
        group.warm_up_time(warm_up_time);
    }

    for case in &cases {
        let nnz = case.properties.nnz;
        group.throughput(Throughput::Elements(nnz as u64));

        let rhs: Vec<f64> = vec![1.0; case.matrix.nrows()];

        // Check if roundtrip is available
        let available = solver.bench_roundtrip(&case.matrix, &rhs).is_some();
        if !available {
            eprintln!("  SKIP {}/roundtrip: phase not implemented", case.name);
            continue;
        }

        group.bench_with_input(BenchmarkId::from_parameter(&case.name), &case, |b, case| {
            let rhs: Vec<f64> = vec![1.0; case.matrix.nrows()];
            b.iter(|| solver.bench_roundtrip(&case.matrix, &rhs));
        });
    }
    group.finish();
}

fn bench_components(c: &mut Criterion) {
    let rss_before = read_peak_rss_kb();

    let config = BenchmarkConfig::default_components();
    let solver = MockBenchmarkable::new();
    run_component_benchmarks(c, &config, &solver);

    let rss_after = read_peak_rss_kb();
    if let (Some(before), Some(after)) = (rss_before, rss_after) {
        eprintln!(
            "\nPeak RSS: {} KB -> {} KB (delta: {} KB)",
            before,
            after,
            after.saturating_sub(before)
        );
    } else if let Some(after) = rss_after {
        eprintln!("\nPeak RSS: {} KB", after);
    }
}

fn bench_e2e(c: &mut Criterion) {
    let config = BenchmarkConfig::default_components().with_phases(vec![BenchmarkPhase::Roundtrip]);
    let solver = MockBenchmarkable::new();
    run_e2e_benchmarks(c, &config, &solver);
}

fn bench_ci_subset(c: &mut Criterion) {
    let config = BenchmarkConfig::default_components().with_filter(TestCaseFilter::ci_subset());
    let solver = MockBenchmarkable::new();
    run_component_benchmarks(c, &config, &solver);
}

criterion_group!(component_benches, bench_components);
criterion_group!(e2e_benches, bench_e2e);
criterion_group!(ci_benches, bench_ci_subset);
criterion_main!(component_benches, e2e_benches, ci_benches);
