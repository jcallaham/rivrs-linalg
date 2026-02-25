//! Parallel scaling benchmark for the SSIDS solver.
//!
//! Measures factorization and solve performance across multiple thread counts
//! and produces a structured report showing speedup and efficiency per matrix
//! and per workload class (IntraNode, TreeLevel, Mixed).
//!
//! Usage:
//!   cargo run --example parallel_scaling --release -- --ci-only
//!   cargo run --example parallel_scaling --release -- --ci-only --threads 1,2,4
//!   cargo run --example parallel_scaling --release

use std::time::Instant;

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::{Col, Par};
use serde::Serialize;

use rivrs_sparse::aptp::{AnalyzeOptions, FactorOptions, SparseLDLT};
use rivrs_sparse::io::registry;
use rivrs_sparse::validate::sparse_backward_error;

/// Per-thread-count timing for a single matrix.
#[derive(Debug, Clone, Serialize)]
struct ThreadTiming {
    threads: usize,
    factor_ms: f64,
    solve_ms: f64,
    backward_error: f64,
    factor_speedup: f64,
    solve_speedup: f64,
    factor_efficiency: f64,
    solve_efficiency: f64,
}

/// Per-matrix scaling result.
#[derive(Debug, Clone, Serialize)]
struct MatrixScaling {
    name: String,
    n: usize,
    nnz: usize,
    max_front_size: usize,
    num_supernodes: usize,
    timings: Vec<ThreadTiming>,
}

/// Overall scaling report.
#[derive(Debug, Clone, Serialize)]
struct ScalingReport {
    timestamp: String,
    thread_counts: Vec<usize>,
    matrices: Vec<MatrixScaling>,
}

fn parse_args() -> (bool, Vec<usize>) {
    let args: Vec<String> = std::env::args().collect();
    let ci_only = args.iter().any(|a| a == "--ci-only");

    let threads = if let Some(pos) = args.iter().position(|a| a == "--threads") {
        args.get(pos + 1)
            .expect("--threads requires a comma-separated list")
            .split(',')
            .map(|s| s.trim().parse::<usize>().expect("invalid thread count"))
            .collect()
    } else {
        vec![1, 2, 4]
    };

    (ci_only, threads)
}

fn main() {
    let (ci_only, thread_counts) = parse_args();

    eprintln!("Parallel Scaling Benchmark");
    eprintln!("Thread counts: {:?}", thread_counts);

    let reg = registry::load_registry().expect("load registry");
    let matrices: Vec<_> = if ci_only {
        reg.iter().filter(|m| m.ci_subset).collect()
    } else {
        reg.iter().collect()
    };

    eprintln!("Matrices: {} (ci_only={})", matrices.len(), ci_only);
    eprintln!();

    let analyze_opts = AnalyzeOptions::default();
    let mut results = Vec::new();

    for meta in &matrices {
        let test = match registry::load_test_matrix(&meta.name) {
            Ok(Some(t)) => t,
            Ok(None) => {
                eprintln!("  SKIP {} (not found)", meta.name);
                continue;
            }
            Err(e) => {
                eprintln!("  SKIP {} ({})", meta.name, e);
                continue;
            }
        };
        let matrix = &test.matrix;
        let n = matrix.nrows();
        let nnz = matrix.compute_nnz();

        eprint!("  {} (n={}, nnz={}) ", meta.name, n, nnz);

        let mut timings = Vec::new();
        let mut base_factor_ms = 0.0;
        let mut base_solve_ms = 0.0;
        let mut max_front = 0;
        let mut num_supernodes = 0;

        let rss_before = rivrs_sparse::benchmarking::read_peak_rss_kb();
        for &nthreads in &thread_counts {
            let par = if nthreads <= 1 {
                Par::Seq
            } else {
                Par::rayon(nthreads)
            };
            let factor_opts = FactorOptions {
                par,
                ..FactorOptions::default()
            };

            let mut solver =
                SparseLDLT::analyze_with_matrix(matrix, &analyze_opts).expect("analyze");

            // Factor (3 runs, take median)
            let mut factor_times = Vec::new();
            for _ in 0..3 {
                let start = Instant::now();
                solver.factor(matrix, &factor_opts).expect("factor");
                factor_times.push(start.elapsed().as_secs_f64() * 1000.0);
            }
            factor_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let factor_ms = factor_times[1]; // median

            if let Some(stats) = solver.stats() {
                max_front = stats.max_front_size;
            }
            if let Some(sn_stats) = solver.per_supernode_stats() {
                num_supernodes = sn_stats.len();
            }

            // Solve
            let b = Col::from_fn(n, |i| ((i + 1) as f64).sin());
            let scratch_req = solver.solve_scratch(1);

            let mut solve_times = Vec::new();
            for _ in 0..3 {
                let mut mem = MemBuffer::new(scratch_req);
                let start = Instant::now();
                solver
                    .solve(&b, MemStack::new(&mut mem), par)
                    .expect("solve");
                solve_times.push(start.elapsed().as_secs_f64() * 1000.0);
            }
            solve_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let solve_ms = solve_times[1]; // median

            // Backward error (last solve result)
            let mut mem = MemBuffer::new(scratch_req);
            let x = solver
                .solve(&b, MemStack::new(&mut mem), par)
                .expect("solve");
            let be = sparse_backward_error(matrix, &x, &b);

            if nthreads == thread_counts[0] || (nthreads == 1) {
                base_factor_ms = factor_ms;
                base_solve_ms = solve_ms;
            }

            let factor_speedup = if base_factor_ms > 0.0 {
                base_factor_ms / factor_ms
            } else {
                1.0
            };
            let solve_speedup = if base_solve_ms > 0.0 {
                base_solve_ms / solve_ms
            } else {
                1.0
            };

            let effective_threads = nthreads.max(1) as f64;
            timings.push(ThreadTiming {
                threads: nthreads,
                factor_ms,
                solve_ms,
                backward_error: be,
                factor_speedup,
                solve_speedup,
                factor_efficiency: factor_speedup / effective_threads,
                solve_efficiency: solve_speedup / effective_threads,
            });

            // Log current RSS before dropping solver
            let rss_with_solver = rivrs_sparse::benchmarking::read_current_rss_kb().unwrap_or(0);
            eprint!("T{}:{:.1}ms ", nthreads, factor_ms);

            // Drop solver explicitly and force glibc to return freed pages to OS.
            // Without malloc_trim, glibc holds freed pages in its arena — for H2O
            // (5.9 GB peak), only ~1.3 GB is returned on drop, leaving 4.5 GB of
            // unreturned-but-freed memory that causes OOM on the next thread count.
            drop(solver);
            #[cfg(target_os = "linux")]
            {
                unsafe extern "C" {
                    fn malloc_trim(pad: usize) -> i32;
                }
                unsafe {
                    malloc_trim(0);
                }
            }
            let rss_after_drop = rivrs_sparse::benchmarking::read_current_rss_kb().unwrap_or(0);
            if rss_with_solver > 100_000 {
                eprint!(
                    "[rss:{}→{}MB] ",
                    rss_with_solver / 1024,
                    rss_after_drop / 1024
                );
            }
        }
        let rss_end = rivrs_sparse::benchmarking::read_current_rss_kb().unwrap_or(0);
        if let Some(before) = rss_before {
            if rss_end > 100_000 || before > 100_000 {
                eprint!("[matrix rss:{}→{}MB] ", before / 1024, rss_end / 1024);
            }
        }
        eprintln!();

        results.push(MatrixScaling {
            name: meta.name.clone(),
            n,
            nnz,
            max_front_size: max_front,
            num_supernodes,
            timings,
        });
    }

    // Human-readable summary table
    eprintln!();
    eprintln!("{:<25} {:>8} {:>8} {:>8}", "Matrix", "n", "max_front", "sn");
    for tc in &thread_counts {
        eprint!("  T{:>2}_ms spdup", tc);
    }
    eprintln!();
    eprintln!("{}", "-".repeat(25 + 8 + 8 + 8 + thread_counts.len() * 14));

    for ms in &results {
        eprint!(
            "{:<25} {:>8} {:>8} {:>8}",
            ms.name, ms.n, ms.max_front_size, ms.num_supernodes
        );
        for t in &ms.timings {
            eprint!("  {:>6.1} {:>5.2}x", t.factor_ms, t.factor_speedup);
        }
        eprintln!();
    }

    // JSON output to stdout
    let report = ScalingReport {
        timestamp: chrono_lite_timestamp(),
        thread_counts: thread_counts.clone(),
        matrices: results,
    };

    // Write JSON to file
    let out_dir = std::path::Path::new("target/benchmarks/parallel");
    std::fs::create_dir_all(out_dir).expect("create output dir");
    let filename = format!("scaling-{}.json", chrono_lite_timestamp());
    let path = out_dir.join(&filename);
    let json = serde_json::to_string_pretty(&report).expect("serialize");
    std::fs::write(&path, &json).expect("write json");
    eprintln!();
    eprintln!("JSON output: {}", path.display());
}

/// Simple timestamp without chrono dependency.
fn chrono_lite_timestamp() -> String {
    use std::time::SystemTime;
    let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap();
    format!("{}", d.as_secs())
}
