//! MA27 benchmark suite — runs HSL MA27 on SuiteSparse test matrices
//! for timing and accuracy comparison against rivrs-sparse.
//!
//! MA27 (Duff & Reid, 1982) is the classic multifrontal solver. It uses
//! its own built-in minimum degree ordering (not METIS). This means the
//! ordering quality differs from SPRAL/rivrs (which use METIS-based orderings),
//! and timing comparisons should note this difference.
//!
//! # Prerequisites
//!
//! MA27 source is NOT redistributable. Obtain from HSL:
//!   https://www.hsl.rl.ac.uk/catalogue/ma27.html
//!
//! Then build the driver:
//! ```sh
//! MA27_SRC=/path/to/ma27 comparisons/drivers/build_ma27.sh
//! ```
//!
//! # Usage
//!
//! ```sh
//! # MA27-only on CI subset
//! cargo run --bin ma27-comparison --release -- --ci-only
//!
//! # Side-by-side comparison with rivrs
//! cargo run --bin ma27-comparison --release -- --ci-only --rivrs
//!
//! # Compare against a previously collected rivrs baseline
//! cargo run --bin ma27-comparison --release -- --ci-only \
//!   --compare target/benchmarks/baselines/baseline-latest.json
//!
//! # Filter by category
//! cargo run --bin ma27-comparison --release -- --category hard-indefinite
//! ```

#[path = "common.rs"]
mod common;

use std::path::Path;

use rivrs_sparse::io::registry;

const DEFAULT_MA27_BIN: &str = "/tmp/ma27_benchmark";

// ---------------------------------------------------------------------------
// CLI argument parsing
// ---------------------------------------------------------------------------

struct CliArgs {
    common: common::CommonCliArgs,
    ma27_binary: String,
}

fn parse_args() -> CliArgs {
    let mut ma27_binary = DEFAULT_MA27_BIN.to_string();

    let common = common::parse_common_args(|args, i| match args[i].as_str() {
        "--ma27-binary" => {
            if let Some(path) = args.get(i + 1) {
                ma27_binary = path.clone();
            }
            Ok(1)
        }
        other => {
            eprintln!("Unknown argument: {other}");
            eprintln!(
                "Usage: ma27-comparison [--ci-only] [--threads N] [--category CAT] \
                 [--rivrs] [--compare FILE] [--ma27-binary PATH]"
            );
            std::process::exit(1);
        }
    })
    .unwrap_or_else(|e| {
        eprintln!("Error: {e}");
        std::process::exit(1);
    });

    CliArgs {
        common,
        ma27_binary,
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let cli = parse_args();

    if !Path::new(&cli.ma27_binary).exists() {
        eprintln!(
            "Error: MA27 benchmark driver not found at {}",
            cli.ma27_binary
        );
        eprintln!("MA27 source is not redistributable. To set up:");
        eprintln!("  1. Obtain MA27 from https://www.hsl.rl.ac.uk/catalogue/ma27.html");
        eprintln!("  2. Run: MA27_SRC=/path/to/ma27 comparisons/drivers/build_ma27.sh");
        std::process::exit(1);
    }

    let matrices = common::load_matrix_entries(&cli.common);

    eprintln!(
        "MA27 benchmark: {} matrices{}",
        matrices.len(),
        if cli.common.ci_only {
            " (CI subset)"
        } else {
            ""
        },
    );

    let mut results = Vec::new();
    let mut pass = 0usize;
    let mut fail = 0usize;

    for meta in &matrices {
        let test = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(t)) => t,
            Ok(None) => {
                eprintln!("{:<35} SKIP (not found on disk)", meta.name);
                continue;
            }
            Err(e) => {
                eprintln!("{:<35} LOAD FAILED: {}", meta.name, e);
                fail += 1;
                continue;
            }
        };

        let matrix = &test.matrix;
        let n = matrix.nrows();

        // Format matrix as COO for MA27
        let input = common::format_lower_coo_text(matrix);
        let nnz_lower = input.lines().count().saturating_sub(1);

        let use_file_mode = nnz_lower > common::FILE_MODE_NNZ_THRESHOLD;

        // Run MA27
        let ma27_result = common::run_solver_subprocess(
            &cli.ma27_binary,
            &input,
            "MA27_BENCHMARK_BEGIN",
            "MA27_BENCHMARK_END",
            &[],
            &[],
            use_file_mode,
            "ma27_benchmark",
        );

        let solver_result = match ma27_result {
            Ok(r) => Some(r),
            Err(e) => {
                eprintln!("{:<35} MA27 FAILED: {}", meta.name, e);
                fail += 1;
                None
            }
        };

        // Run rivrs if requested
        let rivrs_result = if cli.common.rivrs {
            let r = match common::run_rivrs_solver(matrix, &meta.name, cli.common.threads) {
                Ok(r) => Some(r),
                Err(e) => {
                    eprintln!("{:<35} RIVRS FAILED: {}", meta.name, e);
                    None
                }
            };
            common::maybe_trim_memory();
            r
        } else {
            None
        };

        // Print progress row
        if let Some(ref sr) = solver_result {
            if cli.common.rivrs {
                if let Some(ref rr) = rivrs_result {
                    let ratio = if sr.factor_s > 0.0 {
                        rr.factor_s / sr.factor_s
                    } else {
                        f64::NAN
                    };
                    eprintln!(
                        "{:<35} {:>8} ma27={:.3}s rivrs={:.3}s ratio={:.2} be={:.1e}/{:.1e}",
                        meta.name,
                        n,
                        sr.factor_s,
                        rr.factor_s,
                        ratio,
                        sr.backward_error,
                        rr.backward_error
                    );
                }
            } else {
                eprintln!(
                    "{:<35} {:>8} fac={:.3}s be={:.1e}",
                    meta.name, n, sr.factor_s, sr.backward_error
                );
            }
            pass += 1;
        }

        results.push(common::ComparisonRecord {
            matrix_name: meta.name.clone(),
            category: meta.category.clone(),
            n,
            nnz: meta.nnz,
            solver: solver_result,
            rivrs: rivrs_result,
        });
    }

    // Print summary tables
    let categories = ["positive-definite", "easy-indefinite", "hard-indefinite"];

    if cli.common.rivrs {
        common::print_solver_comparison_header("MA27", cli.common.threads);
        for rec in &results {
            if let Some(ref sr) = rec.solver {
                common::print_solver_comparison_row(
                    &rec.matrix_name,
                    rec.n,
                    sr.factor_s,
                    sr.backward_error,
                    rec.rivrs.as_ref(),
                );
            }
        }
    } else {
        for cat in &categories {
            let in_cat: Vec<_> = results.iter().filter(|r| r.category == *cat).collect();
            if in_cat.is_empty() {
                continue;
            }
            eprintln!("\n=== {} ===", cat);
            common::print_solver_only_header("MA27");
            for rec in &in_cat {
                if let Some(ref sr) = rec.solver {
                    common::print_solver_only_row(
                        &rec.matrix_name,
                        rec.n,
                        rec.nnz,
                        sr.analyse_s,
                        sr.factor_s,
                        sr.solve_s,
                        sr.backward_error,
                    );
                }
            }
        }
    }

    // Print comparison with baseline if requested
    if let Some(ref compare_path) = cli.common.compare {
        match std::fs::read_to_string(compare_path) {
            Ok(json) => match serde_json::from_str::<common::BaselineSuite>(&json) {
                Ok(baseline) => {
                    common::print_baseline_comparison("MA27", &results, &baseline);
                }
                Err(e) => eprintln!("\nFailed to parse baseline JSON: {e}"),
            },
            Err(e) => eprintln!("\nFailed to read baseline file: {e}"),
        }
    }

    eprintln!("\n{pass}/{} completed successfully", pass + fail);

    // Write JSON output
    let mode = if cli.common.rivrs {
        "comparison"
    } else {
        "ma27-only"
    };

    let suite = common::BenchmarkSuite {
        timestamp: common::timestamp(),
        platform: format!("{} {}", std::env::consts::OS, std::env::consts::ARCH),
        solver_name: "MA27".to_string(),
        solver_version: "HSL-MA27".to_string(),
        rivrs_version: common::git_hash(),
        threads: cli.common.threads,
        mode: mode.to_string(),
        results,
    };

    let out_path = common::write_suite_json(&suite, "target/benchmarks/ma27");
    eprintln!("JSON written to: {}", out_path.display());
}
