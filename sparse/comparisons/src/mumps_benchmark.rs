//! MUMPS benchmark suite — runs MUMPS on SuiteSparse test matrices
//! for timing and accuracy comparison against rivrs-sparse.
//!
//! MUMPS uses its own ordering (auto-select by default, configurable
//! via `MUMPS_ORDERING` environment variable).
//!
//! # Prerequisites
//!
//! Build the MUMPS benchmark driver:
//! ```sh
//! comparisons/drivers/build_mumps.sh
//! ```
//!
//! # Usage
//!
//! ```sh
//! # MUMPS-only on CI subset
//! cargo run --bin mumps-comparison --release -- --ci-only
//!
//! # Side-by-side comparison with rivrs
//! cargo run --bin mumps-comparison --release -- --ci-only --rivrs
//!
//! # Set thread count for rivrs (MUMPS sequential version is single-threaded)
//! cargo run --bin mumps-comparison --release -- --ci-only --rivrs --threads 4
//!
//! # Compare against a previously collected rivrs baseline
//! cargo run --bin mumps-comparison --release -- --ci-only \
//!   --compare target/benchmarks/baselines/baseline-latest.json
//!
//! # Filter by category
//! cargo run --bin mumps-comparison --release -- --category hard-indefinite
//! ```

#[path = "common.rs"]
mod common;

use std::path::Path;

use rivrs_sparse::io::registry;

const DEFAULT_MUMPS_BIN: &str = "/tmp/mumps_benchmark";

// ---------------------------------------------------------------------------
// CLI argument parsing
// ---------------------------------------------------------------------------

struct CliArgs {
    common: common::CommonCliArgs,
    mumps_binary: String,
}

fn parse_args() -> CliArgs {
    let mut mumps_binary = DEFAULT_MUMPS_BIN.to_string();

    let common = common::parse_common_args(|args, i| match args[i].as_str() {
        "--mumps-binary" => {
            if let Some(path) = args.get(i + 1) {
                mumps_binary = path.clone();
            }
            Ok(1)
        }
        other => {
            eprintln!("Unknown argument: {other}");
            eprintln!(
                "Usage: mumps-comparison [--ci-only] [--threads N] [--category CAT] \
                 [--rivrs] [--compare FILE] [--mumps-binary PATH]"
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
        mumps_binary,
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let cli = parse_args();

    if !Path::new(&cli.mumps_binary).exists() {
        eprintln!(
            "Error: MUMPS benchmark driver not found at {}",
            cli.mumps_binary
        );
        eprintln!("Run comparisons/drivers/build_mumps.sh to compile the MUMPS driver.");
        std::process::exit(1);
    }

    let matrices = common::load_matrix_entries(&cli.common);

    eprintln!(
        "MUMPS benchmark: {} matrices{}{}",
        matrices.len(),
        if cli.common.ci_only {
            " (CI subset)"
        } else {
            ""
        },
        cli.common
            .threads
            .map(|t| format!(", threads={t}"))
            .unwrap_or_default()
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

        // Format matrix as COO for MUMPS
        let input = common::format_lower_coo_text(matrix);
        let nnz_lower = input.lines().count().saturating_sub(1); // subtract header line

        let use_file_mode = nnz_lower > common::FILE_MODE_NNZ_THRESHOLD;

        // Run MUMPS
        let mumps_result = common::run_solver_subprocess(
            &cli.mumps_binary,
            &input,
            "MUMPS_BENCHMARK_BEGIN",
            "MUMPS_BENCHMARK_END",
            &[], // MUMPS sequential doesn't use OMP
            &[],
            use_file_mode,
            "mumps_benchmark",
        );

        let solver_result = match mumps_result {
            Ok(r) => Some(r),
            Err(e) => {
                eprintln!("{:<35} MUMPS FAILED: {}", meta.name, e);
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
                        "{:<35} {:>8} mumps={:.3}s rivrs={:.3}s ratio={:.2} be={:.1e}/{:.1e}",
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
        common::print_solver_comparison_header("MUMPS", cli.common.threads);
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
            common::print_solver_only_header("MUMPS");
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
                    common::print_baseline_comparison("MUMPS", &results, &baseline);
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
        "mumps-only"
    };

    let suite = common::BenchmarkSuite {
        timestamp: common::timestamp(),
        platform: format!("{} {}", std::env::consts::OS, std::env::consts::ARCH),
        solver_name: "MUMPS".to_string(),
        solver_version: "MUMPS-seq".to_string(),
        rivrs_version: common::git_hash(),
        threads: cli.common.threads,
        mode: mode.to_string(),
        results,
    };

    let out_path = common::write_suite_json(&suite, "target/benchmarks/mumps");
    eprintln!("JSON written to: {}", out_path.display());
}
