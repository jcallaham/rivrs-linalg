//! SPRAL benchmark suite — runs SPRAL's SSIDS on SuiteSparse test matrices
//! to produce timing and accuracy tables similar to Duff et al (2020) Tables 3 & 4.
//!
//! SPRAL uses its own ordering pipeline (MC64 + METIS) internally.
//! When `--rivrs` is passed, the rivrs solver also runs with its own
//! independent ordering for side-by-side comparison.
//!
//! # Prerequisites
//!
//! Build SPRAL library and driver binaries:
//! ```sh
//! comparisons/drivers/build_spral.sh
//! ```
//!
//! # Usage
//!
//! ```sh
//! # SPRAL-only on CI subset
//! cargo run --bin spral-comparison --release -- --ci-only
//!
//! # Side-by-side comparison with rivrs
//! cargo run --bin spral-comparison --release -- --ci-only --rivrs
//!
//! # Set thread count for both SPRAL (OMP_NUM_THREADS) and rivrs (Par::rayon)
//! cargo run --bin spral-comparison --release -- --ci-only --threads 4
//!
//! # Compare against a previously collected rivrs baseline
//! cargo run --bin spral-comparison --release -- --ci-only \
//!   --compare target/benchmarks/baselines/baseline-latest.json
//!
//! # Filter by category
//! cargo run --bin spral-comparison --release -- --category hard-indefinite
//! ```

#[path = "common.rs"]
mod common;

use std::path::Path;

use rivrs_sparse::io::registry;

use common::extract_lower_triangle;

const DEFAULT_SPRAL_BIN: &str = "/tmp/spral_benchmark";

// ---------------------------------------------------------------------------
// SPRAL-specific stats extracted from ExternalSolverResult.extra
// ---------------------------------------------------------------------------

/// SPRAL-specific statistics beyond the common analyse/factor/solve/backward_error.
/// Additional SPRAL fields (num_neg, factor_flag, analyse_flag, etc.) are preserved
/// in `ExternalSolverResult.extra` and serialized to JSON automatically.
struct SpralStats {
    num_delay: i64,
    num_two: i64,
    maxfront: i64,
}

impl SpralStats {
    /// Extract typed SPRAL stats from the generic extra map.
    fn from_extra(extra: &std::collections::BTreeMap<String, serde_json::Value>) -> Self {
        Self {
            num_delay: extra.get("num_delay").and_then(|v| v.as_i64()).unwrap_or(0),
            num_two: extra.get("num_two").and_then(|v| v.as_i64()).unwrap_or(0),
            maxfront: extra.get("maxfront").and_then(|v| v.as_i64()).unwrap_or(0),
        }
    }
}

// ---------------------------------------------------------------------------
// CLI argument parsing
// ---------------------------------------------------------------------------

struct CliArgs {
    common: common::CommonCliArgs,
    spral_binary: String,
}

fn parse_args() -> CliArgs {
    let mut spral_binary = DEFAULT_SPRAL_BIN.to_string();

    let common = common::parse_common_args(|args, i| match args[i].as_str() {
        "--spral-binary" => {
            if let Some(path) = args.get(i + 1) {
                spral_binary = path.clone();
            }
            Ok(1) // consumed one extra arg
        }
        other => {
            eprintln!("Unknown argument: {other}");
            eprintln!(
                "Usage: spral-comparison [--ci-only] [--threads N] [--category CAT] \
                 [--rivrs] [--compare FILE] [--spral-binary PATH]"
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
        spral_binary,
    }
}

// ---------------------------------------------------------------------------
// Table printing (SPRAL-specific columns: delays, 2x2, maxfrt, sl_st, sl_nd)
// ---------------------------------------------------------------------------

fn print_spral_table_header() {
    eprintln!(
        "{:<35} {:>8} {:>10} {:>7} {:>7} {:>7} {:>9} {:>7} {:>6} {:>6}",
        "Matrix", "n", "nnz", "ana_s", "fac_s", "slv_s", "bwd_err", "delays", "2x2", "maxfrt"
    );
    eprintln!("{}", "-".repeat(110));
}

fn print_spral_row(
    name: &str,
    n: usize,
    nnz: usize,
    r: &common::ExternalSolverResult,
    stats: &SpralStats,
) {
    eprintln!(
        "{:<35} {:>8} {:>10} {:>7.3} {:>7.3} {:>7.3} {:>9.1e} {:>7} {:>6} {:>6}",
        name,
        n,
        nnz,
        r.analyse_s,
        r.factor_s,
        r.solve_s,
        r.backward_error,
        stats.num_delay,
        stats.num_two,
        stats.maxfront
    );
}

fn print_comparison_header(threads: Option<usize>) {
    let thread_str = threads
        .map(|t| format!("threads={t}"))
        .unwrap_or_else(|| "threads=default".to_string());
    eprintln!("\n=== Comparison: SPRAL vs rivrs ({thread_str}) ===");
    eprintln!(
        "{:<35} {:>8} {:>10} {:>10} {:>7} {:>10} {:>10} {:>6} {:>6}",
        "Matrix", "n", "spral_fac", "rivrs_fac", "ratio", "spral_be", "rivrs_be", "sl_st", "sl_nd"
    );
    eprintln!("{}", "-".repeat(114));
}

fn print_comparison_row(rec: &common::ComparisonRecord) {
    let (spral_fac, spral_be) = rec
        .solver
        .as_ref()
        .map(|s| (s.factor_s, s.backward_error))
        .unwrap_or((f64::NAN, f64::NAN));
    let (rivrs_fac, rivrs_be, sl_st, sl_nd) = rec
        .rivrs
        .as_ref()
        .map(|r| {
            (
                r.factor_s,
                r.backward_error,
                r.small_leaf_subtrees,
                r.small_leaf_nodes,
            )
        })
        .unwrap_or((f64::NAN, f64::NAN, 0, 0));
    let ratio = if spral_fac > 0.0 {
        rivrs_fac / spral_fac
    } else {
        f64::NAN
    };
    eprintln!(
        "{:<35} {:>8} {:>10.3} {:>10.3} {:>7.2} {:>10.1e} {:>10.1e} {:>6} {:>6}",
        rec.matrix_name, rec.n, spral_fac, rivrs_fac, ratio, spral_be, rivrs_be, sl_st, sl_nd
    );
}

fn print_category_tables(
    results: &[common::ComparisonRecord],
    has_rivrs: bool,
    threads: Option<usize>,
) {
    let categories = ["positive-definite", "easy-indefinite", "hard-indefinite"];

    if has_rivrs {
        print_comparison_header(threads);
        for rec in results {
            print_comparison_row(rec);
        }
    } else {
        for cat in &categories {
            let in_cat: Vec<_> = results.iter().filter(|r| r.category == *cat).collect();
            if in_cat.is_empty() {
                continue;
            }
            eprintln!("\n=== {} ===", cat);
            print_spral_table_header();
            for rec in &in_cat {
                if let Some(ref solver) = rec.solver {
                    let stats = SpralStats::from_extra(&solver.extra);
                    print_spral_row(&rec.matrix_name, rec.n, rec.nnz, solver, &stats);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let cli = parse_args();

    if !Path::new(&cli.spral_binary).exists() {
        eprintln!(
            "Error: SPRAL benchmark driver not found at {}",
            cli.spral_binary
        );
        eprintln!(
            "Run comparisons/drivers/build_spral.sh to build SPRAL and compile driver binaries."
        );
        std::process::exit(1);
    }

    let matrices = common::load_matrix_entries(&cli.common);

    eprintln!(
        "SPRAL benchmark: {} matrices{}{}",
        matrices.len(),
        if cli.common.ci_only {
            " (CI subset)"
        } else {
            ""
        },
        cli.common
            .threads
            .map(|t| format!(", OMP_NUM_THREADS={t}"))
            .unwrap_or_default()
    );

    let mut results = Vec::new();
    let mut pass = 0usize;
    let mut fail = 0usize;

    // Build SPRAL-specific env vars
    let mut env_vars: Vec<(&str, String)> = vec![
        ("OMP_CANCELLATION", "TRUE".to_string()),
        ("OMP_PROC_BIND", "TRUE".to_string()),
    ];
    if let Some(t) = cli.common.threads {
        env_vars.push(("OMP_NUM_THREADS", t.to_string()));
    }

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

        // Extract lower triangle to compute nnz_lower for threshold
        let (lower_ptr, _, _) = extract_lower_triangle(matrix);
        let nnz_lower = (lower_ptr[n] - 1) as usize;

        let input = common::format_spral_input(matrix);
        let use_file_mode = nnz_lower > common::FILE_MODE_NNZ_THRESHOLD;

        // Run SPRAL via shared subprocess infrastructure
        let solver_result = match common::run_solver_subprocess(
            &cli.spral_binary,
            &input,
            "SPRAL_BENCHMARK_BEGIN",
            "SPRAL_BENCHMARK_END",
            &env_vars,
            &[],
            use_file_mode,
            "spral_benchmark",
        ) {
            Ok(r) => Some(r),
            Err(e) => {
                eprintln!("{:<35} SPRAL FAILED: {}", meta.name, e);
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
            let stats = SpralStats::from_extra(&sr.extra);
            if cli.common.rivrs {
                if let Some(ref rr) = rivrs_result {
                    let ratio = if sr.factor_s > 0.0 {
                        rr.factor_s / sr.factor_s
                    } else {
                        f64::NAN
                    };
                    eprintln!(
                        "{:<35} {:>8} spral={:.3}s rivrs={:.3}s ratio={:.2} be={:.1e}/{:.1e}",
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
                    "{:<35} {:>8} fac={:.3}s be={:.1e} delays={} 2x2={} maxfrt={}",
                    meta.name,
                    n,
                    sr.factor_s,
                    sr.backward_error,
                    stats.num_delay,
                    stats.num_two,
                    stats.maxfront
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
    print_category_tables(&results, cli.common.rivrs, cli.common.threads);

    // Print comparison with baseline if requested
    if let Some(ref compare_path) = cli.common.compare {
        match std::fs::read_to_string(compare_path) {
            Ok(json) => match serde_json::from_str::<common::BaselineSuite>(&json) {
                Ok(baseline) => {
                    common::print_baseline_comparison("SPRAL", &results, &baseline);
                }
                Err(e) => eprintln!("\nFailed to parse baseline JSON: {e}"),
            },
            Err(e) => eprintln!("\nFailed to read baseline file: {e}"),
        }
    }

    eprintln!("\n{pass}/{} completed successfully", pass + fail);

    // Write JSON output via shared infrastructure
    let mode = if cli.common.rivrs {
        "comparison"
    } else {
        "spral-only"
    };

    let suite = common::BenchmarkSuite {
        timestamp: common::timestamp(),
        platform: format!("{} {}", std::env::consts::OS, std::env::consts::ARCH),
        solver_name: "SPRAL".to_string(),
        solver_version: "SPRAL-git".to_string(),
        rivrs_version: common::git_hash(),
        threads: cli.common.threads,
        mode: mode.to_string(),
        results,
    };

    let out_path = common::write_suite_json(&suite, "target/benchmarks/spral");
    eprintln!("JSON written to: {}", out_path.display());
}
