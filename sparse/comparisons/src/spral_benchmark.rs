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

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use faer::sparse::SparseColMat;
use serde::{Deserialize, Serialize};

use rivrs_sparse::io::registry;

use common::{RivrsResult, extract_lower_triangle, format_spral_input, write_temp_matrix_file};

const DEFAULT_SPRAL_BIN: &str = "/tmp/spral_benchmark";

// ---------------------------------------------------------------------------
// SPRAL-specific data types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SpralResult {
    matrix_name: String,
    category: String,
    n: usize,
    nnz: usize,
    analyse_s: f64,
    factor_s: f64,
    solve_s: f64,
    total_s: f64,
    backward_error: f64,
    num_delay: i64,
    num_two: i64,
    num_neg: i64,
    maxfront: i64,
    maxsupernode: i64,
    num_sup: i64,
    num_factor: i64,
    num_flops: i64,
    not_first_pass: i64,
    not_second_pass: i64,
    matrix_rank: i64,
    status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ComparisonRecord {
    matrix_name: String,
    category: String,
    n: usize,
    nnz: usize,
    spral: Option<SpralResult>,
    rivrs: Option<RivrsResult>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SpralBenchmarkSuite {
    timestamp: String,
    platform: String,
    omp_num_threads: Option<usize>,
    spral_version: String,
    rivrs_version: String,
    mode: String,
    results: Vec<ComparisonRecord>,
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
// SPRAL subprocess invocation
// ---------------------------------------------------------------------------

fn run_spral_benchmark(
    matrix: &SparseColMat<usize, f64>,
    spral_bin: &str,
    threads: Option<usize>,
    nnz_lower: usize,
) -> Result<SpralParsed, String> {
    let use_file_mode = nnz_lower > common::FILE_MODE_NNZ_THRESHOLD;

    let mut cmd = Command::new(spral_bin);
    cmd.stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .env("OMP_CANCELLATION", "TRUE")
        .env("OMP_PROC_BIND", "TRUE");

    if let Some(t) = threads {
        cmd.env("OMP_NUM_THREADS", t.to_string());
    }

    if use_file_mode {
        let input = format_spral_input(matrix);
        let temp_path = write_temp_matrix_file(&input, "spral_benchmark")
            .map_err(|e| format!("temp file: {e}"))?;
        cmd.arg(
            temp_path
                .to_str()
                .unwrap_or("/tmp/spral_benchmark_input/matrix.txt"),
        );

        let output = cmd
            .spawn()
            .map_err(|e| format!("spawn: {e}"))?
            .wait_with_output()
            .map_err(|e| format!("wait: {e}"))?;

        let _ = std::fs::remove_file(&temp_path);
        parse_spral_result(&output)
    } else {
        cmd.stdin(std::process::Stdio::piped());

        let mut child = cmd.spawn().map_err(|e| format!("spawn: {e}"))?;

        let input = format_spral_input(matrix);
        if let Some(stdin) = child.stdin.as_mut() {
            stdin
                .write_all(input.as_bytes())
                .map_err(|e| format!("stdin write: {e}"))?;
        }
        drop(child.stdin.take());

        let output = child.wait_with_output().map_err(|e| format!("wait: {e}"))?;
        parse_spral_result(&output)
    }
}

struct SpralParsed {
    analyse_s: f64,
    factor_s: f64,
    solve_s: f64,
    backward_error: f64,
    num_delay: i64,
    num_two: i64,
    num_neg: i64,
    maxfront: i64,
    maxsupernode: i64,
    num_sup: i64,
    num_factor: i64,
    num_flops: i64,
    not_first_pass: i64,
    not_second_pass: i64,
    matrix_rank: i64,
    analyse_flag: i64,
    factor_flag: i64,
    solve_flag: i64,
}

fn parse_spral_result(output: &std::process::Output) -> Result<SpralParsed, String> {
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("exit {}: {}", output.status, stderr));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);

    let begin_marker = "SPRAL_BENCHMARK_BEGIN";
    let end_marker = "SPRAL_BENCHMARK_END";

    let begin = stdout
        .find(begin_marker)
        .ok_or("missing SPRAL_BENCHMARK_BEGIN")?;
    let end = stdout
        .find(end_marker)
        .ok_or("missing SPRAL_BENCHMARK_END")?;

    let block = &stdout[begin + begin_marker.len()..end];

    let mut result = SpralParsed {
        analyse_s: 0.0,
        factor_s: 0.0,
        solve_s: 0.0,
        backward_error: f64::NAN,
        num_delay: 0,
        num_two: 0,
        num_neg: 0,
        maxfront: 0,
        maxsupernode: 0,
        num_sup: 0,
        num_factor: 0,
        num_flops: 0,
        not_first_pass: 0,
        not_second_pass: 0,
        matrix_rank: 0,
        analyse_flag: 0,
        factor_flag: 0,
        solve_flag: 0,
    };

    for line in block.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut parts = line.splitn(2, char::is_whitespace);
        let key = match parts.next() {
            Some(k) => k,
            None => continue,
        };
        let val = match parts.next() {
            Some(v) => v.trim(),
            None => continue,
        };

        match key {
            "analyse_s" => result.analyse_s = val.parse().unwrap_or(0.0),
            "factor_s" => result.factor_s = val.parse().unwrap_or(0.0),
            "solve_s" => result.solve_s = val.parse().unwrap_or(0.0),
            "backward_error" => result.backward_error = val.parse().unwrap_or(f64::NAN),
            "num_delay" => result.num_delay = val.parse().unwrap_or(0),
            "num_two" => result.num_two = val.parse().unwrap_or(0),
            "num_neg" => result.num_neg = val.parse().unwrap_or(0),
            "maxfront" => result.maxfront = val.parse().unwrap_or(0),
            "maxsupernode" => result.maxsupernode = val.parse().unwrap_or(0),
            "num_sup" => result.num_sup = val.parse().unwrap_or(0),
            "num_factor" => result.num_factor = val.parse().unwrap_or(0),
            "num_flops" => result.num_flops = val.parse().unwrap_or(0),
            "not_first_pass" => result.not_first_pass = val.parse().unwrap_or(0),
            "not_second_pass" => result.not_second_pass = val.parse().unwrap_or(0),
            "matrix_rank" => result.matrix_rank = val.parse().unwrap_or(0),
            "analyse_flag" => result.analyse_flag = val.parse().unwrap_or(0),
            "factor_flag" => result.factor_flag = val.parse().unwrap_or(0),
            "solve_flag" => result.solve_flag = val.parse().unwrap_or(0),
            _ => {}
        }
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// Table printing
// ---------------------------------------------------------------------------

fn print_spral_table_header() {
    eprintln!(
        "{:<35} {:>8} {:>10} {:>7} {:>7} {:>7} {:>9} {:>7} {:>6} {:>6}",
        "Matrix", "n", "nnz", "ana_s", "fac_s", "slv_s", "bwd_err", "delays", "2x2", "maxfrt"
    );
    eprintln!("{}", "-".repeat(110));
}

fn print_spral_row(name: &str, n: usize, nnz: usize, r: &SpralResult) {
    eprintln!(
        "{:<35} {:>8} {:>10} {:>7.3} {:>7.3} {:>7.3} {:>9.1e} {:>7} {:>6} {:>6}",
        name,
        n,
        nnz,
        r.analyse_s,
        r.factor_s,
        r.solve_s,
        r.backward_error,
        r.num_delay,
        r.num_two,
        r.maxfront
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

fn print_comparison_row(rec: &ComparisonRecord) {
    let (spral_fac, spral_be) = rec
        .spral
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

fn print_category_tables(results: &[ComparisonRecord], has_rivrs: bool, threads: Option<usize>) {
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
                if let Some(ref spral) = rec.spral {
                    print_spral_row(&rec.matrix_name, rec.n, rec.nnz, spral);
                }
            }
        }
    }
}

fn print_baseline_comparison(results: &[ComparisonRecord], baseline: &common::BaselineSuite) {
    eprintln!("\n=== Comparison: SPRAL vs rivrs baseline ===");
    eprintln!(
        "{:<35} {:>8} {:>10} {:>10} {:>7} {:>10} {:>10}",
        "Matrix", "n", "spral_fac", "rivrs_fac", "ratio", "spral_be", "rivrs_be"
    );
    eprintln!("{}", "-".repeat(100));

    for rec in results {
        let spral = match rec.spral.as_ref() {
            Some(s) => s,
            None => continue,
        };

        let bl = match baseline
            .baselines
            .iter()
            .find(|b| b.matrix_name == rec.matrix_name)
        {
            Some(b) => b,
            None => continue,
        };

        let rivrs_fac_s = bl.factor_ms / 1000.0;
        let ratio = if spral.factor_s > 0.0 {
            rivrs_fac_s / spral.factor_s
        } else {
            f64::NAN
        };

        eprintln!(
            "{:<35} {:>8} {:>10.3} {:>10.3} {:>7.2} {:>10.1e} {:>10.1e}",
            rec.matrix_name,
            rec.n,
            spral.factor_s,
            rivrs_fac_s,
            ratio,
            spral.backward_error,
            bl.backward_error
        );
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

        // Run SPRAL
        let spral_result =
            match run_spral_benchmark(matrix, &cli.spral_binary, cli.common.threads, nnz_lower) {
                Ok(parsed) => {
                    let status = if parsed.factor_flag < 0 {
                        format!("FACTOR_FAIL({})", parsed.factor_flag)
                    } else if parsed.analyse_flag < 0 {
                        format!("ANALYSE_FAIL({})", parsed.analyse_flag)
                    } else {
                        "OK".to_string()
                    };

                    Some(SpralResult {
                        matrix_name: meta.name.clone(),
                        category: meta.category.clone(),
                        n,
                        nnz: meta.nnz,
                        analyse_s: parsed.analyse_s,
                        factor_s: parsed.factor_s,
                        solve_s: parsed.solve_s,
                        total_s: parsed.analyse_s + parsed.factor_s + parsed.solve_s,
                        backward_error: parsed.backward_error,
                        num_delay: parsed.num_delay,
                        num_two: parsed.num_two,
                        num_neg: parsed.num_neg,
                        maxfront: parsed.maxfront,
                        maxsupernode: parsed.maxsupernode,
                        num_sup: parsed.num_sup,
                        num_factor: parsed.num_factor,
                        num_flops: parsed.num_flops,
                        not_first_pass: parsed.not_first_pass,
                        not_second_pass: parsed.not_second_pass,
                        matrix_rank: parsed.matrix_rank,
                        status,
                    })
                }
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
        if let Some(ref sr) = spral_result {
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
                    sr.num_delay,
                    sr.num_two,
                    sr.maxfront
                );
            }
            pass += 1;
        }

        results.push(ComparisonRecord {
            matrix_name: meta.name.clone(),
            category: meta.category.clone(),
            n,
            nnz: meta.nnz,
            spral: spral_result,
            rivrs: rivrs_result,
        });
    }

    // Print summary tables
    print_category_tables(&results, cli.common.rivrs, cli.common.threads);

    // Print comparison with baseline if requested
    if let Some(ref compare_path) = cli.common.compare {
        match std::fs::read_to_string(compare_path) {
            Ok(json) => match serde_json::from_str::<common::BaselineSuite>(&json) {
                Ok(baseline) => print_baseline_comparison(&results, &baseline),
                Err(e) => eprintln!("\nFailed to parse baseline JSON: {e}"),
            },
            Err(e) => eprintln!("\nFailed to read baseline file: {e}"),
        }
    }

    eprintln!("\n{pass}/{} completed successfully", pass + fail);

    // Build output suite
    let mode = if cli.common.rivrs {
        "comparison"
    } else {
        "spral-only"
    };

    let suite = SpralBenchmarkSuite {
        timestamp: common::timestamp(),
        platform: format!("{} {}", std::env::consts::OS, std::env::consts::ARCH),
        omp_num_threads: cli.common.threads,
        spral_version: "SPRAL-git".to_string(),
        rivrs_version: common::git_hash(),
        mode: mode.to_string(),
        results,
    };

    // Write JSON
    let out_dir = PathBuf::from("target/benchmarks/spral");
    std::fs::create_dir_all(&out_dir).expect("create output dir");
    let timestamp_slug = suite.timestamp.replace(':', "-").replace(' ', "_");
    let out_path = out_dir.join(format!("spral-benchmark-{timestamp_slug}.json"));
    let json = serde_json::to_string_pretty(&suite).expect("serialize");
    std::fs::write(&out_path, &json).expect("write JSON");
    eprintln!("JSON written to: {}", out_path.display());
}
