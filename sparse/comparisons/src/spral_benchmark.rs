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
//! cargo run --example spral_benchmark --release -- --ci-only
//!
//! # Side-by-side comparison with rivrs
//! cargo run --example spral_benchmark --release -- --ci-only --rivrs
//!
//! # Set OMP thread count for SPRAL
//! cargo run --example spral_benchmark --release -- --ci-only --threads 4
//!
//! # Compare against a previously collected rivrs baseline
//! cargo run --example spral_benchmark --release -- --ci-only \
//!   --compare target/benchmarks/baselines/baseline-latest.json
//!
//! # Filter by category
//! cargo run --example spral_benchmark --release -- --category hard-indefinite
//! ```

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use faer::Col;
use faer::Par;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::SparseColMat;
use serde::{Deserialize, Serialize};

use rivrs_sparse::io::registry;
use rivrs_sparse::symmetric::{AnalyzeOptions, FactorOptions, SparseLDLT};
use rivrs_sparse::validate::sparse_backward_error;

const DEFAULT_SPRAL_BIN: &str = "/tmp/spral_benchmark";

/// Threshold (in lower-triangle nnz) above which we use file mode instead of stdin.
const FILE_MODE_NNZ_THRESHOLD: usize = 1_000_000;

// ---------------------------------------------------------------------------
// Data types
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
struct RivrsResult {
    matrix_name: String,
    analyse_s: f64,
    factor_s: f64,
    solve_s: f64,
    total_s: f64,
    backward_error: f64,
    num_delayed: usize,
    num_2x2: usize,
    max_front_size: usize,
    small_leaf_subtrees: usize,
    small_leaf_nodes: usize,
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

/// Rivrs baseline format (from baseline_collection.rs).
#[derive(Debug, Clone, Deserialize)]
struct BaselineSuite {
    #[allow(dead_code)]
    timestamp: String,
    #[allow(dead_code)]
    platform: String,
    #[allow(dead_code)]
    solver_version: String,
    baselines: Vec<BaselineEntry>,
}

#[derive(Debug, Clone, Deserialize)]
struct BaselineEntry {
    matrix_name: String,
    #[allow(dead_code)]
    matrix_dim: usize,
    #[allow(dead_code)]
    matrix_nnz: usize,
    #[allow(dead_code)]
    ordering_ms: f64,
    #[allow(dead_code)]
    symbolic_ms: f64,
    factor_ms: f64,
    #[allow(dead_code)]
    solve_ms: f64,
    #[allow(dead_code)]
    total_ms: f64,
    backward_error: f64,
    #[allow(dead_code)]
    num_supernodes: usize,
    #[allow(dead_code)]
    max_front_size: usize,
    #[allow(dead_code)]
    factorization_stats: BaselineFactStats,
}

#[derive(Debug, Clone, Deserialize)]
struct BaselineFactStats {
    #[allow(dead_code)]
    total_1x1_pivots: usize,
    #[allow(dead_code)]
    total_2x2_pivots: usize,
    #[allow(dead_code)]
    total_delayed: usize,
    #[allow(dead_code)]
    zero_pivots: usize,
    #[allow(dead_code)]
    max_front_size: usize,
}

// ---------------------------------------------------------------------------
// CLI argument parsing
// ---------------------------------------------------------------------------

struct CliArgs {
    ci_only: bool,
    threads: Option<usize>,
    category: Option<String>,
    rivrs: bool,
    compare: Option<String>,
    spral_binary: String,
}

fn parse_args() -> CliArgs {
    let args: Vec<String> = std::env::args().collect();
    let mut cli = CliArgs {
        ci_only: false,
        threads: None,
        category: None,
        rivrs: false,
        compare: None,
        spral_binary: DEFAULT_SPRAL_BIN.to_string(),
    };

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--ci-only" => cli.ci_only = true,
            "--rivrs" => cli.rivrs = true,
            "--threads" => {
                i += 1;
                cli.threads = args.get(i).and_then(|s| s.parse().ok());
            }
            "--category" => {
                i += 1;
                cli.category = args.get(i).cloned();
            }
            "--compare" => {
                i += 1;
                cli.compare = args.get(i).cloned();
            }
            "--spral-binary" => {
                i += 1;
                if let Some(path) = args.get(i) {
                    cli.spral_binary = path.clone();
                }
            }
            other => {
                eprintln!("Unknown argument: {other}");
                eprintln!(
                    "Usage: spral_benchmark [--ci-only] [--threads N] [--category CAT] \
                     [--rivrs] [--compare FILE] [--spral-binary PATH]"
                );
                std::process::exit(1);
            }
        }
        i += 1;
    }

    cli
}

// ---------------------------------------------------------------------------
// Matrix formatting for SPRAL
// ---------------------------------------------------------------------------

/// Extract lower triangle from a full symmetric CSC matrix.
/// SPRAL expects lower-triangle-only CSC input.
fn extract_lower_triangle(matrix: &SparseColMat<usize, f64>) -> (Vec<i64>, Vec<i64>, Vec<f64>) {
    let n = matrix.nrows();
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();

    let mut lower_ptr = vec![0i64; n + 1];
    let mut lower_row = Vec::new();
    let mut lower_val = Vec::new();

    for j in 0..n {
        lower_ptr[j] = lower_row.len() as i64 + 1; // 1-indexed
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for idx in start..end {
            let i = row_indices[idx];
            if i >= j {
                lower_row.push(i as i64 + 1); // 1-indexed
                lower_val.push(values[idx]);
            }
        }
    }
    lower_ptr[n] = lower_row.len() as i64 + 1;

    (lower_ptr, lower_row, lower_val)
}

/// Format lower-triangle CSC as text for the SPRAL benchmark driver.
fn format_spral_input(matrix: &SparseColMat<usize, f64>) -> String {
    let n = matrix.nrows();
    let (lower_ptr, lower_row, lower_val) = extract_lower_triangle(matrix);
    let nnz = lower_row.len();

    let mut out = String::new();
    out.push_str(&format!("{} {}\n", n, nnz));

    for &ptr in &lower_ptr[..=n] {
        out.push_str(&format!("{}\n", ptr));
    }

    for k in 0..nnz {
        out.push_str(&format!("{} {:.17e}\n", lower_row[k], lower_val[k]));
    }

    out
}

/// Write formatted matrix to a temporary file, returning the path.
fn write_temp_matrix_file(matrix: &SparseColMat<usize, f64>) -> std::io::Result<PathBuf> {
    let input = format_spral_input(matrix);
    let dir = PathBuf::from("/tmp/spral_benchmark_input");
    std::fs::create_dir_all(&dir)?;
    let path = dir.join("matrix.txt");
    std::fs::write(&path, input)?;
    Ok(path)
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
    let use_file_mode = nnz_lower > FILE_MODE_NNZ_THRESHOLD;

    let mut cmd = Command::new(spral_bin);
    cmd.stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .env("OMP_CANCELLATION", "TRUE")
        .env("OMP_PROC_BIND", "TRUE");

    if let Some(t) = threads {
        cmd.env("OMP_NUM_THREADS", t.to_string());
    }

    if use_file_mode {
        // Write matrix to temp file and pass path as argument
        let temp_path = write_temp_matrix_file(matrix).map_err(|e| format!("temp file: {e}"))?;
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

        // Clean up temp file
        let _ = std::fs::remove_file(&temp_path);

        parse_spral_result(&output)
    } else {
        // Use stdin mode
        cmd.stdin(std::process::Stdio::piped());

        let mut child = cmd.spawn().map_err(|e| format!("spawn: {e}"))?;

        let input = format_spral_input(matrix);
        if let Some(stdin) = child.stdin.as_mut() {
            stdin
                .write_all(input.as_bytes())
                .map_err(|e| format!("stdin write: {e}"))?;
        }
        drop(child.stdin.take()); // Close stdin to signal EOF

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

    // Find content between sentinels
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
        // Split on first whitespace
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
// Rivrs solver
// ---------------------------------------------------------------------------

fn run_rivrs_solver(matrix: &SparseColMat<usize, f64>, name: &str) -> Result<RivrsResult, String> {
    let n = matrix.nrows();

    // Generate RHS: b = A * ones(n) (using full symmetric CSC)
    let symbolic = matrix.symbolic();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();
    let values = matrix.val();

    let mut b_vec = vec![0.0f64; n];
    for j in 0..n {
        for idx in col_ptrs[j]..col_ptrs[j + 1] {
            b_vec[row_indices[idx]] += values[idx];
        }
    }
    let b = Col::from_fn(n, |i| b_vec[i]);

    let t_total = Instant::now();

    // Analyse
    let t_analyse = Instant::now();
    let opts = AnalyzeOptions::default();
    let mut solver =
        SparseLDLT::analyze_with_matrix(matrix, &opts).map_err(|e| format!("analyse: {e}"))?;
    let analyse_s = t_analyse.elapsed().as_secs_f64();

    // Factor
    let t_factor = Instant::now();
    let factor_opts = FactorOptions::default();
    solver
        .factor(matrix, &factor_opts)
        .map_err(|e| format!("factor: {e}"))?;
    let factor_s = t_factor.elapsed().as_secs_f64();

    // Solve
    let t_solve = Instant::now();
    let scratch = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);
    let stack = MemStack::new(&mut mem);
    let x = solver
        .solve(&b, stack, Par::Seq)
        .map_err(|e| format!("solve: {e}"))?;
    let solve_s = t_solve.elapsed().as_secs_f64();

    let total_s = t_total.elapsed().as_secs_f64();

    let be = sparse_backward_error(matrix, &x, &b);

    let stats = solver.stats().ok_or("no factorization stats")?;

    Ok(RivrsResult {
        matrix_name: name.to_string(),
        analyse_s,
        factor_s,
        solve_s,
        total_s,
        backward_error: be,
        num_delayed: stats.total_delayed,
        num_2x2: stats.total_2x2_pivots,
        max_front_size: stats.max_front_size,
        small_leaf_subtrees: stats.small_leaf_subtrees,
        small_leaf_nodes: stats.small_leaf_nodes,
    })
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

fn print_baseline_comparison(results: &[ComparisonRecord], baseline: &BaselineSuite) {
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

    // Load matrix registry
    let all = registry::load_registry().expect("Failed to load registry");
    let matrices: Vec<_> = all
        .iter()
        .filter(|m| {
            if m.source != "suitesparse" {
                return false;
            }
            if cli.ci_only && !m.ci_subset {
                return false;
            }
            if let Some(ref cat) = cli.category {
                if m.category != *cat {
                    return false;
                }
            }
            true
        })
        .collect();

    eprintln!(
        "SPRAL benchmark: {} matrices{}{}",
        matrices.len(),
        if cli.ci_only { " (CI subset)" } else { "" },
        cli.threads
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
            match run_spral_benchmark(matrix, &cli.spral_binary, cli.threads, nnz_lower) {
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
        let rivrs_result = if cli.rivrs {
            match run_rivrs_solver(matrix, &meta.name) {
                Ok(r) => Some(r),
                Err(e) => {
                    eprintln!("{:<35} RIVRS FAILED: {}", meta.name, e);
                    None
                }
            }
        } else {
            None
        };

        // Print progress row
        if let Some(ref sr) = spral_result {
            if cli.rivrs {
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
    print_category_tables(&results, cli.rivrs, cli.threads);

    // Print comparison with baseline if requested
    if let Some(ref compare_path) = cli.compare {
        match std::fs::read_to_string(compare_path) {
            Ok(json) => match serde_json::from_str::<BaselineSuite>(&json) {
                Ok(baseline) => print_baseline_comparison(&results, &baseline),
                Err(e) => eprintln!("\nFailed to parse baseline JSON: {e}"),
            },
            Err(e) => eprintln!("\nFailed to read baseline file: {e}"),
        }
    }

    eprintln!("\n{pass}/{} completed successfully", pass + fail);

    // Build output suite
    let git_hash = Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .unwrap_or_else(|| "unknown".to_string())
        .trim()
        .to_string();

    let mode = if cli.rivrs {
        "comparison"
    } else {
        "spral-only"
    };

    let suite = SpralBenchmarkSuite {
        timestamp: chrono_timestamp(),
        platform: format!("{} {}", std::env::consts::OS, std::env::consts::ARCH),
        omp_num_threads: cli.threads,
        spral_version: "SPRAL-git".to_string(),
        rivrs_version: git_hash,
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

    // Also print JSON to stdout for piping
    //     println!("{json}");
}

fn chrono_timestamp() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    format!("{secs}")
}
