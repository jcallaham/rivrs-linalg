//! Shared infrastructure for external solver comparison binaries.
//!
//! Provides matrix format conversion, rivrs solver invocation, subprocess
//! management, and table printing helpers used by all comparison drivers.
//!
//! Included via `#[path = "common.rs"] mod common;` in each binary — not all
//! functions are used by every binary, so suppress dead code warnings.
#![allow(dead_code)]

use std::collections::BTreeMap;
use std::io::Write;
use std::path::PathBuf;
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

/// Threshold (in lower-triangle nnz) above which we use file mode instead of stdin.
pub const FILE_MODE_NNZ_THRESHOLD: usize = 1_000_000;

// ---------------------------------------------------------------------------
// Matrix format conversion
// ---------------------------------------------------------------------------

/// Extract lower triangle from a full symmetric CSC matrix.
/// Returns 1-indexed CSC arrays (col_ptr, row_idx, values).
pub fn extract_lower_triangle(matrix: &SparseColMat<usize, f64>) -> (Vec<i64>, Vec<i64>, Vec<f64>) {
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
/// Format: `n nnz\n` then `col_ptr[i]\n` lines, then `row[k] val[k]\n` lines.
pub fn format_spral_input(matrix: &SparseColMat<usize, f64>) -> String {
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

/// Format lower-triangle as COO text for MUMPS drivers.
/// Format: `n nnz\n` then `row col val\n` lines (1-indexed).
pub fn format_lower_coo_text(matrix: &SparseColMat<usize, f64>) -> String {
    let n = matrix.nrows();
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();

    let mut entries: Vec<(i64, i64, f64)> = Vec::new();

    for j in 0..n {
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for idx in start..end {
            let i = row_indices[idx];
            if i >= j {
                entries.push((i as i64 + 1, j as i64 + 1, values[idx])); // 1-indexed
            }
        }
    }

    let nnz = entries.len();
    let mut out = String::with_capacity(nnz * 40 + 32);
    out.push_str(&format!("{} {}\n", n, nnz));

    for &(row, col, val) in &entries {
        out.push_str(&format!("{} {} {:.17e}\n", row, col, val));
    }

    out
}

/// Write formatted matrix data to a temporary file, returning the path.
pub fn write_temp_matrix_file(data: &str, dir_name: &str) -> std::io::Result<PathBuf> {
    let dir = PathBuf::from(format!("/tmp/{dir_name}_input"));
    std::fs::create_dir_all(&dir)?;
    let path = dir.join("matrix.txt");
    std::fs::write(&path, data)?;
    Ok(path)
}

// ---------------------------------------------------------------------------
// Subprocess protocol
// ---------------------------------------------------------------------------

/// Solver-agnostic result parsed from subprocess output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExternalSolverResult {
    pub analyse_s: f64,
    pub factor_s: f64,
    pub solve_s: f64,
    pub backward_error: f64,
    /// Solver-specific stats (e.g., num_delay, num_neg, ordering, etc.)
    pub extra: BTreeMap<String, serde_json::Value>,
}

/// Run an external solver subprocess with sentinel-delimited output parsing.
///
/// - `binary`: path to the solver executable
/// - `input`: matrix data to feed (via stdin or temp file)
/// - `begin_sentinel` / `end_sentinel`: markers delimiting the output block
/// - `env_vars`: additional environment variables to set
/// - `args`: additional command-line arguments
/// - `use_file_mode`: if true, write input to a temp file and pass as argument
/// - `temp_dir_name`: directory name for temp file (used with file mode)
#[allow(clippy::too_many_arguments)]
pub fn run_solver_subprocess(
    binary: &str,
    input: &str,
    begin_sentinel: &str,
    end_sentinel: &str,
    env_vars: &[(&str, String)],
    args: &[String],
    use_file_mode: bool,
    temp_dir_name: &str,
) -> Result<ExternalSolverResult, String> {
    let mut cmd = Command::new(binary);
    cmd.stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());

    for &(key, ref val) in env_vars {
        cmd.env(key, val);
    }
    for arg in args {
        cmd.arg(arg);
    }

    if use_file_mode {
        let temp_path =
            write_temp_matrix_file(input, temp_dir_name).map_err(|e| format!("temp file: {e}"))?;
        cmd.arg(temp_path.to_str().unwrap_or("/tmp/solver_input/matrix.txt"));

        let output = cmd
            .spawn()
            .map_err(|e| format!("spawn: {e}"))?
            .wait_with_output()
            .map_err(|e| format!("wait: {e}"))?;

        let _ = std::fs::remove_file(&temp_path);
        parse_solver_output(&output, begin_sentinel, end_sentinel)
    } else {
        cmd.stdin(std::process::Stdio::piped());

        let mut child = cmd.spawn().map_err(|e| format!("spawn: {e}"))?;

        if let Some(stdin) = child.stdin.as_mut() {
            stdin
                .write_all(input.as_bytes())
                .map_err(|e| format!("stdin write: {e}"))?;
        }
        drop(child.stdin.take());

        let output = child.wait_with_output().map_err(|e| format!("wait: {e}"))?;
        parse_solver_output(&output, begin_sentinel, end_sentinel)
    }
}

/// Parse sentinel-delimited key-value output from a solver subprocess.
fn parse_solver_output(
    output: &std::process::Output,
    begin_sentinel: &str,
    end_sentinel: &str,
) -> Result<ExternalSolverResult, String> {
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("exit {}: {}", output.status, stderr));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);

    let begin = stdout
        .find(begin_sentinel)
        .ok_or_else(|| format!("missing {begin_sentinel}"))?;
    let end = stdout
        .find(end_sentinel)
        .ok_or_else(|| format!("missing {end_sentinel}"))?;

    let block = &stdout[begin + begin_sentinel.len()..end];

    let mut result = ExternalSolverResult {
        analyse_s: 0.0,
        factor_s: 0.0,
        solve_s: 0.0,
        backward_error: f64::NAN,
        extra: BTreeMap::new(),
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
            _ => {
                // Store everything else as extra stats
                if let Ok(i) = val.parse::<i64>() {
                    result.extra.insert(key.to_string(), serde_json::json!(i));
                } else if let Ok(f) = val.parse::<f64>() {
                    result.extra.insert(key.to_string(), serde_json::json!(f));
                } else {
                    result.extra.insert(key.to_string(), serde_json::json!(val));
                }
            }
        }
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// Rivrs solver invocation
// ---------------------------------------------------------------------------

/// Result from running the rivrs solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RivrsResult {
    pub matrix_name: String,
    pub analyse_s: f64,
    pub factor_s: f64,
    pub solve_s: f64,
    pub total_s: f64,
    pub backward_error: f64,
    pub num_delayed: usize,
    pub num_2x2: usize,
    pub max_front_size: usize,
    pub small_leaf_subtrees: usize,
    pub small_leaf_nodes: usize,
}

/// Run the rivrs SparseLDLT solver on a matrix and return timing/accuracy results.
pub fn run_rivrs_solver(
    matrix: &SparseColMat<usize, f64>,
    name: &str,
    threads: Option<usize>,
) -> Result<RivrsResult, String> {
    let n = matrix.nrows();

    let par = match threads {
        Some(t) if t > 1 => Par::rayon(t),
        _ => Par::Seq,
    };

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
    let factor_opts = FactorOptions {
        par,
        ..FactorOptions::default()
    };
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
        .solve(&b, stack, par)
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

/// Trim glibc memory on Linux when memory pressure is high.
/// Call after each matrix factorization to prevent OOM in containers.
pub fn maybe_trim_memory() {
    #[cfg(target_os = "linux")]
    {
        unsafe extern "C" {
            fn malloc_trim(pad: usize) -> i32;
        }
        let should_trim = std::fs::read_to_string("/sys/fs/cgroup/memory.current")
            .ok()
            .and_then(|c| c.trim().parse::<u64>().ok())
            .zip(
                std::fs::read_to_string("/sys/fs/cgroup/memory.max")
                    .ok()
                    .and_then(|m| m.trim().parse::<u64>().ok()),
            )
            .is_some_and(|(current, max)| current > max * 3 / 4);
        if should_trim {
            unsafe {
                malloc_trim(0);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// CLI argument parsing helpers
// ---------------------------------------------------------------------------

/// Common CLI arguments shared across comparison binaries.
pub struct CommonCliArgs {
    pub ci_only: bool,
    pub threads: Option<usize>,
    pub category: Option<String>,
    pub rivrs: bool,
    pub compare: Option<String>,
}

/// Parse common CLI arguments from argv, returning unrecognized args for
/// solver-specific parsing. The callback `extra_handler` is called for
/// unrecognized flags; return `Ok(skip)` with the number of additional
/// args consumed (0 for flags, 1 for --key value), or `Err(msg)` to abort.
pub fn parse_common_args<F>(mut extra_handler: F) -> Result<CommonCliArgs, String>
where
    F: FnMut(&[String], usize) -> Result<usize, String>,
{
    let args: Vec<String> = std::env::args().collect();
    let mut cli = CommonCliArgs {
        ci_only: false,
        threads: None,
        category: None,
        rivrs: false,
        compare: None,
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
            _ => {
                let skip = extra_handler(&args, i)?;
                i += skip;
            }
        }
        i += 1;
    }

    Ok(cli)
}

// ---------------------------------------------------------------------------
// Matrix loading
// ---------------------------------------------------------------------------

/// Load matrices from the registry, filtered by CLI args.
pub fn load_matrix_entries(cli: &CommonCliArgs) -> Vec<registry::MatrixMetadata> {
    let all = registry::load_registry().expect("Failed to load registry");
    all.into_iter()
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
        .collect()
}

// ---------------------------------------------------------------------------
// Table printing helpers
// ---------------------------------------------------------------------------

/// Print a comparison table header for solver vs rivrs.
pub fn print_solver_comparison_header(solver_name: &str, threads: Option<usize>) {
    let thread_str = threads
        .map(|t| format!("threads={t}"))
        .unwrap_or_else(|| "threads=default".to_string());
    eprintln!("\n=== Comparison: {solver_name} vs rivrs ({thread_str}) ===");
    eprintln!(
        "{:<35} {:>8} {:>10} {:>10} {:>7} {:>10} {:>10}",
        "Matrix",
        "n",
        &format!("{}_fac", solver_name.to_lowercase()),
        "rivrs_fac",
        "ratio",
        &format!("{}_be", solver_name.to_lowercase()),
        "rivrs_be"
    );
    eprintln!("{}", "-".repeat(100));
}

/// Print a comparison row for solver vs rivrs.
pub fn print_solver_comparison_row(
    name: &str,
    n: usize,
    solver_fac: f64,
    solver_be: f64,
    rivrs: Option<&RivrsResult>,
) {
    let (rivrs_fac, rivrs_be) = rivrs
        .map(|r| (r.factor_s, r.backward_error))
        .unwrap_or((f64::NAN, f64::NAN));
    let ratio = if solver_fac > 0.0 {
        rivrs_fac / solver_fac
    } else {
        f64::NAN
    };
    eprintln!(
        "{:<35} {:>8} {:>10.3} {:>10.3} {:>7.2} {:>10.1e} {:>10.1e}",
        name, n, solver_fac, rivrs_fac, ratio, solver_be, rivrs_be
    );
}

/// Print a solver-only table header.
pub fn print_solver_only_header() {
    eprintln!(
        "{:<35} {:>8} {:>10} {:>7} {:>7} {:>7} {:>9}",
        "Matrix", "n", "nnz", "ana_s", "fac_s", "slv_s", "bwd_err"
    );
    eprintln!("{}", "-".repeat(90));
}

/// Print a solver-only result row.
pub fn print_solver_only_row(
    name: &str,
    n: usize,
    nnz: usize,
    analyse_s: f64,
    factor_s: f64,
    solve_s: f64,
    backward_error: f64,
) {
    eprintln!(
        "{:<35} {:>8} {:>10} {:>7.3} {:>7.3} {:>7.3} {:>9.1e}",
        name, n, nnz, analyse_s, factor_s, solve_s, backward_error
    );
}

// ---------------------------------------------------------------------------
// JSON output
// ---------------------------------------------------------------------------

/// Generic comparison record for any solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComparisonRecord {
    pub matrix_name: String,
    pub category: String,
    pub n: usize,
    pub nnz: usize,
    pub solver: Option<ExternalSolverResult>,
    pub rivrs: Option<RivrsResult>,
}

/// Benchmark suite output envelope.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkSuite {
    pub timestamp: String,
    pub platform: String,
    pub solver_name: String,
    pub solver_version: String,
    pub rivrs_version: String,
    pub threads: Option<usize>,
    pub mode: String,
    pub results: Vec<ComparisonRecord>,
}

/// Get a simple timestamp string (Unix seconds).
pub fn timestamp() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    format!("{secs}")
}

/// Get the short git hash of HEAD.
pub fn git_hash() -> String {
    Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .unwrap_or_else(|| "unknown".to_string())
        .trim()
        .to_string()
}

/// Write a benchmark suite to JSON in the given output directory.
pub fn write_suite_json(suite: &BenchmarkSuite, output_dir: &str) -> PathBuf {
    let out_dir = PathBuf::from(output_dir);
    std::fs::create_dir_all(&out_dir).expect("create output dir");
    let timestamp_slug = suite.timestamp.replace(':', "-").replace(' ', "_");
    let solver_slug = suite.solver_name.to_lowercase();
    let out_path = out_dir.join(format!("{solver_slug}-benchmark-{timestamp_slug}.json"));
    let json = serde_json::to_string_pretty(suite).expect("serialize");
    std::fs::write(&out_path, &json).expect("write JSON");
    out_path
}

// ---------------------------------------------------------------------------
// Baseline comparison (from baseline_collection.rs format)
// ---------------------------------------------------------------------------

/// Rivrs baseline format (from baseline_collection.rs).
#[derive(Debug, Clone, Deserialize)]
pub struct BaselineSuite {
    #[allow(dead_code)]
    pub timestamp: String,
    #[allow(dead_code)]
    pub platform: String,
    #[allow(dead_code)]
    pub solver_version: String,
    pub baselines: Vec<BaselineEntry>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct BaselineEntry {
    pub matrix_name: String,
    #[allow(dead_code)]
    pub matrix_dim: usize,
    #[allow(dead_code)]
    pub matrix_nnz: usize,
    #[allow(dead_code)]
    pub ordering_ms: f64,
    #[allow(dead_code)]
    pub symbolic_ms: f64,
    pub factor_ms: f64,
    #[allow(dead_code)]
    pub solve_ms: f64,
    #[allow(dead_code)]
    pub total_ms: f64,
    pub backward_error: f64,
    #[allow(dead_code)]
    pub num_supernodes: usize,
    #[allow(dead_code)]
    pub max_front_size: usize,
    #[allow(dead_code)]
    pub factorization_stats: BaselineFactStats,
}

#[derive(Debug, Clone, Deserialize)]
pub struct BaselineFactStats {
    #[allow(dead_code)]
    pub total_1x1_pivots: usize,
    #[allow(dead_code)]
    pub total_2x2_pivots: usize,
    #[allow(dead_code)]
    pub total_delayed: usize,
    #[allow(dead_code)]
    pub zero_pivots: usize,
    #[allow(dead_code)]
    pub max_front_size: usize,
}

/// Print a comparison table between an external solver and a rivrs baseline file.
pub fn print_baseline_comparison(
    solver_name: &str,
    results: &[ComparisonRecord],
    baseline: &BaselineSuite,
) {
    eprintln!("\n=== Comparison: {solver_name} vs rivrs baseline ===");
    eprintln!(
        "{:<35} {:>8} {:>10} {:>10} {:>7} {:>10} {:>10}",
        "Matrix",
        "n",
        &format!("{}_fac", solver_name.to_lowercase()),
        "rivrs_fac",
        "ratio",
        &format!("{}_be", solver_name.to_lowercase()),
        "rivrs_be"
    );
    eprintln!("{}", "-".repeat(100));

    for rec in results {
        let solver = match rec.solver.as_ref() {
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
        let ratio = if solver.factor_s > 0.0 {
            rivrs_fac_s / solver.factor_s
        } else {
            f64::NAN
        };

        eprintln!(
            "{:<35} {:>8} {:>10.3} {:>10.3} {:>7.2} {:>10.1e} {:>10.1e}",
            rec.matrix_name,
            rec.n,
            solver.factor_s,
            rivrs_fac_s,
            ratio,
            solver.backward_error,
            bl.backward_error
        );
    }
}
