//! Systematic performance baseline collection for the SuiteSparse test suite.
//!
//! Records per-phase timing, per-supernode statistics, peak RSS, and backward
//! error for each SuiteSparse matrix. Outputs structured JSON for regression
//! tracking and a human-readable summary to stderr.
//!
//! Usage:
//!   cargo run --example baseline_collection --features diagnostic --release
//!   cargo run --example baseline_collection --features diagnostic --release -- --ci-only
//!   cargo run --example baseline_collection --features diagnostic --release -- --compare <prev.json>

use std::time::Instant;

use faer::Col;
use serde::{Deserialize, Serialize};

use rivrs_sparse::aptp::{AnalyzeOptions, FactorOptions, SparseLDLT};
use rivrs_sparse::benchmarking::read_peak_rss_kb;
use rivrs_sparse::io::registry;
use rivrs_sparse::validate::sparse_backward_error;

/// Per-supernode stats snapshot (serializable subset of PerSupernodeStats).
#[derive(Debug, Clone, Serialize, Deserialize)]
struct SupernodeRecord {
    snode_id: usize,
    front_size: usize,
    num_fully_summed: usize,
    num_eliminated: usize,
    num_delayed: usize,
    num_1x1: usize,
    num_2x2: usize,
    max_l_entry: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    assembly_us: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    kernel_us: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    extraction_us: Option<f64>,
}

/// Aggregate factorization stats (serializable).
#[derive(Debug, Clone, Serialize, Deserialize)]
struct FactorizationRecord {
    total_1x1_pivots: usize,
    total_2x2_pivots: usize,
    total_delayed: usize,
    zero_pivots: usize,
    max_front_size: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    total_assembly_us: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    total_kernel_us: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    total_extraction_us: Option<f64>,
}

/// Per-matrix performance baseline.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct PerformanceBaseline {
    matrix_name: String,
    matrix_dim: usize,
    matrix_nnz: usize,
    ordering_ms: f64,
    symbolic_ms: f64,
    factor_ms: f64,
    solve_ms: f64,
    total_ms: f64,
    peak_rss_kb: Option<u64>,
    backward_error: f64,
    num_supernodes: usize,
    max_front_size: usize,
    factorization_stats: FactorizationRecord,
    per_supernode_stats: Vec<SupernodeRecord>,
}

/// Collection of baselines with metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct BaselineSuite {
    timestamp: String,
    platform: String,
    solver_version: String,
    ordering_strategy: String,
    baselines: Vec<PerformanceBaseline>,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let ci_only = args.iter().any(|a| a == "--ci-only");
    let compare_path = args
        .iter()
        .position(|a| a == "--compare")
        .and_then(|i| args.get(i + 1));

    let all = registry::load_registry().expect("Failed to load registry");
    let matrices: Vec<_> = all
        .iter()
        .filter(|m| m.source == "suitesparse" && (!ci_only || m.ci_subset))
        .collect();

    eprintln!(
        "Collecting baselines for {} SuiteSparse matrices{}",
        matrices.len(),
        if ci_only { " (CI subset)" } else { "" }
    );

    eprintln!(
        "{:<30} {:>8} {:>10} {:>8} {:>8} {:>8} {:>8} {:>12} {:>10}",
        "Matrix", "n", "nnz", "ord_ms", "sym_ms", "fac_ms", "slv_ms", "backward_err", "rss_kb"
    );
    eprintln!("{}", "-".repeat(120));

    let mut baselines = Vec::new();

    for meta in &matrices {
        let test = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(t)) => t,
            _ => {
                eprintln!("{:<30} SKIP (not found)", meta.name);
                continue;
            }
        };
        let a = &test.matrix;
        let n = a.nrows();

        // Create RHS: b = A * x_true
        let x_true: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
        let sym = a.symbolic();
        let cp = sym.col_ptr();
        let ri = sym.row_idx();
        let vals = a.val();
        let mut b_vec = vec![0.0f64; n];
        for j in 0..n {
            for idx in cp[j]..cp[j + 1] {
                b_vec[ri[idx]] += vals[idx] * x_true[j];
            }
        }
        let b = Col::from_fn(n, |i| b_vec[i]);

        // Phase 1: Ordering (analyze includes ordering)
        let t_total = Instant::now();
        let t_ord = Instant::now();
        let opts = AnalyzeOptions::default();
        let mut solver = match SparseLDLT::analyze_with_matrix(a, &opts) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("{:<30} ANALYZE FAILED: {}", meta.name, e);
                continue;
            }
        };
        let ordering_ms = t_ord.elapsed().as_secs_f64() * 1000.0;

        // Phase 2: Symbolic (already done in analyze_with_matrix, so record ~0)
        let symbolic_ms = 0.0;

        // Phase 3: Factor
        let t_fac = Instant::now();
        let factor_opts = FactorOptions::default();
        match solver.factor(a, &factor_opts) {
            Ok(()) => {}
            Err(e) => {
                eprintln!("{:<30} FACTOR FAILED: {}", meta.name, e);
                continue;
            }
        }
        let factor_ms = t_fac.elapsed().as_secs_f64() * 1000.0;

        // Phase 4: Solve
        let t_slv = Instant::now();
        let scratch = solver.solve_scratch(1);
        let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
        let stack = faer::dyn_stack::MemStack::new(&mut mem);
        let x = match solver.solve(&b, stack) {
            Ok(x) => x,
            Err(e) => {
                eprintln!("{:<30} SOLVE FAILED: {}", meta.name, e);
                continue;
            }
        };
        let solve_ms = t_slv.elapsed().as_secs_f64() * 1000.0;
        let total_ms = t_total.elapsed().as_secs_f64() * 1000.0;

        let be = sparse_backward_error(a, &x, &b);
        let peak_rss = read_peak_rss_kb();

        let stats = solver.stats().unwrap();
        let per_sn = solver.per_supernode_stats().unwrap();

        let factorization_record = FactorizationRecord {
            total_1x1_pivots: stats.total_1x1_pivots,
            total_2x2_pivots: stats.total_2x2_pivots,
            total_delayed: stats.total_delayed,
            zero_pivots: stats.zero_pivots,
            max_front_size: stats.max_front_size,
            #[cfg(feature = "diagnostic")]
            total_assembly_us: Some(stats.total_assembly_time.as_secs_f64() * 1_000_000.0),
            #[cfg(feature = "diagnostic")]
            total_kernel_us: Some(stats.total_kernel_time.as_secs_f64() * 1_000_000.0),
            #[cfg(feature = "diagnostic")]
            total_extraction_us: Some(stats.total_extraction_time.as_secs_f64() * 1_000_000.0),
            #[cfg(not(feature = "diagnostic"))]
            total_assembly_us: None,
            #[cfg(not(feature = "diagnostic"))]
            total_kernel_us: None,
            #[cfg(not(feature = "diagnostic"))]
            total_extraction_us: None,
        };

        let sn_records: Vec<SupernodeRecord> = per_sn
            .iter()
            .map(|s| SupernodeRecord {
                snode_id: s.snode_id,
                front_size: s.front_size,
                num_fully_summed: s.num_fully_summed,
                num_eliminated: s.num_eliminated,
                num_delayed: s.num_delayed,
                num_1x1: s.num_1x1,
                num_2x2: s.num_2x2,
                max_l_entry: s.max_l_entry,
                #[cfg(feature = "diagnostic")]
                assembly_us: Some(s.assembly_time.as_secs_f64() * 1_000_000.0),
                #[cfg(feature = "diagnostic")]
                kernel_us: Some(s.kernel_time.as_secs_f64() * 1_000_000.0),
                #[cfg(feature = "diagnostic")]
                extraction_us: Some(s.extraction_time.as_secs_f64() * 1_000_000.0),
                #[cfg(not(feature = "diagnostic"))]
                assembly_us: None,
                #[cfg(not(feature = "diagnostic"))]
                kernel_us: None,
                #[cfg(not(feature = "diagnostic"))]
                extraction_us: None,
            })
            .collect();

        eprintln!(
            "{:<30} {:>8} {:>10} {:>8.1} {:>8.1} {:>8.1} {:>8.1} {:>12.2e} {:>10}",
            meta.name,
            n,
            meta.nnz,
            ordering_ms,
            symbolic_ms,
            factor_ms,
            solve_ms,
            be,
            peak_rss.map_or("N/A".to_string(), |v| v.to_string())
        );

        baselines.push(PerformanceBaseline {
            matrix_name: meta.name.clone(),
            matrix_dim: n,
            matrix_nnz: meta.nnz,
            ordering_ms,
            symbolic_ms,
            factor_ms,
            solve_ms,
            total_ms,
            peak_rss_kb: peak_rss,
            backward_error: be,
            num_supernodes: per_sn.len(),
            max_front_size: stats.max_front_size,
            factorization_stats: factorization_record,
            per_supernode_stats: sn_records,
        });
    }

    // Build suite
    let git_hash = std::process::Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .unwrap_or_else(|| "unknown".to_string())
        .trim()
        .to_string();

    let suite = BaselineSuite {
        timestamp: chrono_timestamp(),
        platform: format!("{} {}", std::env::consts::OS, std::env::consts::ARCH),
        solver_version: git_hash,
        ordering_strategy: "MatchOrderMetis".to_string(),
        baselines,
    };

    // Write JSON
    let baseline_dir = std::path::PathBuf::from("target/benchmarks/baselines");
    std::fs::create_dir_all(&baseline_dir).expect("create baseline dir");
    let timestamp_slug = suite.timestamp.replace(':', "-").replace(' ', "_");
    let out_path = baseline_dir.join(format!("baseline-{timestamp_slug}.json"));
    let json = serde_json::to_string_pretty(&suite).expect("serialize");
    std::fs::write(&out_path, &json).expect("write baseline");
    eprintln!("\nBaseline written to: {}", out_path.display());
    eprintln!(
        "{} matrices, {} baselines",
        matrices.len(),
        suite.baselines.len()
    );

    // Also print JSON to stdout for piping
    println!("{json}");

    // Compare with previous baseline if requested
    if let Some(prev_path) = compare_path {
        let prev_json = std::fs::read_to_string(prev_path).expect("read previous baseline");
        let prev: BaselineSuite =
            serde_json::from_str(&prev_json).expect("parse previous baseline");
        compare_baselines(&prev, &suite);
    }
}

fn compare_baselines(prev: &BaselineSuite, curr: &BaselineSuite) {
    eprintln!("\n--- Comparison with previous baseline ---");
    eprintln!("Previous: {} ({})", prev.solver_version, prev.timestamp);
    eprintln!("Current:  {} ({})", curr.solver_version, curr.timestamp);
    eprintln!(
        "{:<30} {:>10} {:>10} {:>10} {:>10}",
        "Matrix", "prev_ms", "curr_ms", "change%", "status"
    );
    eprintln!("{}", "-".repeat(80));

    for cb in &curr.baselines {
        if let Some(pb) = prev
            .baselines
            .iter()
            .find(|b| b.matrix_name == cb.matrix_name)
        {
            let change = if pb.factor_ms > 0.001 {
                (cb.factor_ms - pb.factor_ms) / pb.factor_ms * 100.0
            } else {
                0.0
            };
            let status = if change > 10.0 {
                "REGRESS"
            } else if change < -10.0 {
                "IMPROVE"
            } else {
                "OK"
            };
            eprintln!(
                "{:<30} {:>10.1} {:>10.1} {:>+9.1}% {:>10}",
                cb.matrix_name, pb.factor_ms, cb.factor_ms, change, status
            );
        }
    }
}

fn chrono_timestamp() -> String {
    // Simple ISO-8601 timestamp without chrono dependency
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    format!("{secs}")
}
