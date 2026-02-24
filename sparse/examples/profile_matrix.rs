//! Profile the full factorization pipeline for a single matrix.
//!
//! Runs analyze → factor → solve with `diagnostic` instrumentation, then
//! prints a per-supernode timing breakdown and optionally exports Chrome
//! Trace JSON (viewable in Perfetto or chrome://tracing).
//!
//! # Usage
//!
//! ```bash
//! # By registry name (from metadata.json):
//! cargo run --example profile_matrix --features diagnostic --release -- d_pretok
//!
//! # By file path:
//! cargo run --example profile_matrix --features diagnostic --release -- --path test-data/suitesparse/hard-indefinite/stokes128/stokes128.mtx
//!
//! # Export Chrome Trace JSON:
//! cargo run --example profile_matrix --features diagnostic --release -- d_pretok --trace /tmp/trace.json
//!
//! # List available matrices:
//! cargo run --example profile_matrix --features diagnostic --release -- --list
//! ```

use std::path::Path;
use std::time::Instant;

use faer::Col;
use faer::Par;

use rivrs_sparse::aptp::{AnalyzeOptions, FactorOptions, PerSupernodeStats, SparseLDLT};
use rivrs_sparse::io::mtx::load_mtx;
use rivrs_sparse::io::registry;
use rivrs_sparse::validate::sparse_backward_error;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.iter().any(|a| a == "--list") {
        list_matrices();
        return;
    }

    if args.iter().any(|a| a == "--help" || a == "-h") {
        print_usage();
        return;
    }

    // Parse arguments
    let trace_path = args
        .iter()
        .position(|a| a == "--trace")
        .and_then(|i| args.get(i + 1))
        .cloned();

    let matrix_path = args
        .iter()
        .position(|a| a == "--path")
        .and_then(|i| args.get(i + 1))
        .cloned();

    // Matrix name is the first positional arg (not a flag or flag value)
    let matrix_name = args
        .iter()
        .skip(1)
        .find(|a| {
            !a.starts_with("--") && {
                // Skip values that follow --trace or --path
                let prev_idx = args.iter().position(|x| std::ptr::eq(x, *a)).unwrap();
                prev_idx == 0 || (args[prev_idx - 1] != "--trace" && args[prev_idx - 1] != "--path")
            }
        })
        .cloned();

    // Load matrix
    let (name, matrix) = if let Some(ref path) = matrix_path {
        let p = Path::new(path);
        let name = p
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();
        let m = load_mtx(p).unwrap_or_else(|e| {
            eprintln!("Failed to load {}: {}", path, e);
            std::process::exit(1);
        });
        (name, m)
    } else if let Some(ref name) = matrix_name {
        match registry::load_test_matrix(name) {
            Ok(Some(t)) => (name.clone(), t.matrix),
            Ok(None) => {
                eprintln!("Matrix '{}' found in registry but not on disk.", name);
                eprintln!("Run with --list to see available matrices.");
                std::process::exit(1);
            }
            Err(e) => {
                eprintln!("Failed to load '{}': {}", name, e);
                eprintln!("Run with --list to see available matrices.");
                std::process::exit(1);
            }
        }
    } else {
        eprintln!("No matrix specified. Provide a name or --path.");
        eprintln!();
        print_usage();
        std::process::exit(1);
    };

    let n = matrix.nrows();
    let nnz = matrix.val().len();
    eprintln!("Matrix: {} (n={}, nnz={})", name, n, nnz);
    eprintln!();

    // --- Analyze ---
    let t0 = Instant::now();
    let opts = AnalyzeOptions::default();
    let mut solver = SparseLDLT::analyze_with_matrix(&matrix, &opts).unwrap_or_else(|e| {
        eprintln!("Analyze failed: {}", e);
        std::process::exit(1);
    });
    let analyze_ms = t0.elapsed().as_secs_f64() * 1000.0;

    // --- Factor ---
    let t1 = Instant::now();
    let factor_opts = FactorOptions::default();
    solver.factor(&matrix, &factor_opts).unwrap_or_else(|e| {
        eprintln!("Factor failed: {}", e);
        std::process::exit(1);
    });
    let factor_ms = t1.elapsed().as_secs_f64() * 1000.0;

    // --- Solve ---
    let x_true: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
    let sym = matrix.symbolic();
    let cp = sym.col_ptr();
    let ri = sym.row_idx();
    let vals = matrix.val();
    let mut b_vec = vec![0.0f64; n];
    for j in 0..n {
        for idx in cp[j]..cp[j + 1] {
            b_vec[ri[idx]] += vals[idx] * x_true[j];
        }
    }
    let b = Col::from_fn(n, |i| b_vec[i]);

    let t2 = Instant::now();
    let scratch = solver.solve_scratch(1);
    let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
    let stack = faer::dyn_stack::MemStack::new(&mut mem);
    let x = solver.solve(&b, stack, Par::Seq).unwrap_or_else(|e| {
        eprintln!("Solve failed: {}", e);
        std::process::exit(1);
    });
    let solve_ms = t2.elapsed().as_secs_f64() * 1000.0;

    let be = sparse_backward_error(&matrix, &x, &b);
    let total_ms = t0.elapsed().as_secs_f64() * 1000.0;

    // --- Phase summary ---
    let stats = solver.stats().unwrap();
    eprintln!("=== Phase Timing ===");
    eprintln!("  Analyze:  {:>10.2} ms", analyze_ms);
    eprintln!("  Factor:   {:>10.2} ms", factor_ms);
    eprintln!("  Solve:    {:>10.2} ms", solve_ms);
    eprintln!("  Total:    {:>10.2} ms", total_ms);
    eprintln!();

    // --- Factorization summary ---
    eprintln!("=== Factorization Summary ===");
    eprintln!(
        "  Supernodes:     {} (before amalg: {}, merges: {})",
        stats.supernodes_after_amalgamation,
        stats.supernodes_before_amalgamation,
        stats.merges_performed,
    );
    eprintln!("  Max front size: {}", stats.max_front_size);
    eprintln!("  1x1 pivots:    {}", stats.total_1x1_pivots);
    eprintln!("  2x2 pivots:    {}", stats.total_2x2_pivots);
    eprintln!("  Delayed:        {}", stats.total_delayed);
    eprintln!("  Zero pivots:    {}", stats.zero_pivots);
    eprintln!("  Backward error: {:.2e}", be);

    #[cfg(feature = "diagnostic")]
    {
        let asm_ms = stats.total_assembly_time.as_secs_f64() * 1000.0;
        let kern_ms = stats.total_kernel_time.as_secs_f64() * 1000.0;
        let ext_ms = stats.total_extraction_time.as_secs_f64() * 1000.0;
        let accounted = asm_ms + kern_ms + ext_ms;
        eprintln!();
        eprintln!("=== Factor Time Breakdown ===");
        eprintln!(
            "  Assembly:    {:>10.2} ms  ({:>5.1}%)",
            asm_ms,
            asm_ms / factor_ms * 100.0
        );
        eprintln!(
            "  Kernel:      {:>10.2} ms  ({:>5.1}%)",
            kern_ms,
            kern_ms / factor_ms * 100.0
        );
        eprintln!(
            "  Extraction:  {:>10.2} ms  ({:>5.1}%)",
            ext_ms,
            ext_ms / factor_ms * 100.0
        );
        eprintln!(
            "  Unaccounted: {:>10.2} ms  ({:>5.1}%)",
            factor_ms - accounted,
            (factor_ms - accounted) / factor_ms * 100.0,
        );

        // Sub-phase breakdown
        let zero_ms = stats.total_zero_time.as_secs_f64() * 1000.0;
        let g2l_ms = stats.total_g2l_time.as_secs_f64() * 1000.0;
        let scatter_ms = stats.total_scatter_time.as_secs_f64() * 1000.0;
        let ea_ms = stats.total_extend_add_time.as_secs_f64() * 1000.0;
        let exf_ms = stats.total_extract_factors_time.as_secs_f64() * 1000.0;
        let exc_ms = stats.total_extract_contrib_time.as_secs_f64() * 1000.0;
        eprintln!();
        eprintln!("=== Sub-Phase Breakdown ===");
        eprintln!(
            "  Zeroing:      {:>10.2} ms  ({:>5.1}%)",
            zero_ms,
            zero_ms / factor_ms * 100.0
        );
        eprintln!(
            "  G2L setup:    {:>10.2} ms  ({:>5.1}%)",
            g2l_ms,
            g2l_ms / factor_ms * 100.0
        );
        eprintln!(
            "  Scatter:      {:>10.2} ms  ({:>5.1}%)",
            scatter_ms,
            scatter_ms / factor_ms * 100.0
        );
        eprintln!(
            "  Extend-add:   {:>10.2} ms  ({:>5.1}%)",
            ea_ms,
            ea_ms / factor_ms * 100.0
        );
        eprintln!(
            "  ExtractFactr: {:>10.2} ms  ({:>5.1}%)",
            exf_ms,
            exf_ms / factor_ms * 100.0
        );
        eprintln!(
            "  ExtractContr: {:>10.2} ms  ({:>5.1}%)",
            exc_ms,
            exc_ms / factor_ms * 100.0
        );
        let sub_accounted = zero_ms + g2l_ms + scatter_ms + ea_ms + exf_ms + exc_ms + kern_ms;
        eprintln!(
            "  Other:        {:>10.2} ms  ({:>5.1}%)",
            factor_ms - sub_accounted,
            (factor_ms - sub_accounted) / factor_ms * 100.0
        );
    }

    // --- Per-supernode top-N ---
    let per_sn = solver.per_supernode_stats().unwrap();
    print_top_supernodes(per_sn, 20, factor_ms);

    // --- Chrome Trace export ---
    if let Some(ref path) = trace_path {
        export_trace(per_sn, &name, path);
    }

    #[cfg(not(feature = "diagnostic"))]
    {
        eprintln!();
        eprintln!("NOTE: Build with --features diagnostic for per-supernode timing.");
    }
}

fn print_top_supernodes(per_sn: &[PerSupernodeStats], top_n: usize, _factor_ms: f64) {
    eprintln!();
    eprintln!(
        "=== Top {} Supernodes by Front Size ===",
        top_n.min(per_sn.len())
    );

    #[cfg(feature = "diagnostic")]
    {
        // Sort by total time (assembly + kernel + extraction)
        let mut ranked: Vec<(usize, &PerSupernodeStats, f64)> = per_sn
            .iter()
            .enumerate()
            .map(|(i, s)| {
                let total_us = s.assembly_time.as_secs_f64() * 1e6
                    + s.kernel_time.as_secs_f64() * 1e6
                    + s.extraction_time.as_secs_f64() * 1e6;
                (i, s, total_us)
            })
            .collect();
        ranked.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());

        let total_us: f64 = ranked.iter().map(|(_, _, t)| t).sum();

        eprintln!(
            "  {:>5} {:>6} {:>5} {:>4} {:>4} {:>10} {:>10} {:>10} {:>10} {:>6}",
            "snode",
            "front",
            "elim",
            "dlyd",
            "2x2",
            "asm_ms",
            "kern_ms",
            "ext_ms",
            "total_ms",
            "cum%",
        );
        eprintln!("  {}", "-".repeat(85));

        let mut cumulative_us = 0.0;
        for (_, s, t_us) in ranked.iter().take(top_n) {
            cumulative_us += t_us;
            let cum_pct = if total_us > 0.0 {
                cumulative_us / total_us * 100.0
            } else {
                0.0
            };
            eprintln!(
                "  {:>5} {:>6} {:>5} {:>4} {:>4} {:>10.3} {:>10.3} {:>10.3} {:>10.3} {:>5.1}%",
                s.snode_id,
                s.front_size,
                s.num_eliminated,
                s.num_delayed,
                s.num_2x2,
                s.assembly_time.as_secs_f64() * 1000.0,
                s.kernel_time.as_secs_f64() * 1000.0,
                s.extraction_time.as_secs_f64() * 1000.0,
                t_us / 1000.0,
                cum_pct,
            );
        }

        // Show how much the top N covers
        let top_n_us: f64 = ranked.iter().take(top_n).map(|(_, _, t)| t).sum();
        if total_us > 0.0 {
            eprintln!();
            eprintln!(
                "  Top {} of {} supernodes account for {:.1}% of factor time.",
                top_n.min(per_sn.len()),
                per_sn.len(),
                top_n_us / total_us * 100.0,
            );
        }
    }

    #[cfg(not(feature = "diagnostic"))]
    {
        // Without diagnostic, sort by front size
        let mut ranked: Vec<&PerSupernodeStats> = per_sn.iter().collect();
        ranked.sort_by(|a, b| b.front_size.cmp(&a.front_size));

        eprintln!(
            "  {:>5} {:>6} {:>5} {:>4} {:>4} {:>8}",
            "snode", "front", "elim", "dlyd", "2x2", "max|L|",
        );
        eprintln!("  {}", "-".repeat(45));

        for s in ranked.iter().take(top_n) {
            eprintln!(
                "  {:>5} {:>6} {:>5} {:>4} {:>4} {:>8.2e}",
                s.snode_id, s.front_size, s.num_eliminated, s.num_delayed, s.num_2x2, s.max_l_entry,
            );
        }
    }
}

#[cfg(feature = "diagnostic")]
fn export_trace(per_sn: &[PerSupernodeStats], matrix_name: &str, path: &str) {
    use serde_json::json;

    let pid = std::process::id();
    let mut events = Vec::new();
    let mut offset_us = 0.0;

    for s in per_sn {
        // Assembly
        let asm_us = s.assembly_time.as_secs_f64() * 1e6;
        if asm_us > 0.0 {
            events.push(json!({
                "name": format!("s{}_asm", s.snode_id),
                "cat": "assembly",
                "ph": "X",
                "ts": offset_us,
                "dur": asm_us,
                "pid": pid,
                "tid": 0,
                "args": { "snode": s.snode_id, "front_size": s.front_size }
            }));
            offset_us += asm_us;
        }

        // Kernel
        let kern_us = s.kernel_time.as_secs_f64() * 1e6;
        if kern_us > 0.0 {
            events.push(json!({
                "name": format!("s{}_kern", s.snode_id),
                "cat": "kernel",
                "ph": "X",
                "ts": offset_us,
                "dur": kern_us,
                "pid": pid,
                "tid": 0,
                "args": {
                    "snode": s.snode_id,
                    "front_size": s.front_size,
                    "eliminated": s.num_eliminated,
                    "delayed": s.num_delayed,
                }
            }));
            offset_us += kern_us;
        }

        // Extraction
        let ext_us = s.extraction_time.as_secs_f64() * 1e6;
        if ext_us > 0.0 {
            events.push(json!({
                "name": format!("s{}_ext", s.snode_id),
                "cat": "extraction",
                "ph": "X",
                "ts": offset_us,
                "dur": ext_us,
                "pid": pid,
                "tid": 0,
                "args": { "snode": s.snode_id, "front_size": s.front_size }
            }));
            offset_us += ext_us;
        }
    }

    let trace = json!({
        "traceEvents": events,
        "metadata": { "matrix": matrix_name }
    });
    let json_str = serde_json::to_string_pretty(&trace).unwrap();
    std::fs::write(path, &json_str).unwrap_or_else(|e| {
        eprintln!("Failed to write trace to {}: {}", path, e);
        std::process::exit(1);
    });
    eprintln!();
    eprintln!("Chrome Trace written to: {}", path);
    eprintln!("  Open in https://ui.perfetto.dev or chrome://tracing");
}

#[cfg(not(feature = "diagnostic"))]
fn export_trace(_per_sn: &[PerSupernodeStats], _matrix_name: &str, _path: &str) {
    eprintln!("Chrome Trace export requires --features diagnostic");
    std::process::exit(1);
}

fn list_matrices() {
    let all = registry::load_registry().unwrap_or_else(|e| {
        eprintln!("Failed to load registry: {}", e);
        std::process::exit(1);
    });

    eprintln!(
        "{:<30} {:>8} {:>10} {:>6} {:>20}",
        "Name", "n", "nnz", "CI?", "Category"
    );
    eprintln!("{}", "-".repeat(80));
    for m in &all {
        eprintln!(
            "{:<30} {:>8} {:>10} {:>6} {:>20}",
            m.name,
            m.size,
            m.nnz,
            if m.ci_subset { "yes" } else { "" },
            m.category,
        );
    }
    eprintln!("\n{} matrices total", all.len());
}

fn print_usage() {
    eprintln!("Usage: profile_matrix [OPTIONS] <MATRIX_NAME>");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  <MATRIX_NAME>        Registry name (from metadata.json)");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  --path <FILE>        Load matrix from .mtx file path instead");
    eprintln!("  --trace <FILE>       Export Chrome Trace JSON to file");
    eprintln!("  --list               List all available matrices");
    eprintln!("  -h, --help           Show this help");
}
