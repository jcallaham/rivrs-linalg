//! Workload distribution analysis for Phase 8.2 parallelism strategy.
//!
//! Analyzes per-supernode timing data from the SuiteSparse suite to classify
//! matrices by workload distribution and recommend parallelism strategies.
//!
//! Usage:
//!   cargo run --example workload_analysis --features diagnostic --release

use faer::Col;
use serde::{Deserialize, Serialize};

use rivrs_sparse::aptp::{AnalyzeOptions, FactorOptions, SparseLDLT};
use rivrs_sparse::io::registry;

/// Parallelism strategy classification.
#[derive(Debug, Clone, Serialize, Deserialize)]
enum ParallelismClass {
    /// Many small independent fronts; tree-level parallelism most beneficial.
    TreeLevel,
    /// Few large fronts dominate; intra-node BLAS-3 parallelism most beneficial.
    IntraNode,
    /// Significant time in both small and large fronts; both strategies needed.
    Mixed,
}

impl std::fmt::Display for ParallelismClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParallelismClass::TreeLevel => write!(f, "TreeLevel"),
            ParallelismClass::IntraNode => write!(f, "IntraNode"),
            ParallelismClass::Mixed => write!(f, "Mixed"),
        }
    }
}

/// Per-matrix workload profile.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct WorkloadProfile {
    matrix_name: String,
    num_supernodes: usize,
    max_front_size: usize,
    total_factor_us: f64,
    top_10pct_time_fraction: f64,
    top_1_front_time_fraction: f64,
    front_size_histogram: Vec<(String, usize)>,
    time_by_front_size: Vec<(String, f64)>,
    parallelism_recommendation: ParallelismClass,
}

fn classify(top_10pct_fraction: f64) -> ParallelismClass {
    if top_10pct_fraction > 0.80 {
        ParallelismClass::IntraNode
    } else if top_10pct_fraction < 0.30 {
        ParallelismClass::TreeLevel
    } else {
        ParallelismClass::Mixed
    }
}

fn front_size_bucket(size: usize) -> &'static str {
    match size {
        0..=10 => "1-10",
        11..=50 => "11-50",
        51..=100 => "51-100",
        101..=500 => "101-500",
        501..=1000 => "501-1k",
        _ => "1k+",
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let ci_only = args.iter().any(|a| a == "--ci-only");

    let all = registry::load_registry().expect("Failed to load registry");
    let matrices: Vec<_> = all
        .iter()
        .filter(|m| m.source == "suitesparse" && (!ci_only || m.ci_subset))
        .collect();

    eprintln!(
        "Workload analysis for {} SuiteSparse matrices{}",
        matrices.len(),
        if ci_only { " (CI subset)" } else { "" }
    );
    eprintln!();

    println!(
        "{:<30} {:>8} {:>6} {:>10} {:>10} {:>10} {:>12}",
        "Matrix", "n", "snodes", "max_front", "top10%_t", "top1_t", "class"
    );
    println!("{}", "-".repeat(100));

    let mut profiles = Vec::new();
    let mut class_counts = [0usize; 3]; // TreeLevel, IntraNode, Mixed

    for meta in &matrices {
        let test = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(t)) => t,
            _ => continue,
        };
        let a = &test.matrix;
        let n = a.nrows();

        // Create RHS
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

        // Analyze + factor
        let opts = AnalyzeOptions::default();
        let mut solver = match SparseLDLT::analyze_with_matrix(a, &opts) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let factor_opts = FactorOptions::default();
        if solver.factor(a, &factor_opts).is_err() {
            continue;
        }

        // Solve for backward error
        let scratch = solver.solve_scratch(1);
        let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
        let stack = faer::dyn_stack::MemStack::new(&mut mem);
        let _x = solver.solve(&b, stack).ok();

        let stats = solver.stats().unwrap();
        let per_sn = solver.per_supernode_stats().unwrap();

        // Compute per-supernode timing (from diagnostic fields)
        let mut sn_times: Vec<(usize, f64)> = per_sn
            .iter()
            .map(|s| {
                #[cfg(feature = "diagnostic")]
                let time_us = s.assembly_time.as_secs_f64() * 1e6
                    + s.kernel_time.as_secs_f64() * 1e6
                    + s.extraction_time.as_secs_f64() * 1e6;
                #[cfg(not(feature = "diagnostic"))]
                let time_us = 0.0;
                (s.front_size, time_us)
            })
            .collect();

        let total_time_us: f64 = sn_times.iter().map(|(_, t)| t).sum();

        // Sort by front_size descending
        sn_times.sort_by(|a, b| b.0.cmp(&a.0));

        // Top 10% by size
        let top_10_count = (per_sn.len() as f64 * 0.10).ceil() as usize;
        let top_10_time: f64 = sn_times.iter().take(top_10_count).map(|(_, t)| t).sum();
        let top_10_frac = if total_time_us > 0.0 {
            top_10_time / total_time_us
        } else {
            0.0
        };

        // Top 1 front
        let top_1_time = sn_times.first().map_or(0.0, |(_, t)| *t);
        let top_1_frac = if total_time_us > 0.0 {
            top_1_time / total_time_us
        } else {
            0.0
        };

        // Build histogram
        let buckets = ["1-10", "11-50", "51-100", "101-500", "501-1k", "1k+"];
        let mut hist_count = vec![0usize; buckets.len()];
        let mut hist_time = vec![0.0f64; buckets.len()];
        for (fs, t) in &sn_times {
            let bucket = front_size_bucket(*fs);
            let idx = buckets.iter().position(|&b| b == bucket).unwrap_or(0);
            hist_count[idx] += 1;
            hist_time[idx] += t;
        }

        let front_size_histogram: Vec<(String, usize)> = buckets
            .iter()
            .zip(&hist_count)
            .map(|(&b, &c)| (b.to_string(), c))
            .collect();
        let time_by_front_size: Vec<(String, f64)> = buckets
            .iter()
            .zip(&hist_time)
            .map(|(&b, t)| (b.to_string(), *t))
            .collect();

        let class = classify(top_10_frac);
        match &class {
            ParallelismClass::TreeLevel => class_counts[0] += 1,
            ParallelismClass::IntraNode => class_counts[1] += 1,
            ParallelismClass::Mixed => class_counts[2] += 1,
        }

        println!(
            "{:<30} {:>8} {:>6} {:>10} {:>9.1}% {:>9.1}% {:>12}",
            meta.name,
            n,
            per_sn.len(),
            stats.max_front_size,
            top_10_frac * 100.0,
            top_1_frac * 100.0,
            class
        );

        profiles.push(WorkloadProfile {
            matrix_name: meta.name.clone(),
            num_supernodes: per_sn.len(),
            max_front_size: stats.max_front_size,
            total_factor_us: total_time_us,
            top_10pct_time_fraction: top_10_frac,
            top_1_front_time_fraction: top_1_frac,
            front_size_histogram,
            time_by_front_size,
            parallelism_recommendation: class,
        });
    }

    // Summary
    println!("\n--- Classification Summary ---");
    println!("TreeLevel: {} matrices", class_counts[0]);
    println!("IntraNode: {} matrices", class_counts[1]);
    println!("Mixed:     {} matrices", class_counts[2]);
    println!("Total:     {} matrices", profiles.len());

    println!("\n--- Parallelism Recommendation for Phase 8.2 ---");
    if class_counts[1] > class_counts[0] && class_counts[1] > class_counts[2] {
        println!(
            "Primary strategy: Intra-node BLAS-3 parallelism (majority of matrices have concentrated workload)"
        );
    } else if class_counts[0] > class_counts[1] && class_counts[0] > class_counts[2] {
        println!(
            "Primary strategy: Tree-level parallelism (majority of matrices have distributed workload)"
        );
    } else {
        println!("Mixed strategy: Both tree-level and intra-node parallelism needed");
    }

    // Export profiles as JSON to stdout
    let json = serde_json::to_string_pretty(&profiles).expect("serialize");
    eprintln!("\nWorkload profiles JSON:");
    eprintln!("{json}");
}
