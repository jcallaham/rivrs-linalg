//! Per-supernode diagnostic statistics for sparse LDL^T factorization.
//!
//! Loads a matrix, factors with MatchOrderMetis, and prints per-supernode
//! stats (front size, delays, pivot types, max L entry) to help identify
//! factorization behavior and compare with reference solvers.
//!
//! Usage:
//!   cargo run --release --example supernode_stats -- <path/to/matrix.mtx>
//!
//! If no matrix path is given, defaults to d_pretok.

use rivrs_sparse::aptp::{AnalyzeOptions, FactorOptions, OrderingStrategy, SparseLDLT};
use rivrs_sparse::io::mtx::load_mtx;

fn main() {
    let path = std::env::args().nth(1).unwrap_or_else(|| {
        "test-data/suitesparse/hard-indefinite/d_pretok/d_pretok.mtx".to_string()
    });
    eprintln!("Loading {}...", path);
    let matrix = load_mtx(std::path::Path::new(&path)).unwrap();
    let n = matrix.nrows();
    let nnz = matrix.val().len();
    eprintln!("n={}, nnz={}", n, nnz);

    // Analyze
    let analyze_opts = AnalyzeOptions {
        ordering: OrderingStrategy::MatchOrderMetis,
    };
    eprintln!("Analyzing (MatchOrderMetis)...");
    let mut solver = SparseLDLT::analyze_with_matrix(&matrix, &analyze_opts).unwrap();

    // Factor
    let factor_opts = FactorOptions::default();
    eprintln!("Factoring...");
    solver.factor(&matrix, &factor_opts).unwrap();

    // Retrieve stats
    let agg_stats = solver.stats().unwrap();
    let per_sn = solver.per_supernode_stats().unwrap();

    // Print TSV header (stdout for piping)
    println!(
        "snode_id\tfront_size\tnum_fully_summed\tnum_eliminated\tnum_delayed\tnum_1x1\tnum_2x2\tmax_l_entry"
    );

    for s in per_sn {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6e}",
            s.snode_id,
            s.front_size,
            s.num_fully_summed,
            s.num_eliminated,
            s.num_delayed,
            s.num_1x1,
            s.num_2x2,
            s.max_l_entry,
        );
    }

    // Summary (stderr so it doesn't mix with TSV)
    eprintln!();
    eprintln!("=== Summary ===");
    eprintln!("Total supernodes:     {}", per_sn.len());
    eprintln!("Total 1x1 pivots:     {}", agg_stats.total_1x1_pivots);
    eprintln!("Total 2x2 pivots:     {}", agg_stats.total_2x2_pivots);
    eprintln!("Total delayed:        {}", agg_stats.total_delayed);
    eprintln!("Zero pivots:          {}", agg_stats.zero_pivots);
    eprintln!("Max front size:       {}", agg_stats.max_front_size);

    // Supernodes with delays
    let mut delayed: Vec<_> = per_sn.iter().filter(|s| s.num_delayed > 0).collect();
    delayed.sort_by(|a, b| b.num_delayed.cmp(&a.num_delayed));

    eprintln!();
    eprintln!("=== Supernodes with delays (sorted by delay count) ===");
    eprintln!(
        "{:>8}  {:>10}  {:>6}  {:>6}  {:>7}  {:>5}  {:>5}  {:>12}",
        "snode_id", "front_size", "k", "ne", "delayed", "1x1", "2x2", "max_l_entry"
    );
    for s in &delayed {
        eprintln!(
            "{:>8}  {:>10}  {:>6}  {:>6}  {:>7}  {:>5}  {:>5}  {:>12.4e}",
            s.snode_id,
            s.front_size,
            s.num_fully_summed,
            s.num_eliminated,
            s.num_delayed,
            s.num_1x1,
            s.num_2x2,
            s.max_l_entry,
        );
    }
    eprintln!("Total supernodes with delays: {}", delayed.len());

    // Largest fronts
    let mut by_front_size: Vec<_> = per_sn.iter().collect();
    by_front_size.sort_by(|a, b| b.front_size.cmp(&a.front_size));

    eprintln!();
    eprintln!("=== Top 20 largest fronts ===");
    eprintln!(
        "{:>8}  {:>10}  {:>6}  {:>6}  {:>7}  {:>5}  {:>5}  {:>12}",
        "snode_id", "front_size", "k", "ne", "delayed", "1x1", "2x2", "max_l_entry"
    );
    for s in by_front_size.iter().take(20) {
        eprintln!(
            "{:>8}  {:>10}  {:>6}  {:>6}  {:>7}  {:>5}  {:>5}  {:>12.4e}",
            s.snode_id,
            s.front_size,
            s.num_fully_summed,
            s.num_eliminated,
            s.num_delayed,
            s.num_1x1,
            s.num_2x2,
            s.max_l_entry,
        );
    }

    // Front size distribution
    eprintln!();
    eprintln!("=== Front size distribution ===");
    let mut buckets = std::collections::BTreeMap::new();
    for s in per_sn {
        let bucket = if s.front_size <= 10 {
            "1-10"
        } else if s.front_size <= 50 {
            "11-50"
        } else if s.front_size <= 100 {
            "51-100"
        } else if s.front_size <= 256 {
            "101-256"
        } else if s.front_size <= 512 {
            "257-512"
        } else if s.front_size <= 1024 {
            "513-1024"
        } else {
            "1025+"
        };
        *buckets.entry(bucket).or_insert(0usize) += 1;
    }
    for (bucket, count) in &buckets {
        eprintln!("  {:>10}: {}", bucket, count);
    }
}
