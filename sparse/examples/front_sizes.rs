//! Report front size statistics for all SuiteSparse matrices.
//!
//! For each matrix, performs symbolic analysis with METIS ordering
//! and reports supernode statistics including front sizes. The "front size"
//! for a supernode is the total dimension of its frontal matrix:
//! `supernode_width + off_diagonal_pattern_length`.
//!
//! Usage: cargo run --example front_sizes --release

use faer::sparse::linalg::cholesky::SymmetricOrdering;
use rivrs_sparse::aptp::ordering::metis_ordering;
use rivrs_sparse::aptp::AptpSymbolic;
use rivrs_sparse::io::registry;

fn percentile(sorted: &[usize], p: f64) -> usize {
    if sorted.is_empty() {
        return 0;
    }
    let idx = ((sorted.len() as f64 - 1.0) * p) as usize;
    sorted[idx.min(sorted.len() - 1)]
}

fn main() {
    let all = registry::load_registry().expect("failed to load registry");
    let mut matrices: Vec<_> = all.iter().collect();
    matrices.sort_by_key(|m| m.size);

    println!(
        "{:<40} {:>8} {:>10} {:>7} {:>7} {:>7} {:>7} {:>7} {:>6}",
        "Matrix", "n", "nnz", "n_sn", "med_fr", "p90_fr", "p99_fr", "max_fr", ">500"
    );
    println!("{}", "-".repeat(105));

    let mut all_max_fronts: Vec<(String, usize, usize)> = Vec::new(); // (name, n, max_front)

    for meta in &matrices {
        let test = match registry::load_test_matrix(&meta.name) {
            Ok(Some(t)) => t,
            _ => {
                continue; // silently skip missing matrices
            }
        };
        let a = &test.matrix;
        let n = a.nrows();

        let perm = match metis_ordering(a.symbolic()) {
            Ok(p) => p,
            Err(e) => {
                println!(
                    "{:<40} {:>8} {:>10} METIS FAILED: {}",
                    meta.name, n, meta.nnz, e
                );
                continue;
            }
        };

        let symbolic =
            match AptpSymbolic::analyze(a.symbolic(), SymmetricOrdering::Custom(perm.as_ref())) {
                Ok(s) => s,
                Err(e) => {
                    println!(
                        "{:<40} {:>8} {:>10} ANALYZE FAILED: {}",
                        meta.name, n, meta.nnz, e
                    );
                    continue;
                }
            };

        let n_supernodes = match symbolic.n_supernodes() {
            Some(ns) => ns,
            None => {
                println!(
                    "{:<40} {:>8} {:>10} simplicial",
                    meta.name, n, meta.nnz
                );
                all_max_fronts.push((meta.name.clone(), n, 0));
                continue;
            }
        };

        let begin = symbolic.supernode_begin().unwrap();
        let end = symbolic.supernode_end().unwrap();

        let mut front_sizes: Vec<usize> = (0..n_supernodes)
            .map(|s| {
                let sn_width = end[s] - begin[s];
                let pattern_len = symbolic.supernode_pattern(s).unwrap().len();
                sn_width + pattern_len
            })
            .collect();

        front_sizes.sort();
        let max_front = front_sizes.last().copied().unwrap_or(0);
        let med_front = percentile(&front_sizes, 0.5);
        let p90_front = percentile(&front_sizes, 0.9);
        let p99_front = percentile(&front_sizes, 0.99);
        let large_fronts = front_sizes.iter().filter(|&&f| f > 500).count();

        println!(
            "{:<40} {:>8} {:>10} {:>7} {:>7} {:>7} {:>7} {:>7} {:>6}",
            meta.name, n, meta.nnz, n_supernodes, med_front, p90_front, p99_front, max_front, large_fronts
        );

        all_max_fronts.push((meta.name.clone(), n, max_front));
    }

    // Summary statistics
    println!("\n=== Summary ===");
    let total = all_max_fronts.len();
    let simplicial = all_max_fronts.iter().filter(|(_, _, mf)| *mf == 0).count();
    let small = all_max_fronts.iter().filter(|(_, _, mf)| *mf > 0 && *mf <= 256).count();
    let medium = all_max_fronts.iter().filter(|(_, _, mf)| *mf > 256 && *mf <= 1000).count();
    let large = all_max_fronts.iter().filter(|(_, _, mf)| *mf > 1000 && *mf <= 5000).count();
    let huge = all_max_fronts.iter().filter(|(_, _, mf)| *mf > 5000).count();

    println!("Total matrices analyzed: {}", total);
    println!("  Simplicial (no supernodes):  {:>3}", simplicial);
    println!("  max_front ≤ 256:             {:>3}  (single-level fine)", small);
    println!("  max_front 257–1000:          {:>3}  (borderline)", medium);
    println!("  max_front 1001–5000:         {:>3}  (needs two-level)", large);
    println!("  max_front > 5000:            {:>3}  (critical for two-level)", huge);

    println!("\nTop 10 by max front size:");
    let mut sorted = all_max_fronts.clone();
    sorted.sort_by(|a, b| b.2.cmp(&a.2));
    for (name, n, mf) in sorted.iter().take(10) {
        println!("  {:<40} n={:>8}  max_front={:>6}", name, n, mf);
    }
}
