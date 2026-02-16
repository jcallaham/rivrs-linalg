//! Report front size statistics for each CI matrix.
//!
//! For each matrix in the CI subset, performs symbolic analysis with METIS ordering
//! and reports supernode statistics including front sizes. The "front size"
//! for a supernode is the total dimension of its frontal matrix:
//! `supernode_width + off_diagonal_pattern_length`.
//!
//! Usage: cargo run --example front_sizes --release

use faer::sparse::linalg::cholesky::SymmetricOrdering;
use rivrs_sparse::aptp::ordering::metis_ordering;
use rivrs_sparse::aptp::AptpSymbolic;
use rivrs_sparse::io::registry;

fn main() {
    let all = registry::load_registry().expect("failed to load registry");
    let ci: Vec<_> = all.iter().filter(|m| m.ci_subset).collect();

    println!(
        "{:<30} {:>8} {:>10} {:>8} {:>8} {:>8} {:>8}",
        "Matrix", "n", "nnz", "n_sn", "max_fr", "med_fr", "lg_fronts"
    );
    println!("{}", "-".repeat(85));

    for meta in &ci {
        let test = match registry::load_test_matrix(&meta.name) {
            Ok(Some(t)) => t,
            _ => {
                println!("{:<30} SKIP", meta.name);
                continue;
            }
        };
        let a = &test.matrix;
        let n = a.nrows();

        // Compute METIS ordering
        let perm = match metis_ordering(a.symbolic()) {
            Ok(p) => p,
            Err(e) => {
                println!(
                    "{:<30} {:>8} {:>10} METIS FAILED: {}",
                    meta.name, n, meta.nnz, e
                );
                continue;
            }
        };

        // Perform symbolic analysis with METIS ordering
        let symbolic =
            match AptpSymbolic::analyze(a.symbolic(), SymmetricOrdering::Custom(perm.as_ref())) {
                Ok(s) => s,
                Err(e) => {
                    println!(
                        "{:<30} {:>8} {:>10} ANALYZE FAILED: {}",
                        meta.name, n, meta.nnz, e
                    );
                    continue;
                }
            };

        let n_supernodes = match symbolic.n_supernodes() {
            Some(ns) => ns,
            None => {
                println!(
                    "{:<30} {:>8} {:>10} simplicial (no supernodes)",
                    meta.name, n, meta.nnz
                );
                continue;
            }
        };

        let begin = symbolic.supernode_begin().unwrap();
        let end = symbolic.supernode_end().unwrap();

        // Compute front sizes: supernode_width + off-diagonal pattern length
        let mut front_sizes: Vec<usize> = (0..n_supernodes)
            .map(|s| {
                let sn_width = end[s] - begin[s];
                let pattern_len = symbolic.supernode_pattern(s).unwrap().len();
                sn_width + pattern_len
            })
            .collect();

        front_sizes.sort();
        let max_front = front_sizes.last().copied().unwrap_or(0);
        let median_front = if front_sizes.is_empty() {
            0
        } else {
            front_sizes[front_sizes.len() / 2]
        };
        let large_fronts = front_sizes.iter().filter(|&&f| f > 500).count();

        println!(
            "{:<30} {:>8} {:>10} {:>8} {:>8} {:>8} {:>8}",
            meta.name, n, meta.nnz, n_supernodes, max_front, median_front, large_fronts
        );
    }
}
