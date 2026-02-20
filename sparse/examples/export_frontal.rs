//! Export an assembled frontal matrix (pre-factorization) for SPRAL comparison.
//!
//! Loads a sparse matrix, performs symbolic analysis with MC64+METIS ordering,
//! assembles the frontal matrix for a target supernode (default: largest front),
//! and writes it to a file that can be read by the SPRAL dense APTP driver.
//!
//! # Output Format
//!
//! ```text
//! m k
//! a[0,0] a[0,1] ... a[0,m-1]
//! a[1,0] a[1,1] ... a[1,m-1]
//! ...
//! a[m-1,0] a[m-1,1] ... a[m-1,m-1]
//! ```
//!
//! The matrix is stored as full m x m (both triangles), row by row.
//! Only the lower triangle is assembled; the upper triangle is mirrored.
//!
//! # Usage
//!
//! ```bash
//! # Export largest front for d_pretok:
//! cargo run --example export_frontal --release
//!
//! # Export specific supernode:
//! SNODE=42 cargo run --example export_frontal --release
//!
//! # Use a different matrix:
//! MATRIX=test-data/suitesparse/hard-indefinite/stokes128/stokes128.mtx \
//!   cargo run --example export_frontal --release
//! ```

use std::env;
use std::fs::File;
use std::io::Write;

use faer::sparse::linalg::cholesky::SymmetricOrdering;

use rivrs_sparse::aptp::ordering::metis_ordering;
use rivrs_sparse::aptp::symbolic::AptpSymbolic;
use rivrs_sparse::aptp::{AptpNumeric, AptpOptions, match_order_metis};
use rivrs_sparse::io::mtx::load_mtx;

fn main() {
    // --- Configuration from environment variables ---
    let matrix_path = env::var("MATRIX").unwrap_or_else(|_| {
        "test-data/suitesparse/hard-indefinite/d_pretok/d_pretok.mtx".to_string()
    });
    let output_path = env::var("OUTPUT").unwrap_or_else(|_| "/tmp/frontal_matrix.txt".to_string());
    let target_snode: Option<usize> = env::var("SNODE").ok().and_then(|s| s.parse().ok());
    let use_matching = env::var("NO_MATCHING").is_err(); // default: use MC64+METIS

    println!("Matrix:       {}", matrix_path);
    println!("Output:       {}", output_path);
    println!(
        "Target snode: {}",
        target_snode
            .map(|s| s.to_string())
            .unwrap_or_else(|| "largest front".to_string())
    );
    println!(
        "Ordering:     {}",
        if use_matching {
            "MatchOrderMetis (MC64+METIS)"
        } else {
            "METIS only"
        }
    );

    // --- Load matrix ---
    println!("\nLoading matrix...");
    let matrix = load_mtx(std::path::Path::new(&matrix_path)).expect("Failed to load matrix");
    let n = matrix.nrows();
    let nnz = matrix.val().len();
    println!("  n = {}, nnz = {}", n, nnz);

    // --- Symbolic analysis ---
    println!("Running symbolic analysis...");
    let (symbolic, scaling) = if use_matching {
        let mo_result = match_order_metis(&matrix).expect("MC64+METIS failed");
        let ordering_perm = mo_result.ordering;
        let sym = AptpSymbolic::analyze(
            matrix.symbolic(),
            SymmetricOrdering::Custom(ordering_perm.as_ref()),
        )
        .expect("Symbolic analysis failed");

        // Transform scaling to elimination order
        let (perm_fwd, _) = sym.perm_vecs();
        let elim_scaling: Vec<f64> = (0..n).map(|i| mo_result.scaling[perm_fwd[i]]).collect();
        (sym, Some(elim_scaling))
    } else {
        let perm = metis_ordering(matrix.symbolic()).expect("METIS ordering failed");
        let sym =
            AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Custom(perm.as_ref()))
                .expect("Symbolic analysis failed");
        (sym, None)
    };

    // Print supernode statistics
    if let Some(ns) = symbolic.n_supernodes() {
        let begin = symbolic.supernode_begin().unwrap();
        let end = symbolic.supernode_end().unwrap();
        let mut front_sizes: Vec<(usize, usize)> = (0..ns)
            .map(|s| {
                let sn_width = end[s] - begin[s];
                let pattern_len = symbolic.supernode_pattern(s).unwrap().len();
                (s, sn_width + pattern_len)
            })
            .collect();
        front_sizes.sort_by(|a, b| b.1.cmp(&a.1));

        println!("  {} supernodes", ns);
        println!("  Top 5 fronts by size:");
        for (s, sz) in front_sizes.iter().take(5) {
            let sn_width = end[*s] - begin[*s];
            let pat_len = symbolic.supernode_pattern(*s).unwrap().len();
            println!(
                "    snode {:>5}: front_size={:>5} (cols={:>4}, pattern={:>4})",
                s, sz, sn_width, pat_len
            );
        }
    }

    // --- Export assembled frontal matrix ---
    println!("\nAssembling frontal matrix for target supernode...");
    let options = AptpOptions::default();
    let (frontal, m, k, row_indices) = AptpNumeric::export_assembled_frontal(
        &symbolic,
        &matrix,
        &options,
        scaling.as_deref(),
        target_snode,
    )
    .expect("Failed to export assembled frontal matrix");

    println!("  Frontal matrix: m={}, k={}", m, k);
    println!(
        "  Row indices (first 20): {:?}",
        &row_indices[..20.min(row_indices.len())]
    );

    // Mirror lower triangle to get full symmetric matrix
    let mut full = frontal.clone();
    for i in 0..m {
        for j in (i + 1)..m {
            full[(i, j)] = full[(j, i)];
        }
    }

    // Some diagnostic stats
    let mut max_val = 0.0f64;
    let mut min_diag = f64::INFINITY;
    let mut max_diag = 0.0f64;
    let mut zero_diag = 0usize;
    for i in 0..m {
        let d = full[(i, i)].abs();
        if d == 0.0 {
            zero_diag += 1;
        }
        min_diag = min_diag.min(d);
        max_diag = max_diag.max(d);
        for j in 0..m {
            max_val = max_val.max(full[(i, j)].abs());
        }
    }
    println!("  max |a_ij| = {:.6e}", max_val);
    println!(
        "  diagonal: min |a_ii| = {:.6e}, max |a_ii| = {:.6e}, zeros = {}",
        min_diag, max_diag, zero_diag
    );

    // Norm of the matrix (Frobenius)
    let mut fro_sq = 0.0f64;
    for i in 0..m {
        for j in 0..m {
            fro_sq += full[(i, j)] * full[(i, j)];
        }
    }
    let fro_norm = fro_sq.sqrt();
    println!("  ||F||_F = {:.6e}", fro_norm);

    // --- Write to file ---
    println!("\nWriting to {}...", output_path);
    let mut f = File::create(&output_path).expect("Failed to create output file");

    // Header: m k
    writeln!(f, "{} {}", m, k).expect("Write failed");

    // Matrix entries row by row (full symmetric)
    for i in 0..m {
        let row: Vec<String> = (0..m).map(|j| format!("{:.17e}", full[(i, j)])).collect();
        writeln!(f, "{}", row.join(" ")).expect("Write failed");
    }

    let file_size = std::fs::metadata(&output_path)
        .map(|md| md.len())
        .unwrap_or(0);
    println!(
        "  Done. File size: {:.1} MB",
        file_size as f64 / (1024.0 * 1024.0)
    );

    // Also write metadata
    let meta_path = format!("{}.meta", output_path);
    let mut mf = File::create(&meta_path).expect("Failed to create metadata file");
    writeln!(mf, "matrix: {}", matrix_path).unwrap();
    writeln!(mf, "m: {}", m).unwrap();
    writeln!(mf, "k: {}", k).unwrap();
    writeln!(
        mf,
        "ordering: {}",
        if use_matching {
            "MatchOrderMetis"
        } else {
            "METIS"
        }
    )
    .unwrap();
    writeln!(
        mf,
        "target_snode: {}",
        target_snode
            .map(|s| s.to_string())
            .unwrap_or_else(|| "auto (largest front)".to_string())
    )
    .unwrap();
    writeln!(mf, "fro_norm: {:.17e}", fro_norm).unwrap();
    writeln!(mf, "max_abs: {:.17e}", max_val).unwrap();
    writeln!(mf, "zero_diags: {}", zero_diag).unwrap();
    println!("  Metadata: {}", meta_path);

    // --- Also factor with our APTP kernel for comparison ---
    if env::var("NO_FACTOR").is_err() {
        println!("\n=== Factoring with our APTP kernel for comparison ===");

        // Clone the assembled matrix (lower triangle only, same as what we pass to factor)
        let mut our_frontal = frontal.clone();
        // Mirror to get full lower triangle populated (our factor expects lower triangle)
        for i in 0..m {
            for j in (i + 1)..m {
                our_frontal[(i, j)] = our_frontal[(j, i)];
            }
        }

        let result = rivrs_sparse::aptp::aptp_factor_in_place(our_frontal.as_mut(), k, &options);
        match result {
            Ok(res) => {
                println!("  num_eliminated: {}", res.num_eliminated);
                println!("  num_1x1: {}", res.stats.num_1x1);
                println!("  num_2x2: {}", res.stats.num_2x2);
                println!("  num_delayed: {}", res.stats.num_delayed);
                println!("  max_l_entry: {:.6e}", res.stats.max_l_entry);

                if res.num_eliminated < k {
                    println!(
                        "  WARNING: {} columns delayed (SPRAL eliminates all {})",
                        k - res.num_eliminated,
                        k
                    );
                }
            }
            Err(e) => {
                println!("  ERROR: {:?}", e);
            }
        }
    }

    println!("\nTo factor with SPRAL:");
    println!(
        "  /tmp/spral_dense_factor -v < {} 2>&1 | tail -20",
        output_path
    );
}
