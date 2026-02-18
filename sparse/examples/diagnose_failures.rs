//! Diagnose backward error failures on SuiteSparse matrices.
//!
//! Usage: cargo run --example diagnose_failures --release

use faer::Col;
use rivrs_sparse::aptp::{
    AnalyzeOptions, FactorOptions, FactorizationStats, OrderingStrategy, SparseLDLT,
};
use rivrs_sparse::io::registry;
use rivrs_sparse::validate::sparse_backward_error;

fn sparse_matvec(a: &faer::sparse::SparseColMat<usize, f64>, x: &[f64]) -> Vec<f64> {
    let n = a.nrows();
    let sym = a.symbolic();
    let cp = sym.col_ptr();
    let ri = sym.row_idx();
    let vals = a.val();
    let mut result = vec![0.0f64; n];
    for j in 0..n {
        for idx in cp[j]..cp[j + 1] {
            result[ri[idx]] += vals[idx] * x[j];
        }
    }
    result
}

fn diagnose(name: &str, a: &faer::sparse::SparseColMat<usize, f64>) {
    let n = a.nrows();
    let nnz = a.compute_nnz();

    // Check diagonal properties
    let sym = a.symbolic();
    let cp = sym.col_ptr();
    let ri = sym.row_idx();
    let vals = a.val();

    let mut zero_diag = 0usize;
    let mut neg_diag = 0usize;
    let mut pos_diag = 0usize;
    let mut min_abs_diag = f64::INFINITY;
    let mut max_abs_diag = 0.0f64;

    for j in 0..n {
        let mut found_diag = false;
        for idx in cp[j]..cp[j + 1] {
            if ri[idx] == j {
                let v = vals[idx];
                found_diag = true;
                if v == 0.0 {
                    zero_diag += 1;
                } else if v < 0.0 {
                    neg_diag += 1;
                } else {
                    pos_diag += 1;
                }
                let abs_v = v.abs();
                if abs_v > 0.0 && abs_v < min_abs_diag {
                    min_abs_diag = abs_v;
                }
                if abs_v > max_abs_diag {
                    max_abs_diag = abs_v;
                }
                break;
            }
        }
        if !found_diag {
            zero_diag += 1; // missing diagonal = zero
        }
    }

    println!("\n{} (n={}, nnz={})", name, n, nnz);
    println!(
        "  Diagonal: {} pos, {} neg, {} zero, min_abs={:.2e}, max_abs={:.2e}, ratio={:.2e}",
        pos_diag,
        neg_diag,
        zero_diag,
        min_abs_diag,
        max_abs_diag,
        max_abs_diag / min_abs_diag
    );

    // Build RHS
    let x_exact: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
    let b_vec = sparse_matvec(a, &x_exact);
    let b = Col::from_fn(n, |i| b_vec[i]);

    // Try MatchOrderMetis (default)
    solve_and_report("MatchOrderMetis", a, &b, OrderingStrategy::MatchOrderMetis);

    // Try plain Metis
    solve_and_report("Metis", a, &b, OrderingStrategy::Metis);

    // Try AMD
    solve_and_report("Amd", a, &b, OrderingStrategy::Amd);
}

fn solve_and_report(
    ordering_name: &str,
    a: &faer::sparse::SparseColMat<usize, f64>,
    b: &Col<f64>,
    ordering: OrderingStrategy,
) {
    let opts = AnalyzeOptions { ordering };
    let mut solver = match SparseLDLT::analyze_with_matrix(a, &opts) {
        Ok(s) => s,
        Err(e) => {
            println!("  {:<20} analyze error: {}", ordering_name, e);
            return;
        }
    };

    // Try default block size and single-level
    for &nb in &[256usize, usize::MAX] {
        let factor_opts = FactorOptions {
            outer_block_size: nb,
            ..FactorOptions::default()
        };
        if let Err(e) = solver.factor(a, &factor_opts) {
            println!(
                "  {:<20} nb={:<6} factor error: {}",
                ordering_name,
                if nb == usize::MAX {
                    "MAX".to_string()
                } else {
                    nb.to_string()
                },
                e
            );
            continue;
        }

        let scratch = solver.solve_scratch(1);
        let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
        let stack = faer::dyn_stack::MemStack::new(&mut mem);
        match solver.solve(b, stack) {
            Ok(x) => {
                let be = sparse_backward_error(a, &x, b);
                let stats = solver.stats().unwrap();
                print_stats(ordering_name, nb, be, &stats);
            }
            Err(e) => {
                println!(
                    "  {:<20} nb={:<6} solve error: {}",
                    ordering_name,
                    if nb == usize::MAX {
                        "MAX".to_string()
                    } else {
                        nb.to_string()
                    },
                    e
                );
            }
        }
    }
}

fn print_stats(ordering_name: &str, nb: usize, be: f64, stats: &FactorizationStats) {
    let nb_str = if nb == usize::MAX {
        "MAX".to_string()
    } else {
        nb.to_string()
    };
    println!(
        "  {:<20} nb={:<6} BE={:.2e}  1x1={:<6} 2x2={:<6} delay={:<6} zero={:<4} max_front={}",
        ordering_name,
        nb_str,
        be,
        stats.total_1x1_pivots,
        stats.total_2x2_pivots,
        stats.total_delayed,
        stats.zero_pivots,
        stats.max_front_size,
    );
}

fn main() {
    let all = registry::load_registry().expect("registry");

    // Failing matrices from the full SuiteSparse run, sorted by n
    let targets = [
        "GHS_indef/linverse",  // n=11999
        "Newman/astro-ph",     // n=16706
        "GHS_indef/helm3d01",  // n=32226
        "GHS_indef/dawson5",   // n=51537
        "GHS_indef/copter2",   // n=55476
    ];

    for target in &targets {
        let meta = all.iter().find(|m| m.name == *target);
        let meta = match meta {
            Some(m) => m,
            None => {
                println!("{}: not found in registry", target);
                continue;
            }
        };

        let test = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(t)) => t,
            Ok(None) => {
                println!("{}: not found on disk", target);
                continue;
            }
            Err(e) => {
                println!("{}: load error: {}", target, e);
                continue;
            }
        };

        diagnose(target, &test.matrix);
    }
}
