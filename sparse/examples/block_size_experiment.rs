//! Compare backward error with different block sizes on specific matrices.
//!
//! Tests whether the BLAS-3 two-level factorization path introduces precision loss
//! compared to the single-level (factor_inner only) path.
//!
//! Usage:
//!   cargo run --release --example block_size_experiment -- <matrix.mtx>

use faer::Col;
use faer::Par;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::SparseColMat;
use std::path::Path;

use rivrs_sparse::aptp::{AnalyzeOptions, FactorOptions, OrderingStrategy, SparseLDLT};
use rivrs_sparse::io::mtx::load_mtx;

fn backward_error(matrix: &SparseColMat<usize, f64>, x: &Col<f64>, b: &Col<f64>) -> f64 {
    let n = matrix.nrows();
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();

    let mut ax = vec![0.0f64; n];
    for j in 0..n {
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for idx in start..end {
            let i = row_indices[idx];
            ax[i] += values[idx] * x[j];
        }
    }
    let norm_r: f64 = ax
        .iter()
        .zip(b.as_ref().iter())
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0, f64::max);
    let norm_a: f64 = values.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let norm_x: f64 = (0..n).map(|i| x[i].abs()).fold(0.0, f64::max);
    let norm_b: f64 = b.as_ref().iter().map(|v| v.abs()).fold(0.0, f64::max);
    let denom = norm_a * norm_x + norm_b;
    if denom > 0.0 { norm_r / denom } else { 0.0 }
}

fn solve_with_options(
    matrix: &SparseColMat<usize, f64>,
    rhs: &Col<f64>,
    outer_bs: usize,
    inner_bs: usize,
) -> Result<f64, Box<dyn std::error::Error>> {
    let _n = matrix.nrows();
    let analyze_opts = AnalyzeOptions {
        ordering: OrderingStrategy::MatchOrderMetis,
    };
    let factor_opts = FactorOptions {
        outer_block_size: outer_bs,
        inner_block_size: inner_bs,
        ..FactorOptions::default()
    };

    let mut solver = SparseLDLT::analyze_with_matrix(matrix, &analyze_opts)?;
    solver.factor(matrix, &factor_opts)?;

    let scratch = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);
    let stack = MemStack::new(&mut mem);
    let x = solver.solve(rhs, stack, Par::Seq)?;

    Ok(backward_error(matrix, &x, rhs))
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        eprintln!("Usage: block_size_experiment <matrix.mtx>");
        std::process::exit(1);
    }

    let path = Path::new(&args[0]);
    let matrix = load_mtx(path).expect("Failed to load matrix");
    let n = matrix.nrows();

    // Generate RHS: A * ones(n)
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();
    let mut rhs_vec = vec![0.0f64; n];
    for j in 0..n {
        for idx in col_ptrs[j]..col_ptrs[j + 1] {
            rhs_vec[row_indices[idx]] += values[idx];
        }
    }
    let rhs = Col::from_fn(n, |i| rhs_vec[i]);

    let name = path.file_stem().unwrap().to_str().unwrap();
    eprintln!("Matrix: {} (n={})", name, n);
    eprintln!();

    // Test configurations: (outer_block_size, inner_block_size, description)
    let configs: Vec<(usize, usize, &str)> = vec![
        (256, 32, "default (ob=256, ib=32)"),
        (256, 256, "ib=ob (ob=256, ib=256)"),
        (100000, 100000, "unblocked (ob=huge, ib=huge)"),
    ];

    eprintln!("  {:<45} {:>12}", "Configuration", "Backward Error");
    eprintln!("  {:-<60}", "");

    for (obs, ibs, desc) in &configs {
        match solve_with_options(&matrix, &rhs, *obs, *ibs) {
            Ok(be) => {
                eprintln!("  {:<45} {:>12.2e}", desc, be);
            }
            Err(e) => {
                eprintln!("  {:<45} FAILED: {}", desc, e);
            }
        }
    }
}
