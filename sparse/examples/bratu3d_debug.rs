//! Compare single-level vs two-level APTP on SuiteSparse matrices.
//!
//! Uses the same solver configuration as `test_solve_suitesparse_full`:
//! - Category-dependent ordering (hard-indefinite → MatchOrderMetis, else → Metis)
//! - Same RHS construction: x_exact[i] = ((i % 7) - 3) / 3
//! - Same backward error thresholds (SPRAL 5e-11, relaxed 1e-8)
//!
//! Usage: cargo run --example bratu3d_debug --release

use faer::Col;
use rivrs_sparse::aptp::{
    AnalyzeOptions, FactorOptions, OrderingStrategy, SparseLDLT,
};
use rivrs_sparse::io::registry::{self, MatrixMetadata};
use rivrs_sparse::validate::sparse_backward_error;

const SPRAL_BE_THRESHOLD: f64 = 5e-11;
const RELAXED_BE_THRESHOLD: f64 = 1e-8;

/// Compute b = A*x for a full symmetric CSC matrix (both triangles stored).
fn sparse_matvec(a: &faer::sparse::SparseColMat<usize, f64>, x: &[f64]) -> Vec<f64> {
    let n = a.nrows();
    let sym = a.symbolic();
    let col_ptrs = sym.col_ptr();
    let row_indices = sym.row_idx();
    let values = a.val();

    let mut result = vec![0.0f64; n];
    for j in 0..n {
        for idx in col_ptrs[j]..col_ptrs[j + 1] {
            result[row_indices[idx]] += values[idx] * x[j];
        }
    }
    result
}

/// Select ordering strategy based on matrix category, matching test_solve_suitesparse_full.
fn ordering_for_category(category: &str) -> OrderingStrategy {
    if category == "hard-indefinite" {
        OrderingStrategy::MatchOrderMetis
    } else {
        OrderingStrategy::Metis
    }
}

fn solve_with_block_size(
    name: &str,
    a: &faer::sparse::SparseColMat<usize, f64>,
    b: &Col<f64>,
    ordering: &OrderingStrategy,
    outer_block_size: usize,
) {
    let analyze_opts = AnalyzeOptions {
        ordering: ordering.clone(),
    };
    let mut solver = SparseLDLT::analyze_with_matrix(a, &analyze_opts).expect("analyze");

    let factor_opts = FactorOptions {
        outer_block_size,
        ..FactorOptions::default()
    };
    solver.factor(a, &factor_opts).expect("factor");

    let scratch = solver.solve_scratch(1);
    let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
    let stack = faer::dyn_stack::MemStack::new(&mut mem);
    let x = solver.solve(b, stack).expect("solve");

    let be = sparse_backward_error(a, &x, b);
    let stats = solver.stats().unwrap();

    let status = if be < SPRAL_BE_THRESHOLD {
        "PASS"
    } else if be < RELAXED_BE_THRESHOLD {
        "RELAXED"
    } else {
        "FAIL"
    };

    println!(
        "  {:<20} nb={:<6} BE={:.2e} ({:<7})  1x1={:<6} 2x2={:<6} delay={:<6} max_front={}",
        name,
        if outer_block_size == usize::MAX {
            "MAX".to_string()
        } else {
            outer_block_size.to_string()
        },
        be,
        status,
        stats.total_1x1_pivots,
        stats.total_2x2_pivots,
        stats.total_delayed,
        stats.max_front_size,
    );
}

fn run_matrix(meta: &MatrixMetadata) {
    let test = match registry::load_test_matrix_from_entry(meta) {
        Ok(Some(t)) => t,
        Ok(None) => {
            println!("{}: MISSING", meta.name);
            return;
        }
        Err(e) => {
            println!("{}: LOAD ERROR: {}", meta.name, e);
            return;
        }
    };
    let a = &test.matrix;
    let n = a.nrows();

    // Same RHS construction as test_solve_suitesparse_full
    let x_vec: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
    let b_vec = sparse_matvec(a, &x_vec);
    let b = Col::from_fn(n, |i| b_vec[i]);

    // Same category-dependent ordering as test_solve_suitesparse_full
    let ordering = ordering_for_category(&meta.category);

    println!(
        "{}  (n={}, nnz={}, category={}, ordering={:?}):",
        meta.name, n, meta.nnz, meta.category, ordering
    );

    // Single-level (force via large outer_block_size)
    solve_with_block_size("single-level", a, &b, &ordering, usize::MAX);

    // Two-level with default block size (matches solve_full via FactorOptions::default)
    solve_with_block_size("two-level (default)", a, &b, &ordering, 256);

    // Two-level with other block sizes for comparison
    for &nb in &[128, 512, 1024] {
        solve_with_block_size("two-level", a, &b, &ordering, nb);
    }

    println!();
}

fn main() {
    let all = registry::load_registry().expect("registry");

    let targets = [
        "GHS_indef/bratu3d",
        "GHS_indef/stokes128",
        "GHS_indef/bloweybq",
        "GHS_indef/sparsine",
    ];

    for target in &targets {
        let meta = match all.iter().find(|m| m.name == *target) {
            Some(m) => m,
            None => {
                println!("{}: not found in registry", target);
                continue;
            }
        };
        run_matrix(meta);
    }
}
