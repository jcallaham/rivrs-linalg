//! Basic usage of the rivrs-sparse SSIDS solver.
//!
//! Demonstrates the three-phase API (analyze → factor → solve) on a small
//! symmetric indefinite matrix constructed from triplets. No external data
//! files or feature flags are required.
//!
//! ```sh
//! cargo run --example basic_usage
//! ```

use faer::sparse::{SparseColMat, Triplet};
use faer::{Col, Par};

use rivrs_sparse::symmetric::{AnalyzeOptions, FactorOptions, OrderingStrategy, SparseLDLT};
use rivrs_sparse::validate::sparse_backward_error;

fn main() {
    // -----------------------------------------------------------------------
    // 1. Build a 6×6 symmetric indefinite matrix from triplets.
    //
    //    A = [ 2  -1   .   .   3   . ]
    //        [-1  -4   1   .   .   . ]
    //        [ .   1   3  -2   .   . ]
    //        [ .   .  -2   1   .   4 ]
    //        [ 3   .   .   .  -5   2 ]
    //        [ .   .   .   4   2   6 ]
    //
    //    This matrix has both positive and negative eigenvalues (indefinite).
    //    Full symmetric storage: both upper and lower triangles are stored.
    // -----------------------------------------------------------------------
    let n = 6;
    #[rustfmt::skip]
    let triplets = vec![
        // Diagonal
        Triplet::new(0, 0,  2.0),
        Triplet::new(1, 1, -4.0),
        Triplet::new(2, 2,  3.0),
        Triplet::new(3, 3,  1.0),
        Triplet::new(4, 4, -5.0),
        Triplet::new(5, 5,  6.0),
        // Off-diagonal (both triangles for full symmetric storage)
        Triplet::new(0, 1, -1.0), Triplet::new(1, 0, -1.0),
        Triplet::new(1, 2,  1.0), Triplet::new(2, 1,  1.0),
        Triplet::new(2, 3, -2.0), Triplet::new(3, 2, -2.0),
        Triplet::new(0, 4,  3.0), Triplet::new(4, 0,  3.0),
        Triplet::new(3, 5,  4.0), Triplet::new(5, 3,  4.0),
        Triplet::new(4, 5,  2.0), Triplet::new(5, 4,  2.0),
    ];

    let a =
        SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets).expect("valid triplets");

    // -----------------------------------------------------------------------
    // 2. Create a right-hand side b = A * x_true for a known solution.
    // -----------------------------------------------------------------------
    let x_true = vec![1.0, -2.0, 3.0, -1.0, 0.5, 2.0];
    let b_vec = sparse_matvec(&a, &x_true);
    let b = Col::from_fn(n, |i| b_vec[i]);

    // -----------------------------------------------------------------------
    // 3. Analyze: symbolic factorization (reusable for same sparsity pattern).
    // -----------------------------------------------------------------------
    let analyze_opts = AnalyzeOptions {
        ordering: OrderingStrategy::Amd, // AMD is fine for small matrices
    };
    let mut solver = SparseLDLT::analyze_with_matrix(&a, &analyze_opts).expect("symbolic analysis");

    // -----------------------------------------------------------------------
    // 4. Factor: numeric LDL^T factorization with APTP pivoting.
    // -----------------------------------------------------------------------
    let factor_opts = FactorOptions::default();
    solver.factor(&a, &factor_opts).expect("factorization");

    // Print factorization statistics.
    if let Some(stats) = solver.stats() {
        println!("Factorization statistics:");
        println!("  1×1 pivots:     {}", stats.total_1x1_pivots);
        println!("  2×2 pivots:     {}", stats.total_2x2_pivots);
        println!("  Delayed pivots: {}", stats.total_delayed);
        println!("  Max front size: {}", stats.max_front_size);
    }

    // Print inertia (eigenvalue sign counts).
    if let Some(inertia) = solver.inertia() {
        println!(
            "  Inertia: ({} positive, {} negative, {} zero)",
            inertia.positive, inertia.negative, inertia.zero
        );
    }

    // -----------------------------------------------------------------------
    // 5. Solve: forward/diagonal/backward substitution.
    // -----------------------------------------------------------------------
    let scratch = solver.solve_scratch(1);
    let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
    let stack = faer::dyn_stack::MemStack::new(&mut mem);
    let x = solver.solve(&b, stack, Par::Seq).expect("solve");

    // -----------------------------------------------------------------------
    // 6. Validate: check backward error ||Ax - b|| / (||A|| ||x|| + ||b||).
    // -----------------------------------------------------------------------
    let be = sparse_backward_error(&a, &x, &b);
    println!("\nSolution:");
    for i in 0..n {
        println!("  x[{i}] = {:+.12}  (true: {:+.1})", x[i], x_true[i]);
    }
    println!("\nBackward error: {be:.2e}");

    if be < 1e-12 {
        println!("PASS (backward error < 1e-12)");
    } else {
        println!("WARN (backward error = {be:.2e})");
    }
}

/// Sparse matrix-vector product b = A * x using CSC storage.
fn sparse_matvec(a: &SparseColMat<usize, f64>, x: &[f64]) -> Vec<f64> {
    let n = a.nrows();
    let sym = a.symbolic();
    let cp = sym.col_ptr();
    let ri = sym.row_idx();
    let vals = a.val();
    let mut b = vec![0.0; n];
    for j in 0..n {
        for idx in cp[j]..cp[j + 1] {
            b[ri[idx]] += vals[idx] * x[j];
        }
    }
    b
}
