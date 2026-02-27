//! Solving for multiple right-hand sides with workspace reuse.
//!
//! Demonstrates factoring a matrix once and solving for several different
//! right-hand sides, reusing both the factorization and the solve workspace.
//! This is the typical pattern for problems with multiple load cases or
//! source terms.
//!
//! ```sh
//! cargo run --example multiple_rhs
//! ```

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::{SparseColMat, Triplet};
use faer::{Col, Par};

use rivrs_sparse::symmetric::{AnalyzeOptions, FactorOptions, OrderingStrategy, SparseLDLT};
use rivrs_sparse::validate::sparse_backward_error;

fn main() {
    // -----------------------------------------------------------------------
    // 1. Build a 5×5 symmetric indefinite matrix (full storage).
    //
    //    A = [ 4  -1   .   .   2 ]
    //        [-1   3   1   .   . ]
    //        [ .   1  -2   1   . ]
    //        [ .   .   1   5  -1 ]
    //        [ 2   .   .  -1   3 ]
    // -----------------------------------------------------------------------
    let n = 5;
    #[rustfmt::skip]
    let triplets = vec![
        Triplet::new(0, 0,  4.0),
        Triplet::new(1, 1,  3.0),
        Triplet::new(2, 2, -2.0),
        Triplet::new(3, 3,  5.0),
        Triplet::new(4, 4,  3.0),
        Triplet::new(0, 1, -1.0), Triplet::new(1, 0, -1.0),
        Triplet::new(1, 2,  1.0), Triplet::new(2, 1,  1.0),
        Triplet::new(2, 3,  1.0), Triplet::new(3, 2,  1.0),
        Triplet::new(3, 4, -1.0), Triplet::new(4, 3, -1.0),
        Triplet::new(0, 4,  2.0), Triplet::new(4, 0,  2.0),
    ];
    let a =
        SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets).expect("valid triplets");

    // -----------------------------------------------------------------------
    // 2. Analyze and factor once.
    // -----------------------------------------------------------------------
    let analyze_opts = AnalyzeOptions {
        ordering: OrderingStrategy::Amd,
    };
    let mut solver =
        SparseLDLT::analyze_with_matrix(&a, &analyze_opts).expect("symbolic analysis");
    solver
        .factor(&a, &FactorOptions::default())
        .expect("factorization");

    // -----------------------------------------------------------------------
    // 3. Allocate solve workspace once, reuse for all RHS vectors.
    // -----------------------------------------------------------------------
    let scratch = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);

    // Four different right-hand sides.
    let rhs_vectors: Vec<Col<f64>> = vec![
        Col::from_fn(n, |i| [1.0, 0.0, 0.0, 0.0, 0.0][i]),
        Col::from_fn(n, |i| [0.0, 1.0, 0.0, 0.0, 0.0][i]),
        Col::from_fn(n, |i| [1.0, 1.0, 1.0, 1.0, 1.0][i]),
        Col::from_fn(n, |i| [2.0, -1.0, 3.0, 0.5, -2.0][i]),
    ];

    println!("Solving Ax = b for {} different right-hand sides\n", rhs_vectors.len());

    for (k, b) in rhs_vectors.iter().enumerate() {
        // Reuse the same MemStack for each solve.
        let stack = MemStack::new(&mut mem);
        let x = solver.solve(b, stack, Par::Seq).expect("solve");

        let be = sparse_backward_error(&a, &x, b);
        print!("RHS {}: x = [", k + 1);
        for i in 0..n {
            if i > 0 {
                print!(", ");
            }
            print!("{:+.6}", x[i]);
        }
        println!("]  backward error = {be:.2e}");

        assert!(be < 1e-12, "backward error too large: {be:.2e}");
    }

    println!("\nAll solves passed (backward error < 1e-12).");
}
