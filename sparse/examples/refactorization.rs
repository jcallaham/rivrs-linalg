//! Refactorization: same sparsity pattern, different numeric values.
//!
//! Demonstrates analyzing a matrix once, then calling `factor()` repeatedly
//! as the matrix values change. This is the typical pattern for interior-point
//! methods, continuation/homotopy solvers, and time-stepping schemes where
//! the sparsity structure is fixed but the entries evolve.
//!
//! ```sh
//! cargo run --example refactorization
//! ```

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::{SparseColMat, Triplet};
use faer::{Col, Par};

use rivrs_sparse::symmetric::{AnalyzeOptions, FactorOptions, OrderingStrategy, SparseLDLT};
use rivrs_sparse::validate::sparse_backward_error;

/// Build a parametric 4×4 symmetric indefinite matrix A(t).
///
///    A(t) = [ 2+t   -1     .     1  ]
///           [ -1    3-t     2     .  ]
///           [  .     2    -1+t   -1  ]
///           [  1     .    -1    4-2t  ]
///
/// The sparsity pattern is constant; only the values depend on `t`.
fn build_matrix(t: f64) -> SparseColMat<usize, f64> {
    #[rustfmt::skip]
    let triplets = vec![
        Triplet::new(0, 0, 2.0 + t),
        Triplet::new(1, 1, 3.0 - t),
        Triplet::new(2, 2, -1.0 + t),
        Triplet::new(3, 3, 4.0 - 2.0 * t),
        Triplet::new(0, 1, -1.0), Triplet::new(1, 0, -1.0),
        Triplet::new(1, 2,  2.0), Triplet::new(2, 1,  2.0),
        Triplet::new(2, 3, -1.0), Triplet::new(3, 2, -1.0),
        Triplet::new(0, 3,  1.0), Triplet::new(3, 0,  1.0),
    ];
    SparseColMat::<usize, f64>::try_new_from_triplets(4, 4, &triplets).expect("valid triplets")
}

fn main() {
    let n = 4;
    let b = Col::from_fn(n, |i| [1.0, -1.0, 2.0, 0.5][i]);

    // -----------------------------------------------------------------------
    // 1. Analyze once using the sparsity pattern at t = 0.
    // -----------------------------------------------------------------------
    let a0 = build_matrix(0.0);
    let analyze_opts = AnalyzeOptions {
        ordering: OrderingStrategy::Amd,
    };
    let mut solver =
        SparseLDLT::analyze_with_matrix(&a0, &analyze_opts).expect("symbolic analysis");

    println!("Symbolic analysis done (reused for all parameter values)\n");
    println!("{:<8} {:>40}  {:>12}  {:>10}", "t", "solution", "bwd_err", "inertia");
    println!("{}", "-".repeat(78));

    // -----------------------------------------------------------------------
    // 2. Sweep parameter t, refactoring each time.
    // -----------------------------------------------------------------------
    let factor_opts = FactorOptions::default();
    let mut mem: Option<MemBuffer> = None;

    for &t in &[0.0, 0.5, 1.0, 1.5, 2.0] {
        let a = build_matrix(t);

        // Refactor with the new values (reuses symbolic analysis).
        solver.factor(&a, &factor_opts).expect("factorization");

        // Allocate workspace on first iteration, reuse thereafter.
        let scratch = solver.solve_scratch(1);
        let buf = mem.get_or_insert_with(|| MemBuffer::new(scratch));
        let stack = MemStack::new(buf);

        let x = solver.solve(&b, stack, Par::Seq).expect("solve");
        let be = sparse_backward_error(&a, &x, &b);

        let inertia = solver.inertia().expect("inertia available");

        let x_str = format!(
            "[{:+.4}, {:+.4}, {:+.4}, {:+.4}]",
            x[0], x[1], x[2], x[3]
        );
        println!(
            "{:<8.1} {:>40}  {:>12.2e}  ({}, {}, {})",
            t, x_str, be, inertia.positive, inertia.negative, inertia.zero
        );

        assert!(be < 1e-12, "backward error too large at t={t}: {be:.2e}");
    }

    println!("\nAll refactorizations passed (backward error < 1e-12).");
}
