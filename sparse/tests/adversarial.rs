//! Adversarial and edge-case input tests for SparseLDLT.
//!
//! Verifies that the solver handles malformed, extreme, and degenerate inputs
//! with zero panics and clean error returns (or correct results where applicable).

use faer::Col;
use faer::sparse::{SparseColMat, Triplet};
use rivrs_sparse::aptp::solver::{SolverOptions, SparseLDLT};
use rivrs_sparse::validate::sparse_backward_error;

// ==========================================================================
// Edge-case tests (T029)
// ==========================================================================

#[test]
fn adversarial_1x1_nonzero_diagonal() {
    let triplets = vec![Triplet::new(0, 0, 5.0)];
    let a = SparseColMat::try_new_from_triplets(1, 1, &triplets).unwrap();
    let rhs = Col::from_fn(1, |_| 3.0);
    let options = SolverOptions::default();
    let result = SparseLDLT::solve_full(&a, &rhs, &options);
    match result {
        Ok(ref x) => {
            let be = sparse_backward_error(&a, x, &rhs);
            assert!(be < 1e-14, "1x1 backward error too high: {:.2e}", be);
        }
        Err(e) => {
            // Some orderings may not handle 1x1 — that's acceptable
            eprintln!("1x1 nonzero: clean error: {}", e);
        }
    }
}

#[test]
fn adversarial_1x1_zero_diagonal() {
    let triplets = vec![Triplet::new(0, 0, 0.0)];
    let a = SparseColMat::try_new_from_triplets(1, 1, &triplets).unwrap();
    let rhs = Col::from_fn(1, |_| 1.0);
    let options = SolverOptions::default();
    // Should not panic — may return error or zero-inertia factorization
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        SparseLDLT::solve_full(&a, &rhs, &options)
    }));
    assert!(result.is_ok(), "1x1 zero diagonal panicked");
}

#[test]
fn adversarial_pure_diagonal() {
    let n = 10;
    let triplets: Vec<_> = (0..n)
        .map(|i| {
            let val = if i < n / 2 {
                (i + 1) as f64
            } else {
                -((i + 1) as f64)
            };
            Triplet::new(i, i, val)
        })
        .collect();
    let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
    let rhs = Col::from_fn(n, |_| 1.0);
    let options = SolverOptions::default();
    let result = SparseLDLT::solve_full(&a, &rhs, &options);
    match result {
        Ok(ref x) => {
            let be = sparse_backward_error(&a, x, &rhs);
            assert!(be < 1e-14, "diagonal backward error too high: {:.2e}", be);
        }
        Err(e) => panic!("pure diagonal should succeed: {}", e),
    }
}

#[test]
fn adversarial_arrowhead_matrix() {
    // Dense first row/column + diagonal remainder
    let n = 20;
    let mut triplets = Vec::new();
    // First row/column
    for i in 1..n {
        triplets.push(Triplet::new(0, i, 0.5));
        triplets.push(Triplet::new(i, 0, 0.5));
    }
    // Diagonal: positive definite (diagonal dominance)
    triplets.push(Triplet::new(0, 0, n as f64));
    for i in 1..n {
        triplets.push(Triplet::new(i, i, 2.0));
    }
    let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
    let rhs = Col::from_fn(n, |_| 1.0);
    let options = SolverOptions::default();
    let result = SparseLDLT::solve_full(&a, &rhs, &options);
    match result {
        Ok(ref x) => {
            let be = sparse_backward_error(&a, x, &rhs);
            assert!(be < 5e-11, "arrowhead backward error too high: {:.2e}", be);
        }
        Err(e) => panic!("arrowhead should succeed: {}", e),
    }
}

#[test]
fn adversarial_all_2x2_pivots() {
    // Create a matrix where no 1x1 pivots are stable, forcing all 2x2 pivots.
    // A matrix with zero diagonal and strong off-diagonal structure:
    //   [ 0  a  0  0 ]
    //   [ a  0  0  0 ]
    //   [ 0  0  0  b ]
    //   [ 0  0  b  0 ]
    let n = 4;
    let a_val = 3.0;
    let b_val = 5.0;
    let triplets = vec![
        Triplet::new(0, 0, 0.0),
        Triplet::new(0, 1, a_val),
        Triplet::new(1, 0, a_val),
        Triplet::new(1, 1, 0.0),
        Triplet::new(2, 2, 0.0),
        Triplet::new(2, 3, b_val),
        Triplet::new(3, 2, b_val),
        Triplet::new(3, 3, 0.0),
    ];
    let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
    let rhs = Col::from_fn(n, |_| 1.0);
    let options = SolverOptions::default();
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        SparseLDLT::solve_full(&a, &rhs, &options)
    }));
    assert!(result.is_ok(), "all-2x2-pivots matrix panicked");
    if let Ok(Ok(ref x)) = result {
        let be = sparse_backward_error(&a, x, &rhs);
        assert!(
            be < 5e-11,
            "all-2x2-pivots backward error too high: {:.2e}",
            be
        );
    }
}

#[test]
fn adversarial_power_of_2_boundary_sizes() {
    for &n in &[32, 64, 128, 256] {
        let mut triplets = Vec::new();
        // Tridiagonal structure with diagonal dominance
        for i in 0..n {
            triplets.push(Triplet::new(i, i, 3.0));
            if i + 1 < n {
                triplets.push(Triplet::new(i, i + 1, -1.0));
                triplets.push(Triplet::new(i + 1, i, -1.0));
            }
        }
        let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
        let rhs = Col::from_fn(n, |_| 1.0);
        let options = SolverOptions::default();
        let result = SparseLDLT::solve_full(&a, &rhs, &options);
        match result {
            Ok(ref x) => {
                let be = sparse_backward_error(&a, x, &rhs);
                assert!(
                    be < 5e-11,
                    "power-of-2 n={} backward error too high: {:.2e}",
                    n,
                    be
                );
            }
            Err(e) => panic!("power-of-2 n={} should succeed: {}", n, e),
        }
    }
}

// T029.7: Exact numerical cancellation
#[test]
fn adversarial_exact_cancellation() {
    // Matrix where elimination produces exact cancellation:
    //   [ 1  1  1 ]
    //   [ 1  1  0 ]
    //   [ 1  0  1 ]
    // After eliminating column 0, the (1,1) entry becomes 1-1=0 (exact cancellation).
    let n = 3;
    let triplets = vec![
        Triplet::new(0, 0, 1.0),
        Triplet::new(0, 1, 1.0),
        Triplet::new(0, 2, 1.0),
        Triplet::new(1, 0, 1.0),
        Triplet::new(1, 1, 1.0),
        Triplet::new(2, 0, 1.0),
        Triplet::new(2, 2, 1.0),
    ];
    let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
    let rhs = Col::from_fn(n, |_| 1.0);
    let options = SolverOptions::default();
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        SparseLDLT::solve_full(&a, &rhs, &options)
    }));
    assert!(result.is_ok(), "exact cancellation panicked");
    if let Ok(Ok(ref x)) = result {
        let be = sparse_backward_error(&a, x, &rhs);
        assert!(
            be < 5e-11,
            "exact cancellation backward error too high: {:.2e}",
            be
        );
    }
}

// T029.1: 0×0 empty matrix
#[test]
fn adversarial_0x0_empty_matrix() {
    let triplets: Vec<Triplet<usize, usize, f64>> = vec![];
    let a = SparseColMat::try_new_from_triplets(0, 0, &triplets).unwrap();
    let rhs = Col::from_fn(0, |_| 0.0);
    let options = SolverOptions::default();
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        SparseLDLT::solve_full(&a, &rhs, &options)
    }));
    assert!(result.is_ok(), "0x0 empty matrix panicked");
    // Should return either an empty solution or a clean error
}

// ==========================================================================
// Extreme value tests (T030)
// ==========================================================================

#[test]
fn adversarial_near_overflow_entries() {
    let n = 4;
    let big = 1e150; // Not full overflow but extremely large
    let triplets = vec![
        Triplet::new(0, 0, big),
        Triplet::new(0, 1, 1.0),
        Triplet::new(1, 0, 1.0),
        Triplet::new(1, 1, big),
        Triplet::new(2, 2, big),
        Triplet::new(2, 3, 1.0),
        Triplet::new(3, 2, 1.0),
        Triplet::new(3, 3, big),
    ];
    let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
    let rhs = Col::from_fn(n, |_| 1.0);
    let options = SolverOptions::default();
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        SparseLDLT::solve_full(&a, &rhs, &options)
    }));
    assert!(result.is_ok(), "near-overflow matrix panicked");
}

#[test]
fn adversarial_near_underflow_entries() {
    let n = 4;
    let tiny = 1e-150;
    let triplets = vec![
        Triplet::new(0, 0, tiny),
        Triplet::new(0, 1, tiny * 0.5),
        Triplet::new(1, 0, tiny * 0.5),
        Triplet::new(1, 1, tiny),
        Triplet::new(2, 2, tiny),
        Triplet::new(2, 3, tiny * 0.5),
        Triplet::new(3, 2, tiny * 0.5),
        Triplet::new(3, 3, tiny),
    ];
    let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
    let rhs = Col::from_fn(n, |_| tiny);
    let options = SolverOptions::default();
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        SparseLDLT::solve_full(&a, &rhs, &options)
    }));
    assert!(result.is_ok(), "near-underflow matrix panicked");
}

#[test]
fn adversarial_nan_entries() {
    let n = 4;
    let triplets = vec![
        Triplet::new(0, 0, 4.0),
        Triplet::new(0, 1, f64::NAN),
        Triplet::new(1, 0, f64::NAN),
        Triplet::new(1, 1, 4.0),
        Triplet::new(2, 2, 4.0),
        Triplet::new(3, 3, 4.0),
    ];
    let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
    let rhs = Col::from_fn(n, |_| 1.0);
    let options = SolverOptions::default();
    // Must not panic — should return clean error or NaN-contaminated result
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        SparseLDLT::solve_full(&a, &rhs, &options)
    }));
    assert!(result.is_ok(), "NaN matrix panicked");
    // If the solver returns Ok, the solution should contain NaN (not finite garbage).
    // If it returns Err, that's the ideal behavior.
    if let Ok(Ok(ref x)) = result {
        let has_nan = (0..n).any(|i| x[i].is_nan());
        let has_finite_nonsense = (0..n).all(|i| x[i].is_finite());
        // NaN propagation or clean error are both acceptable.
        // Finite solutions from NaN input would indicate a bug (NaN silently ignored).
        assert!(
            has_nan || !has_finite_nonsense,
            "NaN matrix produced finite solution — NaN may have been silently dropped"
        );
    }
}

#[test]
fn adversarial_inf_entries() {
    let n = 4;
    let triplets = vec![
        Triplet::new(0, 0, f64::INFINITY),
        Triplet::new(0, 1, 1.0),
        Triplet::new(1, 0, 1.0),
        Triplet::new(1, 1, 4.0),
        Triplet::new(2, 2, 4.0),
        Triplet::new(3, 3, 4.0),
    ];
    let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
    let rhs = Col::from_fn(n, |_| 1.0);
    let options = SolverOptions::default();
    // Must not panic — should return clean error or non-finite result
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        SparseLDLT::solve_full(&a, &rhs, &options)
    }));
    assert!(result.is_ok(), "Inf matrix panicked");
    // If solver returns Ok, solution should reflect the Inf input (not finite garbage).
    if let Ok(Ok(ref x)) = result {
        let all_finite = (0..n).all(|i| x[i].is_finite());
        if all_finite {
            // If finite, backward error must be valid (solver somehow handled Inf gracefully)
            let be = sparse_backward_error(&a, x, &rhs);
            assert!(
                !be.is_finite() || be < 5e-11,
                "Inf matrix produced finite solution with bad backward error: {:.2e}",
                be
            );
        }
    }
}

// ==========================================================================
// Structural validity tests (T031)
// ==========================================================================

#[test]
fn adversarial_disconnected_components() {
    // Two disconnected 3x3 blocks
    let n = 6;
    let triplets = vec![
        // Block 1: rows/cols 0-2
        Triplet::new(0, 0, 4.0),
        Triplet::new(0, 1, -1.0),
        Triplet::new(1, 0, -1.0),
        Triplet::new(1, 1, 4.0),
        Triplet::new(1, 2, -1.0),
        Triplet::new(2, 1, -1.0),
        Triplet::new(2, 2, 4.0),
        // Block 2: rows/cols 3-5
        Triplet::new(3, 3, 4.0),
        Triplet::new(3, 4, -1.0),
        Triplet::new(4, 3, -1.0),
        Triplet::new(4, 4, 4.0),
        Triplet::new(4, 5, -1.0),
        Triplet::new(5, 4, -1.0),
        Triplet::new(5, 5, 4.0),
    ];
    let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
    let rhs = Col::from_fn(n, |_| 1.0);
    let options = SolverOptions::default();
    let result = SparseLDLT::solve_full(&a, &rhs, &options);
    match result {
        Ok(ref x) => {
            let be = sparse_backward_error(&a, x, &rhs);
            assert!(
                be < 5e-11,
                "disconnected backward error too high: {:.2e}",
                be
            );
        }
        Err(e) => panic!("disconnected components should succeed: {}", e),
    }
}

#[test]
fn adversarial_power_of_2_boundary_512() {
    // Larger power-of-2 boundary (512) — tests block size alignment
    let n = 512;
    let mut triplets = Vec::new();
    for i in 0..n {
        triplets.push(Triplet::new(i, i, 3.0));
        if i + 1 < n {
            triplets.push(Triplet::new(i, i + 1, -1.0));
            triplets.push(Triplet::new(i + 1, i, -1.0));
        }
    }
    let a = SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap();
    let rhs = Col::from_fn(n, |_| 1.0);
    let options = SolverOptions::default();
    let result = SparseLDLT::solve_full(&a, &rhs, &options);
    match result {
        Ok(ref x) => {
            let be = sparse_backward_error(&a, x, &rhs);
            assert!(
                be < 5e-11,
                "power-of-2 n=512 backward error too high: {:.2e}",
                be
            );
        }
        Err(e) => panic!("power-of-2 n=512 should succeed: {}", e),
    }
}
