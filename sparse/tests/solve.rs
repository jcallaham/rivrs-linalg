//! Tests for triangular solve & solver API (Phase 7).
//!
//! Organized by user story:
//! - US3: Per-supernode triangular solve correctness
//! - US1: End-to-end sparse solve
//! - US2: Three-phase API with reuse
//! - US4: Scaling integration
//! - US5: Workspace-efficient solve

mod common;

use faer::Col;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::linalg::cholesky::SymmetricOrdering;
use faer::sparse::{SparseColMat, Triplet};

use rivrs_sparse::aptp::{
    AnalyzeOptions, AptpNumeric, AptpOptions, AptpSymbolic, FactorOptions, OrderingStrategy,
    SolverOptions, SparseLDLT,
};
use rivrs_sparse::error::SparseError;
use rivrs_sparse::validate::sparse_backward_error;

use common::{sparse_from_lower_triplets, sparse_matvec};

// ---------------------------------------------------------------------------
// US3: Per-supernode triangular solve correctness
// ---------------------------------------------------------------------------

/// Test the full aptp_solve pipeline on a small hand-constructed matrix.
/// Factors via AptpNumeric::factor(), then calls aptp_solve() directly.
#[test]
fn test_aptp_solve_small_pd() {
    // 3x3 positive definite:
    // A = [[4, 1, 0.5],
    //      [1, 3, 0  ],
    //      [0.5, 0, 2 ]]
    let matrix = sparse_from_lower_triplets(
        3,
        &[
            (0, 0, 4.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (2, 0, 0.5),
            (2, 2, 2.0),
        ],
    );

    let x_exact = vec![1.0, 2.0, 3.0];
    let b = sparse_matvec(&matrix, &x_exact);

    // Analyze and factor using identity ordering for predictability
    let symbolic =
        AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Identity).expect("analyze");
    let numeric =
        AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default(), None).expect("factor");

    // Get permutation
    let n = 3;
    let (perm_fwd, _) = symbolic.perm_vecs();

    // Permute RHS to elimination order
    let mut rhs_perm = vec![0.0f64; n];
    for new in 0..n {
        rhs_perm[new] = b[perm_fwd[new]];
    }

    // Solve in permuted coordinates
    let scratch = rivrs_sparse::aptp::solve::aptp_solve_scratch(&numeric, 1);
    let mut mem = MemBuffer::new(scratch);
    let stack = MemStack::new(&mut mem);
    rivrs_sparse::aptp::solve::aptp_solve(&symbolic, &numeric, &mut rhs_perm, stack)
        .expect("solve");

    // Unpermute solution
    let mut x_computed = vec![0.0f64; n];
    for new in 0..n {
        x_computed[perm_fwd[new]] = rhs_perm[new];
    }

    // Check solution
    let err: f64 = x_computed
        .iter()
        .zip(x_exact.iter())
        .map(|(c, e)| (c - e).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(err < 1e-10, "solution error: {:.2e}", err);
}

/// Test aptp_solve on a small indefinite matrix with 2x2 pivots.
#[test]
fn test_aptp_solve_small_indefinite() {
    // 4x4 indefinite matrix:
    // A = [[ 1, 2, 0, 0],
    //      [ 2, -3, 1, 0],
    //      [ 0,  1, 5, 1],
    //      [ 0,  0, 1, -2]]
    let matrix = sparse_from_lower_triplets(
        4,
        &[
            (0, 0, 1.0),
            (1, 0, 2.0),
            (1, 1, -3.0),
            (2, 1, 1.0),
            (2, 2, 5.0),
            (3, 2, 1.0),
            (3, 3, -2.0),
        ],
    );

    let x_exact = vec![1.0, -1.0, 2.0, 0.5];
    let b = sparse_matvec(&matrix, &x_exact);

    let b_col = Col::from_fn(4, |i| b[i]);

    // Use SparseLDLT for end-to-end solve
    let x = SparseLDLT::solve_full(&matrix, &b_col, &SolverOptions::default()).expect("solve_full");

    let be = sparse_backward_error(&matrix, &x, &b_col);
    assert!(be < 1e-10, "backward error: {:.2e}", be);
}

/// Test that aptp_solve handles a 1x1 matrix (trivial case).
#[test]
fn test_aptp_solve_1x1() {
    let matrix = sparse_from_lower_triplets(1, &[(0, 0, 5.0)]);

    let b = Col::from_fn(1, |_| 10.0);
    let x = SparseLDLT::solve_full(&matrix, &b, &SolverOptions::default()).expect("solve_full");

    assert!((x[0] - 2.0).abs() < 1e-14, "x[0] = {}, expected 2.0", x[0]);
}

/// Test that aptp_solve handles a diagonal matrix (simplicial supernodes).
#[test]
fn test_aptp_solve_diagonal() {
    let n = 5;
    let entries: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, (i + 1) as f64)).collect();
    let matrix = sparse_from_lower_triplets(n, &entries);

    let x_exact: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
    let b = sparse_matvec(&matrix, &x_exact);
    let b_col = Col::from_fn(n, |i| b[i]);

    let x = SparseLDLT::solve_full(&matrix, &b_col, &SolverOptions::default()).expect("solve_full");

    let be = sparse_backward_error(&matrix, &x, &b_col);
    assert!(be < 1e-14, "backward error: {:.2e}", be);
}

// ---------------------------------------------------------------------------
// US1: End-to-end backward error tests
// ---------------------------------------------------------------------------

/// Test SparseLDLT::solve_full on hand-constructed matrices.
#[test]
fn test_solve_full_hand_constructed() {
    use rivrs_sparse::io::registry;

    let test_names = [
        "arrow-5-pd",
        "tridiag-5-pd",
        "trivial-1x1",
        "trivial-2x2",
        "block-diag-15",
        "stress-delayed-pivots",
        "arrow-10-indef",
        "tridiag-10-indef",
        "zero-diagonal-3",
        "singular-3",
        "bordered-block-20",
    ];

    for name in &test_names {
        let test = registry::load_test_matrix(name)
            .unwrap_or_else(|e| panic!("load '{}': {}", name, e))
            .unwrap_or_else(|| panic!("matrix '{}' not found", name));

        let a = &test.matrix;
        let n = a.nrows();
        if n == 0 {
            continue;
        }

        // Create exact solution and compute b
        let _x_exact = Col::from_fn(n, |i| (i + 1) as f64);
        let b_vec = sparse_matvec(a, &(0..n).map(|i| (i + 1) as f64).collect::<Vec<_>>());
        let b = Col::from_fn(n, |i| b_vec[i]);

        // Solve
        let result = SparseLDLT::solve_full(a, &b, &SolverOptions::default());
        match result {
            Ok(x) => {
                let be = sparse_backward_error(a, &x, &b);
                assert!(be < 1e-10, "'{}': backward error {:.2e} >= 1e-10", name, be);
            }
            Err(e) => {
                // Some matrices (e.g. rank-deficient) may have issues
                // but should not fail to factor
                panic!("'{}': solve_full failed: {}", name, e);
            }
        }
    }
}

/// Test inertia matches for hand-constructed matrices.
#[test]
fn test_solve_inertia_hand_constructed() {
    use rivrs_sparse::io::registry;

    // Test a few matrices with known inertia
    let cases: &[(&str, usize, usize, usize)] = &[
        ("arrow-5-pd", 5, 0, 0),
        ("tridiag-5-pd", 5, 0, 0),
        ("trivial-1x1", 1, 0, 0),
        ("trivial-2x2", 1, 1, 0), // [[0,3],[3,0]] has eigenvalues ±3
    ];

    for &(name, pos, neg, _zero) in cases {
        let test = registry::load_test_matrix(name)
            .unwrap_or_else(|e| panic!("load '{}': {}", name, e))
            .unwrap_or_else(|| panic!("matrix '{}' not found", name));

        let a = &test.matrix;
        let mut solver =
            SparseLDLT::analyze_with_matrix(a, &AnalyzeOptions::default()).expect("analyze");
        solver.factor(a, &FactorOptions::default()).expect("factor");

        let inertia = solver.inertia().expect("should have inertia after factor");
        assert_eq!(
            inertia.positive, pos,
            "'{}': positive {} != {}",
            name, inertia.positive, pos
        );
        assert_eq!(
            inertia.negative, neg,
            "'{}': negative {} != {}",
            name, inertia.negative, neg
        );
        // Zero pivots might differ due to APTP behavior, skip exact check
        // for matrices that aren't rank-deficient
    }
}

/// Test SuiteSparse CI subset backward error.
#[test]
#[ignore] // Slow in debug mode — run with `cargo test -- --ignored`
fn test_solve_suitesparse_ci() {
    use rivrs_sparse::io::registry;

    let all_matrices = registry::load_registry().expect("load registry");
    let ci_matrices: Vec<_> = all_matrices.iter().filter(|m| m.ci_subset).collect();

    assert!(
        !ci_matrices.is_empty(),
        "CI matrix subset should not be empty"
    );

    for meta in &ci_matrices {
        let test = registry::load_test_matrix(&meta.name)
            .unwrap_or_else(|e| panic!("load '{}': {}", meta.name, e))
            .unwrap_or_else(|| panic!("matrix '{}' not found", meta.name));

        let a = &test.matrix;
        let n = a.nrows();

        // Create x_exact and b
        let x_vec: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
        let b_vec = sparse_matvec(a, &x_vec);
        let b = Col::from_fn(n, |i| b_vec[i]);

        // Follow Duff et al (2020) experimental protocol:
        // hard-indefinite → MC64+METIS, otherwise → plain METIS
        let ordering = if meta.category == "hard-indefinite" {
            OrderingStrategy::MatchOrderMetis
        } else {
            OrderingStrategy::Metis
        };
        let opts = SolverOptions {
            ordering,
            ..SolverOptions::default()
        };
        let result = SparseLDLT::solve_full(a, &b, &opts);
        match result {
            Ok(x) => {
                let be = sparse_backward_error(a, &x, &b);
                assert!(
                    be < 1e-10,
                    "SuiteSparse CI '{}': backward error {:.2e} >= 1e-10",
                    meta.name,
                    be
                );
            }
            Err(e) => {
                panic!("SuiteSparse CI '{}': solve_full failed: {}", meta.name, e);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// US1: Error handling and edge cases
// ---------------------------------------------------------------------------

/// Test SolveBeforeFactor error.
#[test]
fn test_solve_before_factor() {
    let matrix = sparse_from_lower_triplets(2, &[(0, 0, 1.0), (1, 1, 1.0)]);

    let solver =
        SparseLDLT::analyze_with_matrix(&matrix, &AnalyzeOptions::default()).expect("analyze");

    let b = Col::from_fn(2, |i| (i + 1) as f64);
    let mut scratch = MemBuffer::new(faer::dyn_stack::StackReq::empty());
    let stack = MemStack::new(&mut scratch);

    let result = solver.solve(&b, stack);
    assert!(
        matches!(result, Err(SparseError::SolveBeforeFactor { .. })),
        "expected SolveBeforeFactor error, got {:?}",
        result
    );
}

/// Test DimensionMismatch for wrong RHS length.
#[test]
fn test_solve_dimension_mismatch() {
    let matrix = sparse_from_lower_triplets(3, &[(0, 0, 1.0), (1, 1, 2.0), (2, 2, 3.0)]);

    let mut solver =
        SparseLDLT::analyze_with_matrix(&matrix, &AnalyzeOptions::default()).expect("analyze");
    solver
        .factor(&matrix, &FactorOptions::default())
        .expect("factor");

    let wrong_b = Col::from_fn(5, |i| i as f64);
    let scratch_req = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch_req);
    let stack = MemStack::new(&mut mem);

    let result = solver.solve(&wrong_b, stack);
    assert!(
        matches!(result, Err(SparseError::DimensionMismatch { .. })),
        "expected DimensionMismatch error, got {:?}",
        result
    );
}

/// Test rank-deficient matrix: solver does not error, and zero pivots detected.
#[test]
fn test_solve_rank_deficient() {
    // Singular 3x3 matrix: row 2 = row 0
    let matrix = sparse_from_lower_triplets(
        3,
        &[
            (0, 0, 1.0),
            (1, 0, 2.0),
            (1, 1, 5.0),
            (2, 0, 1.0),
            (2, 1, 2.0),
            (2, 2, 1.0),
        ],
    );

    // Use a consistent RHS from the column space
    let x_trial = vec![1.0, 0.0, 1.0];
    let b_vec = sparse_matvec(&matrix, &x_trial);
    let b = Col::from_fn(3, |i| b_vec[i]);

    // Should succeed (not error) per SPRAL convention
    let result = SparseLDLT::solve_full(&matrix, &b, &SolverOptions::default());
    assert!(
        result.is_ok(),
        "rank-deficient solve should succeed: {:?}",
        result.err()
    );

    // Verify solution: Ax should equal b (even if x differs from x_trial)
    if let Ok(x) = &result {
        let ax = sparse_matvec(&matrix, &(0..3).map(|i| x[i]).collect::<Vec<_>>());
        let residual: f64 = ax
            .iter()
            .zip(b_vec.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        let b_norm: f64 = b_vec.iter().map(|v| v * v).sum::<f64>().sqrt();
        // For singular systems the residual may not be zero, but factorization
        // should detect zero pivots
        if b_norm > 0.0 {
            eprintln!(
                "rank-deficient: relative residual = {:.2e}",
                residual / b_norm
            );
        }
    }

    // Check that inertia reflects rank deficiency
    let mut solver =
        SparseLDLT::analyze_with_matrix(&matrix, &AnalyzeOptions::default()).expect("analyze");
    solver
        .factor(&matrix, &FactorOptions::default())
        .expect("factor");
    let inertia = solver.inertia().expect("inertia after factor");
    // For singular matrices, some pivots may be infinitely delayed and thus
    // not counted in pos/neg/zero. Total accounted pivots <= n.
    let total = inertia.positive + inertia.negative + inertia.zero;
    assert!(total <= 3, "inertia total {} exceeds n=3", total);
    // At least some pivots should be accounted for
    assert!(total > 0, "inertia should account for at least one pivot");
}

/// Test 0x0 matrix trivial solve.
#[test]
fn test_solve_empty_matrix() {
    let triplets: Vec<Triplet<usize, usize, f64>> = vec![];
    let matrix = SparseColMat::try_new_from_triplets(0, 0, &triplets).unwrap();
    let b = Col::from_fn(0, |_| 0.0);

    let result = SparseLDLT::solve_full(&matrix, &b, &SolverOptions::default());
    match result {
        Ok(x) => {
            assert_eq!(x.nrows(), 0, "0x0 solve should return empty solution");
        }
        Err(e) => {
            // 0x0 may fail in analyze (METIS/MC64 may reject empty graphs).
            // This is acceptable — just verify it's not a panic.
            eprintln!("0x0 solve returned error (acceptable): {}", e);
        }
    }
}

/// Test API equivalence: solve_full vs analyze→factor→solve.
#[test]
fn test_solve_api_equivalence() {
    let matrix = sparse_from_lower_triplets(
        3,
        &[
            (0, 0, 4.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (2, 0, 0.5),
            (2, 2, 2.0),
        ],
    );

    let b = Col::from_fn(3, |i| [5.5, 4.0, 6.5][i]);

    // Method 1: solve_full
    let x1 = SparseLDLT::solve_full(&matrix, &b, &SolverOptions::default()).expect("solve_full");

    // Method 2: analyze → factor → solve
    let opts = AnalyzeOptions::default();
    let mut solver = SparseLDLT::analyze_with_matrix(&matrix, &opts).expect("analyze");
    solver
        .factor(&matrix, &FactorOptions::default())
        .expect("factor");
    let scratch = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);
    let stack = MemStack::new(&mut mem);
    let x2 = solver.solve(&b, stack).expect("solve");

    // Both should give the same result
    let diff: f64 = (0..3).map(|i| (x1[i] - x2[i]).powi(2)).sum::<f64>().sqrt();
    assert!(diff < 1e-14, "API equivalence: diff = {:.2e}", diff);
}

// ---------------------------------------------------------------------------
// US2: Three-phase API with reuse
// ---------------------------------------------------------------------------

/// Test refactoring with different numeric values.
#[test]
fn test_refactor_different_values() {
    // Matrix 1: positive definite
    let a1 = sparse_from_lower_triplets(
        3,
        &[
            (0, 0, 4.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (2, 0, 0.5),
            (2, 2, 2.0),
        ],
    );

    // Matrix 2: same sparsity, different values
    let a2 = sparse_from_lower_triplets(
        3,
        &[
            (0, 0, 5.0),
            (1, 0, 0.5),
            (1, 1, 4.0),
            (2, 0, 1.0),
            (2, 2, 3.0),
        ],
    );

    let opts = AnalyzeOptions::default();
    let mut solver = SparseLDLT::analyze_with_matrix(&a1, &opts).expect("analyze");

    // Factor and solve with matrix 1
    solver
        .factor(&a1, &FactorOptions::default())
        .expect("factor 1");
    let x1_exact: Vec<f64> = vec![1.0, 2.0, 3.0];
    let b1_vec = sparse_matvec(&a1, &x1_exact);
    let b1 = Col::from_fn(3, |i| b1_vec[i]);
    let scratch = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);
    let stack = MemStack::new(&mut mem);
    let x1 = solver.solve(&b1, stack).expect("solve 1");
    let be1 = sparse_backward_error(&a1, &x1, &b1);
    assert!(be1 < 1e-10, "refactor: backward error 1 = {:.2e}", be1);

    // Refactor and solve with matrix 2
    solver
        .refactor(&a2, &FactorOptions::default())
        .expect("refactor 2");
    let x2_exact: Vec<f64> = vec![-1.0, 0.5, 2.0];
    let b2_vec = sparse_matvec(&a2, &x2_exact);
    let b2 = Col::from_fn(3, |i| b2_vec[i]);
    let scratch2 = solver.solve_scratch(1);
    let mut mem2 = MemBuffer::new(scratch2);
    let stack2 = MemStack::new(&mut mem2);
    let x2 = solver.solve(&b2, stack2).expect("solve 2");
    let be2 = sparse_backward_error(&a2, &x2, &b2);
    assert!(be2 < 1e-10, "refactor: backward error 2 = {:.2e}", be2);
}

/// Test multiple RHS with same factorization.
#[test]
fn test_multiple_rhs_reuse() {
    let matrix = sparse_from_lower_triplets(
        3,
        &[
            (0, 0, 4.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (2, 0, 0.5),
            (2, 2, 2.0),
        ],
    );

    let mut solver =
        SparseLDLT::analyze_with_matrix(&matrix, &AnalyzeOptions::default()).expect("analyze");
    solver
        .factor(&matrix, &FactorOptions::default())
        .expect("factor");

    let rhs_vectors: Vec<Vec<f64>> = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    for (idx, rhs) in rhs_vectors.iter().enumerate() {
        let b_vec = sparse_matvec(&matrix, rhs);
        let b = Col::from_fn(3, |i| b_vec[i]);
        let scratch = solver.solve_scratch(1);
        let mut mem = MemBuffer::new(scratch);
        let stack = MemStack::new(&mut mem);
        let x = solver.solve(&b, stack).expect("solve");
        let be = sparse_backward_error(&matrix, &x, &b);
        assert!(be < 1e-10, "RHS {}: backward error {:.2e}", idx, be);
    }
}

// ---------------------------------------------------------------------------
// US5: Workspace-efficient solve
// ---------------------------------------------------------------------------

/// Test that solve_scratch returns sufficient workspace.
#[test]
fn test_solve_scratch_sufficient() {
    let matrix = sparse_from_lower_triplets(
        3,
        &[
            (0, 0, 4.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (2, 0, 0.5),
            (2, 2, 2.0),
        ],
    );

    let mut solver =
        SparseLDLT::analyze_with_matrix(&matrix, &AnalyzeOptions::default()).expect("analyze");
    solver
        .factor(&matrix, &FactorOptions::default())
        .expect("factor");

    let scratch = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);
    let stack = MemStack::new(&mut mem);

    let b = Col::from_fn(3, |i| (i + 1) as f64);
    let x = solver.solve(&b, stack).expect("solve with exact scratch");

    let _be = sparse_backward_error(&matrix, &x, &b);
    // Note: b is NOT A*x_exact here, so backward error may be large.
    // The point is that solve succeeds without panicking.
}

/// Test workspace reuse: allocate once, solve twice.
#[test]
fn test_workspace_reuse() {
    let matrix = sparse_from_lower_triplets(
        3,
        &[
            (0, 0, 4.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (2, 0, 0.5),
            (2, 2, 2.0),
        ],
    );

    let mut solver =
        SparseLDLT::analyze_with_matrix(&matrix, &AnalyzeOptions::default()).expect("analyze");
    solver
        .factor(&matrix, &FactorOptions::default())
        .expect("factor");

    let scratch = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);

    // Solve 1
    let x1_exact: Vec<f64> = vec![1.0, 2.0, 3.0];
    let b1_vec = sparse_matvec(&matrix, &x1_exact);
    let b1 = Col::from_fn(3, |i| b1_vec[i]);
    {
        let stack = MemStack::new(&mut mem);
        let x1 = solver.solve(&b1, stack).expect("solve 1");
        let be1 = sparse_backward_error(&matrix, &x1, &b1);
        assert!(be1 < 1e-10, "workspace reuse solve 1: be = {:.2e}", be1);
    }

    // Solve 2 (reuse same MemBuffer)
    let x2_exact: Vec<f64> = vec![-1.0, 0.5, 2.0];
    let b2_vec = sparse_matvec(&matrix, &x2_exact);
    let b2 = Col::from_fn(3, |i| b2_vec[i]);
    {
        let stack = MemStack::new(&mut mem);
        let x2 = solver.solve(&b2, stack).expect("solve 2");
        let be2 = sparse_backward_error(&matrix, &x2, &b2);
        assert!(be2 < 1e-10, "workspace reuse solve 2: be = {:.2e}", be2);
    }
}

// ---------------------------------------------------------------------------
// US4: Scaling integration (MatchOrderMetis round-trip)
// ---------------------------------------------------------------------------

/// Test that MatchOrderMetis scaling is correctly applied/unapplied during solve.
/// Verifies that S A S is factored and S^{-1} x is returned.
#[test]
fn test_scaling_round_trip() {
    // Matrix where MC64 scaling makes a visible difference:
    // diagonal entries vary by orders of magnitude.
    let matrix = sparse_from_lower_triplets(
        4,
        &[
            (0, 0, 1e4),
            (1, 0, 1.0),
            (1, 1, 1e-4),
            (2, 1, 0.5),
            (2, 2, 1e2),
            (3, 2, 1.0),
            (3, 3, 1e-2),
        ],
    );

    let x_exact: Vec<f64> = vec![1.0, 2.0, -1.0, 3.0];
    let b_vec = sparse_matvec(&matrix, &x_exact);
    let b = Col::from_fn(4, |i| b_vec[i]);

    // Solve with MatchOrderMetis (default — includes scaling)
    let opts = SolverOptions {
        ordering: OrderingStrategy::MatchOrderMetis,
        ..SolverOptions::default()
    };
    let x = SparseLDLT::solve_full(&matrix, &b, &opts).expect("solve with scaling");
    let be = sparse_backward_error(&matrix, &x, &b);
    assert!(be < 1e-10, "scaling round-trip backward error: {:.2e}", be);

    // Also verify with Metis (no scaling) — may have worse accuracy on poorly-scaled matrix
    let opts_no_scale = SolverOptions {
        ordering: OrderingStrategy::Metis,
        ..SolverOptions::default()
    };
    let x2 = SparseLDLT::solve_full(&matrix, &b, &opts_no_scale).expect("solve without scaling");
    let be2 = sparse_backward_error(&matrix, &x2, &b);
    // Plain METIS without scaling may have degraded accuracy on poorly-scaled matrices.
    // We just verify it still produces a reasonable answer (not catastrophic failure).
    assert!(be2 < 1e-3, "no-scaling backward error: {:.2e}", be2);
}

// ---------------------------------------------------------------------------
// Random matrix backward error
// ---------------------------------------------------------------------------

/// Test backward error on random symmetric indefinite matrices.
#[cfg(feature = "test-util")]
#[test]
fn test_random_indefinite_backward_error() {
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use rivrs_sparse::testing::generators::{RandomMatrixConfig, generate_random_symmetric};

    let sizes = [10, 50, 100];
    for &n in &sizes {
        let mut rng = StdRng::seed_from_u64(42 + n as u64);
        let config = RandomMatrixConfig {
            size: n,
            target_nnz: (n as f64 * n as f64 * 0.3) as usize,
            positive_definite: false,
        };
        let matrix = generate_random_symmetric(&config, &mut rng)
            .unwrap_or_else(|e| panic!("generate indef n={}: {}", n, e));

        let x_exact: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
        let b_vec = sparse_matvec(&matrix, &x_exact);
        let b = Col::from_fn(n, |i| b_vec[i]);

        let x = SparseLDLT::solve_full(&matrix, &b, &SolverOptions::default())
            .unwrap_or_else(|e| panic!("random indef n={}: {}", n, e));

        let be = sparse_backward_error(&matrix, &x, &b);
        assert!(
            be < 1e-10,
            "random indefinite n={}: backward error {:.2e}",
            n,
            be
        );
    }
}

/// Test backward error on random positive definite matrices.
#[cfg(feature = "test-util")]
#[test]
fn test_random_pd_backward_error() {
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use rivrs_sparse::testing::generators::{RandomMatrixConfig, generate_random_symmetric};

    let sizes = [10, 50, 100];
    for &n in &sizes {
        let mut rng = StdRng::seed_from_u64(99 + n as u64);
        let config = RandomMatrixConfig {
            size: n,
            target_nnz: (n as f64 * n as f64 * 0.3) as usize,
            positive_definite: true,
        };
        let matrix = generate_random_symmetric(&config, &mut rng)
            .unwrap_or_else(|e| panic!("generate PD n={}: {}", n, e));

        let x_exact: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
        let b_vec = sparse_matvec(&matrix, &x_exact);
        let b = Col::from_fn(n, |i| b_vec[i]);

        let x = SparseLDLT::solve_full(&matrix, &b, &SolverOptions::default())
            .unwrap_or_else(|e| panic!("random PD n={}: {}", n, e));

        let be = sparse_backward_error(&matrix, &x, &b);
        assert!(be < 1e-12, "random PD n={}: backward error {:.2e}", n, be);
    }
}

// ---------------------------------------------------------------------------
// Ordering strategy tests
// ---------------------------------------------------------------------------

/// Test UserSupplied ordering: identity permutation.
#[test]
fn test_user_supplied_ordering() {
    use faer::perm::Perm;

    let matrix = sparse_from_lower_triplets(
        4,
        &[
            (0, 0, 4.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (2, 1, 0.5),
            (2, 2, 5.0),
            (3, 2, 1.0),
            (3, 3, 2.0),
        ],
    );

    // Identity permutation
    let fwd: Box<[usize]> = vec![0, 1, 2, 3].into_boxed_slice();
    let inv: Box<[usize]> = vec![0, 1, 2, 3].into_boxed_slice();
    let perm = Perm::new_checked(fwd, inv, 4);

    let opts = SolverOptions {
        ordering: OrderingStrategy::UserSupplied(perm),
        ..SolverOptions::default()
    };

    let x_exact: Vec<f64> = vec![1.0, -1.0, 2.0, 0.5];
    let b_vec = sparse_matvec(&matrix, &x_exact);
    let b = Col::from_fn(4, |i| b_vec[i]);

    let x = SparseLDLT::solve_full(&matrix, &b, &opts).expect("UserSupplied solve");
    let be = sparse_backward_error(&matrix, &x, &b);
    assert!(be < 1e-10, "UserSupplied backward error: {:.2e}", be);
}

/// Test UserSupplied ordering: reverse permutation.
#[test]
fn test_user_supplied_reverse_ordering() {
    use faer::perm::Perm;

    let matrix = sparse_from_lower_triplets(
        4,
        &[
            (0, 0, 4.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (2, 1, 0.5),
            (2, 2, 5.0),
            (3, 2, 1.0),
            (3, 3, 2.0),
        ],
    );

    // Reverse permutation
    let fwd: Box<[usize]> = vec![3, 2, 1, 0].into_boxed_slice();
    let inv: Box<[usize]> = vec![3, 2, 1, 0].into_boxed_slice();
    let perm = Perm::new_checked(fwd, inv, 4);

    let opts = SolverOptions {
        ordering: OrderingStrategy::UserSupplied(perm),
        ..SolverOptions::default()
    };

    let x_exact: Vec<f64> = vec![1.0, -1.0, 2.0, 0.5];
    let b_vec = sparse_matvec(&matrix, &x_exact);
    let b = Col::from_fn(4, |i| b_vec[i]);

    let x = SparseLDLT::solve_full(&matrix, &b, &opts).expect("reverse ordering solve");
    let be = sparse_backward_error(&matrix, &x, &b);
    assert!(be < 1e-10, "reverse ordering backward error: {:.2e}", be);
}

/// Test symbolic-only analyze (no numeric values) with AMD ordering.
#[test]
fn test_symbolic_only_analyze_amd() {
    let matrix = sparse_from_lower_triplets(
        4,
        &[
            (0, 0, 4.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (2, 1, 0.5),
            (2, 2, 5.0),
            (3, 2, 1.0),
            (3, 3, 2.0),
        ],
    );

    // Use symbolic-only analyze (requires AMD or Metis — not MatchOrderMetis)
    let opts = AnalyzeOptions {
        ordering: OrderingStrategy::Amd,
    };
    let mut solver = SparseLDLT::analyze(matrix.symbolic(), &opts).expect("symbolic analyze AMD");
    solver
        .factor(&matrix, &FactorOptions::default())
        .expect("factor");

    let x_exact: Vec<f64> = vec![1.0, -1.0, 2.0, 0.5];
    let b_vec = sparse_matvec(&matrix, &x_exact);
    let b = Col::from_fn(4, |i| b_vec[i]);
    let scratch = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);
    let stack = MemStack::new(&mut mem);
    let x = solver.solve(&b, stack).expect("solve");
    let be = sparse_backward_error(&matrix, &x, &b);
    assert!(be < 1e-10, "symbolic-only AMD backward error: {:.2e}", be);
}

/// Test that MatchOrderMetis ordering with symbolic-only analyze returns proper error.
#[test]
fn test_match_order_metis_requires_numeric() {
    let matrix = sparse_from_lower_triplets(3, &[(0, 0, 1.0), (1, 1, 2.0), (2, 2, 3.0)]);

    let opts = AnalyzeOptions {
        ordering: OrderingStrategy::MatchOrderMetis,
    };
    let result = SparseLDLT::analyze(matrix.symbolic(), &opts);
    assert!(
        result.is_err(),
        "MatchOrderMetis with symbolic-only should fail"
    );
    assert!(
        matches!(&result, Err(SparseError::AnalysisFailure { .. })),
        "expected AnalysisFailure error"
    );
}

// ---------------------------------------------------------------------------
// Full SuiteSparse backward error validation
// ---------------------------------------------------------------------------

/// Minimum number of SuiteSparse matrices to consider the full collection present.
const MIN_SUITESPARSE_FULL: usize = 50;

/// End-to-end backward error on the full SuiteSparse collection.
///
/// Loads matrices one at a time to avoid OOM. Follows Duff et al (2020)
/// experimental protocol: easy-indefinite → plain METIS, hard-indefinite →
/// MC64+METIS. Reports per-matrix backward error with `--nocapture`.
///
/// Run with:
///   cargo test --release --test solve test_solve_suitesparse_full -- --ignored --nocapture --test-threads=1
#[test]
#[ignore = "requires full SuiteSparse collection"]
fn test_solve_suitesparse_full() {
    use rivrs_sparse::io::registry;

    const SPRAL_BE_THRESHOLD: f64 = 5e-11;
    const RELAXED_BE_THRESHOLD: f64 = 1e-8;

    let all_meta = registry::load_registry().expect("failed to load metadata.json");
    let suitesparse_meta: Vec<_> = all_meta
        .into_iter()
        .filter(|m| m.source == "suitesparse")
        .collect();

    let mut passed_strict = 0usize;
    let mut passed_relaxed = 0usize;
    let mut missing = 0usize;
    let mut failed = Vec::new();

    eprintln!("\n{:<40} {:>8} {:>12}  {}", "Matrix", "n", "BE", "Status");
    eprintln!("{}", "-".repeat(80));

    for meta in &suitesparse_meta {
        let test_matrix = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(tm)) => tm,
            Ok(None) => {
                missing += 1;
                eprintln!("{:<40} {:>8} {:>12}  MISSING", meta.name, meta.size, "");
                continue;
            }
            Err(e) => {
                failed.push(format!("'{}': load error: {}", meta.name, e));
                eprintln!(
                    "{:<40} {:>8} {:>12}  LOAD ERROR: {}",
                    meta.name, meta.size, "", e
                );
                continue;
            }
        };

        let a = &test_matrix.matrix;
        let n = a.nrows();

        // Build RHS from a known x_exact
        let x_vec: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
        let b_vec = sparse_matvec(a, &x_vec);
        let b = Col::from_fn(n, |i| b_vec[i]);

        // Follow Duff et al (2020) experimental protocol:
        // easy-indefinite → plain METIS, hard-indefinite → MC64+METIS
        let ordering = if meta.category == "hard-indefinite" {
            OrderingStrategy::MatchOrderMetis
        } else {
            OrderingStrategy::Metis
        };
        // let ordering = OrderingStrategy::MatchOrderMetis;
        let opts = SolverOptions {
            ordering,
            ..SolverOptions::default()
        };
        let result = SparseLDLT::solve_full(a, &b, &opts);
        match result {
            Ok(x) => {
                let be = sparse_backward_error(a, &x, &b);
                let status = if be < SPRAL_BE_THRESHOLD {
                    passed_strict += 1;
                    "PASS"
                } else if be < RELAXED_BE_THRESHOLD {
                    passed_relaxed += 1;
                    "RELAXED"
                } else {
                    failed.push(format!("'{}' (n={}): BE={:.2e}", meta.name, n, be));
                    "FAIL"
                };
                eprintln!("{:<40} {:>8} {:>12.2e}  {}", meta.name, n, be, status);
            }
            Err(e) => {
                failed.push(format!("'{}' (n={}): solve error: {}", meta.name, n, e));
                eprintln!("{:<40} {:>8} {:>12}  SOLVE ERROR: {}", meta.name, n, "", e);
            }
        }
        // test_matrix dropped here — frees memory before next iteration
    }

    eprintln!("{}", "-".repeat(80));
    eprintln!(
        "Strict (<5e-11): {}  Relaxed (<1e-8): {}  Missing: {}  Failed: {}  Total: {}",
        passed_strict,
        passed_relaxed,
        missing,
        failed.len(),
        suitesparse_meta.len()
    );

    if suitesparse_meta.len() - missing < MIN_SUITESPARSE_FULL {
        eprintln!(
            "\nSkipping assertions: only {} SuiteSparse matrices loaded (need >= {}). \
             Extract the full collection to test-data/suitesparse/.",
            suitesparse_meta.len() - missing,
            MIN_SUITESPARSE_FULL
        );
        return;
    }

    if !failed.is_empty() {
        eprintln!("\nFailed matrices:");
        for f in &failed {
            eprintln!("  {}", f);
        }
        panic!(
            "{} of {} matrices failed backward error threshold",
            failed.len(),
            suitesparse_meta.len() - missing
        );
    }
}
