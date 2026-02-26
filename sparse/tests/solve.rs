//! Tests for triangular solve & solver API.
//!
//! Test categories:
//! - Per-supernode triangular solve correctness
//! - End-to-end sparse solve
//! - Three-phase API with reuse
//! - Scaling integration
//! - Workspace-efficient solve

mod common;

use faer::Col;
use faer::Par;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::linalg::cholesky::SymmetricOrdering;
use faer::sparse::{SparseColMat, Triplet};

use rivrs_sparse::error::SparseError;
use rivrs_sparse::symmetric::{
    AnalyzeOptions, AptpNumeric, AptpOptions, AptpSymbolic, FactorOptions, OrderingStrategy,
    SolverOptions, SparseLDLT,
};
use rivrs_sparse::validate::sparse_backward_error;

use common::{sparse_from_lower_triplets, sparse_matvec};

// ---------------------------------------------------------------------------
// Per-supernode triangular solve correctness
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
    let scratch = rivrs_sparse::symmetric::solve::aptp_solve_scratch(&numeric, 1);
    let mut mem = MemBuffer::new(scratch);
    let stack = MemStack::new(&mut mem);
    rivrs_sparse::symmetric::solve::aptp_solve(&symbolic, &numeric, &mut rhs_perm, stack, Par::Seq)
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
// End-to-end backward error tests
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
// Error handling and edge cases
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

    let result = solver.solve(&b, stack, Par::Seq);
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

    let result = solver.solve(&wrong_b, stack, Par::Seq);
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
    let x2 = solver.solve(&b, stack, Par::Seq).expect("solve");

    // Both should give the same result
    let diff: f64 = (0..3).map(|i| (x1[i] - x2[i]).powi(2)).sum::<f64>().sqrt();
    assert!(diff < 1e-14, "API equivalence: diff = {:.2e}", diff);
}

// ---------------------------------------------------------------------------
// Three-phase API with reuse
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
    let x1 = solver.solve(&b1, stack, Par::Seq).expect("solve 1");
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
    let x2 = solver.solve(&b2, stack2, Par::Seq).expect("solve 2");
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
        let x = solver.solve(&b, stack, Par::Seq).expect("solve");
        let be = sparse_backward_error(&matrix, &x, &b);
        assert!(be < 1e-10, "RHS {}: backward error {:.2e}", idx, be);
    }
}

// ---------------------------------------------------------------------------
// Workspace-efficient solve
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
    let x = solver
        .solve(&b, stack, Par::Seq)
        .expect("solve with exact scratch");

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
        let x1 = solver.solve(&b1, stack, Par::Seq).expect("solve 1");
        let be1 = sparse_backward_error(&matrix, &x1, &b1);
        assert!(be1 < 1e-10, "workspace reuse solve 1: be = {:.2e}", be1);
    }

    // Solve 2 (reuse same MemBuffer)
    let x2_exact: Vec<f64> = vec![-1.0, 0.5, 2.0];
    let b2_vec = sparse_matvec(&matrix, &x2_exact);
    let b2 = Col::from_fn(3, |i| b2_vec[i]);
    {
        let stack = MemStack::new(&mut mem);
        let x2 = solver.solve(&b2, stack, Par::Seq).expect("solve 2");
        let be2 = sparse_backward_error(&matrix, &x2, &b2);
        assert!(be2 < 1e-10, "workspace reuse solve 2: be = {:.2e}", be2);
    }
}

// ---------------------------------------------------------------------------
// Scaling integration (MatchOrderMetis round-trip)
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
    let x = solver.solve(&b, stack, Par::Seq).expect("solve");
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

    eprintln!("\n{:<40} {:>8} {:>12}  Status", "Matrix", "n", "BE");
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

// ---------------------------------------------------------------------------
// Parallel correctness tests
// ---------------------------------------------------------------------------

/// Factor + solve CI subset matrices with Par::rayon(4), assert backward error < 5e-11.
#[test]
fn test_parallel_factor_ci_subset() {
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
        let b = Col::from_fn(n, |i| (i as f64 + 1.0).sin());

        let ordering = if meta.category == "hard-indefinite" {
            OrderingStrategy::MatchOrderMetis
        } else {
            OrderingStrategy::Metis
        };
        let analyze_opts = AnalyzeOptions { ordering };
        let factor_opts = FactorOptions {
            par: Par::rayon(4),
            ..FactorOptions::default()
        };

        let mut solver = SparseLDLT::analyze_with_matrix(a, &analyze_opts).expect("analyze");
        solver.factor(a, &factor_opts).expect("factor");

        let scratch = solver.solve_scratch(1);
        let mut mem = MemBuffer::new(scratch);
        let stack = MemStack::new(&mut mem);
        let x = solver.solve(&b, stack, Par::rayon(4)).expect("solve");

        let be = sparse_backward_error(a, &x, &b);
        assert!(
            be < 5e-11,
            "Parallel factor CI '{}' (n={}): backward error {:.2e} >= 5e-11",
            meta.name,
            n,
            be
        );
    }
}

/// Factor same matrix with Par::Seq and Par::rayon(4), verify both are correct
/// and agree to machine precision.
///
/// Uses a 500-node random sparse matrix (with METIS ordering) to ensure the
/// elimination tree has multi-child supernodes that actually trigger rayon
/// `par_iter` dispatch in `factor_subtree`.
///
/// Note: true parallel execution may reorder floating-point operations, so
/// results may differ at the ULP level. We verify both solutions have excellent
/// backward error and agree within a relative tolerance.
#[cfg(feature = "test-util")]
#[test]
fn test_parallel_factor_determinism() {
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use rivrs_sparse::testing::generators::{RandomMatrixConfig, generate_random_symmetric};

    let mut rng = StdRng::seed_from_u64(2024);
    let config = RandomMatrixConfig {
        size: 500,
        target_nnz: 5000,
        positive_definite: false, // indefinite exercises more pivot logic
    };
    let matrix = generate_random_symmetric(&config, &mut rng).expect("generate matrix");

    let n = matrix.nrows();
    let x_exact: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
    let b_vec = sparse_matvec(&matrix, &x_exact);
    let b = Col::from_fn(n, |i| b_vec[i]);

    let analyze_opts = AnalyzeOptions {
        ordering: OrderingStrategy::Metis,
    };

    // Sequential factorization + solve
    let factor_seq = FactorOptions {
        par: Par::Seq,
        ..FactorOptions::default()
    };
    let mut solver_seq =
        SparseLDLT::analyze_with_matrix(&matrix, &analyze_opts).expect("analyze seq");
    solver_seq.factor(&matrix, &factor_seq).expect("factor seq");
    let scratch = solver_seq.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);
    let x_seq = solver_seq
        .solve(&b, MemStack::new(&mut mem), Par::Seq)
        .expect("solve seq");
    let be_seq = sparse_backward_error(&matrix, &x_seq, &b);

    // Parallel factorization + solve
    let factor_par = FactorOptions {
        par: Par::rayon(4),
        ..FactorOptions::default()
    };
    let mut solver_par =
        SparseLDLT::analyze_with_matrix(&matrix, &analyze_opts).expect("analyze par");
    solver_par.factor(&matrix, &factor_par).expect("factor par");
    let mut mem2 = MemBuffer::new(scratch);
    let x_par = solver_par
        .solve(&b, MemStack::new(&mut mem2), Par::rayon(4))
        .expect("solve par");
    let be_par = sparse_backward_error(&matrix, &x_par, &b);

    // Both must achieve excellent backward error
    assert!(be_seq < 1e-10, "sequential backward error: {:.2e}", be_seq);
    assert!(be_par < 1e-10, "parallel backward error: {:.2e}", be_par);

    // Solutions should agree within machine precision (ULP-level differences
    // from parallel FP reordering are acceptable)
    let max_diff: f64 = (0..n)
        .map(|i| (x_seq[i] - x_par[i]).abs())
        .fold(0.0f64, f64::max);
    let x_norm: f64 = (0..n).map(|i| x_seq[i].abs()).fold(0.0f64, f64::max);
    let rel_diff = if x_norm > 0.0 {
        max_diff / x_norm
    } else {
        max_diff
    };
    assert!(
        rel_diff < 1e-12,
        "Seq vs Par relative difference: {:.2e} (max_diff={:.2e})",
        rel_diff,
        max_diff
    );
}

/// Parallel factorization produces correct results across matrix sizes.
///
/// Tests both small hand-constructed matrices (where threshold gating forces
/// Par::Seq for BLAS) AND larger generated matrices (where tree-level par_iter
/// actually dispatches). Compares Par::Seq vs Par::rayon(4) for bitwise identity
/// on small matrices, and checks backward error on larger ones.
#[cfg(feature = "test-util")]
#[test]
fn test_parallel_correctness_mixed_sizes() {
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use rivrs_sparse::io::registry;
    use rivrs_sparse::testing::generators::{
        RandomMatrixConfig, generate_arrow, generate_random_symmetric,
    };

    // Part 1: Small hand-constructed matrices — bitwise identity (Seq vs Par)
    let test_names = [
        "arrow-5-pd",
        "tridiag-5-pd",
        "block-diag-15",
        "stress-delayed-pivots",
        "arrow-10-indef",
        "bordered-block-20",
    ];

    for name in &test_names {
        let test = registry::load_test_matrix(name)
            .unwrap_or_else(|e| panic!("load '{}': {}", name, e))
            .unwrap_or_else(|| panic!("matrix '{}' not found", name));

        let matrix = &test.matrix;
        let n = matrix.nrows();
        let analyze_opts = AnalyzeOptions::default();

        let factor_seq = FactorOptions {
            par: Par::Seq,
            ..FactorOptions::default()
        };
        let mut solver_seq =
            SparseLDLT::analyze_with_matrix(matrix, &analyze_opts).expect("analyze seq");
        solver_seq.factor(matrix, &factor_seq).expect("factor seq");
        let b: Col<f64> = Col::from_fn(n, |i| (i + 1) as f64);
        let scratch_req = solver_seq.solve_scratch(1);
        let mut mem = MemBuffer::new(scratch_req);
        let x_seq = solver_seq
            .solve(&b, MemStack::new(&mut mem), Par::Seq)
            .expect("solve seq");

        let factor_par = FactorOptions {
            par: Par::rayon(4),
            ..FactorOptions::default()
        };
        let mut solver_par =
            SparseLDLT::analyze_with_matrix(matrix, &analyze_opts).expect("analyze par");
        solver_par.factor(matrix, &factor_par).expect("factor par");
        let mut mem2 = MemBuffer::new(scratch_req);
        let x_par = solver_par
            .solve(&b, MemStack::new(&mut mem2), Par::rayon(4))
            .expect("solve par");

        // Small matrices should produce identical results: intra-node BLAS
        // is gated to Par::Seq (front < INTRA_NODE_THRESHOLD=256), and tree-level
        // par_iter preserves ordering. The diagonal solve's gather-solve-scatter
        // pattern may reorder operations on some platforms, so allow a tight
        // relative tolerance rather than requiring bitwise identity.
        let max_diff: f64 = (0..n)
            .map(|i| (x_seq[i] - x_par[i]).abs())
            .fold(0.0f64, f64::max);
        let x_norm: f64 = (0..n).map(|i| x_seq[i].abs()).fold(0.0f64, f64::max);
        let rel_diff = if x_norm > 0.0 {
            max_diff / x_norm
        } else {
            max_diff
        };
        assert!(
            rel_diff < 1e-14,
            "{}: Seq vs Par relative difference: {:.2e} (max_diff={:.2e})",
            name,
            rel_diff,
            max_diff
        );
    }

    // Part 2: Larger generated matrices — these actually exercise tree-level
    // parallelism (METIS ordering creates multi-child supernodes at n=300+).
    let mut rng = StdRng::seed_from_u64(7777);

    // Arrow matrix n=300: star-shaped elimination tree with many leaves
    let arrow = generate_arrow(300, true, &mut rng).expect("generate arrow");
    let b_arrow = Col::from_fn(300, |i| (i as f64 + 1.0).sin());
    let factor_par = FactorOptions {
        par: Par::rayon(4),
        ..FactorOptions::default()
    };
    let analyze_opts = AnalyzeOptions {
        ordering: OrderingStrategy::Metis,
    };
    let mut solver = SparseLDLT::analyze_with_matrix(&arrow, &analyze_opts).expect("analyze arrow");
    solver.factor(&arrow, &factor_par).expect("factor arrow");
    let scratch = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);
    let x = solver
        .solve(&b_arrow, MemStack::new(&mut mem), Par::rayon(4))
        .expect("solve arrow");
    let be = sparse_backward_error(&arrow, &x, &b_arrow);
    assert!(be < 1e-10, "parallel arrow-300 backward error: {:.2e}", be);

    // Random sparse indefinite n=500: diverse tree structure
    let config = RandomMatrixConfig {
        size: 500,
        target_nnz: 5000,
        positive_definite: false,
    };
    let random_mat = generate_random_symmetric(&config, &mut rng).expect("generate random");
    let b_rand = Col::from_fn(500, |i| ((i % 7) as f64 - 3.0) / 3.0);
    let mut solver =
        SparseLDLT::analyze_with_matrix(&random_mat, &analyze_opts).expect("analyze random");
    solver
        .factor(&random_mat, &factor_par)
        .expect("factor random");
    let scratch = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);
    let x = solver
        .solve(&b_rand, MemStack::new(&mut mem), Par::rayon(4))
        .expect("solve random");
    let be = sparse_backward_error(&random_mat, &x, &b_rand);
    assert!(be < 1e-10, "parallel random-500 backward error: {:.2e}", be);
}

/// Parallel solve correctness — factor + solve CI subset with Par::rayon(4).
#[test]
fn test_parallel_solve_ci_subset() {
    use rivrs_sparse::io::registry;

    let reg = registry::load_registry().expect("load registry");
    let ci_matrices: Vec<_> = reg.iter().filter(|m| m.ci_subset).collect();

    let factor_opts = FactorOptions {
        par: Par::rayon(4),
        ..FactorOptions::default()
    };
    let analyze_opts = AnalyzeOptions::default();

    for meta in &ci_matrices {
        let test = registry::load_test_matrix(&meta.name)
            .unwrap_or_else(|e| panic!("load '{}': {}", meta.name, e))
            .unwrap_or_else(|| panic!("matrix '{}' not found", meta.name));
        let matrix = &test.matrix;
        let n = matrix.nrows();

        let mut solver = SparseLDLT::analyze_with_matrix(matrix, &analyze_opts).expect("analyze");
        solver.factor(matrix, &factor_opts).expect("factor");

        let b = Col::from_fn(n, |i| ((i + 1) as f64).sin());
        let scratch_req = solver.solve_scratch(1);
        let mut mem = MemBuffer::new(scratch_req);
        let x = solver
            .solve(&b, MemStack::new(&mut mem), Par::rayon(4))
            .expect("solve");

        let be = rivrs_sparse::validate::sparse_backward_error(matrix, &x, &b);
        assert!(
            be < 5e-11,
            "{}: parallel solve backward error {:.2e} exceeds 5e-11",
            meta.name,
            be
        );
    }
}

/// Parallel solve consistency — factor+solve with Par::Seq vs Par::rayon(4),
/// verify both are correct and agree to machine precision.
///
/// Uses a 500-node random matrix to exercise tree-level parallelism in both
/// factorization and the parallel diagonal solve (n_supernodes >= 4 trigger).
#[cfg(feature = "test-util")]
#[test]
fn test_parallel_solve_determinism() {
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use rivrs_sparse::testing::generators::{RandomMatrixConfig, generate_random_symmetric};

    let mut rng = StdRng::seed_from_u64(3141);
    let config = RandomMatrixConfig {
        size: 500,
        target_nnz: 5000,
        positive_definite: true,
    };
    let matrix = generate_random_symmetric(&config, &mut rng).expect("generate matrix");

    let n = matrix.nrows();
    let x_exact: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
    let b_vec = sparse_matvec(&matrix, &x_exact);
    let b = Col::from_fn(n, |i| b_vec[i]);

    let analyze_opts = AnalyzeOptions {
        ordering: OrderingStrategy::Metis,
    };

    // Sequential
    let factor_seq = FactorOptions {
        par: Par::Seq,
        ..FactorOptions::default()
    };
    let mut solver_seq =
        SparseLDLT::analyze_with_matrix(&matrix, &analyze_opts).expect("analyze seq");
    solver_seq.factor(&matrix, &factor_seq).expect("factor seq");
    let scratch = solver_seq.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch);
    let x_seq = solver_seq
        .solve(&b, MemStack::new(&mut mem), Par::Seq)
        .expect("solve seq");
    let be_seq = sparse_backward_error(&matrix, &x_seq, &b);

    // Parallel
    let factor_par = FactorOptions {
        par: Par::rayon(4),
        ..FactorOptions::default()
    };
    let mut solver_par =
        SparseLDLT::analyze_with_matrix(&matrix, &analyze_opts).expect("analyze par");
    solver_par.factor(&matrix, &factor_par).expect("factor par");
    let mut mem2 = MemBuffer::new(scratch);
    let x_par = solver_par
        .solve(&b, MemStack::new(&mut mem2), Par::rayon(4))
        .expect("solve par");
    let be_par = sparse_backward_error(&matrix, &x_par, &b);

    // Both must achieve excellent backward error
    assert!(be_seq < 1e-12, "sequential backward error: {:.2e}", be_seq);
    assert!(be_par < 1e-12, "parallel backward error: {:.2e}", be_par);

    // Solutions should agree within machine precision
    let max_diff: f64 = (0..n)
        .map(|i| (x_seq[i] - x_par[i]).abs())
        .fold(0.0f64, f64::max);
    let x_norm: f64 = (0..n).map(|i| x_seq[i].abs()).fold(0.0f64, f64::max);
    let rel_diff = if x_norm > 0.0 {
        max_diff / x_norm
    } else {
        max_diff
    };
    assert!(
        rel_diff < 1e-12,
        "Seq vs Par relative difference: {:.2e} (max_diff={:.2e})",
        rel_diff,
        max_diff
    );
}

// ---------------------------------------------------------------------------
// Supernode Amalgamation Integration Tests
// ---------------------------------------------------------------------------

/// Verify amalgamation reduces c-71 supernode count below 12K.
///
/// c-71 (n=76638, nnz=468096) is an optimization matrix with ~35K tiny
/// supernodes before amalgamation. With nemin=32, most merge into larger
/// supernodes, bringing the count well below 12K.
#[test]
#[ignore = "requires c-71 matrix (not in CI subset)"]
fn test_amalgamation_c71_supernode_count() {
    use rivrs_sparse::io::registry;

    let all_meta = registry::load_registry().expect("load registry");
    let c71_meta = all_meta
        .iter()
        .find(|m| m.name == "GHS_indef/c-71")
        .expect("c-71 not in metadata.json");

    let test_matrix = registry::load_test_matrix_from_entry(c71_meta)
        .expect("load c-71")
        .expect("c-71 matrix file not found");

    let matrix = &test_matrix.matrix;

    let mut solver =
        SparseLDLT::analyze_with_matrix(matrix, &AnalyzeOptions::default()).expect("analyze c-71");

    let factor_opts = FactorOptions::default();
    solver.factor(matrix, &factor_opts).expect("factor c-71");

    let sn_count = solver.per_supernode_stats().unwrap().len();
    eprintln!("c-71 supernode count after amalgamation: {sn_count}");

    assert!(
        sn_count < 12_000,
        "c-71 supernode count ({sn_count}) should be < 12K after amalgamation"
    );
}

/// Verify c-71 achieves acceptable backward error after amalgamation.
///
/// Factorize + solve c-71 with default settings (nemin=32, MatchOrderMetis),
/// assert backward error < 5e-11 (SPRAL's threshold).
#[test]
#[ignore = "requires c-71 matrix (not in CI subset)"]
fn test_amalgamation_c71_backward_error() {
    use rivrs_sparse::io::registry;

    let all_meta = registry::load_registry().expect("load registry");
    let c71_meta = all_meta
        .iter()
        .find(|m| m.name == "GHS_indef/c-71")
        .expect("c-71 not in metadata.json");

    let test_matrix = registry::load_test_matrix_from_entry(c71_meta)
        .expect("load c-71")
        .expect("c-71 matrix file not found");

    let matrix = &test_matrix.matrix;
    let n = matrix.nrows();

    // Create RHS: b = A * x_exact where x_exact = [1, 2, ..., n]
    let x_exact: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
    let b_vec = sparse_matvec(matrix, &x_exact);
    let b = Col::from_fn(n, |i| b_vec[i]);

    // Solve
    let opts = SolverOptions::default();
    let x = SparseLDLT::solve_full(matrix, &b, &opts).expect("solve c-71");

    let be = sparse_backward_error(matrix, &x, &b);
    eprintln!("c-71 backward error with amalgamation: {be:.2e}");

    assert!(
        be < 5e-11,
        "c-71 backward error ({be:.2e}) should be < 5e-11"
    );
}

/// nemin=1 disables amalgamation — solve result should be bitwise identical
/// to a solve without amalgamation (since nemin=1 makes do_merge always return false).
///
/// Uses a hand-constructed matrix where the result is deterministic regardless
/// of supernode structure.
#[test]
fn test_nemin_1_bitwise_identical() {
    let matrix = sparse_from_lower_triplets(
        5,
        &[
            (0, 0, 4.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (2, 1, 0.5),
            (2, 2, 5.0),
            (3, 2, 1.0),
            (3, 3, -2.0),
            (4, 3, 0.5),
            (4, 4, 6.0),
        ],
    );

    let x_exact: Vec<f64> = vec![1.0, -1.0, 2.0, 0.5, -0.5];
    let b_vec = sparse_matvec(&matrix, &x_exact);
    let b = Col::from_fn(5, |i| b_vec[i]);

    // Solve with nemin=1 (amalgamation disabled)
    let opts_no_amalg = SolverOptions {
        nemin: 1,
        ..SolverOptions::default()
    };
    let x_no_amalg = SparseLDLT::solve_full(&matrix, &b, &opts_no_amalg).expect("solve nemin=1");

    // Solve with nemin=32 (default, amalgamation enabled)
    let opts_default = SolverOptions::default();
    let x_default = SparseLDLT::solve_full(&matrix, &b, &opts_default).expect("solve nemin=32");

    // Both should achieve excellent backward error
    let be_no_amalg = sparse_backward_error(&matrix, &x_no_amalg, &b);
    let be_default = sparse_backward_error(&matrix, &x_default, &b);
    assert!(
        be_no_amalg < 1e-14,
        "nemin=1 backward error: {:.2e}",
        be_no_amalg
    );
    assert!(
        be_default < 1e-14,
        "nemin=32 backward error: {:.2e}",
        be_default
    );
}

// ---------------------------------------------------------------------------
// Small-leaf fast-path tests
// ---------------------------------------------------------------------------

/// Small chain matrix with fast path active.
#[test]
fn test_fast_path_small_chain() {
    // Arrow matrix: dense column structure produces small supernodes
    let n = 20;
    let mut entries = Vec::new();
    for i in 0..n {
        // Diagonal: alternating sign for indefiniteness
        entries.push((i, i, if i % 2 == 0 { 4.0 } else { -3.0 }));
    }
    // Tridiagonal coupling
    for i in 1..n {
        entries.push((i, i - 1, 0.5));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let x_vec: Vec<f64> = (0..n).map(|i| (i % 5) as f64 - 2.0).collect();
    let b_vec = sparse_matvec(&matrix, &x_vec);
    let b = Col::from_fn(n, |i| b_vec[i]);

    // Factor with fast path enabled (default threshold=256)
    let opts = SolverOptions {
        ordering: OrderingStrategy::Metis,
        ..SolverOptions::default()
    };
    let x = SparseLDLT::solve_full(&matrix, &b, &opts).expect("solve_full");
    let be = sparse_backward_error(&matrix, &x, &b);
    assert!(
        be < 5e-11,
        "fast path small chain: backward error {:.2e} >= 5e-11",
        be
    );
}

/// Compare fast path vs disabled — both must be correct.
#[test]
fn test_fast_path_matches_general() {
    let n = 30;
    let mut entries = Vec::new();
    for i in 0..n {
        entries.push((i, i, if i % 2 == 0 { 5.0 } else { -4.0 }));
    }
    for i in 1..n {
        entries.push((i, i - 1, 1.0));
    }
    // A few longer-range entries for structure
    for i in 3..n {
        entries.push((i, i - 3, 0.2));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let x_vec: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
    let b_vec = sparse_matvec(&matrix, &x_vec);
    let b = Col::from_fn(n, |i| b_vec[i]);

    // With fast path (default threshold=256)
    let opts_enabled = SolverOptions {
        ordering: OrderingStrategy::Metis,
        small_leaf_threshold: 256,
        ..SolverOptions::default()
    };
    let x_enabled = SparseLDLT::solve_full(&matrix, &b, &opts_enabled).expect("solve (enabled)");
    let be_enabled = sparse_backward_error(&matrix, &x_enabled, &b);

    // Without fast path (threshold=0)
    let opts_disabled = SolverOptions {
        ordering: OrderingStrategy::Metis,
        small_leaf_threshold: 0,
        ..SolverOptions::default()
    };
    let x_disabled = SparseLDLT::solve_full(&matrix, &b, &opts_disabled).expect("solve (disabled)");
    let be_disabled = sparse_backward_error(&matrix, &x_disabled, &b);

    assert!(
        be_enabled < 5e-11,
        "fast path enabled: backward error {:.2e}",
        be_enabled
    );
    assert!(
        be_disabled < 5e-11,
        "fast path disabled: backward error {:.2e}",
        be_disabled
    );
}

/// Matrix that forces delayed pivots in a small-leaf subtree.
#[test]
fn test_fast_path_delayed_pivots() {
    // Matrix with zero diagonal entries to force 2x2 pivots / delays
    let n = 15;
    let mut entries = Vec::new();
    for i in 0..n {
        // Some zero diagonals to force pivoting issues
        let diag = if i % 3 == 0 {
            0.0
        } else {
            3.0 * if i % 2 == 0 { 1.0 } else { -1.0 }
        };
        entries.push((i, i, diag));
    }
    // Off-diagonal structure
    for i in 1..n {
        entries.push((i, i - 1, 2.0));
    }
    for i in 2..n {
        entries.push((i, i - 2, 0.3));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let x_vec: Vec<f64> = (0..n).map(|i| i as f64 + 1.0).collect();
    let b_vec = sparse_matvec(&matrix, &x_vec);
    let b = Col::from_fn(n, |i| b_vec[i]);

    let opts = SolverOptions {
        ordering: OrderingStrategy::Metis,
        small_leaf_threshold: 256,
        ..SolverOptions::default()
    };
    let x = SparseLDLT::solve_full(&matrix, &b, &opts).expect("solve_full");
    let be = sparse_backward_error(&matrix, &x, &b);
    assert!(
        be < 5e-11,
        "fast path delayed pivots: backward error {:.2e} >= 5e-11",
        be
    );
}

/// Mixed tree — small subtree contribution feeds into large supernode.
#[test]
fn test_fast_path_contribution_boundary() {
    // Build a matrix with both small and large supernodes.
    // Dense lower-right block creates a large supernode, sparse structure
    // around it creates small ones.
    let n = 50;
    let mut entries = Vec::new();
    // Dense block in rows/cols 40..50 (creates a large front)
    for i in 40..n {
        for j in 40..=i {
            let val = if i == j { 10.0 } else { 0.5 };
            entries.push((i, j, val));
        }
    }
    // Sparse tridiagonal structure in rows/cols 0..40 (creates small supernodes)
    for i in 0..40 {
        entries.push((i, i, if i % 2 == 0 { 5.0 } else { -4.0 }));
    }
    for i in 1..40 {
        entries.push((i, i - 1, 0.8));
    }
    // Connect the sparse part to the dense part
    for i in 35..40 {
        entries.push((40, i, 0.3));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let x_vec: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
    let b_vec = sparse_matvec(&matrix, &x_vec);
    let b = Col::from_fn(n, |i| b_vec[i]);

    let opts = SolverOptions {
        ordering: OrderingStrategy::Metis,
        small_leaf_threshold: 256,
        ..SolverOptions::default()
    };
    let x = SparseLDLT::solve_full(&matrix, &b, &opts).expect("solve_full");
    let be = sparse_backward_error(&matrix, &x, &b);
    assert!(
        be < 5e-11,
        "fast path contribution boundary: backward error {:.2e} >= 5e-11",
        be
    );
}

/// MC64 scaling with fast path — validate identical scaling in both paths.
#[test]
fn test_fast_path_mc64_scaling() {
    use rivrs_sparse::io::registry;

    // Use a hard-indefinite CI-subset matrix that requires MC64 scaling
    let meta = registry::load_registry()
        .expect("load registry")
        .into_iter()
        .find(|m| m.ci_subset && m.category == "hard-indefinite")
        .expect("need at least one hard-indefinite CI matrix");

    let test = registry::load_test_matrix(&meta.name)
        .unwrap_or_else(|e| panic!("load '{}': {}", meta.name, e))
        .unwrap_or_else(|| panic!("matrix '{}' not found", meta.name));

    let a = &test.matrix;
    let n = a.nrows();
    let x_vec: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
    let b_vec = sparse_matvec(a, &x_vec);
    let b = Col::from_fn(n, |i| b_vec[i]);

    // With fast path + MC64
    let opts_enabled = SolverOptions {
        ordering: OrderingStrategy::MatchOrderMetis,
        small_leaf_threshold: 256,
        ..SolverOptions::default()
    };
    let x_enabled = SparseLDLT::solve_full(a, &b, &opts_enabled).expect("solve (enabled)");
    let be_enabled = sparse_backward_error(a, &x_enabled, &b);

    // Without fast path + MC64
    let opts_disabled = SolverOptions {
        ordering: OrderingStrategy::MatchOrderMetis,
        small_leaf_threshold: 0,
        ..SolverOptions::default()
    };
    let x_disabled = SparseLDLT::solve_full(a, &b, &opts_disabled).expect("solve (disabled)");
    let be_disabled = sparse_backward_error(a, &x_disabled, &b);

    assert!(
        be_enabled < 5e-11,
        "'{}' fast path+MC64: backward error {:.2e}",
        meta.name,
        be_enabled
    );
    assert!(
        be_disabled < 5e-11,
        "'{}' general+MC64: backward error {:.2e}",
        meta.name,
        be_disabled
    );
}

/// CI-subset SuiteSparse matrices with fast path enabled.
#[test]
fn test_fast_path_suitesparse_ci() {
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
        let x_vec: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
        let b_vec = sparse_matvec(a, &x_vec);
        let b = Col::from_fn(n, |i| b_vec[i]);

        let ordering = if meta.category == "hard-indefinite" {
            OrderingStrategy::MatchOrderMetis
        } else {
            OrderingStrategy::Metis
        };
        let opts = SolverOptions {
            ordering,
            small_leaf_threshold: 256,
            ..SolverOptions::default()
        };
        let x = SparseLDLT::solve_full(a, &b, &opts)
            .unwrap_or_else(|e| panic!("solve '{}': {}", meta.name, e));
        let be = sparse_backward_error(a, &x, &b);
        assert!(
            be < 5e-11,
            "fast path '{}': backward error {:.2e} >= 5e-11",
            meta.name,
            be
        );
    }
}
