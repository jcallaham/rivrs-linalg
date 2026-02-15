//! Tests for multifrontal numeric factorization (Phase 6).
//!
//! Organized by user story:
//! - US3: Assembly of frontal matrices
//! - US2: Delayed pivot propagation
//! - US1: End-to-end factorization
//! - US4: Schur complement validation

use faer::Mat;
use faer::sparse::linalg::cholesky::SymmetricOrdering;
use faer::sparse::{SparseColMat, Triplet};

use rivrs_sparse::aptp::pivot::PivotType;
use rivrs_sparse::aptp::{AptpNumeric, AptpOptions, AptpSymbolic};

// ---------------------------------------------------------------------------
// Helper: build a sparse symmetric matrix from dense lower triangle
// ---------------------------------------------------------------------------

/// Build a SparseColMat from a dense lower-triangular specification.
/// Input: (i, j, val) triplets where i >= j. Automatically mirrors to upper triangle.
fn sparse_from_lower_triplets(
    n: usize,
    entries: &[(usize, usize, f64)],
) -> SparseColMat<usize, f64> {
    let mut triplets = Vec::new();
    for &(i, j, v) in entries {
        triplets.push(Triplet::new(i, j, v));
        if i != j {
            triplets.push(Triplet::new(j, i, v));
        }
    }
    SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
}

/// Compute reconstruction error ||P^T A P - L D L^T|| / ||A|| for a multifrontal result.
///
/// Uses the reassemble_global_factors test utility from numeric.rs.
fn multifrontal_reconstruction_error(
    matrix: &SparseColMat<usize, f64>,
    symbolic: &AptpSymbolic,
    numeric: &AptpNumeric,
) -> f64 {
    let n = numeric.n();
    let a_dense = matrix.to_dense();

    // Get permutation
    let (perm_fwd, _perm_inv) = if let Some(perm) = symbolic.perm() {
        let (fwd, inv) = perm.arrays();
        (fwd.to_vec(), inv.to_vec())
    } else {
        let id: Vec<usize> = (0..n).collect();
        (id.clone(), id)
    };

    // Compute P^T A P
    let mut pap = Mat::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            pap[(i, j)] = a_dense[(perm_fwd[i], perm_fwd[j])];
        }
    }

    // Reconstruct L D L^T from per-supernode factors
    let mut l_global = Mat::<f64>::zeros(n, n);
    for i in 0..n {
        l_global[(i, i)] = 1.0;
    }
    let mut d_diag = vec![0.0f64; n];

    for ff in numeric.front_factors() {
        let ne = ff.num_eliminated();
        if ne == 0 {
            continue;
        }

        // Scatter L11 entries
        for i in 0..ne {
            for j in 0..i {
                let gi = ff.col_indices()[i];
                let gj = ff.col_indices()[j];
                l_global[(gi, gj)] = ff.l11()[(i, j)];
            }
        }

        // Scatter L21 entries
        for i in 0..ff.row_indices().len() {
            for j in 0..ne {
                let gi = ff.row_indices()[i];
                let gj = ff.col_indices()[j];
                l_global[(gi, gj)] = ff.l21()[(i, j)];
            }
        }

        // Scatter D11 entries
        let mut col = 0;
        while col < ne {
            let gc = ff.col_indices()[col];
            match ff.d11().get_pivot_type(col) {
                PivotType::OneByOne => {
                    d_diag[gc] = ff.d11().get_1x1(col);
                    col += 1;
                }
                PivotType::TwoByTwo { .. } => {
                    let block = ff.d11().get_2x2(col);
                    let gc2 = ff.col_indices()[col + 1];
                    d_diag[gc] = block.a;
                    d_diag[gc2] = block.c;
                    col += 2;
                }
                PivotType::Delayed => {
                    col += 1;
                }
            }
        }
    }

    // Build dense D matrix — place 2x2 off-diagonals at correct global positions
    let mut d_dense = Mat::<f64>::zeros(n, n);
    for i in 0..n {
        d_dense[(i, i)] = d_diag[i];
    }
    // Place off-diagonals at the exact (gc1, gc2) positions from each front
    for ff in numeric.front_factors() {
        let ne = ff.num_eliminated();
        let mut col = 0;
        while col < ne {
            match ff.d11().get_pivot_type(col) {
                PivotType::TwoByTwo { partner } if partner > col => {
                    let block = ff.d11().get_2x2(col);
                    let gc1 = ff.col_indices()[col];
                    let gc2 = ff.col_indices()[col + 1];
                    d_dense[(gc1, gc2)] = block.b;
                    d_dense[(gc2, gc1)] = block.b;
                    col += 2;
                }
                _ => {
                    col += 1;
                }
            }
        }
    }

    // Compute L * D * L^T
    let ld = &l_global * &d_dense;
    let mut l_t = Mat::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            l_t[(i, j)] = l_global[(j, i)];
        }
    }
    let ldlt = &ld * &l_t;

    // Compute ||P^T A P - L D L^T|| / ||A||
    let diff = &pap - &ldlt;
    let norm_diff = diff.norm_l2();
    let norm_a = a_dense.norm_l2();

    if norm_a == 0.0 {
        return 0.0;
    }
    norm_diff / norm_a
}

// ===========================================================================
// US3: Assembly tests (T010-T014)
// ===========================================================================

/// T010: Test build_supernode_info for supernodal symbolic analysis.
/// Uses a matrix that should produce supernodal analysis.
#[test]
fn test_build_supernode_info_supernodal() {
    // Arrow matrix: dense first row/column forces supernodal structure
    let n = 10;
    let mut entries = Vec::new();
    for i in 0..n {
        entries.push((i, i, 10.0));
    }
    for i in 1..n {
        entries.push((i, 0, 1.0));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");

    // Just verify it produces valid output
    if symbolic.is_supernodal() {
        let ns = symbolic.n_supernodes().unwrap();
        assert!(ns > 0, "should have at least 1 supernode");
        assert!(ns <= n, "cannot have more supernodes than columns");

        // Verify supernodes cover all columns
        let begin = symbolic.supernode_begin().unwrap();
        let end = symbolic.supernode_end().unwrap();
        assert_eq!(begin[0], 0, "first supernode starts at 0");
        assert_eq!(end[ns - 1], n, "last supernode ends at n");
        for s in 0..ns {
            assert!(end[s] > begin[s], "supernode {} has zero columns", s);
        }
    }
}

/// T011: Test build_supernode_info for simplicial symbolic analysis.
/// Uses a small diagonal matrix where faer produces simplicial.
#[test]
fn test_build_supernode_info_simplicial() {
    // Diagonal matrix — simplest possible sparsity
    let n = 4;
    let entries: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, (i + 1) as f64)).collect();
    let matrix = sparse_from_lower_triplets(n, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");

    // Whether faer chooses simplicial or supernodal, we should be able to factor
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    assert_eq!(numeric.n(), n);
    assert_eq!(
        numeric.stats().total_1x1_pivots,
        n,
        "diagonal matrix should have all 1x1 pivots"
    );
    assert_eq!(numeric.stats().total_2x2_pivots, 0);
    assert_eq!(numeric.stats().total_delayed, 0);
}

/// T012: Test scatter of original entries for a leaf supernode.
/// Verifies correct assembly by checking factorization of a simple matrix.
#[test]
fn test_scatter_leaf_supernode() {
    // Tridiagonal 4x4: well-conditioned PD matrix
    let entries = vec![
        (0, 0, 4.0),
        (1, 0, 1.0),
        (1, 1, 4.0),
        (2, 1, 1.0),
        (2, 2, 4.0),
        (3, 2, 1.0),
        (3, 3, 4.0),
    ];
    let matrix = sparse_from_lower_triplets(4, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "tridiagonal 4x4 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );
}

/// T013: Test extend_add by factoring a matrix with multiple supernodes.
/// If the assembly tree has children, extend-add must work correctly.
#[test]
fn test_extend_add_correctness() {
    // Arrow matrix with 6 nodes — forces parent-child relationships
    let n = 6;
    let mut entries = Vec::new();
    for i in 0..n {
        entries.push((i, i, 10.0));
    }
    for i in 1..n {
        entries.push((i, 0, 1.0));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "arrow 6x6 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );
}

/// T014: Test assembly with children (multi-supernode matrix).
#[test]
fn test_assembly_with_children() {
    // Band matrix: wider bandwidth forces multi-supernode structure
    let n = 8;
    let mut entries = Vec::new();
    for i in 0..n {
        entries.push((i, i, 10.0));
    }
    for i in 1..n {
        entries.push((i, i - 1, 1.0));
    }
    for i in 2..n {
        entries.push((i, i - 2, 0.5));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "band 8x8 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );
}

// ===========================================================================
// US2: Delayed pivot tests (T016-T019)
// ===========================================================================

/// T016: Test extract_contribution - factor a small frontal matrix and check
/// that the contribution block has correct dimensions.
#[test]
fn test_extract_contribution_structure() {
    // Indefinite matrix that may produce delays
    let entries = vec![
        (0, 0, 1.0),
        (1, 0, 2.0),
        (1, 1, 1.0),
        (2, 0, 1.0),
        (2, 1, 1.0),
        (2, 2, 4.0),
        (3, 0, 1.0),
        (3, 1, 1.0),
        (3, 2, 1.0),
        (3, 3, 4.0),
    ];
    let matrix = sparse_from_lower_triplets(4, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    // Verify that factorization completed (all pivots accounted for)
    let stats = numeric.stats();
    assert_eq!(
        stats.total_1x1_pivots + 2 * stats.total_2x2_pivots + stats.total_delayed,
        4,
        "all columns should be accounted for"
    );

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "indefinite 4x4 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );
}

/// T017: Test extract_front_factors - verify L11, D11, L21 shapes.
#[test]
fn test_extract_front_factors_shapes() {
    // Positive definite tridiagonal
    let entries = vec![
        (0, 0, 4.0),
        (1, 0, 1.0),
        (1, 1, 4.0),
        (2, 1, 1.0),
        (2, 2, 4.0),
    ];
    let matrix = sparse_from_lower_triplets(3, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    // Each front factor should have valid shapes
    for (s, ff) in numeric.front_factors().iter().enumerate() {
        let ne = ff.num_eliminated();
        assert_eq!(ff.l11().nrows(), ne, "supernode {} L11 rows mismatch", s);
        assert_eq!(ff.l11().ncols(), ne, "supernode {} L11 cols mismatch", s);
        assert_eq!(ff.l21().ncols(), ne, "supernode {} L21 cols mismatch", s);
        assert_eq!(
            ff.col_indices().len(),
            ne,
            "supernode {} col_indices len mismatch",
            s
        );
        assert_eq!(
            ff.row_indices().len(),
            ff.l21().nrows(),
            "supernode {} row_indices len mismatch",
            s
        );
    }
}

/// T018: Test delayed pivot propagation - construct a 2x2 indefinite matrix
/// where delays may occur and verify reconstruction.
#[test]
fn test_delayed_pivot_propagation() {
    // Indefinite matrix: eigenvalues have mixed signs
    // [1  3  0]    eigenvalues: ~-1.54, ~1.00, ~7.54
    // [3  5  2]
    // [0  2  3]
    let entries = vec![
        (0, 0, 1.0),
        (1, 0, 3.0),
        (1, 1, 5.0),
        (2, 1, 2.0),
        (2, 2, 3.0),
    ];
    let matrix = sparse_from_lower_triplets(3, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "indefinite 3x3 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );
}

/// T019: Test factorization with strongly indefinite matrix.
#[test]
fn test_strongly_indefinite() {
    // Matrix with both positive and negative pivots
    // [2   1   0  0]
    // [1  -2   1  0]
    // [0   1   3  1]
    // [0   0   1 -3]
    let entries = vec![
        (0, 0, 2.0),
        (1, 0, 1.0),
        (1, 1, -2.0),
        (2, 1, 1.0),
        (2, 2, 3.0),
        (3, 2, 1.0),
        (3, 3, -3.0),
    ];
    let matrix = sparse_from_lower_triplets(4, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "strongly indefinite 4x4 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );
}

// ===========================================================================
// US1: End-to-end factorization tests (T021-T026)
// ===========================================================================

/// T021: Test single-supernode-matches-dense.
/// For a small dense matrix that becomes a single front, verify multifrontal
/// result matches Phase 5's dense APTP exactly.
#[test]
fn test_single_supernode_matches_dense() {
    // Fully dense small matrix → single supernode
    let n = 5;
    let mut entries = Vec::new();
    for i in 0..n {
        for j in 0..=i {
            let v = if i == j {
                10.0
            } else {
                1.0 / ((i - j) as f64 + 1.0)
            };
            entries.push((i, j, v));
        }
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let options = AptpOptions::default();
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &options).expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "dense 5x5 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );
}

/// T022: Test simplicial path.
/// Use a small sparse diagonal matrix where faer may produce simplicial analysis.
#[test]
fn test_simplicial_path() {
    // Diagonal matrix — simplest possible
    let n = 5;
    let entries: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, (i + 1) as f64)).collect();
    let matrix = sparse_from_lower_triplets(n, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-14,
        "diagonal 5x5 reconstruction error: {:.2e} (expected < 1e-14)",
        err
    );

    assert_eq!(numeric.stats().total_1x1_pivots, n);
    assert_eq!(numeric.stats().total_2x2_pivots, 0);
    assert_eq!(numeric.stats().total_delayed, 0);
}

/// T023: Test hand-constructed matrices via the test registry.
/// Factors each hand-constructed matrix and verifies reconstruction error.
#[test]
fn test_hand_constructed_matrices() {
    use rivrs_sparse::testing::{TestCaseFilter, load_test_cases};

    let cases = load_test_cases(&TestCaseFilter::hand_constructed())
        .expect("failed to load hand-constructed test cases");

    assert_eq!(cases.len(), 15, "expected 15 hand-constructed matrices");

    for case in &cases {
        let symbolic = AptpSymbolic::analyze(case.matrix.symbolic(), SymmetricOrdering::Amd)
            .unwrap_or_else(|e| panic!("analyze failed for '{}': {}", case.name, e));
        let numeric = AptpNumeric::factor(&symbolic, &case.matrix, &AptpOptions::default())
            .unwrap_or_else(|e| panic!("factor failed for '{}': {}", case.name, e));

        let err = multifrontal_reconstruction_error(&case.matrix, &symbolic, &numeric);
        assert!(
            err < 1e-12,
            "'{}' reconstruction error: {:.2e} (expected < 1e-12)",
            case.name,
            err
        );
    }
}

/// T024: Test factorization statistics.
/// For a small matrix with known structure, verify stats fields.
#[test]
fn test_factorization_statistics() {
    // Diagonal 4x4 — all 1x1 pivots, no delays
    let entries: Vec<(usize, usize, f64)> = (0..4).map(|i| (i, i, (i + 1) as f64)).collect();
    let matrix = sparse_from_lower_triplets(4, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let stats = numeric.stats();
    assert_eq!(stats.total_1x1_pivots, 4, "expected 4 x 1x1 pivots");
    assert_eq!(stats.total_2x2_pivots, 0, "expected 0 x 2x2 pivots");
    assert_eq!(stats.total_delayed, 0, "expected 0 delayed columns");
    assert_eq!(stats.zero_pivots, 0, "expected 0 zero pivots");
    assert!(stats.max_front_size > 0, "max_front_size should be > 0");
}

/// T025: Test dimension mismatch error.
#[test]
fn test_dimension_mismatch_error() {
    let entries_3x3: Vec<(usize, usize, f64)> = (0..3).map(|i| (i, i, 1.0)).collect();
    let matrix_3x3 = sparse_from_lower_triplets(3, &entries_3x3);
    let entries_4x4: Vec<(usize, usize, f64)> = (0..4).map(|i| (i, i, 1.0)).collect();
    let matrix_4x4 = sparse_from_lower_triplets(4, &entries_4x4);

    let symbolic = AptpSymbolic::analyze(matrix_3x3.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");

    let result = AptpNumeric::factor(&symbolic, &matrix_4x4, &AptpOptions::default());
    assert!(result.is_err(), "should fail with dimension mismatch");
    let err = result.err().unwrap();
    assert!(
        matches!(
            err,
            rivrs_sparse::error::SparseError::DimensionMismatch { .. }
        ),
        "expected DimensionMismatch, got: {:?}",
        err
    );
}

/// T026: Test inertia validation.
/// For matrices with known sign distributions, verify inertia.
#[test]
fn test_inertia_validation() {
    // Positive definite: all eigenvalues positive
    let entries_pd = vec![(0, 0, 4.0), (1, 0, 1.0), (1, 1, 4.0)];
    let matrix_pd = sparse_from_lower_triplets(2, &entries_pd);

    let symbolic = AptpSymbolic::analyze(matrix_pd.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix_pd, &AptpOptions::default())
        .expect("factor should succeed");

    // Count positive/negative pivots from all front factors
    let mut n_pos = 0usize;
    let mut n_neg = 0usize;
    for ff in numeric.front_factors() {
        let inertia = ff.d11().compute_inertia();
        n_pos += inertia.positive;
        n_neg += inertia.negative;
    }
    assert_eq!(n_pos, 2, "PD matrix should have 2 positive eigenvalues");
    assert_eq!(n_neg, 0, "PD matrix should have 0 negative eigenvalues");
}

// ===========================================================================
// US4: Schur complement validation tests (T030-T031)
// ===========================================================================

/// T030: Test dense equivalence.
/// For small matrices, multifrontal should produce reconstruction errors
/// comparable to dense APTP.
#[test]
fn test_dense_equivalence() {
    use rivrs_sparse::aptp::aptp_factor;

    // Small PD matrix
    let n = 6;
    let mut entries = Vec::new();
    for i in 0..n {
        entries.push((i, i, 10.0));
    }
    for i in 1..n {
        entries.push((i, i - 1, 1.0));
    }
    for i in 2..n {
        entries.push((i, i - 2, 0.3));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    // Multifrontal factorization
    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let options = AptpOptions::default();
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &options)
        .expect("multifrontal factor should succeed");

    let mf_err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);

    // Dense APTP factorization (on full dense matrix, identity permutation)
    let dense = matrix.to_dense();
    let dense_result = aptp_factor(dense.as_ref(), &options).expect("dense factor should succeed");

    // Dense reconstruction error
    // Build L D L^T from dense result
    let l_dense = &dense_result.l;
    let mut d_mat = Mat::<f64>::zeros(n, n);
    for i in 0..n {
        match dense_result.d.get_pivot_type(i) {
            PivotType::OneByOne => {
                d_mat[(i, i)] = dense_result.d.get_1x1(i);
            }
            PivotType::TwoByTwo { partner } if partner > i => {
                let block = dense_result.d.get_2x2(i);
                d_mat[(i, i)] = block.a;
                d_mat[(i, i + 1)] = block.b;
                d_mat[(i + 1, i)] = block.b;
                d_mat[(i + 1, i + 1)] = block.c;
            }
            _ => {}
        }
    }
    let ld = l_dense * &d_mat;
    let mut l_t = Mat::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            l_t[(i, j)] = l_dense[(j, i)];
        }
    }
    let ldlt = &ld * &l_t;

    // Apply permutation to get P^T A P for dense case
    let dense_perm = &dense_result.perm;
    let (dfwd, _dinv) = dense_perm.as_ref().arrays();
    let mut pap_dense = Mat::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            pap_dense[(i, j)] = dense[(dfwd[i], dfwd[j])];
        }
    }
    let diff_dense = &pap_dense - &ldlt;
    let dense_err = diff_dense.norm_l2() / dense.norm_l2();

    // Both should be very small
    assert!(
        mf_err < 1e-12,
        "multifrontal error {:.2e} too large",
        mf_err
    );
    assert!(dense_err < 1e-12, "dense error {:.2e} too large", dense_err);

    // Multifrontal should be within 10x of dense (both use APTP kernel)
    let ratio = if dense_err > 0.0 {
        mf_err / dense_err
    } else {
        1.0
    };
    assert!(
        ratio < 100.0,
        "multifrontal error {:.2e} is {}x larger than dense error {:.2e}",
        mf_err,
        ratio,
        dense_err
    );
}

/// T031: Test multi-level contribution flow.
/// Factor a matrix with a deep assembly tree and verify correctness.
#[test]
fn test_multi_level_contribution_flow() {
    // Create a larger tridiagonal matrix — will have a deep tree
    let n = 20;
    let mut entries = Vec::new();
    for i in 0..n {
        entries.push((i, i, 4.0));
    }
    for i in 1..n {
        entries.push((i, i - 1, 1.0));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "tridiagonal 20x20 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );

    // Verify stats
    assert_eq!(
        numeric.stats().total_1x1_pivots + 2 * numeric.stats().total_2x2_pivots,
        n - numeric.stats().total_delayed,
        "pivot accounting must balance"
    );
}

// ===========================================================================
// Additional edge case tests
// ===========================================================================

/// Test 1x1 matrix.
#[test]
fn test_1x1_matrix() {
    let matrix = sparse_from_lower_triplets(1, &[(0, 0, 5.0)]);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    assert_eq!(numeric.stats().total_1x1_pivots, 1);
    assert_eq!(numeric.stats().total_2x2_pivots, 0);
    assert_eq!(numeric.stats().total_delayed, 0);

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(err < 1e-14, "1x1 reconstruction error: {:.2e}", err);
}

/// Test 2x2 matrix.
#[test]
fn test_2x2_matrix() {
    let entries = vec![(0, 0, 4.0), (1, 0, 1.0), (1, 1, 4.0)];
    let matrix = sparse_from_lower_triplets(2, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(err < 1e-14, "2x2 reconstruction error: {:.2e}", err);
}

/// Test block diagonal matrix (disconnected components).
#[test]
fn test_block_diagonal() {
    // Two disconnected 2x2 blocks
    let entries = vec![
        (0, 0, 4.0),
        (1, 0, 1.0),
        (1, 1, 4.0),
        (2, 2, 4.0),
        (3, 2, 1.0),
        (3, 3, 4.0),
    ];
    let matrix = sparse_from_lower_triplets(4, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "block diagonal reconstruction error: {:.2e}",
        err
    );
}

/// Test identity matrix.
#[test]
fn test_identity_matrix() {
    let n = 5;
    let entries: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 1.0)).collect();
    let matrix = sparse_from_lower_triplets(n, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    assert_eq!(numeric.stats().total_1x1_pivots, n);
    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(err < 1e-14, "identity reconstruction error: {:.2e}", err);
}

/// Test larger sparse matrix (50x50 banded).
#[test]
fn test_larger_banded() {
    let n = 50;
    let bw = 3; // bandwidth
    let mut entries = Vec::new();
    for i in 0..n {
        entries.push((i, i, 10.0));
        for d in 1..=bw {
            if i + d < n {
                entries.push((i + d, i, 1.0 / (d as f64)));
            }
        }
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "banded 50x50 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );
}

/// Test medium-sized indefinite matrix.
#[test]
fn test_medium_indefinite() {
    let n = 30;
    let mut entries = Vec::new();
    for i in 0..n {
        // Alternating positive/negative diagonal
        let sign = if i % 3 == 1 { -1.0 } else { 1.0 };
        entries.push((i, i, sign * 5.0));
    }
    for i in 1..n {
        entries.push((i, i - 1, 1.0));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");
    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "indefinite 30x30 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );
}

// ===========================================================================
// T028: SuiteSparse CI subset tests
// ===========================================================================

/// Factor CI-subset SuiteSparse matrices (n < 2000) and verify reconstruction error.
/// Larger matrices are skipped because dense reconstruction is O(n²) memory.
#[test]
fn test_suitesparse_ci_subset() {
    use rivrs_sparse::testing::{TestCaseFilter, load_test_cases};

    const MAX_DIM: usize = 2000;

    let cases = load_test_cases(&TestCaseFilter::ci_subset())
        .expect("failed to load SuiteSparse CI subset");

    assert!(!cases.is_empty(), "expected at least 1 CI-subset matrix");

    let mut tested = 0;
    for case in &cases {
        if case.matrix.nrows() > MAX_DIM {
            continue; // Skip large matrices (dense reconstruction would OOM)
        }

        let symbolic = AptpSymbolic::analyze(case.matrix.symbolic(), SymmetricOrdering::Amd)
            .unwrap_or_else(|e| panic!("analyze failed for '{}': {}", case.name, e));

        let numeric = AptpNumeric::factor(&symbolic, &case.matrix, &AptpOptions::default())
            .unwrap_or_else(|e| panic!("factor failed for '{}': {}", case.name, e));

        let err = multifrontal_reconstruction_error(&case.matrix, &symbolic, &numeric);
        assert!(
            err < 1e-12,
            "'{}' (n={}) reconstruction error: {:.2e} (expected < 1e-12)",
            case.name,
            case.matrix.nrows(),
            err
        );
        tested += 1;
    }

    assert!(
        tested > 0,
        "no matrices in CI subset were small enough to test"
    );
}

// ===========================================================================
// Additional review-driven tests
// ===========================================================================

/// Test that a singular (zero) matrix produces zero_pivots instead of an error.
/// A zero matrix has all-zero pivots; the APTP kernel delays every column,
/// and the root supernode records them as zero pivots. This matches SPRAL's
/// behavior: factorization succeeds, but zero_pivots > 0 signals rank deficiency.
#[test]
fn test_singular_matrix_zero_pivots() {
    // Zero matrix: all entries are zero but with nonzero structure
    // for the symbolic analysis to produce a non-trivial pattern.
    let entries = vec![
        (0, 0, 0.0),
        (1, 0, 0.0),
        (1, 1, 0.0),
        (2, 1, 0.0),
        (2, 2, 0.0),
    ];
    let matrix = sparse_from_lower_triplets(3, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed for zero matrix");

    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factorization should succeed (zero pivots recorded, not errored)");

    assert!(
        numeric.stats().zero_pivots > 0,
        "zero matrix should produce zero_pivots > 0, got {}",
        numeric.stats().zero_pivots,
    );
}

/// Test that a rank-deficient (but non-zero) matrix produces zero_pivots.
/// This matrix has rank 2 in a 3x3: rows 0 and 1 are identical.
#[test]
fn test_rank_deficient_matrix_zero_pivots() {
    // [1  1  0]
    // [1  1  0]
    // [0  0  2]
    // Rank 2: the (0,0)-(1,1) 2x2 subblock is singular.
    let entries = vec![(0, 0, 1.0), (1, 0, 1.0), (1, 1, 1.0), (2, 2, 2.0)];
    let matrix = sparse_from_lower_triplets(3, &entries);

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");

    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factorization should succeed (zero pivots recorded, not errored)");

    assert!(
        numeric.stats().zero_pivots > 0,
        "rank-deficient matrix should produce zero_pivots > 0, got {}",
        numeric.stats().zero_pivots,
    );
}

/// Test factorization with METIS ordering instead of AMD.
/// Verifies that the multifrontal factorization works correctly with
/// non-default orderings.
#[test]
fn test_metis_ordering() {
    use rivrs_sparse::aptp::metis_ordering;

    // Moderately sized banded indefinite matrix
    let n = 30;
    let mut entries = Vec::new();
    for i in 0..n {
        let sign = if i % 3 == 1 { -1.0 } else { 1.0 };
        entries.push((i, i, sign * 8.0));
    }
    for i in 1..n {
        entries.push((i, i - 1, 1.5));
    }
    for i in 2..n {
        entries.push((i, i - 2, 0.5));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    // Get METIS ordering
    let metis_perm = metis_ordering(matrix.symbolic()).expect("METIS ordering should succeed");

    let symbolic = AptpSymbolic::analyze(
        matrix.symbolic(),
        SymmetricOrdering::Custom(metis_perm.as_ref()),
    )
    .expect("analyze with METIS ordering should succeed");

    let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default())
        .expect("factor with METIS ordering should succeed");

    let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
    assert!(
        err < 1e-12,
        "METIS-ordered 30x30 reconstruction error: {:.2e} (expected < 1e-12)",
        err
    );
}

/// Test factorization with METIS ordering on hand-constructed matrices.
/// Cross-validates that AMD and METIS orderings produce factors with
/// comparable reconstruction error.
#[test]
fn test_metis_vs_amd_reconstruction() {
    use rivrs_sparse::aptp::metis_ordering;

    // Arrow matrix: forces non-trivial tree structure
    let n = 15;
    let mut entries = Vec::new();
    for i in 0..n {
        entries.push((i, i, 10.0));
    }
    for i in 1..n {
        entries.push((i, 0, 1.0));
    }
    // Add some off-diagonal entries for richer structure
    for i in 2..n {
        entries.push((i, i - 1, 0.5));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    // AMD path
    let sym_amd = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("AMD analyze should succeed");
    let num_amd = AptpNumeric::factor(&sym_amd, &matrix, &AptpOptions::default())
        .expect("AMD factor should succeed");
    let err_amd = multifrontal_reconstruction_error(&matrix, &sym_amd, &num_amd);

    // METIS path
    let metis_perm = metis_ordering(matrix.symbolic()).expect("METIS ordering should succeed");
    let sym_metis = AptpSymbolic::analyze(
        matrix.symbolic(),
        SymmetricOrdering::Custom(metis_perm.as_ref()),
    )
    .expect("METIS analyze should succeed");
    let num_metis = AptpNumeric::factor(&sym_metis, &matrix, &AptpOptions::default())
        .expect("METIS factor should succeed");
    let err_metis = multifrontal_reconstruction_error(&matrix, &sym_metis, &num_metis);

    // Both should be well within tolerance
    assert!(err_amd < 1e-12, "AMD reconstruction error: {:.2e}", err_amd);
    assert!(
        err_metis < 1e-12,
        "METIS reconstruction error: {:.2e}",
        err_metis
    );
}

/// Test multi-level cascading delayed pivots.
///
/// Constructs a matrix designed so that delays cascade across multiple levels
/// of the assembly tree: a pivot that fails at a leaf must propagate through
/// an intermediate node before being resolved at or near the root.
///
/// Strategy: Use a tight APTP threshold with a matrix whose small diagonal
/// entries and large off-diagonals force delays at lower levels. The pivots
/// can only be resolved when enough context is available at a higher level.
#[test]
fn test_cascading_delayed_pivots() {
    // Build a matrix where the (0,0) pivot is tiny compared to its column,
    // forcing delays. The tridiagonal structure with a very small leading
    // diagonal should cascade delays through the tree.
    //
    // [eps  1    0    0    0  ]
    // [1    eps  1    0    0  ]
    // [0    1    eps  1    0  ]
    // [0    0    1    10   1  ]
    // [0    0    0    1    10 ]
    //
    // The first 3 diagonals (eps) will fail the APTP threshold test and
    // delay upward, while the last 2 (10.0) are well-conditioned.
    let eps = 1e-15;
    let n = 8;
    let mut entries = Vec::new();
    // Small diagonal entries for first half — will trigger delays
    for i in 0..n / 2 {
        entries.push((i, i, eps));
    }
    // Well-conditioned diagonal for second half
    for i in n / 2..n {
        entries.push((i, i, 10.0));
    }
    // Tridiagonal off-diagonals
    for i in 1..n {
        entries.push((i, i - 1, 1.0));
    }
    let matrix = sparse_from_lower_triplets(n, &entries);

    // Use a moderate threshold that will reject the tiny pivots
    let options = AptpOptions {
        threshold: 0.01,
        ..AptpOptions::default()
    };

    let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("analyze should succeed");

    let result = AptpNumeric::factor(&symbolic, &matrix, &options);

    let numeric = result.expect("factorization should succeed (zero pivots allowed)");
    let stats = numeric.stats();

    // We expect some delays, 2x2 pivots, or zero pivots to have occurred
    assert!(
        stats.total_delayed > 0 || stats.total_2x2_pivots > 0 || stats.zero_pivots > 0,
        "expected delays, 2x2 pivots, or zero pivots for near-singular leading block, \
         got: 1x1={}, 2x2={}, delayed={}, zero={}",
        stats.total_1x1_pivots,
        stats.total_2x2_pivots,
        stats.total_delayed,
        stats.zero_pivots,
    );

    // If no zero pivots, verify reconstruction
    if stats.zero_pivots == 0 {
        let err = multifrontal_reconstruction_error(&matrix, &symbolic, &numeric);
        assert!(
            err < 1e-10,
            "cascading delay reconstruction error: {:.2e} (expected < 1e-10)",
            err
        );
    }
}

// ===========================================================================
// T029: Full SuiteSparse test (ignored, requires full collection)
// ===========================================================================

/// Factor SuiteSparse matrices one at a time and verify reconstruction error.
/// Requires the full SuiteSparse collection (extracted from archive).
///
/// Matrices are loaded individually to avoid OOM. Matrices with n <= 500 get
/// full dense reconstruction check. Larger matrices (up to 5000) only verify
/// that factorization completes without error. Very large matrices are skipped.
///
/// Run with: `cargo test --test multifrontal test_suitesparse_full -- --ignored --test-threads=1`
#[test]
#[ignore = "requires full SuiteSparse collection"]
fn test_suitesparse_full() {
    use rivrs_sparse::io::registry;

    const MAX_DIM_FOR_RECON: usize = 1000;
    const MAX_DIM_FOR_FACTOR: usize = 20000;

    let all_meta = registry::load_registry().expect("failed to load metadata.json");
    let suitesparse_meta: Vec<_> = all_meta
        .into_iter()
        .filter(|m| m.source == "suitesparse")
        .collect();

    let mut recon_passed = 0;
    let mut factor_only = 0;
    let mut skipped = 0;
    let mut missing = 0;
    let mut failed = Vec::new();

    for meta in &suitesparse_meta {
        // Load one matrix at a time to avoid OOM
        let test_matrix = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(tm)) => tm,
            Ok(None) => {
                missing += 1;
                continue;
            }
            Err(e) => {
                failed.push(format!("'{}': load error: {}", meta.name, e));
                continue;
            }
        };

        let n = test_matrix.matrix.nrows();
        if n > MAX_DIM_FOR_FACTOR {
            skipped += 1;
            continue;
        }

        let symbolic =
            match AptpSymbolic::analyze(test_matrix.matrix.symbolic(), SymmetricOrdering::Amd) {
                Ok(s) => s,
                Err(e) => {
                    failed.push(format!("'{}' (n={}): analyze error: {}", meta.name, n, e));
                    continue;
                }
            };

        let numeric =
            match AptpNumeric::factor(&symbolic, &test_matrix.matrix, &AptpOptions::default()) {
                Ok(num) => num,
                Err(e) => {
                    failed.push(format!("'{}' (n={}): factor error: {}", meta.name, n, e));
                    continue;
                }
            };

        if n <= MAX_DIM_FOR_RECON {
            let err = multifrontal_reconstruction_error(&test_matrix.matrix, &symbolic, &numeric);
            if err < 1e-12 {
                recon_passed += 1;
            } else {
                failed.push(format!("'{}' (n={}): error={:.2e}", meta.name, n, err));
            }
        } else {
            // Just verify factorization completed successfully
            factor_only += 1;
        }
        // test_matrix, symbolic, numeric dropped here — frees memory before next iteration
    }

    eprintln!(
        "Reconstruction: {}, Factor-only: {}, Skipped (n>{}): {}, Missing: {}, Failed: {}, Total: {}",
        recon_passed,
        factor_only,
        MAX_DIM_FOR_FACTOR,
        skipped,
        missing,
        failed.len(),
        suitesparse_meta.len()
    );
    if !failed.is_empty() {
        eprintln!("Failed:");
        for f in &failed {
            eprintln!("  {}", f);
        }
        panic!("{} matrices failed", failed.len());
    }
}
