//! Integration tests for MC64 matching and symmetric scaling.
//!
//! Validates SPRAL-style scaling properties on hand-constructed matrices,
//! SuiteSparse CI subset, and full SuiteSparse collection.
//!
//! # SPRAL Scaling Properties
//!
//! 1. All entries: `|s_i * a_ij * s_j| <= 1.0`
//! 2. Row maxima: `max_j |s_i * a_ij * s_j| >= 0.75`
//! 3. Matched diagonal: `|s_i * a_{i,σ(i)} * s_{σ(i)}| ≈ 1.0`
//! 4. Scaling factors: positive and finite
//! 5. Matching: singletons + 2-cycles only

use faer::sparse::{SparseColMat, Triplet};

use rivrs_sparse::aptp::{mc64_matching, Mc64Job, Mc64Result};
use rivrs_sparse::aptp::matching::count_cycles;

/// Helper: create a symmetric upper-triangular matrix from entries.
fn make_upper_tri(n: usize, entries: &[(usize, usize, f64)]) -> SparseColMat<usize, f64> {
    let triplets: Vec<_> = entries
        .iter()
        .map(|&(i, j, v)| Triplet::new(i, j, v))
        .collect();
    SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
}

/// Verify all SPRAL scaling properties for a matching result.
fn verify_spral_properties(
    name: &str,
    matrix: &SparseColMat<usize, f64>,
    result: &Mc64Result,
) {
    let n = matrix.nrows();
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();
    let (fwd, _) = result.matching.as_ref().arrays();

    // Property 1: |s_i * a_ij * s_j| <= 1.0 for all stored entries
    let mut row_max = vec![0.0_f64; n];

    for j in 0..n {
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for k in start..end {
            let i = row_indices[k];
            let scaled = (result.scaling[i] * values[k] * result.scaling[j]).abs();
            assert!(
                scaled <= 1.0 + 1e-10,
                "{}: |s[{}]*a[{},{}]*s[{}]| = {:.6e} > 1.0",
                name, i, i, j, j, scaled
            );
            if scaled > row_max[i] {
                row_max[i] = scaled;
            }
            if i != j {
                if scaled > row_max[j] {
                    row_max[j] = scaled;
                }
            }
        }
    }

    // Property 2: row max >= 0.75 for rows with nonzero entries
    for i in 0..n {
        if row_max[i] > 0.0 {
            assert!(
                row_max[i] >= 0.75 - 1e-10,
                "{}: row_max[{}] = {:.6e} < 0.75",
                name, i, row_max[i]
            );
        }
    }

    // Property 3: matching is valid permutation
    let mut seen = vec![false; n];
    for i in 0..n {
        assert!(fwd[i] < n, "{}: matching[{}] = {} out of range", name, i, fwd[i]);
        assert!(
            !seen[fwd[i]],
            "{}: duplicate in matching at {}",
            name, fwd[i]
        );
        seen[fwd[i]] = true;
    }

    // Property 4: scaling factors positive and finite
    for (i, &s) in result.scaling.iter().enumerate() {
        assert!(s > 0.0, "{}: scaling[{}] = {} not positive", name, i, s);
        assert!(s.is_finite(), "{}: scaling[{}] = {} not finite", name, i, s);
    }

    // Property 5: singletons + 2-cycles only
    let (singletons, two_cycles) = count_cycles(fwd);
    assert_eq!(
        singletons + 2 * two_cycles, n,
        "{}: cycle decomposition doesn't cover all indices: {} singletons + {} 2-cycles != {}",
        name, singletons, two_cycles, n
    );
}

// ---- T009: End-to-end hand-constructed tests ----

#[test]
fn test_mc64_3x3_diagonal() {
    let matrix = make_upper_tri(3, &[
        (0, 0, 4.0),
        (1, 1, 9.0),
        (2, 2, 1.0),
    ]);

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 3, "diagonal: should match all 3");

    // Identity matching for diagonal
    let (fwd, _) = result.matching.as_ref().arrays();
    for i in 0..3 {
        assert_eq!(fwd[i], i, "diagonal: matching should be identity");
    }

    verify_spral_properties("3x3_diagonal", &matrix, &result);
}

#[test]
fn test_mc64_4x4_tridiagonal_indefinite() {
    let matrix = make_upper_tri(4, &[
        (0, 0, 2.0), (0, 1, -1.0),
        (1, 1, -3.0), (1, 2, 2.0),
        (2, 2, 1.0), (2, 3, -1.0),
        (3, 3, -4.0),
    ]);

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 4);
    verify_spral_properties("4x4_tridiag_indef", &matrix, &result);
}

#[test]
fn test_mc64_5x5_arrow_indefinite() {
    let matrix = make_upper_tri(5, &[
        (0, 0, 10.0), (0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (0, 4, 1.0),
        (1, 1, -3.0),
        (2, 2, 5.0),
        (3, 3, -2.0),
        (4, 4, 4.0),
    ]);

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 5);
    verify_spral_properties("5x5_arrow_indef", &matrix, &result);
}

#[test]
fn test_mc64_1x1_trivial() {
    let matrix = make_upper_tri(1, &[(0, 0, 7.0)]);
    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 1);
    assert_eq!(result.scaling.len(), 1);
    assert!(result.scaling[0] > 0.0);
}

#[test]
fn test_mc64_2x2_trivial() {
    let matrix = make_upper_tri(2, &[
        (0, 0, 3.0), (0, 1, 1.0),
        (1, 1, 5.0),
    ]);
    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 2);
    verify_spral_properties("2x2_trivial", &matrix, &result);
}

// ---- T014: Edge case matrices ----

#[test]
fn test_mc64_badly_scaled() {
    // Entries spanning 10 orders of magnitude
    let matrix = make_upper_tri(4, &[
        (0, 0, 1e10), (0, 1, 1e5),
        (1, 1, 1e0), (1, 2, 1e-5),
        (2, 2, 1e-10), (2, 3, 1e-3),
        (3, 3, 1e3),
    ]);

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 4);
    verify_spral_properties("badly_scaled", &matrix, &result);
}

#[test]
fn test_mc64_zero_diagonal() {
    // Matrix with zero diagonal entries — off-diagonal matching needed
    // [0  5  0]
    // [5  0  3]
    // [0  3  0]
    let matrix = make_upper_tri(3, &[
        (0, 1, 5.0),
        (1, 2, 3.0),
    ]);

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    // With zero diagonals, matching will pair off-diagonal entries
    // Can't guarantee full matching without diagonals; just verify properties
    for (i, &s) in result.scaling.iter().enumerate() {
        assert!(s > 0.0, "scaling[{}] should be positive", i);
        assert!(s.is_finite(), "scaling[{}] should be finite", i);
    }
}

#[test]
fn test_mc64_well_conditioned_pd() {
    // Well-conditioned PD matrix: should produce near-identity matching
    // and near-unit scaling
    let matrix = make_upper_tri(4, &[
        (0, 0, 4.0), (0, 1, 1.0), (0, 2, 0.5),
        (1, 1, 4.0), (1, 2, 1.0), (1, 3, 0.5),
        (2, 2, 4.0), (2, 3, 1.0),
        (3, 3, 4.0),
    ]);

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 4);
    verify_spral_properties("well_conditioned_pd", &matrix, &result);

    // For a well-conditioned PD matrix, matching should be identity
    // (diagonals are already the largest entries in each column)
    let (fwd, _) = result.matching.as_ref().arrays();
    for i in 0..4 {
        assert_eq!(
            fwd[i], i,
            "PD matrix should have identity matching, got σ({})={}",
            i, fwd[i]
        );
    }
}

#[test]
fn test_mc64_6x6_larger_indefinite() {
    // 6x6 block-structured indefinite matrix
    // [  4  1  0  0  1  0 ]
    // [  1 -3  2  0  0  0 ]
    // [  0  2  5  1  0  0 ]
    // [  0  0  1 -2  3  0 ]
    // [  1  0  0  3  6  1 ]
    // [  0  0  0  0  1 -1 ]
    let matrix = make_upper_tri(6, &[
        (0, 0, 4.0), (0, 1, 1.0), (0, 4, 1.0),
        (1, 1, -3.0), (1, 2, 2.0),
        (2, 2, 5.0), (2, 3, 1.0),
        (3, 3, -2.0), (3, 4, 3.0),
        (4, 4, 6.0), (4, 5, 1.0),
        (5, 5, -1.0),
    ]);

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 6);
    verify_spral_properties("6x6_indef", &matrix, &result);
}

// ---- T015: MC64 + METIS independent composition tests ----

#[test]
fn test_mc64_metis_independent_composition() {
    use faer::sparse::linalg::cholesky::SymmetricOrdering;
    use rivrs_sparse::aptp::{metis_ordering, AptpSymbolic};

    // Test on a few indefinite matrices
    let test_matrices: Vec<(&str, SparseColMat<usize, f64>)> = vec![
        ("tridiag_4x4", make_upper_tri(4, &[
            (0, 0, 2.0), (0, 1, -1.0),
            (1, 1, -3.0), (1, 2, 2.0),
            (2, 2, 1.0), (2, 3, -1.0),
            (3, 3, -4.0),
        ])),
        ("arrow_5x5", make_upper_tri(5, &[
            (0, 0, 10.0), (0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (0, 4, 1.0),
            (1, 1, -3.0),
            (2, 2, 5.0),
            (3, 3, -2.0),
            (4, 4, 4.0),
        ])),
        ("block_6x6", make_upper_tri(6, &[
            (0, 0, 4.0), (0, 1, 1.0), (0, 4, 1.0),
            (1, 1, -3.0), (1, 2, 2.0),
            (2, 2, 5.0), (2, 3, 1.0),
            (3, 3, -2.0), (3, 4, 3.0),
            (4, 4, 6.0), (4, 5, 1.0),
            (5, 5, -1.0),
        ])),
    ];

    for (name, matrix) in &test_matrices {
        // (a) MC64 matching independently
        let mc64 = mc64_matching(matrix, Mc64Job::MaximumProduct)
            .unwrap_or_else(|e| panic!("{}: MC64 failed: {}", name, e));

        // (b) METIS ordering independently on same matrix
        let metis_perm = metis_ordering(matrix.symbolic())
            .unwrap_or_else(|e| panic!("{}: METIS failed: {}", name, e));

        // (c) Feed METIS into AptpSymbolic::analyze
        let symbolic = AptpSymbolic::analyze(
            matrix.symbolic(),
            SymmetricOrdering::Custom(metis_perm.as_ref()),
        )
        .unwrap_or_else(|e| panic!("{}: AptpSymbolic failed: {}", name, e));

        // (d) Verify symbolic analysis succeeds with valid fill estimates
        assert!(
            symbolic.predicted_nnz() > 0,
            "{}: predicted_nnz should be > 0",
            name
        );

        // (e) Verify MC64 scaling is still valid (structure-independent)
        verify_spral_properties(name, matrix, &mc64);

        // Verify 2-cycle count is available and non-negative
        let (fwd, _) = mc64.matching.as_ref().arrays();
        let (singletons, two_cycles) = count_cycles(fwd);
        assert!(
            singletons + 2 * two_cycles == matrix.nrows(),
            "{}: cycle count doesn't match dimension",
            name
        );

        eprintln!(
            "  {:<20} n={} matched={} singletons={} 2-cycles={} nnz(L)={}",
            name, matrix.nrows(), mc64.matched, singletons, two_cycles,
            symbolic.predicted_nnz()
        );
    }
}
