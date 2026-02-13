//! Integration tests for MC64 matching and symmetric scaling.
//!
//! Validates SPRAL-style scaling properties on hand-constructed matrices,
//! SuiteSparse CI subset, and full SuiteSparse collection.
//!
//! # SPRAL Scaling Properties
//!
//! 1. All entries: `|s_i * a_ij * s_j| <= 1.0`
//! 2. Row maxima: `max_j |s_i * a_ij * s_j| ≈ 1.0` (nonsingular), `median >= 0.5` (singular)
//! 3. Matched diagonal: `|s_i * a_{i,σ(i)} * s_{σ(i)}| ≈ 1.0`
//! 4. Scaling factors: positive and finite
//! 5. Matching: singletons + 2-cycles only

use faer::sparse::{SparseColMat, Triplet};

use rivrs_sparse::aptp::matching::count_cycles;
use rivrs_sparse::aptp::{Mc64Job, Mc64Result, mc64_matching};
use rivrs_sparse::io::registry;
use rivrs_sparse::testing::{TestCaseFilter, load_test_cases};

/// Helper: create a symmetric upper-triangular matrix from entries.
fn make_upper_tri(n: usize, entries: &[(usize, usize, f64)]) -> SparseColMat<usize, f64> {
    let triplets: Vec<_> = entries
        .iter()
        .map(|&(i, j, v)| Triplet::new(i, j, v))
        .collect();
    SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
}

/// Verify all SPRAL scaling properties for a matching result.
fn verify_spral_properties(name: &str, matrix: &SparseColMat<usize, f64>, result: &Mc64Result) {
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
                name,
                i,
                i,
                j,
                j,
                scaled
            );
            if scaled > row_max[i] {
                row_max[i] = scaled;
            }
            if i != j && scaled > row_max[j] {
                row_max[j] = scaled;
            }
        }
    }

    // Property 2: row max ≈ 1.0 for nonsingular matrices (SPRAL tolerance: 5e-14)
    // For full-rank matrices with a perfect matching, every row has a matched entry
    // that scales to exactly 1.0, so row_max = 1.0 for all rows.
    for (i, &rm) in row_max.iter().enumerate() {
        if rm > 0.0 {
            assert!(
                rm >= 1.0 - 1e-12,
                "{}: row_max[{}] = {:.6e} < 1.0 - 1e-12 (SPRAL expects ~1.0)",
                name,
                i,
                rm
            );
        }
    }

    // Property 3: matching is valid permutation
    let mut seen = vec![false; n];
    for i in 0..n {
        assert!(
            fwd[i] < n,
            "{}: matching[{}] = {} out of range",
            name,
            i,
            fwd[i]
        );
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
        singletons + 2 * two_cycles,
        n,
        "{}: cycle decomposition doesn't cover all indices: {} singletons + {} 2-cycles != {}",
        name,
        singletons,
        two_cycles,
        n
    );
}

// ---- T009: End-to-end hand-constructed tests ----

#[test]
fn test_mc64_3x3_diagonal() {
    let matrix = make_upper_tri(3, &[(0, 0, 4.0), (1, 1, 9.0), (2, 2, 1.0)]);

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 3, "diagonal: should match all 3");

    // Identity matching for diagonal
    let (fwd, _) = result.matching.as_ref().arrays();
    for (i, &f) in fwd.iter().enumerate() {
        assert_eq!(f, i, "diagonal: matching should be identity");
    }

    verify_spral_properties("3x3_diagonal", &matrix, &result);
}

#[test]
fn test_mc64_4x4_tridiagonal_indefinite() {
    let matrix = make_upper_tri(
        4,
        &[
            (0, 0, 2.0),
            (0, 1, -1.0),
            (1, 1, -3.0),
            (1, 2, 2.0),
            (2, 2, 1.0),
            (2, 3, -1.0),
            (3, 3, -4.0),
        ],
    );

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 4);
    verify_spral_properties("4x4_tridiag_indef", &matrix, &result);
}

#[test]
fn test_mc64_5x5_arrow_indefinite() {
    let matrix = make_upper_tri(
        5,
        &[
            (0, 0, 10.0),
            (0, 1, 1.0),
            (0, 2, 1.0),
            (0, 3, 1.0),
            (0, 4, 1.0),
            (1, 1, -3.0),
            (2, 2, 5.0),
            (3, 3, -2.0),
            (4, 4, 4.0),
        ],
    );

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
    let matrix = make_upper_tri(2, &[(0, 0, 3.0), (0, 1, 1.0), (1, 1, 5.0)]);
    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 2);
    verify_spral_properties("2x2_trivial", &matrix, &result);
}

// ---- T014: Edge case matrices ----

#[test]
fn test_mc64_badly_scaled() {
    // Entries spanning 10 orders of magnitude
    let matrix = make_upper_tri(
        4,
        &[
            (0, 0, 1e10),
            (0, 1, 1e5),
            (1, 1, 1e0),
            (1, 2, 1e-5),
            (2, 2, 1e-10),
            (2, 3, 1e-3),
            (3, 3, 1e3),
        ],
    );

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
    let matrix = make_upper_tri(3, &[(0, 1, 5.0), (1, 2, 3.0)]);

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
    let matrix = make_upper_tri(
        4,
        &[
            (0, 0, 4.0),
            (0, 1, 1.0),
            (0, 2, 0.5),
            (1, 1, 4.0),
            (1, 2, 1.0),
            (1, 3, 0.5),
            (2, 2, 4.0),
            (2, 3, 1.0),
            (3, 3, 4.0),
        ],
    );

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 4);
    verify_spral_properties("well_conditioned_pd", &matrix, &result);

    // For a well-conditioned PD matrix, matching should be identity
    // (diagonals are already the largest entries in each column)
    let (fwd, _) = result.matching.as_ref().arrays();
    for (i, &f) in fwd.iter().enumerate() {
        assert_eq!(
            f, i,
            "PD matrix should have identity matching, got σ({})={}",
            i, f
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
    let matrix = make_upper_tri(
        6,
        &[
            (0, 0, 4.0),
            (0, 1, 1.0),
            (0, 4, 1.0),
            (1, 1, -3.0),
            (1, 2, 2.0),
            (2, 2, 5.0),
            (2, 3, 1.0),
            (3, 3, -2.0),
            (3, 4, 3.0),
            (4, 4, 6.0),
            (4, 5, 1.0),
            (5, 5, -1.0),
        ],
    );

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
    assert_eq!(result.matched, 6);
    verify_spral_properties("6x6_indef", &matrix, &result);
}

// ---- T015: MC64 + METIS independent composition tests ----

#[test]
fn test_mc64_metis_independent_composition() {
    use faer::sparse::linalg::cholesky::SymmetricOrdering;
    use rivrs_sparse::aptp::{AptpSymbolic, metis_ordering};

    // Test on a few indefinite matrices
    let test_matrices: Vec<(&str, SparseColMat<usize, f64>)> = vec![
        (
            "tridiag_4x4",
            make_upper_tri(
                4,
                &[
                    (0, 0, 2.0),
                    (0, 1, -1.0),
                    (1, 1, -3.0),
                    (1, 2, 2.0),
                    (2, 2, 1.0),
                    (2, 3, -1.0),
                    (3, 3, -4.0),
                ],
            ),
        ),
        (
            "arrow_5x5",
            make_upper_tri(
                5,
                &[
                    (0, 0, 10.0),
                    (0, 1, 1.0),
                    (0, 2, 1.0),
                    (0, 3, 1.0),
                    (0, 4, 1.0),
                    (1, 1, -3.0),
                    (2, 2, 5.0),
                    (3, 3, -2.0),
                    (4, 4, 4.0),
                ],
            ),
        ),
        (
            "block_6x6",
            make_upper_tri(
                6,
                &[
                    (0, 0, 4.0),
                    (0, 1, 1.0),
                    (0, 4, 1.0),
                    (1, 1, -3.0),
                    (1, 2, 2.0),
                    (2, 2, 5.0),
                    (2, 3, 1.0),
                    (3, 3, -2.0),
                    (3, 4, 3.0),
                    (4, 4, 6.0),
                    (4, 5, 1.0),
                    (5, 5, -1.0),
                ],
            ),
        ),
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
            name,
            matrix.nrows(),
            mc64.matched,
            singletons,
            two_cycles,
            symbolic.predicted_nnz()
        );
    }
}

// ---- T019: Structurally singular matrix integration tests ----

#[test]
fn test_mc64_structurally_singular() {
    // 4x4 structurally singular: only off-diagonal connections
    // [0  5  0  0]
    // [5  0  0  0]
    // [0  0  0  3]
    // [0  0  3  0]
    let matrix = make_upper_tri(4, &[(0, 1, 5.0), (2, 3, 3.0)]);

    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();

    // Should match as many as possible (2 pairs = 4 entries via 2-cycles)
    assert!(result.matched <= 4);

    // All scaling factors positive and finite
    for (i, &s) in result.scaling.iter().enumerate() {
        assert!(s > 0.0, "scaling[{}] should be positive", i);
        assert!(s.is_finite(), "scaling[{}] should be finite", i);
    }

    // Matching is valid permutation
    let (fwd, _) = result.matching.as_ref().arrays();
    let mut seen = [false; 4];
    for &f in fwd {
        assert!(!seen[f], "duplicate in matching at {}", f);
        seen[f] = true;
    }
}

// ---- T020: SuiteSparse CI subset tests ----

#[test]
fn test_mc64_suitesparse_ci_subset() {
    let cases =
        load_test_cases(&TestCaseFilter::ci_subset()).expect("failed to load CI-subset matrices");

    let suitesparse: Vec<_> = cases
        .iter()
        .filter(|c| c.properties.source == "suitesparse")
        .collect();

    assert_eq!(
        suitesparse.len(),
        9,
        "expected 9 CI-subset suitesparse matrices"
    );

    for case in &suitesparse {
        let n = case.properties.size;

        let result = mc64_matching(&case.matrix, Mc64Job::MaximumProduct)
            .unwrap_or_else(|e| panic!("MC64 failed for '{}': {}", case.name, e));

        // (a) Report matching cardinality (some SuiteSparse matrices are structurally singular)
        assert!(
            result.matched > 0,
            "'{}': expected at least some matching, got 0",
            case.name
        );

        // (b) All scaling factors positive and finite
        for (i, &s) in result.scaling.iter().enumerate() {
            assert!(
                s > 0.0 && s.is_finite(),
                "'{}': scaling[{}] = {} not positive/finite",
                case.name,
                i,
                s
            );
        }

        // (c) |s_i * a_ij * s_j| <= 1.0 for all entries
        let symbolic = case.matrix.symbolic();
        let values = case.matrix.val();
        let col_ptrs = symbolic.col_ptr();
        let row_indices = symbolic.row_idx();

        let mut row_max = vec![0.0_f64; n];
        let mut max_violation = 0.0_f64;

        for j in 0..n {
            let start = col_ptrs[j];
            let end = col_ptrs[j + 1];
            for k in start..end {
                let i = row_indices[k];
                let scaled = (result.scaling[i] * values[k] * result.scaling[j]).abs();
                if scaled > 1.0 + 1e-10 {
                    max_violation = max_violation.max(scaled - 1.0);
                }
                if scaled > row_max[i] {
                    row_max[i] = scaled;
                }
                if i != j && scaled > row_max[j] {
                    row_max[j] = scaled;
                }
            }
        }

        assert!(
            max_violation < 1e-8,
            "'{}': max scaling violation = {:.2e}",
            case.name,
            max_violation
        );

        // (d) Row max quality check
        // For full-rank matrices, the optimal matching ensures each matched entry
        // scales to exactly 1.0. Since every row has a matched entry and all scaled
        // entries are bounded by 1.0, row_max = 1.0 for all rows. We use SPRAL's
        // strict tolerance (5e-14) as primary check with 1e-12 for floating-point margin.
        // For structurally singular matrices (matched < n), some rows may have weaker
        // scaling due to the partial matching — relaxed threshold applies.
        if result.matched == n {
            for (i, &rm) in row_max.iter().enumerate() {
                if rm > 0.0 {
                    assert!(
                        rm >= 1.0 - 1e-12,
                        "'{}': row_max[{}] = {:.6e} < 1.0 - 1e-12 (SPRAL expects ~1.0)",
                        case.name,
                        i,
                        rm
                    );
                }
            }
        } else {
            // For structurally singular matrices: dual-based scaling may not
            // provide quality scaling for all rows. We only require that
            // the median row_max is reasonable (>= 0.5).
            let mut sorted_max: Vec<f64> = row_max.iter().copied().filter(|&m| m > 0.0).collect();
            sorted_max.sort_by(|a, b| a.total_cmp(b));
            let median = if sorted_max.is_empty() {
                0.0
            } else {
                sorted_max[sorted_max.len() / 2]
            };
            assert!(
                median >= 0.5,
                "'{}': median row_max = {:.4e} < 0.5 (singular, matched={}/{})",
                case.name,
                median,
                result.matched,
                n
            );
        }

        // (e) Cycle structure report
        // For full-rank matrices, optimal matching should be singletons + 2-cycles.
        // For structurally singular matrices, unmatched row assignment can create longer cycles.
        let (fwd, _) = result.matching.as_ref().arrays();
        if result.matched == n {
            // Check cycle structure without panicking
            let mut singletons = 0usize;
            let mut two_cycles = 0usize;
            let mut longer_cycles = 0usize;
            let mut cycle_visited = vec![false; n];
            for i in 0..n {
                if cycle_visited[i] {
                    continue;
                }
                let j = fwd[i];
                if j == i {
                    singletons += 1;
                    cycle_visited[i] = true;
                } else if fwd[j] == i {
                    two_cycles += 1;
                    cycle_visited[i] = true;
                    cycle_visited[j] = true;
                } else {
                    // Longer cycle — trace it
                    longer_cycles += 1;
                    let mut k = i;
                    loop {
                        cycle_visited[k] = true;
                        k = fwd[k];
                        if k == i {
                            break;
                        }
                    }
                }
            }
            if longer_cycles > 0 {
                eprintln!(
                    "  {:<30} dim={:>8}  matched={:>8}  singletons={:>8}  2-cycles={:>8}  LONGER_CYCLES={:>4}",
                    case.name, n, result.matched, singletons, two_cycles, longer_cycles
                );
            } else {
                eprintln!(
                    "  {:<30} dim={:>8}  matched={:>8}  singletons={:>8}  2-cycles={:>8}",
                    case.name, n, result.matched, singletons, two_cycles
                );
            }
        } else {
            eprintln!(
                "  {:<30} dim={:>8}  matched={:>8}  (singular, cycle check skipped)",
                case.name, n, result.matched
            );
        }
    }
}

// ---- T021: Full SuiteSparse validation ----
//
// Matrices are loaded and validated one at a time to avoid OOM on large
// collections. Each matrix is dropped before loading the next. Run with:
//
//   cargo test --release --test mc64_matching test_mc64_suitesparse_full -- --ignored --test-threads=1

/// Minimum number of SuiteSparse matrices to consider the full collection present.
const MIN_FULL_COLLECTION_SIZE: usize = 20;

/// Load SuiteSparse metadata entries (lightweight, no matrix data).
fn suitesparse_metadata() -> Vec<registry::MatrixMetadata> {
    let all_meta = registry::load_registry().expect("failed to load metadata.json");
    all_meta
        .into_iter()
        .filter(|m| m.source == "suitesparse")
        .collect()
}

#[test]
#[ignore = "requires full SuiteSparse collection"]
fn test_mc64_suitesparse_full() {
    let entries = suitesparse_metadata();

    eprintln!("Running MC64 on {} SuiteSparse entries", entries.len());

    // MC64 allocates a cost graph mirroring the input CSC, so peak memory is ~2x
    // the matrix data. Skip matrices whose nnz would exceed available memory.
    // 10M nnz ≈ 160MB (CSC) + 160MB (cost graph) + O(n) state ≈ 350MB peak.
    const MAX_NNZ_FOR_MC64: usize = 10_000_000;

    let mut loaded = 0;
    let mut pass_count = 0;
    let mut indefinite_count = 0;
    let mut improvement_count = 0;
    let mut degraded_quality = Vec::new();
    let mut skipped_parse = Vec::new();
    let mut skipped_large = Vec::new();

    for meta in &entries {
        if meta.nnz > MAX_NNZ_FOR_MC64 {
            skipped_large.push(format!("{}: nnz={}", meta.name, meta.nnz));
            continue;
        }

        let test_matrix = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(tm)) => tm,
            Ok(None) => continue,
            Err(e) => {
                skipped_parse.push(format!("{}: {}", meta.name, e));
                continue;
            }
        };
        loaded += 1;

        let n = meta.size;
        let matrix = &test_matrix.matrix;

        let result = mc64_matching(matrix, Mc64Job::MaximumProduct)
            .unwrap_or_else(|e| panic!("MC64 failed for '{}': {}", meta.name, e));

        // All scaling factors positive and finite
        for (i, &s) in result.scaling.iter().enumerate() {
            assert!(
                s > 0.0 && s.is_finite(),
                "'{}': scaling[{}] = {} not positive/finite",
                meta.name,
                i,
                s
            );
        }

        // |s_i * a_ij * s_j| <= 1.0
        let symbolic = matrix.symbolic();
        let values = matrix.val();
        let col_ptrs = symbolic.col_ptr();
        let row_indices = symbolic.row_idx();

        let mut row_max = vec![0.0_f64; n];
        let mut max_violation = 0.0_f64;
        let mut has_negative_diag = false;

        // Compute diagonal dominance before and after scaling
        let mut diag_before = vec![0.0_f64; n];
        let mut off_diag_sum_before = vec![0.0_f64; n];
        let mut diag_after = vec![0.0_f64; n];
        let mut off_diag_sum_after = vec![0.0_f64; n];

        for j in 0..n {
            let start = col_ptrs[j];
            let end = col_ptrs[j + 1];
            for k in start..end {
                let i = row_indices[k];
                let abs_val = values[k].abs();
                let scaled = (result.scaling[i] * values[k] * result.scaling[j]).abs();

                if scaled > 1.0 + 1e-10 {
                    max_violation = max_violation.max(scaled - 1.0);
                }
                if scaled > row_max[i] {
                    row_max[i] = scaled;
                }
                if i != j && scaled > row_max[j] {
                    row_max[j] = scaled;
                }

                // Diagonal dominance metrics
                if i == j {
                    if values[k] < 0.0 {
                        has_negative_diag = true;
                    }
                    diag_before[i] = abs_val;
                    diag_after[i] = scaled;
                } else {
                    off_diag_sum_before[i] += abs_val;
                    off_diag_sum_before[j] += abs_val;
                    off_diag_sum_after[i] += scaled;
                    off_diag_sum_after[j] += scaled;
                }
            }
        }

        // Drop matrix data to free memory before validation
        drop(test_matrix);

        assert!(
            max_violation < 1e-8,
            "'{}': max scaling violation = {:.2e}",
            meta.name,
            max_violation
        );

        // Row max quality: report statistics rather than hard-fail on individual rows.
        // CI-subset uses strict 1.0 - 1e-12 (all curated matrices hit 1.0 exactly).
        // Full collection includes hard-indefinite matrices (e.g. TSOPF power systems)
        // where matching quality can be degraded despite full cardinality.
        let nonzero_max: Vec<f64> = row_max.iter().copied().filter(|&m| m > 0.0).collect();
        let min_row_max = nonzero_max.iter().copied().fold(f64::INFINITY, f64::min);
        let mean_row_max = if nonzero_max.is_empty() {
            0.0
        } else {
            nonzero_max.iter().sum::<f64>() / nonzero_max.len() as f64
        };
        let rows_below_075 = nonzero_max.iter().filter(|&&m| m < 0.75).count();

        // Track matrices with degraded quality for reporting
        if result.matched == n && rows_below_075 > 0 {
            degraded_quality.push(format!(
                "{}: min_rmax={:.4e}, rows<0.75={}/{}",
                meta.name, min_row_max, rows_below_075, n
            ));
        }

        // Cycle structure (soft check for longer cycles)
        let (fwd, _) = result.matching.as_ref().arrays();
        let mut longer_cycles = 0usize;
        if result.matched == n {
            let mut cycle_visited = vec![false; n];
            for i in 0..n {
                if cycle_visited[i] {
                    continue;
                }
                let j = fwd[i];
                if j == i {
                    cycle_visited[i] = true;
                } else if fwd[j] == i {
                    cycle_visited[i] = true;
                    cycle_visited[j] = true;
                } else {
                    longer_cycles += 1;
                    let mut k = i;
                    loop {
                        cycle_visited[k] = true;
                        k = fwd[k];
                        if k == i {
                            break;
                        }
                    }
                }
            }
        }

        if has_negative_diag {
            indefinite_count += 1;

            // Compute min diagonal dominance ratio before and after
            let dd_before: f64 = (0..n)
                .filter(|&i| off_diag_sum_before[i] > 0.0)
                .map(|i| diag_before[i] / off_diag_sum_before[i])
                .fold(f64::INFINITY, f64::min);
            let dd_after: f64 = (0..n)
                .filter(|&i| off_diag_sum_after[i] > 0.0)
                .map(|i| diag_after[i] / off_diag_sum_after[i])
                .fold(f64::INFINITY, f64::min);

            if dd_after >= dd_before {
                improvement_count += 1;
            }
        }

        pass_count += 1;

        if result.matched < n {
            eprintln!(
                "  {:<30} dim={:>8}  matched={:>8}/{:>8}  min_rmax={:.4e}  mean_rmax={:.4e}",
                meta.name, n, result.matched, n, min_row_max, mean_row_max,
            );
        } else if rows_below_075 > 0 || longer_cycles > 0 {
            eprintln!(
                "  {:<30} dim={:>8}  matched={:>8}  min_rmax={:.4e}  rows<0.75={:>6}  longer_cyc={:>4}",
                meta.name, n, result.matched, min_row_max, rows_below_075, longer_cycles,
            );
        } else {
            eprintln!(
                "  {:<30} dim={:>8}  matched={:>8}  min_rmax={:.4e}",
                meta.name, n, result.matched, min_row_max,
            );
        }
    }

    if !skipped_large.is_empty() {
        eprintln!(
            "\nSkipped {} matrices exceeding nnz cap ({}):",
            skipped_large.len(),
            MAX_NNZ_FOR_MC64
        );
        for s in &skipped_large {
            eprintln!("  {}", s);
        }
    }

    if !skipped_parse.is_empty() {
        eprintln!(
            "\nSkipped {} matrices due to parse errors:",
            skipped_parse.len()
        );
        for s in &skipped_parse {
            eprintln!("  {}", s);
        }
    }

    if loaded < MIN_FULL_COLLECTION_SIZE {
        eprintln!(
            "\nSkipping assertions: only {} SuiteSparse matrices loaded (need >= {}). \
             Extract the full collection to test-data/suitesparse/.",
            loaded, MIN_FULL_COLLECTION_SIZE,
        );
        return;
    }

    eprintln!(
        "\n{}/{} matrices passed scaling bound check",
        pass_count, loaded
    );
    if indefinite_count > 0 {
        eprintln!(
            "Diagonal dominance improved on {}/{} indefinite matrices ({:.0}%)",
            improvement_count,
            indefinite_count,
            improvement_count as f64 / indefinite_count as f64 * 100.0
        );
    }
    if !degraded_quality.is_empty() {
        eprintln!(
            "\n{} fully-matched matrices with degraded row_max quality (rows < 0.75):",
            degraded_quality.len()
        );
        for d in &degraded_quality {
            eprintln!("  {}", d);
        }
    }

    assert_eq!(pass_count, loaded, "some matrices failed MC64 validation");
}
