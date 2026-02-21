//! Integration tests for AptpSymbolic analysis on the full test matrix suite.
//!
//! Tests symbolic analysis with AMD ordering across all hand-constructed and
//! SuiteSparse CI-subset matrices, verifying valid statistics and structural
//! invariants.

use faer::sparse::linalg::cholesky::SymmetricOrdering;

use rivrs_sparse::aptp::AptpSymbolic;
use rivrs_sparse::testing::{TestCaseFilter, load_test_cases};

#[test]
fn test_analyze_all_hand_constructed() {
    let cases = load_test_cases(&TestCaseFilter::hand_constructed())
        .expect("failed to load hand-constructed test cases");

    assert_eq!(cases.len(), 15, "expected 15 hand-constructed matrices");

    for case in &cases {
        let result = AptpSymbolic::analyze(case.matrix.symbolic(), SymmetricOrdering::Amd);
        assert!(
            result.is_ok(),
            "symbolic analysis failed for '{}': {:?}",
            case.name,
            result.err()
        );

        let sym = result.unwrap();

        // Dimension should match the test matrix
        assert_eq!(
            sym.nrows(),
            case.properties.size,
            "dimension mismatch for '{}'",
            case.name
        );

        // Predicted NNZ should be positive
        assert!(
            sym.predicted_nnz() > 0,
            "predicted_nnz should be > 0 for '{}'",
            case.name
        );

        // Statistics should be valid
        let stats = sym.statistics();
        assert_eq!(stats.dimension, case.properties.size);
        assert!(
            stats.average_col_count >= 1.0,
            "average_col_count should be >= 1.0 for '{}', got {}",
            case.name,
            stats.average_col_count
        );

        // Etree should have correct length
        assert_eq!(
            sym.etree().len(),
            case.properties.size,
            "etree length mismatch for '{}'",
            case.name
        );

        // Col counts should have correct length
        assert_eq!(
            sym.col_counts().len(),
            case.properties.size,
            "col_counts length mismatch for '{}'",
            case.name
        );
    }
}

#[test]
fn test_analyze_suitesparse_ci_subset() {
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
        let result = AptpSymbolic::analyze(case.matrix.symbolic(), SymmetricOrdering::Amd);
        assert!(
            result.is_ok(),
            "symbolic analysis failed for '{}': {:?}",
            case.name,
            result.err()
        );

        let sym = result.unwrap();

        // Dimension should match
        assert_eq!(
            sym.nrows(),
            case.properties.size,
            "dimension mismatch for '{}'",
            case.name
        );

        // Predicted NNZ should be positive
        assert!(
            sym.predicted_nnz() > 0,
            "predicted_nnz should be > 0 for '{}'",
            case.name
        );

        // Statistics should be valid
        let stats = sym.statistics();
        assert_eq!(stats.dimension, case.properties.size);
        assert!(
            stats.average_col_count >= 1.0,
            "average_col_count should be >= 1.0 for '{}'",
            case.name
        );

        // Pivot buffer should exist and be non-empty
        assert!(
            !sym.pivot_buffer_estimates().is_empty(),
            "pivot_buffer_estimates should not be empty for '{}'",
            case.name
        );
        assert!(
            sym.total_pivot_buffer() > 0,
            "total_pivot_buffer should be > 0 for '{}'",
            case.name
        );
    }
}

#[test]
fn test_custom_ordering_on_hand_constructed() {
    let cases = load_test_cases(&TestCaseFilter::hand_constructed())
        .expect("failed to load hand-constructed test cases");

    // Pick the first case with size > 1 as representative
    let case = cases
        .iter()
        .find(|c| c.properties.size > 1)
        .expect("need at least one matrix with size > 1");

    // AMD ordering
    let sym_amd = AptpSymbolic::analyze(case.matrix.symbolic(), SymmetricOrdering::Amd)
        .expect("AMD analysis should succeed");

    // Identity ordering
    let n = case.properties.size;
    let fwd: Vec<usize> = (0..n).collect();
    let perm = rivrs_sparse::aptp::perm_from_forward(fwd).unwrap();
    let sym_identity = AptpSymbolic::analyze(
        case.matrix.symbolic(),
        SymmetricOrdering::Custom(perm.as_ref()),
    )
    .expect("identity ordering analysis should succeed");

    // Both should produce valid results
    assert_eq!(sym_amd.nrows(), n);
    assert_eq!(sym_identity.nrows(), n);
    assert!(sym_amd.predicted_nnz() > 0);
    assert!(sym_identity.predicted_nnz() > 0);

    // Both should have valid etrees
    assert_eq!(sym_amd.etree().len(), n);
    assert_eq!(sym_identity.etree().len(), n);
}

#[test]
fn test_supernodal_structure_suitesparse() {
    let cases =
        load_test_cases(&TestCaseFilter::ci_subset()).expect("failed to load CI-subset matrices");

    let suitesparse: Vec<_> = cases
        .iter()
        .filter(|c| c.properties.source == "suitesparse")
        .collect();

    for case in &suitesparse {
        let sym = AptpSymbolic::analyze(case.matrix.symbolic(), SymmetricOrdering::Amd)
            .unwrap_or_else(|e| panic!("analysis should succeed for '{}': {}", case.name, e));

        if sym.is_supernodal() {
            let ns = sym.n_supernodes().unwrap();
            assert!(ns >= 1, "'{}' should have at least 1 supernode", case.name);

            // Supernode ranges should be valid
            let begin = sym.supernode_begin().unwrap();
            let end = sym.supernode_end().unwrap();

            assert_eq!(
                begin[0], 0,
                "'{}' first supernode should start at 0",
                case.name
            );
            assert_eq!(
                end[ns - 1],
                sym.nrows(),
                "'{}' last supernode should end at dimension",
                case.name
            );

            // Each supernode should have at least one column
            for s in 0..ns {
                assert!(
                    begin[s] < end[s],
                    "'{}' supernode {} should have at least one column",
                    case.name,
                    s
                );
            }

            // Assembly tree should be valid
            let mut root_count = 0;
            for s in 0..ns {
                match sym.supernode_parent(s) {
                    None => root_count += 1,
                    Some(parent) => {
                        assert!(
                            parent > s,
                            "'{}' parent of supernode {} should be > {} (postorder), got {}",
                            case.name,
                            s,
                            s,
                            parent
                        );
                    }
                }
            }
            assert!(
                root_count >= 1,
                "'{}' assembly tree should have at least one root",
                case.name
            );

            // Row patterns should be accessible
            for s in 0..ns {
                let pattern = sym.supernode_pattern(s);
                assert!(
                    pattern.is_some(),
                    "'{}' supernode {} pattern should be Some",
                    case.name,
                    s
                );
            }
        }
    }
}
