//! Full SuiteSparse collection integration tests for AptpSymbolic analysis.
//!
//! These tests require the full 65-matrix SuiteSparse collection extracted at
//! `sparse/test-data/suitesparse/`. They are marked `#[ignore]` so that
//! `cargo test` skips them (CI-safe), but they can be run with:
//!
//!   cargo test --test symbolic_analysis_full -- --ignored --test-threads=1
//!
//! The `--test-threads=1` flag is **required** to avoid OOM when multiple tests
//! load large matrices concurrently.
//!
//! All tests use METIS ordering via `metis_ordering()` + `SymmetricOrdering::Custom`.
//! METIS handles all matrix sizes well with no dimension cap needed.
//!
//! # Extracting the full SuiteSparse collection
//!
//! The full collection is stored in `references/ssids/suitesparse.tar.gz`.
//! To extract into the expected location:
//!
//!   cd sparse/test-data
//!   tar xzf ../../references/ssids/suitesparse.tar.gz
//!
//! After extraction, `sparse/test-data/suitesparse/` should contain ~65 `.mtx`
//! files across `easy-indefinite/`, `hard-indefinite/`, and `positive-definite/`
//! subdirectories.
//!
//! # Memory considerations
//!
//! Matrices are loaded and analyzed one at a time to avoid OOM on large
//! collections. Each matrix is dropped before loading the next.
//!
//! # Unsupported / corrupt matrix files
//!
//! Matrices that fail to parse (e.g., truncated files, unsupported formats)
//! are silently skipped and reported in the summary.

use faer::sparse::linalg::cholesky::SymmetricOrdering;

use rivrs_sparse::io::registry;
use rivrs_sparse::ordering::metis_ordering;
use rivrs_sparse::symmetric::AptpSymbolic;

/// Minimum number of SuiteSparse matrices to consider the full collection present.
const MIN_FULL_COLLECTION_SIZE: usize = 20;

/// Load the SuiteSparse metadata entries (lightweight, no matrix data).
fn suitesparse_metadata() -> Vec<registry::MatrixMetadata> {
    let all_meta = registry::load_registry().expect("failed to load metadata.json");
    all_meta
        .into_iter()
        .filter(|m| m.source == "suitesparse")
        .collect()
}

/// Validate supernode column range partitioning and assembly tree validity.
///
/// For each supernodal matrix: verifies that supernode ranges partition [0, n),
/// assembly tree has parent > child (postorder) with at least one root, and
/// row patterns are accessible.
#[test]
#[ignore = "requires full SuiteSparse collection"]
fn test_supernodal_structure_full_suitesparse() {
    let entries = suitesparse_metadata();

    let mut loaded = 0usize;
    let mut supernodal_count = 0usize;

    for meta in &entries {
        let test_matrix = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(tm)) => tm,
            Ok(None) | Err(_) => continue,
        };
        loaded += 1;

        let perm = metis_ordering(test_matrix.matrix.symbolic())
            .unwrap_or_else(|e| panic!("METIS ordering failed for '{}': {}", meta.name, e));
        let sym = AptpSymbolic::analyze(
            test_matrix.matrix.symbolic(),
            SymmetricOrdering::Custom(perm.as_ref()),
        )
        .unwrap_or_else(|e| panic!("analysis should succeed for '{}': {}", meta.name, e));
        drop(test_matrix); // Free matrix memory before validation

        if !sym.is_supernodal() {
            continue;
        }
        supernodal_count += 1;

        let ns = sym.n_supernodes().unwrap();
        let begin = sym.supernode_begin().unwrap();
        let end = sym.supernode_end().unwrap();

        // Ranges partition [0, n)
        assert_eq!(
            begin[0], 0,
            "'{}' first supernode should start at 0",
            meta.name
        );
        assert_eq!(
            end[ns - 1],
            sym.nrows(),
            "'{}' last supernode should end at dimension",
            meta.name
        );
        for s in 0..ns - 1 {
            assert_eq!(
                end[s],
                begin[s + 1],
                "'{}' supernode {} end should equal supernode {} begin",
                meta.name,
                s,
                s + 1,
            );
        }
        for s in 0..ns {
            assert!(
                begin[s] < end[s],
                "'{}' supernode {} should have at least one column",
                meta.name,
                s,
            );
        }

        // Assembly tree: parent > child (postorder), at least one root
        let mut root_count = 0;
        for s in 0..ns {
            match sym.supernode_parent(s) {
                None => root_count += 1,
                Some(parent) => {
                    assert!(
                        parent > s,
                        "'{}' parent of supernode {} should be > {} (postorder), got {}",
                        meta.name,
                        s,
                        s,
                        parent,
                    );
                    assert!(
                        parent < ns,
                        "'{}' parent index {} out of range for {} supernodes",
                        meta.name,
                        parent,
                        ns,
                    );
                }
            }
        }
        assert!(
            root_count >= 1,
            "'{}' assembly tree should have at least one root",
            meta.name,
        );

        // Row patterns accessible
        for s in 0..ns {
            assert!(
                sym.supernode_pattern(s).is_some(),
                "'{}' supernode {} pattern should be accessible",
                meta.name,
                s,
            );
        }
    }

    if loaded < MIN_FULL_COLLECTION_SIZE {
        eprintln!(
            "Skipping: only {} SuiteSparse matrices loaded (need >= {}).",
            loaded, MIN_FULL_COLLECTION_SIZE,
        );
        return;
    }

    eprintln!(
        "Validated supernodal structure for {} out of {} loaded matrices",
        supernodal_count, loaded,
    );
}

/// Verify pivot buffer estimates are non-negative and proportional.
///
/// For each SuiteSparse matrix, checks that buffer estimates are non-negative,
/// the buffer/nnz ratio is within reasonable bounds (0.1% to 200%), and
/// buffer length matches the structural mode (n_supernodes for supernodal,
/// nrows for simplicial).
#[test]
#[ignore = "requires full SuiteSparse collection"]
fn test_pivot_buffer_sanity_full_suitesparse() {
    let entries = suitesparse_metadata();

    let mut loaded = 0usize;

    for meta in &entries {
        let test_matrix = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(tm)) => tm,
            Ok(None) | Err(_) => continue,
        };
        loaded += 1;

        let perm = metis_ordering(test_matrix.matrix.symbolic())
            .unwrap_or_else(|e| panic!("METIS ordering failed for '{}': {}", meta.name, e));
        let sym = AptpSymbolic::analyze(
            test_matrix.matrix.symbolic(),
            SymmetricOrdering::Custom(perm.as_ref()),
        )
        .unwrap_or_else(|e| panic!("analysis failed for '{}': {}", meta.name, e));
        drop(test_matrix); // Free matrix memory

        // Non-negative buffer values (usize is always >= 0, but check non-empty)
        assert!(
            !sym.pivot_buffer_estimates().is_empty(),
            "'{}' should have non-empty pivot buffer",
            meta.name,
        );

        // Buffer length matches structural mode
        if sym.is_supernodal() {
            assert_eq!(
                sym.pivot_buffer_estimates().len(),
                sym.n_supernodes().unwrap(),
                "'{}' supernodal buffer length should match n_supernodes",
                meta.name,
            );
        } else {
            assert_eq!(
                sym.pivot_buffer_estimates().len(),
                sym.nrows(),
                "'{}' simplicial buffer length should match nrows",
                meta.name,
            );
        }

        // Buffer/nnz ratio should be in reasonable bounds
        let total_buffer = sym.total_pivot_buffer();
        let predicted_nnz = sym.predicted_nnz();
        if predicted_nnz > 0 {
            let ratio = total_buffer as f64 / predicted_nnz as f64;
            assert!(
                (0.001..=2.0).contains(&ratio),
                "'{}' buffer/nnz ratio ({:.4}) should be between 0.1% and 200%",
                meta.name,
                ratio,
            );
        }
    }

    if loaded < MIN_FULL_COLLECTION_SIZE {
        eprintln!(
            "Skipping: only {} SuiteSparse matrices loaded (need >= {}).",
            loaded, MIN_FULL_COLLECTION_SIZE,
        );
    }
}

/// Analyze all SuiteSparse matrices with METIS ordering.
///
/// Verifies that every matrix in the full collection completes symbolic analysis
/// successfully with METIS ordering. No dimension cap needed.
#[test]
#[ignore = "requires full SuiteSparse collection"]
fn test_analyze_full_suitesparse_metis() {
    let entries = suitesparse_metadata();
    let start = std::time::Instant::now();

    let mut loaded = 0usize;
    let mut success_count = 0usize;
    let mut fail_count = 0usize;
    let mut skipped_parse = Vec::new();

    eprintln!(
        "\n{:<30} {:>8} {:>14} {:>10} {:>12}",
        "Matrix", "Dim", "Pred NNZ", "Mode", "Supernodes"
    );
    eprintln!("{}", "-".repeat(80));

    for meta in &entries {
        let test_matrix = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(tm)) => tm,
            Ok(None) => continue,
            Err(e) => {
                skipped_parse.push(format!("{}: {}", meta.name, e));
                continue;
            }
        };
        loaded += 1;

        // METIS ordering
        let perm = match metis_ordering(test_matrix.matrix.symbolic()) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("  METIS FAIL: {}: {}", meta.name, e);
                fail_count += 1;
                continue;
            }
        };

        match AptpSymbolic::analyze(
            test_matrix.matrix.symbolic(),
            SymmetricOrdering::Custom(perm.as_ref()),
        ) {
            Ok(sym) => {
                assert_eq!(
                    sym.nrows(),
                    meta.size,
                    "dimension mismatch for '{}'",
                    meta.name
                );
                assert!(
                    sym.predicted_nnz() > 0,
                    "predicted_nnz should be > 0 for '{}'",
                    meta.name
                );
                assert_eq!(
                    sym.etree().len(),
                    meta.size,
                    "etree length mismatch for '{}'",
                    meta.name
                );

                let mode = if sym.is_supernodal() {
                    "supernodal"
                } else {
                    "simplicial"
                };
                let sn = sym
                    .n_supernodes()
                    .map(|n| n.to_string())
                    .unwrap_or_else(|| "-".to_string());

                eprintln!(
                    "{:<30} {:>8} {:>14} {:>10} {:>12}",
                    meta.name,
                    sym.nrows(),
                    sym.predicted_nnz(),
                    mode,
                    sn,
                );

                success_count += 1;
            }
            Err(e) => {
                eprintln!("  ANALYZE FAIL: {}: {}", meta.name, e);
                fail_count += 1;
            }
        }
    }

    let elapsed = start.elapsed();

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
        "\nSummary: {} passed, {} failed, {} parse-skipped out of {} loaded ({:.1}s)",
        success_count,
        fail_count,
        skipped_parse.len(),
        loaded,
        elapsed.as_secs_f64(),
    );

    // Wall-clock timing: full collection with METIS should complete in reasonable time.
    // Target: < 2 minutes on production hardware. Allow 5 minutes in test environments
    // (Docker containers, CI) where CPUs may be slower or throttled.
    assert!(
        elapsed.as_secs() < 300,
        "Full SuiteSparse METIS analysis took {:?}, expected < 300s",
        elapsed
    );

    assert_eq!(
        fail_count, 0,
        "all loadable SuiteSparse matrices should analyze successfully with METIS"
    );
}
