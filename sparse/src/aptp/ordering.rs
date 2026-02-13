//! METIS nested dissection ordering for symmetric sparse matrices.
//!
//! This module provides fill-reducing permutations computed by METIS's multilevel
//! nested dissection algorithm. The resulting permutations integrate with
//! [`AptpSymbolic::analyze()`](super::AptpSymbolic::analyze) via
//! [`SymmetricOrdering::Custom`](faer::sparse::linalg::cholesky::SymmetricOrdering::Custom).
//!
//! # Algorithm References
//!
//! - Karypis & Kumar (1998), "A Fast and High Quality Multilevel Scheme for
//!   Partitioning Irregular Graphs", SIAM J. Sci. Comput. 20(1)
//! - George (1973), "Nested Dissection of a Regular Finite Element Mesh",
//!   SIAM J. Numer. Anal. 10(2)
//!
//! # Implementation Notes
//!
//! Uses `metis-sys` for raw FFI bindings to METIS 5.x (vendored C source,
//! compiled via `cc`). The `METIS_NodeND` function computes nested dissection
//! orderings on undirected graphs extracted from the matrix sparsity pattern.

use faer::perm::Perm;
use faer::sparse::SymbolicSparseColMatRef;

use crate::error::SparseError;

/// Compute a METIS nested dissection fill-reducing ordering for a symmetric sparse matrix.
///
/// Extracts the adjacency graph from the matrix's sparsity pattern and calls METIS
/// to compute a fill-reducing permutation via multilevel nested dissection.
///
/// The returned permutation is suitable for use with
/// [`AptpSymbolic::analyze()`](super::AptpSymbolic::analyze) via
/// [`SymmetricOrdering::Custom`](faer::sparse::linalg::cholesky::SymmetricOrdering::Custom).
///
/// # Arguments
///
/// * `matrix` — Symbolic sparsity pattern of a symmetric sparse matrix.
///   Only the structural pattern is used; numerical values are ignored.
///   May contain upper triangle, lower triangle, or full symmetric structure.
///
/// # Returns
///
/// A fill-reducing permutation as `Perm<usize>`.
///
/// # Errors
///
/// Returns [`SparseError::InvalidInput`] if the matrix dimension exceeds `i32::MAX`.
/// Returns [`SparseError::AnalysisFailure`] if METIS reports an internal error.
///
/// # Algorithm References
///
/// - Karypis & Kumar (1998), "A Fast and High Quality Multilevel Scheme for
///   Partitioning Irregular Graphs", SIAM J. Sci. Comput.
/// - George (1973), "Nested Dissection of a Regular Finite Element Mesh"
///
/// # Examples
///
/// ```no_run
/// use rivrs_sparse::aptp::{AptpSymbolic, metis_ordering};
/// use faer::sparse::linalg::cholesky::SymmetricOrdering;
/// use faer::sparse::{SparseColMat, Triplet};
///
/// let triplets = vec![
///     Triplet::new(0, 0, 4.0),
///     Triplet::new(0, 1, 1.0),
///     Triplet::new(1, 0, 1.0),
///     Triplet::new(1, 1, 4.0),
/// ];
/// let matrix = SparseColMat::try_new_from_triplets(2, 2, &triplets).unwrap();
///
/// let perm = metis_ordering(matrix.symbolic()).expect("METIS ordering failed");
/// let symbolic = AptpSymbolic::analyze(
///     matrix.symbolic(),
///     SymmetricOrdering::Custom(perm.as_ref()),
/// ).expect("symbolic analysis failed");
/// ```
pub fn metis_ordering(
    matrix: SymbolicSparseColMatRef<'_, usize>,
) -> Result<Perm<usize>, SparseError> {
    let n = matrix.nrows();

    // Validate dimension fits in i32 (METIS uses idx_t = i32).
    if n > i32::MAX as usize {
        return Err(SparseError::InvalidInput {
            reason: format!(
                "Matrix dimension {} exceeds i32::MAX ({}); METIS requires 32-bit indices",
                n,
                i32::MAX
            ),
        });
    }

    // Trivial cases: dim 0 or 1 — return identity permutation without calling METIS.
    if n <= 1 {
        return Ok(identity_perm(n));
    }

    // Extract adjacency structure.
    let (mut xadj, mut adjncy) = extract_adjacency(matrix)?;

    // If no edges (diagonal matrix), return identity permutation.
    if adjncy.is_empty() {
        return Ok(identity_perm(n));
    }

    // Set up METIS options: 0-based numbering, defaults for everything else.
    let mut options = [0i32; metis_sys::METIS_NOPTIONS as usize];
    unsafe {
        metis_sys::METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    options[metis_sys::moptions_et_METIS_OPTION_NUMBERING as usize] = 0;

    // Allocate output arrays.
    let mut nvtxs = n as i32;
    let mut perm_i32 = vec![0i32; n]; // METIS perm: old → new
    let mut iperm_i32 = vec![0i32; n]; // METIS iperm: new → old

    // Call METIS_NodeND.
    let ret = unsafe {
        metis_sys::METIS_NodeND(
            &mut nvtxs,
            xadj.as_mut_ptr(),
            adjncy.as_mut_ptr(),
            std::ptr::null_mut(), // vwgt: no vertex weights
            options.as_mut_ptr(),
            perm_i32.as_mut_ptr(),
            iperm_i32.as_mut_ptr(),
        )
    };

    // Check return code.
    if ret != metis_sys::rstatus_et_METIS_OK {
        let msg = match ret {
            x if x == metis_sys::rstatus_et_METIS_ERROR_INPUT => "METIS_ERROR_INPUT: invalid input",
            x if x == metis_sys::rstatus_et_METIS_ERROR_MEMORY => {
                "METIS_ERROR_MEMORY: allocation failure"
            }
            x if x == metis_sys::rstatus_et_METIS_ERROR => "METIS_ERROR: general error",
            _ => "METIS: unknown error code",
        };
        return Err(SparseError::AnalysisFailure {
            reason: format!("{} (return code: {})", msg, ret),
        });
    }

    // Convert i32 → usize and construct Perm.
    // METIS manual: A' = A(perm, perm) — perm maps new→old (forward permutation).
    // METIS perm[new_pos] = old_idx  → faer forward array
    // METIS iperm[old_idx] = new_pos → faer inverse array
    let fwd: Box<[usize]> = perm_i32.iter().map(|&v| v as usize).collect();
    let inv: Box<[usize]> = iperm_i32.iter().map(|&v| v as usize).collect();

    Ok(Perm::new_checked(fwd, inv, n))
}

/// Construct an identity permutation of dimension `n`.
fn identity_perm(n: usize) -> Perm<usize> {
    let id: Vec<usize> = (0..n).collect();
    Perm::new_checked(id.clone().into_boxed_slice(), id.into_boxed_slice(), n)
}

/// Extract CSR adjacency arrays (xadj, adjncy) from a symmetric sparse matrix.
///
/// Converts the CSC symbolic sparsity pattern into METIS's CSR graph format.
/// Diagonal entries are excluded (no self-loops). The structural pattern is
/// symmetrized: if entry (i, j) exists, entry (j, i) is also included.
///
/// # Arguments
///
/// * `matrix` — Symbolic CSC sparsity pattern of a symmetric sparse matrix.
///
/// # Returns
///
/// A tuple `(xadj, adjncy)` where:
/// - `xadj` has length `n + 1` (row pointers)
/// - `adjncy` has length `xadj[n]` (column indices of neighbors)
///
/// Both arrays use `i32` indices as required by METIS (`idx_t`).
///
/// # Errors
///
/// Returns [`SparseError::InvalidInput`] if the total number of adjacency entries
/// exceeds `i32::MAX`.
fn extract_adjacency(
    matrix: SymbolicSparseColMatRef<'_, usize>,
) -> Result<(Vec<i32>, Vec<i32>), SparseError> {
    let n = matrix.nrows();
    let col_ptrs = matrix.col_ptr();
    let row_indices = matrix.row_idx();

    // Step 1: Build a symmetric neighbor set for each vertex, excluding diagonal.
    // Use Vec<Vec<usize>> rather than HashSet; sorted and deduped in step 2.
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];

    // The CSC structure gives us column j's row indices.
    // For symmetric adjacency, entry (i, j) with i != j means both i→j and j→i.
    for j in 0..n {
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for &i in &row_indices[start..end] {
            if i != j {
                neighbors[j].push(i);
                neighbors[i].push(j);
            }
        }
    }

    // Step 2: Sort and deduplicate each neighbor list.
    for nbrs in &mut neighbors {
        nbrs.sort_unstable();
        nbrs.dedup();
    }

    // Step 3: Validate total edge count fits in i32.
    let total_edges: usize = neighbors.iter().map(|v| v.len()).sum();
    if total_edges > i32::MAX as usize {
        return Err(SparseError::InvalidInput {
            reason: format!(
                "Total adjacency entries {} exceeds i32::MAX ({}); METIS requires 32-bit indices",
                total_edges,
                i32::MAX
            ),
        });
    }

    // Step 4: Build CSR arrays (xadj, adjncy) with i32 indices.
    let mut xadj = Vec::with_capacity(n + 1);
    xadj.push(0i32);

    let mut adjncy = Vec::with_capacity(total_edges);

    for nbrs in &neighbors {
        for &neighbor in nbrs {
            adjncy.push(neighbor as i32);
        }
        xadj.push(adjncy.len() as i32);
    }

    Ok((xadj, adjncy))
}

#[cfg(test)]
mod tests {
    use super::*;

    use faer::sparse::{SparseColMat, Triplet};

    /// Helper: create an n×n symmetric tridiagonal sparse matrix.
    fn make_tridiagonal(n: usize) -> SparseColMat<usize, f64> {
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push(Triplet::new(i, i, 4.0));
        }
        for i in 0..n - 1 {
            triplets.push(Triplet::new(i, i + 1, 1.0));
            triplets.push(Triplet::new(i + 1, i, 1.0));
        }
        SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
    }

    /// Helper: create an n×n diagonal-only sparse matrix.
    fn make_diagonal(n: usize) -> SparseColMat<usize, f64> {
        let triplets: Vec<_> = (0..n).map(|i| Triplet::new(i, i, (i + 1) as f64)).collect();
        SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
    }

    /// Helper: create an n×n arrow matrix (dense first row/column + diagonal).
    fn make_arrow(n: usize) -> SparseColMat<usize, f64> {
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push(Triplet::new(i, i, 10.0));
        }
        for i in 1..n {
            triplets.push(Triplet::new(0, i, 1.0));
            triplets.push(Triplet::new(i, 0, 1.0));
        }
        SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
    }

    /// Helper: create a matrix stored as upper triangle only.
    fn make_upper_triangle_only(n: usize) -> SparseColMat<usize, f64> {
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push(Triplet::new(i, i, 4.0));
        }
        for i in 0..n - 1 {
            // Only upper triangle: row < col
            triplets.push(Triplet::new(i, i + 1, 1.0));
        }
        SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
    }

    // ---- T008: metis_ordering unit tests ----

    #[test]
    fn test_metis_ordering_tridiagonal_5() {
        let matrix = make_tridiagonal(5);
        let perm = metis_ordering(matrix.symbolic()).expect("METIS ordering should succeed");

        // Verify valid permutation: every index 0..4 appears exactly once
        let (fwd, inv) = perm.as_ref().arrays();
        assert_eq!(fwd.len(), 5);
        assert_eq!(inv.len(), 5);

        let mut seen_fwd = [false; 5];
        let mut seen_inv = [false; 5];
        for i in 0..5 {
            assert!(fwd[i] < 5, "fwd[{}] = {} out of range", i, fwd[i]);
            assert!(inv[i] < 5, "inv[{}] = {} out of range", i, inv[i]);
            seen_fwd[fwd[i]] = true;
            seen_inv[inv[i]] = true;
        }
        assert!(
            seen_fwd.iter().all(|&s| s),
            "fwd is not a valid permutation"
        );
        assert!(
            seen_inv.iter().all(|&s| s),
            "inv is not a valid permutation"
        );
    }

    #[test]
    fn test_metis_ordering_fwd_inv_consistency() {
        let matrix = make_tridiagonal(10);
        let perm = metis_ordering(matrix.symbolic()).expect("METIS ordering should succeed");

        let (fwd, inv) = perm.as_ref().arrays();
        for i in 0..10 {
            assert_eq!(fwd[inv[i]], i, "fwd[inv[{}]] = {} != {}", i, fwd[inv[i]], i);
            assert_eq!(inv[fwd[i]], i, "inv[fwd[{}]] = {} != {}", i, inv[fwd[i]], i);
        }
    }

    #[test]
    fn test_metis_ordering_dim_0() {
        let triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
        let matrix = SparseColMat::try_new_from_triplets(0, 0, &triplets).unwrap();
        let perm = metis_ordering(matrix.symbolic()).expect("dim-0 should succeed");
        let (fwd, inv) = perm.as_ref().arrays();
        assert_eq!(fwd.len(), 0);
        assert_eq!(inv.len(), 0);
    }

    #[test]
    fn test_metis_ordering_dim_1() {
        let matrix = make_diagonal(1);
        let perm = metis_ordering(matrix.symbolic()).expect("dim-1 should succeed");
        let (fwd, inv) = perm.as_ref().arrays();
        assert_eq!(fwd, &[0]);
        assert_eq!(inv, &[0]);
    }

    #[test]
    fn test_metis_ordering_diagonal_identity() {
        // Diagonal matrix has no off-diagonal entries — should return identity perm
        let matrix = make_diagonal(5);
        let perm = metis_ordering(matrix.symbolic()).expect("diagonal should succeed");
        let (fwd, inv) = perm.as_ref().arrays();
        for i in 0..5 {
            assert_eq!(
                fwd[i], i,
                "diagonal perm should be identity: fwd[{}] = {}",
                i, fwd[i]
            );
            assert_eq!(
                inv[i], i,
                "diagonal perm should be identity: inv[{}] = {}",
                i, inv[i]
            );
        }
    }

    // ---- T005: extract_adjacency tests ----

    #[test]
    fn test_adjacency_tridiagonal_3x3() {
        // 3x3 tridiagonal: edges 0-1, 1-2
        let matrix = make_tridiagonal(3);
        let (xadj, adjncy) = extract_adjacency(matrix.symbolic()).unwrap();

        // Expected CSR adjacency (no diagonal):
        // vertex 0: neighbors [1]
        // vertex 1: neighbors [0, 2]
        // vertex 2: neighbors [1]
        assert_eq!(xadj, vec![0, 1, 3, 4]);
        assert_eq!(adjncy.len(), 4);

        // Check neighbors for each vertex
        let n0: Vec<i32> = adjncy[xadj[0] as usize..xadj[1] as usize].to_vec();
        let n1: Vec<i32> = adjncy[xadj[1] as usize..xadj[2] as usize].to_vec();
        let n2: Vec<i32> = adjncy[xadj[2] as usize..xadj[3] as usize].to_vec();
        assert_eq!(n0, vec![1]);
        assert_eq!(n1, vec![0, 2]);
        assert_eq!(n2, vec![1]);
    }

    #[test]
    fn test_adjacency_arrow_5x5() {
        // 5x5 arrow: vertex 0 connected to all others (star graph)
        let matrix = make_arrow(5);
        let (xadj, adjncy) = extract_adjacency(matrix.symbolic()).unwrap();

        // vertex 0: neighbors [1, 2, 3, 4]
        // vertex 1: neighbors [0]
        // vertex 2: neighbors [0]
        // vertex 3: neighbors [0]
        // vertex 4: neighbors [0]
        assert_eq!(xadj, vec![0, 4, 5, 6, 7, 8]);
        assert_eq!(adjncy.len(), 8);

        let n0: Vec<i32> = adjncy[xadj[0] as usize..xadj[1] as usize].to_vec();
        assert_eq!(n0, vec![1, 2, 3, 4]);
        for v in 1..5i32 {
            let start = xadj[v as usize] as usize;
            let end = xadj[v as usize + 1] as usize;
            assert_eq!(
                &adjncy[start..end],
                &[0],
                "vertex {} should only neighbor 0",
                v
            );
        }
    }

    #[test]
    fn test_adjacency_diagonal_only() {
        // Diagonal matrix: no off-diagonal entries, empty adjacency
        let matrix = make_diagonal(4);
        let (xadj, adjncy) = extract_adjacency(matrix.symbolic()).unwrap();

        assert_eq!(xadj, vec![0, 0, 0, 0, 0]);
        assert!(adjncy.is_empty());
    }

    #[test]
    fn test_adjacency_upper_triangle_only() {
        // Upper-triangle-only input should produce symmetric adjacency output
        let matrix = make_upper_triangle_only(4);
        let (xadj, adjncy) = extract_adjacency(matrix.symbolic()).unwrap();

        // Tridiag structure (edges 0-1, 1-2, 2-3) stored as upper-only
        // Result should still be symmetric:
        // vertex 0: [1]
        // vertex 1: [0, 2]
        // vertex 2: [1, 3]
        // vertex 3: [2]
        assert_eq!(xadj, vec![0, 1, 3, 5, 6]);
        assert_eq!(adjncy.len(), 6);

        // Verify symmetry: for each (i, j) in adjncy, also (j, i)
        let n = 4;
        for i in 0..n {
            let start = xadj[i] as usize;
            let end = xadj[i + 1] as usize;
            for &j in &adjncy[start..end] {
                let j_start = xadj[j as usize] as usize;
                let j_end = xadj[j as usize + 1] as usize;
                assert!(
                    adjncy[j_start..j_end].contains(&(i as i32)),
                    "symmetry violated: ({}, {}) exists but ({}, {}) does not",
                    i,
                    j,
                    j,
                    i
                );
            }
        }
    }

    #[test]
    fn test_adjacency_invariants() {
        // Test structural invariants on various matrices
        for (name, matrix) in [
            ("tridiag3", make_tridiagonal(3)),
            ("arrow5", make_arrow(5)),
            ("diag4", make_diagonal(4)),
            ("upper4", make_upper_triangle_only(4)),
        ] {
            let n = matrix.nrows();
            let (xadj, adjncy) = extract_adjacency(matrix.symbolic()).unwrap();

            // xadj[0] == 0
            assert_eq!(xadj[0], 0, "{}: xadj[0] should be 0", name);

            // xadj has length n+1
            assert_eq!(xadj.len(), n + 1, "{}: xadj length should be n+1", name);

            // xadj is monotonically non-decreasing
            for i in 0..n {
                assert!(
                    xadj[i + 1] >= xadj[i],
                    "{}: xadj not monotonic at {}: {} > {}",
                    name,
                    i,
                    xadj[i],
                    xadj[i + 1]
                );
            }

            // No self-loops
            for i in 0..n {
                let start = xadj[i] as usize;
                let end = xadj[i + 1] as usize;
                for &j in &adjncy[start..end] {
                    assert_ne!(j, i as i32, "{}: self-loop at vertex {}", name, i);
                }
            }

            // All indices in [0, n)
            for &j in &adjncy {
                assert!(
                    j >= 0 && (j as usize) < n,
                    "{}: adjncy index {} out of range [0, {})",
                    name,
                    j,
                    n
                );
            }

            // Symmetry
            for i in 0..n {
                let start = xadj[i] as usize;
                let end = xadj[i + 1] as usize;
                for &j in &adjncy[start..end] {
                    let j_start = xadj[j as usize] as usize;
                    let j_end = xadj[j as usize + 1] as usize;
                    assert!(
                        adjncy[j_start..j_end].contains(&(i as i32)),
                        "{}: symmetry violated: ({}, {}) without ({}, {})",
                        name,
                        i,
                        j,
                        j,
                        i
                    );
                }
            }
        }
    }
}
