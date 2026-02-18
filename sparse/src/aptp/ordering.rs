//! Ordering algorithms for symmetric sparse matrices.
//!
//! This module provides fill-reducing permutations for symmetric sparse matrices:
//!
//! - [`metis_ordering()`] — METIS nested dissection on the raw sparsity pattern
//! - [`match_order_metis()`] — Combined MC64 matching + METIS ordering with cycle
//!   condensation, guaranteeing matched 2-cycle pairs are adjacent in the
//!   elimination order (SPRAL `ordering=2` mode)
//!
//! The resulting permutations integrate with
//! [`AptpSymbolic::analyze()`](super::AptpSymbolic::analyze) via
//! [`SymmetricOrdering::Custom`](faer::sparse::linalg::cholesky::SymmetricOrdering::Custom).
//!
//! # Algorithm References
//!
//! - Karypis & Kumar (1998), "A Fast and High Quality Multilevel Scheme for
//!   Partitioning Irregular Graphs", SIAM J. Sci. Comput. 20(1)
//! - George (1973), "Nested Dissection of a Regular Finite Element Mesh",
//!   SIAM J. Numer. Anal. 10(2)
//! - Hogg & Scott (2013), HSL_MC80 — cycle condensation approach
//! - Duff & Koster (2001), Algorithm MPD — MC64 matching
//! - Duff & Pralet (2005) — MC64SYM symmetric scaling
//! - SPRAL `match_order.f90` (BSD-3) — reference implementation
//!
//! # Implementation Notes
//!
//! Uses `metis-sys` for raw FFI bindings to METIS 5.x (vendored C source,
//! compiled via `cc`). The `METIS_NodeND` function computes nested dissection
//! orderings on undirected graphs extracted from the matrix sparsity pattern.

use faer::perm::Perm;
use faer::sparse::{SparseColMat, SymbolicSparseColMatRef};

use super::matching::{Mc64Job, mc64_matching};
use crate::error::SparseError;

/// Result of the combined matching-ordering pipeline.
///
/// Contains both the fill-reducing ordering (for symbolic analysis) and
/// the MC64 scaling factors (for numeric factorization), plus diagnostics
/// about the condensation process.
///
/// # Usage
///
/// ```no_run
/// use rivrs_sparse::aptp::match_order_metis;
/// use faer::sparse::linalg::cholesky::SymmetricOrdering;
/// # let matrix = todo!();
///
/// let result = match_order_metis(&matrix).unwrap();
/// // Use ordering for symbolic analysis:
/// // SymmetricOrdering::Custom(result.ordering.as_ref())
/// ```
pub struct MatchOrderResult {
    /// Fill-reducing permutation with matched pair adjacency guarantee.
    /// Use with `SymmetricOrdering::Custom(result.ordering.as_ref())`.
    pub ordering: Perm<usize>,

    /// MC64 symmetric scaling factors (linear domain).
    /// Apply as `A_scaled[i,j] = scaling[i] * A[i,j] * scaling[j]`.
    pub scaling: Vec<f64>,

    /// Number of matched entries from MC64.
    /// Equals n for structurally nonsingular matrices.
    pub matched: usize,

    /// Dimension of the condensed graph passed to METIS.
    /// Strictly less than n when 2-cycles exist.
    pub condensed_dim: usize,

    /// Number of singleton nodes (self-matched).
    pub singletons: usize,

    /// Number of 2-cycle pairs in the decomposed matching.
    pub two_cycles: usize,
}

/// Compute a combined MC64 matching + METIS ordering with cycle condensation.
///
/// This implements SPRAL's `match_order_metis` pipeline:
/// 1. MC64 matching → scaling + matching permutation
/// 2. Cycle splitting → decompose matching into singletons and 2-cycles
/// 3. Condensed graph → fuse 2-cycle pairs into single super-nodes
/// 4. METIS ordering → fill-reducing ordering on condensed graph
/// 5. Expansion → map back to original indices with pair adjacency
///
/// Matched 2-cycle pairs are guaranteed to occupy consecutive positions
/// in the output ordering, making them natural candidates for 2x2 pivots
/// in APTP factorization.
///
/// # Arguments
///
/// * `matrix` — Sparse symmetric matrix in upper-triangular CSC format.
///   Must be square. Numeric values required (used by MC64).
///
/// # Returns
///
/// * `Ok(MatchOrderResult)` — Ordering, scaling, and diagnostics.
/// * `Err(SparseError::NotSquare)` — Matrix is not square.
/// * `Err(SparseError::InvalidInput)` — Zero dimension or non-finite entries.
/// * `Err(SparseError::AnalysisFailure)` — METIS or MC64 internal failure.
///
/// # Algorithm References
///
/// - Hogg & Scott (2013), HSL_MC80 — cycle condensation approach
/// - Duff & Koster (2001), Algorithm MPD — MC64 matching
/// - Duff & Pralet (2005) — MC64SYM symmetric scaling
/// - SPRAL `match_order.f90` (BSD-3) — reference implementation
///
/// # SPRAL Equivalent
///
/// Corresponds to `spral_match_order::match_order_metis` with
/// `options%ordering = 2` in SSIDS.
pub fn match_order_metis(
    matrix: &SparseColMat<usize, f64>,
) -> Result<MatchOrderResult, SparseError> {
    let n = matrix.nrows();

    // Trivial case: dim 0 or 1 — return identity ordering with empty matching info
    if n <= 1 {
        let scaling = if n == 1 { vec![1.0] } else { vec![] };
        return Ok(MatchOrderResult {
            ordering: identity_perm(n),
            scaling,
            matched: n,
            condensed_dim: n,
            singletons: n,
            two_cycles: 0,
        });
    }

    // Step 1: MC64 matching
    let mc64_result = mc64_matching(matrix, Mc64Job::MaximumProduct)?;

    // Step 2: Cycle decomposition
    let matching_fwd = mc64_result.matching.as_ref().arrays().0;
    let decomp = split_matching_cycles(matching_fwd, &mc64_result.is_matched, n);

    // Step 3: Build condensed adjacency graph
    let (mut xadj, mut adjncy) = build_condensed_adjacency(matrix.symbolic(), &decomp)?;

    // Step 4: METIS on condensed graph
    let cdim = decomp.condensed_dim;
    let condensed_order = if cdim <= 1 || adjncy.is_empty() {
        // Trivial condensed graph — identity ordering
        (0..cdim as i32).collect::<Vec<_>>()
    } else {
        if cdim > i32::MAX as usize {
            return Err(SparseError::InvalidInput {
                reason: format!(
                    "Condensed dimension {} exceeds i32::MAX; METIS requires 32-bit indices",
                    cdim
                ),
            });
        }

        // METIS iperm: iperm[old_idx] = new_pos (this is the "order" we need)
        let (_perm, iperm) = call_metis_node_nd(cdim, &mut xadj, &mut adjncy)?;
        iperm
    };

    // Step 5: Expand ordering
    let ordering = expand_ordering(&condensed_order, &decomp, n);

    Ok(MatchOrderResult {
        ordering,
        scaling: mc64_result.scaling,
        matched: mc64_result.matched,
        condensed_dim: decomp.condensed_dim,
        singletons: decomp.singletons,
        two_cycles: decomp.two_cycles,
    })
}

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

    // Call METIS_NodeND via shared helper.
    // METIS perm[new_pos] = old_idx  → faer forward array
    // METIS iperm[old_idx] = new_pos → faer inverse array
    let (perm_i32, iperm_i32) = call_metis_node_nd(n, &mut xadj, &mut adjncy)?;

    let fwd: Box<[usize]> = perm_i32.iter().map(|&v| v as usize).collect();
    let inv: Box<[usize]> = iperm_i32.iter().map(|&v| v as usize).collect();

    Ok(Perm::new_checked(fwd, inv, n))
}

/// Result of decomposing an MC64 matching permutation into singletons and 2-cycles.
///
/// Longer cycles in the matching are split into 2-cycles plus at most one singleton
/// (for odd-length cycles), following SPRAL's `mo_split` algorithm.
///
/// Each matched index is mapped to a condensed super-node index, where 2-cycle pairs
/// share the same super-node. Unmatched indices are excluded from the condensed graph.
struct CycleDecomposition {
    /// Per-node classification:
    /// - `-2`: unmatched by MC64
    /// - `-1`: singleton (matched to self)
    /// - `>= 0`: partner index (forms a 2-cycle pair)
    partner: Vec<isize>,

    /// Original index → condensed super-node index.
    /// Both members of a 2-cycle pair map to the same condensed index.
    /// Unmatched indices have undefined values (not used).
    old_to_new: Vec<usize>,

    /// Condensed super-node index → first original index.
    /// For singletons, this is the singleton index.
    /// For 2-cycle pairs, this is one member (the other is found via `partner`).
    new_to_old: Vec<usize>,

    /// Number of condensed super-nodes (singletons + 2-cycle pairs, matched only).
    condensed_dim: usize,

    /// Count of singleton nodes (matched to self).
    singletons: usize,

    /// Count of 2-cycle pairs.
    two_cycles: usize,
}

/// Decompose an MC64 matching into singletons, 2-cycles, and unmatched indices.
///
/// Follows SPRAL's `mo_split` algorithm: walk each cycle in the matching permutation,
/// pairing consecutive members into 2-cycles. Odd-length cycles produce one extra
/// singleton. The `is_matched` slice distinguishes true singletons (matched, `fwd[i]==i`)
/// from unmatched indices (assigned to arbitrary free columns by `build_singular_permutation`).
///
/// # Arguments
///
/// * `matching_fwd` — Forward permutation from MC64 matching (`Mc64Result.matching`)
/// * `is_matched` — Per-index matched status from `Mc64Result.is_matched`
/// * `n` — Matrix dimension
fn split_matching_cycles(
    matching_fwd: &[usize],
    is_matched: &[bool],
    n: usize,
) -> CycleDecomposition {
    let mut partner = vec![isize::MIN; n];
    let mut old_to_new = vec![0usize; n];
    let mut new_to_old = Vec::new();
    let mut visited = vec![false; n];
    let mut singletons = 0usize;
    let mut two_cycles = 0usize;

    // First pass: mark unmatched indices (visited, so cycle walker skips them)
    for i in 0..n {
        if !is_matched[i] {
            partner[i] = -2;
            visited[i] = true;
        }
    }

    // Second pass: walk cycles, split into singletons + 2-cycles.
    // Following SPRAL's mo_split: ALL nodes get condensed indices (including
    // unmatched, which become singleton condensed nodes for METIS).
    let mut condensed_idx = 0usize;
    for i in 0..n {
        if visited[i] && partner[i] != -2 {
            // Already processed as part of a matched cycle
            continue;
        }

        // Unmatched node: assign singleton condensed index
        if partner[i] == -2 {
            old_to_new[i] = condensed_idx;
            new_to_old.push(i);
            condensed_idx += 1;
            singletons += 1;
            continue;
        }

        // Check for singleton (self-matched)
        if matching_fwd[i] == i {
            partner[i] = -1;
            old_to_new[i] = condensed_idx;
            new_to_old.push(i);
            condensed_idx += 1;
            singletons += 1;
            visited[i] = true;
            continue;
        }

        // Walk the cycle starting at i, pairing consecutive members
        let mut j = i;
        loop {
            let k = matching_fwd[j];
            if visited[k] || k == i {
                // Odd cycle: j is the leftover — becomes singleton
                if !visited[j] {
                    partner[j] = -1;
                    old_to_new[j] = condensed_idx;
                    new_to_old.push(j);
                    condensed_idx += 1;
                    singletons += 1;
                    visited[j] = true;
                }
                break;
            }

            // Pair j and k as a 2-cycle
            partner[j] = k as isize;
            partner[k] = j as isize;
            old_to_new[j] = condensed_idx;
            old_to_new[k] = condensed_idx;
            new_to_old.push(j);
            condensed_idx += 1;
            two_cycles += 1;
            visited[j] = true;
            visited[k] = true;

            // Advance past k
            let next = matching_fwd[k];
            if next == i {
                break;
            }
            j = next;
            // If j was already visited (e.g. an unmatched node from the first pass),
            // the permutation cycle passes through a "hole" — stop walking here.
            // The remaining matched indices beyond the hole will be picked up when
            // the outer loop reaches them.
            if visited[j] {
                break;
            }
        }
    }

    CycleDecomposition {
        partner,
        old_to_new,
        new_to_old,
        condensed_dim: condensed_idx,
        singletons,
        two_cycles,
    }
}

/// Build condensed CSR adjacency arrays (xadj, adjncy) for METIS from the
/// original matrix sparsity pattern and cycle decomposition.
///
/// Each condensed super-node absorbs all edges from its constituent original
/// indices. Duplicate edges and self-loops are removed using a marker array.
/// The output is a full symmetric graph in METIS CSR format.
///
/// # Arguments
///
/// * `matrix` — Symbolic sparsity pattern of the original matrix
/// * `decomp` — Cycle decomposition from `split_matching_cycles`
fn build_condensed_adjacency(
    matrix: SymbolicSparseColMatRef<'_, usize>,
    decomp: &CycleDecomposition,
) -> Result<(Vec<i32>, Vec<i32>), SparseError> {
    let n = matrix.nrows();
    let col_ptrs = matrix.col_ptr();
    let row_indices = matrix.row_idx();
    let cdim = decomp.condensed_dim;

    // Build symmetric neighbor lists for condensed nodes using marker dedup
    let mut neighbors: Vec<Vec<i32>> = vec![Vec::new(); cdim];
    let mut marker = vec![usize::MAX; cdim]; // marker[cnode] = current source cnode

    for col in 0..n {
        let c_col = decomp.old_to_new[col];

        // For 2-cycle second members: old_to_new maps to same condensed index
        // as the first member. The marker dedup handles this — we process both
        // columns of a pair, and all edges merge into the same condensed node.

        // Mark self to avoid self-loops
        marker[c_col] = c_col;

        let start = col_ptrs[col];
        let end = col_ptrs[col + 1];
        for &row in &row_indices[start..end] {
            let c_row = decomp.old_to_new[row];
            if marker[c_row] != c_col {
                marker[c_row] = c_col;
                // Add edge in both directions (symmetric)
                neighbors[c_col].push(c_row as i32);
                neighbors[c_row].push(c_col as i32);
            }
        }
    }

    // Deduplicate (the bidirectional insertion can create duplicates when
    // both col→row and row→col are visited from the original CSC)
    for nbrs in &mut neighbors {
        nbrs.sort_unstable();
        nbrs.dedup();
    }

    // Validate total edges fit in i32
    let total_edges: usize = neighbors.iter().map(|v| v.len()).sum();
    if total_edges > i32::MAX as usize {
        return Err(SparseError::InvalidInput {
            reason: format!(
                "Condensed adjacency entries {} exceeds i32::MAX",
                total_edges
            ),
        });
    }

    // Build CSR arrays
    let mut xadj = Vec::with_capacity(cdim + 1);
    xadj.push(0i32);
    let mut adjncy = Vec::with_capacity(total_edges);

    for nbrs in &neighbors {
        adjncy.extend_from_slice(nbrs);
        xadj.push(adjncy.len() as i32);
    }

    Ok((xadj, adjncy))
}

/// Expand a condensed METIS ordering back to full-size permutation.
///
/// Maps each condensed super-node position to its original index(es).
/// 2-cycle pairs get consecutive positions, preserving pair adjacency.
/// Unmatched indices are appended at the end.
///
/// # Arguments
///
/// * `condensed_order` — METIS output: `condensed_order[i]` = position of condensed node `i`
/// * `decomp` — Cycle decomposition from `split_matching_cycles`
/// * `n` — Original matrix dimension
fn expand_ordering(condensed_order: &[i32], decomp: &CycleDecomposition, n: usize) -> Perm<usize> {
    let cdim = decomp.condensed_dim;

    // Build inverse of METIS output: inv_order[position] = condensed_node
    let mut inv_order = vec![0usize; cdim];
    for (cnode, &pos) in condensed_order.iter().enumerate() {
        inv_order[pos as usize] = cnode;
    }

    // Walk positions in order, emitting original indices
    let mut fwd = Vec::with_capacity(n); // fwd[new_pos] = old_idx
    for &cnode in &inv_order {
        let orig = decomp.new_to_old[cnode];
        fwd.push(orig);

        // If this condensed node is a 2-cycle pair, emit partner next
        if decomp.partner[orig] >= 0 {
            fwd.push(decomp.partner[orig] as usize);
        }
    }

    // All nodes (matched + unmatched) have condensed indices and are emitted
    // by the loop above. Unmatched nodes are singleton condensed nodes placed
    // at METIS-chosen positions, not appended at the end.
    debug_assert_eq!(fwd.len(), n);

    // Build inverse: inv[old_idx] = new_pos
    let mut inv = vec![0usize; n];
    for (new_pos, &old_idx) in fwd.iter().enumerate() {
        inv[old_idx] = new_pos;
    }

    Perm::new_checked(fwd.into_boxed_slice(), inv.into_boxed_slice(), n)
}

/// Construct an identity permutation of dimension `n`.
fn identity_perm(n: usize) -> Perm<usize> {
    let id: Vec<usize> = (0..n).collect();
    Perm::new_checked(id.clone().into_boxed_slice(), id.into_boxed_slice(), n)
}

/// Call METIS_NodeND on a CSR adjacency graph and return (perm, iperm) as i32 vectors.
///
/// Handles options setup, FFI call, and error mapping. The caller is responsible
/// for validating that `n > 1`, `adjncy` is non-empty, and `n <= i32::MAX`.
///
/// Returns `(perm, iperm)` where:
/// - `perm[new_pos] = old_idx` (forward permutation)
/// - `iperm[old_idx] = new_pos` (inverse permutation)
fn call_metis_node_nd(
    n: usize,
    xadj: &mut [i32],
    adjncy: &mut [i32],
) -> Result<(Vec<i32>, Vec<i32>), SparseError> {
    let mut options = [0i32; metis_sys::METIS_NOPTIONS as usize];
    unsafe {
        metis_sys::METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    options[metis_sys::moptions_et_METIS_OPTION_NUMBERING as usize] = 0;

    let mut nvtxs = n as i32;
    let mut perm_i32 = vec![0i32; n];
    let mut iperm_i32 = vec![0i32; n];

    let ret = unsafe {
        metis_sys::METIS_NodeND(
            &mut nvtxs,
            xadj.as_mut_ptr(),
            adjncy.as_mut_ptr(),
            std::ptr::null_mut(),
            options.as_mut_ptr(),
            perm_i32.as_mut_ptr(),
            iperm_i32.as_mut_ptr(),
        )
    };

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
            reason: format!("{} (return code: {}, dim: {})", msg, ret, n),
        });
    }

    Ok((perm_i32, iperm_i32))
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

    // ---- T003: split_matching_cycles unit tests ----

    #[test]
    fn test_split_all_singletons() {
        // Identity matching: fwd[i] = i, all matched
        let fwd = vec![0, 1, 2, 3];
        let is_matched = vec![true, true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 4);

        assert_eq!(decomp.singletons, 4);
        assert_eq!(decomp.two_cycles, 0);
        assert_eq!(decomp.condensed_dim, 4);
        for i in 0..4 {
            assert_eq!(decomp.partner[i], -1, "index {} should be singleton", i);
        }
    }

    #[test]
    fn test_split_pure_2_cycles() {
        // Two 2-cycles: 0↔1, 2↔3
        let fwd = vec![1, 0, 3, 2];
        let is_matched = vec![true, true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 4);

        assert_eq!(decomp.singletons, 0);
        assert_eq!(decomp.two_cycles, 2);
        assert_eq!(decomp.condensed_dim, 2);
        assert_eq!(decomp.partner[0], 1);
        assert_eq!(decomp.partner[1], 0);
        assert_eq!(decomp.partner[2], 3);
        assert_eq!(decomp.partner[3], 2);
        // Partners share condensed index
        assert_eq!(decomp.old_to_new[0], decomp.old_to_new[1]);
        assert_eq!(decomp.old_to_new[2], decomp.old_to_new[3]);
        assert_ne!(decomp.old_to_new[0], decomp.old_to_new[2]);
    }

    #[test]
    fn test_split_3_cycle() {
        // 3-cycle: 0→1→2→0
        let fwd = vec![1, 2, 0];
        let is_matched = vec![true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 3);

        // Should split into 1 two-cycle + 1 singleton
        assert_eq!(decomp.two_cycles, 1);
        assert_eq!(decomp.singletons, 1);
        assert_eq!(decomp.condensed_dim, 2);

        // Count partner types
        let pairs: Vec<_> = (0..3).filter(|&i| decomp.partner[i] >= 0).collect();
        let sings: Vec<_> = (0..3).filter(|&i| decomp.partner[i] == -1).collect();
        assert_eq!(pairs.len(), 2);
        assert_eq!(sings.len(), 1);
    }

    #[test]
    fn test_split_4_cycle() {
        // 4-cycle: 0→1→2→3→0
        let fwd = vec![1, 2, 3, 0];
        let is_matched = vec![true, true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 4);

        // Should split into 2 two-cycles
        assert_eq!(decomp.two_cycles, 2);
        assert_eq!(decomp.singletons, 0);
        assert_eq!(decomp.condensed_dim, 2);
    }

    #[test]
    fn test_split_5_cycle() {
        // 5-cycle: 0→1→2→3→4→0
        let fwd = vec![1, 2, 3, 4, 0];
        let is_matched = vec![true, true, true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 5);

        // Should split into 2 two-cycles + 1 singleton
        assert_eq!(decomp.two_cycles, 2);
        assert_eq!(decomp.singletons, 1);
        assert_eq!(decomp.condensed_dim, 3);
    }

    #[test]
    fn test_split_mixed_with_unmatched() {
        // 6 nodes: 0↔1 (2-cycle), 2 singleton, 3↔4 (2-cycle), 5 unmatched
        let fwd = vec![1, 0, 2, 4, 3, 5]; // 5 maps to self but unmatched
        let is_matched = vec![true, true, true, true, true, false];
        let decomp = split_matching_cycles(&fwd, &is_matched, 6);

        assert_eq!(decomp.two_cycles, 2);
        // Singletons include both matched singletons and unmatched-as-singletons
        assert_eq!(decomp.singletons, 2); // {2} matched singleton + {5} unmatched singleton
        assert_eq!(decomp.condensed_dim, 4); // All nodes: 2 pairs + 2 singletons
        assert_eq!(decomp.partner[5], -2); // Still tagged as unmatched
        assert_eq!(decomp.partner[2], -1); // Matched singleton
        // Unmatched node 5 has a condensed index
        assert!(
            decomp.old_to_new[5] < decomp.condensed_dim,
            "unmatched node should have valid condensed index"
        );
    }

    #[test]
    fn test_split_trivial_n0() {
        let decomp = split_matching_cycles(&[], &[], 0);
        assert_eq!(decomp.condensed_dim, 0);
        assert_eq!(decomp.singletons, 0);
        assert_eq!(decomp.two_cycles, 0);
    }

    #[test]
    fn test_split_trivial_n1() {
        let decomp = split_matching_cycles(&[0], &[true], 1);
        assert_eq!(decomp.condensed_dim, 1);
        assert_eq!(decomp.singletons, 1);
        assert_eq!(decomp.two_cycles, 0);
        assert_eq!(decomp.partner[0], -1);
    }

    // ---- T005: build_condensed_adjacency unit tests ----

    /// Helper: create a symmetric matrix from upper-triangular entries.
    fn make_upper_tri(n: usize, entries: &[(usize, usize, f64)]) -> SparseColMat<usize, f64> {
        let triplets: Vec<_> = entries
            .iter()
            .map(|&(i, j, v)| Triplet::new(i, j, v))
            .collect();
        SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
    }

    #[test]
    fn test_condensed_adjacency_arrow_with_2cycle() {
        // 4x4 arrow: 0 connected to all, plus 2-cycle {0,1}
        // After condensation with pair {0,1}: condensed has 3 nodes
        //   cnode 0 = {0,1} (pair), cnode 1 = {2}, cnode 2 = {3}
        // Edges: {0,1} connects to {2}, {3}; {2} connects to {0,1}; {3} connects to {0,1}
        let matrix = make_arrow(4);

        // Simulate matching: 0↔1, 2=singleton, 3=singleton
        let fwd = vec![1, 0, 2, 3];
        let is_matched = vec![true, true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 4);

        assert_eq!(decomp.condensed_dim, 3);

        let (xadj, adjncy) = build_condensed_adjacency(matrix.symbolic(), &decomp).unwrap();

        // Condensed should have 3 nodes
        assert_eq!(xadj.len(), 4); // 3+1

        // No self-loops
        for i in 0..3 {
            let start = xadj[i] as usize;
            let end = xadj[i + 1] as usize;
            for &j in &adjncy[start..end] {
                assert_ne!(j, i as i32, "self-loop at condensed vertex {}", i);
            }
        }

        // Total edges should be deduplicated
        // The pair node {0,1} absorbs edges from both 0 and 1
        // No duplicates should remain
        let total = adjncy.len();
        assert!(total > 0, "condensed graph should have edges");
    }

    #[test]
    fn test_condensed_adjacency_diagonal() {
        // Diagonal matrix: no off-diagonal edges
        let matrix = make_diagonal(4);
        let fwd = vec![0, 1, 2, 3];
        let is_matched = vec![true, true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 4);

        let (xadj, adjncy) = build_condensed_adjacency(matrix.symbolic(), &decomp).unwrap();
        assert!(
            adjncy.is_empty(),
            "diagonal should have zero condensed edges"
        );
        assert_eq!(xadj, vec![0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_condensed_adjacency_full_4x4_two_pairs() {
        // 4x4 fully connected (upper tri), 2 pairs: {0,1}, {2,3}
        let matrix = make_upper_tri(
            4,
            &[
                (0, 0, 1.0),
                (0, 1, 1.0),
                (0, 2, 1.0),
                (0, 3, 1.0),
                (1, 1, 1.0),
                (1, 2, 1.0),
                (1, 3, 1.0),
                (2, 2, 1.0),
                (2, 3, 1.0),
                (3, 3, 1.0),
            ],
        );
        let fwd = vec![1, 0, 3, 2];
        let is_matched = vec![true, true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 4);

        assert_eq!(decomp.condensed_dim, 2);

        let (xadj, adjncy) = build_condensed_adjacency(matrix.symbolic(), &decomp).unwrap();

        // 2 condensed nodes, connected to each other
        assert_eq!(xadj.len(), 3);
        // Each should neighbor the other (symmetric)
        assert_eq!(adjncy.len(), 2);

        // Verify symmetry
        let n0_start = xadj[0] as usize;
        let n0_end = xadj[1] as usize;
        let n1_start = xadj[1] as usize;
        let n1_end = xadj[2] as usize;
        assert!(adjncy[n0_start..n0_end].contains(&1));
        assert!(adjncy[n1_start..n1_end].contains(&0));
    }

    // ---- T007: expand_ordering unit tests ----

    #[test]
    fn test_expand_ordering_with_pairs() {
        // 4 nodes, 2 pairs: {0,1} and {2,3}
        // condensed_dim = 2, condensed_order = [1, 0] (swap order)
        let fwd = vec![1, 0, 3, 2];
        let is_matched = vec![true, true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 4);

        // condensed_order[cnode] = position
        // If cnode 0 ({0,1}) gets position 1 and cnode 1 ({2,3}) gets position 0:
        let condensed_order = vec![1i32, 0i32];
        let perm = expand_ordering(&condensed_order, &decomp, 4);

        let (expanded_fwd, expanded_inv) = perm.as_ref().arrays();
        assert_eq!(expanded_fwd.len(), 4);

        // Verify valid permutation
        let mut seen = [false; 4];
        for &v in expanded_fwd {
            assert!(v < 4);
            seen[v] = true;
        }
        assert!(seen.iter().all(|&s| s), "not a valid permutation");

        // Verify pair adjacency: partners must be consecutive
        let pos_0 = expanded_inv[0];
        let pos_1 = expanded_inv[1];
        assert_eq!(pos_0.abs_diff(pos_1), 1, "pair (0,1) not consecutive");

        let pos_2 = expanded_inv[2];
        let pos_3 = expanded_inv[3];
        assert_eq!(pos_2.abs_diff(pos_3), 1, "pair (2,3) not consecutive");
    }

    #[test]
    fn test_expand_ordering_mixed_singletons_pairs() {
        // 5 nodes: pair {0,1}, singleton {2}, pair {3,4}
        let fwd = vec![1, 0, 2, 4, 3];
        let is_matched = vec![true, true, true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 5);

        assert_eq!(decomp.condensed_dim, 3);

        // Identity condensed ordering
        let condensed_order = vec![0i32, 1i32, 2i32];
        let perm = expand_ordering(&condensed_order, &decomp, 5);

        let (expanded_fwd, expanded_inv) = perm.as_ref().arrays();
        assert_eq!(expanded_fwd.len(), 5);

        // Valid permutation
        let mut seen = [false; 5];
        for &v in expanded_fwd {
            seen[v] = true;
        }
        assert!(seen.iter().all(|&s| s));

        // Pair {0,1} consecutive
        let diff = (expanded_inv[0] as isize - expanded_inv[1] as isize).unsigned_abs();
        assert_eq!(diff, 1, "pair (0,1) not consecutive");

        // Pair {3,4} consecutive
        let diff = (expanded_inv[3] as isize - expanded_inv[4] as isize).unsigned_abs();
        assert_eq!(diff, 1, "pair (3,4) not consecutive");
    }

    #[test]
    fn test_expand_ordering_with_unmatched() {
        // 5 nodes: pair {0,1}, singleton {2}, unmatched {3}, unmatched {4}
        let fwd = vec![1, 0, 2, 3, 4];
        let is_matched = vec![true, true, true, false, false];
        let decomp = split_matching_cycles(&fwd, &is_matched, 5);

        // Now includes unmatched as singleton condensed nodes
        assert_eq!(decomp.condensed_dim, 4); // pair + singleton + 2 unmatched

        let condensed_order = vec![0i32, 1i32, 2i32, 3i32];
        let perm = expand_ordering(&condensed_order, &decomp, 5);

        let (expanded_fwd, expanded_inv) = perm.as_ref().arrays();
        assert_eq!(expanded_fwd.len(), 5);

        // All nodes should have valid positions
        let mut seen = [false; 5];
        for &v in expanded_fwd {
            assert!(v < 5);
            seen[v] = true;
        }
        assert!(seen.iter().all(|&s| s), "not a valid permutation");

        // Pair {0,1} still consecutive
        let diff = (expanded_inv[0] as isize - expanded_inv[1] as isize).unsigned_abs();
        assert_eq!(diff, 1, "pair (0,1) not consecutive");

        // Unmatched nodes are now at METIS-determined positions (not forced to end)
        // With identity condensed ordering, they get positions based on their
        // condensed index, not appended
    }

    // ---- Unmatched node condensation fix tests ----

    #[test]
    fn test_split_unmatched_get_condensed_indices() {
        // 4 nodes: pair {0,1}, unmatched {2}, unmatched {3}
        let fwd = vec![1, 0, 2, 3]; // 2,3 map to self but unmatched
        let is_matched = vec![true, true, false, false];
        let decomp = split_matching_cycles(&fwd, &is_matched, 4);

        // Both unmatched nodes should have condensed indices
        assert_eq!(decomp.condensed_dim, 3); // 1 pair + 2 unmatched singletons
        assert!(decomp.old_to_new[2] < decomp.condensed_dim);
        assert!(decomp.old_to_new[3] < decomp.condensed_dim);

        // Unmatched nodes have distinct condensed indices
        assert_ne!(decomp.old_to_new[2], decomp.old_to_new[3]);

        // Pair members share condensed index
        assert_eq!(decomp.old_to_new[0], decomp.old_to_new[1]);
    }

    #[test]
    fn test_condensed_adjacency_includes_unmatched_edges() {
        // 4x4: pair {0,1} connected to unmatched {2}. {3} isolated unmatched.
        // [1  1  1  0]
        // [1  1  0  0]
        // [1  0  1  0]
        // [0  0  0  1]
        let matrix = make_upper_tri(
            4,
            &[
                (0, 0, 1.0),
                (0, 1, 1.0),
                (0, 2, 1.0),
                (1, 1, 1.0),
                (2, 2, 1.0),
                (3, 3, 1.0),
            ],
        );

        let fwd = vec![1, 0, 2, 3];
        let is_matched = vec![true, true, false, false];
        let decomp = split_matching_cycles(&fwd, &is_matched, 4);

        let (xadj, adjncy) = build_condensed_adjacency(matrix.symbolic(), &decomp).unwrap();

        // Should have 3 condensed nodes: pair{0,1}, singleton{2}, singleton{3}
        assert_eq!(xadj.len(), 4); // 3 + 1

        // The pair {0,1} should have edge to unmatched singleton {2}
        // (since original node 0 connects to node 2)
        let c_pair = decomp.old_to_new[0];
        let c_unmatched_2 = decomp.old_to_new[2];

        let pair_start = xadj[c_pair] as usize;
        let pair_end = xadj[c_pair + 1] as usize;
        let pair_neighbors: Vec<i32> = adjncy[pair_start..pair_end].to_vec();
        assert!(
            pair_neighbors.contains(&(c_unmatched_2 as i32)),
            "pair should neighbor unmatched node 2; neighbors = {:?}",
            pair_neighbors
        );

        // The unmatched singleton {3} should be isolated (no edges except diagonal)
        let c_unmatched_3 = decomp.old_to_new[3];
        let n3_start = xadj[c_unmatched_3] as usize;
        let n3_end = xadj[c_unmatched_3 + 1] as usize;
        assert_eq!(
            n3_end - n3_start,
            0,
            "isolated unmatched node 3 should have no neighbors"
        );
    }

    #[test]
    fn test_expand_unmatched_not_at_end() {
        // 6 nodes: pair {0,1}, unmatched {2}, singleton {3}, pair {4,5}
        // With condensed_order that places unmatched first, it should appear first
        let fwd = vec![1, 0, 2, 3, 5, 4];
        let is_matched = vec![true, true, false, true, true, true];
        let decomp = split_matching_cycles(&fwd, &is_matched, 6);

        // condensed_dim should include unmatched
        let cdim = decomp.condensed_dim;

        // Create a condensed ordering that puts the unmatched node first
        let c_unmatched = decomp.old_to_new[2];
        let mut condensed_order = vec![0i32; cdim];
        // Put unmatched at position 0
        condensed_order[c_unmatched] = 0;
        // Fill remaining positions
        let mut pos = 1i32;
        for cnode in 0..cdim {
            if cnode == c_unmatched {
                continue;
            }
            condensed_order[cnode] = pos;
            pos += 1;
        }

        let perm = expand_ordering(&condensed_order, &decomp, 6);
        let (expanded_fwd, expanded_inv) = perm.as_ref().arrays();

        // Unmatched node 2 should be at position 0 (not at the end)
        assert_eq!(
            expanded_inv[2], 0,
            "unmatched node should be at position 0, got {}",
            expanded_inv[2]
        );

        // Valid permutation
        let mut seen = [false; 6];
        for &v in expanded_fwd {
            assert!(v < 6);
            seen[v] = true;
        }
        assert!(seen.iter().all(|&s| s));
    }

    #[test]
    fn test_all_unmatched_gets_ordering() {
        // Edge case: all nodes unmatched (empty matching)
        let fwd = vec![0, 1, 2, 3];
        let is_matched = vec![false, false, false, false];
        let decomp = split_matching_cycles(&fwd, &is_matched, 4);

        // All should be singleton condensed nodes
        assert_eq!(decomp.condensed_dim, 4);
        assert_eq!(decomp.singletons, 4);
        assert_eq!(decomp.two_cycles, 0);

        // All have valid condensed indices
        for i in 0..4 {
            assert!(decomp.old_to_new[i] < 4);
        }

        // Identity condensed ordering should produce valid expansion
        let condensed_order = vec![0i32, 1, 2, 3];
        let perm = expand_ordering(&condensed_order, &decomp, 4);
        let (fwd_expanded, _) = perm.as_ref().arrays();
        assert_eq!(fwd_expanded.len(), 4);

        let mut seen = [false; 4];
        for &v in fwd_expanded {
            seen[v] = true;
        }
        assert!(seen.iter().all(|&s| s));
    }
}
