//! Multifrontal numeric factorization for sparse symmetric indefinite matrices.
//!
//! Implements the multifrontal method with A Posteriori Threshold Pivoting (APTP)
//! for computing P^T A P = L D L^T factorizations of sparse symmetric indefinite
//! matrices. The factorization traverses the assembly tree in postorder, assembling
//! dense frontal matrices and factoring them using Phase 5's [`aptp_factor_in_place`]
//! kernel.
//!
//! # Algorithm Overview
//!
//! For each supernode in assembly-tree postorder:
//! 1. **Assemble** a dense frontal matrix by scattering original sparse entries
//!    and extend-adding child contribution blocks
//! 2. **Factor** the fully-summed portion using APTP (1x1/2x2 pivoting with delays)
//! 3. **Extract** per-supernode factors (L11, D11, L21) and contribution block
//! 4. **Propagate** contribution block (including delayed columns) to parent
//!
//! # Key Types
//!
//! - [`AptpNumeric`] — complete factorization result (public API)
//! - [`FrontFactors`] — per-supernode factors (public, used by Phase 7 solve)
//! - [`FactorizationStats`] — aggregate statistics (public)
//! - `SupernodeInfo` — unified supernode descriptor (internal)
//! - `FrontalMatrix` — temporary dense assembly matrix (internal)
//! - `ContributionBlock` — Schur complement passed to parent (internal)
//!
//! # References
//!
//! - Duff & Reid (1983), "The multifrontal solution of indefinite sparse symmetric
//!   linear equations" — foundational multifrontal method
//! - Liu (1992), "The Multifrontal Method for Sparse Matrix Solution: Theory and
//!   Practice" — frontal matrices, elimination trees, extend-add
//! - Hogg, Duff & Lopez (2020), "A New Sparse LDL^T Solver Using A Posteriori
//!   Threshold Pivoting" — APTP algorithm integration with multifrontal framework

use faer::Mat;
use faer::sparse::SparseColMat;
use faer::sparse::linalg::cholesky::SymbolicCholeskyRaw;

use super::diagonal::MixedDiagonal;
use super::factor::{AptpFactorResult, AptpOptions, aptp_factor_in_place};
use super::pivot::{Block2x2, PivotType};
use super::symbolic::AptpSymbolic;
use crate::error::SparseError;

// ---------------------------------------------------------------------------
// Internal types
// ---------------------------------------------------------------------------

/// Precomputed supernode descriptor unifying faer's supernodal and simplicial
/// decompositions.
///
/// Built once from [`AptpSymbolic`] before the factorization loop via
/// [`build_supernode_info`]. For supernodal decompositions, maps directly to
/// faer's supernode structure. For simplicial decompositions, each column
/// becomes a trivial 1-column supernode.
///
/// # Invariants
///
/// - `col_end > col_begin`
/// - `parent.map_or(true, |p| p > s)` for supernode index `s` (postorder)
///
/// # References
///
/// - Liu (1992), Section 3: supernode definitions and assembly trees
pub(crate) struct SupernodeInfo {
    /// First column index of this supernode (inclusive).
    pub col_begin: usize,
    /// Past-the-end column index (exclusive).
    pub col_end: usize,
    /// Off-diagonal row indices (global permuted column indices of non-fully-summed rows).
    pub pattern: Vec<usize>,
    /// Parent supernode index in assembly tree, or `None` for root.
    pub parent: Option<usize>,
}

// No methods needed — fields accessed directly within this module.

/// Temporary dense matrix for assembling and factoring one supernode.
///
/// Allocated, populated (scatter + extend-add), factored, and deallocated within
/// a single iteration of the factorization loop.
///
/// The matrix is partitioned into:
/// - `F11 = data[0..k, 0..k]` — fully-summed block (factored by APTP)
/// - `F21 = data[k..m, 0..k]` — subdiagonal block
/// - `F22 = data[k..m, k..m]` — contribution block (Schur complement)
///
/// where `k = num_fully_summed` and `m = data.nrows()`.
pub(crate) struct FrontalMatrix {
    /// Dense m × m storage (lower triangle used).
    pub data: Mat<f64>,
    /// Global permuted column indices for each local row (length m).
    pub row_indices: Vec<usize>,
    /// Number of fully-summed rows/columns (supernode cols + delayed from children).
    pub num_fully_summed: usize,
}

/// Schur complement and delayed columns from a factored supernode.
///
/// Created by [`extract_contribution`] and consumed by [`extend_add`] into
/// the parent supernode's frontal matrix.
///
/// # Structure
///
/// - Positions `0..num_delayed`: delayed columns (become fully-summed at parent)
/// - Positions `num_delayed..size`: non-fully-summed rows with Schur complement
pub(crate) struct ContributionBlock {
    /// Dense (m - ne) × (m - ne) trailing submatrix from factored frontal matrix.
    pub data: Mat<f64>,
    /// Global permuted column indices for rows/columns of the contribution.
    pub row_indices: Vec<usize>,
    /// Number of delayed columns at the start of this block.
    pub num_delayed: usize,
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Per-supernode factorization result.
///
/// Stores the L11, D11, and L21 blocks extracted from the factored frontal
/// matrix, along with permutation and index information needed by Phase 7
/// (triangular solve).
///
/// # Invariants
///
/// - `l11.nrows() == l11.ncols() == num_eliminated`
/// - `d11.dimension() == num_eliminated`
/// - `l21.ncols() == num_eliminated`
/// - `l21.nrows() == row_indices.len()`
/// - `col_indices.len() == num_eliminated`
#[derive(Debug)]
pub struct FrontFactors {
    /// Unit lower triangular factor (ne × ne).
    l11: Mat<f64>,
    /// Block diagonal with 1x1 and 2x2 pivots (ne entries).
    d11: MixedDiagonal,
    /// Subdiagonal factor block (r × ne).
    l21: Mat<f64>,
    /// APTP local pivot permutation (length k). Maps factored position to
    /// original front-local column.
    local_perm: Vec<usize>,
    /// Number of columns successfully eliminated (ne ≤ k).
    num_eliminated: usize,
    /// Global permuted column indices for the eliminated columns (length ne).
    col_indices: Vec<usize>,
    /// Global permuted row indices for L21 rows (length r).
    row_indices: Vec<usize>,
}

impl FrontFactors {
    /// Unit lower triangular factor L11 (ne × ne).
    pub fn l11(&self) -> &Mat<f64> {
        &self.l11
    }

    /// Block diagonal D11 with mixed 1x1/2x2 pivots.
    pub fn d11(&self) -> &MixedDiagonal {
        &self.d11
    }

    /// Subdiagonal factor block L21 (r × ne).
    pub fn l21(&self) -> &Mat<f64> {
        &self.l21
    }

    /// APTP local pivot permutation within this front.
    pub fn local_perm(&self) -> &[usize] {
        &self.local_perm
    }

    /// Number of successfully eliminated columns.
    pub fn num_eliminated(&self) -> usize {
        self.num_eliminated
    }

    /// Global permuted column indices for eliminated columns.
    pub fn col_indices(&self) -> &[usize] {
        &self.col_indices
    }

    /// Global permuted row indices for L21 rows.
    pub fn row_indices(&self) -> &[usize] {
        &self.row_indices
    }
}

/// Aggregate statistics from multifrontal factorization.
///
/// Reports pivot counts, delay events, and front sizes across all supernodes.
#[derive(Debug, Clone)]
pub struct FactorizationStats {
    /// Total 1x1 pivots across all supernodes.
    pub total_1x1_pivots: usize,
    /// Total 2x2 pivot pairs across all supernodes.
    pub total_2x2_pivots: usize,
    /// Total delay events across all supernodes.
    pub total_delayed: usize,
    /// Columns that could not be eliminated at root supernodes (zero pivots).
    /// These represent rank deficiency in the matrix. A nonzero value means
    /// the matrix is numerically singular (or nearly so) — the solve phase
    /// must handle this appropriately.
    pub zero_pivots: usize,
    /// Largest frontal matrix dimension encountered.
    pub max_front_size: usize,
}

/// Complete multifrontal numeric factorization result.
///
/// Contains per-supernode factors and aggregate statistics. Created by
/// [`AptpNumeric::factor`] and used by Phase 7 (triangular solve) to
/// perform forward/backward substitution.
///
/// # Usage
///
/// ```no_run
/// use faer::sparse::{SparseColMat, Triplet};
/// use faer::sparse::linalg::cholesky::SymmetricOrdering;
/// use rivrs_sparse::aptp::{AptpSymbolic, AptpOptions, AptpNumeric};
///
/// // Analyze
/// # let triplets = vec![Triplet::new(0, 0, 1.0)];
/// # let matrix = SparseColMat::try_new_from_triplets(1, 1, &triplets).unwrap();
/// let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd).unwrap();
///
/// // Factor
/// let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default()).unwrap();
/// println!("Stats: {:?}", numeric.stats());
/// ```
#[derive(Debug)]
pub struct AptpNumeric {
    /// Per-supernode factors, indexed by supernode ID.
    front_factors: Vec<FrontFactors>,
    /// Aggregate factorization statistics.
    stats: FactorizationStats,
    /// Matrix dimension.
    n: usize,
}

impl AptpNumeric {
    /// Per-supernode factors, indexed by supernode ID.
    pub fn front_factors(&self) -> &[FrontFactors] {
        &self.front_factors
    }

    /// Aggregate factorization statistics.
    pub fn stats(&self) -> &FactorizationStats {
        &self.stats
    }

    /// Matrix dimension.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Factor a sparse symmetric matrix using the multifrontal method with APTP.
    ///
    /// Traverses the assembly tree in postorder, assembling and factoring dense
    /// frontal matrices at each supernode using Phase 5's APTP kernel.
    ///
    /// # Arguments
    ///
    /// - `symbolic`: Symbolic analysis result from [`AptpSymbolic::analyze`]
    /// - `matrix`: Sparse symmetric matrix (lower triangle stored). Dimensions
    ///   must match `symbolic.nrows()`
    /// - `options`: APTP configuration (threshold, fallback strategy)
    ///
    /// # Errors
    ///
    /// - [`SparseError::DimensionMismatch`] if matrix dimensions don't match symbolic
    /// - [`SparseError::AnalysisFailure`] if symbolic analysis is inconsistent
    ///
    /// # Zero Pivots
    ///
    /// If a root supernode has columns that cannot be eliminated (all pivots
    /// delayed to root but still fail the threshold), these are recorded as
    /// zero pivots in [`FactorizationStats::zero_pivots`] rather than returning
    /// an error. This matches SPRAL's behavior: the factorization succeeds but
    /// the solve phase must handle the rank deficiency.
    ///
    /// # References
    ///
    /// - Duff & Reid (1983), "The multifrontal solution of indefinite sparse
    ///   symmetric linear equations"
    /// - Liu (1992), "The Multifrontal Method for Sparse Matrix Solution"
    /// - Hogg, Duff & Lopez (2020), "A New Sparse LDL^T Solver Using A Posteriori
    ///   Threshold Pivoting"
    pub fn factor(
        symbolic: &AptpSymbolic,
        matrix: &SparseColMat<usize, f64>,
        options: &AptpOptions,
    ) -> Result<Self, SparseError> {
        let n = symbolic.nrows();

        // Validate dimensions
        if matrix.nrows() != n || matrix.ncols() != n {
            return Err(SparseError::DimensionMismatch {
                expected: (n, n),
                got: (matrix.nrows(), matrix.ncols()),
                context: "Matrix dimensions must match symbolic analysis".to_string(),
            });
        }

        // Build supernode info (unified abstraction over supernodal/simplicial)
        let supernodes = build_supernode_info(symbolic);
        let n_supernodes = supernodes.len();
        let children = build_children_map(&supernodes);

        // Get fill-reducing permutation
        let (perm_fwd, perm_inv) = if let Some(perm) = symbolic.perm() {
            let (fwd, inv) = perm.arrays();
            (fwd.to_vec(), inv.to_vec())
        } else {
            let id: Vec<usize> = (0..n).collect();
            (id.clone(), id)
        };

        // Allocate shared data structures
        let mut global_to_local = vec![usize::MAX; n];
        let mut contributions: Vec<Option<ContributionBlock>> =
            (0..n_supernodes).map(|_| None).collect();
        let mut front_factors_vec = Vec::with_capacity(n_supernodes);
        let mut stats = FactorizationStats {
            total_1x1_pivots: 0,
            total_2x2_pivots: 0,
            total_delayed: 0,
            zero_pivots: 0,
            max_front_size: 0,
        };

        // Postorder traversal (index order = postorder for faer's supernodal layout)
        for s in 0..n_supernodes {
            let sn = &supernodes[s];

            // 1. Collect delayed columns from children
            let mut delayed_cols: Vec<usize> = Vec::new();
            for &c in &children[s] {
                if let Some(ref cb) = contributions[c] {
                    delayed_cols.extend_from_slice(&cb.row_indices[..cb.num_delayed]);
                }
            }

            // 2. Compute frontal matrix structure
            let sn_cols: Vec<usize> = (sn.col_begin..sn.col_end).collect();
            let k = sn_cols.len() + delayed_cols.len(); // num_fully_summed

            // Row indices: fully-summed (sn columns + delayed) + off-diagonal pattern
            let mut frontal_rows: Vec<usize> = Vec::with_capacity(k + sn.pattern.len());
            frontal_rows.extend_from_slice(&sn_cols);
            frontal_rows.extend_from_slice(&delayed_cols);
            frontal_rows.extend_from_slice(&sn.pattern);
            let m = frontal_rows.len();

            // Track max front size
            if m > stats.max_front_size {
                stats.max_front_size = m;
            }

            // 3. Build global-to-local mapping
            for (i, &global) in frontal_rows.iter().enumerate() {
                global_to_local[global] = i;
            }

            // 4. Allocate and assemble frontal matrix
            let mut frontal = FrontalMatrix {
                data: Mat::zeros(m, m),
                row_indices: frontal_rows.clone(),
                num_fully_summed: k,
            };

            // Scatter original sparse entries (only for supernode's own columns)
            scatter_original_entries(
                &mut frontal,
                matrix,
                &perm_fwd,
                &perm_inv,
                &global_to_local,
                sn.col_begin,
                sn.col_end,
            );

            // Extend-add child contributions
            for &c in &children[s] {
                if let Some(cb) = contributions[c].take() {
                    extend_add(&mut frontal, &cb, &global_to_local);
                }
            }

            // 5. Factor the frontal matrix
            let result = aptp_factor_in_place(frontal.data.as_mut(), k, options)?;
            let ne = result.num_eliminated;

            // 6. Accumulate statistics
            stats.total_1x1_pivots += result.stats.num_1x1;
            stats.total_2x2_pivots += result.stats.num_2x2;
            stats.total_delayed += result.stats.num_delayed;

            // 7. Extract and store front factors
            let ff = extract_front_factors(&frontal, &result);
            front_factors_vec.push(ff);

            // 8. Prepare contribution for parent (if not root and not fully eliminated)
            if sn.parent.is_some() && ne < m {
                contributions[s] = Some(extract_contribution(&frontal, &result));
            } else if sn.parent.is_none() && ne < m {
                // Root supernode with unresolved delayed columns — record as zero
                // pivots rather than erroring. This matches SPRAL's behavior:
                // the factorization succeeds but the matrix is rank-deficient.
                // The solve phase must handle zero pivots appropriately.
                let n_unresolved = k - ne;
                stats.zero_pivots += n_unresolved;
            }

            // 9. Cleanup global-to-local mapping
            for &global in &frontal_rows {
                global_to_local[global] = usize::MAX;
            }
        }

        Ok(AptpNumeric {
            front_factors: front_factors_vec,
            stats,
            n,
        })
    }
}

// ---------------------------------------------------------------------------
// Internal functions
// ---------------------------------------------------------------------------

/// Build unified supernode descriptors from symbolic analysis.
///
/// For supernodal decompositions, maps directly from faer's supernode structure.
/// For simplicial decompositions, each column becomes a trivial 1-column supernode
/// with pattern derived from the elimination tree and L structure.
///
/// # Postcondition
///
/// `info[s].parent.map_or(true, |p| p > s)` for all `s` (postorder).
pub(crate) fn build_supernode_info(symbolic: &AptpSymbolic) -> Vec<SupernodeInfo> {
    match symbolic.raw() {
        SymbolicCholeskyRaw::Supernodal(sn) => {
            let ns = sn.n_supernodes();
            let begin = sn.supernode_begin();
            let end = sn.supernode_end();
            (0..ns)
                .map(|s| {
                    let pattern = sn.supernode(s).pattern().to_vec();
                    let parent = symbolic.supernode_parent(s);
                    SupernodeInfo {
                        col_begin: begin[s],
                        col_end: end[s],
                        pattern,
                        parent,
                    }
                })
                .collect()
        }
        SymbolicCholeskyRaw::Simplicial(simp) => {
            let n = symbolic.nrows();
            let etree = symbolic.etree();
            // Get L structure from simplicial symbolic factorization
            let l_symbolic = simp.factor();
            let col_ptr = l_symbolic.col_ptr();
            let row_idx = l_symbolic.row_idx();
            (0..n)
                .map(|j| {
                    // Pattern = row indices in column j of L that are > j
                    // (off-diagonal entries below the diagonal)
                    let start = col_ptr[j];
                    let end = col_ptr[j + 1];
                    let pattern: Vec<usize> = row_idx[start..end]
                        .iter()
                        .copied()
                        .filter(|&r| r > j)
                        .collect();
                    // Parent from etree
                    let parent = if etree[j] < 0 {
                        None
                    } else {
                        Some(etree[j] as usize)
                    };
                    SupernodeInfo {
                        col_begin: j,
                        col_end: j + 1,
                        pattern,
                        parent,
                    }
                })
                .collect()
        }
    }
}

/// Build a map from each supernode to its children in the assembly tree.
pub(crate) fn build_children_map(infos: &[SupernodeInfo]) -> Vec<Vec<usize>> {
    let n = infos.len();
    let mut children = vec![Vec::new(); n];
    for (s, info) in infos.iter().enumerate() {
        if let Some(p) = info.parent {
            children[p].push(s);
        }
    }
    children
}

/// Scatter original sparse matrix entries into a frontal matrix.
///
/// Only scatters entries from the supernode's own columns `[col_begin, col_end)`.
/// Delayed columns get their values from child contribution blocks via
/// extend-add, not from scatter.
///
/// # Double-counting avoidance
///
/// Two mechanisms prevent double-counting:
///
/// 1. **Between two supernode columns**: For a full-stored symmetric CSC matrix,
///    each off-diagonal entry (i,j) appears in both column i and column j. We
///    skip the upper-triangle instance (`orig_row < orig_col`) when both
///    endpoints are supernode columns, counting only the lower-triangle instance.
///
/// 2. **Between supernode column and delayed column**: Entries connecting a
///    supernode's own column to a delayed column from a child are already in
///    the child's contribution block (as part of F21 or the unfactored frontal
///    data). Scatter skips these to avoid adding them twice.
///
/// # References
///
/// - Liu (1992), Section 4: original entries are scattered only for
///   fully-summed columns, not for update (F22) entries
fn scatter_original_entries(
    frontal: &mut FrontalMatrix,
    matrix: &SparseColMat<usize, f64>,
    perm_fwd: &[usize],
    perm_inv: &[usize],
    global_to_local: &[usize],
    col_begin: usize,
    col_end: usize,
) {
    let symbolic = matrix.symbolic();
    let col_ptrs = symbolic.col_ptr();
    let row_indices_csc = symbolic.row_idx();
    let values = matrix.val();
    let sn_ncols = col_end - col_begin;
    let k = frontal.num_fully_summed;

    // Only iterate over the supernode's own columns (permuted indices)
    for pj in col_begin..col_end {
        let orig_col = perm_fwd[pj];
        let local_col = global_to_local[pj];

        // Iterate nonzeros in this column of the original matrix
        let start = col_ptrs[orig_col];
        let end = col_ptrs[orig_col + 1];
        for idx in start..end {
            let orig_row = row_indices_csc[idx];
            // Skip upper-triangle entries ONLY if both endpoints are supernode
            // columns (meaning we'll see this entry from the other column too).
            // If orig_row maps to a non-supernode row (pattern or other front),
            // we must include it here — it won't be seen from the other side.
            if orig_row < orig_col {
                let perm_row = perm_inv[orig_row];
                if perm_row >= col_begin && perm_row < col_end {
                    continue; // both endpoints are supernode cols; avoid double-count
                }
            }
            let global_row = perm_inv[orig_row];
            let local_row = global_to_local[global_row];
            if local_row == usize::MAX {
                continue; // not in this front
            }
            // Skip entries where the row is a delayed column from a child.
            // These entries are already in the child's contribution block
            // and will be placed by extend-add.
            if local_row >= sn_ncols && local_row < k {
                continue;
            }
            let val = values[idx];
            // Place in lower triangle of frontal matrix
            if local_row >= local_col {
                frontal.data[(local_row, local_col)] += val;
            } else {
                frontal.data[(local_col, local_row)] += val;
            }
        }
    }
}

/// Perform extend-add: merge a child's contribution block into the parent frontal matrix.
///
/// For each lower-triangle entry in the child contribution, maps global indices
/// to parent local positions and adds the value.
///
/// # References
///
/// - Liu (1992), Section 4.2: extend-add operator
pub(crate) fn extend_add(
    parent: &mut FrontalMatrix,
    child: &ContributionBlock,
    global_to_local: &[usize],
) {
    let cb_size = child.data.nrows();
    for i in 0..cb_size {
        let gi = child.row_indices[i];
        let li = global_to_local[gi];
        debug_assert!(
            li != usize::MAX,
            "extend_add: child row {} not in parent",
            gi
        );
        for j in 0..=i {
            // Lower triangle only
            let gj = child.row_indices[j];
            let lj = global_to_local[gj];
            debug_assert!(
                lj != usize::MAX,
                "extend_add: child col {} not in parent",
                gj
            );
            let val = child.data[(i, j)];
            if val != 0.0 {
                // Map to parent: ensure lower triangle (li >= lj)
                if li >= lj {
                    parent.data[(li, lj)] += val;
                } else {
                    parent.data[(lj, li)] += val;
                }
            }
        }
    }
}

/// Extract per-supernode factors from a factored frontal matrix.
///
/// Copies L11, D11, L21, and permutation information from the in-place
/// factored frontal matrix into a persistent [`FrontFactors`] struct.
pub(crate) fn extract_front_factors(
    frontal: &FrontalMatrix,
    result: &AptpFactorResult,
) -> FrontFactors {
    let m = frontal.data.nrows();
    let k = frontal.num_fully_summed;
    let ne = result.num_eliminated;

    // Extract L11 (ne × ne) — unit lower triangular part
    // For 2x2 pivots, a[(col+1, col)] is the D off-diagonal, NOT an L entry.
    let l11 = if ne > 0 {
        let mut l = Mat::zeros(ne, ne);
        let mut col = 0;
        while col < ne {
            l[(col, col)] = 1.0; // unit diagonal
            match result.d.get_pivot_type(col) {
                PivotType::OneByOne => {
                    // L entries start at row col+1
                    for i in (col + 1)..ne {
                        l[(i, col)] = frontal.data[(i, col)];
                    }
                    col += 1;
                }
                PivotType::TwoByTwo { partner } if partner > col => {
                    l[(col + 1, col + 1)] = 1.0; // unit diagonal for partner
                    // a[(col+1, col)] is D off-diagonal, skip it
                    // L entries for both columns start at row col+2
                    for i in (col + 2)..ne {
                        l[(i, col)] = frontal.data[(i, col)];
                        l[(i, col + 1)] = frontal.data[(i, col + 1)];
                    }
                    col += 2;
                }
                PivotType::TwoByTwo { .. } => {
                    // Second column of a 2x2 pair — already handled above
                    col += 1;
                }
                PivotType::Delayed => {
                    // Should not appear for columns 0..ne (they were eliminated)
                    debug_assert!(false, "unexpected Delayed pivot at col {} in 0..ne", col);
                    col += 1;
                }
            }
        }
        l
    } else {
        Mat::zeros(0, 0)
    };

    // Build truncated D11 (first ne pivots from result.d)
    let mut d11 = MixedDiagonal::new(ne);
    let mut col = 0;
    while col < ne {
        match result.d.get_pivot_type(col) {
            PivotType::OneByOne => {
                d11.set_1x1(col, result.d.get_1x1(col));
                col += 1;
            }
            PivotType::TwoByTwo { partner: _ } => {
                let block = result.d.get_2x2(col);
                d11.set_2x2(Block2x2 {
                    first_col: col,
                    a: block.a,
                    b: block.b,
                    c: block.c,
                });
                col += 2;
            }
            PivotType::Delayed => {
                // Should not happen for columns 0..ne (they were eliminated)
                col += 1;
            }
        }
    }

    // Extract L21 (r × ne) where r = m - k
    let r = m - k;
    let l21 = if ne > 0 && r > 0 {
        let mut l = Mat::zeros(r, ne);
        for i in 0..r {
            for j in 0..ne {
                l[(i, j)] = frontal.data[(k + i, j)];
            }
        }
        l
    } else {
        Mat::zeros(r, ne)
    };

    // Local permutation (maps factored position to original front-local column)
    let local_perm = result.perm[..k].to_vec();

    // Column indices: the global permuted indices of the eliminated columns
    // local_perm[0..ne] gives the front-local columns that were eliminated,
    // mapped through frontal.row_indices to get global permuted indices
    let col_indices: Vec<usize> = local_perm[..ne]
        .iter()
        .map(|&lp| frontal.row_indices[lp])
        .collect();

    // Row indices: global permuted indices for L21 rows
    let row_indices = frontal.row_indices[k..].to_vec();

    FrontFactors {
        l11,
        d11,
        l21,
        local_perm,
        num_eliminated: ne,
        col_indices,
        row_indices,
    }
}

/// Extract the contribution block from a factored frontal matrix.
///
/// The contribution block is the trailing `(m - ne) × (m - ne)` submatrix
/// of the factored frontal matrix, containing:
/// - Delayed columns (positions 0..num_delayed)
/// - Schur complement entries (positions num_delayed..size)
pub(crate) fn extract_contribution(
    frontal: &FrontalMatrix,
    result: &AptpFactorResult,
) -> ContributionBlock {
    let m = frontal.data.nrows();
    let ne = result.num_eliminated;
    let k = frontal.num_fully_summed;
    let size = m - ne;
    let num_delayed = k - ne;

    // Copy trailing submatrix
    let mut data = Mat::zeros(size, size);
    for i in 0..size {
        for j in 0..=i {
            data[(i, j)] = frontal.data[(ne + i, ne + j)];
        }
    }

    // Build row indices:
    // - First num_delayed entries: delayed columns mapped through result.perm and frontal.row_indices
    // - Remaining entries: non-fully-summed rows from frontal.row_indices[k..m]
    let mut row_indices = Vec::with_capacity(size);

    // Delayed columns: these are at positions ne..k in the APTP-permuted frontal matrix
    // result.perm[ne..k] gives the original front-local indices of delayed columns
    for &lp in &result.perm[ne..k] {
        row_indices.push(frontal.row_indices[lp]);
    }

    // Non-fully-summed rows
    row_indices.extend_from_slice(&frontal.row_indices[k..]);

    ContributionBlock {
        data,
        row_indices,
        num_delayed,
    }
}
