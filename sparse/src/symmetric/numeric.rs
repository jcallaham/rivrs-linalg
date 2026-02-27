//! Multifrontal numeric factorization for sparse symmetric indefinite matrices.
//!
//! Implements the multifrontal method with A Posteriori Threshold Pivoting (APTP)
//! for computing P^T A P = L D L^T factorizations of sparse symmetric indefinite
//! matrices. The factorization traverses the assembly tree in postorder, assembling
//! dense frontal matrices and factoring them using the [`aptp_factor_in_place`]
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
//! - [`FrontFactors`] — per-supernode factors (public, used by triangular solve)
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

use faer::sparse::SparseColMat;
use faer::sparse::linalg::cholesky::SymbolicCholeskyRaw;
use faer::{Mat, MatMut, Par};

use super::diagonal::{MixedDiagonal, PivotEntry};
use super::factor::{
    AptpFactorResult, AptpOptions, aptp_factor_in_place, compute_contribution_gemm,
    compute_contribution_gemm_rect, extract_contribution_rect, extract_front_factors_rect,
    tpp_factor_rect,
};
use super::pivot::{Block2x2, PivotType};
use super::symbolic::AptpSymbolic;
use crate::error::SparseError;

#[cfg(feature = "diagnostic")]
use crate::profiling::{FinishedSession, ProfileSession};

/// Sentinel value indicating a global column is not part of the current front.
const NOT_IN_FRONT: usize = usize::MAX;

// ---------------------------------------------------------------------------
// Workspace types
// ---------------------------------------------------------------------------

/// Pre-allocated reusable buffers for the per-supernode factorization loop.
///
/// Eliminates per-supernode heap allocations for the frontal matrix, row indices,
/// delayed columns buffer, and global-to-local mapping. Sized to the maximum
/// front dimension from symbolic analysis and reused across supernodes.
///
/// # Invariant
///
/// Never shared between concurrent supernodes. In sequential mode, passed as
/// `&mut`. In parallel mode, each rayon worker has its own instance via
/// `Cell`-based move semantics.
///
/// # References
///
/// - Duff, Hogg & Lopez (2020), Section 5: workspace reuse strategy
/// - Liu (1992), Section 4: multifrontal workspace management
pub(crate) struct FactorizationWorkspace {
    /// Dense buffer, max_front × max_front. Reused as each supernode's frontal matrix.
    frontal_data: Mat<f64>,
    /// Buffer for frontal matrix row indices, capacity max_front.
    frontal_row_indices: Vec<usize>,
    /// Buffer for collecting delayed columns from children.
    delayed_cols_buf: Vec<usize>,
    /// Global-to-local row index mapping, length n (matrix dimension).
    global_to_local: Vec<usize>,
    /// Pre-allocated contribution buffer, max_front × max_front. Receives the
    /// deferred contribution GEMM output. Reused across supernodes via swap protocol:
    /// final GEMM writes here → buffer moved into ContributionBlock → parent's
    /// extend-add recycles buffer back.
    contrib_buffer: Mat<f64>,
    /// Pre-allocated L·D product workspace for `compute_contribution_gemm`.
    /// Sized nfs × ne (grows on demand). Eliminates per-supernode allocation
    /// in the hot loop.
    ld_workspace: Mat<f64>,
}

impl Default for FactorizationWorkspace {
    fn default() -> Self {
        Self::empty()
    }
}

impl FactorizationWorkspace {
    /// Create a new workspace sized for the given maximum front dimension.
    ///
    /// # Arguments
    ///
    /// - `max_front`: Maximum frontal matrix dimension across all supernodes
    /// - `n`: Matrix dimension (for global-to-local mapping)
    pub(crate) fn new(max_front: usize, n: usize) -> Self {
        if max_front == 0 {
            return Self {
                frontal_data: Mat::zeros(0, 0),
                frontal_row_indices: Vec::new(),
                delayed_cols_buf: Vec::new(),
                global_to_local: vec![NOT_IN_FRONT; n],
                contrib_buffer: Mat::zeros(0, 0),
                ld_workspace: Mat::new(),
            };
        }
        Self {
            frontal_data: Mat::zeros(max_front, max_front),
            frontal_row_indices: Vec::with_capacity(max_front),
            delayed_cols_buf: Vec::with_capacity(max_front),
            global_to_local: vec![NOT_IN_FRONT; n],
            // Lazily allocated on first use in factor_single_supernode.
            // Recycled via extend_add buffer return, so typically allocates once.
            contrib_buffer: Mat::new(),
            ld_workspace: Mat::new(),
        }
    }

    /// Empty workspace for `Cell` default in thread-local storage.
    const fn empty() -> Self {
        Self {
            frontal_data: Mat::new(),
            frontal_row_indices: Vec::new(),
            delayed_cols_buf: Vec::new(),
            global_to_local: Vec::new(),
            contrib_buffer: Mat::new(),
            ld_workspace: Mat::new(),
        }
    }

    /// Zero the lower triangle of the m × m subregion of `frontal_data`.
    ///
    /// Called after the front size `m` is known and row indices are populated.
    /// Does NOT clear `frontal_row_indices` or `delayed_cols_buf` — those
    /// are managed by the caller before determining `m`.
    ///
    /// # Panics
    ///
    /// Panics if `m > self.frontal_data.nrows()`.
    fn zero_frontal(&mut self, m: usize) {
        assert!(
            m <= self.frontal_data.nrows(),
            "front size {} exceeds workspace capacity {}",
            m,
            self.frontal_data.nrows()
        );
        // Zero the lower triangle of the m × m subregion via per-column fill.
        // col_as_slice_mut gives a contiguous &mut [f64] for each column of the
        // owned Mat, so fill(0.0) compiles to memset — much faster than indexed writes.
        for j in 0..m {
            self.frontal_data.col_as_slice_mut(j)[j..m].fill(0.0);
        }
    }

    /// Ensure `global_to_local` mapping is at least `n` entries.
    ///
    /// Used by the parallel path where frontal_data grows on demand in
    /// `factor_single_supernode` (avoiding eager max_front² allocation
    /// per thread that causes OOM on large-front matrices like H2O).
    fn ensure_g2l(&mut self, n: usize) {
        if self.global_to_local.len() < n {
            self.global_to_local.resize(n, NOT_IN_FRONT);
        }
    }
}

/// Front dimension below which intra-node BLAS uses `Par::Seq` regardless of
/// the user-supplied parallelism setting. Fronts smaller than this threshold
/// do not benefit from parallel BLAS (overhead exceeds computation).
pub(crate) const INTRA_NODE_THRESHOLD: usize = 256;

/// Compute an upper bound on the maximum front size across all supernodes.
///
/// The front size for a supernode is `owned_cols + delayed_cols + pattern.len()`.
/// Since delayed columns depend on the actual factorization, we estimate using
/// `owned_cols + pattern.len()` (front size without delays) and add a safety
/// factor. If delayed columns cause a front to exceed this estimate, the
/// workspace will be resized at that point.
fn estimate_max_front_size(supernodes: &[SupernodeInfo]) -> usize {
    supernodes
        .iter()
        .map(|sn| {
            let owned: usize = sn.owned_ranges.iter().map(|r| r.len()).sum();
            owned + sn.pattern.len()
        })
        .max()
        .unwrap_or(0)
}

// ---------------------------------------------------------------------------
// Assembly maps (precomputed scatter indices)
// ---------------------------------------------------------------------------

/// Precomputed index mappings for assembly operations.
///
/// Eliminates per-entry index arithmetic in `scatter_original_entries_multi`
/// and `extend_add` by storing source→destination index pairs computed from
/// the symbolic structure. Computed once after amalgamation at the start of
/// factorization and read-only during the assembly loop.
///
/// Uses `u32` indices to halve memory (frontal matrices are bounded by
/// max_front_size, well within u32 range).
///
/// # Layout
///
/// Both `amap_entries` and `ea_map` use CSC-style compressed storage with
/// per-supernode/per-child offset arrays for O(1) slice access.
///
/// # References
///
/// - SPRAL `assemble.hxx:51-80` (`add_a_block`): amap-based assembly pattern
pub(crate) struct AssemblyMaps {
    /// Flattened tuples for original matrix scatter. Each entry is 4 u32 values:
    /// `[src_csc_index, dest_frontal_linear_index, scale_row, scale_col]` where
    /// `src_csc_index` is the position in the CSC values array, `dest_frontal_linear_index`
    /// is `col * m + row` in column-major layout, and `scale_row`/`scale_col` are
    /// the permuted indices for MC64 scaling lookup (`scaling[scale_row] * scaling[scale_col]`).
    /// Length: 4 * total_entries across all supernodes.
    pub amap_entries: Vec<u32>,
    /// Per-supernode start offset into `amap_entries` (in entry units; each entry
    /// is 4 u32 values). `amap_entries[amap_offsets[s]*4 .. amap_offsets[s+1]*4]`
    /// gives the flattened entries for supernode s.
    /// Length: num_supernodes + 1.
    pub amap_offsets: Vec<usize>,

    /// Flattened child→parent row index mappings for extend-add.
    /// For each child's contribution row, stores the local row index in the
    /// parent's frontal matrix (assuming zero delays).
    pub ea_map: Vec<u32>,
    /// Per-child start offset into `ea_map`. Indexed by the flattened child
    /// enumeration (all children of all supernodes in tree order).
    /// Length: total_children + 1.
    pub ea_offsets: Vec<usize>,
    /// Per-supernode start offset into the children enumeration.
    /// `ea_snode_child_begin[s]..ea_snode_child_begin[s+1]` gives the range
    /// in `ea_offsets` for supernode s's children.
    /// Length: num_supernodes + 1.
    pub ea_snode_child_begin: Vec<usize>,
}

/// Build precomputed assembly maps from post-amalgamation supernode structure.
///
/// # Arguments
///
/// - `supernodes`: Post-amalgamation supernode descriptors.
/// - `children`: Children map (children[s] = list of child supernode indices).
/// - `matrix`: Original sparse matrix (for CSC structure).
/// - `perm_fwd`: Forward permutation (perm_fwd[new] = old).
/// - `perm_inv`: Inverse permutation (perm_inv[old] = new).
/// - `n`: Matrix dimension.
fn build_assembly_maps(
    supernodes: &[SupernodeInfo],
    children: &[Vec<usize>],
    matrix: &SparseColMat<usize, f64>,
    perm_fwd: &[usize],
    perm_inv: &[usize],
    n: usize,
) -> AssemblyMaps {
    let n_supernodes = supernodes.len();
    let symbolic = matrix.symbolic();
    let col_ptrs = symbolic.col_ptr();
    let row_indices_csc = symbolic.row_idx();

    // ---- Build amap (original entry scatter maps) ----
    // Each entry is 4 u32: [src_csc_index, dest_frontal_linear, scale_row, scale_col]
    // where scale_row = perm_inv[orig_row], scale_col = perm_inv[orig_col]
    // (precomputed for efficient MC64 scaling at assembly time).

    let mut amap_entries = Vec::new();
    let mut amap_offsets = vec![0usize; n_supernodes + 1];

    // Temporary global-to-local mapping, reused per supernode
    let mut g2l = vec![NOT_IN_FRONT; n];

    for (s, sn) in supernodes.iter().enumerate() {
        // Compute front structure for this supernode (without delays)
        let sn_ncols: usize = sn.owned_ranges.iter().map(|r| r.len()).sum();
        let mut frontal_rows: Vec<usize> = Vec::with_capacity(sn_ncols + sn.pattern.len());
        for range in &sn.owned_ranges {
            frontal_rows.extend(range.clone());
        }
        frontal_rows.extend_from_slice(&sn.pattern);
        let m = frontal_rows.len();

        // Build global-to-local
        for (i, &global) in frontal_rows.iter().enumerate() {
            g2l[global] = i;
        }

        let total_owned = sn_ncols;

        // Iterate owned columns and build amap entries
        for range in &sn.owned_ranges {
            for pj in range.clone() {
                let orig_col = perm_fwd[pj];
                let local_col = g2l[pj];

                let start = col_ptrs[orig_col];
                let end = col_ptrs[orig_col + 1];
                for (idx, &orig_row) in row_indices_csc.iter().enumerate().take(end).skip(start) {
                    // Upper-triangle dedup: skip if both are owned supernode cols
                    if orig_row < orig_col {
                        let perm_row = perm_inv[orig_row];
                        let local_peer = g2l[perm_row];
                        if local_peer != NOT_IN_FRONT && local_peer < total_owned {
                            continue;
                        }
                    }
                    let global_row = perm_inv[orig_row];
                    let local_row = g2l[global_row];
                    if local_row == NOT_IN_FRONT {
                        continue;
                    }

                    // Determine destination in lower triangle
                    let (dest_row, dest_col) = if local_row >= local_col {
                        (local_row, local_col)
                    } else {
                        (local_col, local_row)
                    };
                    let dest_linear = dest_col * m + dest_row;

                    amap_entries.push(idx as u32);
                    amap_entries.push(dest_linear as u32);
                    amap_entries.push(global_row as u32); // perm_inv[orig_row]
                    amap_entries.push(pj as u32); // perm_inv[orig_col] = pj
                }
            }
        }

        // Reset g2l
        for &global in &frontal_rows {
            g2l[global] = NOT_IN_FRONT;
        }

        let entry_end = amap_entries.len() / 4;
        amap_offsets[s + 1] = entry_end;
    }

    // ---- Build extend-add maps ----

    let mut ea_map = Vec::new();
    let mut ea_offsets = vec![0usize];
    let mut ea_snode_child_begin = vec![0usize; n_supernodes + 1];

    for (s, sn) in supernodes.iter().enumerate() {
        // Build parent's frontal row structure (without delays)
        let sn_ncols: usize = sn.owned_ranges.iter().map(|r| r.len()).sum();
        let mut parent_rows: Vec<usize> = Vec::with_capacity(sn_ncols + sn.pattern.len());
        for range in &sn.owned_ranges {
            parent_rows.extend(range.clone());
        }
        parent_rows.extend_from_slice(&sn.pattern);

        // Build g2l for parent
        for (i, &global) in parent_rows.iter().enumerate() {
            g2l[global] = i;
        }

        for &c in &children[s] {
            let child_sn = &supernodes[c];

            // Child's contribution rows = child's pattern (off-diagonal rows).
            // Without delays, the contribution block rows are exactly the pattern.
            // With delays, additional delayed rows appear — the precomputed map
            // won't cover those (fallback to g2l at factorization time).
            for &child_row in &child_sn.pattern {
                let parent_local = g2l[child_row];
                debug_assert!(
                    parent_local != NOT_IN_FRONT,
                    "child pattern row {} not in parent supernode {}",
                    child_row,
                    s
                );
                ea_map.push(parent_local as u32);
            }
            ea_offsets.push(ea_map.len());
        }

        // Reset g2l
        for &global in &parent_rows {
            g2l[global] = NOT_IN_FRONT;
        }

        ea_snode_child_begin[s + 1] = ea_offsets.len() - 1;
    }

    AssemblyMaps {
        amap_entries,
        amap_offsets,
        ea_map,
        ea_offsets,
        ea_snode_child_begin,
    }
}

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
#[derive(Clone)]
pub(crate) struct SupernodeInfo {
    /// First column index of this supernode (inclusive).
    pub col_begin: usize,
    /// Past-the-end column index (exclusive).
    pub col_end: usize,
    /// Off-diagonal row indices (global permuted column indices of non-fully-summed rows).
    pub pattern: Vec<usize>,
    /// Parent supernode index in assembly tree, or `None` for root.
    pub parent: Option<usize>,
    /// Column ranges actually owned by this supernode. After amalgamation of
    /// non-adjacent supernodes, `col_begin..col_end` may span columns belonging
    /// to other supernodes. `owned_ranges` tracks the true owned columns for
    /// scatter. Default: `[col_begin..col_end]`.
    pub owned_ranges: Vec<std::ops::Range<usize>>,
    /// True if this supernode belongs to a classified small-leaf subtree.
    /// Set by `classify_small_leaf_subtrees()` after amalgamation.
    pub in_small_leaf: bool,
}

/// A classified small-leaf subtree eligible for the fast-path factorization.
///
/// Produced by [`classify_small_leaf_subtrees`] after amalgamation and consumed
/// by `factor_small_leaf_subtree` during the pre-pass before the main level-set
/// loop. All supernodes in the subtree have front_size < `small_leaf_threshold`,
/// enabling a streamlined code path with a single small reusable workspace.
///
/// # Invariants
///
/// - `nodes.len() >= 2` (single-node "subtrees" are excluded)
/// - `*nodes.last().unwrap() == root`
/// - All nodes have `in_small_leaf = true`
/// - `max_front_size < small_leaf_threshold`
/// - Nodes are in postorder: children appear before parents
///
/// # References
///
/// - SPRAL `SymbolicSubtree.hxx:57-84` (BSD-3): subtree classification
/// - SPRAL `SmallLeafNumericSubtree.hxx:187-446` (BSD-3): fast-path factorization
pub(crate) struct SmallLeafSubtree {
    /// Supernode ID of the subtree root (topmost node in the subtree).
    /// Used by tests for invariant verification; always equals `*nodes.last()`.
    #[allow(dead_code)]
    pub root: usize,
    /// Supernode IDs in postorder (leaves first, root last).
    pub nodes: Vec<usize>,
    /// Maximum front size across all nodes in the subtree.
    pub max_front_size: usize,
    /// Parent supernode outside the subtree that receives the root's contribution,
    /// or `None` if the subtree root is also a tree root.
    pub parent_of_root: Option<usize>,
}

// No methods needed — fields accessed directly within this module.

/// Temporary dense matrix for assembling and factoring one supernode.
///
/// References workspace buffers during the factorization loop, or owns its
/// data in diagnostic contexts (`export_assembled_frontal`).
///
/// The matrix is partitioned into:
/// - `F11 = data[0..k, 0..k]` — fully-summed block (factored by APTP)
/// - `F21 = data[k..m, 0..k]` — subdiagonal block
/// - `F22 = data[k..m, k..m]` — contribution block (Schur complement)
///
/// where `k = num_fully_summed` and `m = data.nrows()`.
pub(crate) struct FrontalMatrix<'a> {
    /// Dense m × m storage (lower triangle used). Borrows from workspace.
    pub data: MatMut<'a, f64>,
    /// Global permuted column indices for each local row (length m).
    /// Used by callers that need index mapping (e.g., scatter maps).
    #[allow(dead_code)]
    pub row_indices: &'a [usize],
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
/// matrix, along with permutation and index information needed by the
/// triangular solve.
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
    pub(crate) l11: Mat<f64>,
    /// Block diagonal with 1x1 and 2x2 pivots (ne entries).
    pub(crate) d11: MixedDiagonal,
    /// Subdiagonal factor block (r × ne).
    pub(crate) l21: Mat<f64>,
    /// APTP local pivot permutation (length k). Maps factored position to
    /// original front-local column.
    pub(crate) local_perm: Vec<usize>,
    /// Number of columns successfully eliminated (ne ≤ k).
    pub(crate) num_eliminated: usize,
    /// Global permuted column indices for the eliminated columns (length ne).
    pub(crate) col_indices: Vec<usize>,
    /// Global permuted row indices for L21 rows (length r).
    pub(crate) row_indices: Vec<usize>,
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

/// Per-supernode diagnostic statistics.
///
/// Collected during [`AptpNumeric::factor`] for each supernode processed.
/// Provides visibility into per-front behavior for analysis and comparison
/// with reference solvers (e.g., SPRAL).
#[derive(Debug, Clone)]
pub struct PerSupernodeStats {
    /// Supernode index in postorder.
    pub snode_id: usize,
    /// Frontal matrix dimension (m).
    pub front_size: usize,
    /// Number of fully-summed rows/columns (k = supernode cols + delayed from children).
    pub num_fully_summed: usize,
    /// Number of columns successfully eliminated (ne <= k).
    pub num_eliminated: usize,
    /// Number of delayed columns (k - ne).
    pub num_delayed: usize,
    /// Number of 1x1 pivots accepted at this supernode.
    pub num_1x1: usize,
    /// Number of 2x2 pivot pairs accepted at this supernode.
    pub num_2x2: usize,
    /// Maximum absolute L entry at this supernode (stability metric).
    pub max_l_entry: f64,
    /// Wall-clock time for scatter + extend-add assembly.
    #[cfg(feature = "diagnostic")]
    pub assembly_time: std::time::Duration,
    /// Wall-clock time for dense APTP kernel.
    #[cfg(feature = "diagnostic")]
    pub kernel_time: std::time::Duration,
    /// Wall-clock time for front factor extraction.
    #[cfg(feature = "diagnostic")]
    pub extraction_time: std::time::Duration,
    // -- Sub-phase timing (finer-grained diagnostic breakdown) --
    /// Wall-clock time for zeroing the frontal matrix (memset-equivalent).
    #[cfg(feature = "diagnostic")]
    pub zero_time: std::time::Duration,
    /// Wall-clock time for building the global-to-local mapping.
    #[cfg(feature = "diagnostic")]
    pub g2l_time: std::time::Duration,
    /// Wall-clock time for scattering original CSC entries into the frontal matrix.
    #[cfg(feature = "diagnostic")]
    pub scatter_time: std::time::Duration,
    /// Wall-clock time for extend-add (merging child contributions).
    #[cfg(feature = "diagnostic")]
    pub extend_add_time: std::time::Duration,
    /// Wall-clock time for extracting L11, D11, L21 from the factored frontal matrix.
    #[cfg(feature = "diagnostic")]
    pub extract_factors_time: std::time::Duration,
    /// Wall-clock time for extracting the contribution block.
    #[cfg(feature = "diagnostic")]
    pub extract_contrib_time: std::time::Duration,
    /// Wall-clock time for the deferred contribution GEMM.
    #[cfg(feature = "diagnostic")]
    pub contrib_gemm_time: std::time::Duration,
    /// Wall-clock time for resetting the global-to-local mapping after factoring.
    #[cfg(feature = "diagnostic")]
    pub g2l_reset_time: std::time::Duration,
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
    /// Number of supernodes before amalgamation.
    pub supernodes_before_amalgamation: usize,
    /// Number of supernodes after amalgamation.
    pub supernodes_after_amalgamation: usize,
    /// Number of merge operations performed during amalgamation.
    pub merges_performed: usize,
    /// Sum of assembly times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_assembly_time: std::time::Duration,
    /// Sum of kernel times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_kernel_time: std::time::Duration,
    /// Sum of extraction times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_extraction_time: std::time::Duration,
    // -- Sub-phase totals --
    /// Sum of zeroing times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_zero_time: std::time::Duration,
    /// Sum of g2l build times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_g2l_time: std::time::Duration,
    /// Sum of scatter times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_scatter_time: std::time::Duration,
    /// Sum of extend-add times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_extend_add_time: std::time::Duration,
    /// Sum of factor extraction times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_extract_factors_time: std::time::Duration,
    /// Sum of contribution extraction times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_extract_contrib_time: std::time::Duration,
    /// Sum of deferred contribution GEMM times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_contrib_gemm_time: std::time::Duration,
    /// Sum of g2l reset times across all supernodes.
    #[cfg(feature = "diagnostic")]
    pub total_g2l_reset_time: std::time::Duration,
    /// Number of small-leaf subtrees identified by classification.
    pub small_leaf_subtrees: usize,
    /// Number of supernodes processed via the small-leaf fast path.
    pub small_leaf_nodes: usize,
}

/// Assembled frontal matrix exported for diagnostic comparison.
#[cfg(feature = "diagnostic")]
pub struct ExportedFrontal {
    /// Dense m x m frontal matrix (lower triangle populated).
    pub data: Mat<f64>,
    /// Total front size.
    pub front_size: usize,
    /// Number of fully-summed columns.
    pub num_fully_summed: usize,
    /// Global permuted column indices for each local row.
    pub row_indices: Vec<usize>,
}

/// Complete multifrontal numeric factorization result.
///
/// Contains per-supernode factors and aggregate statistics. Created by
/// [`AptpNumeric::factor`] and used by the triangular solve to
/// perform forward/backward substitution.
///
/// # Usage
///
/// ```
/// use faer::sparse::{SparseColMat, Triplet};
/// use faer::sparse::linalg::cholesky::SymmetricOrdering;
/// use rivrs_sparse::symmetric::{AptpSymbolic, AptpOptions, AptpNumeric};
///
/// // Analyze
/// # let triplets = vec![Triplet::new(0, 0, 1.0)];
/// # let matrix = SparseColMat::try_new_from_triplets(1, 1, &triplets).unwrap();
/// let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd).unwrap();
///
/// // Factor
/// let numeric = AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default(), None).unwrap();
/// println!("Stats: {:?}", numeric.stats());
/// ```
pub struct AptpNumeric {
    /// Per-supernode factors, indexed by supernode ID.
    front_factors: Vec<FrontFactors>,
    /// Aggregate factorization statistics.
    stats: FactorizationStats,
    /// Per-supernode diagnostic statistics (same order as `front_factors`).
    per_supernode_stats: Vec<PerSupernodeStats>,
    /// Matrix dimension.
    n: usize,
    /// Profiling session from factorization (Chrome Trace export).
    #[cfg(feature = "diagnostic")]
    profile_session: Option<FinishedSession>,
}

impl std::fmt::Debug for AptpNumeric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AptpNumeric")
            .field("n", &self.n)
            .field("stats", &self.stats)
            .field("front_factors_count", &self.front_factors.len())
            .finish()
    }
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

    /// Per-supernode diagnostic statistics.
    ///
    /// One entry per supernode in postorder, matching [`front_factors`](Self::front_factors).
    pub fn per_supernode_stats(&self) -> &[PerSupernodeStats] {
        &self.per_supernode_stats
    }

    /// Matrix dimension.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Profiling session from factorization, for Chrome Trace export.
    ///
    /// Returns `None` if factorization was not profiled.
    /// Only available with `diagnostic` feature.
    #[cfg(feature = "diagnostic")]
    pub fn profile_session(&self) -> Option<&FinishedSession> {
        self.profile_session.as_ref()
    }

    /// Factor a sparse symmetric matrix using the multifrontal method with APTP.
    ///
    /// Traverses the assembly tree in postorder, assembling and factoring dense
    /// frontal matrices at each supernode using the dense APTP kernel.
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
    /// an error. The factorization succeeds;
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
        scaling: Option<&[f64]>,
    ) -> Result<Self, SparseError> {
        let nemin = options.nemin;
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
        let n_supernodes_before = supernodes.len();
        let mut supernodes = super::amalgamation::amalgamate(supernodes, nemin);
        let n_supernodes = supernodes.len();
        let merges_performed = n_supernodes_before - n_supernodes;
        let children = build_children_map(&supernodes);

        // Classify small-leaf subtrees for the fast-path pre-pass
        let small_leaf_subtrees =
            classify_small_leaf_subtrees(&mut supernodes, &children, options.small_leaf_threshold);

        // Get fill-reducing permutation
        let (perm_fwd, perm_inv) = symbolic.perm_vecs();

        // Precompute assembly index maps (scatter + extend-add) for the
        // zero-delay case. Supernodes with delayed children fall back to
        // per-entry index arithmetic.
        let assembly_maps =
            build_assembly_maps(&supernodes, &children, matrix, &perm_fwd, &perm_inv, n);

        // Profiling session (diagnostic only)
        #[cfg(feature = "diagnostic")]
        let session = ProfileSession::new();
        #[cfg(feature = "diagnostic")]
        let _factor_loop_guard = session.enter_section("factor_loop");

        // Iterative level-set factorization: process all ready (leaf) nodes first,
        // then propagate upward. No recursion — avoids stack overflow on deep
        // elimination trees (e.g., c-71 with 76K supernodes on rayon's 2MB stacks).
        let all_node_results = factor_tree_levelset(
            &supernodes,
            &children,
            matrix,
            &perm_fwd,
            &perm_inv,
            options,
            scaling,
            n,
            &assembly_maps,
            &small_leaf_subtrees,
        )?;

        #[cfg(feature = "diagnostic")]
        drop(_factor_loop_guard);

        // Scatter returned results into supernode-indexed vectors
        let mut front_factors_vec: Vec<Option<FrontFactors>> =
            (0..n_supernodes).map(|_| None).collect();
        let mut per_sn_stats_vec: Vec<Option<PerSupernodeStats>> =
            (0..n_supernodes).map(|_| None).collect();
        for (idx, ff, stats) in all_node_results {
            front_factors_vec[idx] = Some(ff);
            per_sn_stats_vec[idx] = Some(stats);
        }

        let front_factors: Vec<FrontFactors> = front_factors_vec
            .into_iter()
            .enumerate()
            .map(|(i, opt)| opt.unwrap_or_else(|| panic!("supernode {} not factored", i)))
            .collect();
        let per_sn_stats: Vec<PerSupernodeStats> = per_sn_stats_vec
            .into_iter()
            .enumerate()
            .map(|(i, opt)| opt.unwrap_or_else(|| panic!("supernode {} missing stats", i)))
            .collect();

        // Compute aggregate stats from per-supernode stats
        let mut stats = FactorizationStats {
            total_1x1_pivots: 0,
            total_2x2_pivots: 0,
            total_delayed: 0,
            zero_pivots: 0,
            max_front_size: 0,
            supernodes_before_amalgamation: n_supernodes_before,
            supernodes_after_amalgamation: n_supernodes,
            merges_performed,
            #[cfg(feature = "diagnostic")]
            total_assembly_time: std::time::Duration::ZERO,
            #[cfg(feature = "diagnostic")]
            total_kernel_time: std::time::Duration::ZERO,
            #[cfg(feature = "diagnostic")]
            total_extraction_time: std::time::Duration::ZERO,
            #[cfg(feature = "diagnostic")]
            total_zero_time: std::time::Duration::ZERO,
            #[cfg(feature = "diagnostic")]
            total_g2l_time: std::time::Duration::ZERO,
            #[cfg(feature = "diagnostic")]
            total_scatter_time: std::time::Duration::ZERO,
            #[cfg(feature = "diagnostic")]
            total_extend_add_time: std::time::Duration::ZERO,
            #[cfg(feature = "diagnostic")]
            total_extract_factors_time: std::time::Duration::ZERO,
            #[cfg(feature = "diagnostic")]
            total_extract_contrib_time: std::time::Duration::ZERO,
            #[cfg(feature = "diagnostic")]
            total_contrib_gemm_time: std::time::Duration::ZERO,
            #[cfg(feature = "diagnostic")]
            total_g2l_reset_time: std::time::Duration::ZERO,
            small_leaf_subtrees: small_leaf_subtrees.len(),
            small_leaf_nodes: small_leaf_subtrees.iter().map(|s| s.nodes.len()).sum(),
        };
        for sn_stat in &per_sn_stats {
            stats.total_1x1_pivots += sn_stat.num_1x1;
            stats.total_2x2_pivots += sn_stat.num_2x2;
            stats.total_delayed += sn_stat.num_delayed;
            if sn_stat.front_size > stats.max_front_size {
                stats.max_front_size = sn_stat.front_size;
            }
            #[cfg(feature = "diagnostic")]
            {
                stats.total_assembly_time += sn_stat.assembly_time;
                stats.total_kernel_time += sn_stat.kernel_time;
                stats.total_extraction_time += sn_stat.extraction_time;
                stats.total_zero_time += sn_stat.zero_time;
                stats.total_g2l_time += sn_stat.g2l_time;
                stats.total_scatter_time += sn_stat.scatter_time;
                stats.total_extend_add_time += sn_stat.extend_add_time;
                stats.total_extract_factors_time += sn_stat.extract_factors_time;
                stats.total_extract_contrib_time += sn_stat.extract_contrib_time;
                stats.total_contrib_gemm_time += sn_stat.contrib_gemm_time;
                stats.total_g2l_reset_time += sn_stat.g2l_reset_time;
            }
        }

        // Count zero pivots from front factors
        for ff in &front_factors {
            for (_, entry) in ff.d11().iter_pivots() {
                if let PivotEntry::OneByOne(val) = entry {
                    if val == 0.0 {
                        stats.zero_pivots += 1;
                    }
                }
            }
        }
        // Also count unresolved delays at root supernodes
        for (s, sn) in supernodes.iter().enumerate() {
            if sn.parent.is_none() {
                let sn_stat = &per_sn_stats[s];
                let m = sn_stat.front_size;
                let ne = sn_stat.num_eliminated;
                if ne < m {
                    let k = sn_stat.num_fully_summed;
                    let n_unresolved = k - ne;
                    stats.zero_pivots += n_unresolved;
                }
            }
        }

        // Finish profiling session
        #[cfg(feature = "diagnostic")]
        let finished_session = session.finish();

        Ok(AptpNumeric {
            front_factors,
            stats,
            per_supernode_stats: per_sn_stats,
            n,
            #[cfg(feature = "diagnostic")]
            profile_session: Some(finished_session),
        })
    }

    /// Export the assembled (pre-factorization) frontal matrix for a specific supernode.
    ///
    /// Runs the multifrontal assembly up to the target supernode, factoring all
    /// preceding supernodes to generate their contributions, then returns the
    /// assembled frontal matrix for the target **before** APTP factoring.
    ///
    /// # Selecting the Target Supernode
    ///
    /// If `target_snode` is `None`, the supernode with the largest front size is used.
    /// Otherwise, the specified supernode index is used.
    ///
    /// # Use Case
    ///
    /// This is a diagnostic tool for comparing our APTP dense kernel with SPRAL's
    /// `ldlt_app_factor` on the same assembled input. The exported matrix can be
    /// written to a file and fed to a SPRAL driver.
    #[cfg(feature = "diagnostic")]
    pub fn export_assembled_frontal(
        symbolic: &AptpSymbolic,
        matrix: &SparseColMat<usize, f64>,
        options: &AptpOptions,
        scaling: Option<&[f64]>,
        target_snode: Option<usize>,
    ) -> Result<ExportedFrontal, SparseError> {
        let n = symbolic.nrows();

        if matrix.nrows() != n || matrix.ncols() != n {
            return Err(SparseError::DimensionMismatch {
                expected: (n, n),
                got: (matrix.nrows(), matrix.ncols()),
                context: "Matrix dimensions must match symbolic analysis".to_string(),
            });
        }

        let supernodes = build_supernode_info(symbolic);
        let n_supernodes = supernodes.len();
        let children = build_children_map(&supernodes);
        let (perm_fwd, perm_inv) = symbolic.perm_vecs();

        // Determine target supernode
        let target = match target_snode {
            Some(t) => {
                if t >= n_supernodes {
                    return Err(SparseError::AnalysisFailure {
                        reason: format!(
                            "Target supernode {} out of range (n_supernodes={})",
                            t, n_supernodes
                        ),
                    });
                }
                t
            }
            None => {
                // Find supernode with largest front
                let mut best = 0;
                let mut best_m = 0;
                for (s, sn) in supernodes.iter().enumerate() {
                    let sn_ncols = sn.col_end - sn.col_begin;
                    // Approximate: delayed cols not counted here, just pattern + cols
                    let m_approx = sn_ncols + sn.pattern.len();
                    if m_approx > best_m {
                        best_m = m_approx;
                        best = s;
                    }
                }
                best
            }
        };

        let mut global_to_local = vec![NOT_IN_FRONT; n];
        let mut contributions: Vec<Option<ContributionBlock>> =
            (0..n_supernodes).map(|_| None).collect();

        // Process supernodes in postorder up to and including target
        for s in 0..=target {
            let sn = &supernodes[s];

            // 1. Collect delayed columns from children
            let mut delayed_cols: Vec<usize> = Vec::new();
            for &c in &children[s] {
                if let Some(ref cb) = contributions[c] {
                    delayed_cols.extend_from_slice(&cb.row_indices[..cb.num_delayed]);
                }
            }

            // 2. Compute frontal matrix structure
            let sn_ncols: usize = sn.owned_ranges.iter().map(|r| r.len()).sum();
            let k = sn_ncols + delayed_cols.len();
            let mut frontal_rows: Vec<usize> = Vec::with_capacity(k + sn.pattern.len());
            for range in &sn.owned_ranges {
                frontal_rows.extend(range.clone());
            }
            frontal_rows.extend_from_slice(&delayed_cols);
            frontal_rows.extend_from_slice(&sn.pattern);
            let m = frontal_rows.len();

            // 3. Build global-to-local mapping
            for (i, &global) in frontal_rows.iter().enumerate() {
                global_to_local[global] = i;
            }

            // 4. Assemble frontal matrix (owned data — diagnostic path)
            let mut frontal_data = Mat::zeros(m, m);
            {
                let mut frontal = FrontalMatrix {
                    data: frontal_data.as_mut(),
                    row_indices: &frontal_rows,
                    num_fully_summed: k,
                };

                scatter_original_entries_multi(
                    &mut frontal,
                    matrix,
                    &perm_fwd,
                    &perm_inv,
                    &global_to_local,
                    &sn.owned_ranges,
                    scaling,
                );

                for &c in &children[s] {
                    if let Some(cb) = contributions[c].take() {
                        let _ = extend_add(&mut frontal, cb, &global_to_local);
                    }
                }
            }

            if s == target {
                // Return the assembled frontal matrix BEFORE factoring
                // Clean up global_to_local
                for &global in &frontal_rows {
                    global_to_local[global] = NOT_IN_FRONT;
                }
                return Ok(ExportedFrontal {
                    data: frontal_data,
                    front_size: m,
                    num_fully_summed: k,
                    row_indices: frontal_rows,
                });
            }

            // Factor this supernode (needed for contribution propagation)
            let result = aptp_factor_in_place(frontal_data.as_mut(), k, options)?;
            let ne = result.num_eliminated;

            if sn.parent.is_some() && ne < m {
                let nfs = m - k;
                let mut cb = Mat::zeros(m, m);
                if nfs > 0 {
                    let mut ld_ws = Mat::zeros(m, ne.max(1));
                    compute_contribution_gemm(
                        &frontal_data,
                        k,
                        ne,
                        m,
                        &result.d,
                        &mut cb,
                        &mut ld_ws,
                        Par::Seq,
                    );
                }
                contributions[s] = Some(extract_contribution(
                    &frontal_data,
                    m,
                    k,
                    &frontal_rows,
                    &result,
                    cb,
                ));
            }

            // Cleanup
            for &global in &frontal_rows {
                global_to_local[global] = NOT_IN_FRONT;
            }
        }

        // Should not reach here — we return inside the loop at s == target
        unreachable!("Target supernode {} not reached in postorder loop", target)
    }
}

// ---------------------------------------------------------------------------
// Tree-level factorization functions
// ---------------------------------------------------------------------------

/// Result of factoring a single supernode.
struct SupernodeResult {
    ff: FrontFactors,
    contribution: Option<ContributionBlock>,
    stats: PerSupernodeStats,
}

/// Factor a single supernode given its children's contribution blocks.
///
/// Assembles the frontal matrix from sparse entries and child contributions,
/// factors it with APTP, and extracts the front factors and contribution block.
///
/// Uses the pre-allocated `workspace` for the frontal matrix buffer, row indices,
/// delayed columns, and global-to-local mapping to avoid per-supernode allocations.
#[allow(clippy::too_many_arguments)]
fn factor_single_supernode(
    s: usize,
    sn: &SupernodeInfo,
    child_contributions: Vec<Option<ContributionBlock>>,
    matrix: &SparseColMat<usize, f64>,
    perm_fwd: &[usize],
    perm_inv: &[usize],
    options: &AptpOptions,
    scaling: Option<&[f64]>,
    workspace: &mut FactorizationWorkspace,
    assembly_maps: &AssemblyMaps,
) -> Result<SupernodeResult, SparseError> {
    // 1. Collect delayed columns from children into workspace buffer
    workspace.delayed_cols_buf.clear();
    for cb in child_contributions.iter().flatten() {
        workspace
            .delayed_cols_buf
            .extend_from_slice(&cb.row_indices[..cb.num_delayed]);
    }

    // 2. Compute frontal matrix structure
    // Use owned_ranges (not col_begin..col_end) to enumerate fully-summed columns.
    // After amalgamation of non-adjacent supernodes, col_begin..col_end may span
    // columns belonging to other supernodes. Only owned columns should be in the front.
    let sn_ncols: usize = sn.owned_ranges.iter().map(|r| r.len()).sum();
    let k = sn_ncols + workspace.delayed_cols_buf.len(); // num_fully_summed

    workspace.frontal_row_indices.clear();
    workspace.frontal_row_indices.reserve(k + sn.pattern.len());
    for range in &sn.owned_ranges {
        workspace.frontal_row_indices.extend(range.clone());
    }
    workspace
        .frontal_row_indices
        .extend_from_slice(&workspace.delayed_cols_buf);
    workspace.frontal_row_indices.extend_from_slice(&sn.pattern);
    let m = workspace.frontal_row_indices.len();

    // 3. Ensure workspace capacity and zero frontal data for this supernode
    // If delayed columns make this front larger than estimated, grow the workspace
    if m > workspace.frontal_data.nrows() {
        workspace.frontal_data = Mat::zeros(m, m);
    }

    #[cfg(feature = "diagnostic")]
    let zero_start = std::time::Instant::now();

    workspace.zero_frontal(m);

    #[cfg(feature = "diagnostic")]
    let zero_time = zero_start.elapsed();

    // 4. Build global-to-local mapping
    #[cfg(feature = "diagnostic")]
    let g2l_start = std::time::Instant::now();

    for (i, &global) in workspace.frontal_row_indices.iter().enumerate() {
        workspace.global_to_local[global] = i;
    }

    #[cfg(feature = "diagnostic")]
    let g2l_time = g2l_start.elapsed();

    // 5. Assemble frontal matrix using workspace buffer
    #[cfg(feature = "diagnostic")]
    let assembly_start = std::time::Instant::now();

    let ndelay_in = workspace.delayed_cols_buf.len();

    let mut frontal = FrontalMatrix {
        data: workspace.frontal_data.as_mut().submatrix_mut(0, 0, m, m),
        row_indices: &workspace.frontal_row_indices,
        num_fully_summed: k,
    };

    // 5a. Scatter original CSC entries
    #[cfg(feature = "diagnostic")]
    let scatter_start = std::time::Instant::now();

    // Use precomputed amap when no delays shift row positions
    if ndelay_in == 0 {
        // Fast path: scatter using precomputed index tuples (4 u32 per entry)
        let amap_start = assembly_maps.amap_offsets[s];
        let amap_end = assembly_maps.amap_offsets[s + 1];
        let values = matrix.val();
        let amap = &assembly_maps.amap_entries[amap_start * 4..amap_end * 4];

        if let Some(sc) = scaling {
            for entry in amap.chunks_exact(4) {
                let src_idx = entry[0] as usize;
                let dest_linear = entry[1] as usize;
                let scale_row = entry[2] as usize;
                let scale_col = entry[3] as usize;
                let val = values[src_idx] * sc[scale_row] * sc[scale_col];
                let dest_col = dest_linear / m;
                let dest_row = dest_linear % m;
                frontal.data[(dest_row, dest_col)] += val;
            }
        } else {
            for entry in amap.chunks_exact(4) {
                let src_idx = entry[0] as usize;
                let dest_linear = entry[1] as usize;
                let dest_col = dest_linear / m;
                let dest_row = dest_linear % m;
                frontal.data[(dest_row, dest_col)] += values[src_idx];
            }
        }
    } else {
        // Fallback: delayed columns shift positions, use per-entry assembly
        scatter_original_entries_multi(
            &mut frontal,
            matrix,
            perm_fwd,
            perm_inv,
            &workspace.global_to_local,
            &sn.owned_ranges,
            scaling,
        );
    }

    #[cfg(feature = "diagnostic")]
    let scatter_time = scatter_start.elapsed();

    // 5b. Extend-add: merge child contributions into parent frontal matrix
    #[cfg(feature = "diagnostic")]
    let ea_start_time = std::time::Instant::now();

    let ea_children_begin = assembly_maps.ea_snode_child_begin[s];
    let mut ea_child_idx = ea_children_begin;

    for opt_cb in child_contributions {
        if let Some(cb) = opt_cb {
            // Check if this child had delays — if so, fall back to g2l
            let recycled = if cb.num_delayed > 0 || ndelay_in > 0 {
                extend_add(&mut frontal, cb, &workspace.global_to_local)
            } else {
                // Fast path: use precomputed row mapping
                let ea_start = assembly_maps.ea_offsets[ea_child_idx];
                let ea_end = assembly_maps.ea_offsets[ea_child_idx + 1];
                let ea_row_map = &assembly_maps.ea_map[ea_start..ea_end];
                extend_add_mapped(&mut frontal, cb, ea_row_map)
            };
            // Recycle the returned buffer into workspace if it's larger than
            // what we currently have (or if ours was moved out).
            if recycled.nrows() >= workspace.contrib_buffer.nrows() {
                workspace.contrib_buffer = recycled;
            }
        }
        // Always increment: ea_offsets indexes all children, including those
        // with no contribution (fully eliminated, ne == m).
        ea_child_idx += 1;
    }

    #[cfg(feature = "diagnostic")]
    let extend_add_time = ea_start_time.elapsed();

    #[cfg(feature = "diagnostic")]
    let assembly_time = assembly_start.elapsed();

    // 6. Factor the frontal matrix
    #[cfg(feature = "diagnostic")]
    let kernel_start = std::time::Instant::now();

    let effective_par = if m < INTRA_NODE_THRESHOLD {
        Par::Seq
    } else {
        options.par
    };
    let per_sn_options = AptpOptions {
        par: effective_par,
        ..options.clone()
    };
    let result = aptp_factor_in_place(
        workspace.frontal_data.as_mut().submatrix_mut(0, 0, m, m),
        k,
        &per_sn_options,
    )?;
    let ne = result.num_eliminated;

    #[cfg(feature = "diagnostic")]
    let kernel_time = kernel_start.elapsed();

    // 6b. Deferred contribution GEMM: compute NFS×NFS Schur complement directly
    //     into workspace.contrib_buffer. This replaces the per-block NFS×NFS
    //     updates that were skipped during the blocking loop.
    #[cfg(feature = "diagnostic")]
    let contrib_gemm_start = std::time::Instant::now();

    let nfs = m - k;
    if sn.parent.is_some() && ne < m && nfs > 0 {
        // Ensure contrib_buffer is large enough (rare: delayed cascades may exceed estimate)
        if workspace.contrib_buffer.nrows() < m || workspace.contrib_buffer.ncols() < m {
            workspace.contrib_buffer = Mat::zeros(m, m);
        }
        compute_contribution_gemm(
            &workspace.frontal_data,
            k,
            ne,
            m,
            &result.d,
            &mut workspace.contrib_buffer,
            &mut workspace.ld_workspace,
            effective_par,
        );
    }

    #[cfg(feature = "diagnostic")]
    let contrib_gemm_time = contrib_gemm_start.elapsed();

    // 7. Extract front factors (reads from workspace frontal data)
    #[cfg(feature = "diagnostic")]
    let extraction_start = std::time::Instant::now();
    #[cfg(feature = "diagnostic")]
    let extract_factors_start = std::time::Instant::now();

    // Use owned Mat<f64> directly for col_as_slice-based bulk extraction
    let ff = extract_front_factors(
        &workspace.frontal_data,
        m,
        k,
        &workspace.frontal_row_indices,
        &result,
    );

    #[cfg(feature = "diagnostic")]
    let extract_factors_time = extract_factors_start.elapsed();

    // 8. Prepare contribution for parent (if not root and not fully eliminated)
    //    NFS×NFS data is in workspace.contrib_buffer (from deferred GEMM).
    //    Delayed-column data (if any) is still in workspace.frontal_data.
    #[cfg(feature = "diagnostic")]
    let extract_contrib_start = std::time::Instant::now();

    let contribution = if sn.parent.is_some() && ne < m {
        Some(extract_contribution(
            &workspace.frontal_data,
            m,
            k,
            &workspace.frontal_row_indices,
            &result,
            std::mem::replace(&mut workspace.contrib_buffer, Mat::new()),
        ))
    } else {
        None
    };

    #[cfg(feature = "diagnostic")]
    let extract_contrib_time = extract_contrib_start.elapsed();

    #[cfg(feature = "diagnostic")]
    let extraction_time = extraction_start.elapsed();

    let stats = PerSupernodeStats {
        snode_id: s,
        front_size: m,
        num_fully_summed: k,
        num_eliminated: ne,
        num_delayed: k - ne,
        num_1x1: result.stats.num_1x1,
        num_2x2: result.stats.num_2x2,
        max_l_entry: result.stats.max_l_entry,
        #[cfg(feature = "diagnostic")]
        assembly_time,
        #[cfg(feature = "diagnostic")]
        kernel_time,
        #[cfg(feature = "diagnostic")]
        extraction_time,
        #[cfg(feature = "diagnostic")]
        zero_time,
        #[cfg(feature = "diagnostic")]
        g2l_time,
        #[cfg(feature = "diagnostic")]
        scatter_time,
        #[cfg(feature = "diagnostic")]
        extend_add_time,
        #[cfg(feature = "diagnostic")]
        extract_factors_time,
        #[cfg(feature = "diagnostic")]
        extract_contrib_time,
        #[cfg(feature = "diagnostic")]
        contrib_gemm_time,
        #[cfg(feature = "diagnostic")]
        g2l_reset_time: std::time::Duration::ZERO, // filled below after reset
    };

    // Reset global_to_local entries for reuse (O(m), not O(n))
    #[cfg(feature = "diagnostic")]
    let g2l_reset_start = std::time::Instant::now();

    for &global in &workspace.frontal_row_indices[..m] {
        workspace.global_to_local[global] = NOT_IN_FRONT;
    }

    #[cfg(feature = "diagnostic")]
    let g2l_reset_time = g2l_reset_start.elapsed();
    #[cfg(feature = "diagnostic")]
    let stats = {
        let mut s = stats;
        s.g2l_reset_time = g2l_reset_time;
        s
    };

    Ok(SupernodeResult {
        ff,
        contribution,
        stats,
    })
}

/// Factor a small-leaf subtree using rectangular L storage.
///
/// Replaces the general `factor_single_supernode` path for classified small-leaf
/// subtrees with a streamlined per-node loop that:
/// 1. Uses rectangular m×k L storage instead of square m×m frontal matrix
/// 2. Factors via `tpp_factor_rect` (rectangular TPP kernel)
/// 3. Computes contribution GEMM from rectangular storage
/// 4. Reuses a single small workspace across all nodes in the subtree
///
/// # Arguments
///
/// - `subtree`: Small-leaf subtree descriptor (nodes in postorder).
/// - `supernodes`: All supernode descriptors.
/// - `children`: Children map.
/// - `matrix`: Original sparse matrix.
/// - `perm_fwd`/`perm_inv`: Forward/inverse permutations.
/// - `options`: APTP options.
/// - `scaling`: Optional MC64 scaling factors.
/// - `n`: Matrix dimension.
/// - `contributions`: Global contribution storage (read children, write results).
/// - `assembly_maps`: Precomputed scatter maps.
///
/// # Returns
///
/// Vec of (node_id, FrontFactors, PerSupernodeStats) for each node in the subtree.
#[allow(clippy::too_many_arguments)]
fn factor_small_leaf_subtree(
    subtree: &SmallLeafSubtree,
    supernodes: &[SupernodeInfo],
    children: &[Vec<usize>],
    matrix: &SparseColMat<usize, f64>,
    perm_fwd: &[usize],
    perm_inv: &[usize],
    options: &AptpOptions,
    scaling: Option<&[f64]>,
    n: usize,
    contributions: &mut [Option<ContributionBlock>],
    assembly_maps: &AssemblyMaps,
) -> Result<Vec<(usize, FrontFactors, PerSupernodeStats)>, SparseError> {
    let max_front = subtree.max_front_size;
    let mut results = Vec::with_capacity(subtree.nodes.len());

    // Shared workspace for the subtree:
    // - l_storage: rectangular m×k buffer, reused per node (max_front × max_front covers all)
    // - global_to_local: length n, reused per node
    // - contrib_buffer: for contribution GEMM output, recycled
    // - ld_workspace: for L·D product in contribution GEMM
    // - frontal_row_indices: reusable buffer for building front structure
    // - delayed_cols_buf: reusable buffer for collecting child delays
    let mut l_storage = Mat::<f64>::zeros(max_front, max_front);
    let mut global_to_local = vec![NOT_IN_FRONT; n];
    let mut contrib_buffer = Mat::<f64>::new();
    let mut ld_workspace = Mat::<f64>::new();
    let mut frontal_row_indices = Vec::<usize>::with_capacity(max_front);
    let mut delayed_cols_buf = Vec::<usize>::new();

    let symbolic_matrix = matrix.symbolic();
    let col_ptrs = symbolic_matrix.col_ptr();
    let row_indices_csc = symbolic_matrix.row_idx();
    let values = matrix.val();

    for &node_id in &subtree.nodes {
        let sn = &supernodes[node_id];

        // 1. Collect delayed columns from children
        delayed_cols_buf.clear();
        for &c in &children[node_id] {
            if let Some(ref cb) = contributions[c] {
                delayed_cols_buf.extend_from_slice(&cb.row_indices[..cb.num_delayed]);
            }
        }

        // 2. Compute frontal matrix structure
        let sn_ncols: usize = sn.owned_ranges.iter().map(|r| r.len()).sum();
        let k = sn_ncols + delayed_cols_buf.len(); // num_fully_summed

        frontal_row_indices.clear();
        frontal_row_indices.reserve(k + sn.pattern.len());
        for range in &sn.owned_ranges {
            frontal_row_indices.extend(range.clone());
        }
        frontal_row_indices.extend_from_slice(&delayed_cols_buf);
        frontal_row_indices.extend_from_slice(&sn.pattern);
        let m = frontal_row_indices.len();

        // 3. Ensure l_storage capacity and zero rectangular m×k region
        #[cfg(feature = "diagnostic")]
        let zero_start = std::time::Instant::now();

        if m > l_storage.nrows() || k > l_storage.ncols() {
            l_storage = Mat::zeros(m, k.max(l_storage.ncols()));
        }
        // Zero only the m×k rectangle (smaller than m×m for NFS nodes)
        for j in 0..k {
            l_storage.col_as_slice_mut(j)[0..m].fill(0.0);
        }

        #[cfg(feature = "diagnostic")]
        let zero_time = zero_start.elapsed();

        // 4. Build global-to-local mapping
        #[cfg(feature = "diagnostic")]
        let g2l_start = std::time::Instant::now();

        for (i, &global) in frontal_row_indices.iter().enumerate() {
            global_to_local[global] = i;
        }

        #[cfg(feature = "diagnostic")]
        let g2l_time = g2l_start.elapsed();

        // 5. Split assembly: scatter into l_storage (FS columns) and contrib_buffer (NFS×NFS)
        #[cfg(feature = "diagnostic")]
        let assembly_start = std::time::Instant::now();
        #[cfg(feature = "diagnostic")]
        let scatter_start = std::time::Instant::now();

        let ndelay_in = delayed_cols_buf.len();
        let total_owned = sn_ncols;
        let nfs = m - k;

        // Pre-allocate and zero contrib_buffer if there will be NFS rows.
        // This allows both original matrix scatter and child extend-add to
        // deposit NFS×NFS entries directly, avoiding a re-scatter pass.
        let has_nfs = nfs > 0 && sn.parent.is_some();
        if has_nfs {
            if contrib_buffer.nrows() < nfs || contrib_buffer.ncols() < nfs {
                contrib_buffer = Mat::zeros(
                    nfs.max(contrib_buffer.nrows()),
                    nfs.max(contrib_buffer.ncols()),
                );
            }
            for j in 0..nfs {
                contrib_buffer.col_as_slice_mut(j)[0..nfs].fill(0.0);
            }
        }

        // 5a. Scatter original CSC entries — FS columns → l_storage, NFS×NFS → contrib_buffer
        if ndelay_in == 0 {
            let amap_start = assembly_maps.amap_offsets[node_id];
            let amap_end = assembly_maps.amap_offsets[node_id + 1];
            let amap = &assembly_maps.amap_entries[amap_start * 4..amap_end * 4];

            if let Some(sc) = scaling {
                for entry in amap.chunks_exact(4) {
                    let src_idx = entry[0] as usize;
                    let dest_linear = entry[1] as usize;
                    let scale_row = entry[2] as usize;
                    let scale_col = entry[3] as usize;
                    let val = values[src_idx] * sc[scale_row] * sc[scale_col];
                    let amap_col = dest_linear / m;
                    let amap_row = dest_linear % m;
                    if amap_col < k {
                        l_storage[(amap_row, amap_col)] += val;
                    } else if has_nfs && amap_row >= k {
                        // NFS×NFS entry → contrib_buffer
                        contrib_buffer[(amap_row - k, amap_col - k)] += val;
                    }
                    // FS×NFS cross-terms (amap_col >= k, amap_row < k) are handled
                    // by the symmetric counterpart where amap_col < k.
                }
            } else {
                for entry in amap.chunks_exact(4) {
                    let src_idx = entry[0] as usize;
                    let dest_linear = entry[1] as usize;
                    let amap_col = dest_linear / m;
                    let amap_row = dest_linear % m;
                    if amap_col < k {
                        l_storage[(amap_row, amap_col)] += values[src_idx];
                    } else if has_nfs && amap_row >= k {
                        contrib_buffer[(amap_row - k, amap_col - k)] += values[src_idx];
                    }
                }
            }
        } else {
            // Fallback: delayed columns shift positions, use per-entry g2l scatter
            for range in &sn.owned_ranges {
                for pj in range.clone() {
                    let orig_col = perm_fwd[pj];
                    let local_col = global_to_local[pj];

                    let start = col_ptrs[orig_col];
                    let end = col_ptrs[orig_col + 1];
                    for idx in start..end {
                        let orig_row = row_indices_csc[idx];
                        if orig_row < orig_col {
                            let perm_row = perm_inv[orig_row];
                            let local_peer = global_to_local[perm_row];
                            if local_peer != NOT_IN_FRONT && local_peer < total_owned {
                                continue;
                            }
                        }
                        let global_row = perm_inv[orig_row];
                        let local_row = global_to_local[global_row];
                        if local_row == NOT_IN_FRONT {
                            continue;
                        }
                        if local_row >= total_owned && local_row < k {
                            continue;
                        }
                        let mut val = values[idx];
                        if let Some(s) = scaling {
                            val *= s[perm_inv[orig_row]] * s[perm_inv[orig_col]];
                        }
                        let (r, c) = if local_row >= local_col {
                            (local_row, local_col)
                        } else {
                            (local_col, local_row)
                        };
                        if c < k {
                            l_storage[(r, c)] += val;
                        } else if has_nfs && r >= k {
                            contrib_buffer[(r - k, c - k)] += val;
                        }
                    }
                }
            }
        }

        #[cfg(feature = "diagnostic")]
        let scatter_time = scatter_start.elapsed();

        // 5b. Extend-add child contributions
        //     FS columns → l_storage, NFS×NFS → contrib_buffer
        #[cfg(feature = "diagnostic")]
        let extend_add_start = std::time::Instant::now();

        for &c in &children[node_id] {
            if let Some(cb) = contributions[c].take() {
                let cb_size = cb.row_indices.len();
                for i in 0..cb_size {
                    let gi = cb.row_indices[i];
                    let li = global_to_local[gi];
                    debug_assert!(li != NOT_IN_FRONT, "child row {} not in parent", gi);
                    for j in 0..=i {
                        let gj = cb.row_indices[j];
                        let lj = global_to_local[gj];
                        debug_assert!(lj != NOT_IN_FRONT, "child col {} not in parent", gj);
                        let val = cb.data[(i, j)];
                        if val != 0.0 {
                            let (li, lj) = if li >= lj { (li, lj) } else { (lj, li) };
                            if lj < k {
                                // FS column → l_storage
                                l_storage[(li, lj)] += val;
                            } else if has_nfs {
                                // NFS×NFS → contrib_buffer
                                contrib_buffer[(li - k, lj - k)] += val;
                            }
                        }
                    }
                }
                // Only recycle buffer if we're NOT using contrib_buffer for NFS assembly
                if !has_nfs && cb.data.nrows() >= contrib_buffer.nrows() {
                    contrib_buffer = cb.data;
                }
            }
        }

        #[cfg(feature = "diagnostic")]
        let extend_add_time = extend_add_start.elapsed();
        #[cfg(feature = "diagnostic")]
        let assembly_time = assembly_start.elapsed();

        // 6. Factor using rectangular TPP kernel
        #[cfg(feature = "diagnostic")]
        let kernel_start = std::time::Instant::now();

        let result = tpp_factor_rect(l_storage.as_mut().submatrix_mut(0, 0, m, k), k, options)?;
        let ne = result.num_eliminated;

        #[cfg(feature = "diagnostic")]
        let kernel_time = kernel_start.elapsed();

        // 7. Build contribution block
        #[cfg(feature = "diagnostic")]
        let mut contrib_gemm_time = std::time::Duration::ZERO;
        #[cfg(feature = "diagnostic")]
        let mut extract_contrib_time = std::time::Duration::ZERO;
        let contribution = if sn.parent.is_some() && ne < m {
            // Ensure contrib_buffer is large enough for extract_contribution_rect
            // which may need to assemble delayed + NFS into a single buffer
            if contrib_buffer.nrows() < m || contrib_buffer.ncols() < m {
                // Grow if delayed columns expanded the need
                let new_size = m.max(contrib_buffer.nrows());
                let mut new_buf = Mat::zeros(new_size, new_size);
                // Preserve NFS×NFS data already assembled
                if nfs > 0 {
                    for j in 0..nfs {
                        let col_len = nfs - j;
                        let src = &contrib_buffer.col_as_slice(j)[j..j + col_len];
                        new_buf.col_as_slice_mut(j)[j..j + col_len].copy_from_slice(src);
                    }
                }
                contrib_buffer = new_buf;
            }

            // Apply contribution GEMM: updates NFS×NFS with -L21_NFS * D * L21_NFS^T
            #[cfg(feature = "diagnostic")]
            let gemm_start = std::time::Instant::now();

            if nfs > 0 && ne > 0 {
                compute_contribution_gemm_rect(
                    &l_storage,
                    k,
                    ne,
                    m,
                    &result.d,
                    &mut contrib_buffer,
                    &mut ld_workspace,
                    Par::Seq,
                );
            }

            #[cfg(feature = "diagnostic")]
            {
                contrib_gemm_time = gemm_start.elapsed();
            }

            #[cfg(feature = "diagnostic")]
            let extract_contrib_start = std::time::Instant::now();

            let cb = extract_contribution_rect(
                &l_storage,
                m,
                k,
                &frontal_row_indices,
                &result,
                std::mem::replace(&mut contrib_buffer, Mat::new()),
            );

            #[cfg(feature = "diagnostic")]
            {
                extract_contrib_time = extract_contrib_start.elapsed();
            }

            Some(cb)
        } else {
            None
        };

        // 8. Extract FrontFactors from rectangular L storage
        #[cfg(feature = "diagnostic")]
        let extract_factors_start = std::time::Instant::now();

        let ff = extract_front_factors_rect(&l_storage, m, k, &frontal_row_indices, &result);

        #[cfg(feature = "diagnostic")]
        let extract_factors_time = extract_factors_start.elapsed();

        let stats = PerSupernodeStats {
            snode_id: node_id,
            front_size: m,
            num_fully_summed: k,
            num_eliminated: ne,
            num_delayed: k - ne,
            num_1x1: result.stats.num_1x1,
            num_2x2: result.stats.num_2x2,
            max_l_entry: result.stats.max_l_entry,
            #[cfg(feature = "diagnostic")]
            assembly_time,
            #[cfg(feature = "diagnostic")]
            kernel_time,
            #[cfg(feature = "diagnostic")]
            extraction_time: extract_factors_time + extract_contrib_time,
            #[cfg(feature = "diagnostic")]
            zero_time,
            #[cfg(feature = "diagnostic")]
            g2l_time,
            #[cfg(feature = "diagnostic")]
            scatter_time,
            #[cfg(feature = "diagnostic")]
            extend_add_time,
            #[cfg(feature = "diagnostic")]
            extract_factors_time,
            #[cfg(feature = "diagnostic")]
            extract_contrib_time,
            #[cfg(feature = "diagnostic")]
            contrib_gemm_time,
            #[cfg(feature = "diagnostic")]
            g2l_reset_time: std::time::Duration::ZERO, // filled below
        };

        // 9. Reset g2l
        #[cfg(feature = "diagnostic")]
        let g2l_reset_start = std::time::Instant::now();

        for &global in &frontal_row_indices[..m] {
            global_to_local[global] = NOT_IN_FRONT;
        }

        #[cfg(feature = "diagnostic")]
        let stats = {
            let mut s = stats;
            s.g2l_reset_time = g2l_reset_start.elapsed();
            s
        };

        contributions[node_id] = contribution;
        results.push((node_id, ff, stats));
    }

    Ok(results)
}

/// Iterative level-set factorization of the assembly tree.
///
/// Processes supernodes bottom-up: first all leaves (in parallel if enabled),
/// then all nodes whose children are complete, and so on up to the roots.
/// This avoids deep recursion that overflows rayon's 2MB worker stacks on
/// matrices with deep elimination trees (e.g., c-71 with 76K supernodes).
///
/// The parallelism is "level-set" style: at each wave, all ready nodes are
/// independent (their children are already factored) and can be processed
/// concurrently. This provides the same tree-level parallelism as the
/// recursive approach but without stack depth concerns.
///
/// # References
///
/// - Liu (1992), Section 5: level-set scheduling for multifrontal methods
/// - Duff & Reid (1983): assembly tree parallelism
#[allow(clippy::too_many_arguments)]
fn factor_tree_levelset(
    supernodes: &[SupernodeInfo],
    children: &[Vec<usize>],
    matrix: &SparseColMat<usize, f64>,
    perm_fwd: &[usize],
    perm_inv: &[usize],
    options: &AptpOptions,
    scaling: Option<&[f64]>,
    n: usize,
    assembly_maps: &AssemblyMaps,
    small_leaf_subtrees: &[SmallLeafSubtree],
) -> Result<Vec<(usize, FrontFactors, PerSupernodeStats)>, SparseError> {
    let n_supernodes = supernodes.len();

    // Track contributions from completed children
    let mut contributions: Vec<Option<ContributionBlock>> =
        (0..n_supernodes).map(|_| None).collect();

    // Track how many children remain unprocessed for each supernode
    let mut remaining_children: Vec<usize> = children.iter().map(|c| c.len()).collect();

    // Collect all results
    let mut all_results: Vec<(usize, FrontFactors, PerSupernodeStats)> =
        Vec::with_capacity(n_supernodes);

    // --- Small-leaf subtree pre-pass (rectangular workspace fast path) ---
    // Factor all classified small-leaf subtrees before the main level-set loop.
    // Each subtree uses a rectangular m×k L storage kernel instead of the general
    // m×m frontal matrix, eliminating the large square allocation, frontal zeroing,
    // and extract_front_factors copy for every node.
    for subtree in small_leaf_subtrees {
        let subtree_results = factor_small_leaf_subtree(
            subtree,
            supernodes,
            children,
            matrix,
            perm_fwd,
            perm_inv,
            options,
            scaling,
            n,
            &mut contributions,
            assembly_maps,
        )?;
        all_results.extend(subtree_results);

        // Decrement remaining_children for the subtree root's parent
        if let Some(parent) = subtree.parent_of_root {
            remaining_children[parent] -= 1;
        }
    }

    let is_parallel = !matches!(options.par, Par::Seq);
    let max_front = estimate_max_front_size(supernodes);

    // Sequential mode: pre-allocate a single workspace at max_front to avoid
    // repeated reallocations as the postorder traversal encounters progressively
    // larger fronts. No pool needed — only one workspace exists.
    //
    // Parallel mode: use a pool of grow-on-demand workspaces. Pre-allocating
    // a separate shared_workspace alongside the pool would waste max_front²×8
    // bytes (685 MB for H2O) during parallel waves when only the pool is used.
    // Single-node batches (sequential fallback near the tree root) pop from
    // the pool instead.
    //
    // The pool is dropped when factor_tree_levelset returns, freeing all
    // workspace memory (unlike thread_local! which leaks for the lifetime
    // of the rayon thread pool).
    let mut seq_workspace = if !is_parallel {
        FactorizationWorkspace::new(max_front, n)
    } else {
        FactorizationWorkspace::empty()
    };
    let workspace_pool: std::sync::Mutex<Vec<FactorizationWorkspace>> =
        std::sync::Mutex::new(Vec::new());

    // Initial ready set: all leaf supernodes (no children) that are NOT in small-leaf subtrees
    let mut ready: Vec<usize> = (0..n_supernodes)
        .filter(|&s| remaining_children[s] == 0 && !supernodes[s].in_small_leaf)
        .collect();

    while !ready.is_empty() {
        // Collect child contributions for each ready node. Children of different
        // ready nodes are disjoint (tree structure), so each contribution is
        // taken exactly once.
        let mut batch_inputs: Vec<(usize, Vec<Option<ContributionBlock>>)> =
            Vec::with_capacity(ready.len());
        for &s in &ready {
            let child_contribs: Vec<Option<ContributionBlock>> = children[s]
                .iter()
                .map(|&c| contributions[c].take())
                .collect();
            batch_inputs.push((s, child_contribs));
        }

        // Factor all ready nodes (parallel if enabled and batch > 1)
        let batch_results: Vec<SupernodeResult> = if !is_parallel {
            // Pure sequential: use the pre-allocated workspace (no pool
            // overhead, no reallocations during postorder traversal).
            batch_inputs
                .into_iter()
                .map(|(s, contribs)| {
                    factor_single_supernode(
                        s,
                        &supernodes[s],
                        contribs,
                        matrix,
                        perm_fwd,
                        perm_inv,
                        options,
                        scaling,
                        &mut seq_workspace,
                        assembly_maps,
                    )
                })
                .collect::<Result<_, _>>()?
        } else if batch_inputs.len() == 1 {
            // Parallel mode, single-node batch (near tree root): pop a
            // workspace from the pool. Drop remaining pool workspaces to
            // free memory before processing large supernodes.
            let mut pool = workspace_pool.lock().unwrap();
            let mut ws = pool.pop().unwrap_or_else(FactorizationWorkspace::empty);
            pool.clear();
            drop(pool);
            ws.ensure_g2l(n);
            let results = batch_inputs
                .into_iter()
                .map(|(s, contribs)| {
                    factor_single_supernode(
                        s,
                        &supernodes[s],
                        contribs,
                        matrix,
                        perm_fwd,
                        perm_inv,
                        options,
                        scaling,
                        &mut ws,
                        assembly_maps,
                    )
                })
                .collect::<Result<_, _>>();
            workspace_pool.lock().unwrap().push(ws);
            results?
        } else {
            // Parallel: take workspaces from the caller-owned pool so each
            // rayon worker gets exclusive ownership during factorization.
            // The pool is reused across level-set waves and dropped when
            // factor_tree_levelset returns (unlike thread_local! which leaks
            // for the lifetime of the rayon thread pool).
            //
            // If rayon work-steals another par_iter task onto the same thread
            // (due to nested Par::rayon BLAS), that task takes a separate
            // workspace from the pool — no aliasing.
            use rayon::iter::{IntoParallelIterator, ParallelIterator};
            batch_inputs
                .into_par_iter()
                .map(|(s, contribs)| {
                    // Poisoned mutex means a worker panicked — already fatal.
                    let mut ws = workspace_pool
                        .lock()
                        .unwrap()
                        .pop()
                        .unwrap_or_else(FactorizationWorkspace::empty);
                    // Only ensure global_to_local is sized; frontal_data
                    // grows on demand in factor_single_supernode. Avoids
                    // eagerly allocating max_front² per thread (OOM on
                    // large-front matrices like H2O: 9258² × 8 = 685 MB).
                    ws.ensure_g2l(n);
                    let result = factor_single_supernode(
                        s,
                        &supernodes[s],
                        contribs,
                        matrix,
                        perm_fwd,
                        perm_inv,
                        options,
                        scaling,
                        &mut ws,
                        assembly_maps,
                    );
                    // Poisoned mutex means a worker panicked — already fatal.
                    workspace_pool.lock().unwrap().push(ws);
                    result
                })
                .collect::<Result<_, _>>()?
        };

        // Store results and determine next ready set
        let mut next_ready = Vec::new();
        for (result, &s) in batch_results.into_iter().zip(ready.iter()) {
            contributions[s] = result.contribution;
            all_results.push((s, result.ff, result.stats));

            if let Some(parent) = supernodes[s].parent {
                remaining_children[parent] -= 1;
                if remaining_children[parent] == 0 {
                    next_ready.push(parent);
                }
            }
        }

        ready = next_ready;
    }

    Ok(all_results)
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
#[allow(clippy::single_range_in_vec_init)]
fn build_supernode_info(symbolic: &AptpSymbolic) -> Vec<SupernodeInfo> {
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
                        owned_ranges: vec![begin[s]..end[s]],
                        in_small_leaf: false,
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
                        owned_ranges: vec![j..j + 1],
                        in_small_leaf: false,
                    }
                })
                .collect()
        }
    }
}

/// Build a map from each supernode to its children in the assembly tree.
fn build_children_map(infos: &[SupernodeInfo]) -> Vec<Vec<usize>> {
    let n = infos.len();
    let mut children = vec![Vec::new(); n];
    for (s, info) in infos.iter().enumerate() {
        if let Some(p) = info.parent {
            children[p].push(s);
        }
    }
    children
}

/// Classify post-amalgamation supernodes into small-leaf subtrees for the fast path.
///
/// Performs a single O(n_supernodes) bottom-up pass to identify contiguous leaf
/// subtrees where every supernode has `front_size < threshold`. Sets `in_small_leaf`
/// on qualifying supernodes and returns a `Vec<SmallLeafSubtree>` describing each
/// identified subtree (with at least 2 nodes).
///
/// # Arguments
///
/// - `supernodes`: Post-amalgamation supernode descriptors (mutated: `in_small_leaf` set).
/// - `children_map`: Children map (`children[s]` = list of child supernode indices).
/// - `threshold`: Front-size threshold. Supernodes with `front_size >= threshold` are
///   excluded. If `threshold == 0`, returns empty (fast path disabled).
///
/// # Algorithm
///
/// 1. Compute `front_size[s] = sum(owned_ranges[s].len()) + pattern[s].len()`.
/// 2. Bottom-up (postorder): mark `in_small_leaf[s] = true` iff `front_size[s] < threshold`
///    AND (s is a leaf OR all children have `in_small_leaf = true`).
/// 3. Identify subtree roots: `in_small_leaf[s] = true` and parent has `in_small_leaf = false`
///    (or no parent).
/// 4. Collect descendant nodes in postorder via iterative stack-based traversal.
/// 5. Filter subtrees with < 2 nodes.
///
/// # References
///
/// - SPRAL `SymbolicSubtree.hxx:57-84` (BSD-3): bottom-up classification
/// - Duff, Hogg & Lopez (2020), Algorithm 6.1: find_subtree_partition
fn classify_small_leaf_subtrees(
    supernodes: &mut [SupernodeInfo],
    children_map: &[Vec<usize>],
    threshold: usize,
) -> Vec<SmallLeafSubtree> {
    if threshold == 0 {
        return Vec::new();
    }

    let n_supernodes = supernodes.len();

    // Step 1-2: Bottom-up classification (supernodes are in postorder: s=0 is leftmost leaf)
    for s in 0..n_supernodes {
        let owned_cols: usize = supernodes[s].owned_ranges.iter().map(|r| r.len()).sum();
        let front_size = owned_cols + supernodes[s].pattern.len();

        if front_size >= threshold {
            supernodes[s].in_small_leaf = false;
            continue;
        }

        if children_map[s].is_empty() {
            // Leaf node with small front
            supernodes[s].in_small_leaf = true;
        } else if children_map[s].iter().all(|&c| supernodes[c].in_small_leaf) {
            // Interior node: all children are in_small_leaf and this node is small
            supernodes[s].in_small_leaf = true;
        } else {
            supernodes[s].in_small_leaf = false;
        }
    }

    // Step 3: Identify subtree roots
    let mut subtrees = Vec::new();
    for s in 0..n_supernodes {
        if !supernodes[s].in_small_leaf {
            continue;
        }
        let is_root = match supernodes[s].parent {
            None => true,
            Some(p) => !supernodes[p].in_small_leaf,
        };
        if !is_root {
            continue;
        }

        // Step 4: Collect descendant nodes in postorder via iterative DFS.
        // Two-stack approach: first stack drives DFS, second collects in reverse postorder.
        let mut nodes = Vec::new();
        let mut stack = vec![s];
        let mut visit_stack = Vec::new();
        while let Some(node) = stack.pop() {
            visit_stack.push(node);
            // Push children (they'll be processed before parent due to stack LIFO)
            for &c in &children_map[node] {
                if supernodes[c].in_small_leaf {
                    stack.push(c);
                }
            }
        }
        // visit_stack is in reverse postorder — reverse to get postorder
        while let Some(node) = visit_stack.pop() {
            nodes.push(node);
        }

        // Step 5: Filter subtrees with < 2 nodes
        if nodes.len() < 2 {
            // Single-node "subtree" — not worth the fast path overhead
            supernodes[s].in_small_leaf = false;
            continue;
        }

        // Compute max_front_size across all nodes in the subtree
        let max_front_size = nodes
            .iter()
            .map(|&node| {
                let owned: usize = supernodes[node].owned_ranges.iter().map(|r| r.len()).sum();
                owned + supernodes[node].pattern.len()
            })
            .max()
            .unwrap_or(0);

        subtrees.push(SmallLeafSubtree {
            root: s,
            nodes,
            max_front_size,
            parent_of_root: supernodes[s].parent,
        });
    }

    subtrees
}

/// Scatter original sparse matrix entries into a frontal matrix from multiple owned ranges.
///
/// After amalgamation, a supernode may own multiple non-contiguous column ranges.
/// This function iterates all owned columns and correctly handles upper-triangle
/// deduplication across ranges.
fn scatter_original_entries_multi(
    frontal: &mut FrontalMatrix<'_>,
    matrix: &SparseColMat<usize, f64>,
    perm_fwd: &[usize],
    perm_inv: &[usize],
    global_to_local: &[usize],
    owned_ranges: &[std::ops::Range<usize>],
    scaling: Option<&[f64]>,
) {
    let symbolic = matrix.symbolic();
    let col_ptrs = symbolic.col_ptr();
    let row_indices_csc = symbolic.row_idx();
    let values = matrix.val();
    let total_owned: usize = owned_ranges.iter().map(|r| r.len()).sum();
    let k = frontal.num_fully_summed;

    // Build a fast check for "is this permuted column an owned supernode column?"
    // We reuse global_to_local: if g2l[pj] < total_owned, it's an owned column
    // (since frontal_rows puts owned columns first). This works because owned
    // columns always occupy the first `total_owned` positions in frontal_rows.

    for range in owned_ranges {
        for pj in range.clone() {
            let orig_col = perm_fwd[pj];
            let local_col = global_to_local[pj];

            let start = col_ptrs[orig_col];
            let end = col_ptrs[orig_col + 1];
            for idx in start..end {
                let orig_row = row_indices_csc[idx];
                // Skip upper-triangle if both are owned supernode cols
                if orig_row < orig_col {
                    let perm_row = perm_inv[orig_row];
                    let local_peer = global_to_local[perm_row];
                    if local_peer != NOT_IN_FRONT && local_peer < total_owned {
                        continue; // both owned cols — avoid double-count
                    }
                }
                let global_row = perm_inv[orig_row];
                let local_row = global_to_local[global_row];
                if local_row == NOT_IN_FRONT {
                    continue;
                }
                // Skip delayed columns from children (indices total_owned..k)
                if local_row >= total_owned && local_row < k {
                    continue;
                }
                let mut val = values[idx];
                if let Some(s) = scaling {
                    val *= s[perm_inv[orig_row]] * s[perm_inv[orig_col]];
                }
                if local_row >= local_col {
                    frontal.data[(local_row, local_col)] += val;
                } else {
                    frontal.data[(local_col, local_row)] += val;
                }
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
    parent: &mut FrontalMatrix<'_>,
    child: ContributionBlock,
    global_to_local: &[usize],
) -> Mat<f64> {
    let cb_size = child.row_indices.len();
    // Column-oriented iteration: for each child column j, scatter the
    // lower-triangle entries (rows j..cb_size) into the parent. Source
    // reads are contiguous in column-major layout. We keep the zero-check
    // and symmetry swap since delayed columns can break monotonicity.
    for j in 0..cb_size {
        let gj = child.row_indices[j];
        let lj = global_to_local[gj];
        debug_assert!(
            lj != NOT_IN_FRONT,
            "extend_add: child col {} not in parent",
            gj
        );
        let child_col = &child.data.col_as_slice(j)[j..cb_size];
        let row_indices = &child.row_indices[j..cb_size];

        for (&val, &gi) in child_col.iter().zip(row_indices) {
            if val != 0.0 {
                let li = global_to_local[gi];
                debug_assert!(
                    li != NOT_IN_FRONT,
                    "extend_add: child row {} not in parent",
                    gi
                );
                // Map to parent: ensure lower triangle (li >= lj)
                if li >= lj {
                    parent.data[(li, lj)] += val;
                } else {
                    parent.data[(lj, li)] += val;
                }
            }
        }
    }
    // Return the consumed data buffer for recycling into workspace.contrib_buffer
    child.data
}

/// Extend-add using a precomputed child→parent row mapping.
///
/// Used when the child had zero delayed columns, so the contribution block
/// rows correspond exactly to the child's symbolic pattern and the precomputed
/// `ea_row_map` provides direct local index lookup without `global_to_local`.
///
/// Uses column-oriented processing: for each child
/// column `j`, scatters the lower-triangle entries `(j..cb_size, j)` into
/// the parent. Source reads are contiguous in column-major layout.
///
/// # References
///
/// - SPRAL `assemble.hxx:27-38` — `asm_col` column-oriented scatter (BSD-3)
///
/// # Arguments
///
/// - `parent`: Parent frontal matrix to accumulate into (lower triangle).
/// - `child`: Child contribution block.
/// - `ea_row_map`: Precomputed mapping from child contribution row index to
///   parent frontal local row index. Length equals the child's pattern size.
fn extend_add_mapped(
    parent: &mut FrontalMatrix<'_>,
    child: ContributionBlock,
    ea_row_map: &[u32],
) -> Mat<f64> {
    let cb_size = child.row_indices.len();
    debug_assert!(
        ea_row_map.len() >= cb_size,
        "extend_add_mapped: ea_row_map len {} < cb_size {}",
        ea_row_map.len(),
        cb_size
    );

    // Column-oriented scatter.
    // For each child column j, scatter lower-triangle entries (rows j..cb_size)
    // into the parent. Source reads via col_as_slice are contiguous in column-
    // major layout; writes target a single parent column when li >= lj (common
    // case), with a symmetry swap for the rare li < lj case.
    for j in 0..cb_size {
        let lj = ea_row_map[j] as usize;
        let child_col = child.data.col_as_slice(j);
        let row_map = &ea_row_map[j..cb_size];
        let src = &child_col[j..cb_size];

        // 4× unrolled inner loop for ILP
        let n = cb_size - j;
        let n4 = n / 4 * 4;
        let mut k = 0;
        while k < n4 {
            ea_scatter_one(&mut parent.data, row_map[k] as usize, lj, src[k]);
            ea_scatter_one(&mut parent.data, row_map[k + 1] as usize, lj, src[k + 1]);
            ea_scatter_one(&mut parent.data, row_map[k + 2] as usize, lj, src[k + 2]);
            ea_scatter_one(&mut parent.data, row_map[k + 3] as usize, lj, src[k + 3]);
            k += 4;
        }
        while k < n {
            ea_scatter_one(&mut parent.data, row_map[k] as usize, lj, src[k]);
            k += 1;
        }
    }
    // Return the consumed data buffer for recycling into workspace.contrib_buffer
    child.data
}

/// Scatter a single extend-add value into the parent's lower triangle.
#[inline(always)]
fn ea_scatter_one(parent: &mut MatMut<'_, f64>, li: usize, lj: usize, val: f64) {
    if li >= lj {
        parent[(li, lj)] += val;
    } else {
        parent[(lj, li)] += val;
    }
}

/// Extract per-supernode factors from a factored frontal matrix.
///
/// Copies L11, D11, L21, and permutation information from the in-place
/// factored frontal matrix into a persistent [`FrontFactors`] struct.
///
/// Uses per-column `col_as_slice` + `copy_from_slice` for L21 extraction and
/// per-column slice copies for L11 1x1-pivot columns, replacing element-by-element
/// indexing with bulk memory operations.
///
/// # Arguments
///
/// - `frontal_data`: The factored frontal matrix data (owned `Mat<f64>`).
///   Columns must be contiguous (guaranteed by `Mat` layout). Only the first
///   `m` rows and columns are used.
/// - `m`: Front size (number of rows/columns used in `frontal_data`).
/// - `k`: Number of fully-summed rows/columns.
/// - `frontal_row_indices`: Global permuted column indices for each local row (length >= m).
/// - `result`: APTP factorization result (pivot types, permutation, etc.).
pub(crate) fn extract_front_factors(
    frontal_data: &Mat<f64>,
    m: usize,
    k: usize,
    frontal_row_indices: &[usize],
    result: &AptpFactorResult,
) -> FrontFactors {
    let ne = result.num_eliminated;

    // Extract L11 (ne × ne) — unit lower triangular part
    // For 2x2 pivots, a[(col+1, col)] is the D off-diagonal, NOT an L entry.
    // Uses per-column slice copies for 1x1 pivot columns where the L entries
    // form a contiguous segment.
    let l11 = if ne > 0 {
        let mut l = Mat::zeros(ne, ne);
        let mut col = 0;
        while col < ne {
            l[(col, col)] = 1.0; // unit diagonal
            match result.d.pivot_type(col) {
                PivotType::OneByOne => {
                    // L entries: rows (col+1)..ne — contiguous in both source and dest
                    let n_entries = ne - (col + 1);
                    if n_entries > 0 {
                        let src = &frontal_data.col_as_slice(col)[col + 1..ne];
                        l.col_as_slice_mut(col)[col + 1..ne].copy_from_slice(src);
                    }
                    col += 1;
                }
                PivotType::TwoByTwo { partner } if partner > col => {
                    l[(col + 1, col + 1)] = 1.0; // unit diagonal for partner
                    // a[(col+1, col)] is D off-diagonal, skip it
                    // L entries for both columns start at row col+2
                    let n_entries = ne - (col + 2);
                    if n_entries > 0 {
                        let src0 = &frontal_data.col_as_slice(col)[col + 2..ne];
                        l.col_as_slice_mut(col)[col + 2..ne].copy_from_slice(src0);
                        let src1 = &frontal_data.col_as_slice(col + 1)[col + 2..ne];
                        l.col_as_slice_mut(col + 1)[col + 2..ne].copy_from_slice(src1);
                    }
                    col += 2;
                }
                PivotType::TwoByTwo { .. } => {
                    // Second column of a 2x2 pair — already handled above
                    col += 1;
                }
                PivotType::Delayed => {
                    // Should not appear for columns 0..ne (they were eliminated)
                    unreachable!("unexpected Delayed pivot at col {} in 0..ne", col);
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
        match result.d.pivot_type(col) {
            PivotType::OneByOne => {
                d11.set_1x1(col, result.d.diagonal_1x1(col));
                col += 1;
            }
            PivotType::TwoByTwo { partner: _ } => {
                let block = result.d.diagonal_2x2(col);
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
                unreachable!("unexpected Delayed pivot at col {} in 0..ne", col);
            }
        }
    }

    // Extract L21 (r × ne) where r = m - ne (includes delayed rows ne..k AND
    // non-fully-summed rows k..m). The factored frontal has valid L entries at rows
    // ne..m for columns 0..ne: rows ne..k are delayed fully-summed rows that received
    // TRSM updates before being delayed; rows k..m are non-fully-summed rows.
    //
    // Uses per-column copy_from_slice for bulk extraction.
    //
    // NOTE: extract_contribution also includes delayed rows (see lines below).
    // These two functions MUST be consistent in their row treatment.
    let r = m - ne;
    let l21 = if ne > 0 && r > 0 {
        let mut l = Mat::zeros(r, ne);
        for j in 0..ne {
            let src = &frontal_data.col_as_slice(j)[ne..m];
            l.col_as_slice_mut(j)[..r].copy_from_slice(src);
        }
        l
    } else {
        Mat::zeros(r, ne)
    };

    // Local permutation (maps factored position to original front-local column)
    let local_perm = result.perm[..k].to_vec();

    // Column indices: the global permuted indices of the eliminated columns
    // local_perm[0..ne] gives the front-local columns that were eliminated,
    // mapped through frontal_row_indices to get global permuted indices
    let col_indices: Vec<usize> = local_perm[..ne]
        .iter()
        .map(|&lp| frontal_row_indices[lp])
        .collect();

    // Row indices: global permuted indices for L21 rows.
    // First num_delayed entries: delayed fully-summed columns (ne..k), mapped through
    // result.perm and frontal_row_indices to get global permuted indices.
    // Remaining entries: non-fully-summed rows (k..m).
    // This matches extract_contribution's row_indices construction.
    let mut row_indices = Vec::with_capacity(m - ne);
    for &lp in &result.perm[ne..k] {
        row_indices.push(frontal_row_indices[lp]);
    }
    row_indices.extend_from_slice(&frontal_row_indices[k..m]);

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
/// The NFS×NFS Schur complement data is already in `contrib_buffer` (written
/// by [`compute_contribution_gemm`]). This function adds any delayed-column
/// data from the workspace and builds the index metadata.
///
/// # Layout within returned `ContributionBlock.data`
///
/// ```text
/// [0..num_delayed, 0..num_delayed]     → delayed × delayed (from workspace)
/// [num_delayed..size, 0..num_delayed]  → NFS × delayed cross-terms (from workspace)
/// [num_delayed..size, num_delayed..size] → NFS × NFS Schur complement (from contrib_buffer)
/// ```
///
/// # Arguments
///
/// - `frontal_data`: The factored frontal matrix data (owned `Mat<f64>`).
/// - `m`: Front size (number of rows/columns used).
/// - `k`: Number of fully-summed rows/columns.
/// - `frontal_row_indices`: Global permuted column indices for each local row.
/// - `result`: APTP factorization result.
/// - `contrib_buffer`: Pre-allocated buffer already containing NFS×NFS data
///   in `[0..nfs, 0..nfs]` (from deferred GEMM). Moved in; becomes the
///   `ContributionBlock.data`.
pub(crate) fn extract_contribution(
    frontal_data: &Mat<f64>,
    m: usize,
    k: usize,
    frontal_row_indices: &[usize],
    result: &AptpFactorResult,
    mut contrib_buffer: Mat<f64>,
) -> ContributionBlock {
    let ne = result.num_eliminated;
    let size = m - ne;
    let num_delayed = k - ne;
    let nfs = m - k; // = size - num_delayed

    // The NFS×NFS region (contrib_buffer[0..nfs, 0..nfs]) is already filled by
    // compute_contribution_gemm. We need to handle delayed columns (if any) by
    // shifting the NFS×NFS data to make room for the delayed portion.
    if num_delayed > 0 {
        // Delayed columns exist: we need to assemble the full (size × size) contribution.
        // Layout: [delayed × delayed | NFS × delayed | NFS × NFS]
        //
        // Strategy: allocate a new buffer of the correct size, copy delayed regions
        // from workspace, and NFS×NFS from contrib_buffer.
        let mut data = Mat::zeros(size, size);

        // Copy delayed × delayed from workspace: frontal_data[ne..k, ne..k]
        // These positions are in the APTP-permuted order. The delayed columns
        // are at positions ne..k in the factored frontal matrix.
        for j in 0..num_delayed {
            let col_len = num_delayed - j;
            let src = &frontal_data.col_as_slice(ne + j)[ne + j..ne + j + col_len];
            data.col_as_slice_mut(j)[j..j + col_len].copy_from_slice(src);
        }

        // Copy NFS × delayed cross-terms from workspace: frontal_data[k..m, ne..k]
        for j in 0..num_delayed {
            let src = &frontal_data.col_as_slice(ne + j)[k..m];
            data.col_as_slice_mut(j)[num_delayed..size].copy_from_slice(src);
        }

        // Copy NFS × NFS from contrib_buffer[0..nfs, 0..nfs] into data[num_delayed..size, num_delayed..size]
        for j in 0..nfs {
            let col_len = nfs - j;
            let src = &contrib_buffer.col_as_slice(j)[j..j + col_len];
            data.col_as_slice_mut(num_delayed + j)[num_delayed + j..size].copy_from_slice(src);
        }

        contrib_buffer = data;
    }
    // else: num_delayed == 0 → contrib_buffer already has the full contribution
    // (NFS×NFS = entire contribution), no copy needed. True zero-copy path.

    // Build row indices:
    // - First num_delayed entries: delayed columns mapped through result.perm and frontal_row_indices
    // - Remaining entries: non-fully-summed rows from frontal_row_indices[k..m]
    let mut row_indices = Vec::with_capacity(size);

    // Delayed columns: these are at positions ne..k in the APTP-permuted frontal matrix
    // result.perm[ne..k] gives the original front-local indices of delayed columns
    for &lp in &result.perm[ne..k] {
        row_indices.push(frontal_row_indices[lp]);
    }

    // Non-fully-summed rows
    row_indices.extend_from_slice(&frontal_row_indices[k..m]);

    ContributionBlock {
        data: contrib_buffer,
        row_indices,
        num_delayed,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symmetric::factor::AptpOptions;

    /// Helper: build frontal matrix data from a dense lower-triangle matrix.
    ///
    /// Returns owned `(data, row_indices)` — construct a `FrontalMatrix` view from these.
    /// `k` = number of fully-summed columns.
    /// `row_indices` = global permuted indices for each local row (length m).
    fn make_frontal_data(
        lower: &Mat<f64>,
        k: usize,
        row_indices: Vec<usize>,
    ) -> (Mat<f64>, Vec<usize>, usize) {
        let m = lower.nrows();
        assert_eq!(lower.ncols(), m);
        assert_eq!(row_indices.len(), m);
        assert!(k <= m);

        // Build full symmetric from lower triangle
        let mut data = Mat::zeros(m, m);
        for i in 0..m {
            for j in 0..=i {
                data[(i, j)] = lower[(i, j)];
                data[(j, i)] = lower[(i, j)];
            }
        }

        (data, row_indices, k)
    }

    /// Helper: build a symmetric indefinite matrix (lower triangle) that forces
    /// delays. Strategy from SPRAL's `cause_delays`: scale selected rows/cols
    /// by a large factor, making their off-diagonal entries dominate the diagonal
    /// so the threshold test fails.
    ///
    /// Returns (matrix, expected_to_have_delays).
    fn make_delay_matrix(m: usize, k: usize, delay_rows: &[usize]) -> Mat<f64> {
        let mut a = Mat::zeros(m, m);

        // Start with a diagonally dominant indefinite matrix
        for i in 0..m {
            // Alternate positive/negative diagonal for indefiniteness
            a[(i, i)] = if i % 2 == 0 { 4.0 } else { -4.0 };
            // Add some off-diagonal coupling
            for j in 0..i {
                let val = 0.5 / ((i - j) as f64 + 1.0);
                a[(i, j)] = val;
            }
        }

        // Scale delay_rows to create threshold failures:
        // Multiply row/col by 1000 but NOT the diagonal, so off-diag >> diag
        for &r in delay_rows {
            if r < k {
                // Scale off-diagonal entries in row r
                for j in 0..r {
                    a[(r, j)] *= 1000.0;
                }
                for i in (r + 1)..m {
                    a[(i, r)] *= 1000.0;
                }
                // DON'T scale diagonal — this makes |a_rr| << |a_ir| => threshold failure
            }
        }

        a
    }

    /// Test that extract_front_factors produces correct L21 dimensions when
    /// delays are present: L21 should have (m - ne) rows, NOT (m - k) rows.
    /// This is the regression test for the bug fixed in this commit.
    #[test]
    fn test_extract_front_factors_l21_includes_delayed_rows() {
        let m = 12;
        let k = 8;
        let delay_rows = &[1, 5]; // Columns 1 and 5 should be delayed

        let lower = make_delay_matrix(m, k, delay_rows);
        let row_indices: Vec<usize> = (100..100 + m).collect(); // arbitrary global indices

        let (mut data, row_indices, k) = make_frontal_data(&lower, k, row_indices);

        let options = AptpOptions::default();
        let result =
            aptp_factor_in_place(data.as_mut(), k, &options).expect("factor should succeed");

        let ne = result.num_eliminated;
        let num_delayed = k - ne;

        // The bug: if delays occurred, the old code would make L21 have (m - k) rows
        // instead of (m - ne) rows, missing the delayed-row entries.
        if num_delayed > 0 {
            let ff = extract_front_factors(&data, m, k, &row_indices, &result);

            // L21 must have (m - ne) rows = delayed_rows + non-fully-summed rows
            assert_eq!(
                ff.l21.nrows(),
                m - ne,
                "L21 should have m-ne={} rows (not m-k={}), ne={}, delays={}",
                m - ne,
                m - k,
                ne,
                num_delayed
            );

            // row_indices must match L21 rows
            assert_eq!(
                ff.row_indices.len(),
                ff.l21.nrows(),
                "row_indices length must match L21 rows"
            );

            // col_indices must have ne entries
            assert_eq!(ff.col_indices.len(), ne);
        } else {
            // No delays — still verify basic shapes
            let ff = extract_front_factors(&data, m, k, &row_indices, &result);
            assert_eq!(ff.l21.nrows(), m - ne);
            assert_eq!(ff.row_indices.len(), ff.l21.nrows());
        }
    }

    /// Test that delayed rows in extract_front_factors and extract_contribution
    /// have consistent row_indices. Both functions MUST agree on which global
    /// indices correspond to delayed rows.
    #[test]
    fn test_extract_front_factors_contribution_row_indices_consistent() {
        let m = 10;
        let k = 6;
        let delay_rows = &[2]; // Force column 2 to be delayed

        let lower = make_delay_matrix(m, k, delay_rows);
        let row_indices: Vec<usize> = (200..200 + m).collect();

        let (mut data, row_indices, k) = make_frontal_data(&lower, k, row_indices);

        let options = AptpOptions::default();
        let result =
            aptp_factor_in_place(data.as_mut(), k, &options).expect("factor should succeed");

        let ne = result.num_eliminated;
        let num_delayed = k - ne;

        if num_delayed > 0 {
            let ff = extract_front_factors(&data, m, k, &row_indices, &result);
            let nfs = m - k;
            let mut contrib_buf = Mat::zeros(m, m);
            if nfs > 0 {
                let mut ld_ws = Mat::new();
                compute_contribution_gemm(
                    &data,
                    k,
                    ne,
                    m,
                    &result.d,
                    &mut contrib_buf,
                    &mut ld_ws,
                    Par::Seq,
                );
            }
            let cb = extract_contribution(&data, m, k, &row_indices, &result, contrib_buf);

            // Both must have the same number of delayed rows
            assert_eq!(cb.num_delayed, num_delayed);

            // The first num_delayed entries of row_indices must match
            let ff_delayed_indices = &ff.row_indices[..num_delayed];
            let cb_delayed_indices = &cb.row_indices[..num_delayed];
            assert_eq!(
                ff_delayed_indices, cb_delayed_indices,
                "Delayed row indices must be identical between extract_front_factors \
                 and extract_contribution.\n  front_factors: {:?}\n  contribution:  {:?}",
                ff_delayed_indices, cb_delayed_indices
            );

            // The non-fully-summed row indices must also match
            let ff_nfs_indices = &ff.row_indices[num_delayed..];
            let cb_nfs_indices = &cb.row_indices[num_delayed..];
            assert_eq!(
                ff_nfs_indices, cb_nfs_indices,
                "Non-fully-summed row indices must match"
            );
        }
    }

    /// Test that L21 entries at delayed row positions are non-zero.
    /// When APTP delays a column, the TRSM (apply_and_check) has already
    /// computed valid L entries at those row positions. These must be preserved
    /// in the extracted L21.
    #[test]
    fn test_extract_front_factors_delayed_l21_entries_populated() {
        let m = 12;
        let k = 8;
        // Use multiple delay candidates to increase chance of delays
        let delay_rows = &[0, 3, 5];

        let lower = make_delay_matrix(m, k, delay_rows);
        let row_indices: Vec<usize> = (0..m).collect();

        let (mut data, row_indices, k) = make_frontal_data(&lower, k, row_indices);

        let options = AptpOptions::default();
        let result =
            aptp_factor_in_place(data.as_mut(), k, &options).expect("factor should succeed");

        let ne = result.num_eliminated;
        let num_delayed = k - ne;

        if num_delayed > 0 && ne > 0 {
            let ff = extract_front_factors(&data, m, k, &row_indices, &result);

            // Check that L21 entries in the delayed-row portion (rows 0..num_delayed)
            // are not all zero. These rows received TRSM updates before being delayed.
            let delayed_l21 = ff.l21.as_ref().subrows(0, num_delayed);
            let mut has_nonzero = false;
            for i in 0..delayed_l21.nrows() {
                for j in 0..delayed_l21.ncols() {
                    if delayed_l21[(i, j)].abs() > 1e-16 {
                        has_nonzero = true;
                    }
                }
            }
            assert!(
                has_nonzero,
                "L21 delayed rows should have non-zero entries from TRSM. \
                 ne={}, delays={}, L21 shape=({},{})",
                ne,
                num_delayed,
                ff.l21.nrows(),
                ff.l21.ncols()
            );
        }
    }

    /// Reconstruction test for extract_front_factors with delays.
    /// Verify that the extracted L11, D11, L21 correctly reconstruct the
    /// factored portion of the frontal matrix: PAP^T[0..ne, 0..ne] = L11 * D11 * L11^T.
    #[test]
    fn test_extract_front_factors_reconstruction_with_delays() {
        let m = 16;
        let k = 10;
        let delay_rows = &[1, 4, 7];

        let lower = make_delay_matrix(m, k, delay_rows);
        let row_indices: Vec<usize> = (0..m).collect();

        // Save a copy of the original symmetric matrix before factoring
        let mut original = Mat::zeros(m, m);
        for i in 0..m {
            for j in 0..=i {
                original[(i, j)] = lower[(i, j)];
                original[(j, i)] = lower[(i, j)];
            }
        }

        let (mut data, row_indices, k) = make_frontal_data(&lower, k, row_indices);

        let options = AptpOptions::default();
        let result =
            aptp_factor_in_place(data.as_mut(), k, &options).expect("factor should succeed");

        let ne = result.num_eliminated;
        if ne == 0 {
            return; // Nothing to reconstruct
        }

        let ff = extract_front_factors(&data, m, k, &row_indices, &result);

        // Compute L11 * D * L11^T manually
        // Step 1: LD = L11 * D (column-by-column, handling 1x1 and 2x2 pivots)
        let mut ld = Mat::zeros(ne, ne);
        let mut col = 0;
        while col < ne {
            match ff.d11.pivot_type(col) {
                PivotType::OneByOne => {
                    let d = ff.d11.diagonal_1x1(col);
                    for i in 0..ne {
                        ld[(i, col)] = ff.l11[(i, col)] * d;
                    }
                    col += 1;
                }
                PivotType::TwoByTwo { .. } => {
                    let block = ff.d11.diagonal_2x2(col);
                    for i in 0..ne {
                        let l0 = ff.l11[(i, col)];
                        let l1 = ff.l11[(i, col + 1)];
                        ld[(i, col)] = l0 * block.a + l1 * block.b;
                        ld[(i, col + 1)] = l0 * block.b + l1 * block.c;
                    }
                    col += 2;
                }
                PivotType::Delayed => {
                    col += 1;
                }
            }
        }

        // Step 2: reconstructed = LD * L11^T
        let mut reconstructed = Mat::zeros(ne, ne);
        for i in 0..ne {
            for j in 0..=i {
                let mut val = 0.0;
                for p in 0..ne {
                    val += ld[(i, p)] * ff.l11[(j, p)];
                }
                reconstructed[(i, j)] = val;
                reconstructed[(j, i)] = val;
            }
        }

        // Build the permuted original submatrix for comparison
        let mut perm_original = Mat::zeros(ne, ne);
        for i in 0..ne {
            for j in 0..=i {
                let oi = result.perm[i];
                let oj = result.perm[j];
                let val = if oi >= oj {
                    original[(oi, oj)]
                } else {
                    original[(oj, oi)]
                };
                perm_original[(i, j)] = val;
                perm_original[(j, i)] = val;
            }
        }

        // Check reconstruction: ||PAP^T[0..ne,0..ne] - L11*D11*L11^T|| / ||PAP^T||
        let mut diff_norm = 0.0f64;
        let mut orig_norm = 0.0f64;
        for i in 0..ne {
            for j in 0..=i {
                let d = reconstructed[(i, j)] - perm_original[(i, j)];
                diff_norm += d * d;
                orig_norm += perm_original[(i, j)] * perm_original[(i, j)];
            }
        }
        let rel_err = diff_norm.sqrt() / orig_norm.sqrt().max(1e-15);

        assert!(
            rel_err < 1e-10,
            "L11*D11*L11^T reconstruction error: {:.2e} (ne={}, delays={})",
            rel_err,
            ne,
            k - ne
        );
    }

    /// Verify that the solve path produces correct results when delays are
    /// present. Uses a single dense frontal matrix (simulating one supernode),
    /// factors with delays, extracts front_factors, and runs forward/diagonal/
    /// backward solve steps to verify Ax = b.
    ///
    /// This directly tests the code path that the extract_front_factors bug
    /// would break: delayed-row L21 entries participate in forward and backward
    /// solve via gather/scatter operations.
    #[test]
    fn test_extract_front_factors_solve_roundtrip_with_delays() {
        let m = 12;
        let k = 8; // fully-summed columns

        let lower = make_delay_matrix(m, k, &[0, 3, 5]);
        let row_indices: Vec<usize> = (0..m).collect();

        // Build full symmetric matrix for reference
        let mut a_full = Mat::zeros(m, m);
        for i in 0..m {
            for j in 0..=i {
                a_full[(i, j)] = lower[(i, j)];
                a_full[(j, i)] = lower[(i, j)];
            }
        }

        let (mut data, row_indices, k) = make_frontal_data(&lower, k, row_indices);

        let options = AptpOptions::default();
        let result =
            aptp_factor_in_place(data.as_mut(), k, &options).expect("factor should succeed");

        let ne = result.num_eliminated;
        let num_delayed = k - ne;

        // Only meaningful if we actually got some delays and some eliminations
        if num_delayed == 0 || ne == 0 {
            return;
        }

        let ff = extract_front_factors(&data, m, k, &row_indices, &result);

        // Check that L11, D11, L21 are internally consistent:
        // L * D * L^T (where L = [L11; L21]) should reconstruct the factored
        // portion of A when permuted.

        // Build full L = [L11; L21] and check L * D * L^T = P * A[0..m, 0..ne] * P^T
        let r = ff.l21.nrows();
        let full_l_rows = ne + r;
        assert_eq!(full_l_rows, m, "L should cover all rows");

        // Compute L*D*L^T for the eliminated portion
        let mut ld = Mat::zeros(full_l_rows, ne);

        // L = [L11; L21]
        let mut full_l = Mat::zeros(full_l_rows, ne);
        for i in 0..ne {
            for j in 0..ne {
                full_l[(i, j)] = ff.l11[(i, j)];
            }
        }
        for i in 0..r {
            for j in 0..ne {
                full_l[(ne + i, j)] = ff.l21[(i, j)];
            }
        }

        // LD = L * D (handling 1x1 and 2x2)
        let mut col = 0;
        while col < ne {
            match ff.d11.pivot_type(col) {
                PivotType::OneByOne => {
                    let d = ff.d11.diagonal_1x1(col);
                    for i in 0..full_l_rows {
                        ld[(i, col)] = full_l[(i, col)] * d;
                    }
                    col += 1;
                }
                PivotType::TwoByTwo { .. } => {
                    let block = ff.d11.diagonal_2x2(col);
                    for i in 0..full_l_rows {
                        let l0 = full_l[(i, col)];
                        let l1 = full_l[(i, col + 1)];
                        ld[(i, col)] = l0 * block.a + l1 * block.b;
                        ld[(i, col + 1)] = l0 * block.b + l1 * block.c;
                    }
                    col += 2;
                }
                PivotType::Delayed => {
                    col += 1;
                }
            }
        }

        // reconstructed = LD * L^T
        let mut reconstructed = Mat::zeros(full_l_rows, full_l_rows);
        for i in 0..full_l_rows {
            for j in 0..=i {
                let mut val = 0.0;
                for p in 0..ne {
                    val += ld[(i, p)] * full_l[(j, p)];
                }
                reconstructed[(i, j)] = val;
                reconstructed[(j, i)] = val;
            }
        }

        // Build permuted original: row_order = col_indices ++ row_indices
        let mut perm_rows = Vec::with_capacity(full_l_rows);
        perm_rows.extend_from_slice(ff.col_indices());
        perm_rows.extend_from_slice(ff.row_indices());

        // Check that L*D*L^T[0..ne, 0..ne] ≈ A_perm[0..ne, 0..ne]
        let mut diff_norm = 0.0f64;
        let mut orig_norm = 0.0f64;
        for i in 0..ne {
            for j in 0..=i {
                let gi = perm_rows[i];
                let gj = perm_rows[j];
                let orig_val = if gi >= gj {
                    a_full[(gi, gj)]
                } else {
                    a_full[(gj, gi)]
                };
                let d = reconstructed[(i, j)] - orig_val;
                diff_norm += d * d;
                orig_norm += orig_val * orig_val;
            }
        }
        let rel_err = diff_norm.sqrt() / orig_norm.sqrt().max(1e-15);
        assert!(
            rel_err < 1e-10,
            "Full L*D*L^T reconstruction error: {:.2e} (ne={}, delays={}, L21_rows={})",
            rel_err,
            ne,
            num_delayed,
            r
        );
    }

    /// Test scatter_original_entries_multi with non-contiguous owned_ranges.
    ///
    /// After amalgamation, a merged supernode may own columns from multiple
    /// non-contiguous ranges (e.g., [0..2, 4..6] when child columns 2,3 were
    /// NOT merged). This test verifies:
    ///
    /// 1. Entries from all owned ranges are scattered into the frontal matrix
    /// 2. Upper-triangle deduplication works correctly across non-contiguous
    ///    ranges (an entry with both row and col in owned ranges is only counted once)
    /// 3. Entries where one index is owned and the other is a non-fully-summed
    ///    row are scattered correctly
    #[test]
    fn test_scatter_original_entries_multi_non_contiguous() {
        use faer::sparse::{SparseColMat, Triplet};

        // Build a small 6×6 symmetric matrix with identity permutation.
        // Owned ranges: [0..2, 4..6] (columns 0,1,4,5 owned; 2,3 not owned).
        //
        // The frontal has 6 rows: locals 0..4 are owned (4 total_owned),
        // no delayed children, rows 4..6 are non-fully-summed (rows for
        // global indices 2 and 3).
        //
        // Matrix (lower triangle values, symmetric):
        //     col: 0    1    2    3    4    5
        // 0: [10.0                          ]
        // 1: [ 1.0  11.0                    ]
        // 2: [ 2.0   0.0  12.0              ]
        // 3: [ 0.0   3.0   0.0  13.0        ]
        // 4: [ 4.0   0.0   0.0   0.0  14.0  ]
        // 5: [ 0.0   5.0   6.0   0.0   7.0  15.0]
        let triplets = vec![
            // Diagonal
            Triplet::new(0, 0, 10.0),
            Triplet::new(1, 1, 11.0),
            Triplet::new(2, 2, 12.0),
            Triplet::new(3, 3, 13.0),
            Triplet::new(4, 4, 14.0),
            Triplet::new(5, 5, 15.0),
            // Off-diagonal (both triangles for full CSC)
            Triplet::new(1, 0, 1.0),
            Triplet::new(0, 1, 1.0), // mirror
            Triplet::new(2, 0, 2.0),
            Triplet::new(0, 2, 2.0), // mirror
            Triplet::new(3, 1, 3.0),
            Triplet::new(1, 3, 3.0), // mirror
            Triplet::new(4, 0, 4.0),
            Triplet::new(0, 4, 4.0), // mirror
            Triplet::new(5, 1, 5.0),
            Triplet::new(1, 5, 5.0), // mirror
            Triplet::new(5, 2, 6.0),
            Triplet::new(2, 5, 6.0), // mirror
            Triplet::new(5, 4, 7.0),
            Triplet::new(4, 5, 7.0), // mirror
        ];
        let matrix = SparseColMat::try_new_from_triplets(6, 6, &triplets).expect("valid CSC");

        // Identity permutation
        let perm_fwd: Vec<usize> = (0..6).collect();
        let perm_inv: Vec<usize> = (0..6).collect();

        // Frontal rows: owned columns first (0,1,4,5), then non-owned rows (2,3)
        let frontal_rows = vec![0, 1, 4, 5, 2, 3];
        let total_owned = 4; // columns 0,1,4,5
        let m = frontal_rows.len();
        let k = total_owned; // num_fully_summed = total_owned (no delayed children)

        // Build global_to_local mapping
        let mut global_to_local = vec![NOT_IN_FRONT; 6];
        for (local, &global) in frontal_rows.iter().enumerate() {
            global_to_local[global] = local;
        }

        let mut frontal_data = Mat::zeros(m, m);

        // Non-contiguous owned ranges (the key case)
        let owned_ranges = vec![0..2, 4..6];

        {
            let mut frontal = FrontalMatrix {
                data: frontal_data.as_mut(),
                row_indices: &frontal_rows,
                num_fully_summed: k,
            };

            scatter_original_entries_multi(
                &mut frontal,
                &matrix,
                &perm_fwd,
                &perm_inv,
                &global_to_local,
                &owned_ranges,
                None,
            );
        }

        // Verify scattered values.
        // Local indices: 0→global 0, 1→global 1, 2→global 4, 3→global 5, 4→global 2, 5→global 3
        //
        // Expected entries in the frontal (lower triangle):
        //   (0,0)=10.0  diagonal of global 0
        //   (1,0)=1.0   global(1,0) → local(1,0)
        //   (1,1)=11.0  diagonal of global 1
        //   (2,0)=4.0   global(4,0) → local(2,0)
        //   (2,2)=14.0  diagonal of global 4
        //   (3,1)=5.0   global(5,1) → local(3,1)
        //   (3,2)=7.0   global(5,4) → local(3,2)
        //   (3,3)=15.0  diagonal of global 5
        //   (4,0)=2.0   global(2,0) → local(4,0)  [non-owned row, owned col]
        //   (4,3)=6.0   global(2,5) → local(4,3)  [non-owned row, owned col]
        //   (5,1)=3.0   global(3,1) → local(5,1)  [non-owned row, owned col]
        //
        // Upper-triangle dedup: global(0,1) should be skipped because both
        // 0 and 1 are owned. Similarly global(0,4), global(1,5), global(4,5)
        // are upper-triangle entries between owned cols and should be skipped.

        let f = &frontal_data;

        // Diagonal entries
        assert_eq!(f[(0, 0)], 10.0, "diag global 0");
        assert_eq!(f[(1, 1)], 11.0, "diag global 1");
        assert_eq!(f[(2, 2)], 14.0, "diag global 4");
        assert_eq!(f[(3, 3)], 15.0, "diag global 5");

        // Off-diagonal between owned cols (lower triangle only)
        assert_eq!(f[(1, 0)], 1.0, "global(1,0) → local(1,0)");
        assert_eq!(f[(2, 0)], 4.0, "global(4,0) → local(2,0)");
        assert_eq!(f[(3, 1)], 5.0, "global(5,1) → local(3,1)");
        assert_eq!(f[(3, 2)], 7.0, "global(5,4) → local(3,2)");

        // Non-owned rows (global 2 → local 4, global 3 → local 5)
        assert_eq!(f[(4, 0)], 2.0, "global(2,0) → local(4,0)");
        assert_eq!(f[(4, 3)], 6.0, "global(2,5) → local(4,3)");
        assert_eq!(f[(5, 1)], 3.0, "global(3,1) → local(5,1)");

        // Verify no spurious double-counting: the (0,1) position should be 1.0,
        // not 2.0 (which would happen if upper-triangle entry was also scattered)
        assert_eq!(
            f[(1, 0)],
            1.0,
            "upper-triangle dedup: (1,0) should be 1.0, not 2.0"
        );

        // Cross-range upper-triangle dedup: global(4,0) upper entry global(0,4)
        // should NOT be scattered separately, so local(2,0) should be 4.0
        assert_eq!(
            f[(2, 0)],
            4.0,
            "cross-range upper-triangle dedup: (2,0) should be 4.0, not 8.0"
        );

        // Verify zeros where no entries exist
        assert_eq!(f[(4, 1)], 0.0, "no global(2,1) entry");
        assert_eq!(f[(4, 2)], 0.0, "no global(2,4) entry");
        assert_eq!(f[(5, 0)], 0.0, "no global(3,0) entry");
        assert_eq!(f[(5, 2)], 0.0, "no global(3,4) entry");
        assert_eq!(f[(5, 3)], 0.0, "no global(3,5) entry");
    }

    // -----------------------------------------------------------------------
    // extend_add_mapped: precomputed-map scatter (zero-delay fast path)
    // -----------------------------------------------------------------------

    /// Test that `extend_add_mapped` correctly scatters a child contribution
    /// into a parent frontal matrix using a precomputed u32 row mapping,
    /// handling the lower-triangle symmetry swap when child ordering differs
    /// from parent ordering (non-monotonic map from amalgamated supernodes).
    #[test]
    fn test_extend_add_mapped_hand_constructed() {
        // Parent frontal: 5×5, initially zero (represents assembled parent).
        let mut parent_data = Mat::<f64>::zeros(5, 5);
        let parent_row_indices = vec![10, 20, 30, 40, 50]; // global indices (unused by mapped path)
        let mut parent = FrontalMatrix {
            data: parent_data.as_mut(),
            row_indices: &parent_row_indices,
            num_fully_summed: 3,
        };

        // Child contribution: 3×3 lower triangle.
        //   [1.0          ]
        //   [2.0  3.0     ]
        //   [4.0  5.0  6.0]
        let mut child_data = Mat::<f64>::zeros(3, 3);
        child_data[(0, 0)] = 1.0;
        child_data[(1, 0)] = 2.0;
        child_data[(1, 1)] = 3.0;
        child_data[(2, 0)] = 4.0;
        child_data[(2, 1)] = 5.0;
        child_data[(2, 2)] = 6.0;

        let child = ContributionBlock {
            data: child_data,
            row_indices: vec![30, 10, 40], // global indices
            num_delayed: 0,
        };

        // Precomputed map: child row 0 → parent local 2 (global 30),
        //                  child row 1 → parent local 0 (global 10),
        //                  child row 2 → parent local 3 (global 40).
        // Note: child row 1 maps to a LOWER parent index than child row 0,
        // so extend_add_mapped must handle the li < lj swap case.
        let ea_row_map: Vec<u32> = vec![2, 0, 3];

        let recycled = extend_add_mapped(&mut parent, child, &ea_row_map);

        // Verify buffer is returned for recycling
        assert_eq!(recycled.nrows(), 3);
        assert_eq!(recycled.ncols(), 3);

        // Expected scatter (lower triangle of parent):
        //
        // child(0,0)=1.0: li=2, lj=2 → parent(2,2) += 1.0
        // child(1,0)=2.0: li=0, lj=2 → li < lj, swap → parent(2,0) += 2.0
        // child(1,1)=3.0: li=0, lj=0 → parent(0,0) += 3.0
        // child(2,0)=4.0: li=3, lj=2 → parent(3,2) += 4.0
        // child(2,1)=5.0: li=3, lj=0 → parent(3,0) += 5.0
        // child(2,2)=6.0: li=3, lj=3 → parent(3,3) += 6.0

        let p = parent.data.as_ref();
        assert_eq!(p[(0, 0)], 3.0, "child(1,1)=3.0 → parent(0,0)");
        assert_eq!(p[(2, 0)], 2.0, "child(1,0)=2.0 swapped → parent(2,0)");
        assert_eq!(p[(2, 2)], 1.0, "child(0,0)=1.0 → parent(2,2)");
        assert_eq!(p[(3, 0)], 5.0, "child(2,1)=5.0 → parent(3,0)");
        assert_eq!(p[(3, 2)], 4.0, "child(2,0)=4.0 → parent(3,2)");
        assert_eq!(p[(3, 3)], 6.0, "child(2,2)=6.0 → parent(3,3)");

        // Verify no spurious entries (untouched positions remain zero)
        assert_eq!(p[(1, 0)], 0.0, "row 1 untouched");
        assert_eq!(p[(1, 1)], 0.0, "row 1 untouched");
        assert_eq!(p[(4, 0)], 0.0, "row 4 untouched");
        assert_eq!(p[(4, 4)], 0.0, "row 4 untouched");
    }

    /// Test column-oriented `extend_add_mapped` with a larger 10×10
    /// contribution (monotonic map), comparing against a naive reference.
    #[test]
    fn test_extend_add_mapped_10x10_vs_reference() {
        let cb_size = 10;
        let parent_size = 15;

        // Monotonically increasing map: child rows map to
        // parent rows [1, 3, 4, 5, 7, 8, 10, 11, 12, 14].
        let ea_row_map: Vec<u32> = vec![1, 3, 4, 5, 7, 8, 10, 11, 12, 14];

        // Build a child contribution with known lower-triangle values.
        // Use value = (i+1)*100 + (j+1) for easy identification.
        let mut child_data = Mat::<f64>::zeros(cb_size, cb_size);
        for j in 0..cb_size {
            for i in j..cb_size {
                child_data[(i, j)] = (i as f64 + 1.0) * 100.0 + (j as f64 + 1.0);
            }
        }

        // Reference: naive element-by-element scatter with lower-triangle swap
        let mut ref_parent = Mat::<f64>::zeros(parent_size, parent_size);
        for j in 0..cb_size {
            let lj = ea_row_map[j] as usize;
            for i in j..cb_size {
                let li = ea_row_map[i] as usize;
                if li >= lj {
                    ref_parent[(li, lj)] += child_data[(i, j)];
                } else {
                    ref_parent[(lj, li)] += child_data[(i, j)];
                }
            }
        }

        // Actual: column-oriented extend_add_mapped
        let mut parent_data = Mat::<f64>::zeros(parent_size, parent_size);
        let parent_row_indices: Vec<usize> = (0..parent_size).collect();
        let mut parent = FrontalMatrix {
            data: parent_data.as_mut(),
            row_indices: &parent_row_indices,
            num_fully_summed: 5,
        };
        let child = ContributionBlock {
            data: child_data,
            row_indices: (0..cb_size).collect(),
            num_delayed: 0,
        };
        let _ = extend_add_mapped(&mut parent, child, &ea_row_map);

        // Compare every element
        for j in 0..parent_size {
            for i in 0..parent_size {
                assert_eq!(
                    parent.data[(i, j)],
                    ref_parent[(i, j)],
                    "mismatch at ({}, {})",
                    i,
                    j
                );
            }
        }
    }

    /// Test column-oriented `extend_add_mapped` with a non-monotonic map
    /// (simulating post-amalgamation parent with non-contiguous owned ranges).
    #[test]
    fn test_extend_add_mapped_non_monotonic() {
        let cb_size = 4;
        let parent_size = 8;

        // Non-monotonic map: simulates child pattern rows mapping into a
        // parent whose owned_ranges are interleaved with pattern rows.
        let ea_row_map: Vec<u32> = vec![5, 2, 7, 1];

        // Build child contribution with lower-triangle values
        let mut child_data = Mat::<f64>::zeros(cb_size, cb_size);
        for j in 0..cb_size {
            for i in j..cb_size {
                child_data[(i, j)] = (i * 10 + j + 1) as f64;
            }
        }

        // Reference: naive scatter with swap
        let mut ref_parent = Mat::<f64>::zeros(parent_size, parent_size);
        for j in 0..cb_size {
            let lj = ea_row_map[j] as usize;
            for i in j..cb_size {
                let li = ea_row_map[i] as usize;
                if li >= lj {
                    ref_parent[(li, lj)] += child_data[(i, j)];
                } else {
                    ref_parent[(lj, li)] += child_data[(i, j)];
                }
            }
        }

        // Actual
        let mut parent_data = Mat::<f64>::zeros(parent_size, parent_size);
        let parent_row_indices: Vec<usize> = (0..parent_size).collect();
        let mut parent = FrontalMatrix {
            data: parent_data.as_mut(),
            row_indices: &parent_row_indices,
            num_fully_summed: 4,
        };
        let child = ContributionBlock {
            data: child_data,
            row_indices: (0..cb_size).collect(),
            num_delayed: 0,
        };
        let _ = extend_add_mapped(&mut parent, child, &ea_row_map);

        // Compare every element
        for j in 0..parent_size {
            for i in 0..parent_size {
                assert_eq!(
                    parent.data[(i, j)],
                    ref_parent[(i, j)],
                    "mismatch at ({}, {})",
                    i,
                    j
                );
            }
        }
    }

    /// Test that `extend_add_mapped` accumulates contributions from multiple
    /// children into the same parent.
    #[test]
    fn test_extend_add_mapped_accumulates() {
        let mut parent_data = Mat::<f64>::zeros(4, 4);
        let parent_row_indices = vec![0, 1, 2, 3];
        let mut parent = FrontalMatrix {
            data: parent_data.as_mut(),
            row_indices: &parent_row_indices,
            num_fully_summed: 2,
        };

        // First child: 2×2, maps to parent rows [1, 3]
        let mut c1_data = Mat::<f64>::zeros(2, 2);
        c1_data[(0, 0)] = 10.0;
        c1_data[(1, 0)] = 20.0;
        c1_data[(1, 1)] = 30.0;
        let child1 = ContributionBlock {
            data: c1_data,
            row_indices: vec![1, 3],
            num_delayed: 0,
        };
        let map1: Vec<u32> = vec![1, 3];

        let _ = extend_add_mapped(&mut parent, child1, &map1);

        // Second child: 2×2, also maps to parent rows [1, 3]
        let mut c2_data = Mat::<f64>::zeros(2, 2);
        c2_data[(0, 0)] = 5.0;
        c2_data[(1, 0)] = 7.0;
        c2_data[(1, 1)] = 9.0;
        let child2 = ContributionBlock {
            data: c2_data,
            row_indices: vec![1, 3],
            num_delayed: 0,
        };
        let map2: Vec<u32> = vec![1, 3];

        let _ = extend_add_mapped(&mut parent, child2, &map2);

        let p = parent.data.as_ref();
        assert_eq!(p[(1, 1)], 15.0, "10 + 5 accumulated");
        assert_eq!(p[(3, 1)], 27.0, "20 + 7 accumulated");
        assert_eq!(p[(3, 3)], 39.0, "30 + 9 accumulated");
    }

    // -----------------------------------------------------------------------
    // Small-leaf subtree classification tests
    // -----------------------------------------------------------------------

    /// Helper: create a SupernodeInfo with given parameters for classification tests.
    /// `ncols` = number of owned columns, `pattern_len` = number of off-diagonal rows.
    /// front_size = ncols + pattern_len.
    #[allow(clippy::single_range_in_vec_init)]
    fn make_snode(
        col_begin: usize,
        ncols: usize,
        pattern_len: usize,
        parent: Option<usize>,
    ) -> SupernodeInfo {
        SupernodeInfo {
            col_begin,
            col_end: col_begin + ncols,
            pattern: (0..pattern_len).map(|i| 1000 + col_begin + i).collect(),
            parent,
            owned_ranges: vec![col_begin..col_begin + ncols],
            in_small_leaf: false,
        }
    }

    /// Linear chain of 5 small supernodes — all should be in_small_leaf.
    #[test]
    fn test_classify_all_small() {
        // Chain: 0 → 1 → 2 → 3 → 4 (each has 2 owned cols, 3 pattern → front_size=5)
        let mut supernodes = vec![
            make_snode(0, 2, 3, Some(1)),
            make_snode(2, 2, 3, Some(2)),
            make_snode(4, 2, 3, Some(3)),
            make_snode(6, 2, 3, Some(4)),
            make_snode(8, 2, 3, None),
        ];
        let children = build_children_map(&supernodes);

        let subtrees = classify_small_leaf_subtrees(&mut supernodes, &children, 256);

        // All 5 nodes should be in_small_leaf
        for sn in &supernodes {
            assert!(
                sn.in_small_leaf,
                "snode {} should be in_small_leaf",
                sn.col_begin
            );
        }
        // One subtree with all 5 nodes in postorder
        assert_eq!(subtrees.len(), 1);
        assert_eq!(subtrees[0].nodes.len(), 5);
        assert_eq!(subtrees[0].nodes, vec![0, 1, 2, 3, 4]);
        assert_eq!(subtrees[0].root, 4);
        assert_eq!(subtrees[0].parent_of_root, None);
    }

    /// Mixed tree — small leaves, large root.
    #[test]
    fn test_classify_mixed_tree() {
        // Tree: leaves 0,1 → parent 2 (large)
        // Nodes 0,1: front_size=5 (small)
        // Node 2: front_size=300 (large, >= 256)
        let mut supernodes = vec![
            make_snode(0, 2, 3, Some(2)),  // front=5, small
            make_snode(2, 2, 3, Some(2)),  // front=5, small
            make_snode(4, 100, 200, None), // front=300, large
        ];
        let children = build_children_map(&supernodes);

        let subtrees = classify_small_leaf_subtrees(&mut supernodes, &children, 256);

        // Node 2 is large → not in_small_leaf
        assert!(!supernodes[2].in_small_leaf);
        // Nodes 0,1 are small leaves but they don't form a subtree with ≥2 nodes
        // that has a common root — they are separate single-leaf subtrees
        // Each leaf alone is 1 node which is < 2, so no subtrees
        assert!(!supernodes[0].in_small_leaf, "single node excluded");
        assert!(!supernodes[1].in_small_leaf, "single node excluded");
        assert_eq!(subtrees.len(), 0);
    }

    /// Single isolated small leaves — no subtree formed (min 2 nodes).
    #[test]
    fn test_classify_single_node_excluded() {
        // Three independent small leaves (no parent linkage between them)
        let mut supernodes = vec![
            make_snode(0, 2, 3, None), // root, small, standalone
            make_snode(2, 2, 3, None), // root, small, standalone
            make_snode(4, 2, 3, None), // root, small, standalone
        ];
        let children = build_children_map(&supernodes);

        let subtrees = classify_small_leaf_subtrees(&mut supernodes, &children, 256);

        // Each is a single-node subtree → excluded
        assert_eq!(subtrees.len(), 0);
        // in_small_leaf should be false for single-node subtrees (cleared by filter)
        for sn in &supernodes {
            assert!(!sn.in_small_leaf);
        }
    }

    /// Threshold boundary — front_size == 256 excluded, 255 included.
    #[test]
    fn test_classify_threshold_boundary() {
        // Node 0: front_size=255 (qualifies), child of node 1
        // Node 1: front_size=256 (does NOT qualify — strict less-than)
        let mut supernodes = vec![
            make_snode(0, 100, 155, Some(1)), // front=255
            make_snode(100, 128, 128, None),  // front=256
        ];
        let children = build_children_map(&supernodes);

        let subtrees = classify_small_leaf_subtrees(&mut supernodes, &children, 256);

        assert!(
            !supernodes[1].in_small_leaf,
            "front_size=256 must NOT qualify"
        );
        // Node 0 is a small leaf but its parent is NOT small → single-node subtree → excluded
        assert!(!supernodes[0].in_small_leaf, "single node excluded");
        assert_eq!(subtrees.len(), 0);

        // Now with front_size=255 for the parent too
        let mut supernodes2 = vec![
            make_snode(0, 100, 155, Some(1)), // front=255
            make_snode(100, 127, 128, None),  // front=255
        ];
        let children2 = build_children_map(&supernodes2);
        let subtrees2 = classify_small_leaf_subtrees(&mut supernodes2, &children2, 256);

        assert!(supernodes2[0].in_small_leaf);
        assert!(supernodes2[1].in_small_leaf);
        assert_eq!(subtrees2.len(), 1);
        assert_eq!(subtrees2[0].nodes, vec![0, 1]);
    }

    /// Disabled when threshold=0.
    #[test]
    fn test_classify_disabled() {
        let mut supernodes = vec![make_snode(0, 2, 3, Some(1)), make_snode(2, 2, 3, None)];
        let children = build_children_map(&supernodes);

        let subtrees = classify_small_leaf_subtrees(&mut supernodes, &children, 0);

        assert_eq!(subtrees.len(), 0);
        assert!(!supernodes[0].in_small_leaf);
        assert!(!supernodes[1].in_small_leaf);
    }

    /// Two independent small-leaf subtrees.
    #[test]
    fn test_classify_multiple_subtrees() {
        // Tree:
        //   0 → 2 (subtree A)
        //   1 → 2 (subtree A)
        //   2 → 5 (large root)
        //   3 → 4 (subtree B)
        //   4 → 5 (subtree B is just 3→4)
        //   5 = large root
        let mut supernodes = vec![
            make_snode(0, 2, 3, Some(2)),   // 0: front=5, small
            make_snode(2, 2, 3, Some(2)),   // 1: front=5, small
            make_snode(4, 2, 3, Some(5)),   // 2: front=5, small (parent is large)
            make_snode(6, 2, 3, Some(4)),   // 3: front=5, small
            make_snode(8, 2, 3, Some(5)),   // 4: front=5, small (parent is large)
            make_snode(10, 100, 200, None), // 5: front=300, large
        ];
        let children = build_children_map(&supernodes);

        let subtrees = classify_small_leaf_subtrees(&mut supernodes, &children, 256);

        // Two subtrees: {0,1,2} and {3,4}
        assert_eq!(subtrees.len(), 2);

        // Subtree A: root=2, nodes=[0,1,2] in postorder
        assert_eq!(subtrees[0].root, 2);
        assert_eq!(subtrees[0].nodes, vec![0, 1, 2]);
        assert_eq!(subtrees[0].parent_of_root, Some(5));

        // Subtree B: root=4, nodes=[3,4] in postorder
        assert_eq!(subtrees[1].root, 4);
        assert_eq!(subtrees[1].nodes, vec![3, 4]);
        assert_eq!(subtrees[1].parent_of_root, Some(5));

        // Large root is not in_small_leaf
        assert!(!supernodes[5].in_small_leaf);
        // All others are
        for (s, sn) in supernodes.iter().enumerate().take(5) {
            assert!(sn.in_small_leaf, "snode {} should be in_small_leaf", s);
        }
    }
}
