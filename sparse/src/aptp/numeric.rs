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

use faer::sparse::SparseColMat;
use faer::sparse::linalg::cholesky::SymbolicCholeskyRaw;
use faer::{Mat, MatMut, Par};

use super::diagonal::{MixedDiagonal, PivotEntry};
use super::factor::{AptpFactorResult, AptpOptions, aptp_factor_in_place};
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
            };
        }
        Self {
            frontal_data: Mat::zeros(max_front, max_front),
            frontal_row_indices: Vec::with_capacity(max_front),
            delayed_cols_buf: Vec::with_capacity(max_front),
            global_to_local: vec![NOT_IN_FRONT; n],
        }
    }

    /// Empty workspace for `Cell` default in thread-local storage.
    const fn empty() -> Self {
        Self {
            frontal_data: Mat::new(),
            frontal_row_indices: Vec::new(),
            delayed_cols_buf: Vec::new(),
            global_to_local: Vec::new(),
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
        // Zero the lower triangle of the m × m subregion (column-wise)
        let mut sub = self.frontal_data.as_mut().submatrix_mut(0, 0, m, m);
        for j in 0..m {
            for i in j..m {
                sub[(i, j)] = 0.0;
            }
        }
    }

    /// Ensure workspace capacity is at least `max_front × max_front` for frontal
    /// data and `n` for global-to-local mapping. Grows buffers if needed.
    fn ensure_capacity(&mut self, max_front: usize, n: usize) {
        if self.frontal_data.nrows() < max_front || self.frontal_data.ncols() < max_front {
            self.frontal_data = Mat::zeros(max_front, max_front);
        }
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
        let supernodes = super::amalgamation::amalgamate(supernodes, nemin);
        let n_supernodes = supernodes.len();
        let merges_performed = n_supernodes_before - n_supernodes;
        let children = build_children_map(&supernodes);

        // Get fill-reducing permutation
        let (perm_fwd, perm_inv) = symbolic.perm_vecs();

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
                        extend_add(&mut frontal, &cb, &global_to_local);
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
                let frontal_view = FrontalMatrix {
                    data: frontal_data.as_mut(),
                    row_indices: &frontal_rows,
                    num_fully_summed: k,
                };
                contributions[s] = Some(extract_contribution(&frontal_view, &result));
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
    child_contributions: Vec<ContributionBlock>,
    matrix: &SparseColMat<usize, f64>,
    perm_fwd: &[usize],
    perm_inv: &[usize],
    options: &AptpOptions,
    scaling: Option<&[f64]>,
    workspace: &mut FactorizationWorkspace,
) -> Result<SupernodeResult, SparseError> {
    // 1. Collect delayed columns from children into workspace buffer
    workspace.delayed_cols_buf.clear();
    for cb in &child_contributions {
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
    workspace.zero_frontal(m);

    // 4. Build global-to-local mapping
    for (i, &global) in workspace.frontal_row_indices.iter().enumerate() {
        workspace.global_to_local[global] = i;
    }

    // 5. Assemble frontal matrix using workspace buffer
    #[cfg(feature = "diagnostic")]
    let assembly_start = std::time::Instant::now();

    let mut frontal = FrontalMatrix {
        data: workspace.frontal_data.as_mut().submatrix_mut(0, 0, m, m),
        row_indices: &workspace.frontal_row_indices,
        num_fully_summed: k,
    };

    scatter_original_entries_multi(
        &mut frontal,
        matrix,
        perm_fwd,
        perm_inv,
        &workspace.global_to_local,
        &sn.owned_ranges,
        scaling,
    );

    for cb in child_contributions {
        extend_add(&mut frontal, &cb, &workspace.global_to_local);
    }

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

    // 7. Extract front factors (reads from workspace frontal data)
    #[cfg(feature = "diagnostic")]
    let extraction_start = std::time::Instant::now();

    // Build a temporary FrontalMatrix view for extraction functions
    let frontal_view = FrontalMatrix {
        data: workspace.frontal_data.as_mut().submatrix_mut(0, 0, m, m),
        row_indices: &workspace.frontal_row_indices,
        num_fully_summed: k,
    };
    let ff = extract_front_factors(&frontal_view, &result);

    // 8. Prepare contribution for parent (if not root and not fully eliminated)
    let contribution = if sn.parent.is_some() && ne < m {
        Some(extract_contribution(&frontal_view, &result))
    } else {
        None
    };

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
    };

    // Reset global_to_local entries for reuse (O(m), not O(n))
    for &global in &workspace.frontal_row_indices[..m] {
        workspace.global_to_local[global] = NOT_IN_FRONT;
    }

    Ok(SupernodeResult {
        ff,
        contribution,
        stats,
    })
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

    // Pre-allocate workspace for the sequential path. Sized to the estimated
    // maximum front dimension. Reused across all supernodes.
    let max_front = estimate_max_front_size(supernodes);
    let mut shared_workspace = FactorizationWorkspace::new(max_front, n);

    // Initial ready set: all leaf supernodes (no children)
    let mut ready: Vec<usize> = (0..n_supernodes)
        .filter(|&s| remaining_children[s] == 0)
        .collect();

    while !ready.is_empty() {
        // Collect child contributions for each ready node. Children of different
        // ready nodes are disjoint (tree structure), so each contribution is
        // taken exactly once.
        let mut batch_inputs: Vec<(usize, Vec<ContributionBlock>)> =
            Vec::with_capacity(ready.len());
        for &s in &ready {
            let child_contribs: Vec<ContributionBlock> = children[s]
                .iter()
                .filter_map(|&c| contributions[c].take())
                .collect();
            batch_inputs.push((s, child_contribs));
        }

        // Factor all ready nodes (parallel if enabled and batch > 1)
        let batch_results: Vec<SupernodeResult> =
            if matches!(options.par, Par::Seq) || batch_inputs.len() == 1 {
                // Sequential: reuse shared workspace (frontal buffer + g2l)
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
                            &mut shared_workspace,
                        )
                    })
                    .collect::<Result<_, _>>()?
            } else {
                // Parallel: use thread-local workspace so each rayon worker
                // allocates once and reuses across all waves/supernodes.
                // We use Cell::take/set (move semantics) instead of RefCell
                // because factor_single_supernode may invoke rayon-parallel
                // BLAS (Par::rayon) which can work-steal other tasks from
                // this par_iter onto the same thread, causing re-entrant
                // borrow panics with RefCell.
                use rayon::iter::{IntoParallelIterator, ParallelIterator};
                use std::cell::Cell;
                thread_local! {
                    static FACTOR_WORKSPACE: Cell<FactorizationWorkspace> =
                        const { Cell::new(FactorizationWorkspace::empty()) };
                }
                batch_inputs
                    .into_par_iter()
                    .map(|(s, contribs)| {
                        // Take the workspace out of thread-local (leaves empty behind).
                        // If another task on this thread runs while we hold the workspace
                        // (rayon work-stealing during nested parallelism), it will get
                        // an empty workspace, resize it, and put it back when done.
                        let mut ws = FACTOR_WORKSPACE.take();
                        ws.ensure_capacity(max_front, n);
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
                        );
                        // Put workspace back for reuse by the next task on this thread.
                        FACTOR_WORKSPACE.set(ws);
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
                        owned_ranges: vec![begin[s]..end[s]],
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
    child: &ContributionBlock,
    global_to_local: &[usize],
) {
    let cb_size = child.data.nrows();
    for i in 0..cb_size {
        let gi = child.row_indices[i];
        let li = global_to_local[gi];
        debug_assert!(
            li != NOT_IN_FRONT,
            "extend_add: child row {} not in parent",
            gi
        );
        for j in 0..=i {
            // Lower triangle only
            let gj = child.row_indices[j];
            let lj = global_to_local[gj];
            debug_assert!(
                lj != NOT_IN_FRONT,
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
    frontal: &FrontalMatrix<'_>,
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
                unreachable!("unexpected Delayed pivot at col {} in 0..ne", col);
            }
        }
    }

    // Extract L21 (r × ne) where r = m - ne (includes delayed rows ne..k AND
    // non-fully-summed rows k..m). The factored frontal has valid L entries at rows
    // ne..m for columns 0..ne: rows ne..k are delayed fully-summed rows that received
    // TRSM updates before being delayed; rows k..m are non-fully-summed rows.
    //
    // NOTE: extract_contribution also includes delayed rows (see lines below).
    // These two functions MUST be consistent in their row treatment.
    let r = m - ne;
    let l21 = if ne > 0 && r > 0 {
        let mut l = Mat::zeros(r, ne);
        for i in 0..r {
            for j in 0..ne {
                l[(i, j)] = frontal.data[(ne + i, j)];
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

    // Row indices: global permuted indices for L21 rows.
    // First num_delayed entries: delayed fully-summed columns (ne..k), mapped through
    // result.perm and frontal.row_indices to get global permuted indices.
    // Remaining entries: non-fully-summed rows (k..m).
    // This matches extract_contribution's row_indices construction.
    let mut row_indices = Vec::with_capacity(m - ne);
    for &lp in &result.perm[ne..k] {
        row_indices.push(frontal.row_indices[lp]);
    }
    row_indices.extend_from_slice(&frontal.row_indices[k..]);

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
    frontal: &FrontalMatrix<'_>,
    result: &AptpFactorResult,
) -> ContributionBlock {
    let m = frontal.data.nrows();
    let ne = result.num_eliminated;
    let k = frontal.num_fully_summed;
    let size = m - ne;
    let num_delayed = k - ne;

    // Copy trailing submatrix (lower triangle only) via column-wise bulk copy.
    // Each column j copies rows j..size from the frontal submatrix.
    let mut data = Mat::zeros(size, size);
    for j in 0..size {
        let col_len = size - j;
        for i in 0..col_len {
            data[(j + i, j)] = frontal.data[(ne + j + i, ne + j)];
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

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aptp::factor::AptpOptions;

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

        let frontal = FrontalMatrix {
            data: data.as_mut(),
            row_indices: &row_indices,
            num_fully_summed: k,
        };

        // The bug: if delays occurred, the old code would make L21 have (m - k) rows
        // instead of (m - ne) rows, missing the delayed-row entries.
        if num_delayed > 0 {
            let ff = extract_front_factors(&frontal, &result);

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
            let ff = extract_front_factors(&frontal, &result);
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
            let frontal = FrontalMatrix {
                data: data.as_mut(),
                row_indices: &row_indices,
                num_fully_summed: k,
            };
            let ff = extract_front_factors(&frontal, &result);
            let cb = extract_contribution(&frontal, &result);

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
            let frontal = FrontalMatrix {
                data: data.as_mut(),
                row_indices: &row_indices,
                num_fully_summed: k,
            };
            let ff = extract_front_factors(&frontal, &result);

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

        let frontal = FrontalMatrix {
            data: data.as_mut(),
            row_indices: &row_indices,
            num_fully_summed: k,
        };
        let ff = extract_front_factors(&frontal, &result);

        // Compute L11 * D * L11^T manually
        // Step 1: LD = L11 * D (column-by-column, handling 1x1 and 2x2 pivots)
        let mut ld = Mat::zeros(ne, ne);
        let mut col = 0;
        while col < ne {
            match ff.d11.get_pivot_type(col) {
                PivotType::OneByOne => {
                    let d = ff.d11.get_1x1(col);
                    for i in 0..ne {
                        ld[(i, col)] = ff.l11[(i, col)] * d;
                    }
                    col += 1;
                }
                PivotType::TwoByTwo { .. } => {
                    let block = ff.d11.get_2x2(col);
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

        let frontal = FrontalMatrix {
            data: data.as_mut(),
            row_indices: &row_indices,
            num_fully_summed: k,
        };
        let ff = extract_front_factors(&frontal, &result);

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
            match ff.d11.get_pivot_type(col) {
                PivotType::OneByOne => {
                    let d = ff.d11.get_1x1(col);
                    for i in 0..full_l_rows {
                        ld[(i, col)] = full_l[(i, col)] * d;
                    }
                    col += 1;
                }
                PivotType::TwoByTwo { .. } => {
                    let block = ff.d11.get_2x2(col);
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
}
