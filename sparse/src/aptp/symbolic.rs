//! Symbolic analysis for APTP (A Posteriori Threshold Pivoting) factorization.
//!
//! This module provides [`AptpSymbolic`], the symbolic analysis result for the
//! APTP-based sparse symmetric indefinite LDL^T solver. It composes faer's
//! [`SymbolicCholesky`] with APTP-specific metadata: elimination tree parent
//! pointers, per-column nonzero counts, and heuristic pivot buffer estimates.
//!
//! The symbolic analysis is the first step of the three-phase
//! analyze → factorize → solve pipeline. The result is reusable across multiple
//! numeric factorizations sharing the same sparsity pattern.
//!
//! # Algorithm References
//!
//! - Liu (1990), "The role of elimination trees in sparse factorization" —
//!   foundational theory for elimination trees
//! - Gilbert, Ng & Peyton (1992), "An efficient algorithm to compute row and
//!   column counts for sparse Cholesky factorization" — column count prediction
//! - Gilbert, Ng & Peyton (1994), "An efficient algorithm to compute row and
//!   column counts for sparse Cholesky factorization" (journal version)
//! - Hogg, Scott & Sherwood-Jones (2016), "A New Sparse LDLT Solver using
//!   A Posteriori Threshold Pivoting" — SSIDS analyze phase, delayed pivot
//!   buffer estimation (Section 2.2–2.4)

use std::fmt;

use faer::Side;
use faer::dyn_stack::MemStack;
use faer::sparse::linalg::cholesky::simplicial::{
    prefactorize_symbolic_cholesky, prefactorize_symbolic_cholesky_scratch,
};
use faer::sparse::linalg::cholesky::{
    CholeskySymbolicParams, SymbolicCholesky, SymbolicCholeskyRaw, SymmetricOrdering,
    factorize_symbolic_cholesky,
};
use faer::sparse::{SparseColMat, SymbolicSparseColMatRef, Triplet};

use crate::error::SparseError;

/// Build the permuted upper-triangular symbolic structure P^T A P.
///
/// For each column `j_new` in the permuted matrix, the original column is `perm_fwd[j_new]`.
/// Row `i` in the original becomes `perm_inv[i]` in the permuted matrix.
/// Only entries where `permuted_row <= j_new` (upper triangle) are kept.
///
/// Returns a `SparseColMat<usize, f64>` with dummy 1.0 values — only the symbolic
/// structure matters for downstream use in `prefactorize_symbolic_cholesky`.
///
/// # Arguments
///
/// * `matrix` — Symbolic CSC sparsity pattern of the original matrix.
/// * `perm_fwd` — Forward permutation: `perm_fwd[new_index] = old_index`.
/// * `perm_inv` — Inverse permutation: `perm_inv[old_index] = new_index`.
///
/// Reference: Gilbert et al. (1992), Section 3.3 for permuted structure computation.
pub(crate) fn permute_symbolic_upper_triangle(
    matrix: SymbolicSparseColMatRef<'_, usize>,
    perm_fwd: &[usize],
    perm_inv: &[usize],
) -> Result<SparseColMat<usize, f64>, SparseError> {
    let n = matrix.nrows();
    let col_ptrs = matrix.col_ptr();
    let row_indices = matrix.row_idx();

    // Count nonzeros per column in the permuted upper triangle
    let mut permuted_col_nnz = vec![0usize; n];
    for j_new in 0..n {
        let j_orig = perm_fwd[j_new];
        let start = col_ptrs[j_orig];
        let end = col_ptrs[j_orig + 1];
        for &i_orig in &row_indices[start..end] {
            let i_new = perm_inv[i_orig];
            if i_new <= j_new {
                // Upper triangle: row <= col
                permuted_col_nnz[j_new] += 1;
            }
        }
    }

    // Build CSC column pointers
    let mut permuted_col_ptrs = vec![0usize; n + 1];
    for j in 0..n {
        permuted_col_ptrs[j + 1] = permuted_col_ptrs[j] + permuted_col_nnz[j];
    }
    let total_nnz = permuted_col_ptrs[n];

    // Fill row indices (sorted within each column)
    let mut permuted_row_indices = vec![0usize; total_nnz];
    let mut write_pos = permuted_col_ptrs.clone();
    for j_new in 0..n {
        let j_orig = perm_fwd[j_new];
        let start = col_ptrs[j_orig];
        let end = col_ptrs[j_orig + 1];
        for &i_orig in &row_indices[start..end] {
            let i_new = perm_inv[i_orig];
            if i_new <= j_new {
                permuted_row_indices[write_pos[j_new]] = i_new;
                write_pos[j_new] += 1;
            }
        }
        // Sort row indices within this column for CSC validity
        let col_start = permuted_col_ptrs[j_new];
        let col_end = permuted_col_ptrs[j_new + 1];
        permuted_row_indices[col_start..col_end].sort_unstable();
    }

    // Build a SparseColMat from the permuted structure (using dummy f64 values).
    let mut triplets = Vec::with_capacity(total_nnz);
    for j in 0..n {
        let start = permuted_col_ptrs[j];
        let end = permuted_col_ptrs[j + 1];
        for &row in &permuted_row_indices[start..end] {
            triplets.push(Triplet::new(row, j, 1.0f64));
        }
    }
    SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets).map_err(|e| {
        SparseError::AnalysisFailure {
            reason: format!("Failed to construct permuted symbolic structure: {}", e),
        }
    })
}

/// Symbolic analysis result for APTP factorization.
///
/// Wraps faer's [`SymbolicCholesky`] and supplements it with APTP-specific
/// metadata that faer does not expose: the elimination tree parent pointers,
/// per-column nonzero counts, and heuristic delayed-pivot buffer estimates.
///
/// Created via [`AptpSymbolic::analyze`], which is deterministic: the same
/// matrix and ordering always produce identical results.
///
/// # Examples
///
/// ```
/// use faer::sparse::{SparseColMat, Triplet};
/// use faer::sparse::linalg::cholesky::SymmetricOrdering;
/// use rivrs_sparse::aptp::AptpSymbolic;
///
/// // Build a small symmetric matrix
/// let triplets = vec![
///     Triplet::new(0, 0, 4.0),
///     Triplet::new(1, 1, 4.0),
///     Triplet::new(0, 1, 1.0),
///     Triplet::new(1, 0, 1.0),
/// ];
/// let matrix = SparseColMat::try_new_from_triplets(2, 2, &triplets).unwrap();
///
/// let symbolic = AptpSymbolic::analyze(
///     matrix.symbolic(),
///     SymmetricOrdering::Amd,
/// ).expect("symbolic analysis failed");
///
/// println!("{}", symbolic.statistics());
/// ```
pub struct AptpSymbolic {
    /// faer's symbolic factorization result (permutation, L structure, simplicial/supernodal).
    inner: SymbolicCholesky<usize>,

    /// Elimination tree parent pointers for the permuted matrix.
    /// `etree[j]` is the parent of column `j` in the elimination tree, or `-1` if `j` is a root.
    /// Computed via `prefactorize_symbolic_cholesky` on the permuted structure.
    etree: Vec<isize>,

    /// Predicted nonzero count per column of L (from `prefactorize_symbolic_cholesky`).
    col_counts: Vec<usize>,

    /// Per-supernode (supernodal mode) or per-column (simplicial mode) extra space
    /// estimates for APTP delayed pivots.
    pivot_buffer: Vec<usize>,
}

/// Diagnostic summary of a symbolic analysis.
///
/// Provides human-readable statistics about the symbolic factorization result,
/// useful for diagnostics, benchmarking, and debugging.
#[derive(Debug, Clone)]
pub struct SymbolicStatistics {
    /// Matrix dimension (n for n×n).
    pub dimension: usize,
    /// Predicted number of nonzeros in the L factor.
    pub predicted_nnz: usize,
    /// Mean nonzeros per column in L (`predicted_nnz / dimension`).
    pub average_col_count: f64,
    /// Whether faer chose supernodal over simplicial.
    pub is_supernodal: bool,
    /// Number of supernodes (`Some` if supernodal, `None` if simplicial).
    pub n_supernodes: Option<usize>,
    /// Sum of all pivot buffer estimates.
    pub total_pivot_buffer: usize,
}

impl fmt::Display for SymbolicStatistics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Symbolic Analysis Statistics:")?;
        writeln!(f, "  Dimension:          {}", self.dimension)?;
        writeln!(f, "  Predicted NNZ:      {}", self.predicted_nnz)?;
        writeln!(f, "  Avg col count:      {:.2}", self.average_col_count)?;
        writeln!(
            f,
            "  Mode:               {}",
            if self.is_supernodal {
                "supernodal"
            } else {
                "simplicial"
            }
        )?;
        if let Some(ns) = self.n_supernodes {
            writeln!(f, "  Supernodes:         {}", ns)?;
        }
        write!(f, "  Total pivot buffer: {}", self.total_pivot_buffer)
    }
}

impl AptpSymbolic {
    /// Perform symbolic analysis on a sparse symmetric matrix.
    ///
    /// Computes the fill-reducing ordering, elimination tree, predicted factor
    /// structure, and APTP-specific pivot buffer estimates.
    ///
    /// # Arguments
    ///
    /// * `matrix` — Symbolic sparsity pattern of the sparse symmetric matrix (CSC format).
    ///   Only the upper triangle is used.
    /// * `ordering` — Fill-reducing ordering strategy. Use `SymmetricOrdering::Amd` for
    ///   the default approximate minimum degree ordering.
    ///
    /// # Errors
    ///
    /// * [`SparseError::NotSquare`] — if the matrix is not square
    /// * [`SparseError::DimensionMismatch`] — if a custom permutation has wrong dimension
    /// * [`SparseError::AnalysisFailure`] — if the internal symbolic factorization fails
    ///
    /// # Examples
    ///
    /// ```
    /// use faer::sparse::{SparseColMat, Triplet};
    /// use faer::sparse::linalg::cholesky::SymmetricOrdering;
    /// use rivrs_sparse::aptp::AptpSymbolic;
    ///
    /// let triplets = vec![
    ///     Triplet::new(0, 0, 1.0),
    ///     Triplet::new(1, 1, 1.0),
    /// ];
    /// let matrix = SparseColMat::try_new_from_triplets(2, 2, &triplets).unwrap();
    /// let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)?;
    /// assert_eq!(symbolic.nrows(), 2);
    /// # Ok::<(), rivrs_sparse::error::SparseError>(())
    /// ```
    pub fn analyze(
        matrix: SymbolicSparseColMatRef<'_, usize>,
        ordering: SymmetricOrdering<'_, usize>,
    ) -> Result<Self, SparseError> {
        let nrows = matrix.nrows();
        let ncols = matrix.ncols();

        // Validate: must be square
        if nrows != ncols {
            return Err(SparseError::NotSquare {
                dims: (nrows, ncols),
            });
        }

        // Validate custom ordering dimension
        if let SymmetricOrdering::Custom(perm) = &ordering {
            let perm_len = perm.len();
            if perm_len != nrows {
                return Err(SparseError::DimensionMismatch {
                    expected: (nrows, nrows),
                    got: (perm_len, perm_len),
                    context: format!(
                        "Custom permutation has dimension {} but matrix is {}×{}",
                        perm_len, nrows, ncols
                    ),
                });
            }
        }

        let n = nrows;

        // Step 1: Call factorize_symbolic_cholesky for the full symbolic result
        let inner = factorize_symbolic_cholesky(
            matrix,
            Side::Upper,
            ordering,
            CholeskySymbolicParams::default(),
        )
        .map_err(|e| SparseError::AnalysisFailure {
            reason: format!("faer symbolic Cholesky failed: {}", e),
        })?;

        // Step 2: Compute etree and col_counts on the permuted structure
        // We need to apply the permutation to get the correct permuted etree.
        let (etree, col_counts) = if n == 0 {
            (Vec::new(), Vec::new())
        } else {
            Self::compute_permuted_etree_and_col_counts(matrix, &inner)?
        };

        // Step 3: Compute pivot buffer estimates using the 10% heuristic
        // Hogg et al. (2016), Section 2.4: delayed pivots propagate entries
        // up the assembly tree. We estimate 10% of the column count per unit
        // (supernode or column) as extra buffer space.
        let pivot_buffer = Self::compute_pivot_buffer(&inner, &col_counts);

        Ok(AptpSymbolic {
            inner,
            etree,
            col_counts,
            pivot_buffer,
        })
    }

    /// Compute the elimination tree and column counts on the permuted matrix structure.
    ///
    /// Extracts the permutation from the inner `SymbolicCholesky`, delegates to
    /// [`permute_symbolic_upper_triangle`] to build P^T A P, then calls
    /// `prefactorize_symbolic_cholesky` to obtain the etree parent pointers and
    /// column counts.
    ///
    /// # Optimization notes
    ///
    /// The triplet round-trip (CSC → triplets → SparseColMat → internal CSC) allocates
    /// the permuted structure data ~3 times. When the permutation is identity, we could
    /// short-circuit and call `prefactorize_symbolic_cholesky` directly on the original
    /// matrix. Currently acceptable: O(nnz · α(n)) `prefactorize` cost is negligible
    /// compared to the full symbolic factorization.
    fn compute_permuted_etree_and_col_counts(
        matrix: SymbolicSparseColMatRef<'_, usize>,
        inner: &SymbolicCholesky<usize>,
    ) -> Result<(Vec<isize>, Vec<usize>), SparseError> {
        let n = matrix.nrows();

        // Extract permutation arrays (or identity if no permutation).
        let (perm_fwd, perm_inv) = if let Some(perm) = inner.perm() {
            let (fwd, inv) = perm.arrays();
            (fwd.to_vec(), inv.to_vec())
        } else {
            let id: Vec<usize> = (0..n).collect();
            (id.clone(), id)
        };

        // Build the permuted upper-triangular symbolic structure P^T A P.
        let permuted_symbolic = permute_symbolic_upper_triangle(matrix, &perm_fwd, &perm_inv)?;

        // Call prefactorize_symbolic_cholesky on the permuted structure
        let mut etree_buf = vec![0isize; n];
        let mut col_counts_buf = vec![0usize; n];

        let scratch_req =
            prefactorize_symbolic_cholesky_scratch::<usize>(n, permuted_symbolic.compute_nnz());
        let mut mem = faer::dyn_stack::MemBuffer::new(scratch_req);
        let stack = MemStack::new(&mut mem);

        let _etree_ref = prefactorize_symbolic_cholesky(
            &mut etree_buf,
            &mut col_counts_buf,
            permuted_symbolic.symbolic(),
            stack,
        );

        Ok((etree_buf, col_counts_buf))
    }

    /// Compute per-unit pivot buffer estimates using the 10% heuristic.
    ///
    /// For supernodal mode: one estimate per supernode, based on the total
    /// column count across the supernode's columns.
    ///
    /// For simplicial mode: one estimate per column, based on each column's
    /// predicted nonzero count.
    ///
    /// The 10% fraction is a conservative starting point; it will be validated
    /// empirically during numeric factorization.
    ///
    /// Reference: Hogg et al. (2016), Section 2.4 — delayed pivot propagation
    /// and workspace estimation.
    fn compute_pivot_buffer(inner: &SymbolicCholesky<usize>, col_counts: &[usize]) -> Vec<usize> {
        const BUFFER_FRACTION: f64 = 0.10;

        match inner.raw() {
            SymbolicCholeskyRaw::Supernodal(sn) => {
                let ns = sn.n_supernodes();
                let begin = sn.supernode_begin();
                let end = sn.supernode_end();
                (0..ns)
                    .map(|s| {
                        let total_col_count: usize = (begin[s]..end[s])
                            .map(|j| col_counts.get(j).copied().unwrap_or(0))
                            .sum();
                        (total_col_count as f64 * BUFFER_FRACTION).ceil() as usize
                    })
                    .collect()
            }
            SymbolicCholeskyRaw::Simplicial(_) => col_counts
                .iter()
                .map(|&c| (c as f64 * BUFFER_FRACTION).ceil() as usize)
                .collect(),
        }
    }

    // ---- faer-delegating accessors ----

    /// Returns the fill-reducing permutation, if one was computed.
    ///
    /// Returns `None` if `SymmetricOrdering::Identity` was used.
    ///
    /// # Examples
    ///
    /// ```
    /// # use faer::sparse::{SparseColMat, Triplet};
    /// # use faer::sparse::linalg::cholesky::SymmetricOrdering;
    /// # use rivrs_sparse::aptp::AptpSymbolic;
    /// # let triplets = vec![Triplet::new(0, 0, 1.0)];
    /// # let matrix = SparseColMat::try_new_from_triplets(1, 1, &triplets).unwrap();
    /// let symbolic = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)?;
    /// if let Some(perm) = symbolic.perm() {
    ///     let (fwd, inv) = perm.arrays();
    ///     // fwd[new_index] = old_index
    /// }
    /// # Ok::<(), rivrs_sparse::error::SparseError>(())
    /// ```
    pub fn perm(&self) -> Option<faer::perm::PermRef<'_, usize>> {
        self.inner.perm()
    }

    /// Returns the forward and inverse permutation arrays as owned vectors.
    ///
    /// If no permutation is applied (identity ordering), returns identity
    /// arrays `(0..n, 0..n)`.
    ///
    /// - `fwd[new_index] = old_index`
    /// - `inv[old_index] = new_index`
    pub fn perm_vecs(&self) -> (Vec<usize>, Vec<usize>) {
        if let Some(perm) = self.inner.perm() {
            let (fwd, inv) = perm.arrays();
            (fwd.to_vec(), inv.to_vec())
        } else {
            let n = self.nrows();
            let id: Vec<usize> = (0..n).collect();
            (id.clone(), id)
        }
    }

    /// Returns the predicted number of nonzeros in the L factor.
    pub fn predicted_nnz(&self) -> usize {
        self.inner.len_val()
    }

    /// Returns the number of rows in the matrix.
    pub fn nrows(&self) -> usize {
        self.inner.nrows()
    }

    /// Returns the number of columns in the matrix.
    pub fn ncols(&self) -> usize {
        self.inner.ncols()
    }

    /// Returns a reference to the inner symbolic structure (simplicial or supernodal).
    ///
    /// This enables downstream code to access faer's full symbolic API, including
    /// workspace estimation for numeric factorization.
    pub fn raw(&self) -> &SymbolicCholeskyRaw<usize> {
        self.inner.raw()
    }

    /// Returns a reference to the inner `SymbolicCholesky`.
    ///
    /// Useful for passing directly to faer's numeric factorization routines.
    pub fn inner(&self) -> &SymbolicCholesky<usize> {
        &self.inner
    }

    // ---- APTP-specific accessors ----

    /// Returns the elimination tree as a parent pointer array.
    ///
    /// `etree[j]` is the parent of column `j` in the elimination tree of the
    /// **permuted** matrix, or `-1` if `j` is a root.
    ///
    /// Reference: Liu (1990), "The role of elimination trees in sparse factorization".
    pub fn etree(&self) -> &[isize] {
        &self.etree
    }

    /// Returns the predicted nonzero count per column of L.
    ///
    /// Reference: Gilbert, Ng & Peyton (1994).
    pub fn col_counts(&self) -> &[usize] {
        &self.col_counts
    }

    /// Returns per-unit delayed pivot buffer estimates.
    ///
    /// For supernodal mode, the length equals `n_supernodes().unwrap()`.
    /// For simplicial mode, the length equals `nrows()`.
    ///
    /// Reference: Hogg et al. (2016), Section 2.4 — delayed pivot propagation.
    pub fn pivot_buffer_estimates(&self) -> &[usize] {
        &self.pivot_buffer
    }

    /// Returns the total pivot buffer estimate (sum of all per-unit estimates).
    pub fn total_pivot_buffer(&self) -> usize {
        self.pivot_buffer.iter().sum()
    }

    /// Returns whether faer chose supernodal over simplicial factorization.
    pub fn is_supernodal(&self) -> bool {
        matches!(self.inner.raw(), SymbolicCholeskyRaw::Supernodal(_))
    }

    /// Returns the number of supernodes, if the result is supernodal.
    ///
    /// Returns `None` for simplicial factorization.
    pub fn n_supernodes(&self) -> Option<usize> {
        match self.inner.raw() {
            SymbolicCholeskyRaw::Supernodal(sn) => Some(sn.n_supernodes()),
            SymbolicCholeskyRaw::Simplicial(_) => None,
        }
    }

    /// Returns the starting column index for each supernode.
    ///
    /// Returns `None` for simplicial factorization.
    pub fn supernode_begin(&self) -> Option<&[usize]> {
        match self.inner.raw() {
            SymbolicCholeskyRaw::Supernodal(sn) => Some(sn.supernode_begin()),
            SymbolicCholeskyRaw::Simplicial(_) => None,
        }
    }

    /// Returns the ending column index (exclusive) for each supernode.
    ///
    /// Returns `None` for simplicial factorization.
    pub fn supernode_end(&self) -> Option<&[usize]> {
        match self.inner.raw() {
            SymbolicCholeskyRaw::Supernodal(sn) => Some(sn.supernode_end()),
            SymbolicCholeskyRaw::Simplicial(_) => None,
        }
    }

    /// Returns the row pattern (off-diagonal row indices) for supernode `s`.
    ///
    /// The pattern contains row indices of off-diagonal entries for the
    /// supernode's columns in the L factor.
    ///
    /// Returns `None` for simplicial factorization.
    ///
    /// # Panics
    ///
    /// Panics if `s >= n_supernodes()` for a supernodal factorization.
    pub fn supernode_pattern(&self, s: usize) -> Option<&[usize]> {
        match self.inner.raw() {
            SymbolicCholeskyRaw::Supernodal(sn) => Some(sn.supernode(s).pattern()),
            SymbolicCholeskyRaw::Simplicial(_) => None,
        }
    }

    /// Returns the parent supernode of supernode `s` in the assembly tree.
    ///
    /// Derived from the column-level elimination tree: for supernode `s` with
    /// column range `[begin, end)`, the parent is determined by `etree[end - 1]`.
    /// Returns `None` if `s` is a root of the assembly tree.
    ///
    /// For simplicial factorization, always returns `None`.
    ///
    /// Reference: Liu (1992), "The Multifrontal Method" — assembly tree derivation
    /// from elimination tree.
    ///
    /// # Panics
    ///
    /// Panics if `s >= n_supernodes()` for a supernodal factorization.
    // Optimization note: The binary search per call is O(log ns).
    // If multifrontal factorization traverses the assembly tree for all supernodes,
    // consider precomputing a parent array in `analyze()` for O(1) per lookup.
    pub fn supernode_parent(&self, s: usize) -> Option<usize> {
        match self.inner.raw() {
            SymbolicCholeskyRaw::Supernodal(sn) => {
                let begin = sn.supernode_begin();
                let end = sn.supernode_end();
                let last_col = end[s] - 1;
                let parent_col = self.etree[last_col];
                if parent_col < 0 {
                    None // root supernode
                } else {
                    // Find which supernode contains the parent column
                    let parent_col = parent_col as usize;
                    // Binary search: find the supernode s' where begin[s'] <= parent_col < end[s']
                    let ns = sn.n_supernodes();
                    let idx = begin[..ns]
                        .partition_point(|&b| b <= parent_col)
                        .saturating_sub(1);
                    Some(idx)
                }
            }
            SymbolicCholeskyRaw::Simplicial(_) => None,
        }
    }

    /// Returns a diagnostic summary of the symbolic analysis.
    ///
    /// The returned [`SymbolicStatistics`] includes dimension, predicted fill-in,
    /// supernodal status, and pivot buffer totals.
    pub fn statistics(&self) -> SymbolicStatistics {
        let dimension = self.nrows();
        let predicted_nnz = self.predicted_nnz();
        let average_col_count = if dimension > 0 {
            predicted_nnz as f64 / dimension as f64
        } else {
            0.0
        };
        let is_supernodal = self.is_supernodal();
        let n_supernodes = self.n_supernodes();
        let total_pivot_buffer = self.total_pivot_buffer();

        SymbolicStatistics {
            dimension,
            predicted_nnz,
            average_col_count,
            is_supernodal,
            n_supernodes,
            total_pivot_buffer,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: create a 1×1 sparse matrix with a single diagonal entry.
    fn make_1x1() -> SparseColMat<usize, f64> {
        let triplets = vec![Triplet::new(0, 0, 5.0)];
        SparseColMat::try_new_from_triplets(1, 1, &triplets).unwrap()
    }

    /// Helper: create an n×n diagonal sparse matrix.
    fn make_diagonal(n: usize) -> SparseColMat<usize, f64> {
        let triplets: Vec<_> = (0..n).map(|i| Triplet::new(i, i, (i + 1) as f64)).collect();
        SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
    }

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

    /// Helper: create a small arrow matrix (dense first row/column + diagonal).
    fn make_arrow(n: usize) -> SparseColMat<usize, f64> {
        let mut triplets = Vec::new();
        // Diagonal
        for i in 0..n {
            triplets.push(Triplet::new(i, i, 10.0));
        }
        // First row/column (off-diagonal)
        for i in 1..n {
            triplets.push(Triplet::new(0, i, 1.0));
            triplets.push(Triplet::new(i, 0, 1.0));
        }
        SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
    }

    /// Helper: create a block-diagonal matrix with two disconnected blocks.
    fn make_block_diagonal() -> SparseColMat<usize, f64> {
        // 4×4 with two 2×2 blocks
        let triplets = vec![
            // Block 1: rows/cols 0-1
            Triplet::new(0, 0, 4.0),
            Triplet::new(0, 1, 1.0),
            Triplet::new(1, 0, 1.0),
            Triplet::new(1, 1, 4.0),
            // Block 2: rows/cols 2-3
            Triplet::new(2, 2, 4.0),
            Triplet::new(2, 3, 1.0),
            Triplet::new(3, 2, 1.0),
            Triplet::new(3, 3, 4.0),
        ];
        SparseColMat::try_new_from_triplets(4, 4, &triplets).unwrap()
    }

    // ---- Symbolic Analysis with Default Ordering ----

    #[test]
    fn test_analyze_1x1_matrix() {
        let matrix = make_1x1();
        let result = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd);
        assert!(result.is_ok(), "analyze should succeed for 1×1 matrix");
        let sym = result.unwrap();
        assert_eq!(sym.nrows(), 1);
        assert!(sym.predicted_nnz() >= 1);
        let stats = sym.statistics();
        assert_eq!(stats.dimension, 1);
    }

    #[test]
    fn test_analyze_diagonal_matrix() {
        let matrix = make_diagonal(5);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed for diagonal matrix");

        // Diagonal matrix: no fill-in, predicted_nnz should be exactly 5
        assert_eq!(sym.predicted_nnz(), 5);

        // Etree: all roots (parent == -1) for a diagonal matrix
        for &parent in sym.etree() {
            assert_eq!(parent, -1, "diagonal matrix should have all roots in etree");
        }
    }

    #[test]
    fn test_analyze_deterministic() {
        let matrix = make_arrow(8);
        let sym1 = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("first analyze");
        let sym2 = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("second analyze");

        assert_eq!(sym1.predicted_nnz(), sym2.predicted_nnz());
        assert_eq!(sym1.etree(), sym2.etree());
        assert_eq!(sym1.col_counts(), sym2.col_counts());
        assert_eq!(sym1.pivot_buffer_estimates(), sym2.pivot_buffer_estimates());

        let stats1 = sym1.statistics();
        let stats2 = sym2.statistics();
        assert_eq!(stats1.dimension, stats2.dimension);
        assert_eq!(stats1.predicted_nnz, stats2.predicted_nnz);
        assert_eq!(stats1.is_supernodal, stats2.is_supernodal);
    }

    #[test]
    fn test_analyze_predicted_nnz_lower_bound() {
        let matrix = make_arrow(6);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        // Count structural nnz in lower triangle of the original matrix
        let col_ptrs = matrix.symbolic().col_ptr();
        let row_indices = matrix.symbolic().row_idx();
        let mut lower_nnz = 0;
        for j in 0..matrix.ncols() {
            for &row in &row_indices[col_ptrs[j]..col_ptrs[j + 1]] {
                if row >= j {
                    lower_nnz += 1;
                }
            }
        }
        assert!(
            sym.predicted_nnz() >= lower_nnz,
            "predicted_nnz ({}) should be >= structural lower nnz ({})",
            sym.predicted_nnz(),
            lower_nnz
        );
    }

    #[test]
    fn test_analyze_statistics_consistency() {
        let matrix = make_arrow(5);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        let stats = sym.statistics();

        assert_eq!(stats.dimension, sym.nrows());
        assert_eq!(stats.predicted_nnz, sym.predicted_nnz());

        // average_col_count = predicted_nnz / dimension
        let expected_avg = sym.predicted_nnz() as f64 / sym.nrows() as f64;
        assert!(
            (stats.average_col_count - expected_avg).abs() < 1e-12,
            "average_col_count mismatch: {} vs {}",
            stats.average_col_count,
            expected_avg
        );

        // total_pivot_buffer should match
        assert_eq!(stats.total_pivot_buffer, sym.total_pivot_buffer());

        // n_supernodes consistency
        assert_eq!(stats.n_supernodes.is_some(), stats.is_supernodal);
    }

    // ---- Custom Ordering Support ----

    #[test]
    fn test_analyze_custom_identity_ordering() {
        let matrix = make_arrow(5);
        let n = matrix.nrows();

        // Construct identity permutation
        let fwd: Vec<usize> = (0..n).collect();
        let perm = crate::aptp::perm_from_forward(fwd).unwrap();

        let sym =
            AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Custom(perm.as_ref()))
                .expect("analyze with identity custom ordering should succeed");

        // Permutation should reflect identity
        if let Some(p) = sym.perm() {
            let (fwd, _) = p.arrays();
            for (i, &f) in fwd.iter().enumerate() {
                assert_eq!(f, i, "identity ordering should preserve indices");
            }
        }
    }

    #[test]
    fn test_analyze_custom_reverse_ordering() {
        let matrix = make_arrow(5);
        let n = matrix.nrows();

        // Construct reverse permutation
        let fwd: Vec<usize> = (0..n).rev().collect();
        let perm = crate::aptp::perm_from_forward(fwd).unwrap();

        let sym_reverse =
            AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Custom(perm.as_ref()))
                .expect("analyze with reverse ordering should succeed");

        let sym_amd = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze with AMD should succeed");

        // Both should produce valid results, potentially with different nnz
        assert_eq!(sym_reverse.nrows(), n);
        assert_eq!(sym_amd.nrows(), n);
        assert!(sym_reverse.predicted_nnz() > 0);
        assert!(sym_amd.predicted_nnz() > 0);
    }

    #[test]
    fn test_analyze_invalid_perm_dimension() {
        let matrix = make_arrow(5);

        // Construct a permutation with wrong dimension (4 instead of 5)
        let fwd: Vec<usize> = (0..4).collect();
        let perm = crate::aptp::perm_from_forward(fwd).unwrap();

        let result =
            AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Custom(perm.as_ref()));
        assert!(
            result.is_err(),
            "should fail for wrong permutation dimension"
        );
        assert!(
            matches!(result, Err(SparseError::DimensionMismatch { .. })),
            "should return DimensionMismatch error"
        );
    }

    // ---- Supernodal Structure Access ----

    #[test]
    fn test_n_supernodes_simplicial() {
        // A tiny 3×3 matrix that faer classifies as simplicial
        let matrix = make_diagonal(3);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        if !sym.is_supernodal() {
            assert!(
                sym.n_supernodes().is_none(),
                "simplicial should have no supernodes"
            );
        }
        // If faer chose supernodal even for a tiny matrix, that's fine — just verify consistency
        assert_eq!(sym.n_supernodes().is_some(), sym.is_supernodal());
    }

    #[test]
    fn test_supernode_column_ranges_partition() {
        // Use a larger matrix more likely to be supernodal
        let matrix = make_arrow(20);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        if let (Some(begin), Some(end)) = (sym.supernode_begin(), sym.supernode_end()) {
            let ns = sym.n_supernodes().unwrap();
            assert_eq!(begin.len(), ns);
            assert_eq!(end.len(), ns);

            // First supernode starts at 0
            assert_eq!(begin[0], 0, "first supernode should start at column 0");

            // Last supernode ends at dimension
            assert_eq!(
                end[ns - 1],
                sym.nrows(),
                "last supernode should end at dimension"
            );

            // No gaps or overlaps
            for s in 0..ns - 1 {
                assert_eq!(
                    end[s],
                    begin[s + 1],
                    "supernode {} end ({}) should equal supernode {} begin ({})",
                    s,
                    end[s],
                    s + 1,
                    begin[s + 1]
                );
            }

            // Each supernode has at least one column
            for s in 0..ns {
                assert!(
                    begin[s] < end[s],
                    "supernode {} should have at least one column",
                    s
                );
            }
        }
    }

    #[test]
    fn test_assembly_tree_valid() {
        let matrix = make_arrow(20);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        if sym.is_supernodal() {
            let ns = sym.n_supernodes().unwrap();
            let mut root_count = 0;

            for s in 0..ns {
                match sym.supernode_parent(s) {
                    None => root_count += 1,
                    Some(parent) => {
                        assert!(
                            parent > s,
                            "parent of supernode {} should have index > {} (postorder), got {}",
                            s,
                            s,
                            parent
                        );
                        assert!(
                            parent < ns,
                            "parent index {} out of range for {} supernodes",
                            parent,
                            ns
                        );
                    }
                }
            }

            assert!(
                root_count >= 1,
                "assembly tree should have at least one root"
            );
        }
    }

    #[test]
    fn test_supernode_row_pattern_nonempty() {
        let matrix = make_arrow(20);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        if sym.is_supernodal() {
            let ns = sym.n_supernodes().unwrap();
            let mut has_nonempty = false;

            for s in 0..ns {
                let pattern = sym.supernode_pattern(s).unwrap();
                if !pattern.is_empty() {
                    has_nonempty = true;
                }
            }

            // For a non-diagonal matrix, at least one supernode should have off-diagonal entries
            assert!(
                has_nonempty,
                "at least one supernode should have a non-empty row pattern"
            );
        }
    }

    // ---- APTP Pivot Buffer Estimation ----

    #[test]
    fn test_pivot_buffer_nonnegative() {
        for n in [3, 5, 10, 20] {
            let matrix = make_arrow(n);
            let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
                .expect("analyze should succeed");

            assert!(
                !sym.pivot_buffer_estimates().is_empty(),
                "pivot buffer should not be empty for n={}",
                n
            );
            assert!(
                sym.total_pivot_buffer() > 0,
                "total pivot buffer should be > 0 for non-trivial matrix n={}",
                n
            );
        }
    }

    #[test]
    fn test_pivot_buffer_proportional() {
        let matrix = make_arrow(15);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        let total_buffer = sym.total_pivot_buffer();
        let predicted_nnz = sym.predicted_nnz();

        assert!(total_buffer > 0, "total buffer should be > 0");
        assert!(predicted_nnz > 0, "predicted_nnz should be > 0");

        // Buffer should be a reasonable fraction of predicted_nnz.
        // The 10% heuristic is applied to col_counts (which may exceed predicted_nnz
        // for small matrices), so the ratio can exceed 10%. Upper bound is generous.
        let ratio = total_buffer as f64 / predicted_nnz as f64;
        assert!(
            (0.01..=1.0).contains(&ratio),
            "buffer/nnz ratio ({:.4}) should be between 1% and 100%",
            ratio
        );
    }

    #[test]
    fn test_pivot_buffer_reproducible() {
        let matrix = make_arrow(10);
        let sym1 = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("first analyze");
        let sym2 = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("second analyze");

        assert_eq!(
            sym1.pivot_buffer_estimates(),
            sym2.pivot_buffer_estimates(),
            "pivot buffer should be identical for same input"
        );
    }

    #[test]
    fn test_pivot_buffer_length_matches_structure() {
        let matrix = make_arrow(15);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        if sym.is_supernodal() {
            assert_eq!(
                sym.pivot_buffer_estimates().len(),
                sym.n_supernodes().unwrap(),
                "supernodal buffer length should match n_supernodes"
            );
        } else {
            assert_eq!(
                sym.pivot_buffer_estimates().len(),
                sym.nrows(),
                "simplicial buffer length should match nrows"
            );
        }
    }

    // ---- Edge case tests ----

    #[test]
    fn test_analyze_empty_matrix() {
        let triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
        let matrix = SparseColMat::try_new_from_triplets(0, 0, &triplets).unwrap();
        let result = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd);

        // Should either succeed with trivial result or return a descriptive error
        match result {
            Ok(sym) => {
                assert_eq!(sym.nrows(), 0);
                assert_eq!(sym.predicted_nnz(), 0);
            }
            Err(e) => {
                // Acceptable if faer doesn't handle 0×0
                assert!(
                    !format!("{}", e).is_empty(),
                    "error should have descriptive message"
                );
            }
        }
    }

    #[test]
    fn test_analyze_block_diagonal() {
        let matrix = make_block_diagonal();
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed for block-diagonal matrix");

        assert_eq!(sym.nrows(), 4);

        // Etree should have a forest structure (multiple roots)
        let roots: Vec<_> = sym.etree().iter().filter(|&&p| p == -1).collect();
        assert!(
            roots.len() >= 2,
            "block-diagonal should have at least 2 roots, got {}",
            roots.len()
        );
    }

    #[test]
    fn test_analyze_not_square_error() {
        // Create a non-square symbolic structure by constructing a rectangular sparse matrix
        let triplets = vec![Triplet::new(0, 0, 1.0), Triplet::new(1, 1, 1.0)];
        let matrix = SparseColMat::try_new_from_triplets(3, 2, &triplets).unwrap();
        let result = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd);

        assert!(result.is_err(), "should fail for non-square matrix");
        assert!(
            matches!(result, Err(SparseError::NotSquare { .. })),
            "should return NotSquare error"
        );
    }

    // ---- Regression tests with exact expected values ----
    //
    // These tests use SymmetricOrdering::Identity (no AMD reordering) so results
    // are fully deterministic and analytically predictable.

    #[test]
    fn test_regression_diagonal_5() {
        let matrix = make_diagonal(5);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Identity)
            .expect("analyze should succeed");

        // Diagonal: no off-diagonal entries, every column is a root, one nnz per column.
        assert_eq!(sym.etree(), &[-1, -1, -1, -1, -1]);
        assert_eq!(sym.col_counts(), &[1, 1, 1, 1, 1]);
        assert_eq!(sym.predicted_nnz(), 5);
    }

    #[test]
    fn test_regression_tridiagonal_5() {
        let matrix = make_tridiagonal(5);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Identity)
            .expect("analyze should succeed");

        // Tridiagonal: each column j has parent j+1, last is root.
        // L has one sub-diagonal entry per column (except last), so col_counts = [2,2,2,2,1].
        assert_eq!(sym.etree(), &[1, 2, 3, 4, -1]);
        assert_eq!(sym.col_counts(), &[2, 2, 2, 2, 1]);
        assert_eq!(sym.predicted_nnz(), 9);
    }

    #[test]
    fn test_regression_arrow_5_identity() {
        let matrix = make_arrow(5);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Identity)
            .expect("analyze should succeed");

        // Arrow (hub at row/col 0): dense first column causes complete fill-in in L.
        // etree is a chain: 0→1→2→3→4(root).
        // col_counts = [5,4,3,2,1] (dense lower triangle).
        assert_eq!(sym.etree(), &[1, 2, 3, 4, -1]);
        assert_eq!(sym.col_counts(), &[5, 4, 3, 2, 1]);
        assert_eq!(sym.predicted_nnz(), 15);
    }

    #[test]
    fn test_regression_block_diagonal_4() {
        let matrix = make_block_diagonal();
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Identity)
            .expect("analyze should succeed");

        // Two disconnected 2x2 blocks: two independent subtrees, each with one child→parent.
        assert_eq!(sym.etree(), &[1, -1, 3, -1]);
        assert_eq!(sym.col_counts(), &[2, 1, 2, 1]);
        assert_eq!(sym.predicted_nnz(), 6);
    }

    #[test]
    fn test_statistics_display() {
        let matrix = make_arrow(5);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        let stats = sym.statistics();
        let display = format!("{}", stats);
        assert!(display.contains("Dimension:"));
        assert!(display.contains("Predicted NNZ:"));
        assert!(display.contains("Avg col count:"));
    }

    // ---- permute_symbolic_upper_triangle tests ----

    /// Extract CSC col_ptrs and sorted row indices from a SparseColMat for assertion.
    fn extract_csc_structure(mat: &SparseColMat<usize, f64>) -> (Vec<usize>, Vec<usize>) {
        let sym = mat.symbolic();
        let col_ptrs: Vec<usize> = sym.col_ptr().to_vec();
        let row_indices: Vec<usize> = sym.row_idx().to_vec();
        (col_ptrs, row_indices)
    }

    #[test]
    fn test_permute_identity() {
        let matrix = make_arrow(5);
        let perm_fwd: Vec<usize> = (0..5).collect();
        let perm_inv: Vec<usize> = (0..5).collect();

        let permuted = permute_symbolic_upper_triangle(matrix.symbolic(), &perm_fwd, &perm_inv)
            .expect("permute should succeed");

        // Identity permutation: permuted upper triangle == original upper triangle.
        // Arrow(5) upper triangle: diagonal (5 entries) + first row off-diagonal (4 entries) = 9
        let (col_ptrs, row_indices) = extract_csc_structure(&permuted);
        let total_nnz: usize = col_ptrs.last().copied().unwrap_or(0);
        assert_eq!(total_nnz, 9, "arrow(5) upper triangle has 9 entries");

        // Col 0: only diagonal (0)
        assert_eq!(&row_indices[col_ptrs[0]..col_ptrs[1]], &[0]);
        // Cols 1-4: diagonal + entry from row 0 (the hub)
        for j in 1..5 {
            let rows = &row_indices[col_ptrs[j]..col_ptrs[j + 1]];
            assert_eq!(rows, &[0, j], "col {} should have rows [0, {}]", j, j);
        }
    }

    #[test]
    fn test_permute_reverse_arrow() {
        let matrix = make_arrow(5);
        // Reverse: perm_fwd = [4,3,2,1,0], perm_inv = [4,3,2,1,0]
        let perm_fwd: Vec<usize> = (0..5).rev().collect();
        let perm_inv: Vec<usize> = (0..5).rev().collect();

        let permuted = permute_symbolic_upper_triangle(matrix.symbolic(), &perm_fwd, &perm_inv)
            .expect("permute should succeed");

        let (col_ptrs, row_indices) = extract_csc_structure(&permuted);
        let total_nnz: usize = col_ptrs.last().copied().unwrap_or(0);
        assert_eq!(total_nnz, 9, "same number of entries after permutation");

        // Under reverse permutation, the hub (original col 0) moves to col 4.
        // Cols 0-3 should have only diagonal entries.
        for j in 0..4 {
            let rows = &row_indices[col_ptrs[j]..col_ptrs[j + 1]];
            assert_eq!(rows, &[j], "col {} should have only diagonal", j);
        }
        // Col 4 (the hub): rows {0,1,2,3,4}
        let rows_4 = &row_indices[col_ptrs[4]..col_ptrs[5]];
        assert_eq!(rows_4, &[0, 1, 2, 3, 4], "col 4 should be the dense hub");
    }

    #[test]
    fn test_permute_swap_tridiag() {
        let matrix = make_tridiagonal(4);
        // Pair-swap: [1,0,3,2]
        let perm_fwd: Vec<usize> = vec![1, 0, 3, 2];
        let perm_inv: Vec<usize> = vec![1, 0, 3, 2]; // self-inverse

        let permuted = permute_symbolic_upper_triangle(matrix.symbolic(), &perm_fwd, &perm_inv)
            .expect("permute should succeed");

        let (col_ptrs, row_indices) = extract_csc_structure(&permuted);
        let total_nnz: usize = col_ptrs.last().copied().unwrap_or(0);

        // Tridiagonal(4) has 4 diagonal + 3 off-diagonal = 7 upper-triangle entries.
        // Permutation preserves total nnz in the upper triangle.
        assert_eq!(total_nnz, 7, "tridiag(4) upper triangle has 7 entries");

        // Verify every entry is in the upper triangle (row <= col)
        for j in 0..4 {
            for &row in &row_indices[col_ptrs[j]..col_ptrs[j + 1]] {
                assert!(row <= j, "row {} > col {} violates upper triangle", row, j);
            }
        }

        // The swapped tridiagonal should have edges 0-1, 1-2, 2-3 mapped to
        // new edges. Verify the diagonal is present in every column.
        for j in 0..4 {
            let rows = &row_indices[col_ptrs[j]..col_ptrs[j + 1]];
            assert!(
                rows.contains(&j),
                "col {} should contain its diagonal entry",
                j
            );
        }
    }

    #[test]
    fn test_permute_diagonal_invariant() {
        let matrix = make_diagonal(5);
        // Arbitrary permutation
        let perm_fwd: Vec<usize> = vec![3, 1, 4, 0, 2];
        // Compute inverse
        let mut perm_inv = vec![0usize; 5];
        for (new, &old) in perm_fwd.iter().enumerate() {
            perm_inv[old] = new;
        }

        let permuted = permute_symbolic_upper_triangle(matrix.symbolic(), &perm_fwd, &perm_inv)
            .expect("permute should succeed");

        let (col_ptrs, row_indices) = extract_csc_structure(&permuted);
        let total_nnz: usize = col_ptrs.last().copied().unwrap_or(0);

        // Diagonal matrix is permutation-invariant: always n entries, all diagonal.
        assert_eq!(total_nnz, 5, "diagonal(5) always has 5 entries");
        for j in 0..5 {
            let rows = &row_indices[col_ptrs[j]..col_ptrs[j + 1]];
            assert_eq!(rows, &[j], "col {} should have only diagonal entry", j);
        }
    }

    // ---- Cross-validation tests ----

    #[test]
    fn test_col_counts_sum_equals_predicted_nnz() {
        // For each helper matrix, with both Identity and Amd ordering,
        // verify that col_counts sum equals predicted_nnz.
        let matrices: Vec<(&str, SparseColMat<usize, f64>)> = vec![
            ("1x1", make_1x1()),
            ("diagonal(5)", make_diagonal(5)),
            ("tridiagonal(5)", make_tridiagonal(5)),
            ("arrow(5)", make_arrow(5)),
            ("block_diagonal", make_block_diagonal()),
        ];

        for (name, matrix) in &matrices {
            for (ordering_name, ordering) in [
                ("Identity", SymmetricOrdering::Identity),
                ("Amd", SymmetricOrdering::Amd),
            ] {
                let sym = AptpSymbolic::analyze(matrix.symbolic(), ordering).unwrap_or_else(|e| {
                    panic!("analyze failed for {} with {}: {}", name, ordering_name, e)
                });

                let col_counts_sum: usize = sym.col_counts().iter().sum();
                assert_eq!(
                    col_counts_sum,
                    sym.predicted_nnz(),
                    "col_counts sum ({}) != predicted_nnz ({}) for {} with {}",
                    col_counts_sum,
                    sym.predicted_nnz(),
                    name,
                    ordering_name,
                );
            }
        }
    }

    // ---- Non-identity permutation regression test ----

    #[test]
    fn test_regression_arrow_5_reverse() {
        let matrix = make_arrow(5);
        let n = matrix.nrows();

        // Reverse permutation: [4,3,2,1,0]
        let fwd: Vec<usize> = (0..n).rev().collect();
        let perm = crate::aptp::perm_from_forward(fwd).unwrap();

        let sym =
            AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Custom(perm.as_ref()))
                .expect("analyze with reverse ordering should succeed");

        // Under reverse permutation, the hub (original col 0) moves to col 4.
        // This creates a star graph: columns 0-3 are leaves connected only to col 4.
        // No fill-in occurs (leaves are independent — no shared off-diagonal entries).
        // etree: all leaves parent to the hub.
        // col_counts: each leaf has diagonal + one off-diagonal entry to hub (2),
        //             hub has only diagonal (1).
        assert_eq!(sym.etree(), &[4, 4, 4, 4, -1]);
        assert_eq!(sym.col_counts(), &[2, 2, 2, 2, 1]);
        assert_eq!(sym.predicted_nnz(), 9);
    }

    // ---- Supernode parent correctness tests ----

    #[test]
    fn test_supernode_parent_consistent_with_etree() {
        let matrix = make_arrow(20);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        if sym.is_supernodal() {
            let begin = sym.supernode_begin().unwrap();
            let end = sym.supernode_end().unwrap();
            let ns = sym.n_supernodes().unwrap();
            let etree = sym.etree();

            for s in 0..ns {
                let last_col = end[s] - 1;
                let parent_col = etree[last_col];

                match sym.supernode_parent(s) {
                    None => {
                        // Root supernode: etree[last_col] should be -1
                        assert_eq!(
                            parent_col, -1,
                            "supernode {} is root but etree[{}] = {}",
                            s, last_col, parent_col
                        );
                    }
                    Some(parent_sn) => {
                        // Non-root: parent_col should fall within the parent supernode's range
                        let pc = parent_col as usize;
                        assert!(
                            pc >= begin[parent_sn] && pc < end[parent_sn],
                            "supernode {} parent col {} not in supernode {} range [{}, {})",
                            s,
                            pc,
                            parent_sn,
                            begin[parent_sn],
                            end[parent_sn]
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_supernode_parent_child_round_trip() {
        let matrix = make_arrow(20);
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)
            .expect("analyze should succeed");

        if sym.is_supernodal() {
            let ns = sym.n_supernodes().unwrap();

            // Build children map from parent pointers
            let mut children: Vec<Vec<usize>> = vec![Vec::new(); ns];
            let mut root_count = 0;
            for s in 0..ns {
                match sym.supernode_parent(s) {
                    None => root_count += 1,
                    Some(parent) => children[parent].push(s),
                }
            }

            // Every child's parent should point back correctly
            for (parent, kids) in children.iter().enumerate() {
                for &child in kids {
                    assert_eq!(
                        sym.supernode_parent(child),
                        Some(parent),
                        "child {} parent should be {}",
                        child,
                        parent
                    );
                }
            }

            // Total children + roots = total supernodes
            let total_children: usize = children.iter().map(|c| c.len()).sum();
            assert_eq!(
                total_children + root_count,
                ns,
                "children ({}) + roots ({}) should equal total supernodes ({})",
                total_children,
                root_count,
                ns
            );
        }
    }

    #[test]
    fn test_supernode_parent_block_diagonal_roots() {
        let matrix = make_block_diagonal();
        let sym = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Identity)
            .expect("analyze should succeed");

        // Block-diagonal: etree has ≥2 roots (verified in test_regression_block_diagonal_4).
        let etree_roots = sym.etree().iter().filter(|&&p| p == -1).count();
        assert!(
            etree_roots >= 2,
            "block_diagonal should have ≥2 etree roots"
        );

        if sym.is_supernodal() {
            let ns = sym.n_supernodes().unwrap();
            let sn_roots = (0..ns)
                .filter(|&s| sym.supernode_parent(s).is_none())
                .count();
            assert!(
                sn_roots >= 2,
                "supernodal block-diagonal should have ≥2 root supernodes, got {}",
                sn_roots
            );
        }
    }
}
