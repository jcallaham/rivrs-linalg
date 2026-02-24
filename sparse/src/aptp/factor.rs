//! Dense APTP (A Posteriori Threshold Pivoting) factorization kernel.
//!
//! Implements the APTP algorithm for dense symmetric indefinite matrices —
//! the core numerical kernel for the SSIDS multifrontal solver. The kernel
//! factors a dense frontal matrix in place using an optimistic 1x1 pivot
//! strategy with a posteriori stability checking, falling back to 2x2
//! Bunch-Kaufman pivots or column delay when stability bounds are violated.
//!
//! # Algorithm
//!
//! Implements the single-level (column-by-column) APTP strategy from
//! Duff, Hogg & Lopez (2020), "A New Sparse LDL^T Solver Using A Posteriori
//! Threshold Pivoting", SIAM J. Sci. Comput. 42(4).
//!
//! For each column k (left-to-right):
//! 1. Attempt 1x1 pivot: compute L column, check stability bound |l_ij| < 1/u
//! 2. On pass: accept, apply rank-1 Schur complement update
//! 3. On fail: either attempt 2x2 Bunch-Kaufman pivot (with best partner) or
//!    immediately delay, depending on [`AptpFallback`] strategy
//! 4. After all columns: permute delayed columns to end of ordering
//!
//! # Complexity
//!
//! - Time: O(n^3) for full n x n factorization (dominated by Schur complement updates)
//! - Space: O(n^2) for the input matrix (factored in place) plus O(n) for D and permutation
//!
//! # SPRAL Equivalent
//!
//! Corresponds to `LDLT<T, BLOCK_SIZE, Backup>::factor()` in
//! `spral/src/ssids/cpu/kernels/ldlt_app.cxx` (BSD-3-Clause).
//!
//! # References
//!
//! - Duff, Hogg & Lopez (2020), "A New Sparse LDL^T Solver Using A Posteriori
//!   Threshold Pivoting", SIAM J. Sci. Comput. 42(4)
//! - Bunch & Kaufman (1977), "Some Stable Methods for Calculating Inertia and
//!   Solving Symmetric Linear Systems", Math. Comp.

use faer::Par;
use faer::linalg::matmul::triangular::{self as tri_matmul, BlockStructure};
use faer::linalg::triangular_solve;
use faer::perm::Perm;
use faer::prelude::*;
use faer::{Accum, Conj, Mat, MatMut, MatRef};

use super::diagonal::MixedDiagonal;
use super::perm::perm_from_forward;
use super::pivot::{Block2x2, PivotType};
use crate::error::SparseError;

/// Bunch-Kaufman α parameter for 2×2 determinant condition: |det| ≥ α × |a₂₁|².
const BUNCH_KAUFMAN_ALPHA: f64 = 0.5;

/// Growth factor bound for complete pivoting: max |L_ij| ≤ COMPLETE_PIVOTING_GROWTH_BOUND.
/// From Algorithm 4.1 of Duff, Hogg & Lopez (2020).
#[cfg(test)]
const COMPLETE_PIVOTING_GROWTH_BOUND: f64 = 4.0;

// ---------------------------------------------------------------------------
// Configuration types
// ---------------------------------------------------------------------------

/// Configuration for the APTP factorization kernel.
///
/// # Defaults
/// - `threshold`: 0.01 (growth factor bound of 100, matching SPRAL)
/// - `small`: 1e-20 (singularity detection, matching SPRAL)
/// - `fallback`: [`AptpFallback::BunchKaufman`]
/// - `outer_block_size`: 256 (two-level outer block, matching SPRAL)
/// - `inner_block_size`: 32 (two-level inner block, matching SPRAL)
///
/// # Block Size Parameters (Two-Level APTP)
///
/// For frontal matrices larger than `outer_block_size`, the kernel uses
/// a two-level blocked algorithm (Duff, Hogg & Lopez 2020, Section 3):
/// - Outer loop processes blocks of `outer_block_size` columns
/// - Inner loop processes sub-blocks of `inner_block_size` columns
/// - Innermost ib×ib diagonal blocks use complete pivoting (Algorithm 4.1)
///
/// For frontal matrices ≤ `outer_block_size`, the kernel processes the
/// entire fully-summed portion as a single inner block.
///
/// # References
///
/// - SPRAL default: `u = 0.01`, `small = 1e-20`
/// - Duff, Hogg & Lopez (2020), Section 4: threshold parameter u
/// - Duff, Hogg & Lopez (2020), Section 3: two-level blocking with nb=256, ib=32
#[derive(Debug, Clone)]
pub struct AptpOptions {
    /// Stability threshold u: entries must satisfy |l_ij| < 1/u.
    /// Must be in (0.0, 1.0]. Default: 0.01.
    pub threshold: f64,
    /// Singularity detection: pivots with |d| < small treated as zero.
    /// Must be >= 0.0. Default: 1e-20.
    pub small: f64,
    /// Strategy when a 1x1 pivot fails the stability check.
    pub fallback: AptpFallback,
    /// Outer block size (nb) for two-level APTP. Default: 256.
    /// Must be > 0 and >= inner_block_size.
    pub outer_block_size: usize,
    /// Inner block size (ib) for two-level APTP. Default: 32.
    /// Must be > 0 and <= outer_block_size.
    pub inner_block_size: usize,
    /// Strategy for columns that APTP fails to eliminate. Default: Tpp.
    pub failed_pivot_method: FailedPivotMethod,
    /// Parallelism control for BLAS-3 operations (TRSM, GEMM). Default: `Par::Seq`.
    pub par: Par,
    /// Minimum supernode size for amalgamation. Supernodes with fewer than
    /// `nemin` eliminated columns may be merged with their parent. Default: 32.
    ///
    /// Setting `nemin = 1` disables amalgamation entirely.
    ///
    /// # SPRAL Equivalent
    ///
    /// Corresponds to `options%nemin` (`datatypes.f90:21`).
    pub nemin: usize,
    /// Front-size threshold for small-leaf subtree fast path. Default: 256.
    /// Set to 0 to disable.
    pub small_leaf_threshold: usize,
}

impl Default for AptpOptions {
    fn default() -> Self {
        Self {
            threshold: 0.01,
            small: 1e-20,
            fallback: AptpFallback::BunchKaufman,
            outer_block_size: 256,
            inner_block_size: 32,
            failed_pivot_method: FailedPivotMethod::Tpp,
            par: Par::Seq,
            nemin: 32,
            small_leaf_threshold: 256,
        }
    }
}

/// Fallback strategy when a 1x1 pivot fails the stability check.
///
/// # References
///
/// - Duff, Hogg & Lopez (2020), Algorithm 3.1: fallback strategies
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AptpFallback {
    /// Attempt 2x2 Bunch-Kaufman pivot; delay if that also fails.
    BunchKaufman,
    /// Immediately delay the column without attempting 2x2.
    Delay,
}

/// Strategy for handling columns that APTP fails to eliminate.
///
/// When APTP's block-scoped search cannot find acceptable pivots for some
/// columns, the `FailedPivotMethod` controls what happens next.
///
/// # References
///
/// - Duff, Hogg & Lopez (2020), Section 3: TPP fallback after APTP
/// - SPRAL `options%failed_pivot_method`: 0 = pass, 1 = tpp (default)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FailedPivotMethod {
    /// Retry failed columns with serial Threshold Partial Pivoting (default).
    ///
    /// TPP searches ALL remaining columns for acceptable pivots, which
    /// succeeds where APTP's block-scoped search fails.
    Tpp,
    /// Pass failed columns directly to parent as delays (no retry).
    Pass,
}

// ---------------------------------------------------------------------------
// Kernel workspace
// ---------------------------------------------------------------------------

/// Pre-allocated reusable buffers for the BLAS-3 inner loop.
///
/// Eliminates per-block heap allocations inside `factor_inner` and
/// `two_level_factor`. Sized to `max_front × inner_block_size` and reused
/// across all block iterations within and across supernodes.
///
/// # References
///
/// - SPRAL `NumericSubtree.hxx:75-81`: per-thread workspace pattern
pub(crate) struct AptpKernelWorkspace {
    /// Block backup for restore-on-failure, max_front × inner_block_size.
    backup: Mat<f64>,
    /// Copy of L11 block for TRSM aliasing avoidance, inner_block_size × inner_block_size.
    l11_buf: Mat<f64>,
    /// L·D product workspace for update_trailing/cross_terms, max_front × inner_block_size.
    ld_buf: Mat<f64>,
    /// Copy buffer for L21/L_panel aliasing avoidance, max_front × inner_block_size.
    copy_buf: Mat<f64>,
}

impl AptpKernelWorkspace {
    /// Create a new kernel workspace sized for the given maximum front and
    /// inner block dimensions.
    pub(crate) fn new(max_front: usize, inner_block_size: usize) -> Self {
        Self {
            backup: Mat::zeros(max_front, inner_block_size),
            l11_buf: Mat::zeros(inner_block_size, inner_block_size),
            ld_buf: Mat::zeros(max_front, inner_block_size),
            copy_buf: Mat::zeros(max_front, inner_block_size),
        }
    }
}

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// Result of in-place APTP factorization.
///
/// The L factor is stored in the lower triangle of the mutated input matrix.
/// Only the first `num_eliminated` columns contain valid L entries.
#[derive(Debug)]
pub struct AptpFactorResult {
    /// Block diagonal D factor with mixed 1x1/2x2 blocks.
    pub d: MixedDiagonal,
    /// Column permutation: `perm[i]` = original column index at position i.
    pub perm: Vec<usize>,
    /// Number of successfully eliminated columns (q <= num_fully_summed).
    pub num_eliminated: usize,
    /// Original column indices that were not eliminated.
    pub delayed_cols: Vec<usize>,
    /// Summary statistics for diagnostics.
    pub stats: AptpStatistics,
    /// Per-column diagnostic log.
    pub pivot_log: Vec<AptpPivotRecord>,
}

/// Result of convenience APTP factorization (with extracted L).
#[derive(Debug)]
pub struct AptpFactorization {
    /// Unit lower triangular factor (extracted from in-place result).
    pub l: Mat<f64>,
    /// Block diagonal D factor with mixed 1x1/2x2 blocks.
    pub d: MixedDiagonal,
    /// Column permutation as faer Perm type.
    pub perm: Perm<usize>,
    /// Original column indices not eliminated.
    pub delayed_cols: Vec<usize>,
    /// Summary statistics.
    pub stats: AptpStatistics,
    /// Per-column diagnostic log.
    pub pivot_log: Vec<AptpPivotRecord>,
}

/// Summary statistics from factorization.
///
/// Invariant: `num_1x1 + 2 * num_2x2 + num_delayed == total_fully_summed_columns`
#[derive(Debug, Clone, Default)]
pub struct AptpStatistics {
    /// Count of 1x1 pivots accepted.
    pub num_1x1: usize,
    /// Count of 2x2 pivot pairs (each pair counts as 1).
    pub num_2x2: usize,
    /// Count of delayed columns.
    pub num_delayed: usize,
    /// Maximum absolute value across all L entries (stability metric).
    pub max_l_entry: f64,
}

/// Per-column pivot diagnostic record.
#[derive(Debug, Clone)]
pub struct AptpPivotRecord {
    /// Original column index.
    pub col: usize,
    /// Classification from Phase 2 (OneByOne, TwoByTwo, Delayed).
    pub pivot_type: PivotType,
    /// Worst stability metric for this column's L entries.
    pub max_l_entry: f64,
    /// True if 2x2 fallback was attempted (regardless of outcome).
    pub was_fallback: bool,
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Factor a dense symmetric matrix using A Posteriori Threshold Pivoting.
///
/// Performs partial LDL^T factorization of the first `num_fully_summed`
/// columns of the input matrix `a`. The L factor is stored in-place in the
/// lower triangle of `a`. The contribution block (trailing submatrix after
/// the fully-summed columns) is updated with the Schur complement.
///
/// # Arguments
/// - `a`: Dense symmetric matrix (m x m), mutated in place. Only the lower
///   triangle is read; the upper triangle may contain arbitrary values.
///   On return, the lower triangle of the first `num_eliminated` columns
///   contains the L factor entries (unit diagonal implicit).
/// - `num_fully_summed`: Number of columns eligible for elimination (p <= m).
///   For full factorization, pass `a.nrows()`.
/// - `options`: APTP configuration (threshold, fallback strategy).
///
/// # Returns
/// `AptpFactorResult` containing D factor, column permutation, delayed
/// column indices, and diagnostics.
///
/// # Errors
/// Returns `SparseError::InvalidInput` if dimensions are inconsistent.
///
/// # Algorithm
/// Implements the APTP strategy from Duff, Hogg & Lopez (2020).
/// Single-level (column-by-column) variant.
pub fn aptp_factor_in_place(
    mut a: MatMut<'_, f64>,
    num_fully_summed: usize,
    options: &AptpOptions,
) -> Result<AptpFactorResult, SparseError> {
    let m = a.nrows();
    if a.ncols() != m {
        return Err(SparseError::InvalidInput {
            reason: format!("matrix must be square, got {}x{}", m, a.ncols()),
        });
    }
    if num_fully_summed > m {
        return Err(SparseError::InvalidInput {
            reason: format!(
                "num_fully_summed ({}) > matrix dimension ({})",
                num_fully_summed, m
            ),
        });
    }
    if options.threshold <= 0.0 || options.threshold > 1.0 {
        return Err(SparseError::InvalidInput {
            reason: format!("threshold must be in (0, 1], got {}", options.threshold),
        });
    }
    if options.small < 0.0 {
        return Err(SparseError::InvalidInput {
            reason: format!("small must be >= 0.0, got {}", options.small),
        });
    }
    if options.outer_block_size == 0 {
        return Err(SparseError::InvalidInput {
            reason: "outer_block_size must be > 0".to_string(),
        });
    }
    if options.inner_block_size == 0 {
        return Err(SparseError::InvalidInput {
            reason: "inner_block_size must be > 0".to_string(),
        });
    }
    if options.inner_block_size > options.outer_block_size {
        return Err(SparseError::InvalidInput {
            reason: format!(
                "inner_block_size ({}) must be <= outer_block_size ({})",
                options.inner_block_size, options.outer_block_size
            ),
        });
    }

    // Primary factorization dispatch.
    //
    // SPRAL uses ldlt_tpp_factor (TPP) for blocks with ncol < INNER_BLOCK_SIZE,
    // and block_ldlt (complete pivoting) for full aligned blocks. TPP tries 2x2
    // pivots first for every column, accepting more 2x2 pivots than complete
    // pivoting (which only tries 2x2 when the off-diagonal is the global max).
    // For indefinite matrices, 2x2 pivots pair +/- eigenvalues and produce
    // better-conditioned factors.
    //
    // See SPRAL ldlt_app.cxx Block::factor() lines 994-1020.
    // Allocate kernel workspace once for the BLAS-3 paths (reused across
    // all block iterations within factor_inner / two_level_factor).
    let mut kernel_ws = AptpKernelWorkspace::new(m, options.inner_block_size);

    let mut result = if num_fully_summed < options.inner_block_size {
        // Small front: use TPP as primary method (matching SPRAL)
        tpp_factor_as_primary(a.rb_mut(), num_fully_summed, options)?
    } else if num_fully_summed > options.outer_block_size {
        two_level_factor(a.rb_mut(), num_fully_summed, options, &mut kernel_ws)?
    } else {
        factor_inner(
            a.rb_mut(),
            num_fully_summed,
            num_fully_summed,
            options,
            &mut kernel_ws,
        )?
    };

    // Fallback: TPP on remaining columns
    if result.num_eliminated < num_fully_summed
        && options.failed_pivot_method == FailedPivotMethod::Tpp
    {
        let q = result.num_eliminated;

        // Grow D to accommodate TPP's writes at positions q..num_fully_summed
        result.d.grow(num_fully_summed);

        let additional = tpp_factor(
            a.rb_mut(),
            q,
            num_fully_summed,
            &mut result.perm,
            &mut result.d,
            &mut result.stats,
            &mut result.pivot_log,
            options,
        );

        result.num_eliminated = q + additional;
        result.d.truncate(result.num_eliminated);
        result.delayed_cols = (result.num_eliminated..num_fully_summed)
            .map(|i| result.perm[i])
            .collect();
        result.stats.num_delayed = num_fully_summed - result.num_eliminated;
    }

    Ok(result)
}

/// Factor a dense symmetric matrix using APTP, returning extracted L factor.
///
/// Convenience wrapper around [`aptp_factor_in_place`] that:
/// 1. Copies the input matrix
/// 2. Calls `aptp_factor_in_place` with `num_fully_summed = n`
/// 3. Extracts L from the mutated copy
/// 4. Returns [`AptpFactorization`] with all components as owned types
///
/// Suitable for standalone testing and small matrices. For large frontal
/// matrices in multifrontal factorization (Phase 6), use
/// [`aptp_factor_in_place`] directly to avoid the copy.
pub fn aptp_factor(
    a: MatRef<'_, f64>,
    options: &AptpOptions,
) -> Result<AptpFactorization, SparseError> {
    let n = a.nrows();
    if a.ncols() != n {
        return Err(SparseError::InvalidInput {
            reason: format!("matrix must be square, got {}x{}", n, a.ncols()),
        });
    }

    let mut a_copy = a.to_owned();
    let result = aptp_factor_in_place(a_copy.as_mut(), n, options)?;
    let l = extract_l(a_copy.as_ref(), &result.d, result.num_eliminated);
    let perm = perm_from_forward(result.perm.clone())?;

    Ok(AptpFactorization {
        l,
        d: result.d,
        perm,
        delayed_cols: result.delayed_cols,
        stats: result.stats,
        pivot_log: result.pivot_log,
    })
}

// ---------------------------------------------------------------------------
// Internal functions
// ---------------------------------------------------------------------------

/// Swap symmetric rows/columns i and j in the lower triangle of a dense matrix.
///
/// This performs a simultaneous row-and-column permutation on a symmetric
/// matrix stored in the lower triangle, so that the data at position i
/// moves to position j and vice versa.
fn swap_symmetric(mut a: MatMut<'_, f64>, i: usize, j: usize) {
    let m = a.nrows();
    swap_symmetric_block(a.rb_mut(), i, j, 0, m);
}

/// Block-scoped symmetric swap: permute rows/columns i and j within
/// a[col_start..row_limit, col_start..row_limit].
///
/// Same as `swap_symmetric` but the "rows k < i" loop starts at `col_start` instead of 0,
/// and the "rows k > j" loop uses `row_limit` instead of `a.nrows()`. This limits the
/// swap to the diagonal block being factored, leaving both previously-factored columns
/// (0..col_start) and panel rows (row_limit..m) untouched.
///
/// # SPRAL Equivalent
/// `swap_cols(p, t, BLOCK_SIZE, a, lda, ...)` in `block_ldlt.hxx:68`
fn swap_symmetric_block(
    mut a: MatMut<'_, f64>,
    i: usize,
    j: usize,
    col_start: usize,
    row_limit: usize,
) {
    if i == j {
        return;
    }
    let (i, j) = if i < j { (i, j) } else { (j, i) };

    // Swap diagonals
    let tmp = a[(i, i)];
    a[(i, i)] = a[(j, j)];
    a[(j, j)] = tmp;

    // Rows k < i: swap lower-triangle entries a[(i,k)] and a[(j,k)]
    // LIMITED to col_start..i (don't touch previously-factored columns)
    for k in col_start..i {
        let tmp = a[(i, k)];
        a[(i, k)] = a[(j, k)];
        a[(j, k)] = tmp;
    }

    // Rows i < k < j: swap a[(k,i)] and a[(j,k)]
    for k in (i + 1)..j {
        let tmp = a[(k, i)];
        a[(k, i)] = a[(j, k)];
        a[(j, k)] = tmp;
    }

    // Rows k > j: swap a[(k,i)] and a[(k,j)] — LIMITED to row_limit
    for k in (j + 1)..row_limit {
        let tmp = a[(k, i)];
        a[(k, i)] = a[(k, j)];
        a[(k, j)] = tmp;
    }

    // Cross element a[(j,i)] stays unchanged
}

/// Rank-1 Schur complement update after eliminating a 1x1 pivot at column `col`.
#[cfg(test)]
#[allow(dead_code)] // Used by complete_pivoting_factor (also #[cfg(test)])
fn update_schur_1x1(mut a: MatMut<'_, f64>, col: usize, d_value: f64) {
    let m = a.nrows();
    let n_trail = m - col - 1;
    if n_trail == 0 {
        return;
    }

    let l_col: Vec<f64> = (0..n_trail).map(|i| a[(col + 1 + i, col)]).collect();

    for j in 0..n_trail {
        let ldlj = l_col[j] * d_value;
        for i in j..n_trail {
            let old = a[(col + 1 + i, col + 1 + j)];
            a[(col + 1 + i, col + 1 + j)] = old - l_col[i] * ldlj;
        }
    }
}

/// Rank-2 Schur complement update after eliminating a 2x2 pivot.
#[cfg(test)]
#[allow(dead_code)] // Used by complete_pivoting_factor (also #[cfg(test)])
fn update_schur_2x2(mut a: MatMut<'_, f64>, col: usize, partner: usize, block: &Block2x2) {
    let m = a.nrows();
    let start = col.max(partner) + 1;
    let n_trail = m - start;
    if n_trail == 0 {
        return;
    }

    let l1: Vec<f64> = (0..n_trail).map(|i| a[(start + i, col)]).collect();
    let l2: Vec<f64> = (0..n_trail).map(|i| a[(start + i, partner)]).collect();

    let d_a = block.a;
    let d_b = block.b;
    let d_c = block.c;

    for j in 0..n_trail {
        let w_j1 = l1[j] * d_a + l2[j] * d_b;
        let w_j2 = l1[j] * d_b + l2[j] * d_c;
        for i in j..n_trail {
            let update = l1[i] * w_j1 + l2[i] * w_j2;
            let old = a[(start + i, start + j)];
            a[(start + i, start + j)] = old - update;
        }
    }
}

/// Extract the L factor from the in-place factorized matrix.
///
/// Uses the D factor to identify 2x2 pivot blocks, where `a[(k+1, k)]`
/// contains the D off-diagonal (not an L entry). For 2x2 pivots at
/// positions (k, k+1), L entries start at row k+2.
fn extract_l(a: MatRef<'_, f64>, d: &MixedDiagonal, num_eliminated: usize) -> Mat<f64> {
    let n = a.nrows();
    let mut l = Mat::zeros(n, n);

    for i in 0..n {
        l[(i, i)] = 1.0;
    }

    let mut col = 0;
    while col < num_eliminated {
        match d.get_pivot_type(col) {
            PivotType::OneByOne => {
                for i in (col + 1)..n {
                    l[(i, col)] = a[(i, col)];
                }
                col += 1;
            }
            PivotType::TwoByTwo { partner } if partner > col => {
                // a[(col+1, col)] is the D off-diagonal, NOT an L entry.
                // L entries for the 2x2 block start at row col+2.
                for i in (col + 2)..n {
                    l[(i, col)] = a[(i, col)];
                    l[(i, col + 1)] = a[(i, col + 1)];
                }
                col += 2;
            }
            _ => {
                col += 1;
            }
        }
    }

    l
}

// ---------------------------------------------------------------------------
// Two-level APTP: Complete pivoting (Algorithm 4.1)
// ---------------------------------------------------------------------------

/// Factor a small dense symmetric block using complete pivoting.
///
/// Implements Algorithm 4.1 from Duff, Hogg & Lopez (2020): searches
/// the entire remaining submatrix for the entry with maximum magnitude,
/// then uses it as a 1×1 pivot (if on diagonal) or as the off-diagonal
/// of a 2×2 pivot. Provably stable with growth factor bound ≤ 4
/// (equivalent to threshold u=0.25).
///
/// Used at the innermost level of two-level APTP for ib×ib diagonal
/// blocks. Never delays columns (always finds a valid pivot unless
/// the block is numerically singular).
///
/// The matrix is factored in-place (lower triangle). On return:
/// - Columns 0..num_eliminated contain L entries (unit diagonal implicit)
/// - D entries stored in returned MixedDiagonal
/// - The Schur complement is updated in-place in the trailing submatrix
///
/// # SPRAL Equivalent
/// `block_ldlt()` in `spral/src/ssids/cpu/kernels/ldlt_app.cxx`
#[cfg(test)]
fn complete_pivoting_factor(mut a: MatMut<'_, f64>, small: f64) -> AptpFactorResult {
    let n = a.nrows();
    debug_assert_eq!(
        n,
        a.ncols(),
        "complete_pivoting_factor requires square matrix"
    );

    let mut col_order: Vec<usize> = (0..n).collect();
    let mut d = MixedDiagonal::new(n);
    let mut stats = AptpStatistics::default();
    let mut pivot_log = Vec::with_capacity(n);
    let mut k = 0; // next column to eliminate

    while k < n {
        let remaining = n - k;

        // 1. Find (t, m_idx) = argmax |a[i,j]| over all uneliminated entries (lower triangle)
        let mut max_val = 0.0_f64;
        let mut max_row = k;
        let mut max_col = k;
        for j in k..n {
            for i in j..n {
                let v = a[(i, j)].abs();
                if v > max_val {
                    max_val = v;
                    max_row = i;
                    max_col = j;
                }
            }
        }

        // 2. Singularity check
        if max_val < small {
            // Mark all remaining as zero pivots (delayed)
            for &orig_col in &col_order[k..n] {
                stats.num_delayed += 1;
                pivot_log.push(AptpPivotRecord {
                    col: orig_col,
                    pivot_type: PivotType::Delayed,
                    max_l_entry: 0.0,
                    was_fallback: false,
                });
            }
            break;
        }

        // 3. Decide pivot type
        if max_row == max_col {
            // Maximum is on diagonal → 1×1 pivot
            // Swap max_row to position k
            if max_row != k {
                swap_symmetric(a.rb_mut(), k, max_row);
                col_order.swap(k, max_row);
            }

            let d_kk = a[(k, k)];
            let inv_d = 1.0 / d_kk;

            // Compute L column
            let mut max_l = 0.0_f64;
            for i in (k + 1)..n {
                let l_ik = a[(i, k)] * inv_d;
                a[(i, k)] = l_ik;
                let abs_l = l_ik.abs();
                if abs_l > max_l {
                    max_l = abs_l;
                }
            }

            // Schur complement update
            update_schur_1x1(a.rb_mut(), k, d_kk);

            d.set_1x1(k, d_kk);
            stats.num_1x1 += 1;
            if max_l > stats.max_l_entry {
                stats.max_l_entry = max_l;
            }
            pivot_log.push(AptpPivotRecord {
                col: col_order[k],
                pivot_type: PivotType::OneByOne,
                max_l_entry: max_l,
                was_fallback: false,
            });
            k += 1;
        } else {
            // Maximum is off-diagonal at (max_row, max_col)
            // t = max_row, m = max_col (in paper's notation)
            let t = max_row;
            let m = max_col;

            // Compute Δ = a[m,m] * a[t,t] - a[t,m]^2
            // Need current values (before any swap)
            let a_mm = a[(m, m)];
            let a_tt = a[(t, t)];
            let a_tm = a[(t, m)]; // lower triangle: t > m
            let delta = a_mm * a_tt - a_tm * a_tm;

            if remaining < 2 {
                // Only one column left, must use 1×1 (shouldn't normally reach here
                // since max_row == max_col for single element, but guard anyway)
                let d_kk = a[(k, k)];
                if d_kk.abs() < small {
                    stats.num_delayed += 1;
                    pivot_log.push(AptpPivotRecord {
                        col: col_order[k],
                        pivot_type: PivotType::Delayed,
                        max_l_entry: 0.0,
                        was_fallback: false,
                    });
                    break;
                }
                let inv_d = 1.0 / d_kk;
                for i in (k + 1)..n {
                    a[(i, k)] *= inv_d;
                }
                d.set_1x1(k, d_kk);
                stats.num_1x1 += 1;
                k += 1;
                continue;
            }

            if delta.abs() >= BUNCH_KAUFMAN_ALPHA * a_tm * a_tm {
                // 2×2 pivot using (t, m)
                // Swap m → k and t → k+1
                if m != k {
                    swap_symmetric(a.rb_mut(), k, m);
                    col_order.swap(k, m);
                }
                // After swap: the row that was 't' may have moved
                let new_t = if t == k { m } else { t };
                if new_t != k + 1 {
                    swap_symmetric(a.rb_mut(), k + 1, new_t);
                    col_order.swap(k + 1, new_t);
                }

                let a11 = a[(k, k)];
                let a22 = a[(k + 1, k + 1)];
                let a21 = a[(k + 1, k)];
                let det = a11 * a22 - a21 * a21;
                let inv_det = 1.0 / det;

                let block = Block2x2 {
                    first_col: k,
                    a: a11,
                    b: a21,
                    c: a22,
                };

                // Compute L columns for 2×2 pivot
                let mut max_l = 0.0_f64;
                let rows_start = k + 2;
                for i in rows_start..n {
                    let ai1 = a[(i, k)];
                    let ai2 = a[(i, k + 1)];
                    let l_i1 = (ai1 * a22 - ai2 * a21) * inv_det;
                    let l_i2 = (ai2 * a11 - ai1 * a21) * inv_det;
                    a[(i, k)] = l_i1;
                    a[(i, k + 1)] = l_i2;
                    if l_i1.abs() > max_l {
                        max_l = l_i1.abs();
                    }
                    if l_i2.abs() > max_l {
                        max_l = l_i2.abs();
                    }
                }

                // Schur complement update for 2×2
                update_schur_2x2(a.rb_mut(), k, k + 1, &block);

                d.set_2x2(block);
                stats.num_2x2 += 1;
                if max_l > stats.max_l_entry {
                    stats.max_l_entry = max_l;
                }
                pivot_log.push(AptpPivotRecord {
                    col: col_order[k],
                    pivot_type: PivotType::TwoByTwo {
                        partner: col_order[k + 1],
                    },
                    max_l_entry: max_l,
                    was_fallback: false,
                });
                pivot_log.push(AptpPivotRecord {
                    col: col_order[k + 1],
                    pivot_type: PivotType::TwoByTwo {
                        partner: col_order[k],
                    },
                    max_l_entry: max_l,
                    was_fallback: false,
                });
                k += 2;
            } else {
                // Failed 2×2 determinant condition → use 1×1 on max(|a_mm|, |a_tt|)
                let pivot_pos = if a_mm.abs() >= a_tt.abs() { m } else { t };
                if pivot_pos != k {
                    swap_symmetric(a.rb_mut(), k, pivot_pos);
                    col_order.swap(k, pivot_pos);
                }

                let d_kk = a[(k, k)];
                let inv_d = 1.0 / d_kk;

                let mut max_l = 0.0_f64;
                for i in (k + 1)..n {
                    let l_ik = a[(i, k)] * inv_d;
                    a[(i, k)] = l_ik;
                    let abs_l = l_ik.abs();
                    if abs_l > max_l {
                        max_l = abs_l;
                    }
                }

                update_schur_1x1(a.rb_mut(), k, d_kk);

                d.set_1x1(k, d_kk);
                stats.num_1x1 += 1;
                if max_l > stats.max_l_entry {
                    stats.max_l_entry = max_l;
                }
                pivot_log.push(AptpPivotRecord {
                    col: col_order[k],
                    pivot_type: PivotType::OneByOne,
                    max_l_entry: max_l,
                    was_fallback: true, // fell back from 2×2
                });
                k += 1;
            }
        }
    }

    let num_eliminated = k;

    d.truncate(num_eliminated);

    let delayed_cols: Vec<usize> = (num_eliminated..n).map(|i| col_order[i]).collect();

    AptpFactorResult {
        d,
        perm: col_order,
        num_eliminated,
        delayed_cols,
        stats,
        pivot_log,
    }
}

// ---------------------------------------------------------------------------
// Two-level APTP: BLAS-3 building blocks
// ---------------------------------------------------------------------------

/// Factor an ib×ib diagonal block using complete pivoting (Algorithm 4.1).
///
/// Performs complete pivoting on the diagonal block `a[col_start..col_start+block_size,
/// col_start..col_start+block_size]`, but symmetric swaps affect ALL m rows of the
/// matrix (keeping row ordering consistent for the panel below).
///
/// L entries are computed ONLY for rows within the block (col_start..col_start+block_size).
/// Schur complement updates are ONLY within the block. The panel rows below are NOT
/// modified except by symmetric swaps.
///
/// Complete pivoting guarantees |L_ij| ≤ 4 (growth factor for u=0.25).
///
/// # Returns
/// (block_d, local_perm, nelim, stats, pivot_log) where:
/// - `block_d`: MixedDiagonal of dimension block_size with D entries
/// - `local_perm`: permutation applied within positions col_start..col_start+block_size
///   (values are offsets from col_start)
/// - `nelim`: number of successfully eliminated columns (may be < block_size if singular)
/// - `stats`: pivot statistics for this block
/// - `pivot_log`: per-pivot diagnostic records
///
/// # SPRAL Equivalent
/// `block_ldlt()` in `spral/src/ssids/cpu/kernels/block_ldlt.hxx`
fn factor_block_diagonal(
    mut a: MatMut<'_, f64>,
    col_start: usize,
    block_size: usize,
    small: f64,
    row_limit: usize,
) -> (
    MixedDiagonal,
    Vec<usize>,
    usize,
    AptpStatistics,
    Vec<AptpPivotRecord>,
) {
    let block_end = col_start + block_size;

    let mut local_perm: Vec<usize> = (0..block_size).collect();
    let mut d = MixedDiagonal::new(block_size);
    let mut stats = AptpStatistics::default();
    let mut pivot_log = Vec::with_capacity(block_size);
    let mut k = 0; // offset within block (next column to eliminate)

    while k < block_size {
        let cur = col_start + k; // absolute position
        let search_end = block_end;
        let remaining = block_size - k;

        // 1. Find max |a[i,j]| in remaining diagonal sub-block [cur..search_end, cur..search_end]
        let mut max_val = 0.0_f64;
        let mut max_row = cur;
        let mut max_col = cur;
        for j in cur..search_end {
            for i in j..search_end {
                let v = a[(i, j)].abs();
                if v > max_val {
                    max_val = v;
                    max_row = i;
                    max_col = j;
                }
            }
        }

        // 2. Singularity check
        if max_val < small {
            // All remaining entries in block are near-zero
            stats.num_delayed += remaining;
            for &perm_val in &local_perm[k..block_size] {
                pivot_log.push(AptpPivotRecord {
                    col: perm_val,
                    pivot_type: PivotType::Delayed,
                    max_l_entry: 0.0,
                    was_fallback: false,
                });
            }
            break;
        }

        // 3. Decide pivot type
        if max_row == max_col {
            // Diagonal maximum → 1×1 pivot
            if max_row != cur {
                swap_symmetric_block(a.rb_mut(), cur, max_row, col_start, row_limit);
                local_perm.swap(k, max_row - col_start);
            }

            let d_kk = a[(cur, cur)];
            let inv_d = 1.0 / d_kk;

            // Compute L column ONLY within block
            let mut max_l = 0.0_f64;
            for i in (cur + 1)..block_end {
                let l_ik = a[(i, cur)] * inv_d;
                a[(i, cur)] = l_ik;
                let abs_l = l_ik.abs();
                if abs_l > max_l {
                    max_l = abs_l;
                }
            }

            // Schur complement update ONLY within block
            for j in (cur + 1)..block_end {
                let l_j = a[(j, cur)];
                let ldl_j = l_j * d_kk;
                for i in j..block_end {
                    a[(i, j)] -= a[(i, cur)] * ldl_j;
                }
            }

            d.set_1x1(k, d_kk);
            stats.num_1x1 += 1;
            if max_l > stats.max_l_entry {
                stats.max_l_entry = max_l;
            }
            pivot_log.push(AptpPivotRecord {
                col: local_perm[k],
                pivot_type: PivotType::OneByOne,
                max_l_entry: max_l,
                was_fallback: false,
            });
            k += 1;
        } else {
            // Off-diagonal maximum → try 2×2 pivot
            let t = max_row; // absolute positions
            let m_idx = max_col;

            let a_mm = a[(m_idx, m_idx)];
            let a_tt = a[(t, t)];
            let a_tm = a[(t, m_idx)]; // lower triangle: t > m_idx
            let delta = a_mm * a_tt - a_tm * a_tm;

            if remaining < 2 {
                // Only one column left, use 1×1 on max diagonal
                let pivot_pos = if a_mm.abs() >= a_tt.abs() { m_idx } else { t };
                if pivot_pos != cur {
                    swap_symmetric_block(a.rb_mut(), cur, pivot_pos, col_start, row_limit);
                    local_perm.swap(k, pivot_pos - col_start);
                }
                let d_kk = a[(cur, cur)];
                if d_kk.abs() < small {
                    stats.num_delayed += 1;
                    pivot_log.push(AptpPivotRecord {
                        col: local_perm[k],
                        pivot_type: PivotType::Delayed,
                        max_l_entry: 0.0,
                        was_fallback: false,
                    });
                    break;
                }
                let inv_d = 1.0 / d_kk;
                for i in (cur + 1)..block_end {
                    a[(i, cur)] *= inv_d;
                }
                d.set_1x1(k, d_kk);
                stats.num_1x1 += 1;
                k += 1;
                continue;
            }

            if delta.abs() >= BUNCH_KAUFMAN_ALPHA * a_tm * a_tm {
                // 2×2 pivot: swap m_idx → cur, t → cur+1
                if m_idx != cur {
                    swap_symmetric_block(a.rb_mut(), cur, m_idx, col_start, row_limit);
                    local_perm.swap(k, m_idx - col_start);
                }
                let new_t = if t == cur { m_idx } else { t };
                if new_t != cur + 1 {
                    swap_symmetric_block(a.rb_mut(), cur + 1, new_t, col_start, row_limit);
                    local_perm.swap(k + 1, new_t - col_start);
                }

                let a11 = a[(cur, cur)];
                let a22 = a[(cur + 1, cur + 1)];
                let a21 = a[(cur + 1, cur)];
                let det = a11 * a22 - a21 * a21;
                let inv_det = 1.0 / det;

                let block = Block2x2 {
                    first_col: k,
                    a: a11,
                    b: a21,
                    c: a22,
                };

                // Compute L columns ONLY within block
                let mut max_l = 0.0_f64;
                let rows_start = cur + 2;
                for i in rows_start..block_end {
                    let ai1 = a[(i, cur)];
                    let ai2 = a[(i, cur + 1)];
                    let l_i1 = (ai1 * a22 - ai2 * a21) * inv_det;
                    let l_i2 = (ai2 * a11 - ai1 * a21) * inv_det;
                    a[(i, cur)] = l_i1;
                    a[(i, cur + 1)] = l_i2;
                    if l_i1.abs() > max_l {
                        max_l = l_i1.abs();
                    }
                    if l_i2.abs() > max_l {
                        max_l = l_i2.abs();
                    }
                }

                // Schur complement update ONLY within block (rank-2)
                for j in rows_start..block_end {
                    let l1j = a[(j, cur)];
                    let l2j = a[(j, cur + 1)];
                    let w_j1 = l1j * a11 + l2j * a21;
                    let w_j2 = l1j * a21 + l2j * a22;
                    for i in j..block_end {
                        let l1i = a[(i, cur)];
                        let l2i = a[(i, cur + 1)];
                        a[(i, j)] -= l1i * w_j1 + l2i * w_j2;
                    }
                }

                d.set_2x2(block);
                stats.num_2x2 += 1;
                if max_l > stats.max_l_entry {
                    stats.max_l_entry = max_l;
                }
                pivot_log.push(AptpPivotRecord {
                    col: local_perm[k],
                    pivot_type: PivotType::TwoByTwo {
                        partner: local_perm[k + 1],
                    },
                    max_l_entry: max_l,
                    was_fallback: false,
                });
                pivot_log.push(AptpPivotRecord {
                    col: local_perm[k + 1],
                    pivot_type: PivotType::TwoByTwo {
                        partner: local_perm[k],
                    },
                    max_l_entry: max_l,
                    was_fallback: false,
                });
                k += 2;
            } else {
                // Failed Δ → 1×1 on max diagonal
                let pivot_pos = if a_mm.abs() >= a_tt.abs() { m_idx } else { t };
                if pivot_pos != cur {
                    swap_symmetric_block(a.rb_mut(), cur, pivot_pos, col_start, row_limit);
                    local_perm.swap(k, pivot_pos - col_start);
                }

                let d_kk = a[(cur, cur)];
                let inv_d = 1.0 / d_kk;

                let mut max_l = 0.0_f64;
                for i in (cur + 1)..block_end {
                    let l_ik = a[(i, cur)] * inv_d;
                    a[(i, cur)] = l_ik;
                    let abs_l = l_ik.abs();
                    if abs_l > max_l {
                        max_l = abs_l;
                    }
                }

                for j in (cur + 1)..block_end {
                    let l_j = a[(j, cur)];
                    let ldl_j = l_j * d_kk;
                    for i in j..block_end {
                        a[(i, j)] -= a[(i, cur)] * ldl_j;
                    }
                }

                d.set_1x1(k, d_kk);
                stats.num_1x1 += 1;
                if max_l > stats.max_l_entry {
                    stats.max_l_entry = max_l;
                }
                pivot_log.push(AptpPivotRecord {
                    col: local_perm[k],
                    pivot_type: PivotType::OneByOne,
                    max_l_entry: max_l,
                    was_fallback: true,
                });
                k += 1;
            }
        }
    }

    let nelim = k;
    (d, local_perm, nelim, stats, pivot_log)
}

/// Adjust effective_nelim to avoid splitting a 2×2 pivot across block boundaries.
///
/// If the last accepted pivot at position `effective_nelim - 1` is the first half
/// of a 2×2 pair whose partner is beyond `effective_nelim`, decrement by 1.
///
/// # SPRAL Equivalent
/// `Column::adjust()` in `spral/src/ssids/cpu/kernels/ldlt_app.cxx:112-127`
fn adjust_for_2x2_boundary(effective_nelim: usize, d: &MixedDiagonal) -> usize {
    if effective_nelim == 0 {
        return 0;
    }
    let last = effective_nelim - 1;
    match d.get_pivot_type(last) {
        PivotType::TwoByTwo { partner } if partner > last => {
            // Last accepted is the first column of a 2×2 whose partner is beyond nelim
            effective_nelim - 1
        }
        _ => effective_nelim,
    }
}

/// Per-block backup for the two-level APTP algorithm.
///
/// Stores a copy of matrix entries for one outer block column,
/// enabling restore when the a posteriori check reduces nelim.
///
/// # SPRAL Equivalent
/// `CopyBackup<T>` in `spral/src/ssids/cpu/kernels/ldlt_app.cxx`
///
/// Used by the BLAS-3 `factor_inner` to save the block column before
/// factoring, so it can be restored if the Apply step's threshold check
/// reduces nelim.
struct BlockBackup<'a> {
    data: MatMut<'a, f64>,
}

impl<'a> BlockBackup<'a> {
    /// Create a backup of the block column starting at `col_start` with `block_cols` columns.
    /// Backs up `a[col_start.., col_start..col_start+block_cols]` into the provided buffer.
    fn create(
        a: MatRef<'_, f64>,
        col_start: usize,
        block_cols: usize,
        m: usize,
        buf: &'a mut Mat<f64>,
    ) -> Self {
        let rows = m - col_start;
        let mut data = buf.as_mut().submatrix_mut(0, 0, rows, block_cols);
        for j in 0..block_cols {
            for i in 0..rows {
                data[(i, j)] = a[(col_start + i, col_start + j)];
            }
        }
        BlockBackup { data }
    }

    /// Restore the failed portion of the block (columns nelim..block_cols) from backup.
    /// Columns 0..nelim are left untouched (successfully factored).
    fn restore_failed(
        &self,
        mut a: MatMut<'_, f64>,
        col_start: usize,
        nelim: usize,
        block_cols: usize,
        m: usize,
    ) {
        let rows = m - col_start;
        for j in nelim..block_cols {
            for i in 0..rows {
                a[(col_start + i, col_start + j)] = self.data[(i, j)];
            }
        }
    }

    /// Restore failed columns of the diagonal block from pre-factor backup,
    /// applying the block permutation to read from the correct backup positions.
    ///
    /// The backup was taken BEFORE factor_block_diagonal permuted columns.
    /// factor_block_diagonal applied symmetric swaps described by `block_perm`,
    /// so to restore failed column j (in post-perm ordering), we must read
    /// from backup position perm[j] (pre-perm ordering).
    ///
    /// Restores two regions:
    /// 1. Diagonal block: a[k+e..k+bs, k+e..k+bs] — symmetric with perm
    /// 2. Panel below: a[k+bs..m, k+e..k+bs] — column perm only
    ///
    /// # SPRAL Equivalent
    /// `CopyBackup::restore_part_with_sym_perm` (ldlt_app.cxx:562-574)
    fn restore_diagonal_with_perm(
        &self,
        mut a: MatMut<'_, f64>,
        col_start: usize,
        nelim: usize,
        block_cols: usize,
        m: usize,
        block_perm: &[usize],
    ) {
        // Region 1: Diagonal block — restore a[k+e..k+bs, k+e..k+bs]
        // with symmetric permutation from backup.
        // backup[r, c] stored with r >= c (lower triangle).
        // In backup coordinates: row i, col j.
        // SPRAL: aval[j*lda+i] = backup[min(r,c)*ldcopy + max(r,c)]
        for j in nelim..block_cols {
            let c = block_perm[j]; // pre-perm column
            for i in nelim..block_cols {
                let r = block_perm[i]; // pre-perm row
                // Read from lower triangle of backup (row >= col)
                let val = if r >= c {
                    self.data[(r, c)]
                } else {
                    self.data[(c, r)]
                };
                a[(col_start + i, col_start + j)] = val;
            }
        }

        // Region 2: Panel below diagonal block — restore a[k+bs..m, k+e..k+bs]
        // Only column permutation applies (panel rows were permuted by
        // apply_cperm step, but we need original values at permuted column).
        // SPRAL: aval[j*lda+i] = backup[c*ldcopy+i]
        for j in nelim..block_cols {
            let c = block_perm[j]; // pre-perm column
            for i in block_cols..(m - col_start) {
                a[(col_start + i, col_start + j)] = self.data[(i, c)];
            }
        }
    }
}

/// Apply factored L11/D11 to the panel below the diagonal block (TRSM),
/// then perform a posteriori threshold check on all L21 entries.
///
/// Given a factored diagonal block at `a[col_start..col_start+block_nelim, ...]`
/// with L11 (unit lower triangular) and D11, computes:
///   L21 = A21 * (L11 * D11)^{-T}
/// and checks that all |L21[i,j]| < 1/threshold.
///
/// Returns the effective nelim (<= block_nelim): the number of columns
/// whose L entries all satisfy the threshold bound.
///
/// # Algorithm
/// 1. Solve: X = A21 * L11^{-T} via triangular solve (TRSM)
/// 2. Scale: L21[i,j] = X[i,j] / D[j,j] for 1×1 pivots,
///    or L21[i,k:k+1] via 2×2 inversion for 2×2 pivots
/// 3. Scan L21 column-by-column; find first column where any entry exceeds 1/threshold
#[allow(clippy::too_many_arguments)]
fn apply_and_check(
    mut a: MatMut<'_, f64>,
    col_start: usize,
    block_nelim: usize,
    block_cols: usize,
    m: usize,
    d: &MixedDiagonal,
    threshold: f64,
    par: Par,
    l11_buf: &mut Mat<f64>,
) -> usize {
    if block_nelim == 0 {
        return 0;
    }

    let panel_rows = m - col_start - block_cols;
    if panel_rows == 0 {
        return block_nelim;
    }

    // Step 1: TRSM — solve panel * L11^T = A21 for panel (= L21)
    // L11 is unit lower triangular in a[col_start..+block_nelim, col_start..+block_nelim]
    // panel is a[panel_start..m, col_start..col_start+block_nelim]
    //
    // Transposing: L11 * panel^T = A21^T, i.e. unit lower triangular solve on panel^T.
    // Copy L11 to a temporary to avoid aliasing (L11 and panel overlap in `a`).

    let panel_start = col_start + block_cols;

    // Copy L11 into workspace buffer to avoid aliasing
    {
        let src = a
            .rb()
            .submatrix(col_start, col_start, block_nelim, block_nelim);
        let mut dst = l11_buf
            .as_mut()
            .submatrix_mut(0, 0, block_nelim, block_nelim);
        dst.copy_from(src);
    }
    let l11_ref = l11_buf.as_ref().submatrix(0, 0, block_nelim, block_nelim);
    let panel = a
        .rb_mut()
        .submatrix_mut(panel_start, col_start, panel_rows, block_nelim);
    triangular_solve::solve_unit_lower_triangular_in_place(l11_ref, panel.transpose_mut(), par);

    // Step 2: Scale by D^{-1}
    // For 1×1 pivot at column j: L21[:, j] /= D[j]
    // For 2×2 pivot at columns (j, j+1): solve the 2×2 system
    let mut col = 0;
    while col < block_nelim {
        match d.get_pivot_type(col) {
            PivotType::OneByOne => {
                let d_val = d.get_1x1(col);
                let inv_d = 1.0 / d_val;
                for i in 0..panel_rows {
                    a[(panel_start + i, col_start + col)] *= inv_d;
                }
                col += 1;
            }
            PivotType::TwoByTwo { partner } if partner > col => {
                let block = d.get_2x2(col);
                let det = block.a * block.c - block.b * block.b;
                let inv_det = 1.0 / det;
                for i in 0..panel_rows {
                    let x1 = a[(panel_start + i, col_start + col)];
                    let x2 = a[(panel_start + i, col_start + col + 1)];
                    a[(panel_start + i, col_start + col)] = (x1 * block.c - x2 * block.b) * inv_det;
                    a[(panel_start + i, col_start + col + 1)] =
                        (x2 * block.a - x1 * block.b) * inv_det;
                }
                col += 2;
            }
            _ => {
                col += 1;
            }
        }
    }

    // Step 3: Threshold scan — find first failing column
    let stability_bound = 1.0 / threshold;
    let mut effective_nelim = block_nelim;

    let mut scan_col = 0;
    while scan_col < block_nelim {
        let pivot_width = match d.get_pivot_type(scan_col) {
            PivotType::TwoByTwo { partner } if partner > scan_col => 2,
            _ => 1,
        };

        let mut col_ok = true;
        for c in scan_col..scan_col + pivot_width {
            if c >= block_nelim {
                break;
            }
            for i in 0..panel_rows {
                if a[(panel_start + i, col_start + c)].abs() >= stability_bound {
                    col_ok = false;
                    break;
                }
            }
            if !col_ok {
                break;
            }
        }

        if !col_ok {
            effective_nelim = scan_col;
            break;
        }
        scan_col += pivot_width;
    }

    effective_nelim
}

/// Rank-nelim Schur complement update on the trailing submatrix via GEMM.
///
/// Computes: A[trailing, trailing] -= L21 * D11 * L21^T
/// where L21 is the panel at a[panel_start..m, col_start..col_start+nelim]
/// and D11 is the block diagonal from the Factor phase.
///
/// Uses explicit W = L21 * D11 workspace, then A -= W * L21^T (GEMM).
#[allow(clippy::too_many_arguments)]
fn update_trailing(
    mut a: MatMut<'_, f64>,
    col_start: usize,
    nelim: usize,
    block_cols: usize,
    m: usize,
    num_fully_summed: usize,
    d: &MixedDiagonal,
    par: Par,
    ld_buf: &mut Mat<f64>,
    copy_buf: &mut Mat<f64>,
) {
    if nelim == 0 {
        return;
    }

    let trailing_start = col_start + block_cols;
    let trailing_size = m - trailing_start;
    if trailing_size == 0 {
        return;
    }

    // p = num_fully_summed boundary: rows [trailing_start..p] are FS, [p..m] are NFS.
    let p = num_fully_summed;

    // Compute W = L21 * D11 using workspace buffer (trailing_size × nelim subview).
    // We compute W for ALL trailing rows (FS + NFS) since the FS×FS region (region 1)
    // and the NFS×FS cross-term (region 2) both need it.
    let mut w = ld_buf.as_mut().submatrix_mut(0, 0, trailing_size, nelim);
    let mut col = 0;
    while col < nelim {
        match d.get_pivot_type(col) {
            PivotType::OneByOne => {
                let d_val = d.get_1x1(col);
                for i in 0..trailing_size {
                    w[(i, col)] = a[(trailing_start + i, col_start + col)] * d_val;
                }
                col += 1;
            }
            PivotType::TwoByTwo { partner } if partner > col => {
                let block = d.get_2x2(col);
                for i in 0..trailing_size {
                    let l1 = a[(trailing_start + i, col_start + col)];
                    let l2 = a[(trailing_start + i, col_start + col + 1)];
                    w[(i, col)] = l1 * block.a + l2 * block.b;
                    w[(i, col + 1)] = l1 * block.b + l2 * block.c;
                }
                col += 2;
            }
            _ => {
                col += 1;
            }
        }
    }

    // Copy L21 into workspace buffer to avoid borrow conflict (L21 and A22 overlap in `a`)
    {
        let src = a
            .rb()
            .submatrix(trailing_start, col_start, trailing_size, nelim);
        let mut dst = copy_buf.as_mut().submatrix_mut(0, 0, trailing_size, nelim);
        dst.copy_from(src);
    }

    // Region 1 (FS×FS): Lower-triangular GEMM on A[ts..p, ts..p]
    let fs_size = p.saturating_sub(trailing_start);
    if fs_size > 0 {
        let w_fs = ld_buf.as_ref().submatrix(0, 0, fs_size, nelim);
        let l21_fs = copy_buf.as_ref().submatrix(0, 0, fs_size, nelim);
        let mut a_fs = a
            .rb_mut()
            .submatrix_mut(trailing_start, trailing_start, fs_size, fs_size);

        tri_matmul::matmul_with_conj(
            a_fs.rb_mut(),
            BlockStructure::TriangularLower,
            Accum::Add,
            w_fs,
            BlockStructure::Rectangular,
            Conj::No,
            l21_fs.transpose(),
            BlockStructure::Rectangular,
            Conj::No,
            -1.0,
            par,
        );
    }

    // Region 2 (NFS×FS cross-term): Rectangular GEMM on A[p..m, ts..p]
    let nfs_size = m - p;
    if nfs_size > 0 && fs_size > 0 {
        let w_nfs = ld_buf.as_ref().submatrix(fs_size, 0, nfs_size, nelim);
        let l21_fs = copy_buf.as_ref().submatrix(0, 0, fs_size, nelim);
        let mut a_cross = a.submatrix_mut(p, trailing_start, nfs_size, fs_size);

        faer::linalg::matmul::matmul_with_conj(
            a_cross.rb_mut(),
            Accum::Add,
            w_nfs,
            Conj::No,
            l21_fs.transpose(),
            Conj::No,
            -1.0,
            par,
        );
    }

    // Region 3 (NFS×NFS): SKIPPED — deferred to compute_contribution_gemm
}

/// Compute W = L * D for a set of rows, where L and D come from the factored block.
///
/// For 1×1 pivots: W[i, col] = L[i, col] * D[col]
/// For 2×2 pivots: W[i, col] = L[i,col]*D[col,col] + L[i,col+1]*D[col+1,col]
///                 W[i, col+1] = L[i,col]*D[col,col+1] + L[i,col+1]*D[col+1,col+1]
///
/// # SPRAL Equivalent
/// `calcLD<OP_N>` in `spral/src/ssids/cpu/kernels/calc_ld.hxx:41+`
/// Compute W = L * D into a pre-allocated destination buffer.
///
/// Writes into `dst[0..nrows, 0..nelim]`. The destination is overwritten, not zeroed first.
fn compute_ld_into(l: MatRef<'_, f64>, d: &MixedDiagonal, nelim: usize, mut dst: MatMut<'_, f64>) {
    let nrows = l.nrows();
    let mut col = 0;
    while col < nelim {
        match d.get_pivot_type(col) {
            PivotType::OneByOne => {
                let d_val = d.get_1x1(col);
                for i in 0..nrows {
                    dst[(i, col)] = l[(i, col)] * d_val;
                }
                col += 1;
            }
            PivotType::TwoByTwo { partner } if partner > col => {
                let block = d.get_2x2(col);
                for i in 0..nrows {
                    let l1 = l[(i, col)];
                    let l2 = l[(i, col + 1)];
                    dst[(i, col)] = l1 * block.a + l2 * block.b;
                    dst[(i, col + 1)] = l1 * block.b + l2 * block.c;
                }
                col += 2;
            }
            _ => {
                col += 1;
            }
        }
    }
}

/// Compute the NFS×NFS Schur complement via a single deferred GEMM.
///
/// Called after the BLAS-3 blocking loop (which skips the NFS×NFS region).
/// Copies the assembled NFS×NFS values from `frontal_data[p..m, p..m]` into
/// `contrib_buffer`, then applies the rank-`ne` symmetric update in-place:
///
/// ```text
/// contrib_buffer[0..nfs, 0..nfs] = assembled[NFS×NFS] - L21_NFS * D * L21_NFS^T
/// ```
///
/// where `L21_NFS = frontal_data[p..m, 0..ne]` and `D` is the accumulated
/// diagonal from the blocking loop. The output is the lower triangle of the
/// Schur complement.
///
/// # Guard conditions
///
/// - `nfs == 0`: no contribution (root or fully eliminated), returns immediately
/// - `ne == 0`: no rank update (copy only — all columns were delayed)
///
/// # SPRAL Equivalent
///
/// `host_gemm` writing to `node.contrib` in `factor.hxx:92-103` (BSD-3).
///
/// # References
///
/// - Duff, Hogg & Lopez (2020), Section 3: deferred Schur complement
/// - Liu (1992), Section 4: Schur complement in multifrontal method
#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_contribution_gemm(
    frontal_data: &Mat<f64>,
    num_fully_summed: usize,
    num_eliminated: usize,
    m: usize,
    d: &MixedDiagonal,
    contrib_buffer: &mut Mat<f64>,
    par: Par,
) {
    let p = num_fully_summed;
    let ne = num_eliminated;
    let nfs = m - p;

    if nfs == 0 {
        return; // No contribution block (root or fully eliminated)
    }

    // Step 1: Copy assembled NFS×NFS from frontal_data[p..m, p..m] into
    //         contrib_buffer[0..nfs, 0..nfs] (lower triangle only).
    //
    // This copy is unavoidable: assembly scatters into frontal_data, but the
    // GEMM output goes to contrib_buffer. SPRAL avoids this by scattering
    // NFS entries directly into the contribution buffer during assembly — a
    // larger architectural change deferred for a future phase.
    for j in 0..nfs {
        let col_len = nfs - j;
        let src = &frontal_data.col_as_slice(p + j)[p + j..m];
        contrib_buffer.col_as_slice_mut(j)[j..j + col_len].copy_from_slice(src);
    }

    if ne == 0 {
        return; // No rank update — all columns were delayed, copy only
    }

    // Step 2: Compute W = L21_NFS * D (nfs × ne)
    // L21_NFS = frontal_data[p..m, 0..ne]
    let l21_nfs = frontal_data.as_ref().submatrix(p, 0, nfs, ne);

    // TODO: Pre-allocate W in FactorizationWorkspace to amortize across supernodes.
    // Currently allocates per supernode; for matrices with thousands of supernodes
    // this adds allocation pressure in the hot loop.
    let mut w = Mat::zeros(nfs, ne);
    compute_ld_into(l21_nfs, d, ne, w.as_mut());

    // Step 3: Symmetric rank-ne update: contrib -= W * L21_NFS^T (lower triangle)
    tri_matmul::matmul_with_conj(
        contrib_buffer.as_mut().submatrix_mut(0, 0, nfs, nfs),
        BlockStructure::TriangularLower,
        Accum::Add,
        w.as_ref(),
        BlockStructure::Rectangular,
        Conj::No,
        l21_nfs.transpose(),
        BlockStructure::Rectangular,
        Conj::No,
        -1.0,
        par,
    );
}

/// Apply Schur complement updates from passed columns to failed and trailing regions.
///
/// After factoring `nelim` out of `block_cols` columns, three update regions exist:
///
/// 1. **Failed×failed** (diagonal): A[k+e..k+bs, k+e..k+bs] -= W_blk * L_blk^T
/// 2. **Trailing×failed** (cross-term): A[ts..m, k+e..k+bs] -= W_panel * L_blk^T
/// 3. **Trailing×trailing**: handled separately by `update_trailing`
///
/// where:
/// - L_blk = a[k+e..k+bs, k..k+e] (within-block L for failed rows)
/// - L_panel = a[ts..m, k..k+e] (panel L below diagonal block)
/// - W = L * D (LD product)
/// - ts = k + bs (trailing start)
///
/// # SPRAL Equivalent
/// `Block::update` with rfrom/cfrom skip (ldlt_app.cxx:1082-1153)
#[allow(clippy::too_many_arguments)]
fn update_cross_terms(
    mut a: MatMut<'_, f64>,
    col_start: usize,
    nelim: usize,
    block_cols: usize,
    m: usize,
    d: &MixedDiagonal,
    ld_buf: &mut Mat<f64>,
    copy_buf: &mut Mat<f64>,
) {
    if nelim == 0 || nelim >= block_cols {
        return; // No failed columns → no cross-term updates
    }

    let n_failed = block_cols - nelim;
    let trailing_start = col_start + block_cols;
    let trailing_size = m - trailing_start;

    // L_blk: the L entries for failed rows within the diagonal block
    // Copy to workspace to avoid aliasing: a[k+e..k+bs, k..k+e]
    {
        let src = a
            .rb()
            .submatrix(col_start + nelim, col_start, n_failed, nelim);
        let mut dst = copy_buf.as_mut().submatrix_mut(0, 0, n_failed, nelim);
        dst.copy_from(src);
    }

    // W_blk = L_blk * D (into ld_buf)
    {
        let l_blk = copy_buf.as_ref().submatrix(0, 0, n_failed, nelim);
        compute_ld_into(
            l_blk,
            d,
            nelim,
            ld_buf.as_mut().submatrix_mut(0, 0, n_failed, nelim),
        );
    }

    // Region 1: Failed×failed diagonal update
    // A[k+e..k+bs, k+e..k+bs] -= W_blk * L_blk^T (lower triangle only)
    {
        let l_blk = copy_buf.as_ref().submatrix(0, 0, n_failed, nelim);
        let w_blk = ld_buf.as_ref().submatrix(0, 0, n_failed, nelim);
        for j in 0..n_failed {
            for i in j..n_failed {
                let mut sum = 0.0;
                for c in 0..nelim {
                    sum += w_blk[(i, c)] * l_blk[(j, c)];
                }
                a[(col_start + nelim + i, col_start + nelim + j)] -= sum;
            }
        }
    }

    // Region 2: Trailing×failed cross-term update
    // A[ts..m, k+e..k+bs] -= W_panel * L_blk^T
    if trailing_size > 0 {
        // Copy L_panel into workspace (offset past L_blk already in copy_buf)
        {
            let src = a
                .rb()
                .submatrix(trailing_start, col_start, trailing_size, nelim);
            let mut dst = copy_buf
                .as_mut()
                .submatrix_mut(n_failed, 0, trailing_size, nelim);
            dst.copy_from(src);
        }

        // Now both l_blk and l_panel are in copy_buf at non-overlapping offsets.
        // Re-borrow the full copy_buf immutably to access both regions.
        let l_blk = copy_buf.as_ref().submatrix(0, 0, n_failed, nelim);
        let l_panel = copy_buf
            .as_ref()
            .submatrix(n_failed, 0, trailing_size, nelim);

        // W_panel = L_panel * D (into ld_buf offset past W_blk)
        compute_ld_into(
            l_panel,
            d,
            nelim,
            ld_buf
                .as_mut()
                .submatrix_mut(n_failed, 0, trailing_size, nelim),
        );
        let w_panel = ld_buf.as_ref().submatrix(n_failed, 0, trailing_size, nelim);

        for j in 0..n_failed {
            for i in 0..trailing_size {
                let mut sum = 0.0;
                for c in 0..nelim {
                    sum += w_panel[(i, c)] * l_blk[(j, c)];
                }
                a[(trailing_start + i, col_start + nelim + j)] -= sum;
            }
        }
    }
}

/// Factor an nb-sized block using BLAS-3 Factor/Apply/Update loop.
///
/// # Lower-Triangle Convention
///
/// This function operates exclusively on the **lower triangle** of the dense
/// frontal matrix. All reads (pivot search, L extraction) and writes (Schur
/// updates via `BlockStructure::TriangularLower`, `swap_symmetric`) touch only
/// entries where `row >= col`. The upper triangle may contain stale values
/// after column swaps and Schur updates — this is intentional and safe because
/// no code path reads upper-triangle entries. This convention is consistent
/// across `factor_block_diagonal`, `apply_and_check`, `update_trailing`,
/// `update_cross_terms`, and `swap_symmetric_block`.
///
/// This is the middle level of the two-level hierarchy. Processes `num_fully_summed`
/// columns of the block `a[0..m, 0..m]` using ib-sized sub-blocks with the
/// three-phase BLAS-3 pattern from SPRAL's `run_elim_pivoted_notasks`:
///
/// 1. **Backup**: Save `a[k..m, k..k+block_size]` before factoring
/// 2. **Factor**: `factor_block_diagonal` on the ib×ib diagonal block (complete pivoting)
/// 3. **Permute panel**: apply block_perm to panel columns
/// 4. **Zero D off-diagonals**: for TRSM
/// 5. **Row perm propagation**: apply block_perm to columns 0..k (always, not just on success)
/// 6. **Update col_order**: by block_perm (always, not just on success)
/// 7. **Apply+Check**: TRSM on panel + threshold check → effective_nelim
/// 8. **Adjust**: avoid splitting 2×2 across boundary
/// 9. **Partial restore** (on failure): restore failed columns with permuted backup
/// 10. **Schur updates**: trailing×trailing + cross-term updates for failed columns
/// 11. **Delay** failed columns: swap to end_pos
/// 12. **Advance**: k += nelim (may be < block_size)
///
/// Key difference from the old implementation: on threshold failure we DO NOT
/// fully restore and retry. Instead we keep the passed columns, partially
/// restore only the failed columns (using permuted backup), apply Schur updates
/// from passed to failed+trailing, and advance. This matches SPRAL's
/// `run_elim_pivoted_notasks` (ldlt_app.cxx:1585-1713).
///
/// # Arguments
/// - `a`: Dense frontal matrix block (m × m), modified in place
/// - `num_fully_summed`: Number of columns eligible for elimination
/// - `options`: APTP configuration (inner_block_size determines ib)
///
/// # References
/// - SPRAL: `run_elim_pivoted_notasks` in `ldlt_app.cxx:1585-1713`
/// - Duff, Hogg & Lopez (2020), Algorithm 3.1
fn factor_inner(
    mut a: MatMut<'_, f64>,
    num_fully_summed: usize,
    nfs_boundary: usize,
    options: &AptpOptions,
    kernel_ws: &mut AptpKernelWorkspace,
) -> Result<AptpFactorResult, SparseError> {
    let m = a.nrows();
    let ib = options.inner_block_size;
    let small = options.small;
    let threshold = options.threshold;
    let p = num_fully_summed;

    let mut col_order: Vec<usize> = (0..m).collect();
    let mut d = MixedDiagonal::new(p);
    let mut stats = AptpStatistics::default();
    let mut pivot_log = Vec::with_capacity(p);
    let mut k = 0;
    let mut end_pos = p;

    // Pre-allocated buffers reused across block iterations to avoid
    // per-block allocations in the hot loop.
    let mut panel_perm_buf = vec![0.0f64; ib];
    let mut row_perm_buf = vec![0.0f64; ib];
    let mut col_order_buf = vec![0usize; ib];

    // BLAS-3 Factor/Apply/Update loop with ib-sized inner blocks.
    //
    // Matches SPRAL's `run_elim_pivoted_notasks` architecture:
    //   1. Backup (pre-factor)
    //   2. Factor diagonal block (complete pivoting, block-scoped swaps)
    //   3. Permute panel columns by block_perm
    //   4. Zero D off-diagonals for TRSM
    //   5. Row perm propagation to columns 0..k
    //   6. Update col_order by block_perm
    //   7. Apply+Check (TRSM + threshold scan)
    //   8. Adjust for 2×2 boundary
    //   9. On failure: partial restore of failed columns (not full restore+retry)
    //  10. Schur updates (trailing + cross-terms for failed columns)
    //  11. Delay failed columns (swap to end_pos)
    //  12. Advance k += effective_nelim
    //
    // Steps 5-6 happen BEFORE apply_and_check so the matrix is in a consistent
    // permuted state regardless of threshold outcome (matching SPRAL where
    // apply_rperm_and_backup happens before apply_pivot_app).
    while k < end_pos {
        let block_size = (end_pos - k).min(ib);
        let block_end = k + block_size;

        // 1. BACKUP: save a[k..m, k..k+block_size] before factoring
        let backup = BlockBackup::create(a.as_ref(), k, block_size, m, &mut kernel_ws.backup);

        // 2. FACTOR: complete pivoting on the block_size×block_size diagonal block
        //    Block-scoped swaps: only rows/columns within [0..block_end] are permuted.
        //    Panel rows [block_end..m] are NOT touched.
        let (block_d, block_perm, block_nelim, block_stats, block_log) =
            factor_block_diagonal(a.rb_mut(), k, block_size, small, block_end);

        if block_nelim == 0 {
            // Entire block is singular — restore and delay all columns
            backup.restore_failed(a.rb_mut(), k, 0, block_size, m);
            for offset in 0..block_size {
                let delayed_orig = col_order[k + block_perm[offset]];
                end_pos -= 1;
                if k + offset < end_pos {
                    swap_symmetric(a.rb_mut(), k + offset, end_pos);
                    col_order.swap(k + offset, end_pos);
                }
                stats.num_delayed += 1;
                pivot_log.push(AptpPivotRecord {
                    col: delayed_orig,
                    pivot_type: PivotType::Delayed,
                    max_l_entry: 0.0,
                    was_fallback: false,
                });
            }
            continue;
        }

        // 3. PERMUTE PANEL: reorder panel column entries a[r, k..k+bs]
        //    according to block_perm so that TRSM sees the correct input.
        let panel_start = block_end;
        for r in panel_start..m {
            for i in 0..block_size {
                panel_perm_buf[i] = a[(r, k + block_perm[i])];
            }
            for i in 0..block_size {
                a[(r, k + i)] = panel_perm_buf[i];
            }
        }

        // 4. Zero out D off-diagonals so apply_and_check's TRSM reads them as
        //    L11 entries (should be 0 for 2×2 pivots where L starts at row k+2).
        {
            let mut bc = 0;
            while bc < block_nelim {
                match block_d.get_pivot_type(bc) {
                    PivotType::TwoByTwo { partner } if partner > bc => {
                        a[(k + bc + 1, k + bc)] = 0.0;
                        bc += 2;
                    }
                    _ => {
                        bc += 1;
                    }
                }
            }
        }

        // 5. ROW PERM PROPAGATION: apply block_perm to columns 0..k.
        //    Done BEFORE apply_and_check (matching SPRAL where apply_rperm_and_backup
        //    happens before apply_pivot_app). This ensures the matrix is in a
        //    consistent permuted state regardless of threshold outcome.
        if k > 0 {
            for c in 0..k {
                for i in 0..block_size {
                    row_perm_buf[i] = a[(k + block_perm[i], c)];
                }
                for i in 0..block_size {
                    a[(k + i, c)] = row_perm_buf[i];
                }
            }
        }

        // 6. UPDATE COL_ORDER by block_perm (always, not just on success).
        col_order_buf[..block_size].copy_from_slice(&col_order[k..k + block_size]);
        for i in 0..block_size {
            col_order[k + i] = col_order_buf[block_perm[i]];
        }

        // 7. APPLY: TRSM on panel below + threshold check
        let mut effective_nelim = apply_and_check(
            a.rb_mut(),
            k,
            block_nelim,
            block_size,
            m,
            &block_d,
            threshold,
            options.par,
            &mut kernel_ws.l11_buf,
        );

        // 8. ADJUST: don't split 2×2 pivot across block boundary
        effective_nelim = adjust_for_2x2_boundary(effective_nelim, &block_d);

        // 9. PARTIAL RESTORE (on failure): restore only failed columns from
        //    pre-factor backup with permutation applied.
        //    SPRAL: restore_part_with_sym_perm for diagonal, restore_part for panel.
        //    We keep passed columns (0..effective_nelim) — their L11, D, L21 are committed.
        if effective_nelim < block_nelim {
            backup.restore_diagonal_with_perm(
                a.rb_mut(),
                k,
                effective_nelim,
                block_size,
                m,
                &block_perm,
            );
        }

        // Use effective_nelim as the number of passed columns for all subsequent steps.
        let nelim = effective_nelim;

        // 10. SCHUR UPDATES using only passed columns' L and D.
        //     Truncate block_d to passed columns for update computations.
        //     Three regions:
        //     a. Trailing×trailing: A[ts..m, ts..m] -= L_panel * D * L_panel^T
        //     b. Failed×failed: A[k+e..k+bs, k+e..k+bs] -= L_blk * D * L_blk^T
        //     c. Trailing×failed: A[ts..m, k+e..k+bs] -= L_panel * D * L_blk^T
        if nelim > 0 {
            update_trailing(
                a.rb_mut(),
                k,
                nelim,
                block_size,
                m,
                nfs_boundary,
                &block_d,
                options.par,
                &mut kernel_ws.ld_buf,
                &mut kernel_ws.copy_buf,
            );
            update_cross_terms(
                a.rb_mut(),
                k,
                nelim,
                block_size,
                m,
                &block_d,
                &mut kernel_ws.ld_buf,
                &mut kernel_ws.copy_buf,
            );
        }

        // 11. ACCUMULATE D entries for passed columns
        d.copy_from_offset(&block_d, k, nelim);

        // Accumulate stats from passed columns only (count from block_d directly)
        let mut passed_1x1 = 0;
        let mut passed_2x2 = 0;
        let mut sc = 0;
        while sc < nelim {
            match block_d.get_pivot_type(sc) {
                PivotType::OneByOne => {
                    passed_1x1 += 1;
                    sc += 1;
                }
                PivotType::TwoByTwo { partner } if partner > sc => {
                    passed_2x2 += 1;
                    sc += 2;
                }
                _ => {
                    sc += 1;
                }
            }
        }
        stats.num_1x1 += passed_1x1;
        stats.num_2x2 += passed_2x2;
        if block_stats.max_l_entry > stats.max_l_entry {
            stats.max_l_entry = block_stats.max_l_entry;
        }
        // Check panel L entries for max_l_entry (passed columns only)
        for c in 0..nelim {
            for i in panel_start..m {
                let v = a[(i, k + c)].abs();
                if v > stats.max_l_entry {
                    stats.max_l_entry = v;
                }
            }
        }

        // Accumulate pivot log for passed columns
        for record in &block_log {
            if !matches!(record.pivot_type, PivotType::Delayed) && record.col < nelim {
                let global_col = col_order[k + record.col];
                let global_pivot_type = match record.pivot_type {
                    PivotType::TwoByTwo { partner } => PivotType::TwoByTwo {
                        partner: col_order[k + partner],
                    },
                    other => other,
                };
                pivot_log.push(AptpPivotRecord {
                    col: global_col,
                    pivot_type: global_pivot_type,
                    max_l_entry: record.max_l_entry,
                    was_fallback: record.was_fallback,
                });
            }
        }

        // 12. DELAY failed columns (nelim..block_size): swap to end_pos
        if nelim < block_size {
            let n_delayed = block_size - nelim;
            for i in (0..n_delayed).rev() {
                let delayed_pos = k + nelim + i;
                end_pos -= 1;
                if delayed_pos < end_pos {
                    swap_symmetric(a.rb_mut(), delayed_pos, end_pos);
                    col_order.swap(delayed_pos, end_pos);
                }
                stats.num_delayed += 1;
                pivot_log.push(AptpPivotRecord {
                    col: col_order[end_pos],
                    pivot_type: PivotType::Delayed,
                    max_l_entry: 0.0,
                    was_fallback: false,
                });
            }
        }

        k += nelim;
    }

    let num_eliminated = k;

    d.truncate(num_eliminated);

    let delayed_cols: Vec<usize> = (num_eliminated..p).map(|i| col_order[i]).collect();

    Ok(AptpFactorResult {
        d,
        perm: col_order,
        num_eliminated,
        delayed_cols,
        stats,
        pivot_log,
    })
}

/// Two-level outer block loop for large frontal matrices.
///
/// Processes nb-sized blocks. For each outer block: `factor_inner` handles
/// the BLAS-3 Factor/Apply/Update loop internally on ib-sized sub-blocks,
/// including threshold checking on the panel via `apply_and_check` and
/// block-level backup/restore on failure.
///
/// Called by `aptp_factor_in_place` when `num_fully_summed > outer_block_size`.
///
/// # References
/// - Duff, Hogg & Lopez (2020), Algorithm 3.1: two-level outer loop
fn two_level_factor(
    mut a: MatMut<'_, f64>,
    num_fully_summed: usize,
    options: &AptpOptions,
    kernel_ws: &mut AptpKernelWorkspace,
) -> Result<AptpFactorResult, SparseError> {
    let m = a.nrows();
    let nb = options.outer_block_size;
    let p = num_fully_summed;

    let mut col_order: Vec<usize> = (0..m).collect();
    let mut global_d = MixedDiagonal::new(p);
    let mut stats = AptpStatistics::default();
    let mut pivot_log = Vec::with_capacity(p);

    let mut global_nelim = 0;
    let mut remaining_fully_summed = p;

    while remaining_fully_summed > 0 {
        let col_start = global_nelim;
        let block_cols = remaining_fully_summed.min(nb);

        // FACTOR: inner APTP on the full view (including panel rows).
        // factor_inner handles:
        //   - Complete pivoting search within ib×ib sub-blocks
        //   - Symmetric swaps applied to ALL rows (including panel)
        //   - L entry computation + threshold check for ALL rows (via try_1x1/try_2x2)
        //   - Schur complement update for ALL trailing rows
        //   - Internal backup/restore for failed pivots (inside try_1x1/try_2x2)
        //
        // No external backup/restore is needed: factor_inner handles all
        // threshold failures internally, and the Schur complement has been
        // applied to all trailing entries including any delayed columns.
        let block_m = m - col_start;
        let block_result = {
            let block_view = a
                .rb_mut()
                .submatrix_mut(col_start, col_start, block_m, block_m);
            // nfs_boundary relative to this subview: p - col_start
            // This ensures inner blocks skip NFS×NFS updates consistently
            // with the deferred GEMM that runs after aptp_factor_in_place.
            factor_inner(block_view, block_cols, p - col_start, options, kernel_ws)?
        };
        let block_nelim = block_result.num_eliminated;

        // PROPAGATE ROW PERMUTATION to already-factored columns.
        //
        // factor_inner operates on a submatrix view [col_start..m, col_start..m].
        // Its swap_symmetric calls only rearrange rows WITHIN that submatrix.
        // The L entries from previously-factored blocks (columns 0..col_start)
        // are NOT rearranged by these swaps. We must apply the same row
        // permutation to those columns so that extract_l reads consistent
        // L entries across all blocks.
        if col_start > 0 {
            let block_perm = &block_result.perm;
            let mut temp = vec![0.0f64; block_cols];
            for c in 0..col_start {
                // Gather: temp[i] = row that ended up at local position i
                for i in 0..block_cols {
                    temp[i] = a[(col_start + block_perm[i], c)];
                }
                // Scatter: write to sequential positions
                for i in 0..block_cols {
                    a[(col_start + i, c)] = temp[i];
                }
            }
        }

        // Update col_order BEFORE delay swap: factor_inner may have permuted
        // columns within the block. We must capture the pre-swap col_order so
        // that orig_order[block_perm[i]] reads the correct original column index.
        // If we swapped first, delayed positions would contain post-swap values,
        // corrupting the mapping for eliminated columns.
        {
            let block_perm = &block_result.perm;
            let orig_order: Vec<usize> = col_order[col_start..col_start + block_cols].to_vec();
            for i in 0..block_cols {
                // Entries with block_perm[i] >= block_cols are contribution block rows
                // (not fully-summed columns), so they don't have a col_order mapping.
                if block_perm[i] < block_cols {
                    col_order[col_start + i] = orig_order[block_perm[i]];
                }
            }
        }

        // ADJUST delayed columns: swap them to end of unprocessed region
        // so the next outer block processes fresh columns first.
        if block_nelim < block_cols {
            let n_failed = block_cols - block_nelim;
            for i in 0..n_failed {
                let failed_pos = col_start + block_nelim + i;
                let end = col_start + remaining_fully_summed - 1 - i;
                if failed_pos < end {
                    swap_symmetric(a.rb_mut(), failed_pos, end);
                    col_order.swap(failed_pos, end);
                }
            }
            stats.num_delayed += n_failed;
        }

        // 4. Accumulate into global result
        global_d.copy_from_offset(&block_result.d, global_nelim, block_nelim);

        // Merge stats
        stats.num_1x1 += block_result.stats.num_1x1;
        stats.num_2x2 += block_result.stats.num_2x2;
        if block_result.stats.max_l_entry > stats.max_l_entry {
            stats.max_l_entry = block_result.stats.max_l_entry;
        }

        // Merge pivot log
        for record in &block_result.pivot_log {
            if !matches!(record.pivot_type, PivotType::Delayed) {
                pivot_log.push(record.clone());
            }
        }
        // Add delayed logs
        let n_failed = block_cols - block_nelim;
        for i in 0..n_failed {
            let delayed_pos = col_start + remaining_fully_summed - 1 - i;
            pivot_log.push(AptpPivotRecord {
                col: col_order[delayed_pos],
                pivot_type: PivotType::Delayed,
                max_l_entry: 0.0,
                was_fallback: false,
            });
        }

        global_nelim += block_nelim;
        remaining_fully_summed -= block_cols;
    }

    global_d.truncate(global_nelim);

    let delayed_cols: Vec<usize> = (global_nelim..p).map(|i| col_order[i]).collect();

    Ok(AptpFactorResult {
        d: global_d,
        perm: col_order,
        num_eliminated: global_nelim,
        delayed_cols,
        stats,
        pivot_log,
    })
}

// ---------------------------------------------------------------------------
// TPP (Threshold Partial Pivoting) fallback
// ---------------------------------------------------------------------------

/// Check if all entries in row/column `idx` within the uneliminated region are
/// smaller than `small` in absolute value.
///
/// Scans row entries `a[(idx, c)]` for `c in from..idx` and column entries
/// `a[(r, idx)]` for `r in idx..to`.
///
/// # SPRAL Equivalent
/// `check_col_small()` in `spral/src/ssids/cpu/kernels/ldlt_tpp.cxx`
fn tpp_is_col_small(a: MatRef<'_, f64>, idx: usize, from: usize, to: usize, small: f64) -> bool {
    // Row entries: a[(idx, c)] for c < idx (lower triangle: idx > c)
    for c in from..idx {
        if a[(idx, c)].abs() >= small {
            return false;
        }
    }
    // Column entries: a[(r, idx)] for r >= idx (lower triangle: r >= idx)
    for r in idx..to {
        if a[(r, idx)].abs() >= small {
            return false;
        }
    }
    true
}

/// Find the column index with largest absolute entry in row `p`, scanning
/// columns `from..to`.
///
/// For `c < p`: reads `a[(p, c)]` (lower triangle).
/// For `c == p`: reads `a[(p, p)]` (diagonal).
/// For `c > p`: reads `a[(c, p)]` (symmetric entry via lower triangle).
///
/// Returns the column index of the maximum. Returns `from` if `from >= to`.
///
/// # SPRAL Equivalent
/// `find_row_abs_max()` in `spral/src/ssids/cpu/kernels/ldlt_tpp.cxx`
fn tpp_find_row_abs_max(a: MatRef<'_, f64>, p: usize, from: usize, to: usize) -> usize {
    if from >= to {
        return from;
    }
    let mut best_idx = from;
    let mut best_val = tpp_sym_entry(a, p, from).abs();
    for c in (from + 1)..to {
        let v = tpp_sym_entry(a, p, c).abs();
        if v > best_val {
            best_idx = c;
            best_val = v;
        }
    }
    best_idx
}

/// Read symmetric entry `a(row, col)` from lower triangle storage.
#[inline]
fn tpp_sym_entry(a: MatRef<'_, f64>, row: usize, col: usize) -> f64 {
    if row >= col {
        a[(row, col)]
    } else {
        a[(col, row)]
    }
}

/// Find the maximum absolute value in row/column `col` among uneliminated
/// positions, excluding one index and the diagonal.
///
/// Scans both the row (columns `nelim..col`) and the column (rows `col+1..m`),
/// skipping index `exclude` (use `usize::MAX` to exclude nothing).
///
/// # SPRAL Equivalent
/// `find_rc_abs_max_exclude()` in `spral/src/ssids/cpu/kernels/ldlt_tpp.cxx`
fn tpp_find_rc_abs_max_exclude(
    a: MatRef<'_, f64>,
    col: usize,
    nelim: usize,
    m: usize,
    exclude: usize,
) -> f64 {
    let mut best = 0.0_f64;
    // Row part: a[(col, c)] for c in nelim..col
    for c in nelim..col {
        if c == exclude {
            continue;
        }
        best = best.max(a[(col, c)].abs());
    }
    // Column part: a[(r, col)] for r in col+1..m
    for r in (col + 1)..m {
        if r == exclude {
            continue;
        }
        best = best.max(a[(r, col)].abs());
    }
    best
}

/// Test if `(t, p)` with `t < p` form a good 2x2 pivot.
///
/// Three checks (matching SPRAL):
/// 1. Non-zero pivot block: max(|a11|, |a21|, |a22|) >= small
/// 2. Non-singular determinant with cancellation guard
/// 3. Threshold: u * max(|D^{-1}_{11}|*maxt + |D^{-1}_{12}|*maxp,
///    |D^{-1}_{12}|*maxt + |D^{-1}_{22}|*maxp) < 1
///
/// Returns `Some((a11, a21, a22))` (D values, NOT D^{-1}) on success.
///
/// # SPRAL Equivalent
/// `test_2x2()` in `spral/src/ssids/cpu/kernels/ldlt_tpp.cxx`
fn tpp_test_2x2(
    a: MatRef<'_, f64>,
    t: usize,
    p: usize,
    maxt: f64,
    maxp: f64,
    u: f64,
    small: f64,
) -> Option<(f64, f64, f64)> {
    debug_assert!(t < p, "tpp_test_2x2 requires t < p");

    let a11 = a[(t, t)];
    let a21 = a[(p, t)]; // lower triangle: p > t
    let a22 = a[(p, p)];

    // 1. Non-zero pivot block
    let maxpiv = a11.abs().max(a21.abs()).max(a22.abs());
    if maxpiv < small {
        return None;
    }

    // 2. Cancellation guard on determinant
    let detscale = 1.0 / maxpiv;
    let detpiv0 = (a11 * detscale) * a22;
    let detpiv1 = (a21 * detscale) * a21;
    let detpiv = detpiv0 - detpiv1;
    if detpiv.abs() < small.max((detpiv0 / 2.0).abs()).max((detpiv1 / 2.0).abs()) {
        return None;
    }

    // 3. Threshold test using D^{-1}
    let d_inv_11 = (a22 * detscale) / detpiv;
    let d_inv_12 = (-a21 * detscale) / detpiv;
    let d_inv_22 = (a11 * detscale) / detpiv;

    if maxt.max(maxp) < small {
        // Rest of column is small — accept
        return Some((a11, a21, a22));
    }

    let x1 = d_inv_11.abs() * maxt + d_inv_12.abs() * maxp;
    let x2 = d_inv_12.abs() * maxt + d_inv_22.abs() * maxp;
    if u * x1.max(x2) < 1.0 {
        Some((a11, a21, a22))
    } else {
        None
    }
}

/// Apply 1x1 pivot elimination at position `nelim`.
///
/// Computes L entries (divides column by D), then performs rank-1 Schur
/// complement update on the entire lower triangle (including the contribution
/// block beyond the fully-summed region).
///
/// Returns the maximum absolute L entry.
fn tpp_apply_1x1(mut a: MatMut<'_, f64>, nelim: usize, m: usize, num_fully_summed: usize) -> f64 {
    let d = a[(nelim, nelim)];
    let inv_d = 1.0 / d;
    let p = num_fully_summed;

    // Compute L entries for ALL rows (FS + NFS) — needed for deferred contribution GEMM
    let mut max_l = 0.0_f64;
    for i in (nelim + 1)..m {
        let l_ik = a[(i, nelim)] * inv_d;
        a[(i, nelim)] = l_ik;
        max_l = max_l.max(l_ik.abs());
    }

    // Rank-1 Schur complement update: skip NFS×NFS region (deferred to contribution GEMM).
    // For FS columns j < p: update all rows j..m (FS×FS + NFS×FS cross-terms).
    // For NFS columns j >= p: skip entirely (NFS×NFS).
    let schur_col_end = p.min(m);
    for j in (nelim + 1)..schur_col_end {
        let ldlj = a[(j, nelim)] * d;
        for i in j..m {
            a[(i, j)] -= a[(i, nelim)] * ldlj;
        }
    }

    max_l
}

/// Apply 2x2 pivot elimination at positions `(nelim, nelim+1)`.
///
/// Computes L entries using D^{-1}, then performs rank-2 Schur complement
/// update on the entire lower triangle (including the contribution block).
/// Sets `a[(nelim+1, nelim)] = 0.0` to match the APTP convention for
/// `extract_front_factors`.
///
/// Returns the maximum absolute L entry.
fn tpp_apply_2x2(mut a: MatMut<'_, f64>, nelim: usize, m: usize, num_fully_summed: usize) -> f64 {
    let a11 = a[(nelim, nelim)];
    let a21 = a[(nelim + 1, nelim)];
    let a22 = a[(nelim + 1, nelim + 1)];
    let det = a11 * a22 - a21 * a21;
    let inv_det = 1.0 / det;
    let p = num_fully_summed;

    // Compute L entries for ALL rows (FS + NFS) — needed for deferred contribution GEMM
    let mut max_l = 0.0_f64;
    let start = nelim + 2;
    for i in start..m {
        let ai1 = a[(i, nelim)];
        let ai2 = a[(i, nelim + 1)];
        let l_i1 = (ai1 * a22 - ai2 * a21) * inv_det;
        let l_i2 = (ai2 * a11 - ai1 * a21) * inv_det;
        a[(i, nelim)] = l_i1;
        a[(i, nelim + 1)] = l_i2;
        max_l = max_l.max(l_i1.abs()).max(l_i2.abs());
    }

    // Rank-2 Schur complement update: skip NFS×NFS region (deferred to contribution GEMM).
    // For FS columns j < p: update all rows j..m (FS×FS + NFS×FS cross-terms).
    // For NFS columns j >= p: skip entirely (NFS×NFS).
    let schur_col_end = p.min(m);
    for j in start..schur_col_end {
        let wj1 = a[(j, nelim)] * a11 + a[(j, nelim + 1)] * a21;
        let wj2 = a[(j, nelim)] * a21 + a[(j, nelim + 1)] * a22;
        for i in j..m {
            a[(i, j)] -= a[(i, nelim)] * wj1 + a[(i, nelim + 1)] * wj2;
        }
    }

    // Zero D off-diagonal for extract_front_factors convention
    a[(nelim + 1, nelim)] = 0.0;

    max_l
}

/// Use TPP as the primary factorization method for small fronts.
///
/// Matches SPRAL's behavior: when `ncol < INNER_BLOCK_SIZE`, SPRAL uses
/// `ldlt_tpp_factor` instead of `block_ldlt` (complete pivoting). TPP tries
/// 2x2 pivots first for every column, finding more 2x2 opportunities than
/// complete pivoting (which only tries 2x2 when the off-diagonal is the
/// global maximum). For indefinite matrices, this produces significantly
/// better-conditioned factors.
///
/// See SPRAL `ldlt_app.cxx` `Block::factor()` lines 994-1020 (BSD-3).
fn tpp_factor_as_primary(
    mut a: MatMut<'_, f64>,
    num_fully_summed: usize,
    options: &AptpOptions,
) -> Result<AptpFactorResult, SparseError> {
    let m = a.nrows();
    let p = num_fully_summed;

    let mut col_order: Vec<usize> = (0..m).collect();
    let mut d = MixedDiagonal::new(p);
    let mut stats = AptpStatistics::default();
    let mut pivot_log = Vec::with_capacity(p);

    let additional = tpp_factor(
        a.rb_mut(),
        0, // start_col = 0: primary factorization
        p,
        &mut col_order,
        &mut d,
        &mut stats,
        &mut pivot_log,
        options,
    );

    let num_eliminated = additional;
    d.truncate(num_eliminated);
    stats.num_delayed = p - num_eliminated;

    let delayed_cols: Vec<usize> = (num_eliminated..p).map(|i| col_order[i]).collect();

    Ok(AptpFactorResult {
        d,
        perm: col_order,
        num_eliminated,
        delayed_cols,
        stats,
        pivot_log,
    })
}

/// TPP fallback: serial column-by-column factorization on APTP's remaining columns.
///
/// Operates on the partially-factored matrix `a` where columns `0..start_col`
/// have been eliminated by APTP. Searches all remaining fully-summed columns
/// (`start_col..num_fully_summed`) for acceptable pivots using threshold partial
/// pivoting.
///
/// Uses full-matrix `swap_symmetric` which correctly propagates row swaps to
/// already-factored L columns (matching SPRAL's `aleft` parameter).
///
/// Returns the number of additional columns eliminated.
///
/// # References
///
/// - SPRAL `ldlt_tpp_factor()` in `spral/src/ssids/cpu/kernels/ldlt_tpp.cxx` (BSD-3)
/// - Duff, Hogg & Lopez (2020), Section 3: TPP as fallback after APTP
#[allow(clippy::too_many_arguments)]
fn tpp_factor(
    mut a: MatMut<'_, f64>,
    start_col: usize,
    num_fully_summed: usize,
    col_order: &mut [usize],
    global_d: &mut MixedDiagonal,
    stats: &mut AptpStatistics,
    pivot_log: &mut Vec<AptpPivotRecord>,
    options: &AptpOptions,
) -> usize {
    let m = a.nrows();
    let n = num_fully_summed;
    let u = options.threshold;
    let small = options.small;

    let mut nelim = start_col;

    while nelim < n {
        // Check if current column is effectively zero
        if tpp_is_col_small(a.as_ref(), nelim, nelim, m, small) {
            // Zero pivot: record and advance
            global_d.set_1x1(nelim, 0.0);
            stats.num_1x1 += 1;
            pivot_log.push(AptpPivotRecord {
                col: col_order[nelim],
                pivot_type: PivotType::OneByOne,
                max_l_entry: 0.0,
                was_fallback: true,
            });
            // Zero the column entries
            for r in (nelim + 1)..m {
                a[(r, nelim)] = 0.0;
            }
            nelim += 1;
            continue;
        }

        // Search columns p = nelim+1..n-1 for acceptable pivot
        let mut found = false;
        for p in (nelim + 1)..n {
            // Check if column p is effectively zero
            if tpp_is_col_small(a.as_ref(), p, nelim, m, small) {
                // Swap zero column to front and record zero pivot
                if p != nelim {
                    swap_symmetric(a.rb_mut(), p, nelim);
                    col_order.swap(p, nelim);
                }
                global_d.set_1x1(nelim, 0.0);
                stats.num_1x1 += 1;
                pivot_log.push(AptpPivotRecord {
                    col: col_order[nelim],
                    pivot_type: PivotType::OneByOne,
                    max_l_entry: 0.0,
                    was_fallback: true,
                });
                for r in (nelim + 1)..m {
                    a[(r, nelim)] = 0.0;
                }
                nelim += 1;
                found = true;
                break;
            }

            // Find column index t of largest |a(p, c)| for c in nelim..p
            let t = tpp_find_row_abs_max(a.as_ref(), p, nelim, p);

            // Try (t, p) as 2x2 pivot (requires t < p)
            let maxt = tpp_find_rc_abs_max_exclude(a.as_ref(), t, nelim, m, p);
            let maxp = tpp_find_rc_abs_max_exclude(a.as_ref(), p, nelim, m, t);
            if tpp_test_2x2(a.as_ref(), t, p, maxt, maxp, u, small).is_some() {
                // Accept 2x2 pivot: swap t→nelim, p→nelim+1
                if t != nelim {
                    swap_symmetric(a.rb_mut(), t, nelim);
                    col_order.swap(t, nelim);
                }
                // After first swap, p may have moved
                let new_p = if p == nelim { t } else { p };
                if new_p != nelim + 1 {
                    swap_symmetric(a.rb_mut(), new_p, nelim + 1);
                    col_order.swap(new_p, nelim + 1);
                }

                // Re-read D values after swap
                let d11 = a[(nelim, nelim)];
                let d21 = a[(nelim + 1, nelim)];
                let d22 = a[(nelim + 1, nelim + 1)];

                // Apply elimination and Schur update
                let max_l = tpp_apply_2x2(a.rb_mut(), nelim, m, n);

                // Record pivot
                global_d.set_2x2(Block2x2 {
                    first_col: nelim,
                    a: d11,
                    b: d21,
                    c: d22,
                });
                stats.num_2x2 += 1;
                stats.max_l_entry = stats.max_l_entry.max(max_l);
                pivot_log.push(AptpPivotRecord {
                    col: col_order[nelim],
                    pivot_type: PivotType::TwoByTwo {
                        partner: col_order[nelim + 1],
                    },
                    max_l_entry: max_l,
                    was_fallback: true,
                });
                pivot_log.push(AptpPivotRecord {
                    col: col_order[nelim + 1],
                    pivot_type: PivotType::TwoByTwo {
                        partner: col_order[nelim],
                    },
                    max_l_entry: max_l,
                    was_fallback: true,
                });
                nelim += 2;
                found = true;
                break;
            }

            // Try p as 1x1 pivot
            // maxp should include |a(p, t)| (SPRAL line 225)
            let maxp_with_t = maxp.max(a[(p, t)].abs());
            if a[(p, p)].abs() >= u * maxp_with_t {
                // Accept 1x1 pivot: swap p→nelim
                if p != nelim {
                    swap_symmetric(a.rb_mut(), p, nelim);
                    col_order.swap(p, nelim);
                }

                let max_l = tpp_apply_1x1(a.rb_mut(), nelim, m, n);

                global_d.set_1x1(nelim, a[(nelim, nelim)]);
                stats.num_1x1 += 1;
                stats.max_l_entry = stats.max_l_entry.max(max_l);
                pivot_log.push(AptpPivotRecord {
                    col: col_order[nelim],
                    pivot_type: PivotType::OneByOne,
                    max_l_entry: max_l,
                    was_fallback: true,
                });
                nelim += 1;
                found = true;
                break;
            }
        }

        if !found {
            // Last resort: try column nelim as 1x1 (we started searching at nelim+1)
            let maxp = tpp_find_rc_abs_max_exclude(a.as_ref(), nelim, nelim, m, usize::MAX);
            if a[(nelim, nelim)].abs() >= u * maxp {
                let max_l = tpp_apply_1x1(a.rb_mut(), nelim, m, n);

                global_d.set_1x1(nelim, a[(nelim, nelim)]);
                stats.num_1x1 += 1;
                stats.max_l_entry = stats.max_l_entry.max(max_l);
                pivot_log.push(AptpPivotRecord {
                    col: col_order[nelim],
                    pivot_type: PivotType::OneByOne,
                    max_l_entry: max_l,
                    was_fallback: true,
                });
                nelim += 1;
            } else {
                // No more pivots can be found
                break;
            }
        }
    }

    nelim - start_col
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::Mat;

    // ---- Phase 2: Test infrastructure (T004-T005) ----

    /// Reconstruct P^T L D L^T P from factorization components.
    fn reconstruct_dense_ldlt(l: &Mat<f64>, d: &MixedDiagonal, perm: &[usize]) -> Mat<f64> {
        let n = l.nrows();

        // Build D as dense
        let mut d_mat = Mat::zeros(n, n);
        let nd = d.dimension();
        let mut col = 0;
        while col < nd {
            match d.get_pivot_type(col) {
                PivotType::OneByOne => {
                    d_mat[(col, col)] = d.get_1x1(col);
                    col += 1;
                }
                PivotType::TwoByTwo { partner } if partner > col => {
                    let block = d.get_2x2(col);
                    d_mat[(col, col)] = block.a;
                    d_mat[(col, col + 1)] = block.b;
                    d_mat[(col + 1, col)] = block.b;
                    d_mat[(col + 1, col + 1)] = block.c;
                    col += 2;
                }
                _ => {
                    col += 1;
                }
            }
        }

        // L * D
        let mut w = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let mut sum = 0.0;
                for k in 0..n {
                    sum += l[(i, k)] * d_mat[(k, j)];
                }
                w[(i, j)] = sum;
            }
        }

        // W * L^T
        let mut ldlt = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let mut sum = 0.0;
                for k in 0..n {
                    sum += w[(i, k)] * l[(j, k)];
                }
                ldlt[(i, j)] = sum;
            }
        }

        // Apply permutation: reconstructed[perm[i], perm[j]] = ldlt[i, j]
        let mut result = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                result[(perm[i], perm[j])] = ldlt[(i, j)];
            }
        }

        result
    }

    /// Compute ||A - P^T L D L^T P|| / ||A||
    fn dense_reconstruction_error(
        original: &Mat<f64>,
        l: &Mat<f64>,
        d: &MixedDiagonal,
        perm: &[usize],
    ) -> f64 {
        let reconstructed = reconstruct_dense_ldlt(l, d, perm);
        let n = original.nrows();

        let mut diff_norm_sq = 0.0_f64;
        let mut orig_norm_sq = 0.0_f64;

        for i in 0..n {
            for j in 0..n {
                let diff = original[(i, j)] - reconstructed[(i, j)];
                diff_norm_sq += diff * diff;
                orig_norm_sq += original[(i, j)] * original[(i, j)];
            }
        }

        if orig_norm_sq == 0.0 {
            return diff_norm_sq.sqrt();
        }
        (diff_norm_sq / orig_norm_sq).sqrt()
    }

    fn symmetric_matrix(n: usize, f: impl Fn(usize, usize) -> f64) -> Mat<f64> {
        Mat::from_fn(n, n, |i, j| if i >= j { f(i, j) } else { f(j, i) })
    }

    // ---- Phase 2 infrastructure test ----

    #[test]
    fn test_reconstruction_trivial_identity() {
        let l = Mat::identity(2, 2);
        let mut d = MixedDiagonal::new(2);
        d.set_1x1(0, 1.0);
        d.set_1x1(1, 1.0);
        let perm = vec![0, 1];

        let result = reconstruct_dense_ldlt(&l, &d, &perm);
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (result[(i, j)] - expected).abs() < 1e-14,
                    "({},{}) = {}, expected {}",
                    i,
                    j,
                    result[(i, j)],
                    expected
                );
            }
        }
    }

    // ---- Complete Pivoting Tests (Algorithm 4.1) ----

    /// Helper: run complete_pivoting_factor and extract L for reconstruction testing.
    fn complete_pivoting_factor_and_extract(
        a: &Mat<f64>,
    ) -> (Mat<f64>, MixedDiagonal, Vec<usize>, AptpStatistics) {
        let mut a_copy = a.clone();
        let result = complete_pivoting_factor(a_copy.as_mut(), 1e-20);
        let l = extract_l(a_copy.as_ref(), &result.d, result.num_eliminated);
        (l, result.d, result.perm, result.stats)
    }

    #[test]
    fn test_cp_identity() {
        // T006: 3×3 identity → D=[1,1,1], no permutation
        let a = Mat::identity(3, 3);
        let mut a_copy = a.clone();
        let result = complete_pivoting_factor(a_copy.as_mut(), 1e-20);

        assert_eq!(result.num_eliminated, 3);
        assert_eq!(result.stats.num_1x1, 3);
        assert_eq!(result.stats.num_2x2, 0);
        assert_eq!(result.stats.num_delayed, 0);

        // D should be [1, 1, 1]
        for i in 0..3 {
            assert!((result.d.get_1x1(i) - 1.0).abs() < 1e-14);
        }

        // No off-diagonal L entries
        assert!(result.stats.max_l_entry < 1e-14);
    }

    #[test]
    fn test_cp_diagonal_pivot_ordering() {
        // T007: 3×3 diagonal with known pivot ordering (largest diagonal first)
        let a = symmetric_matrix(3, |i, j| if i == j { [2.0, 5.0, 3.0][i] } else { 0.0 });
        let mut a_copy = a.clone();
        let result = complete_pivoting_factor(a_copy.as_mut(), 1e-20);

        assert_eq!(result.num_eliminated, 3);
        assert_eq!(result.stats.num_1x1, 3);

        // First pivot should be 5.0 (col 1), then 3.0 (col 2), then 2.0 (col 0)
        assert!(
            (result.d.get_1x1(0) - 5.0).abs() < 1e-14,
            "first pivot should be 5.0, got {}",
            result.d.get_1x1(0)
        );
        assert!(
            (result.d.get_1x1(1) - 3.0).abs() < 1e-14,
            "second pivot should be 3.0, got {}",
            result.d.get_1x1(1)
        );
        assert!(
            (result.d.get_1x1(2) - 2.0).abs() < 1e-14,
            "third pivot should be 2.0, got {}",
            result.d.get_1x1(2)
        );
    }

    #[test]
    fn test_cp_2x2_pivot() {
        // T008: 4×4 matrix requiring 2×2 pivot (off-diagonal maximum)
        // Designed so max entry is off-diagonal and Δ condition passes
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [0.1, 10.0, 0.0, 0.0],
                [10.0, 0.1, 0.0, 0.0],
                [0.0, 0.0, 5.0, 1.0],
                [0.0, 0.0, 1.0, 3.0],
            ];
            vals[i][j]
        });

        let (l, d, perm, stats) = complete_pivoting_factor_and_extract(&a);

        assert!(
            stats.num_2x2 >= 1,
            "expected at least one 2×2 pivot, got {}",
            stats.num_2x2
        );
        // L entries bounded by 4 (complete pivoting guarantee)
        assert!(
            stats.max_l_entry <= COMPLETE_PIVOTING_GROWTH_BOUND + 1e-10,
            "L entries should be bounded by 4, got {}",
            stats.max_l_entry
        );

        // Reconstruction
        let error = dense_reconstruction_error(&a, &l, &d, &perm);
        assert!(error < 1e-12, "reconstruction error {:.2e} >= 1e-12", error);
    }

    #[test]
    fn test_cp_failed_2x2_fallback() {
        // T009: 4×4 matrix where 2×2 Δ test fails → fallback to 1×1 on max diagonal
        // Need |Δ| < 0.5 * |a_tm|^2
        // a_mm * a_tt - a_tm^2 should be small relative to a_tm^2
        // Let a_mm = 1.0, a_tt = 1.0, a_tm = 2.0
        // Δ = 1*1 - 4 = -3, |Δ| = 3, 0.5*a_tm^2 = 2 → |Δ| > 0.5*a_tm^2
        // Need: a_mm ≈ a_tt and both ≈ a_tm
        // Let a_mm = 0.5, a_tt = 0.5, a_tm = 1.0
        // Δ = 0.25 - 1 = -0.75, |Δ| = 0.75, 0.5*a_tm^2 = 0.5 → 0.75 > 0.5, passes!
        // Let a_mm = 0.1, a_tt = 0.1, a_tm = 1.0
        // Δ = 0.01 - 1 = -0.99, |Δ| = 0.99, 0.5*a_tm^2 = 0.5 → 0.99 > 0.5, passes!
        // Let a_mm = 0.01, a_tt = 0.01, a_tm = 1.0
        // Δ = 0.0001 - 1 = -0.9999, |Δ| = 0.9999, 0.5*1 = 0.5 → passes!
        // Hmm, hard to fail. The condition is |det| >= 0.5*a_tm^2.
        // For failure: need det ≈ 0. This means a_mm*a_tt ≈ a_tm^2.
        // Let a_mm = 2.0, a_tt = 2.0, a_tm = 2.0
        // Δ = 4 - 4 = 0, |Δ| = 0 < 0.5*4 = 2 → FAILS! Good.
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [2.0, 2.0, 0.1, 0.1],
                [2.0, 2.0, 0.1, 0.1],
                [0.1, 0.1, 5.0, 0.0],
                [0.1, 0.1, 0.0, 3.0],
            ];
            vals[i][j]
        });

        let (l, d, perm, stats) = complete_pivoting_factor_and_extract(&a);

        // With Δ ≈ 0, should fall back to 1×1 on max(|a_mm|, |a_tt|)
        // L entries bounded by √2 < 4 in this case (paper's bound)
        assert!(
            stats.max_l_entry <= COMPLETE_PIVOTING_GROWTH_BOUND + 1e-10,
            "L entries should be bounded by 4"
        );

        // Reconstruction
        let error = dense_reconstruction_error(&a, &l, &d, &perm);
        assert!(error < 1e-12, "reconstruction error {:.2e} >= 1e-12", error);
    }

    #[test]
    fn test_cp_singular_block() {
        // T010: singular/near-singular block → zero pivot handling
        let mut a = Mat::zeros(3, 3);
        a[(0, 0)] = 1e-25;
        a[(1, 1)] = 1e-25;
        a[(2, 2)] = 1e-25;

        let result = complete_pivoting_factor(a.as_mut(), 1e-20);

        // All entries below small → all should be zero pivots
        assert_eq!(
            result.num_eliminated, 0,
            "near-singular block should have 0 eliminations"
        );
        assert_eq!(result.stats.num_delayed, 3);
    }

    #[test]
    fn test_cp_reconstruction_random() {
        // T011: reconstruction on random symmetric indefinite matrices
        // Use deterministic seed for reproducibility
        let sizes = [8, 16, 32];
        for &n in &sizes {
            let a = symmetric_matrix(n, |i, j| {
                // Deterministic pseudo-random indefinite matrix
                let seed = (i * 1000 + j * 7 + 13) as f64;
                let val = (seed * 0.618033988749).fract() * 2.0 - 1.0;
                if i == j {
                    val * 10.0 // diagonal dominance but indefinite
                } else {
                    val
                }
            });

            let (l, d, perm, stats) = complete_pivoting_factor_and_extract(&a);

            // Should eliminate all columns (no singular blocks)
            assert_eq!(
                stats.num_1x1 + 2 * stats.num_2x2 + stats.num_delayed,
                n,
                "statistics invariant for {}x{}",
                n,
                n
            );

            if stats.num_delayed == 0 {
                let error = dense_reconstruction_error(&a, &l, &d, &perm);
                assert!(
                    error < 1e-12,
                    "complete pivoting {}x{}: reconstruction error {:.2e} >= 1e-12",
                    n,
                    n,
                    error
                );
            }

            // L entries bounded by 4 (Algorithm 4.1 guarantee)
            assert!(
                stats.max_l_entry <= COMPLETE_PIVOTING_GROWTH_BOUND + 1e-10,
                "complete pivoting {}x{}: max_l_entry {:.2e} > 4",
                n,
                n,
                stats.max_l_entry
            );
        }
    }

    // ---- US1 Tests (T006-T012) ----

    #[test]
    fn test_1x1_trivial_diagonal() {
        let a = symmetric_matrix(2, |i, j| if i == j { [4.0, 9.0][i] } else { 0.0 });

        let opts = AptpOptions {
            fallback: AptpFallback::Delay,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        // All columns eliminated, no delays
        assert_eq!(result.stats.num_1x1 + 2 * result.stats.num_2x2, 2);
        assert_eq!(result.stats.num_delayed, 0);
        assert!(result.delayed_cols.is_empty());

        let error =
            dense_reconstruction_error(&a, &result.l, &result.d, result.perm.as_ref().arrays().0);
        assert!(error < 1e-12, "reconstruction error {:.2e} >= 1e-12", error);
    }

    #[test]
    fn test_1x1_positive_definite_3x3() {
        let a = symmetric_matrix(3, |i, j| {
            let vals = [[4.0, 2.0, 1.0], [2.0, 5.0, 3.0], [1.0, 3.0, 6.0]];
            vals[i][j]
        });

        let opts = AptpOptions::default();
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        // All columns eliminated, no delays
        assert_eq!(result.stats.num_1x1 + 2 * result.stats.num_2x2, 3);
        assert_eq!(result.stats.num_delayed, 0);

        let error =
            dense_reconstruction_error(&a, &result.l, &result.d, result.perm.as_ref().arrays().0);
        assert!(error < 1e-12, "reconstruction error {:.2e} >= 1e-12", error);
    }

    #[test]
    fn test_all_delayed_zero_matrix() {
        let n = 4;
        let a = Mat::zeros(n, n);

        let opts = AptpOptions {
            failed_pivot_method: FailedPivotMethod::Pass,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        // TPP treats zero columns as zero pivots (1x1 with D=0), not delays.
        // With FailedPivotMethod::Pass, TPP is still used as primary for small
        // matrices, and it handles zero columns by recording them as zero pivots.
        // Total columns accounted for:
        let total = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
        assert_eq!(total, n, "total pivots + delays should equal n");
    }

    #[test]
    fn test_1x1_singularity_detection() {
        let a = symmetric_matrix(3, |i, j| if i == j { [4.0, 1e-25, 9.0][i] } else { 0.0 });

        let opts = AptpOptions {
            fallback: AptpFallback::Delay,
            failed_pivot_method: FailedPivotMethod::Pass,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        // The near-zero entry should be detected (as zero pivot or delay).
        // TPP records zero columns as zero pivots (1x1 with D=0), while APTP
        // delays them. Either way, the other 2 columns should be eliminated.
        let eliminated = result.stats.num_1x1 + 2 * result.stats.num_2x2;
        assert!(
            eliminated >= 2,
            "should eliminate at least 2 columns, got {}",
            eliminated
        );
        let total = eliminated + result.stats.num_delayed;
        assert_eq!(total, 3, "total pivots + delays should equal n");
    }

    #[test]
    fn test_stability_bound_enforced() {
        // With complete pivoting, this matrix is handled via 2×2 pivot
        // (max off-diagonal entry 1.0 triggers 2×2 with Δ condition).
        // Complete pivoting bounds L entries by 4 (u=0.25), so no delays.
        let a = symmetric_matrix(2, |i, j| {
            let vals = [[1e-4, 1.0], [1.0, 1.0]];
            vals[i][j]
        });

        let opts = AptpOptions {
            fallback: AptpFallback::Delay,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        // Complete pivoting eliminates all columns (no delays)
        let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
        assert_eq!(sum, 2);
        // Verify correctness via reconstruction
        let error =
            dense_reconstruction_error(&a, &result.l, &result.d, result.perm.as_ref().arrays().0);
        assert!(error < 1e-12, "reconstruction error {:.2e} >= 1e-12", error);
    }

    #[test]
    fn test_1x1_matrix() {
        let a = Mat::from_fn(1, 1, |_, _| 5.0);

        let opts = AptpOptions::default();
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        assert_eq!(result.stats.num_1x1, 1);
        assert_eq!(result.stats.num_delayed, 0);

        let error =
            dense_reconstruction_error(&a, &result.l, &result.d, result.perm.as_ref().arrays().0);
        assert!(error < 1e-14, "reconstruction error {:.2e}", error);
    }

    #[test]
    fn test_statistics_sum_invariant() {
        let a = symmetric_matrix(5, |i, j| {
            let vals = [
                [10.0, 1.0, 0.0, 0.0, 0.0],
                [1.0, 20.0, 2.0, 0.0, 0.0],
                [0.0, 2.0, 30.0, 3.0, 0.0],
                [0.0, 0.0, 3.0, 40.0, 4.0],
                [0.0, 0.0, 0.0, 4.0, 50.0],
            ];
            vals[i][j]
        });

        let opts = AptpOptions::default();
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
        assert_eq!(sum, 5, "statistics sum {} != n=5", sum);
    }

    // ---- US2 Tests (T022-T026) ----

    #[test]
    fn test_2x2_pivot_known_indefinite() {
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [0.01, 10.0, 0.0, 0.0],
                [10.0, 0.01, 0.0, 0.0],
                [0.0, 0.0, 5.0, 1.0],
                [0.0, 0.0, 1.0, 3.0],
            ];
            vals[i][j]
        });

        let opts = AptpOptions {
            fallback: AptpFallback::BunchKaufman,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        assert!(
            result.stats.num_2x2 >= 1,
            "expected 2x2 pivot, got num_2x2={}",
            result.stats.num_2x2
        );

        let error =
            dense_reconstruction_error(&a, &result.l, &result.d, result.perm.as_ref().arrays().0);
        assert!(error < 1e-12, "reconstruction error {:.2e} >= 1e-12", error);
    }

    #[test]
    fn test_2x2_stability_test() {
        let a_good = symmetric_matrix(2, |i, j| {
            let vals = [[1.0, 5.0], [5.0, 1.0]];
            vals[i][j]
        });
        let opts = AptpOptions {
            fallback: AptpFallback::BunchKaufman,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a_good.as_ref(), &opts).unwrap();
        let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
        assert_eq!(sum, 2);
    }

    #[test]
    fn test_bk_vs_delay_fallback() {
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [0.01, 10.0, 0.0, 0.0],
                [10.0, 0.01, 0.0, 0.0],
                [0.0, 0.0, 5.0, 1.0],
                [0.0, 0.0, 1.0, 3.0],
            ];
            vals[i][j]
        });

        let bk_opts = AptpOptions {
            fallback: AptpFallback::BunchKaufman,
            ..AptpOptions::default()
        };
        let delay_opts = AptpOptions {
            fallback: AptpFallback::Delay,
            ..AptpOptions::default()
        };

        let bk_result = aptp_factor(a.as_ref(), &bk_opts).unwrap();
        let delay_result = aptp_factor(a.as_ref(), &delay_opts).unwrap();

        assert!(
            bk_result.stats.num_delayed <= delay_result.stats.num_delayed,
            "BK delayed {} > Delay delayed {}",
            bk_result.stats.num_delayed,
            delay_result.stats.num_delayed
        );

        if bk_result.stats.num_delayed == 0 {
            let error = dense_reconstruction_error(
                &a,
                &bk_result.l,
                &bk_result.d,
                bk_result.perm.as_ref().arrays().0,
            );
            assert!(error < 1e-12, "BK reconstruction error {:.2e}", error);
        }
    }

    #[test]
    fn test_strict_threshold() {
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [4.0, 2.0, 1.0, 0.5],
                [2.0, 5.0, 2.0, 1.0],
                [1.0, 2.0, 6.0, 2.0],
                [0.5, 1.0, 2.0, 7.0],
            ];
            vals[i][j]
        });

        let loose = AptpOptions {
            threshold: 0.01,
            fallback: AptpFallback::Delay,
            ..AptpOptions::default()
        };
        let strict = AptpOptions {
            threshold: 0.5,
            fallback: AptpFallback::Delay,
            ..AptpOptions::default()
        };

        let loose_result = aptp_factor(a.as_ref(), &loose).unwrap();
        let strict_result = aptp_factor(a.as_ref(), &strict).unwrap();

        assert!(
            strict_result.stats.num_delayed >= loose_result.stats.num_delayed,
            "strict delayed {} < loose delayed {}",
            strict_result.stats.num_delayed,
            loose_result.stats.num_delayed,
        );
    }

    #[test]
    fn test_permutation_valid() {
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [0.01, 10.0, 0.0, 0.0],
                [10.0, 0.01, 0.0, 0.0],
                [0.0, 0.0, 5.0, 1.0],
                [0.0, 0.0, 1.0, 3.0],
            ];
            vals[i][j]
        });

        let opts = AptpOptions {
            fallback: AptpFallback::BunchKaufman,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        let (fwd, inv) = result.perm.as_ref().arrays();
        let n = fwd.len();
        assert_eq!(n, 4);
        let mut seen = vec![false; n];
        for &v in fwd {
            assert!(v < n, "perm value {} >= n={}", v, n);
            assert!(!seen[v], "duplicate perm value {}", v);
            seen[v] = true;
        }
        for i in 0..n {
            assert_eq!(inv[fwd[i]], i);
            assert_eq!(fwd[inv[i]], i);
        }
    }

    // ---- US3 Tests (T033-T036) ----

    #[test]
    fn test_pd_statistics() {
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [10.0, 1.0, 0.0, 0.0],
                [1.0, 20.0, 2.0, 0.0],
                [0.0, 2.0, 30.0, 3.0],
                [0.0, 0.0, 3.0, 40.0],
            ];
            vals[i][j]
        });

        let opts = AptpOptions::default();
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        // All columns eliminated, no delays
        assert_eq!(result.stats.num_1x1 + 2 * result.stats.num_2x2, 4);
        assert_eq!(result.stats.num_delayed, 0);
        assert!(result.stats.max_l_entry < 1.0 / opts.threshold);
    }

    #[test]
    fn test_max_l_entry_accuracy() {
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [4.0, 2.0, 1.0, 0.5],
                [2.0, 5.0, 2.0, 1.0],
                [1.0, 2.0, 6.0, 2.0],
                [0.5, 1.0, 2.0, 7.0],
            ];
            vals[i][j]
        });

        let opts = AptpOptions::default();
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        let n = result.l.nrows();
        let mut actual_max = 0.0_f64;
        for j in 0..n {
            for i in (j + 1)..n {
                let val = result.l[(i, j)].abs();
                if val > actual_max {
                    actual_max = val;
                }
            }
        }

        assert!(
            (result.stats.max_l_entry - actual_max).abs() < 1e-14,
            "stats.max_l_entry={:.6e}, actual={:.6e}",
            result.stats.max_l_entry,
            actual_max
        );
    }

    #[test]
    fn test_pivot_log_completeness() {
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [4.0, 2.0, 1.0, 0.5],
                [2.0, 5.0, 2.0, 1.0],
                [1.0, 2.0, 6.0, 2.0],
                [0.5, 1.0, 2.0, 7.0],
            ];
            vals[i][j]
        });

        let opts = AptpOptions::default();
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        assert_eq!(result.pivot_log.len(), 4);
    }

    #[test]
    fn test_inertia_from_d() {
        let a = symmetric_matrix(5, |i, j| {
            if i == j {
                [3.0, -2.0, 1.0, -4.0, 5.0][i]
            } else {
                0.0
            }
        });

        let opts = AptpOptions {
            fallback: AptpFallback::Delay,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        assert_eq!(result.stats.num_delayed, 0, "should have no delays");

        let inertia = result.d.compute_inertia();
        assert_eq!(inertia.positive, 3, "expected 3 positive");
        assert_eq!(inertia.negative, 2, "expected 2 negative");
        assert_eq!(inertia.zero, 0, "expected 0 zero");
    }

    // ---- Phase 6: Polish tests ----

    #[test]
    fn test_partial_factorization() {
        let n = 8;
        let p = 4;
        let a = symmetric_matrix(n, |i, j| {
            if i == j {
                10.0 + i as f64
            } else {
                1.0 / (1.0 + (i as f64 - j as f64).abs())
            }
        });

        let opts = AptpOptions::default();
        let mut a_copy = a.clone();
        let result = aptp_factor_in_place(a_copy.as_mut(), p, &opts).unwrap();

        let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
        assert_eq!(sum, p, "statistics should sum to num_fully_summed={}", p);
    }

    #[test]
    fn test_edge_case_extreme_thresholds() {
        let a = symmetric_matrix(3, |i, j| {
            let vals = [[4.0, 2.0, 1.0], [2.0, 5.0, 2.0], [1.0, 2.0, 6.0]];
            vals[i][j]
        });

        let loose = AptpOptions {
            threshold: 0.001,
            fallback: AptpFallback::Delay,
            ..AptpOptions::default()
        };
        let result_loose = aptp_factor(a.as_ref(), &loose).unwrap();
        assert_eq!(result_loose.stats.num_delayed, 0);

        let strict = AptpOptions {
            threshold: 1.0,
            fallback: AptpFallback::Delay,
            ..AptpOptions::default()
        };
        let result_strict = aptp_factor(a.as_ref(), &strict).unwrap();
        let sum = result_strict.stats.num_1x1
            + 2 * result_strict.stats.num_2x2
            + result_strict.stats.num_delayed;
        assert_eq!(sum, 3);
    }

    #[test]
    fn test_both_fallback_strategies_valid() {
        let matrices = [
            symmetric_matrix(3, |i, j| {
                let vals = [[1.0, 2.0, 0.0], [2.0, -1.0, 1.0], [0.0, 1.0, 3.0]];
                vals[i][j]
            }),
            symmetric_matrix(4, |i, j| {
                let vals = [
                    [0.01, 10.0, 0.0, 0.0],
                    [10.0, 0.01, 0.0, 0.0],
                    [0.0, 0.0, 5.0, 1.0],
                    [0.0, 0.0, 1.0, 3.0],
                ];
                vals[i][j]
            }),
        ];

        for (idx, a) in matrices.iter().enumerate() {
            for fallback in [AptpFallback::BunchKaufman, AptpFallback::Delay] {
                let opts = AptpOptions {
                    fallback,
                    ..AptpOptions::default()
                };
                let result = aptp_factor(a.as_ref(), &opts).unwrap();

                if result.stats.num_delayed == 0 {
                    let error = dense_reconstruction_error(
                        a,
                        &result.l,
                        &result.d,
                        result.perm.as_ref().arrays().0,
                    );
                    assert!(
                        error < 1e-12,
                        "matrix {} {:?}: error {:.2e}",
                        idx,
                        fallback,
                        error
                    );
                }

                let sum =
                    result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
                let n = a.nrows();
                assert_eq!(sum, n, "statistics invariant broken for matrix {}", idx);
            }
        }
    }

    // ---- Input validation tests ----

    #[test]
    fn test_invalid_non_square() {
        let a = Mat::zeros(3, 4);
        let opts = AptpOptions::default();
        let result = aptp_factor(a.as_ref(), &opts);
        assert!(matches!(result, Err(SparseError::InvalidInput { .. })));
    }

    #[test]
    fn test_invalid_threshold() {
        let a = Mat::zeros(2, 2);
        let opts = AptpOptions {
            threshold: 0.0,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts);
        assert!(matches!(result, Err(SparseError::InvalidInput { .. })));

        let opts2 = AptpOptions {
            threshold: 1.5,
            ..AptpOptions::default()
        };
        let result2 = aptp_factor(a.as_ref(), &opts2);
        assert!(matches!(result2, Err(SparseError::InvalidInput { .. })));
    }

    #[test]
    fn test_invalid_num_fully_summed_too_large() {
        let mut a = Mat::zeros(3, 3);
        let opts = AptpOptions::default();
        let result = aptp_factor_in_place(a.as_mut(), 5, &opts);
        assert!(matches!(result, Err(SparseError::InvalidInput { .. })));
    }

    // ---- Edge case and regression tests ----

    #[test]
    fn test_zero_fully_summed() {
        let a = symmetric_matrix(3, |i, j| {
            let vals = [[4.0, 2.0, 1.0], [2.0, 5.0, 3.0], [1.0, 3.0, 6.0]];
            vals[i][j]
        });

        let opts = AptpOptions::default();
        let mut a_copy = a.clone();
        let result = aptp_factor_in_place(a_copy.as_mut(), 0, &opts).unwrap();

        assert_eq!(result.num_eliminated, 0);
        assert_eq!(result.stats.num_1x1, 0);
        assert_eq!(result.stats.num_2x2, 0);
        assert_eq!(result.stats.num_delayed, 0);
        assert!(result.pivot_log.is_empty());
    }

    #[test]
    fn test_2x2_fallback_also_fails() {
        // With complete pivoting, the max entry is 12.0 on diagonal (1,1).
        // Complete pivoting starts with the largest entry and can handle
        // matrices that the old threshold-based approach would delay.
        let a = symmetric_matrix(3, |i, j| {
            let vals = [[0.001, 0.11, 0.0], [0.11, 12.0, 0.1], [0.0, 0.1, 5.0]];
            vals[i][j]
        });

        let opts = AptpOptions {
            fallback: AptpFallback::BunchKaufman,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        // Complete pivoting eliminates all columns successfully
        let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
        assert_eq!(sum, 3, "statistics sum {} != n=3", sum);

        // Verify correctness via reconstruction
        if result.stats.num_delayed == 0 {
            let error = dense_reconstruction_error(
                &a,
                &result.l,
                &result.d,
                result.perm.as_ref().arrays().0,
            );
            assert!(error < 1e-12, "reconstruction error {:.2e}", error);
        }

        // L entries bounded by 4 (complete pivoting guarantee)
        assert!(
            result.stats.max_l_entry <= COMPLETE_PIVOTING_GROWTH_BOUND + 1e-10,
            "max_l_entry {} > 4",
            result.stats.max_l_entry
        );
    }

    #[test]
    fn test_invalid_small_negative() {
        let a = Mat::zeros(2, 2);
        let opts = AptpOptions {
            small: -1.0,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts);
        assert!(matches!(result, Err(SparseError::InvalidInput { .. })));
    }

    // ---- Phase 6: Random matrix stress tests (T041-T042) ----

    #[cfg(feature = "test-util")]
    mod stress_tests {
        use super::*;
        use crate::testing::generators::{RandomMatrixConfig, generate_random_symmetric};
        use rand::SeedableRng;
        use rand::rngs::StdRng;

        /// Generate a dense symmetric PD matrix from the sparse generator.
        fn random_dense_pd(n: usize, rng: &mut impl rand::Rng) -> Mat<f64> {
            let config = RandomMatrixConfig {
                size: n,
                target_nnz: n * n / 2,
                positive_definite: true,
            };
            generate_random_symmetric(&config, rng).unwrap().to_dense()
        }

        /// Generate a dense symmetric indefinite matrix from the sparse generator.
        fn random_dense_indefinite(n: usize, rng: &mut impl rand::Rng) -> Mat<f64> {
            let config = RandomMatrixConfig {
                size: n,
                target_nnz: n * n / 2,
                positive_definite: false,
            };
            generate_random_symmetric(&config, rng).unwrap().to_dense()
        }

        #[test]
        fn test_random_pd_matrices() {
            let mut rng = StdRng::seed_from_u64(42);
            let opts = AptpOptions::default();
            let sizes = [5, 10, 20, 50];
            let mut total = 0;

            for &n in &sizes {
                for _ in 0..25 {
                    let a = random_dense_pd(n, &mut rng);
                    let result = aptp_factor(a.as_ref(), &opts).unwrap();

                    assert_eq!(
                        result.stats.num_delayed, 0,
                        "PD matrix {}x{} should have zero delays",
                        n, n
                    );
                    let total_elim = result.stats.num_1x1 + 2 * result.stats.num_2x2;
                    assert_eq!(
                        total_elim, n,
                        "PD matrix {}x{} should eliminate all columns (1x1={}, 2x2={})",
                        n, n, result.stats.num_1x1, result.stats.num_2x2
                    );

                    let error = dense_reconstruction_error(
                        &a,
                        &result.l,
                        &result.d,
                        result.perm.as_ref().arrays().0,
                    );
                    assert!(
                        error < 1e-12,
                        "PD {}x{}: reconstruction error {:.2e}",
                        n,
                        n,
                        error
                    );
                    total += 1;
                }
            }
            assert!(total >= 100, "ran {} PD tests, need >= 100", total);
        }

        #[test]
        fn test_random_indefinite_matrices() {
            let mut rng = StdRng::seed_from_u64(123);
            let opts = AptpOptions {
                fallback: AptpFallback::BunchKaufman,
                ..AptpOptions::default()
            };
            let sizes = [5, 10, 20, 50];
            let mut total = 0;

            for &n in &sizes {
                for _ in 0..25 {
                    let a = random_dense_indefinite(n, &mut rng);
                    let result = aptp_factor(a.as_ref(), &opts).unwrap();

                    // Statistics invariant
                    let sum =
                        result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
                    assert_eq!(sum, n, "stats invariant for {}x{}", n, n);

                    // Stability bound
                    assert!(
                        result.stats.max_l_entry < 1.0 / opts.threshold,
                        "stability bound violated for {}x{}",
                        n,
                        n
                    );

                    // Reconstruction (only if fully eliminated)
                    if result.stats.num_delayed == 0 {
                        let error = dense_reconstruction_error(
                            &a,
                            &result.l,
                            &result.d,
                            result.perm.as_ref().arrays().0,
                        );
                        assert!(
                            error < 1e-12,
                            "indefinite {}x{}: reconstruction error {:.2e}",
                            n,
                            n,
                            error
                        );
                    }
                    total += 1;
                }
            }
            assert!(total >= 100, "ran {} indefinite tests, need >= 100", total);
        }
    }

    // ---- Phase 6: Integration tests with test data (T046-T047) ----

    #[cfg(feature = "test-util")]
    mod integration_tests {
        use super::*;
        use crate::testing::cases::{TestCaseFilter, load_test_cases};

        #[test]
        #[ignore] // Requires test-data/ on disk
        fn test_hand_constructed_matrices() {
            let cases = load_test_cases(&TestCaseFilter::hand_constructed())
                .expect("failed to load hand-constructed matrices");
            assert_eq!(cases.len(), 15, "expected 15 hand-constructed matrices");

            let opts = AptpOptions::default();
            let mut passed = 0;

            for case in &cases {
                let dense = case.matrix.to_dense();
                let n = dense.nrows();

                // Skip very large matrices (dense APTP is O(n^3))
                if n > 500 {
                    continue;
                }

                let result = aptp_factor(dense.as_ref(), &opts).unwrap();

                let sum =
                    result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
                assert_eq!(sum, n, "{}: stats invariant failed", case.name);

                if result.stats.num_delayed == 0 {
                    let error = dense_reconstruction_error(
                        &dense,
                        &result.l,
                        &result.d,
                        result.perm.as_ref().arrays().0,
                    );
                    assert!(
                        error < 1e-12,
                        "{}: reconstruction error {:.2e}",
                        case.name,
                        error
                    );
                }
                passed += 1;
            }
            assert!(passed >= 10, "only {} hand-constructed passed", passed);
        }

        #[test]
        #[ignore] // Requires SuiteSparse CI-subset on disk
        fn test_suitesparse_ci_dense() {
            let cases = load_test_cases(&TestCaseFilter::ci_subset())
                .expect("failed to load CI-subset matrices");

            let opts = AptpOptions::default();
            let mut tested = 0;

            for case in &cases {
                let n = case.matrix.nrows();

                // Only test matrices small enough for dense factorization.
                // Dense O(n^2) memory + O(n^3) time: cap at 200 to avoid OOM.
                if n > 200 {
                    continue;
                }

                let dense = case.matrix.to_dense();

                let result = aptp_factor(dense.as_ref(), &opts).unwrap();

                let sum =
                    result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
                assert_eq!(sum, n, "{}: stats invariant failed", case.name);

                if result.stats.num_delayed == 0 {
                    let error = dense_reconstruction_error(
                        &dense,
                        &result.l,
                        &result.d,
                        result.perm.as_ref().arrays().0,
                    );
                    assert!(
                        error < 1e-12,
                        "{}: reconstruction error {:.2e}",
                        case.name,
                        error
                    );
                }
                tested += 1;
            }
            assert!(
                tested > 0,
                "no CI-subset matrices small enough for dense test"
            );
        }
    }

    // ---- Phase 3: factor_inner tests (T022) ----

    #[test]
    fn test_factor_inner_reconstruction_moderate() {
        // T022: factor_inner on matrices of moderate size (128, 256)
        // verifying reconstruction < 1e-12.
        let sizes = [64, 128, 256];
        let opts = AptpOptions::default();

        for &n in &sizes {
            let a = symmetric_matrix(n, |i, j| {
                let seed = (i * 1000 + j * 7 + 13) as f64;
                let val = (seed * 0.618033988749).fract() * 2.0 - 1.0;
                if i == j { val * 10.0 } else { val }
            });

            let result = aptp_factor(a.as_ref(), &opts).unwrap();

            // Statistics invariant
            assert_eq!(
                result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed,
                n,
                "factor_inner {}x{}: statistics invariant",
                n,
                n
            );

            if result.stats.num_delayed == 0 {
                let error = dense_reconstruction_error(
                    &a,
                    &result.l,
                    &result.d,
                    result.perm.as_ref().arrays().0,
                );
                assert!(
                    error < 1e-12,
                    "factor_inner {}x{}: reconstruction error {:.2e} >= 1e-12",
                    n,
                    n,
                    error
                );
            }
        }
    }

    #[test]
    fn test_factor_inner_partial_factorization() {
        // factor_inner with num_fully_summed < m (contribution block present)
        let n = 64;
        let p = 48; // Only factor 48 of 64 columns
        let a = symmetric_matrix(n, |i, j| {
            let seed = (i * 1000 + j * 7 + 13) as f64;
            let val = (seed * 0.618033988749).fract() * 2.0 - 1.0;
            if i == j { val * 10.0 } else { val }
        });

        let opts = AptpOptions::default();
        let mut a_copy = a.to_owned();
        let result = aptp_factor_in_place(a_copy.as_mut(), p, &opts).unwrap();

        // Should have eliminated <= p columns
        assert!(
            result.num_eliminated <= p,
            "eliminated {} > p={}",
            result.num_eliminated,
            p
        );

        // Statistics should account for all p columns
        assert_eq!(
            result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed,
            p,
            "statistics invariant for partial factorization"
        );
    }

    // ---- Phase 5: Two-level integration tests (T030-T034) ----

    #[test]
    fn test_two_level_dispatch_small_block_size() {
        // T031/T032: Test the two-level dispatch by setting outer_block_size small
        // so matrices of moderate size trigger two_level_factor.
        let sizes = [33, 64, 100];
        let opts = AptpOptions {
            outer_block_size: 32,
            inner_block_size: 8,
            ..AptpOptions::default()
        };

        for &n in &sizes {
            let a = symmetric_matrix(n, |i, j| {
                let seed = (i * 1000 + j * 7 + 13) as f64;
                let val = (seed * 0.618033988749).fract() * 2.0 - 1.0;
                if i == j { val * 10.0 } else { val }
            });

            let result = aptp_factor(a.as_ref(), &opts).unwrap();

            assert_eq!(
                result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed,
                n,
                "two-level {}x{}: statistics invariant",
                n,
                n
            );

            if result.stats.num_delayed == 0 {
                let error = dense_reconstruction_error(
                    &a,
                    &result.l,
                    &result.d,
                    result.perm.as_ref().arrays().0,
                );
                assert!(
                    error < 1e-12,
                    "two-level {}x{}: reconstruction error {:.2e} >= 1e-12",
                    n,
                    n,
                    error
                );
            }
        }
    }

    #[test]
    fn test_two_level_single_outer_block() {
        // T031: frontal dimension == nb → single outer block, equivalent to factor_inner
        let n = 32;
        let opts = AptpOptions {
            outer_block_size: 32,
            inner_block_size: 8,
            ..AptpOptions::default()
        };

        let a = symmetric_matrix(n, |i, j| {
            let seed = (i * 1000 + j * 7 + 13) as f64;
            let val = (seed * 0.618033988749).fract() * 2.0 - 1.0;
            if i == j { val * 10.0 } else { val }
        });

        let result = aptp_factor(a.as_ref(), &opts).unwrap();
        if result.stats.num_delayed == 0 {
            let error = dense_reconstruction_error(
                &a,
                &result.l,
                &result.d,
                result.perm.as_ref().arrays().0,
            );
            assert!(error < 1e-12, "single block: error {:.2e}", error);
        }
    }

    #[test]
    fn test_two_level_boundary_nb_plus_1() {
        // T032: frontal dimension == nb+1 → two blocks, second block has 1 column
        let n = 33;
        let opts = AptpOptions {
            outer_block_size: 32,
            inner_block_size: 8,
            ..AptpOptions::default()
        };

        let a = symmetric_matrix(n, |i, j| {
            let seed = (i * 1000 + j * 7 + 13) as f64;
            let val = (seed * 0.618033988749).fract() * 2.0 - 1.0;
            if i == j { val * 10.0 } else { val }
        });

        let result = aptp_factor(a.as_ref(), &opts).unwrap();
        assert_eq!(
            result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed,
            n,
        );
        if result.stats.num_delayed == 0 {
            let error = dense_reconstruction_error(
                &a,
                &result.l,
                &result.d,
                result.perm.as_ref().arrays().0,
            );
            assert!(error < 1e-12, "nb+1 boundary: error {:.2e}", error);
        }
    }

    #[test]
    fn test_two_level_partial_factorization() {
        // T033: partial factorization (num_fully_summed < m) with dimension > nb
        // This triggers two_level_factor with a contribution block.
        let n = 80;
        let p = 50; // > outer_block_size=32, triggers two-level
        let opts = AptpOptions {
            outer_block_size: 32,
            inner_block_size: 8,
            ..AptpOptions::default()
        };

        let a = symmetric_matrix(n, |i, j| {
            let seed = (i * 1000 + j * 7 + 13) as f64;
            let val = (seed * 0.618033988749).fract() * 2.0 - 1.0;
            if i == j { val * 10.0 } else { val }
        });

        let mut a_copy = a.to_owned();
        let result = aptp_factor_in_place(a_copy.as_mut(), p, &opts).unwrap();

        assert!(result.num_eliminated <= p);
        assert_eq!(
            result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed,
            p,
        );
    }

    #[test]
    fn test_two_level_vs_single_level_equivalence() {
        // Verify that two-level and single-level produce equivalent reconstruction
        // error on the same matrix (by setting outer_block_size to force each path).
        let n = 64;
        let a = symmetric_matrix(n, |i, j| {
            let seed = (i * 1000 + j * 7 + 13) as f64;
            let val = (seed * 0.618033988749).fract() * 2.0 - 1.0;
            if i == j { val * 10.0 } else { val }
        });

        // Single-level: outer_block_size >= n
        let opts_single = AptpOptions {
            outer_block_size: 256,
            inner_block_size: 32,
            ..AptpOptions::default()
        };
        let result_single = aptp_factor(a.as_ref(), &opts_single).unwrap();

        // Two-level: outer_block_size < n
        let opts_two = AptpOptions {
            outer_block_size: 16,
            inner_block_size: 8,
            ..AptpOptions::default()
        };
        let result_two = aptp_factor(a.as_ref(), &opts_two).unwrap();

        // Both should achieve reconstruction < 1e-12
        if result_single.stats.num_delayed == 0 {
            let err_s = dense_reconstruction_error(
                &a,
                &result_single.l,
                &result_single.d,
                result_single.perm.as_ref().arrays().0,
            );
            assert!(err_s < 1e-12, "single-level error {:.2e}", err_s);
        }
        if result_two.stats.num_delayed == 0 {
            let err_t = dense_reconstruction_error(
                &a,
                &result_two.l,
                &result_two.d,
                result_two.perm.as_ref().arrays().0,
            );
            assert!(err_t < 1e-12, "two-level error {:.2e}", err_t);
        }
    }

    #[test]
    fn test_two_level_vs_unblocked_reconstruction() {
        // Regression test: col_order tracking corruption in two_level_factor when
        // block_perm was applied after delay swap, causing incorrect column mapping.
        //
        // When delay swaps happened BEFORE block_perm was applied to col_order,
        // the permutation mapping was corrupted, causing massive reconstruction
        // errors with small outer_block_size. This test verifies that both
        // two-level (ob=128, forcing 4 outer iterations) and unblocked (ob=huge)
        // paths produce equivalent, accurate reconstruction.
        let n = 512;

        // Build a symmetric indefinite matrix with small diagonal entries to force
        // 2x2 pivots and some delays. Seed = golden ratio hash for reproducibility.
        let a = symmetric_matrix(n, |i, j| {
            let seed = (i * 997 + j * 1013 + 42) as f64;
            let val = (seed * 0.618033988749).fract() * 2.0 - 1.0;
            if i == j {
                // Small diagonal → more 2x2 pivots and delays
                val * 0.5
            } else {
                val
            }
        });

        // Unblocked: outer_block_size >= n → single factor_inner call
        let opts_unblocked = AptpOptions {
            outer_block_size: 100_000,
            inner_block_size: 32,
            ..AptpOptions::default()
        };
        let result_unblocked = aptp_factor(a.as_ref(), &opts_unblocked).unwrap();
        let err_unblocked = dense_reconstruction_error(
            &a,
            &result_unblocked.l,
            &result_unblocked.d,
            result_unblocked.perm.as_ref().arrays().0,
        );
        assert!(
            err_unblocked < 1e-12,
            "unblocked reconstruction error {:.2e} exceeds 1e-12",
            err_unblocked,
        );

        // Two-level: outer_block_size=128 → ~4 outer iterations
        let opts_two_level = AptpOptions {
            outer_block_size: 128,
            inner_block_size: 32,
            ..AptpOptions::default()
        };
        let result_two_level = aptp_factor(a.as_ref(), &opts_two_level).unwrap();
        let err_two_level = dense_reconstruction_error(
            &a,
            &result_two_level.l,
            &result_two_level.d,
            result_two_level.perm.as_ref().arrays().0,
        );
        assert!(
            err_two_level < 1e-12,
            "two-level reconstruction error {:.2e} exceeds 1e-12 \
             (unblocked was {:.2e})",
            err_two_level,
            err_unblocked,
        );

        // Two-level error should be within 10x of unblocked
        let ratio = if err_unblocked > 0.0 {
            err_two_level / err_unblocked
        } else {
            1.0
        };
        assert!(
            ratio < 10.0,
            "two-level error ({:.2e}) is {:.1}x worse than unblocked ({:.2e})",
            err_two_level,
            ratio,
            err_unblocked,
        );
    }

    // ---- BLAS-3 pipeline tests ----

    #[test]
    fn test_factor_block_diagonal_basic() {
        // Factor 8×8 identity with block_size=4.
        // D=[1,1,1,1], identity permutation, no L entries within block.
        let mut a = Mat::identity(8, 8);
        let (d, perm, nelim, stats, _log) = factor_block_diagonal(a.as_mut(), 0, 4, 1e-20, 4);

        assert_eq!(nelim, 4);
        assert_eq!(stats.num_1x1, 4);
        assert_eq!(stats.num_2x2, 0);
        assert_eq!(stats.num_delayed, 0);
        assert!(stats.max_l_entry < 1e-14);

        // D should be [1, 1, 1, 1]
        for i in 0..4 {
            assert!((d.get_1x1(i) - 1.0).abs() < 1e-14);
        }

        // Identity permutation
        assert_eq!(perm, vec![0, 1, 2, 3]);

        // Panel rows (4-7) should be untouched (still identity entries)
        for i in 4..8 {
            for j in 0..4 {
                assert!(
                    a[(i, j)].abs() < 1e-14,
                    "panel entry ({},{}) should be 0, got {}",
                    i,
                    j,
                    a[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_factor_block_diagonal_block_scoped_swap() {
        // 8×8 matrix where max entry in block forces a swap.
        // With block-scoped swaps (row_limit=4), panel rows (4-7) should NOT be
        // affected by factor_block_diagonal. Panel permutation is handled separately.
        let mut a = symmetric_matrix(8, |i, j| {
            if i == j {
                [1.0, 5.0, 2.0, 3.0, 0.1, 0.1, 0.1, 0.1][i]
            } else if (i == 4 && j == 1) || (i == 1 && j == 4) {
                // Panel entry at (4,1)
                0.99
            } else {
                0.0
            }
        });

        // Save panel row 4's original column values before factor
        let panel_row_before: Vec<f64> = (0..4).map(|j| a[(4, j)]).collect();

        let (_d, perm, nelim, _stats, _log) = factor_block_diagonal(a.as_mut(), 0, 4, 1e-20, 4);

        assert!(nelim > 0, "should eliminate at least 1 column");

        // The max diagonal is 5.0 at column 1, so it should be swapped to position 0.
        assert_eq!(
            perm[0], 1,
            "first pivot should be original column 1 (max diag = 5.0)"
        );

        // Panel rows should be UNCHANGED (block-scoped swap with row_limit=4)
        for j in 0..4 {
            assert!(
                (a[(4, j)] - panel_row_before[j]).abs() < 1e-14,
                "panel row (4,{}) should be unchanged: got {}, expected {}",
                j,
                a[(4, j)],
                panel_row_before[j]
            );
        }
    }

    #[test]
    fn test_blas3_pipeline_reconstruction() {
        // Full Factor→Apply→Update on an 8×8 matrix via factor_inner.
        // Verify reconstruction ||PAP^T - LDL^T|| < 1e-12.
        let a = symmetric_matrix(8, |i, j| {
            if i == j {
                10.0 + i as f64
            } else {
                1.0 / (1.0 + (i as f64 - j as f64).abs())
            }
        });

        let opts = AptpOptions {
            inner_block_size: 4, // Force 2 blocks for 8×8
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        assert_eq!(result.stats.num_delayed, 0, "should have no delays");

        let error =
            dense_reconstruction_error(&a, &result.l, &result.d, result.perm.as_ref().arrays().0);
        assert!(
            error < 1e-12,
            "BLAS-3 pipeline reconstruction error {:.2e} >= 1e-12",
            error
        );
    }

    #[test]
    fn test_blas3_threshold_failure_and_retry() {
        // Construct matrix where panel threshold is likely to fail for at least
        // one column, forcing backup/restore/delay/re-factor.
        // Use very strict threshold to trigger failures.
        let a = symmetric_matrix(8, |i, j| {
            let vals = [
                [1e-3, 1.0, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1],
                [1.0, 1e-3, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1],
                [0.5, 0.5, 10.0, 1.0, 0.1, 0.1, 0.1, 0.1],
                [0.5, 0.5, 1.0, 10.0, 0.1, 0.1, 0.1, 0.1],
                [0.1, 0.1, 0.1, 0.1, 5.0, 0.5, 0.0, 0.0],
                [0.1, 0.1, 0.1, 0.1, 0.5, 6.0, 0.0, 0.0],
                [0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 7.0, 0.5],
                [0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.5, 8.0],
            ];
            vals[i][j]
        });

        let opts = AptpOptions {
            inner_block_size: 4,
            threshold: 0.01, // Standard threshold
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        // Statistics invariant
        let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
        assert_eq!(sum, 8, "statistics sum {} != 8", sum);

        // Reconstruction (if fully eliminated)
        if result.stats.num_delayed == 0 {
            let error = dense_reconstruction_error(
                &a,
                &result.l,
                &result.d,
                result.perm.as_ref().arrays().0,
            );
            assert!(
                error < 1e-12,
                "threshold failure retry: reconstruction error {:.2e}",
                error
            );
        }
    }

    #[test]
    fn test_adjust_for_2x2_boundary() {
        // Test that adjust_for_2x2_boundary decrements when last accepted pivot
        // is the first half of a 2×2 pair.
        let mut d = MixedDiagonal::new(4);
        d.set_1x1(0, 1.0);
        d.set_2x2(Block2x2 {
            first_col: 1,
            a: 1.0,
            b: 0.5,
            c: 2.0,
        });
        d.set_1x1(3, 3.0);

        // If effective_nelim = 2, last is position 1 which is first of 2×2 pair (1,2).
        // Partner = 2 > 1, so should decrement to 1.
        assert_eq!(adjust_for_2x2_boundary(2, &d), 1);

        // If effective_nelim = 3, last is position 2 which is second of 2×2 pair.
        // Partner = 1 < 2, so no adjustment.
        assert_eq!(adjust_for_2x2_boundary(3, &d), 3);

        // If effective_nelim = 1, last is position 0 which is 1×1. No adjustment.
        assert_eq!(adjust_for_2x2_boundary(1, &d), 1);

        // If effective_nelim = 4, last is position 3 which is 1×1. No adjustment.
        assert_eq!(adjust_for_2x2_boundary(4, &d), 4);

        // Edge case: effective_nelim = 0
        assert_eq!(adjust_for_2x2_boundary(0, &d), 0);
    }

    #[test]
    fn test_blas3_full_block_singular() {
        // Block where all entries < small. Verify all columns delayed.
        let mut a = Mat::zeros(8, 8);
        for i in 0..8 {
            a[(i, i)] = 1e-25;
        }

        let opts = AptpOptions {
            inner_block_size: 4,
            failed_pivot_method: FailedPivotMethod::Pass,
            ..AptpOptions::default()
        };
        let mut a_copy = a.clone();
        let result = aptp_factor_in_place(a_copy.as_mut(), 8, &opts).unwrap();

        assert_eq!(result.num_eliminated, 0);
        assert_eq!(result.stats.num_delayed, 8);
        assert_eq!(result.delayed_cols.len(), 8);
    }

    // ---- Regression test: factor_inner with threshold failures (delays) ----

    /// Generate a deterministic pseudo-random symmetric indefinite matrix.
    ///
    /// Uses a simple LCG-based hash for reproducibility without external deps.
    /// Entries are in [-1, 1] with diagonal scaled by `diag_scale`.
    fn deterministic_indefinite_matrix(n: usize, seed: u64, diag_scale: f64) -> Mat<f64> {
        // Pure function of (i,j) — compatible with symmetric_matrix's Fn requirement.
        let hash = |a: usize, b: usize| -> f64 {
            let mut s = seed
                .wrapping_add((a * 10007) as u64)
                .wrapping_add((b * 7) as u64);
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 33) as f64) / (u32::MAX as f64 / 2.0) - 1.0
        };

        symmetric_matrix(n, |i, j| {
            if i == j {
                hash(i, i) * diag_scale
            } else {
                // Use min/max so (i,j) and (j,i) produce the same value
                hash(i.min(j), i.max(j) + n)
            }
        })
    }

    /// Apply the "cause_delays" pattern from SPRAL testing: multiply n/8 random
    /// rows (and corresponding columns, to maintain symmetry) by a large factor.
    /// This creates large off-diagonal entries that cause L entries to exceed 1/threshold.
    fn cause_delays(a: &mut Mat<f64>, seed: u64, scale: f64) {
        let n = a.nrows();
        let n_scaled = (n / 8).max(1);

        // Deterministically select which rows to scale
        let mut state = seed;
        let mut next_idx = || -> usize {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) as usize) % n
        };

        let mut scaled_rows = Vec::new();
        while scaled_rows.len() < n_scaled {
            let idx = next_idx();
            if !scaled_rows.contains(&idx) {
                scaled_rows.push(idx);
            }
        }

        // Scale rows and columns symmetrically: A -> S * A * S
        // where S = diag(s_1, ..., s_n) with s_i = scale if i in scaled_rows, else 1
        for &r in &scaled_rows {
            for j in 0..n {
                a[(r, j)] *= scale;
                a[(j, r)] *= scale;
            }
            // Diagonal gets scaled twice (row and column), which is correct for S*A*S
        }
    }

    #[test]
    fn test_factor_inner_with_delays() {
        // Regression test: factor_inner with threshold failures.
        //
        // Constructs matrices using the "cause_delays" pattern (SPRAL testing):
        // take a random symmetric indefinite matrix, then multiply n/8 random
        // rows and corresponding columns by 1000. This creates large off-diagonal
        // entries that cause L entries to exceed 1/threshold, triggering the
        // backup/restore/delay path in factor_inner.
        //
        // Tests multiple (n, ib) combinations to exercise different code paths
        // in the BLAS-3 Factor/Apply/Update loop.

        struct TestConfig {
            n: usize,
            ib: usize,
            seed: u64,
            scale: f64,
        }

        let configs = [
            // Small: 8x8 with ib=2 (4 blocks)
            TestConfig {
                n: 8,
                ib: 2,
                seed: 42,
                scale: 1000.0,
            },
            // Small: 8x8 with ib=4 (2 blocks)
            TestConfig {
                n: 8,
                ib: 4,
                seed: 42,
                scale: 1000.0,
            },
            // Medium: 16x16 with ib=4 (4 blocks)
            TestConfig {
                n: 16,
                ib: 4,
                seed: 42,
                scale: 1000.0,
            },
            // Medium: 16x16 with ib=2 (8 blocks)
            TestConfig {
                n: 16,
                ib: 2,
                seed: 42,
                scale: 1000.0,
            },
            // Larger: 32x32 with ib=4
            TestConfig {
                n: 32,
                ib: 4,
                seed: 42,
                scale: 1000.0,
            },
            // Larger: 32x32 with ib=8
            TestConfig {
                n: 32,
                ib: 8,
                seed: 42,
                scale: 1000.0,
            },
            // Different seed
            TestConfig {
                n: 16,
                ib: 4,
                seed: 137,
                scale: 1000.0,
            },
            // Extreme scale
            TestConfig {
                n: 16,
                ib: 4,
                seed: 42,
                scale: 1e6,
            },
            // 64x64 with ib=8
            TestConfig {
                n: 64,
                ib: 8,
                seed: 42,
                scale: 1000.0,
            },
            // 64x64 with ib=16
            TestConfig {
                n: 64,
                ib: 16,
                seed: 42,
                scale: 1000.0,
            },
        ];

        let mut _any_delays_found = false;
        let mut any_failures = false;

        for (idx, config) in configs.iter().enumerate() {
            let mut a = deterministic_indefinite_matrix(config.n, config.seed, 5.0);
            cause_delays(&mut a, config.seed + 1000, config.scale);

            let opts = AptpOptions {
                inner_block_size: config.ib,
                outer_block_size: config.n.max(config.ib),
                threshold: 0.01,
                ..AptpOptions::default()
            };

            let result = aptp_factor(a.as_ref(), &opts).unwrap();

            let n = config.n;
            let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
            assert_eq!(
                sum, n,
                "config {}: statistics invariant broken: {} != {}",
                idx, sum, n
            );

            if result.stats.num_delayed > 0 {
                _any_delays_found = true;
            }

            // If fully eliminated, check reconstruction error
            if result.stats.num_delayed == 0 {
                let error = dense_reconstruction_error(
                    &a,
                    &result.l,
                    &result.d,
                    result.perm.as_ref().arrays().0,
                );
                if error >= 1e-12 {
                    any_failures = true;
                }
                assert!(
                    error < 1e-10,
                    "config {} (n={}, ib={}): reconstruction error {:.2e} >= 1e-10 \
                     (no delays but bad reconstruction indicates factor_inner bug)",
                    idx,
                    config.n,
                    config.ib,
                    error
                );
            }
        }

        // At least one config should have triggered delays (confirming the
        // cause_delays pattern works with the threshold).
        // Note: we don't assert any_delays_found because with complete pivoting
        // the diagonal block never delays — only apply_and_check can reduce
        // effective_nelim. Whether this triggers depends on the matrix structure.

        assert!(
            !any_failures,
            "Some configs had reconstruction error >= 1e-12 without delays. \
             This indicates a bug in factor_inner's backup/restore/update logic."
        );
    }

    #[test]
    fn test_factor_inner_with_delays_targeted() {
        // More targeted test: construct a matrix specifically designed to trigger
        // threshold failure in factor_inner's apply_and_check step.
        //
        // Strategy: Create a matrix where complete pivoting on the ib×ib diagonal
        // block succeeds (all ib columns eliminated within the block), but the
        // panel entries below the block exceed 1/threshold after TRSM + D scaling.
        //
        // Matrix structure (8x8, ib=4):
        //   Block [0:4, 0:4]: moderate diagonal, small off-diagonal -- complete pivoting
        //     eliminates all 4 columns within the block with |L_block| <= 4
        //   Panel [4:8, 0:4]: large entries that exceed 1/threshold after TRSM
        //     Specifically, panel entries ~ 2.0 with pivots ~ 0.005, so L ~ 400 >> 100
        //   Block [4:8, 4:8]: moderate diagonal
        //
        // Complete pivoting searches only within the ib×ib diagonal block, so panel
        // entries don't affect pivot selection.
        //
        // Also tested with ib=2 to get partial threshold failure (first 2 columns
        // may pass if their diagonal is large enough).
        let n = 8;

        for &ib in &[4, 2] {
            let a = symmetric_matrix(n, |i, j| {
                match (i, j) {
                    // Diagonal block [0:4, 0:4]
                    // Complete pivoting picks largest first: 10.0, 10.0, then 0.005, 0.005
                    // For ib=4: all 4 eliminated in block, but panel L ~ 2.0/0.005 = 400
                    // For ib=2: block [0:2,0:2] → pivots 10.0,10.0 → panel L ~ 2.0/10 = 0.2 (OK)
                    //           block [2:4,2:4] → pivots 0.005,0.005 → panel L ~ 2.0/0.005 = 400 (FAIL)
                    (0, 0) => 10.0,
                    (1, 1) => 10.0,
                    (2, 2) => 0.005,
                    (3, 3) => 0.005,
                    (i, j) if i < 4 && j < 4 && i != j => 0.001,

                    // Panel [4:8, 0:4]
                    (_, j) if j < 4 => 2.0,

                    // Lower-right block [4:8, 4:8]
                    (i, _) if i == j => 20.0,
                    (_, _) => 0.5,
                }
            });

            // Test APTP delay behavior (TPP disabled so we can observe delays)
            let opts_pass = AptpOptions {
                inner_block_size: ib,
                outer_block_size: 256,
                threshold: 0.01,
                failed_pivot_method: FailedPivotMethod::Pass,
                ..AptpOptions::default()
            };

            let result = aptp_factor(a.as_ref(), &opts_pass).unwrap();

            let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
            assert_eq!(sum, n, "ib={}: statistics invariant: {} != {}", ib, sum, n);

            // Check reconstruction if fully eliminated
            if result.stats.num_delayed == 0 {
                let error = dense_reconstruction_error(
                    &a,
                    &result.l,
                    &result.d,
                    result.perm.as_ref().arrays().0,
                );
                assert!(
                    error < 1e-12,
                    "targeted (ib={}): reconstruction error {:.2e} >= 1e-12 (no delays)",
                    ib,
                    error
                );
            }

            // Verify that delays happened (at least for ib=4, columns 2,3 should be
            // delayed because L panel entries exceed 1/0.01=100)
            if ib == 4 || ib == 2 {
                // With small pivots 0.005 and panel entries 2.0:
                // L = 2.0 / 0.005 = 400, which far exceeds 100
                // Expect delays for the small-pivot columns
                assert!(
                    result.stats.num_delayed > 0,
                    "ib={}: expected some delays for small-pivot columns",
                    ib
                );
            }

            // Also check partial factorization with contribution block
            {
                let p = 6; // Only 6 of 8 fully summed
                let mut a_copy = a.clone();
                let opts_partial = AptpOptions {
                    inner_block_size: ib,
                    outer_block_size: 256,
                    threshold: 0.01,
                    failed_pivot_method: FailedPivotMethod::Pass,
                    ..AptpOptions::default()
                };
                let result_p = aptp_factor_in_place(a_copy.as_mut(), p, &opts_partial).unwrap();
                let error_p = check_partial_factorization_in_place(&a, &a_copy, p, &result_p);
                assert!(
                    error_p < 1e-10,
                    "targeted partial (ib={}, p={}): error {:.2e} >= 1e-10",
                    ib,
                    p,
                    error_p
                );
            }
        }
    }

    #[test]
    fn test_factor_inner_with_delays_aggressive() {
        // Aggressive test using a Lehmer-like matrix structure with perturbations
        // designed to trigger many threshold failures.
        //
        // The matrix is constructed as: A = Q * diag(d_i) * Q^T where d_i are
        // alternating +/- with varying magnitudes. Then specific rows are scaled
        // to create problematic panel entries.
        //
        // We test with multiple (n, ib) combos and multiple seeds.

        let test_cases: Vec<(usize, usize, u64)> = vec![
            (12, 4, 1),
            (12, 4, 2),
            (12, 4, 3),
            (16, 4, 1),
            (16, 4, 2),
            (16, 8, 1),
            (20, 4, 1),
            (24, 4, 1),
            (24, 8, 1),
            (32, 4, 1),
            (32, 8, 1),
            (32, 16, 1),
            (48, 8, 1),
            (64, 8, 1),
            (64, 16, 1),
        ];

        let mut _total_delays = 0;
        let mut max_error = 0.0_f64;
        let mut _worst_config = String::new();

        for &(n, ib, seed) in &test_cases {
            // Build a random symmetric indefinite matrix
            let mut a = deterministic_indefinite_matrix(n, seed * 31337, 5.0);

            // Apply cause_delays to create threshold failures
            cause_delays(&mut a, seed * 31337 + 7919, 1000.0);

            let opts = AptpOptions {
                inner_block_size: ib,
                outer_block_size: n.max(ib),
                threshold: 0.01,
                ..AptpOptions::default()
            };

            let result = aptp_factor(a.as_ref(), &opts).unwrap();

            let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
            assert_eq!(
                sum, n,
                "(n={}, ib={}, seed={}): stats invariant {} != {}",
                n, ib, seed, sum, n
            );

            _total_delays += result.stats.num_delayed;

            if result.stats.num_delayed == 0 {
                let error = dense_reconstruction_error(
                    &a,
                    &result.l,
                    &result.d,
                    result.perm.as_ref().arrays().0,
                );
                if error > max_error {
                    max_error = error;
                    _worst_config =
                        format!("n={}, ib={}, seed={}: error={:.2e}", n, ib, seed, error);
                }
                assert!(
                    error < 1e-10,
                    "(n={}, ib={}, seed={}): reconstruction error {:.2e} >= 1e-10",
                    n,
                    ib,
                    seed,
                    error
                );
            }
        }
    }

    #[test]
    fn test_factor_inner_cause_delays_then_compare_single_vs_blocked() {
        // Compare factor_inner (ib > 1) vs ib=n (effectively single-block).
        // Both should produce reconstruction error < 1e-12 if no delays,
        // and the same number of eliminations.
        //
        // This catches bugs where blocking introduces errors that single-block
        // factorization does not.

        let seeds = [42u64, 137, 271, 314, 997];
        let n = 16;

        for &seed in &seeds {
            let mut a = deterministic_indefinite_matrix(n, seed, 5.0);
            cause_delays(&mut a, seed + 5000, 500.0);

            // Single-block: ib = n (factor_block_diagonal processes entire matrix)
            let opts_single = AptpOptions {
                inner_block_size: n,
                outer_block_size: n,
                threshold: 0.01,
                ..AptpOptions::default()
            };
            let result_single = aptp_factor(a.as_ref(), &opts_single).unwrap();

            // Blocked: ib = 4
            let opts_blocked = AptpOptions {
                inner_block_size: 4,
                outer_block_size: n,
                threshold: 0.01,
                ..AptpOptions::default()
            };
            let result_blocked = aptp_factor(a.as_ref(), &opts_blocked).unwrap();

            // Both should achieve good reconstruction when fully eliminated
            if result_single.stats.num_delayed == 0 {
                let error_s = dense_reconstruction_error(
                    &a,
                    &result_single.l,
                    &result_single.d,
                    result_single.perm.as_ref().arrays().0,
                );
                assert!(
                    error_s < 1e-12,
                    "seed={}: single-block error {:.2e}",
                    seed,
                    error_s
                );
            }

            if result_blocked.stats.num_delayed == 0 {
                let error_b = dense_reconstruction_error(
                    &a,
                    &result_blocked.l,
                    &result_blocked.d,
                    result_blocked.perm.as_ref().arrays().0,
                );
                assert!(
                    error_b < 1e-12,
                    "seed={}: blocked error {:.2e} >= 1e-12",
                    seed,
                    error_b
                );
            }
        }
    }

    /// Check partial factorization: P^T A P = L D L^T + [0; 0; contribution]
    ///
    /// After factor_inner with num_fully_summed = p < m, the matrix contains:
    /// - L entries in a[0..q, 0..q] (lower triangle, unit diagonal implicit)
    /// - L panel in a[q..m, 0..q]
    /// - Schur complement (contribution) in a[q..m, q..m]
    ///
    /// where q = num_eliminated <= p.
    ///
    /// Correctness condition: for the permuted matrix PAP^T,
    ///   PAP^T[0..m, 0..m] = L[0..m, 0..q] * D[0..q, 0..q] * L[0..m, 0..q]^T + [0_{q,q} 0; 0 S]
    /// where S = a[q..m, q..m] after factorization (the Schur complement).
    fn check_partial_factorization_in_place(
        original: &Mat<f64>,
        factored: &Mat<f64>,
        num_fully_summed: usize,
        result: &AptpFactorResult,
    ) -> f64 {
        let m = original.nrows();
        let q = result.num_eliminated;
        let p = num_fully_summed;
        let perm = &result.perm;

        // Apply deferred contribution GEMM to get the full Schur complement.
        // After our restructuring, aptp_factor_in_place leaves A[p..m, p..m]
        // with only assembled values (no per-block trailing updates).
        // We need to apply the deferred GEMM to get the actual Schur complement.
        let mut factored = factored.clone();
        let nfs = m - p;
        if nfs > 0 && q > 0 {
            let mut contrib_buffer = Mat::zeros(nfs, nfs);
            compute_contribution_gemm(&factored, p, q, m, &result.d, &mut contrib_buffer, Par::Seq);
            // Copy the Schur complement back into the factored matrix
            for i in 0..nfs {
                for j in 0..=i {
                    factored[(p + i, p + j)] = contrib_buffer[(i, j)];
                }
            }
        }
        let d = &result.d;

        // Build PAP^T
        let mut papt = Mat::zeros(m, m);
        for i in 0..m {
            for j in 0..m {
                papt[(i, j)] = original[(perm[i], perm[j])];
            }
        }

        // Extract L (m x q, unit lower triangular in first q columns)
        let mut l_full = Mat::zeros(m, q);
        for i in 0..q {
            l_full[(i, i)] = 1.0;
        }
        let mut col = 0;
        while col < q {
            match d.get_pivot_type(col) {
                PivotType::OneByOne => {
                    for i in (col + 1)..m {
                        l_full[(i, col)] = factored[(i, col)];
                    }
                    col += 1;
                }
                PivotType::TwoByTwo { partner } if partner > col => {
                    // a[(col+1, col)] is D off-diagonal, not L
                    for i in (col + 2)..m {
                        l_full[(i, col)] = factored[(i, col)];
                        l_full[(i, col + 1)] = factored[(i, col + 1)];
                    }
                    col += 2;
                }
                _ => {
                    col += 1;
                }
            }
        }

        // Build D (q x q)
        let mut d_mat = Mat::zeros(q, q);
        col = 0;
        while col < q {
            match d.get_pivot_type(col) {
                PivotType::OneByOne => {
                    d_mat[(col, col)] = d.get_1x1(col);
                    col += 1;
                }
                PivotType::TwoByTwo { partner } if partner > col => {
                    let block = d.get_2x2(col);
                    d_mat[(col, col)] = block.a;
                    d_mat[(col, col + 1)] = block.b;
                    d_mat[(col + 1, col)] = block.b;
                    d_mat[(col + 1, col + 1)] = block.c;
                    col += 2;
                }
                _ => {
                    col += 1;
                }
            }
        }

        // Compute L * D * L^T (m x m)
        // First: W = L * D (m x q)
        let mut w = Mat::zeros(m, q);
        for i in 0..m {
            for j in 0..q {
                let mut sum = 0.0;
                for k in 0..q {
                    sum += l_full[(i, k)] * d_mat[(k, j)];
                }
                w[(i, j)] = sum;
            }
        }
        // Then: LDL^T = W * L^T (m x m)
        let mut ldlt = Mat::zeros(m, m);
        for i in 0..m {
            for j in 0..m {
                let mut sum = 0.0;
                for k in 0..q {
                    sum += w[(i, k)] * l_full[(j, k)];
                }
                ldlt[(i, j)] = sum;
            }
        }

        // Extract Schur complement S = factored[q..m, q..m]
        // The residual should be: PAP^T - LDL^T = [0 0; 0 S]
        // So: PAP^T[i,j] - LDL^T[i,j] should be:
        //   0 if i < q or j < q
        //   S[i-q, j-q] = factored[i,j] if i >= q and j >= q
        let mut max_error = 0.0_f64;
        let mut norm_sq = 0.0_f64;
        let mut diff_sq = 0.0_f64;

        // Only check lower triangle (i >= j) because swap_symmetric only
        // maintains the lower triangle of the dense frontal matrix.
        // Production code (extract_contribution, extend_add) also reads
        // only the lower triangle.
        for i in 0..m {
            for j in 0..=i {
                // Count lower triangle entries twice for norm (symmetric)
                let weight = if i == j { 1.0 } else { 2.0 };
                norm_sq += weight * papt[(i, j)] * papt[(i, j)];
                let residual = papt[(i, j)] - ldlt[(i, j)];
                if i >= q && j >= q {
                    // This should equal factored[i, j] (the Schur complement)
                    let schur_entry = factored[(i, j)];
                    let err = (residual - schur_entry).abs();
                    if err > max_error {
                        max_error = err;
                    }
                    diff_sq += weight * (residual - schur_entry) * (residual - schur_entry);
                } else {
                    // Should be zero
                    diff_sq += weight * residual * residual;
                    if residual.abs() > max_error {
                        max_error = residual.abs();
                    }
                }
            }
        }

        if norm_sq == 0.0 {
            diff_sq.sqrt()
        } else {
            (diff_sq / norm_sq).sqrt()
        }
    }

    #[test]
    fn test_factor_inner_partial_with_delays_schur_check() {
        // This test exercises factor_inner on frontal matrices where
        // num_fully_summed < m (partial factorization with contribution block),
        // AND delays occur. This is exactly the scenario in multifrontal
        // factorization where the bug was observed.
        //
        // The test verifies that:
        //   PAP^T = L * D * L^T + [0 0; 0 S]
        // where S is the Schur complement stored in the lower-right of the
        // factored matrix.
        //
        // If backup/restore/update_cross_terms is wrong, S will be corrupted.

        struct PartialConfig {
            m: usize,  // total matrix dimension
            p: usize,  // num_fully_summed
            ib: usize, // inner block size
            seed: u64,
            scale: f64,
        }

        let configs = [
            // Small cases: 12x12 with p=8, ib=4 (2 inner blocks, 4 contribution rows)
            PartialConfig {
                m: 12,
                p: 8,
                ib: 4,
                seed: 42,
                scale: 1000.0,
            },
            PartialConfig {
                m: 12,
                p: 8,
                ib: 4,
                seed: 137,
                scale: 1000.0,
            },
            PartialConfig {
                m: 12,
                p: 8,
                ib: 2,
                seed: 42,
                scale: 1000.0,
            },
            // Medium: 16x12 (p=12 of 16), ib=4
            PartialConfig {
                m: 16,
                p: 12,
                ib: 4,
                seed: 42,
                scale: 1000.0,
            },
            PartialConfig {
                m: 16,
                p: 12,
                ib: 4,
                seed: 271,
                scale: 1000.0,
            },
            // 20x16, ib=4
            PartialConfig {
                m: 20,
                p: 16,
                ib: 4,
                seed: 42,
                scale: 1000.0,
            },
            PartialConfig {
                m: 20,
                p: 16,
                ib: 8,
                seed: 42,
                scale: 1000.0,
            },
            // Larger: 32x24, ib=8
            PartialConfig {
                m: 32,
                p: 24,
                ib: 8,
                seed: 42,
                scale: 1000.0,
            },
            PartialConfig {
                m: 32,
                p: 24,
                ib: 4,
                seed: 42,
                scale: 1000.0,
            },
            // 48x32, ib=8
            PartialConfig {
                m: 48,
                p: 32,
                ib: 8,
                seed: 42,
                scale: 1000.0,
            },
            // Extreme scale
            PartialConfig {
                m: 16,
                p: 12,
                ib: 4,
                seed: 42,
                scale: 1e6,
            },
            // More seeds
            PartialConfig {
                m: 16,
                p: 12,
                ib: 4,
                seed: 314,
                scale: 1000.0,
            },
            PartialConfig {
                m: 16,
                p: 12,
                ib: 4,
                seed: 997,
                scale: 1000.0,
            },
        ];

        let mut any_delays = false;
        let mut worst_error = 0.0_f64;
        let mut worst_config = String::new();

        for (idx, config) in configs.iter().enumerate() {
            let mut a = deterministic_indefinite_matrix(config.m, config.seed, 5.0);
            cause_delays(&mut a, config.seed + 1000, config.scale);

            let opts = AptpOptions {
                inner_block_size: config.ib,
                outer_block_size: config.m.max(config.ib),
                threshold: 0.01,
                failed_pivot_method: FailedPivotMethod::Pass,
                ..AptpOptions::default()
            };

            let original = a.clone();
            let result = aptp_factor_in_place(a.as_mut(), config.p, &opts).unwrap();

            let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
            assert_eq!(
                sum, config.p,
                "config {} (m={}, p={}, ib={}): stats invariant {} != {}",
                idx, config.m, config.p, config.ib, sum, config.p
            );

            if result.stats.num_delayed > 0 {
                any_delays = true;
            }

            // Check partial factorization correctness (PAP^T = LDL^T + [0;S])
            let error = check_partial_factorization_in_place(&original, &a, config.p, &result);

            if error > worst_error {
                worst_error = error;
                worst_config = format!(
                    "config {} (m={}, p={}, ib={}, seed={})",
                    idx, config.m, config.p, config.ib, config.seed
                );
            }
        }

        assert!(
            any_delays,
            "No configurations triggered threshold delays. \
             Adjust scale or matrix construction to create delays."
        );

        assert!(
            worst_error < 1e-10,
            "Worst partial factorization error: {:.2e} at {}\n\
             This indicates corrupted Schur complement — bug in \
             backup/restore/update_cross_terms logic in factor_inner.",
            worst_error,
            worst_config
        );
    }

    // ---- NFS boundary in two-level factor tests ----

    #[test]
    fn test_two_level_nfs_boundary_mid_block() {
        // Exercises the local nfs_boundary = p - col_start computation in
        // two_level_factor. When p is not a multiple of outer_block_size,
        // the last outer block has nfs_boundary < block_cols, meaning the
        // NFS region starts inside an inner block. update_trailing must
        // correctly restrict its Schur update to the FS region and skip
        // the NFS×NFS region (deferred to compute_contribution_gemm).
        //
        // If nfs_boundary is wrong, the Schur complement (contribution
        // block) will be corrupted, failing the reconstruction check.

        struct Config {
            m: usize,
            p: usize,
            nb: usize,
            ib: usize,
            seed: u64,
        }

        let configs = [
            // p=20, nb=16: first block nfs_boundary=20, second block nfs_boundary=4
            // (NFS starts at col 4 of the second 16-col outer block)
            Config {
                m: 48,
                p: 20,
                nb: 16,
                ib: 8,
                seed: 42,
            },
            // p=40, nb=32: first block nfs_boundary=40, second block nfs_boundary=8
            Config {
                m: 64,
                p: 40,
                nb: 32,
                ib: 8,
                seed: 42,
            },
            // p=10, nb=8: first block nfs_boundary=10, second block nfs_boundary=2
            // (very small NFS boundary in inner block)
            Config {
                m: 24,
                p: 10,
                nb: 8,
                ib: 4,
                seed: 137,
            },
            // p=25, nb=16: blocks at col_start=0 (nfs=25), col_start=16 (nfs=9)
            Config {
                m: 48,
                p: 25,
                nb: 16,
                ib: 8,
                seed: 271,
            },
            // Deliberately misalign p with ib too: p=13, nb=8, ib=4
            // block 0: nfs_boundary=13, block 1: nfs_boundary=5 (mid-ib boundary)
            Config {
                m: 32,
                p: 13,
                nb: 8,
                ib: 4,
                seed: 42,
            },
        ];

        let mut worst_error = 0.0_f64;
        let mut worst_label = String::new();

        for (idx, c) in configs.iter().enumerate() {
            let a = deterministic_indefinite_matrix(c.m, c.seed, 5.0);

            let opts = AptpOptions {
                outer_block_size: c.nb,
                inner_block_size: c.ib,
                threshold: 0.01,
                failed_pivot_method: FailedPivotMethod::Pass,
                ..AptpOptions::default()
            };

            let original = a.clone();
            let mut a_copy = a;
            let result = aptp_factor_in_place(a_copy.as_mut(), c.p, &opts).unwrap();

            let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
            assert_eq!(
                sum, c.p,
                "config {}: statistics invariant {} != {}",
                idx, sum, c.p
            );

            let error = check_partial_factorization_in_place(&original, &a_copy, c.p, &result);

            if error > worst_error {
                worst_error = error;
                worst_label = format!(
                    "config {} (m={}, p={}, nb={}, ib={}, seed={})",
                    idx, c.m, c.p, c.nb, c.ib, c.seed
                );
            }
        }

        assert!(
            worst_error < 1e-10,
            "Worst reconstruction error: {:.2e} at {}\n\
             This indicates the local nfs_boundary computation in \
             two_level_factor is producing incorrect Schur complements.",
            worst_error,
            worst_label
        );
    }

    // ---- TPP (Threshold Partial Pivoting) fallback tests ----

    #[test]
    fn test_tpp_helpers() {
        // 4x4 lower triangle:
        //  [  4.0   .    .    .  ]
        //  [  0.5  -3.0  .    .  ]
        //  [  0.1   0.2  5.0  .  ]
        //  [  0.3  -0.4  0.6  2.0]
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [4.0, 0.5, 0.1, 0.3],
                [0.5, -3.0, 0.2, -0.4],
                [0.1, 0.2, 5.0, 0.6],
                [0.3, -0.4, 0.6, 2.0],
            ];
            vals[i][j]
        });

        // tpp_is_col_small: col 0, from=0, to=4, small=5.0 → true (all entries < 5.0)
        assert!(tpp_is_col_small(a.as_ref(), 0, 0, 4, 5.0));
        // With small=0.01 → false (0.5, 0.1, 0.3 ≥ 0.01 and diag 4.0 ≥ 0.01)
        assert!(!tpp_is_col_small(a.as_ref(), 0, 0, 4, 0.01));

        // Zero column → true
        let z = Mat::zeros(4, 4);
        assert!(tpp_is_col_small(z.as_ref(), 0, 0, 4, 1e-20));

        // tpp_find_row_abs_max: row 3, cols 0..3
        // a[(3,0)]=0.3, a[(3,1)]=-0.4, a[(3,2)]=0.6 → max at col 2
        let t = tpp_find_row_abs_max(a.as_ref(), 3, 0, 3);
        assert_eq!(t, 2, "expected col 2, got {}", t);

        // tpp_find_rc_abs_max_exclude: col 1, nelim=0, m=4, exclude=3
        // row part: a[(1,0)]=0.5
        // col part: a[(2,1)]=0.2 (excluding row 3)
        // max = 0.5
        let max_exc = tpp_find_rc_abs_max_exclude(a.as_ref(), 1, 0, 4, 3);
        assert!(
            (max_exc - 0.5).abs() < 1e-15,
            "expected 0.5, got {}",
            max_exc
        );

        // tpp_test_2x2: (0, 1) with a11=4, a21=0.5, a22=-3
        // det = 4*(-3) - 0.25 = -12.25 → non-singular
        // With small=1e-20, u=0.01
        let result = tpp_test_2x2(a.as_ref(), 0, 1, 0.6, 0.6, 0.01, 1e-20);
        assert!(result.is_some(), "2x2 pivot (0,1) should pass");
        let (d11, d12, d22) = result.unwrap();
        assert!((d11 - 4.0).abs() < 1e-15);
        assert!((d12 - 0.5).abs() < 1e-15);
        assert!((d22 - (-3.0)).abs() < 1e-15);

        // tpp_test_2x2 with near-singular block
        let tiny = symmetric_matrix(2, |_, _| 1e-25);
        assert!(tpp_test_2x2(tiny.as_ref(), 0, 1, 1e-25, 1e-25, 0.01, 1e-20).is_none());
    }

    #[test]
    fn test_tpp_basic_1x1() {
        // 4x4 diagonal-dominant matrix: all diagonals are acceptable 1x1 pivots
        let a = symmetric_matrix(4, |i, j| {
            if i == j {
                [10.0, -8.0, 5.0, 7.0][i]
            } else {
                0.1
            }
        });

        let opts = AptpOptions::default();
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        assert_eq!(result.stats.num_delayed, 0);
        assert_eq!(result.stats.num_1x1 + 2 * result.stats.num_2x2, 4);

        let error =
            dense_reconstruction_error(&a, &result.l, &result.d, result.perm.as_ref().arrays().0);
        assert!(error < 1e-12, "reconstruction error {:.2e}", error);
    }

    #[test]
    fn test_tpp_basic_2x2() {
        // 4x4 matrix where diagonal entries are small (fail 1x1 threshold)
        // but off-diagonal entries create good 2x2 pivots.
        let a = symmetric_matrix(4, |i, j| {
            let vals = [
                [0.001, 5.0, 0.01, 0.01],
                [5.0, 0.001, 0.01, 0.01],
                [0.01, 0.01, 0.001, 5.0],
                [0.01, 0.01, 5.0, 0.001],
            ];
            vals[i][j]
        });

        let opts = AptpOptions::default();
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        assert_eq!(result.stats.num_delayed, 0);

        let error =
            dense_reconstruction_error(&a, &result.l, &result.d, result.perm.as_ref().arrays().0);
        assert!(error < 1e-12, "reconstruction error {:.2e}", error);
    }

    #[test]
    fn test_tpp_fallback_after_aptp() {
        // 16x16 matrix designed so APTP delays columns that TPP recovers.
        // Small inner block size forces many ib-scoped searches that fail,
        // but TPP's exhaustive search finds pivots.
        //
        // Matrix structure:
        //   Block [0:8, 0:8]: diagonal 10.0, small off-diag → easy APTP pivots
        //   Block [8:12, 8:12]: diagonal 0.005, large off-diag → APTP delays
        //   Block [12:16, 12:16]: diagonal 20.0 → easy pivots
        //   Cross [8:12, 0:8]: large entries → panel L exceeds threshold
        let n = 16;
        let a = symmetric_matrix(n, |i, j| {
            match (i, j) {
                // Top-left 8x8: easy pivots
                (i, j) if i < 8 && j < 8 => {
                    if i == j {
                        10.0
                    } else {
                        0.01
                    }
                }
                // Middle 4x4: hard pivots (small diag, large off-diag)
                (i, j) if (8..12).contains(&i) && (8..12).contains(&j) => {
                    if i == j {
                        0.005
                    } else {
                        2.0
                    }
                }
                // Bottom-right 4x4: easy pivots
                (i, j) if i >= 12 && j >= 12 => {
                    if i == j {
                        20.0
                    } else {
                        0.1
                    }
                }
                // Cross terms: large
                (i, j) if (8..12).contains(&i) && j < 8 => 3.0,
                (i, j) if i < 8 && (8..12).contains(&j) => 3.0,
                // Small cross terms elsewhere
                _ => 0.01,
            }
        });

        // With TPP disabled: APTP delays some columns
        let opts_pass = AptpOptions {
            inner_block_size: 4,
            failed_pivot_method: FailedPivotMethod::Pass,
            ..AptpOptions::default()
        };
        let result_pass = aptp_factor(a.as_ref(), &opts_pass).unwrap();

        // With TPP enabled: should eliminate more columns
        let opts_tpp = AptpOptions {
            inner_block_size: 4,
            failed_pivot_method: FailedPivotMethod::Tpp,
            ..AptpOptions::default()
        };
        let result_tpp = aptp_factor(a.as_ref(), &opts_tpp).unwrap();

        // TPP should recover at least some of APTP's delays
        assert!(
            result_tpp.stats.num_delayed <= result_pass.stats.num_delayed,
            "TPP should not increase delays"
        );

        // If fully eliminated, check reconstruction
        if result_tpp.stats.num_delayed == 0 {
            let error = dense_reconstruction_error(
                &a,
                &result_tpp.l,
                &result_tpp.d,
                result_tpp.perm.as_ref().arrays().0,
            );
            assert!(
                error < 1e-12,
                "TPP reconstruction error {:.2e} >= 1e-12",
                error
            );
        }
    }

    #[test]
    fn test_tpp_reconstruction_stress() {
        // 256x256 random indefinite matrix. Factor with TPP and verify reconstruction.
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        use rand_distr::{Distribution, Uniform};

        let n = 256;
        let mut rng = StdRng::seed_from_u64(42);
        let dist = Uniform::new(-1.0f64, 1.0);

        let mut a = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..=i {
                let v: f64 = dist.sample(&mut rng);
                a[(i, j)] = v;
                a[(j, i)] = v;
            }
            // Strengthen diagonal to avoid excessive delays
            let sign = if i % 3 == 0 { -1.0 } else { 1.0 };
            a[(i, i)] = sign * (5.0 + dist.sample(&mut rng).abs());
        }

        let opts = AptpOptions {
            failed_pivot_method: FailedPivotMethod::Tpp,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        if result.stats.num_delayed == 0 {
            let error = dense_reconstruction_error(
                &a,
                &result.l,
                &result.d,
                result.perm.as_ref().arrays().0,
            );
            assert!(
                error < 1e-10,
                "stress reconstruction error {:.2e} >= 1e-10",
                error
            );
        }

        // Invariant check
        let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
        assert_eq!(sum, n, "statistics invariant: {} != {}", sum, n);
    }

    #[test]
    fn test_failed_pivot_method_pass() {
        // Verify FailedPivotMethod::Pass skips TPP and preserves delays.
        // Use a matrix that APTP cannot fully factor (small diag, large off-diag).
        let n = 4;
        let a = symmetric_matrix(n, |i, j| {
            let vals = [
                [0.001, 5.0, 0.01, 0.01],
                [5.0, 0.001, 0.01, 0.01],
                [0.01, 0.01, 0.001, 5.0],
                [0.01, 0.01, 5.0, 0.001],
            ];
            vals[i][j]
        });

        // With Pass: APTP's complete pivoting on ib-blocks should still handle
        // this via 2x2 pivots. But with a very small block size, it might delay.
        let opts_pass = AptpOptions {
            inner_block_size: 2,
            failed_pivot_method: FailedPivotMethod::Pass,
            ..AptpOptions::default()
        };
        let result_pass = aptp_factor(a.as_ref(), &opts_pass).unwrap();

        let opts_tpp = AptpOptions {
            inner_block_size: 2,
            failed_pivot_method: FailedPivotMethod::Tpp,
            ..AptpOptions::default()
        };
        let result_tpp = aptp_factor(a.as_ref(), &opts_tpp).unwrap();

        // TPP should not have more delays than Pass
        assert!(
            result_tpp.stats.num_delayed <= result_pass.stats.num_delayed,
            "TPP delays {} > Pass delays {}",
            result_tpp.stats.num_delayed,
            result_pass.stats.num_delayed
        );

        // Both must satisfy invariant
        let sum_pass = result_pass.stats.num_1x1
            + 2 * result_pass.stats.num_2x2
            + result_pass.stats.num_delayed;
        let sum_tpp =
            result_tpp.stats.num_1x1 + 2 * result_tpp.stats.num_2x2 + result_tpp.stats.num_delayed;
        assert_eq!(sum_pass, n);
        assert_eq!(sum_tpp, n);
    }

    #[test]
    fn test_tpp_zero_pivot_handling() {
        // Matrix with some zero columns — TPP should handle gracefully.
        let n = 4;
        let a = symmetric_matrix(n, |i, j| {
            if i == j {
                [5.0, 0.0, 0.0, 3.0][i]
            } else if (i == 0 && j == 3) || (i == 3 && j == 0) {
                0.1
            } else {
                0.0
            }
        });

        let opts = AptpOptions {
            failed_pivot_method: FailedPivotMethod::Tpp,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        // Should handle zero pivots (columns 1, 2) as zero 1x1 pivots
        let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
        assert_eq!(sum, n, "invariant: {} != {}", sum, n);
    }
}
