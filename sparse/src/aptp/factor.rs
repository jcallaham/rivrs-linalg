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

use faer::linalg::matmul::triangular::{self as tri_matmul, BlockStructure};
use faer::linalg::triangular_solve;
use faer::perm::Perm;
use faer::prelude::*;
use faer::{Accum, Conj, Mat, MatMut, MatRef};

use super::diagonal::MixedDiagonal;
use super::perm::perm_from_forward;
use super::pivot::{Block2x2, PivotType};
use crate::error::SparseError;

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
}

impl Default for AptpOptions {
    fn default() -> Self {
        Self {
            threshold: 0.01,
            small: 1e-20,
            fallback: AptpFallback::BunchKaufman,
            outer_block_size: 256,
            inner_block_size: 32,
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
#[derive(Debug, Clone)]
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
    a: MatMut<'_, f64>,
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

    // Dispatch: two-level for large fronts, factor_inner for small ones
    if num_fully_summed > options.outer_block_size {
        two_level_factor(a, num_fully_summed, options)
    } else {
        factor_inner(a, num_fully_summed, options)
    }
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

/// Attempt a 1x1 pivot at column `col`.
///
/// Divides column k by diagonal d_kk, checks singularity and stability.
/// On success, L column entries are written in place (a[i, col] = l_ik).
/// On failure, the column is restored to its original state.
///
/// Returns Ok((d_value, max_l_entry)) or Err(max_l_entry).
#[cfg(test)]
#[allow(dead_code)] // Used by complete_pivoting_factor (also #[cfg(test)])
fn try_1x1_pivot(
    mut a: MatMut<'_, f64>,
    col: usize,
    threshold: f64,
    small: f64,
) -> Result<(f64, f64), f64> {
    let m = a.nrows();
    let d_kk = a[(col, col)];

    if d_kk.abs() < small {
        return Err(0.0);
    }

    let inv_d = 1.0 / d_kk;
    let stability_bound = 1.0 / threshold;

    // TODO(Phase 9.1): Replace per-column Vec backup with per-block backup
    // when two-level blocking is implemented.
    let backup: Vec<f64> = ((col + 1)..m).map(|i| a[(i, col)]).collect();

    // Compute L column: l_ik = a_ik / d_kk
    let mut max_l: f64 = 0.0;
    for i in (col + 1)..m {
        let l_ik = a[(i, col)] * inv_d;
        a[(i, col)] = l_ik;
        let abs_l = l_ik.abs();
        if abs_l > max_l {
            max_l = abs_l;
        }
    }

    // Stability check
    if max_l >= stability_bound {
        // Restore column
        for (idx, i) in ((col + 1)..m).enumerate() {
            a[(i, col)] = backup[idx];
        }
        return Err(max_l);
    }

    Ok((d_kk, max_l))
}

/// Attempt a 2x2 Bunch-Kaufman pivot using columns `col` and `partner`.
///
/// Tests the determinant condition from Algorithm 4.1:
/// `|det(D_22)| >= 0.5 * |a_21|^2`
#[cfg(test)]
#[allow(dead_code)] // Used by complete_pivoting_factor (also #[cfg(test)])
fn try_2x2_pivot(
    mut a: MatMut<'_, f64>,
    col: usize,
    partner: usize,
    threshold: f64,
    small: f64,
) -> Result<(Block2x2, f64), ()> {
    let m = a.nrows();
    let a11 = a[(col, col)];
    let a22 = a[(partner, partner)];
    let (r21, c21) = if partner > col {
        (partner, col)
    } else {
        (col, partner)
    };
    let a21 = a[(r21, c21)];

    let det = a11 * a22 - a21 * a21;
    if det.abs() < 0.5 * a21 * a21 {
        return Err(());
    }
    if det.abs() < small {
        return Err(());
    }

    let block = Block2x2 {
        first_col: col,
        a: a11,
        b: a21,
        c: a22,
    };

    let inv_det = 1.0 / det;
    let stability_bound = 1.0 / threshold;
    let mut max_l: f64 = 0.0;

    let rows_start = col.max(partner) + 1;
    let n_rows = m - rows_start;
    // TODO(Phase 9.1): Replace per-column Vec backup with per-block backup
    // when two-level blocking is implemented.
    let mut backup_col: Vec<f64> = Vec::with_capacity(n_rows);
    let mut backup_partner: Vec<f64> = Vec::with_capacity(n_rows);
    for i in rows_start..m {
        backup_col.push(a[(i, col)]);
        backup_partner.push(a[(i, partner)]);
    }

    for i in rows_start..m {
        let ai1 = a[(i, col)];
        let ai2 = a[(i, partner)];

        let l_i1 = (ai1 * a22 - ai2 * a21) * inv_det;
        let l_i2 = (ai2 * a11 - ai1 * a21) * inv_det;

        a[(i, col)] = l_i1;
        a[(i, partner)] = l_i2;

        let abs1 = l_i1.abs();
        let abs2 = l_i2.abs();
        if abs1 > max_l {
            max_l = abs1;
        }
        if abs2 > max_l {
            max_l = abs2;
        }
    }

    if max_l >= stability_bound {
        for (idx, i) in (rows_start..m).enumerate() {
            a[(i, col)] = backup_col[idx];
            a[(i, partner)] = backup_partner[idx];
        }
        return Err(());
    }

    Ok((block, max_l))
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
    let mut stats = AptpStatistics {
        num_1x1: 0,
        num_2x2: 0,
        num_delayed: 0,
        max_l_entry: 0.0,
    };
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

            if delta.abs() >= 0.5 * a_tm * a_tm {
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

    // Build final D of correct dimension
    let mut final_d = MixedDiagonal::new(num_eliminated);
    let mut col = 0;
    while col < num_eliminated {
        match d.get_pivot_type(col) {
            PivotType::OneByOne => {
                final_d.set_1x1(col, d.get_1x1(col));
                col += 1;
            }
            PivotType::TwoByTwo { partner } if partner > col => {
                let block = d.get_2x2(col);
                final_d.set_2x2(block);
                col += 2;
            }
            _ => {
                col += 1;
            }
        }
    }

    let delayed_cols: Vec<usize> = (num_eliminated..n).map(|i| col_order[i]).collect();

    AptpFactorResult {
        d: final_d,
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
    let mut stats = AptpStatistics {
        num_1x1: 0,
        num_2x2: 0,
        num_delayed: 0,
        max_l_entry: 0.0,
    };
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

            if delta.abs() >= 0.5 * a_tm * a_tm {
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
struct BlockBackup {
    data: Mat<f64>,
}

impl BlockBackup {
    /// Create a backup of the block column starting at `col_start` with `block_cols` columns.
    /// Backs up `a[col_start.., col_start..col_start+block_cols]`.
    fn create(a: MatRef<'_, f64>, col_start: usize, block_cols: usize, m: usize) -> Self {
        let rows = m - col_start;
        let mut data = Mat::zeros(rows, block_cols);
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
fn apply_and_check(
    mut a: MatMut<'_, f64>,
    col_start: usize,
    block_nelim: usize,
    block_cols: usize,
    m: usize,
    d: &MixedDiagonal,
    threshold: f64,
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

    let l11_copy = a.rb().submatrix(col_start, col_start, block_nelim, block_nelim).to_owned();
    let panel = a.rb_mut().submatrix_mut(panel_start, col_start, panel_rows, block_nelim);
    triangular_solve::solve_unit_lower_triangular_in_place(
        l11_copy.as_ref(),
        panel.transpose_mut(),
        Par::Seq,
    );

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
fn update_trailing(
    a: MatMut<'_, f64>,
    col_start: usize,
    nelim: usize,
    block_cols: usize,
    m: usize,
    d: &MixedDiagonal,
) {
    if nelim == 0 {
        return;
    }

    let trailing_start = col_start + block_cols;
    let trailing_size = m - trailing_start;
    if trailing_size == 0 {
        return;
    }

    // Compute W = L21 * D11 (trailing_size × nelim workspace)
    let mut w = Mat::zeros(trailing_size, nelim);
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

    // GEMM: A22 -= W * L21^T (lower triangle only)
    // Copy L21 to avoid borrow conflict (L21 and A22 overlap in `a`)
    let l21_copy = a.rb().submatrix(trailing_start, col_start, trailing_size, nelim).to_owned();
    let mut a22 = a.submatrix_mut(trailing_start, trailing_start, trailing_size, trailing_size);

    tri_matmul::matmul_with_conj(
        a22.rb_mut(),
        BlockStructure::TriangularLower,
        Accum::Add,
        w.as_ref(),
        BlockStructure::Rectangular,
        Conj::No,
        l21_copy.as_ref().transpose(),
        BlockStructure::Rectangular,
        Conj::No,
        -1.0,
        Par::Seq,
    );
}

/// Apply updates from newly-factored block to previously-delayed columns.
///
/// Corresponds to UpdateNT/UpdateTN from Algorithm 3.1 (Duff et al. 2020).
/// For each previously-delayed column region, applies the rank-nelim
/// update using the current block's L and D factors.
///
/// delayed_cols: slice of column indices (absolute positions in matrix) that
/// were delayed from earlier blocks and need updating from this block's factors.
#[allow(dead_code)]
fn update_delayed(
    mut a: MatMut<'_, f64>,
    col_start: usize,
    nelim: usize,
    block_cols: usize,
    m: usize,
    delayed_cols: &[usize],
    d: &MixedDiagonal,
) {
    if nelim == 0 || delayed_cols.is_empty() {
        return;
    }

    let panel_start = col_start + block_cols;
    let panel_rows = m - panel_start;

    // For each delayed column d_col, update:
    // a[d_col, d_col2] -= sum_k L_factor[d_col, k] * W[d_col2, k]
    // where L_factor is the L21 panel from this block's factorization
    // and W = L21 * D11
    //
    // But delayed columns are at positions < col_start (from earlier blocks).
    // The L entries for these positions are in the L21 panel of THIS block.
    // Wait — that's not right. Delayed columns are at positions AFTER the
    // eliminated columns (they were swapped to the end). They are in the
    // current frontal matrix at positions that may be in the "still to process"
    // region.
    //
    // In the two-level algorithm, after block j factors nelim_j columns,
    // the trailing submatrix (including delayed columns from earlier blocks)
    // needs the Schur complement update from the newly-factored columns.
    //
    // update_trailing already handles the trailing submatrix below block_cols.
    // But delayed columns from EARLIER blocks are at positions between the
    // current block's eliminated columns and the trailing submatrix.
    //
    // For now, this is handled by the trailing update (which covers everything
    // after col_start + block_cols). Delayed columns from earlier blocks that
    // are still in the "to be processed" region ARE part of the trailing
    // submatrix and will be updated by update_trailing.
    //
    // The explicit update_delayed is needed when delayed columns from block j
    // are at positions within [col_start..col_start+block_cols] but after nelim.
    // These are "failed columns" from THIS block. They were restored from backup
    // and already have correct values — no further update needed from their own
    // block's factors.
    //
    // In practice, for a sequential single-threaded implementation, the trailing
    // update covers all necessary cases. This function is a placeholder for the
    // parallel case (Phase 8.2) where explicit delayed-column updates are needed.
    let _ = (
        a.rb_mut(),
        col_start,
        nelim,
        block_cols,
        m,
        delayed_cols,
        d,
        panel_start,
        panel_rows,
    );
}

/// Factor an nb-sized block using BLAS-3 Factor/Apply/Update loop.
///
/// This is the middle level of the two-level hierarchy. Processes `num_fully_summed`
/// columns of the block `a[0..m, 0..m]` using ib-sized sub-blocks with the
/// three-phase BLAS-3 pattern from SPRAL's `run_elim_pivoted`:
///
/// 1. **Backup**: Save `a[k..m, k..k+block_size]` before factoring
/// 2. **Factor**: `factor_block_diagonal` on the ib×ib diagonal block (complete pivoting)
/// 3. **Apply**: `apply_and_check` — TRSM on panel + threshold check → effective_nelim
/// 4. **Adjust**: `adjust_for_2x2_boundary` — avoid splitting 2×2 across boundaries
/// 5. On failure: full restore from backup, swap failed columns to end, retry
/// 6. **Update**: `update_trailing` — GEMM A22 -= L21*D*L21^T
///
/// # Arguments
/// - `a`: Dense frontal matrix block (m × m), modified in place
/// - `num_fully_summed`: Number of columns eligible for elimination
/// - `options`: APTP configuration (inner_block_size determines ib)
///
/// # References
/// - SPRAL: `run_elim_pivoted` in `ldlt_app.cxx:1273-1579`
/// - Duff, Hogg & Lopez (2020), Algorithm 3.1
fn factor_inner(
    mut a: MatMut<'_, f64>,
    num_fully_summed: usize,
    options: &AptpOptions,
) -> Result<AptpFactorResult, SparseError> {
    let m = a.nrows();
    let ib = options.inner_block_size;
    let small = options.small;
    let threshold = options.threshold;
    let p = num_fully_summed;

    let mut col_order: Vec<usize> = (0..m).collect();
    let mut d = MixedDiagonal::new(p);
    let mut stats = AptpStatistics {
        num_1x1: 0,
        num_2x2: 0,
        num_delayed: 0,
        max_l_entry: 0.0,
    };
    let mut pivot_log = Vec::with_capacity(p);
    let mut k = 0;
    let mut end_pos = p;

    // BLAS-3 Factor/Apply/Update loop with ib-sized inner blocks.
    //
    // Each iteration processes min(end_pos - k, ib) columns:
    //   1. Backup: save a[k..m, k..k+bs]
    //   2. Factor: factor_block_diagonal with block-scoped swaps (row_limit = k+bs)
    //   3. Permute panel column entries by block_perm
    //   4. Zero D off-diagonals for TRSM
    //   5. Apply+Check: TRSM on panel + threshold check
    //   6. Adjust: avoid splitting 2×2 across boundary
    //   7. On fail: restore from backup (columns 0..k untouched), delay, continue
    //   8. On pass: propagate row perm to columns 0..k
    //   9. Update trailing: GEMM A22 -= L21*D*L21^T
    //
    // Block-scoped swaps (step 2) only permute rows/columns within [0..k+bs],
    // leaving panel rows untouched. This matches SPRAL's block_ldlt architecture
    // where panel permutation happens separately in the Apply phase.
    //
    // # References
    // - SPRAL: `block_ldlt::swap_cols` limited to BLOCK_SIZE rows
    // - SPRAL: `apply_rperm_and_backup` / `apply_cperm_and_backup` for panel permutation
    while k < end_pos {
        let block_size = (end_pos - k).min(ib);
        let block_end = k + block_size;

        // 1. BACKUP: save a[k..m, k..k+block_size] before factoring
        let backup = BlockBackup::create(a.as_ref(), k, block_size, m);

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
        //    Block-scoped swaps didn't touch panel rows, so we must gather
        //    the panel entries into the factored column order.
        //    This is SPRAL's `apply_cperm_and_backup` equivalent.
        let panel_start = block_end;
        for r in panel_start..m {
            let orig: Vec<f64> = (0..block_size).map(|i| a[(r, k + block_perm[i])]).collect();
            for i in 0..block_size {
                a[(r, k + i)] = orig[i];
            }
        }

        // 4. Zero out D off-diagonals so apply_and_check's TRSM reads them as
        //    L11 entries (should be 0 for 2×2 pivots where L starts at row k+2).
        //    The D values are safely stored in block_d (MixedDiagonal).
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

        // 5. APPLY: TRSM on panel below + threshold check
        let mut effective_nelim = apply_and_check(
            a.rb_mut(),
            k,
            block_nelim,
            block_size,
            m,
            &block_d,
            threshold,
        );

        // 6. ADJUST: don't split 2×2 pivot across block boundary
        effective_nelim = adjust_for_2x2_boundary(effective_nelim, &block_d);

        if effective_nelim < block_nelim {
            // 7. Threshold failure: full restore, delay failed columns, retry.
            //    Because block-scoped swaps didn't touch columns 0..k, the restore
            //    is clean — no need to un-permute previously-factored columns.

            // a. Full restore from backup (undo factor + apply + panel permute)
            backup.restore_failed(a.rb_mut(), k, 0, block_size, m);

            // b. Identify failed columns via block_perm and delay them.
            //    Failed columns are at block_perm positions effective_nelim..block_nelim.
            //    We also delay any columns that were singular (block_nelim..block_size).
            let n_failed = block_size - effective_nelim;

            // Collect the absolute positions of failed columns (reverse-sorted
            // to swap from end first, avoiding positional conflicts)
            let mut failed_positions: Vec<usize> = (effective_nelim..block_size)
                .map(|i| k + block_perm[i])
                .collect();
            failed_positions.sort_unstable();
            failed_positions.reverse();

            for &failed_abs in &failed_positions {
                end_pos -= 1;
                if failed_abs < end_pos {
                    swap_symmetric(a.rb_mut(), failed_abs, end_pos);
                    col_order.swap(failed_abs, end_pos);
                }
                stats.num_delayed += 1;
                pivot_log.push(AptpPivotRecord {
                    col: col_order[end_pos],
                    pivot_type: PivotType::Delayed,
                    max_l_entry: 0.0,
                    was_fallback: false,
                });
            }

            let _ = n_failed; // used for clarity above
            continue;
        }

        // 8. SUCCESS: propagate row permutation to columns 0..k.
        //    Block-scoped swaps permuted rows within [k..block_end] according to
        //    block_perm, but columns 0..k still have L entries in the old row order.
        //    We must apply the same permutation so extract_l reads consistent rows.
        if k > 0 {
            let mut temp = vec![0.0f64; block_size];
            for c in 0..k {
                // Gather: temp[i] = row that ended up at local position i
                for i in 0..block_size {
                    temp[i] = a[(k + block_perm[i], c)];
                }
                // Scatter: write to sequential positions
                for i in 0..block_size {
                    a[(k + i, c)] = temp[i];
                }
            }
        }

        // Apply the block's permutation to col_order.
        let orig_order: Vec<usize> = col_order[k..k + block_size].to_vec();
        for i in 0..block_size {
            col_order[k + i] = orig_order[block_perm[i]];
        }

        // 9. UPDATE: GEMM A22 -= L21 * D11 * L21^T
        update_trailing(a.rb_mut(), k, block_nelim, block_size, m, &block_d);

        // 10. Accumulate D entries (remap from block-local to global positions)
        let mut bcol = 0;
        while bcol < block_nelim {
            match block_d.get_pivot_type(bcol) {
                PivotType::OneByOne => {
                    d.set_1x1(k + bcol, block_d.get_1x1(bcol));
                    bcol += 1;
                }
                PivotType::TwoByTwo { partner } if partner > bcol => {
                    let blk = block_d.get_2x2(bcol);
                    d.set_2x2(Block2x2 {
                        first_col: k + bcol,
                        a: blk.a,
                        b: blk.b,
                        c: blk.c,
                    });
                    bcol += 2;
                }
                _ => {
                    bcol += 1;
                }
            }
        }

        // Accumulate stats
        stats.num_1x1 += block_stats.num_1x1;
        stats.num_2x2 += block_stats.num_2x2;
        if block_stats.max_l_entry > stats.max_l_entry {
            stats.max_l_entry = block_stats.max_l_entry;
        }
        // Also check panel L entries for max_l_entry
        for c in 0..block_nelim {
            for i in panel_start..m {
                let v = a[(i, k + c)].abs();
                if v > stats.max_l_entry {
                    stats.max_l_entry = v;
                }
            }
        }

        // Accumulate pivot log (skip delayed — we handled delays above)
        for record in &block_log {
            if !matches!(record.pivot_type, PivotType::Delayed) {
                // Remap col indices from block-local to original via col_order
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

        // Handle any columns that were singular within the block (block_nelim < block_size)
        if block_nelim < block_size {
            let n_block_delayed = block_size - block_nelim;
            for i in 0..n_block_delayed {
                let delayed_pos = k + block_nelim + i;
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

        k += block_nelim;
    }

    let num_eliminated = k;

    // Build final D
    let mut final_d = MixedDiagonal::new(num_eliminated);
    let mut col = 0;
    while col < num_eliminated {
        match d.get_pivot_type(col) {
            PivotType::OneByOne => {
                final_d.set_1x1(col, d.get_1x1(col));
                col += 1;
            }
            PivotType::TwoByTwo { partner } if partner > col => {
                let block = d.get_2x2(col);
                final_d.set_2x2(block);
                col += 2;
            }
            _ => {
                col += 1;
            }
        }
    }

    let delayed_cols: Vec<usize> = (num_eliminated..p).map(|i| col_order[i]).collect();

    Ok(AptpFactorResult {
        d: final_d,
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
) -> Result<AptpFactorResult, SparseError> {
    let m = a.nrows();
    let nb = options.outer_block_size;
    let p = num_fully_summed;

    let mut col_order: Vec<usize> = (0..m).collect();
    let mut global_d = MixedDiagonal::new(p);
    let mut stats = AptpStatistics {
        num_1x1: 0,
        num_2x2: 0,
        num_delayed: 0,
        max_l_entry: 0.0,
    };
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
            factor_inner(block_view, block_cols, options)?
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
        // Copy D entries from block_result
        let mut bcol = 0;
        while bcol < block_nelim {
            match block_result.d.get_pivot_type(bcol) {
                PivotType::OneByOne => {
                    global_d.set_1x1(global_nelim + bcol, block_result.d.get_1x1(bcol));
                    bcol += 1;
                }
                PivotType::TwoByTwo { partner } if partner > bcol => {
                    let block = block_result.d.get_2x2(bcol);
                    global_d.set_2x2(Block2x2 {
                        first_col: global_nelim + bcol,
                        a: block.a,
                        b: block.b,
                        c: block.c,
                    });
                    bcol += 2;
                }
                _ => {
                    bcol += 1;
                }
            }
        }

        // Update col_order: factor_inner may have permuted columns within the block
        let block_perm = &block_result.perm;
        let orig_order: Vec<usize> = col_order[col_start..col_start + block_cols].to_vec();
        for i in 0..block_cols {
            if block_perm[i] < block_cols {
                col_order[col_start + i] = orig_order[block_perm[i]];
            }
        }

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

    // Build final D of correct dimension
    let mut final_d = MixedDiagonal::new(global_nelim);
    let mut col = 0;
    while col < global_nelim {
        match global_d.get_pivot_type(col) {
            PivotType::OneByOne => {
                final_d.set_1x1(col, global_d.get_1x1(col));
                col += 1;
            }
            PivotType::TwoByTwo { partner } if partner > col => {
                let block = global_d.get_2x2(col);
                final_d.set_2x2(block);
                col += 2;
            }
            _ => {
                col += 1;
            }
        }
    }

    let delayed_cols: Vec<usize> = (global_nelim..p).map(|i| col_order[i]).collect();

    Ok(AptpFactorResult {
        d: final_d,
        perm: col_order,
        num_eliminated: global_nelim,
        delayed_cols,
        stats,
        pivot_log,
    })
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
            stats.max_l_entry <= 4.0 + 1e-10,
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
            stats.max_l_entry <= 4.0 + 1e-10,
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
                stats.max_l_entry <= 4.0 + 1e-10,
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

        assert_eq!(result.stats.num_1x1, 2);
        assert_eq!(result.stats.num_2x2, 0);
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

        assert_eq!(result.stats.num_1x1, 3);
        assert_eq!(result.stats.num_2x2, 0);
        assert_eq!(result.stats.num_delayed, 0);

        let error =
            dense_reconstruction_error(&a, &result.l, &result.d, result.perm.as_ref().arrays().0);
        assert!(error < 1e-12, "reconstruction error {:.2e} >= 1e-12", error);
    }

    #[test]
    fn test_all_delayed_zero_matrix() {
        let n = 4;
        let a = Mat::zeros(n, n);

        let opts = AptpOptions::default();
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        assert_eq!(result.stats.num_1x1, 0);
        assert_eq!(result.stats.num_2x2, 0);
        assert_eq!(result.stats.num_delayed, n);
        assert_eq!(result.delayed_cols.len(), n);
    }

    #[test]
    fn test_1x1_singularity_detection() {
        let a = symmetric_matrix(3, |i, j| if i == j { [4.0, 1e-25, 9.0][i] } else { 0.0 });

        let opts = AptpOptions {
            fallback: AptpFallback::Delay,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        assert_eq!(result.stats.num_delayed, 1);
        assert_eq!(result.stats.num_1x1, 2);
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

        assert_eq!(result.stats.num_1x1, 4);
        assert_eq!(result.stats.num_2x2, 0);
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
            result.stats.max_l_entry <= 4.0 + 1e-10,
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
                    assert_eq!(
                        result.stats.num_2x2, 0,
                        "PD matrix {}x{} should have zero 2x2 pivots",
                        n, n
                    );
                    assert_eq!(result.stats.num_1x1, n);

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
            ..AptpOptions::default()
        };
        let mut a_copy = a.clone();
        let result = aptp_factor_in_place(a_copy.as_mut(), 8, &opts).unwrap();

        assert_eq!(result.num_eliminated, 0);
        assert_eq!(result.stats.num_delayed, 8);
        assert_eq!(result.delayed_cols.len(), 8);
    }
}
