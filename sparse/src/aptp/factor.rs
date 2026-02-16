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

use faer::perm::Perm;
use faer::prelude::*;
use faer::{Mat, MatMut, MatRef};

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

    let p = num_fully_summed;
    let threshold = options.threshold;
    let small = options.small;

    // col_order[i] = original column index at working position i.
    // This becomes the output permutation.
    let mut col_order: Vec<usize> = (0..m).collect();
    let mut d = MixedDiagonal::new(p);
    let mut stats = AptpStatistics {
        num_1x1: 0,
        num_2x2: 0,
        num_delayed: 0,
        max_l_entry: 0.0,
    };
    let mut pivot_log = Vec::with_capacity(p);

    // Boundary: positions 0..end_pos are unprocessed fully-summed columns.
    // Eliminated columns accumulate at the front (0..k), delayed columns
    // are swapped to positions end_pos..p.
    let mut end_pos = p;
    let mut k = 0;

    while k < end_pos {
        let pivot_result = try_1x1_pivot(a.rb_mut(), k, threshold, small);

        match pivot_result {
            Ok((d_value, max_l)) => {
                // 1x1 pivot accepted at position k
                d.set_1x1(k, d_value);
                update_schur_1x1(a.rb_mut(), k, d_value);

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
            }
            Err(max_l_1x1) => {
                // 1x1 failed — try 2x2 BK fallback if configured
                let bk_accepted = if options.fallback == AptpFallback::BunchKaufman {
                    select_2x2_partner_range(a.rb(), k, k + 1, end_pos).and_then(|partner_pos| {
                        if partner_pos != k + 1 {
                            swap_symmetric(a.rb_mut(), k + 1, partner_pos);
                            col_order.swap(k + 1, partner_pos);
                        }
                        try_2x2_pivot(a.rb_mut(), k, k + 1, threshold, small).ok()
                    })
                } else {
                    None
                };

                if let Some((block, max_l_2x2)) = bk_accepted {
                    // 2x2 pivot accepted
                    let d_block = Block2x2 {
                        first_col: k,
                        a: block.a,
                        b: block.b,
                        c: block.c,
                    };
                    d.set_2x2(d_block);
                    update_schur_2x2(a.rb_mut(), k, k + 1, &d_block);

                    stats.num_2x2 += 1;
                    if max_l_2x2 > stats.max_l_entry {
                        stats.max_l_entry = max_l_2x2;
                    }
                    pivot_log.push(AptpPivotRecord {
                        col: col_order[k],
                        pivot_type: PivotType::TwoByTwo {
                            partner: col_order[k + 1],
                        },
                        max_l_entry: max_l_2x2,
                        was_fallback: true,
                    });
                    pivot_log.push(AptpPivotRecord {
                        col: col_order[k + 1],
                        pivot_type: PivotType::TwoByTwo {
                            partner: col_order[k],
                        },
                        max_l_entry: max_l_2x2,
                        was_fallback: true,
                    });
                    k += 2;
                } else {
                    // Delay column k (unified path for BK-failed and Delay fallback)
                    let was_fallback = options.fallback == AptpFallback::BunchKaufman;
                    stats.num_delayed += 1;
                    end_pos -= 1;
                    pivot_log.push(AptpPivotRecord {
                        col: col_order[k],
                        pivot_type: PivotType::Delayed,
                        max_l_entry: max_l_1x1,
                        was_fallback,
                    });
                    if k != end_pos {
                        swap_symmetric(a.rb_mut(), k, end_pos);
                        col_order.swap(k, end_pos);
                    }
                    // Don't increment k: retry with new column at this position
                }
            }
        }
    }

    let num_eliminated = k; // == end_pos

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
    if i == j {
        return;
    }
    let (i, j) = if i < j { (i, j) } else { (j, i) };
    let m = a.nrows();

    // Swap diagonals
    let tmp = a[(i, i)];
    a[(i, i)] = a[(j, j)];
    a[(j, j)] = tmp;

    // Rows k < i: swap lower-triangle entries a[(i,k)] and a[(j,k)]
    for k in 0..i {
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

    // Rows k > j: swap a[(k,i)] and a[(k,j)]
    for k in (j + 1)..m {
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

/// Select the best partner column for a 2x2 pivot at column `col`.
///
/// Searches positions `range_start..range_end` for the column with the
/// largest absolute off-diagonal entry with column `col`.
fn select_2x2_partner_range(
    a: MatRef<'_, f64>,
    col: usize,
    range_start: usize,
    range_end: usize,
) -> Option<usize> {
    let mut best_idx = None;
    let mut best_val = 0.0_f64;

    for j in range_start..range_end {
        let (r, c) = if j > col { (j, col) } else { (col, j) };
        let val = a[(r, c)].abs();
        if val > best_val {
            best_val = val;
            best_idx = Some(j);
        }
    }

    best_idx
}

/// Attempt a 2x2 Bunch-Kaufman pivot using columns `col` and `partner`.
///
/// Tests the determinant condition from Algorithm 4.1:
/// `|det(D_22)| >= 0.5 * |a_21|^2`
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
fn complete_pivoting_factor(
    mut a: MatMut<'_, f64>,
    small: f64,
) -> AptpFactorResult {
    let n = a.nrows();
    debug_assert_eq!(n, a.ncols(), "complete_pivoting_factor requires square matrix");

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
            for i in k..n {
                stats.num_delayed += 1;
                pivot_log.push(AptpPivotRecord {
                    col: col_order[i],
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
                    a[(i, k)] = a[(i, k)] * inv_d;
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
        let a = symmetric_matrix(3, |i, j| {
            if i == j {
                [2.0, 5.0, 3.0][i]
            } else {
                0.0
            }
        });
        let mut a_copy = a.clone();
        let result = complete_pivoting_factor(a_copy.as_mut(), 1e-20);

        assert_eq!(result.num_eliminated, 3);
        assert_eq!(result.stats.num_1x1, 3);

        // First pivot should be 5.0 (col 1), then 3.0 (col 2), then 2.0 (col 0)
        assert!((result.d.get_1x1(0) - 5.0).abs() < 1e-14, "first pivot should be 5.0, got {}", result.d.get_1x1(0));
        assert!((result.d.get_1x1(1) - 3.0).abs() < 1e-14, "second pivot should be 3.0, got {}", result.d.get_1x1(1));
        assert!((result.d.get_1x1(2) - 2.0).abs() < 1e-14, "third pivot should be 2.0, got {}", result.d.get_1x1(2));
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

        assert!(stats.num_2x2 >= 1, "expected at least one 2×2 pivot, got {}", stats.num_2x2);
        // L entries bounded by 4 (complete pivoting guarantee)
        assert!(stats.max_l_entry <= 4.0 + 1e-10, "L entries should be bounded by 4, got {}", stats.max_l_entry);

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
        assert!(stats.max_l_entry <= 4.0 + 1e-10, "L entries should be bounded by 4");

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
        assert_eq!(result.num_eliminated, 0, "near-singular block should have 0 eliminations");
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
                stats.num_1x1 + 2 * stats.num_2x2 + stats.num_delayed, n,
                "statistics invariant for {}x{}", n, n
            );

            if stats.num_delayed == 0 {
                let error = dense_reconstruction_error(&a, &l, &d, &perm);
                assert!(
                    error < 1e-12,
                    "complete pivoting {}x{}: reconstruction error {:.2e} >= 1e-12",
                    n, n, error
                );
            }

            // L entries bounded by 4 (Algorithm 4.1 guarantee)
            assert!(
                stats.max_l_entry <= 4.0 + 1e-10,
                "complete pivoting {}x{}: max_l_entry {:.2e} > 4",
                n, n, stats.max_l_entry
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
        // l_21 = 1/1e-4 = 10000 > 100 (bound for threshold=0.01)
        let a = symmetric_matrix(2, |i, j| {
            let vals = [[1e-4, 1.0], [1.0, 1.0]];
            vals[i][j]
        });

        let opts = AptpOptions {
            fallback: AptpFallback::Delay,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        assert!(result.stats.num_delayed >= 1);
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
        // Matrix designed so 1x1 fails AND 2x2 also fails → column delayed.
        //
        // 1x1 at col 0: a00 = 0.001, l_10 = 0.11/0.001 = 110 > 1/0.01 = 100 → fails
        // 2x2 at (0,1): det = 0.001*12.0 - 0.11^2 = 0.012 - 0.0121 = -0.0001
        //   0.5 * a21^2 = 0.5 * 0.0121 = 0.00605
        //   |det| = 0.0001 < 0.00605 → determinant condition fails, 2x2 rejected
        let a = symmetric_matrix(3, |i, j| {
            let vals = [[0.001, 0.11, 0.0], [0.11, 12.0, 0.1], [0.0, 0.1, 5.0]];
            vals[i][j]
        });

        let opts = AptpOptions {
            fallback: AptpFallback::BunchKaufman,
            ..AptpOptions::default()
        };
        let result = aptp_factor(a.as_ref(), &opts).unwrap();

        // Column 0 should be delayed (both 1x1 and 2x2 failed)
        assert!(
            result.stats.num_delayed >= 1,
            "expected at least 1 delay, got {}",
            result.stats.num_delayed
        );
        assert_eq!(result.stats.num_2x2, 0, "2x2 should have been rejected");

        // Verify the delayed column's pivot log shows was_fallback = true
        let delayed_entry = result
            .pivot_log
            .iter()
            .find(|p| matches!(p.pivot_type, PivotType::Delayed))
            .expect("should have a delayed pivot log entry");
        assert!(
            delayed_entry.was_fallback,
            "delayed after BK attempt should have was_fallback=true"
        );

        // Statistics invariant still holds
        let sum = result.stats.num_1x1 + 2 * result.stats.num_2x2 + result.stats.num_delayed;
        assert_eq!(sum, 3, "statistics sum {} != n=3", sum);
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
}
