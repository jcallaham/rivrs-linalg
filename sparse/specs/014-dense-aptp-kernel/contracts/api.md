# API Contracts: Dense APTP Factorization Kernel

**Feature Branch**: `014-dense-aptp-kernel`
**Date**: 2026-02-14

## Module: `aptp::factor`

### Types

```rust
/// Configuration for the APTP factorization kernel.
///
/// # Defaults
/// - `threshold`: 0.01 (growth factor bound of 100)
/// - `small`: 1e-20 (singularity detection)
/// - `fallback`: `AptpFallback::BunchKaufman`
pub struct AptpOptions {
    pub threshold: f64,
    pub small: f64,
    pub fallback: AptpFallback,
}

/// Fallback strategy when a 1x1 pivot fails the stability check.
pub enum AptpFallback {
    /// Attempt 2x2 Bunch-Kaufman pivot; delay if that also fails.
    BunchKaufman,
    /// Immediately delay the column without attempting 2x2.
    Delay,
}

/// Result of in-place APTP factorization.
///
/// The L factor is stored in the lower triangle of the mutated input matrix.
/// Only the first `num_eliminated` columns contain valid L entries.
pub struct AptpFactorResult {
    pub d: MixedDiagonal,
    pub perm: Vec<usize>,
    pub num_eliminated: usize,
    pub delayed_cols: Vec<usize>,
    pub stats: AptpStatistics,
    pub pivot_log: Vec<AptpPivotRecord>,
}

/// Result of convenience APTP factorization (with extracted L).
pub struct AptpFactorization {
    pub l: Mat<f64>,
    pub d: MixedDiagonal,
    pub perm: Perm<usize>,
    pub delayed_cols: Vec<usize>,
    pub stats: AptpStatistics,
    pub pivot_log: Vec<AptpPivotRecord>,
}

/// Summary statistics from factorization.
pub struct AptpStatistics {
    pub num_1x1: usize,
    pub num_2x2: usize,
    pub num_delayed: usize,
    pub max_l_entry: f64,
}

/// Per-column pivot diagnostic record.
pub struct AptpPivotRecord {
    pub col: usize,
    pub pivot_type: PivotType,
    pub max_l_entry: f64,
    pub was_fallback: bool,
}
```

### Functions

#### `aptp_factor_in_place` — Core In-Place Factorization

```rust
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
/// Implements the APTP strategy from Duff, Hogg & Lopez (2020),
/// "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting",
/// SIAM J. Sci. Comput. 42(4). Single-level (column-by-column) variant.
///
/// # SPRAL Equivalent
/// Corresponds to `LDLT<T, BLOCK_SIZE, Backup>::factor()` in
/// `spral/src/ssids/cpu/kernels/ldlt_app.cxx` (BSD-3-Clause).
pub fn aptp_factor_in_place(
    a: MatMut<'_, f64>,
    num_fully_summed: usize,
    options: &AptpOptions,
) -> Result<AptpFactorResult, SparseError>;
```

#### `aptp_factor` — Convenience Wrapper

```rust
/// Factor a dense symmetric matrix using APTP, returning extracted L factor.
///
/// Convenience wrapper around `aptp_factor_in_place` that:
/// 1. Copies the input matrix
/// 2. Calls `aptp_factor_in_place` with `num_fully_summed = n`
/// 3. Extracts L from the mutated copy
/// 4. Returns `AptpFactorization` with all components as owned types
///
/// Suitable for standalone testing and small matrices. For large frontal
/// matrices in multifrontal factorization (Phase 6), use
/// `aptp_factor_in_place` directly to avoid the copy.
///
/// # Arguments
/// - `a`: Dense symmetric matrix (n x n), read-only. Full factorization
///   is performed (all n columns are eligible for elimination).
/// - `options`: APTP configuration.
///
/// # Returns
/// `AptpFactorization` containing L, D, permutation, and diagnostics.
pub fn aptp_factor(
    a: MatRef<'_, f64>,
    options: &AptpOptions,
) -> Result<AptpFactorization, SparseError>;
```

### Internal Functions (private, not part of public API)

```rust
/// Attempt a 1x1 pivot at the current column.
///
/// Computes the L column entries and checks the stability bound.
/// Returns Ok(d_value) if the pivot passes, Err(max_l_entry) if it fails.
fn try_1x1_pivot(
    a: MatMut<'_, f64>,
    col: usize,
    d: &MixedDiagonal,
    num_eliminated_before: usize,
    threshold: f64,
    small: f64,
) -> Result<f64, f64>;

/// Attempt a 2x2 Bunch-Kaufman pivot using columns col and partner.
///
/// Computes the 2x2 block factorization and checks stability.
/// Returns Ok(Block2x2) if accepted, Err(()) if rejected.
fn try_2x2_pivot(
    a: MatMut<'_, f64>,
    col: usize,
    partner: usize,
    d: &MixedDiagonal,
    num_eliminated_before: usize,
    threshold: f64,
    small: f64,
) -> Result<Block2x2, ()>;

/// Select the best partner column for a 2x2 pivot at the given column.
///
/// Searches the remaining un-eliminated columns for the one with the
/// largest off-diagonal entry in the current column, following BK strategy.
fn select_2x2_partner(
    a: MatRef<'_, f64>,
    col: usize,
    eliminated: &[bool],
    num_fully_summed: usize,
) -> Option<usize>;

/// Update the trailing submatrix after eliminating a 1x1 pivot.
///
/// Performs: A[i,j] -= l_i * d_kk * l_j for all i,j > k
fn update_schur_1x1(
    a: MatMut<'_, f64>,
    col: usize,
    d_value: f64,
    num_rows: usize,
);

/// Update the trailing submatrix after eliminating a 2x2 pivot.
///
/// Performs: A[i,j] -= [l_i1 l_i2] * D_22 * [l_j1 l_j2]^T for all i,j > k+1
fn update_schur_2x2(
    a: MatMut<'_, f64>,
    col1: usize,
    col2: usize,
    block: &Block2x2,
    num_rows: usize,
);

/// Extract the L factor from the in-place factorized matrix.
///
/// Reads L entries from the lower triangle of the mutated matrix,
/// respecting the column permutation and num_eliminated.
fn extract_l(
    a: MatRef<'_, f64>,
    perm: &[usize],
    num_eliminated: usize,
) -> Mat<f64>;
```

## Error Conditions

| Error                      | Condition                                              |
|----------------------------|--------------------------------------------------------|
| `SparseError::InvalidInput`| Matrix is not square, or `num_fully_summed > nrows`    |
| `SparseError::InvalidInput`| `options.threshold` not in valid range (0, 1]          |

Note: All-columns-delayed is NOT an error — it is a valid result with `num_eliminated == 0`.

## Usage Examples

### Full factorization (testing)

```rust
use rivrs_sparse::aptp::{aptp_factor, AptpOptions, AptpFallback};
use faer::Mat;

let a = Mat::from_fn(5, 5, |i, j| /* symmetric values */);
let options = AptpOptions::default(); // threshold=0.01, BK fallback

let result = aptp_factor(a.as_ref(), &options)?;

// Verify reconstruction: A = P^T L D L^T P
assert!(result.stats.max_l_entry < 1.0 / options.threshold);
assert_eq!(result.stats.num_1x1 + 2 * result.stats.num_2x2
           + result.stats.num_delayed, 5);
```

### Partial factorization (Phase 6 multifrontal)

```rust
use rivrs_sparse::aptp::{aptp_factor_in_place, AptpOptions};
use faer::Mat;

let mut frontal = Mat::zeros(100, 100); // assembled frontal matrix
let num_fully_summed = 60; // first 60 columns are fully summed

let result = aptp_factor_in_place(
    frontal.as_mut(),
    num_fully_summed,
    &AptpOptions::default(),
)?;

// L is in lower triangle of frontal (first result.num_eliminated cols)
// Contribution block (rows/cols 60..100) is updated with Schur complement
// Delayed columns returned for parent node
let delayed = &result.delayed_cols;
```
