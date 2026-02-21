# API Surface: Two-Level APTP Factorization

**Feature**: 017-two-level-aptp
**Date**: 2026-02-16

## Public API Changes

### Modified: `AptpOptions` (factor.rs)

```rust
/// Configuration for the APTP factorization kernel.
///
/// # Block Size Parameters (Two-Level APTP)
///
/// For frontal matrices larger than `outer_block_size`, the kernel uses
/// a two-level blocked algorithm (Duff, Hogg & Lopez 2020, Section 3):
/// - Outer loop processes blocks of `outer_block_size` columns
/// - Inner loop processes sub-blocks of `inner_block_size` columns
/// - Innermost ib×ib diagonal blocks use complete pivoting (Algorithm 4.1)
///
/// For frontal matrices ≤ `outer_block_size`, the kernel uses the
/// single-level column-by-column algorithm (equivalent to one outer block).
pub struct AptpOptions {
    pub threshold: f64,           // unchanged
    pub small: f64,               // unchanged
    pub fallback: AptpFallback,   // unchanged
    pub outer_block_size: usize,  // NEW: default 256
    pub inner_block_size: usize,  // NEW: default 32
}
```

### Modified: `FactorOptions` (solver.rs)

```rust
/// Options for the numeric factorization phase.
pub struct FactorOptions {
    pub threshold: f64,           // unchanged
    pub fallback: AptpFallback,   // unchanged
    pub outer_block_size: usize,  // NEW: default 256
    pub inner_block_size: usize,  // NEW: default 32
}
```

### Unchanged Public API

The following public interfaces are NOT modified:

| Function/Type | Location | Contract |
|---------------|----------|----------|
| `aptp_factor_in_place(MatMut, usize, &AptpOptions) -> Result<AptpFactorResult>` | factor.rs | Same signature, same return type |
| `aptp_factor(MatRef, &AptpOptions) -> Result<AptpFactorization>` | factor.rs | Same signature, same return type |
| `AptpNumeric::factor(symbolic, matrix, options, scaling) -> Result<Self>` | numeric.rs | Same signature |
| `SparseLDLT::factor(matrix, options) -> Result<()>` | solver.rs | Same signature |
| `SparseLDLT::solve_in_place(rhs, stack) -> Result<()>` | solver.rs | Unchanged |
| `AptpFactorResult` | factor.rs | Same fields |
| `AptpStatistics` | factor.rs | Same fields |
| `MixedDiagonal` | diagonal.rs | Same interface |

## Internal (Private) API Additions

These functions are added to `factor.rs` but NOT exported publicly:

### `complete_pivoting_factor`

```rust
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
/// # SPRAL Equivalent
/// `block_ldlt()` in `spral/src/ssids/cpu/kernels/ldlt_app.cxx`
fn complete_pivoting_factor(
    a: MatMut<'_, f64>,
    small: f64,
) -> AptpFactorResult
```

### `factor_inner`

```rust
/// Factor an nb×nb diagonal block using inner APTP with complete
/// pivoting at the ib×ib leaves.
///
/// This is the middle level of the two-level hierarchy:
/// - Processes ib-sized sub-blocks left-to-right
/// - Uses existing try_1x1_pivot/try_2x2_pivot for columns within
///   each sub-block
/// - Calls complete_pivoting_factor for each ib×ib diagonal sub-block
/// - Uses configured fallback strategy (BunchKaufman/Delay) for
///   column failures within sub-blocks
///
/// For fronts ≤ outer_block_size, this is called once on the entire
/// fully-summed portion (single-block case).
fn factor_inner(
    a: MatMut<'_, f64>,
    num_fully_summed: usize,
    options: &AptpOptions,
) -> Result<AptpFactorResult, SparseError>
```

### `apply_and_check`

```rust
/// Apply factored L11/D11 to the panel below the diagonal block (TRSM),
/// then perform a posteriori threshold check on all L21 entries.
///
/// Returns the effective nelim (≤ block_nelim): the number of columns
/// whose L entries all satisfy |l_ij| < 1/threshold.
///
/// # Algorithm
/// 1. Solve L21 = A21 * (L11 * D11)^{-T} via triangular solve + D scaling
/// 2. Scan L21 column-by-column; find first column j where any |l_ij| > 1/u
/// 3. Return min(block_nelim, j) as effective nelim
fn apply_and_check(
    a: MatMut<'_, f64>,
    col_start: usize,
    block_nelim: usize,
    block_cols: usize,
    m: usize,
    d: &MixedDiagonal,
    threshold: f64,
) -> usize
```

### `update_trailing`

```rust
/// Rank-nelim Schur complement update on the trailing submatrix via GEMM.
///
/// Computes: A[trailing, trailing] -= L21 * D11 * L21^T
/// where L21 is the panel below the current block and D11 is the
/// block diagonal from the Factor phase.
///
/// Uses faer's matmul for BLAS-3 performance.
fn update_trailing(
    a: MatMut<'_, f64>,
    col_start: usize,
    nelim: usize,
    m: usize,
    d: &MixedDiagonal,
)
```

### `update_delayed`

```rust
/// Apply updates from newly-factored block to previously-delayed columns.
///
/// Corresponds to UpdateNT/UpdateTN from Algorithm 3.1 (Duff et al. 2020).
/// For each previously-delayed column region, applies the rank-nelim
/// update using the current block's L and D factors.
fn update_delayed(
    a: MatMut<'_, f64>,
    col_start: usize,
    nelim: usize,
    delayed_ranges: &[(usize, usize)],
    d: &MixedDiagonal,
)
```

### `BlockBackup`

```rust
/// Per-block backup for the two-level APTP algorithm.
///
/// Stores a copy of matrix entries for one outer block column,
/// enabling restore when the a posteriori check reduces nelim.
///
/// # SPRAL Equivalent
/// `CopyBackup<T>` in `spral/src/ssids/cpu/kernels/ldlt_app.cxx`
struct BlockBackup {
    col_start: usize,
    block_cols: usize,
    data: Mat<f64>,
}

impl BlockBackup {
    fn create(a: MatRef<'_, f64>, col_start: usize, block_cols: usize, m: usize) -> Self;
    fn restore_failed(
        &self,
        a: MatMut<'_, f64>,
        col_start: usize,
        nelim: usize,
        block_cols: usize,
        m: usize,
    );
}
```

## Backward Compatibility

- **Source compatible**: Existing code using `AptpOptions::default()` and `FactorOptions::default()` continues to work — new fields have defaults (256, 32).
- **Behavior compatible**: For fronts ≤ 256, behavior is equivalent to current single-level (one outer block processed by inner APTP). For fronts > 256, behavior changes to two-level with BLAS-3 operations — numerically equivalent within tolerance, not bitwise identical.
- **No breaking changes**: All existing public function signatures are unchanged. New fields on existing structs use `..Default::default()` pattern.
