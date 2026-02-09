# API Contract: APTP Data Structures

**Feature**: 009-aptp-data-structures
**Date**: 2026-02-08

## Module: `src/aptp/mod.rs`

New top-level module `aptp` (publicly exported from `lib.rs`). Contains pivot types, mixed diagonal, inertia, and permutation utilities.

### Sub-modules

| Module | Visibility | Contents |
|--------|-----------|----------|
| `pivot` | pub | `PivotType`, `Block2x2` |
| `diagonal` | pub | `MixedDiagonal` |
| `inertia` | pub (re-exported) | `Inertia` (relocated from `io/reference`) |
| `perm` | pub | `perm_from_forward` |

## Public API

### `aptp::pivot::PivotType`

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PivotType {
    /// Standard 1x1 scalar pivot.
    OneByOne,
    /// 2x2 Bunch-Kaufman pivot. `partner` identifies the paired column.
    TwoByTwo { partner: usize },
    /// Column failed APTP stability check; deferred to ancestor node.
    Delayed,
}
```

### `aptp::pivot::Block2x2`

```rust
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Block2x2 {
    /// Index of the first (lower-indexed) column.
    pub first_col: usize,
    /// Top-left element D[i, i].
    pub a: f64,
    /// Off-diagonal element D[i, i+1] = D[i+1, i].
    pub b: f64,
    /// Bottom-right element D[i+1, i+1].
    pub c: f64,
}

impl Block2x2 {
    /// Determinant: ac - b².
    pub fn determinant(&self) -> f64;
    /// Trace: a + c.
    pub fn trace(&self) -> f64;
}
```

### `aptp::diagonal::MixedDiagonal`

```rust
pub struct MixedDiagonal { /* private fields */ }

impl MixedDiagonal {
    // -- Construction --

    /// Create a new MixedDiagonal of dimension n.
    /// All columns start as Delayed (unset).
    pub fn new(n: usize) -> Self;

    /// Set column `col` as a 1x1 pivot with the given diagonal value.
    /// Debug-asserts: col < n, column is currently Delayed.
    pub fn set_1x1(&mut self, col: usize, value: f64);

    /// Set a 2x2 pivot block starting at `block.first_col`.
    /// Marks both `first_col` and `first_col + 1` as TwoByTwo.
    /// Debug-asserts: first_col + 1 < n, both columns currently Delayed.
    pub fn set_2x2(&mut self, block: Block2x2);

    // -- Query --

    /// Matrix dimension.
    pub fn dimension(&self) -> usize;

    /// Pivot type for the given column.
    /// Debug-asserts: col < n.
    pub fn get_pivot_type(&self, col: usize) -> PivotType;

    /// Diagonal value for a 1x1 pivot column.
    /// Debug-asserts: col is OneByOne.
    pub fn get_1x1(&self, col: usize) -> f64;

    /// Block data for a 2x2 pivot (by the lower-indexed column).
    /// Debug-asserts: col is TwoByTwo and owns the block.
    pub fn get_2x2(&self, first_col: usize) -> &Block2x2;

    /// Number of columns still marked as Delayed.
    pub fn num_delayed(&self) -> usize;

    /// Number of 1x1 pivots.
    pub fn num_1x1(&self) -> usize;

    /// Number of 2x2 pivot pairs.
    pub fn num_2x2_pairs(&self) -> usize;

    // -- Operations --

    /// Solve D x = b in place.
    /// Debug-asserts: no delayed columns, no singular pivots (zero 1x1
    /// or zero-determinant 2x2), x.len() == n.
    pub fn solve_in_place(&self, x: &mut [f64]);

    /// Compute eigenvalue sign counts from stored pivots.
    /// Debug-asserts: no delayed columns.
    pub fn compute_inertia(&self) -> Inertia;
}
```

### `aptp::inertia::Inertia` (relocated)

```rust
// Relocated from io/reference.rs — same struct, same derives
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Inertia {
    pub positive: usize,
    pub negative: usize,
    pub zero: usize,
}

impl Inertia {
    pub fn dimension(&self) -> usize;
}
```

### `aptp::perm::perm_from_forward`

```rust
use faer::perm::Perm;
use crate::error::SparseError;

/// Construct a `Perm<usize>` from a forward permutation array,
/// computing the inverse automatically.
///
/// # Errors
/// Returns `SparseError::InvalidInput` if the array contains
/// duplicates or out-of-bounds indices.
pub fn perm_from_forward(fwd: Vec<usize>) -> Result<Perm<usize>, SparseError>;
```

## Re-export Changes

### `io/reference.rs`

Add re-export after relocation:
```rust
// Inertia relocated to aptp module; re-export for backward compatibility
pub use crate::aptp::Inertia;
```

### `lib.rs`

Add new module:
```rust
pub mod aptp;
```

## Error Types

No new error variants needed. `perm_from_forward` uses existing `SparseError::InvalidInput`.

## Debug-Assert Summary

All debug-asserts enforce internal invariants (programmer errors, not runtime conditions):

| Function | Assert | Rationale |
|----------|--------|-----------|
| `set_1x1` | col < n | Bounds check |
| `set_1x1` | column is Delayed | Cannot overwrite a set pivot |
| `set_2x2` | first_col + 1 < n | Block must fit |
| `set_2x2` | both columns Delayed | Cannot overwrite set pivots |
| `get_pivot_type` | col < n | Bounds check |
| `get_1x1` | column is OneByOne | Type safety |
| `get_2x2` | column is TwoByTwo owner | Type safety |
| `solve_in_place` | no Delayed columns | Factorization must be complete |
| `solve_in_place` | no singular pivots | Factorization must be valid |
| `solve_in_place` | x.len() == n | Dimension match |
| `compute_inertia` | no Delayed columns | All pivots must be resolved |
