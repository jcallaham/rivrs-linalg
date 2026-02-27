# API Contract: AptpSymbolic

**Feature**: 010-aptp-symbolic
**Date**: 2026-02-10

## Overview

AptpSymbolic provides the symbolic analysis (Phase 1 of analyze→factorize→solve) for APTP-based sparse indefinite LDL^T factorization. It composes faer's `SymbolicCholesky` with APTP-specific metadata.

## Public API

### Construction

#### `AptpSymbolic::analyze`

Perform symbolic analysis on a sparse symmetric matrix.

**Signature** (notional):
```
analyze(matrix, ordering) → Result<AptpSymbolic, SparseError>
```

**Parameters**:
- `matrix` — Symbolic sparsity pattern of the sparse symmetric matrix (CSC format, no values needed). Type: faer's `SymbolicSparseColMatRef<'_, usize>`.
- `ordering` — Fill-reducing ordering strategy. Type: faer's `SymmetricOrdering<'_, usize>`. Default: `Amd`.

**Returns**: `Ok(AptpSymbolic)` on success, `Err(SparseError)` on failure.

**Errors**:
- `SparseError::DimensionMismatch` — custom permutation dimension does not match matrix dimension
- `SparseError::NotSquare` — matrix is not square
- `SparseError::AnalysisFailure` — internal symbolic factorization failed

**Behavior**:
1. Validate input (square matrix, permutation dimension if custom ordering)
2. Compute elimination tree and column counts via faer's `prefactorize_symbolic_cholesky`
3. Compute full symbolic factorization via faer's `factorize_symbolic_cholesky`
4. Compute APTP pivot buffer estimates from column counts
5. Return composed result

**Determinism**: Same inputs always produce identical outputs.

---

### Accessors (delegated to faer)

#### `perm`
```
perm(&self) → Option<PermRef<'_, usize>>
```
Returns the fill-reducing permutation, if one was computed. `None` if `SymmetricOrdering::Identity` was used.

#### `predicted_nnz`
```
predicted_nnz(&self) → usize
```
Returns the predicted number of nonzeros in the L factor (delegated to `inner.len_val()`).

#### `nrows` / `ncols`
```
nrows(&self) → usize
ncols(&self) → usize
```
Returns the matrix dimension.

#### `raw`
```
raw(&self) → &SymbolicCholeskyRaw<usize>
```
Returns a reference to the inner symbolic structure (simplicial or supernodal variant). Enables downstream code to access faer's full symbolic API.

---

### APTP-Specific Accessors

#### `etree`
```
etree(&self) → &[isize]
```
Returns the elimination tree as a parent pointer array. `etree[j]` is the parent of column `j`, or `-1` if `j` is a root.

#### `col_counts`
```
col_counts(&self) → &[usize]
```
Returns the predicted nonzero count per column of L.

#### `pivot_buffer_estimates`
```
pivot_buffer_estimates(&self) → &[usize]
```
Returns per-supernode (supernodal mode) or per-column (simplicial mode) extra space estimates for delayed pivots.

#### `total_pivot_buffer`
```
total_pivot_buffer(&self) → usize
```
Returns the sum of all pivot buffer estimates.

---

### Diagnostics

#### `statistics`
```
statistics(&self) → SymbolicStatistics
```
Returns a summary of the symbolic analysis for diagnostics and debugging.

#### `is_supernodal`
```
is_supernodal(&self) → bool
```
Returns whether faer chose supernodal (vs simplicial) for this matrix.

---

### Supernodal Accessors (available when supernodal)

#### `n_supernodes`
```
n_supernodes(&self) → Option<usize>
```
Returns the number of supernodes if the result is supernodal, `None` if simplicial.

#### `supernode_begin` / `supernode_end`
```
supernode_begin(&self) → Option<&[usize]>
supernode_end(&self) → Option<&[usize]>
```
Returns supernode column range boundaries (supernodal mode only).

---

## SymbolicStatistics

```
struct SymbolicStatistics {
    dimension: usize,
    predicted_nnz: usize,
    average_col_count: f64,
    is_supernodal: bool,
    n_supernodes: Option<usize>,
    total_pivot_buffer: usize,
}
```

Implements `Display` for human-readable output and `Debug` for logging.

---

## Error Contract

| Condition | Error Variant | Message Contains |
|-----------|---------------|------------------|
| Non-square matrix | `SparseError::NotSquare` | matrix dimensions |
| Custom perm wrong size | `SparseError::DimensionMismatch` | expected vs actual size |
| faer internal failure | `SparseError::AnalysisFailure` | faer error message |
| 0×0 matrix | Succeeds | (trivial result) |

---

## Thread Safety

`AptpSymbolic` is `Send + Sync` (all fields are owned data or immutable references). It can be shared across threads for parallel numeric factorizations of different value sets with the same sparsity pattern.

---

## Integration with Downstream Phases

Phase 5-6 (numeric factorization) will call:
```
symbolic.inner.factorize_numeric_ldlt(...)  // or intranode_lblt
```
using the stored `SymbolicCholesky` to drive the numeric factorization, plus `pivot_buffer_estimates()` to pre-allocate APTP-specific workspace.

Phase 4 (MC64) will produce permutations passed as:
```
AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Custom(mc64_perm))
```
