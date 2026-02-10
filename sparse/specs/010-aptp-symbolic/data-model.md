# Data Model: APTP Symbolic Analysis

**Feature**: 010-aptp-symbolic
**Date**: 2026-02-10

## Entities

### AptpSymbolic

The central analysis result. Composes faer's `SymbolicCholesky<usize>` with APTP-specific metadata.

**Fields**:

| Field | Type | Description |
|-------|------|-------------|
| inner | SymbolicCholesky\<usize\> | faer's symbolic factorization result (permutation, L structure, simplicial/supernodal variant) |
| etree | Vec\<isize\> | Elimination tree parent pointers (column j's parent = etree[j]; -1 for root). Computed via `prefactorize_symbolic_cholesky`. |
| col_counts | Vec\<usize\> | Predicted nonzero count per column of L (from prefactorize). Used for buffer estimation. |
| pivot_buffer | Vec\<usize\> | Per-supernode (supernodal mode) or per-column (simplicial mode) extra space estimates for APTP delayed pivots. |

**Invariants**:
- `inner.nrows() == inner.ncols()` (square matrix)
- `etree.len() == inner.nrows()` (one parent per column)
- `col_counts.len() == inner.nrows()` (one count per column)
- `pivot_buffer.len()` equals `n_supernodes()` (supernodal) or `inner.nrows()` (simplicial)
- All `pivot_buffer` entries are non-negative

**Relationships**:
- Produced by `AptpSymbolic::analyze()`
- Consumed by numeric factorization (Phase 5-6) via `inner` for workspace allocation and `pivot_buffer` for APTP-specific sizing
- References faer's `SymbolicCholesky` (composed, not wrapped)

**Lifecycle**:
- Created once per sparsity pattern
- Reusable across multiple numeric factorizations with the same structure but different values
- Immutable after creation

### SymbolicStatistics

Diagnostic summary of the symbolic analysis.

**Fields**:

| Field | Type | Description |
|-------|------|-------------|
| dimension | usize | Matrix dimension (n for n×n) |
| predicted_nnz | usize | Predicted nonzero count in L factor |
| average_col_count | f64 | Mean nonzeros per column in L (predicted_nnz / dimension) |
| is_supernodal | bool | Whether faer chose supernodal over simplicial |
| n_supernodes | Option\<usize\> | Number of supernodes (Some if supernodal, None if simplicial) |
| total_pivot_buffer | usize | Sum of all pivot buffer estimates |

**Invariants**:
- `dimension > 0` (0×0 matrices produce empty statistics)
- `predicted_nnz >= dimension` (at least the diagonal)
- `average_col_count >= 1.0` (at least diagonal entry per column)
- `n_supernodes.is_some() == is_supernodal`

**Relationships**:
- Derived from `AptpSymbolic` (computed on demand, not stored)
- Used for diagnostics, benchmarking, and debugging

## State Transitions

AptpSymbolic has no mutable state — it is created once and read thereafter. The broader solver pipeline has this state flow:

```
[Input Matrix] → analyze() → [AptpSymbolic] → factorize() → [AptpNumeric] → solve() → [Solution]
                  Phase 3          ↑                Phase 5-6                   Phase 7
                                   |
                          reusable for same
                          sparsity pattern
```

## Type Integration with Existing Codebase

| Existing Type | Relationship to AptpSymbolic |
|---------------|------------------------------|
| SparseError | Error type returned by `analyze()` |
| Inertia | Not used at symbolic phase (computed during numeric factorization) |
| MixedDiagonal | Sized using AptpSymbolic's `predicted_nnz` and `pivot_buffer` during numeric phase |
| PivotType | Not used at symbolic phase (determined during numeric factorization) |
| perm_from_forward | Can be used to construct custom permutations passed as `SymmetricOrdering::Custom` |

## faer Types at Boundary

These faer types appear in AptpSymbolic's public API (transparent composition):

| faer Type | Where Used |
|-----------|------------|
| `SymbolicSparseColMatRef<'_, usize>` | Input to `analyze()` — the matrix's sparsity pattern |
| `SymmetricOrdering<'_, usize>` | Input to `analyze()` — ordering strategy |
| `PermRef<'_, usize>` | Returned by `perm()` — the fill-reducing permutation |
| `SymbolicCholeskyRaw<usize>` | Accessible via `raw()` — simplicial/supernodal variant |
| `SymbolicSupernodalCholesky<usize>` | Accessible when supernodal (via `raw()` match) |
| `SymbolicSimplicialCholesky<usize>` | Accessible when simplicial (via `raw()` match) |
| `EliminationTreeRef<'_, usize>` | NOT exposed directly — etree stored as `Vec<isize>` for ownership |
