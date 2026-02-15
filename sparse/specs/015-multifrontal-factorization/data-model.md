# Data Model: Multifrontal Numeric Factorization

**Date**: 2026-02-15
**Feature**: 015-multifrontal-factorization

## Entity Relationship Diagram

```
AptpSymbolic (Phase 3)     SparseColMat (faer)
       │                         │
       └──────────┬───────────────┘
                  │
                  ▼
          SupernodeInfo[]
           (precomputed)
                  │
                  ▼
    ┌─── Factorization Loop ───┐
    │  (postorder traversal)   │
    │                          │
    │  ┌──────────────────┐    │
    │  │  FrontalMatrix    │◄──┤── scatter original entries
    │  │  (dense, per-SN)  │◄──┤── extend-add child ContributionBlocks
    │  └────────┬─────────┘    │
    │           │              │
    │     aptp_factor_in_place │
    │           │              │
    │     ┌─────┴──────┐       │
    │     ▼            ▼       │
    │ FrontFactors  ContributionBlock
    │ (stored)      (passed to parent)
    └──────────────────────────┘
                  │
                  ▼
            AptpNumeric
         (complete result)
```

## Entities

### SupernodeInfo

Unified abstraction over faer's supernodal and simplicial decompositions. Precomputed once before the factorization loop.

| Field | Type | Description |
|-------|------|-------------|
| col_begin | usize | First column index of this supernode (inclusive) |
| col_end | usize | Past-the-end column index (exclusive) |
| pattern | Vec\<usize\> | Off-diagonal row indices (global column indices of non-fully-summed rows) |
| parent | Option\<usize\> | Parent supernode index in assembly tree, or None for root |

**Derived properties**:
- `ncols() = col_end - col_begin` — number of fully-summed columns (supernode size)
- `nrows() = ncols() + pattern.len()` — total frontal matrix dimension
- `is_root() = parent.is_none()`

**Construction**:
- **Supernodal**: From `AptpSymbolic::supernode_begin/end/pattern/parent`
- **Simplicial**: Each column j has `col_begin=j, col_end=j+1, pattern=etree-derived row indices`

### FrontalMatrix

Temporary dense matrix for assembling and factoring one supernode. Allocated, used, and deallocated within a single iteration of the factorization loop.

| Field | Type | Description |
|-------|------|-------------|
| data | Mat\<f64\> | Dense m x m storage (lower triangle used). m = num_fully_summed + num_non_fully_summed |
| row_indices | Vec\<usize\> | Global column indices for each local row (length m). Maps local position to global column index |
| num_fully_summed | usize | First k rows/columns are fully summed (supernode columns + delayed from children) |

**Partitioning** (views, not separate allocations):
- `F11 = data[0..k, 0..k]` — fully-summed block (factored by APTP)
- `F21 = data[k..m, 0..k]` — subdiagonal block (L21 computed by APTP's Schur updates)
- `F22 = data[k..m, k..m]` — contribution block (updated by APTP's Schur updates)

**Assembly**: Populated in two phases:
1. **Scatter**: Original sparse matrix entries at permuted positions
2. **Extend-add**: Child contribution blocks merged using row-index mapping

### FrontFactors

Stored result of factoring one supernode. Persists in memory for Phase 7 (triangular solve).

| Field | Type | Description |
|-------|------|-------------|
| l11 | Mat\<f64\> | Unit lower triangular factor (ne x ne). Extracted from factored frontal matrix |
| d11 | MixedDiagonal | Block diagonal with 1x1 and 2x2 pivots (ne entries) |
| l21 | Mat\<f64\> | Subdiagonal factor block (r x ne). Rows correspond to non-fully-summed pattern |
| local_perm | Vec\<usize\> | APTP local pivot permutation (length k). Maps factored position to original front-local column |
| num_eliminated | usize | Number of columns successfully eliminated (ne <= k) |
| col_indices | Vec\<usize\> | Global column indices for the eliminated columns (length ne). Needed for Phase 7 solve |
| row_indices | Vec\<usize\> | Global row indices for L21 rows (length r). Needed for Phase 7 solve |

**Invariants**:
- `l11.nrows() == l11.ncols() == num_eliminated`
- `d11.dimension() == num_eliminated`
- `l21.ncols() == num_eliminated`
- `l21.nrows() == row_indices.len()`
- `col_indices.len() == num_eliminated`

### ContributionBlock

The Schur complement and delayed columns from a factored supernode. Temporary — consumed by the parent's extend-add step.

| Field | Type | Description |
|-------|------|-------------|
| data | Mat\<f64\> | Dense (m - ne) x (m - ne) trailing submatrix from factored frontal matrix |
| row_indices | Vec\<usize\> | Global column indices for the contribution rows/columns (length m - ne) |
| num_delayed | usize | First num_delayed entries are delayed columns (become fully-summed at parent) |

**Structure within the contribution block**:
- Positions `0..num_delayed`: delayed columns with partially-updated values. Global indices from `AptpFactorResult.delayed_cols`
- Positions `num_delayed..size`: non-fully-summed rows/columns with Schur complement applied. Global indices from original `frontal_row_indices[k..m]`

### AptpNumeric

Complete numeric factorization result. Returned to the caller for use in Phase 7 (triangular solve).

| Field | Type | Description |
|-------|------|-------------|
| front_factors | Vec\<FrontFactors\> | Per-supernode factors, indexed by supernode ID (0..n_supernodes) |
| stats | FactorizationStats | Aggregate factorization statistics |
| n | usize | Matrix dimension |

**Notes**:
- `front_factors[s]` may have `num_eliminated == 0` if all columns of supernode s were delayed
- Phase 7 uses `front_factors` combined with `AptpSymbolic` for the solve path

### FactorizationStats

Aggregate statistics from the factorization.

| Field | Type | Description |
|-------|------|-------------|
| total_1x1_pivots | usize | Total 1x1 pivots across all supernodes |
| total_2x2_pivots | usize | Total 2x2 pivot pairs across all supernodes |
| total_delayed | usize | Total delay events across all supernodes |
| max_front_size | usize | Largest frontal matrix dimension encountered |

## State Transitions

```
FrontalMatrix lifecycle:
  Allocated → Assembled (scatter + extend-add) → Factored (aptp_factor_in_place)
                                                      ↓
                                               Extract FrontFactors + ContributionBlock
                                                      ↓
                                               FrontalMatrix deallocated

ContributionBlock lifecycle:
  Created (from factored FrontalMatrix) → Consumed (extend-add into parent FrontalMatrix) → Deallocated

FrontFactors lifecycle:
  Created (extracted from factored FrontalMatrix) → Stored in AptpNumeric → Used by Phase 7 solve
```

## Global-to-Local Index Mapping

During assembly, each frontal matrix maintains a mapping from global column indices to local positions:

```
global_to_local: HashMap<usize, usize> or dense array of size n

For supernode s with frontal rows [c0, c1, ..., c_{m-1}]:
  global_to_local[c_i] = i

Used for:
  1. Scatter: place A[perm[r], perm[c]] at frontal[local_r, local_c]
  2. Extend-add: map child contribution indices to parent local indices
```

For efficiency, a dense array `global_to_local[0..n]` initialized to a sentinel value (e.g., usize::MAX) is reused across supernodes, with cleanup after each iteration.
