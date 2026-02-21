# Data Model: METIS Nested Dissection Ordering

**Feature**: 011-metis-ordering
**Date**: 2026-02-12

## Entities

### 1. Adjacency Structure (internal, not exposed)

The graph representation extracted from a symmetric sparse matrix for METIS consumption.

**Attributes**:
- `n`: number of vertices (= matrix dimension)
- `xadj`: row pointer array of length `n + 1`, where `xadj[i]..xadj[i+1]` indexes into `adjncy` for vertex `i`'s neighbors
- `adjncy`: flattened neighbor list of length `xadj[n]`, containing column indices (excluding self-loops)

**Invariants**:
- `xadj[0] == 0`
- `xadj` is monotonically non-decreasing
- All values in `adjncy` are in `[0, n)`
- No self-loops: vertex `i` never appears in `adjncy[xadj[i]..xadj[i+1]]`
- Symmetry: if `j` appears in vertex `i`'s neighbor list, then `i` appears in vertex `j`'s neighbor list

**Source**: Extracted from `SymbolicSparseColMatRef` (faer's structural CSC representation)

**Lifecycle**: Temporary — allocated, passed to METIS, then deallocated. Not stored.

### 2. METIS Permutation (output, exposed as `Perm<usize>`)

The fill-reducing ordering produced by METIS, represented as faer's `Perm<usize>`.

**Attributes**:
- Forward array (`fwd`): length `n`, where `fwd[new_pos] = old_idx`
- Inverse array (`inv`): length `n`, where `inv[old_idx] = new_pos`

**Invariants**:
- Both arrays are valid permutations of `{0, 1, ..., n-1}`
- `fwd[inv[i]] == i` for all `i` in `0..n`
- `inv[fwd[j]] == j` for all `j` in `0..n`

**Source**: METIS `METIS_NodeND` output arrays, converted from `i32` to `usize`

**Lifecycle**: Returned to caller. Typically passed to `AptpSymbolic::analyze()` via `SymmetricOrdering::Custom(perm.as_ref())`.

## Relationships

```
SparseColMat/SymbolicSparseColMatRef
         │
         │ extract adjacency (exclude diagonal, ensure symmetry)
         ▼
   Adjacency Structure (xadj, adjncy)
         │
         │ METIS_NodeND
         ▼
   Perm<usize> (forward + inverse arrays)
         │
         │ SymmetricOrdering::Custom(perm.as_ref())
         ▼
   AptpSymbolic::analyze()
```

## Type Conversions at FFI Boundary

| Direction | From | To | Validation |
|-----------|------|----|------------|
| Rust → METIS | `usize` indices | `i32` (`idx_t`) | Dimension fits in `i32` |
| METIS → Rust | `i32` (`idx_t`) arrays | `usize` arrays | Non-negative values |
| METIS → faer | `(perm_i32, iperm_i32)` | `Perm<usize>` | Valid permutation |
