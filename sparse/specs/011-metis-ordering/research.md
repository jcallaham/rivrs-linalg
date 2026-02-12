# Research: METIS Nested Dissection Ordering

**Feature**: 011-metis-ordering
**Date**: 2026-02-12

## Decision 1: Which Rust crate for METIS bindings?

**Decision**: Use `metis-sys` (v0.3.x) directly for raw FFI bindings to `METIS_NodeND`.

**Rationale**: The high-level `metis` crate (v0.2.x) wraps graph *partitioning* functions (`METIS_PartGraphKway`, `METIS_PartGraphRecursive`) but does NOT wrap the node *ordering* function (`METIS_NodeND`), which is what we need for fill-reducing permutations. `metis-sys` exposes `METIS_NodeND` as an unsafe FFI function. We write our own safe Rust wrapper around it.

**Alternatives considered**:
- `metis` crate (high-level): Does not wrap `METIS_NodeND`. Would pull in unnecessary partitioning wrappers.
- Pure Rust reimplementation: METIS is ~15K lines of C with decades of tuning. Reimplementation would be impractical and unlikely to match quality. Does not violate "prefer pure Rust" constitution guideline — that's a preference, not a hard rule, and there is no pure Rust alternative.

## Decision 2: Vendored vs system METIS library

**Decision**: Use vendored METIS source (default behavior of `metis-sys`).

**Rationale**: `metis-sys` vendors METIS 5.x C source and compiles it via `build.rs` using the `cc` crate. This requires only a C compiler (already present in all our build environments). No `apt install`, `brew install`, or system library management needed. Builds are fully reproducible.

**Alternatives considered**:
- System-installed libmetis (`use-system` feature): Requires managing system packages across dev/CI/Docker environments. Adds friction for contributors.
- Both (feature-flagged): Unnecessary complexity. Vendoring works everywhere.

## Decision 3: METIS dependency strategy

**Decision**: Required (non-optional) dependency. No feature flag.

**Rationale**: Clarified with user during `/speckit.clarify`. Since vendoring eliminates installation friction and METIS is essential for the solver's core functionality (AMD is inadequate for real-world matrices), there's no benefit to making it optional. Every build compiles METIS from source (~30s additional compile time).

## Decision 4: idx_t type and index conversion

**Decision**: `metis_sys::idx_t` is `i32`. Convert to/from `usize` at the FFI boundary.

**Rationale**: METIS uses 32-bit signed integers for all indices. Our largest test matrix is ~150K dimension, well within `i32` range (~2.1 billion). Conversion is straightforward: `usize as i32` when calling METIS, `i32 as usize` when reading results. Add validation that matrix dimension fits in `i32` before calling METIS.

**Risk**: Matrices with dimension > 2^31 - 1 (~2.1 billion) cannot use METIS. This is not a practical concern — such matrices would require terabytes of memory.

## Decision 5: Permutation semantics mapping (METIS → faer)

**Decision**: Map METIS output arrays to faer's `Perm<usize>` convention.

**Details**:
- METIS `METIS_NodeND` returns two arrays:
  - `perm[i]` = old vertex at new position `i` (new → old mapping). METIS manual: `A' = A(perm, perm)`.
  - `iperm[j]` = new position for old vertex `j` (old → new mapping). Inverse of perm.
- faer's `Perm<usize>` convention:
  - Forward array: `fwd[new_pos]` = old index (new → old)
  - Inverse array: `inv[old_idx]` = new position (old → new)
- **Mapping**: METIS `perm` = faer forward; METIS `iperm` = faer inverse
- **Implementation**: `Perm::new_checked(perm_box, iperm_box, n)` — no recomputation needed, both arrays provided by METIS.

This avoids using `perm_from_forward()` (which recomputes the inverse), since METIS already provides both arrays.

**Note**: The initial design (pre-implementation) had the mapping inverted. The METIS manual's MATLAB-style `A' = A(perm, perm)` convention clarifies that `perm[new] = old`, making it the forward permutation in faer's convention.

## Decision 6: CSC → adjacency structure extraction

**Decision**: Extract the undirected graph adjacency (CSR format) from the symmetric sparse matrix's CSC structure.

**Details**:
- METIS expects CSR format: `xadj` (row pointer array, length n+1) and `adjncy` (column indices of neighbors)
- The input `SparseColMat` stores the matrix in CSC format (col pointers + row indices)
- For a symmetric matrix, CSC of the full symmetrized structure gives CSR by transposition
- **Key requirements**:
  - Exclude diagonal entries (self-loops) — METIS operates on the graph, not the matrix
  - Ensure symmetry: if entry (i,j) exists, entry (j,i) must also exist
  - Use only structural pattern (no numerical values)
- The extraction function takes `SymbolicSparseColMatRef` (structural pattern only)

## Decision 7: METIS options configuration

**Decision**: Use METIS defaults with 0-based numbering. No user-configurable options in this phase.

**Details**:
- `METIS_OPTION_NUMBERING = 0` (C-style 0-based indexing) — required
- All other options: -1 (METIS defaults)
- Default `NSEPS = 1` (number of separators per level) is adequate for our use case
- METIS parameter tuning is explicitly out of scope per the spec

## Decision 8: Error handling for METIS failures

**Decision**: Map METIS return codes to `SparseError` variants.

**Details**:
- `METIS_OK` (1): success
- `METIS_ERROR_INPUT` (-2): invalid input → `SparseError` with descriptive message
- `METIS_ERROR_MEMORY` (-3): allocation failure → `SparseError` with descriptive message
- `METIS_ERROR` (-4): general error → `SparseError` with descriptive message
- Edge cases (dimension 0, 1, diagonal matrices): handle before calling METIS by returning trivial identity permutation

## Key API Reference: METIS_NodeND

```c
int METIS_NodeND(
    idx_t *nvtxs,     // number of vertices
    idx_t *xadj,      // CSR row pointers (length nvtxs+1)
    idx_t *adjncy,    // CSR column indices (length xadj[nvtxs])
    idx_t *vwgt,      // vertex weights (NULL for unweighted)
    idx_t *options,    // options array (length METIS_NOPTIONS=40)
    idx_t *perm,      // OUTPUT: new position for each vertex
    idx_t *iperm      // OUTPUT: old vertex at each new position
);
```

## License Compatibility

| Component | License | Compatible with Apache-2.0? |
|-----------|---------|---------------------------|
| METIS 5.x source | Apache-2.0 | Yes |
| `metis-sys` crate | MIT / Apache-2.0 | Yes |
| `cc` build dep | MIT / Apache-2.0 | Yes |
