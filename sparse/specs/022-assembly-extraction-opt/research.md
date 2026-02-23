# Research: Assembly & Extraction Optimization (Phase 9.1c)

**Feature**: 022-assembly-extraction-opt
**Date**: 2026-02-22

## R1: What can be precomputed at symbolic time vs factorization time?

**Decision**: Scatter maps (amap) and extend-add row mappings are precomputable for the zero-delay case. Delayed column handling falls back to current per-entry assembly.

**Rationale**: The symbolic analysis provides all the information needed to compute frontal matrix row patterns: supernode column ranges, off-diagonal row patterns (from `supernode_pattern(s)`), and the elimination tree parent-child relationships. The only dynamic element is delayed columns, which shift row positions at factorization time. Since delays are rare (0 on c-71/c-big, 24 on dixmaanl/60K), precomputing for the zero-delay case and falling back on delay covers 99%+ of supernodes.

**Alternatives considered**:
- **Full parameterized maps** (stores offset adjustments for delays): Adds complexity for a rare case. Rejected.
- **Minimal precomputation** (only CSC column boundaries): Saves no per-entry work, only loop setup overhead. Rejected.
- **Lazy map construction** (compute on first factorization, cache): Loses the analyze-once-factor-many benefit. Rejected.

**Key data available at symbolic time**:
- `SupernodeInfo.col_begin/col_end` — which columns each supernode owns
- `SupernodeInfo.owned_ranges` — non-contiguous column ranges (amalgamation)
- `SupernodeInfo.pattern` — off-diagonal row indices (global permuted)
- `SupernodeInfo.parent` — parent supernode in elimination tree
- Forward/inverse permutation arrays
- Original sparse matrix CSC structure (col_ptrs, row_indices)

**Key data NOT available at symbolic time**:
- Number of delayed columns per supernode (determined by APTP kernel)
- Actual contribution block sizes (depend on how many columns are eliminated)
- MC64 scaling values (provided per factorization)

## R2: SPRAL's amap structure and assembly approach

**Decision**: Follow SPRAL's amap pattern — store (source_csc_index, destination_frontal_index) pairs per supernode.

**Rationale**: SPRAL's `assemble.hxx:51-80` (`add_a_block`) demonstrates that a flat array of source→destination index pairs eliminates all per-entry permutation, global-to-local mapping, and deduplication logic. The assembly inner loop becomes a single indexed write: `lcol[dest] = aval[src]`.

**Key SPRAL findings**:
- `amap` is precomputed during Fortran preprocessing (not C++ symbolic phase)
- Format: pairs of 1-based linear indices `[src_in_A, dest_in_frontal]`
- Destination encodes column-major position in the dense frontal storage
- Delayed columns handled via `ndelay_in` offset at assembly time
- Scaling applied separately (not baked into amap)
- SPRAL does NOT explicitly zero frontal matrices — relies on calloc-style allocator

**SPRAL's two-phase assembly**:
- `assemble_pre()`: assembles A entries + child contributions into frontal matrix BEFORE factorization
- `assemble_post()`: assembles child contributions directly into parent's CONTRIBUTION block AFTER factorization (bypasses frontal matrix for Schur complement)
- This two-phase approach eliminates the contribution block copy entirely for the parent's perspective
- Out of scope for 9.1c but noted as a potential follow-up (significant architectural change)

## R3: faer bulk copy API feasibility

**Decision**: Use `col_as_slice()` / `col_as_slice_mut()` for per-column bulk operations. Do NOT use submatrix views for bulk copies (they lose column contiguity).

**Rationale**: faer's `Mat<f64>` uses column-major storage with contiguous columns. The `col_as_slice(j)` method returns a `&[f64]` slice for column j. However, `submatrix_mut(row_start, col_start, nrows, ncols)` at `row_start > 0` creates a view where columns are NOT contiguous (each column starts at an offset within the parent column, with gaps).

**Key findings**:
- `Mat::col_as_slice(j)` → `&[f64]` of length `nrows` (contiguous, column j)
- `Mat::col_as_slice_mut(j)` → `&mut [f64]` of length `nrows` (contiguous, column j)
- For a source range within a column: `col_as_slice(j)[start..end]` gives contiguous `&[f64]`
- `copy_from_slice` is essentially `memcpy` — fastest possible bulk copy
- `slice.fill(0.0)` is essentially `memset` — fastest possible bulk zeroing
- Submatrix views (`submatrix_mut`) do NOT guarantee contiguous columns — must access parent Mat directly

**Extraction applicability**:
- `extract_contribution`: YES — each column's lower triangle is a contiguous slice
- `extract_front_factors` L21: YES — each column from row `ne` to `m` is contiguous
- `extract_front_factors` L11: PARTIAL — 1x1 pivots can use bulk copy; 2x2 pivots need care with diagonal entries
- `zero_frontal`: YES — each column from row `j` to `m` is a contiguous slice

## R4: Memory cost analysis for scatter maps

**Decision**: Accept O(nnz) memory for amap storage using u32 indices. Total cost is ~16 bytes per original matrix entry.

**Rationale**: The amap stores two u32 indices per entry: source position in CSC values array and destination position in the frontal matrix's column-major storage. Using u32 instead of usize halves memory since frontal matrices are bounded by max_front_size (~11K), well within u32 range. The CSC values array length is bounded by nnz(A), also within u32 for practical matrices.

**Memory estimates**:

| Matrix | nnz(A) | amap memory | % of factor storage |
|--------|--------|-------------|---------------------|
| c-71 | 859,520 | ~14 MB | ~2-3% |
| c-big | 2,340,859 | ~37 MB | ~2-3% |
| nd12k | 14,220,946 | ~227 MB | ~5-10% |
| G3_circuit | 7,660,826 | ~122 MB | ~3-5% |

For the extend-add maps, memory is proportional to the sum of contribution block sizes across all parent-child edges, which is bounded by predicted_nnz(L). Estimated at ~4-8 bytes per entry.

**Decision on nd12k**: The 227 MB amap for nd12k is significant but acceptable — nd12k's factor storage itself is several GB. The amap is allocated once during analysis and reused across factorizations.
