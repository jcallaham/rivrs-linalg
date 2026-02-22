# Data Model: Assembly & Extraction Optimization (Phase 9.1c)

**Feature**: 022-assembly-extraction-opt
**Date**: 2026-02-22

## New Entities

### AssemblyMaps

Pre-computed index mappings for assembly operations, stored alongside the symbolic analysis result. One instance per `AptpSymbolic`, computed once during `analyze()` and reused across factorizations.

**Fields:**
- `amap_pairs: Vec<u32>` — Flattened source/destination pairs for all supernodes. Each pair is `[src_csc_index, dest_frontal_index]` where `src_csc_index` is the position in the sparse matrix's CSC values array and `dest_frontal_index` is the column-major linear index in the frontal matrix. Length: 2 * total_original_entries_across_all_supernodes.
- `amap_offsets: Vec<usize>` — Per-supernode start offset into `amap_pairs`. `amap_offsets[s]..amap_offsets[s+1]` gives the range for supernode s. Length: num_supernodes + 1.
- `ea_map: Vec<u32>` — Flattened child contribution row → parent local row mappings for all parent-child edges. For each child's off-diagonal pattern row, stores the corresponding local row index in the parent's frontal matrix (assuming zero delays). Length: sum of all children's off-diagonal pattern sizes.
- `ea_offsets: Vec<usize>` — Per-child start offset into `ea_map`. Indexed by a flattened child enumeration (all children of all supernodes in tree order). Length: total_children + 1.
- `ea_child_snode: Vec<usize>` — Maps each entry in `ea_offsets` to its child supernode index. Used to associate a child offset with the correct supernode when iterating over a parent's children. Length: total_children.

**Lifecycle:**
1. Computed during `AptpSymbolic::analyze()` (or a new `AptpSymbolic::compute_assembly_maps()` called from analyze).
2. Stored as a field on `AptpSymbolic`.
3. Read-only during factorization — multiple factorizations can read concurrently.
4. Dropped when `AptpSymbolic` is dropped.

**Invariants:**
- `amap_offsets` is monotonically non-decreasing with length `num_supernodes + 1`.
- `amap_pairs.len()` == `2 * amap_offsets[num_supernodes]`.
- All `dest_frontal_index` values in a supernode's amap range are valid for a frontal matrix of size `front_size[s]`.
- `ea_map` entries are valid local row indices in the parent supernode's frontal matrix (assuming zero delays).

### FallbackFlag

Not a separate entity but a per-supernode boolean: at factorization time, if a child has delays (ne < k), the extend-add for that child falls back to the current `global_to_local` approach instead of using the precomputed map.

## Modified Entities

### AptpSymbolic (existing)

**Current:** Wraps faer's `SymbolicCholesky` + APTP-specific pivot buffer estimates. Stores supernode info, elimination tree, permutations.

**After change:** Gains an `AssemblyMaps` field (or the individual arrays directly) computed during analysis. The analysis phase is extended to iterate over the sparse matrix structure and build scatter maps.

**Key constraint:** The scatter maps must be computed AFTER amalgamation (which modifies supernode structure) but BEFORE factorization. This fits naturally into the existing `analyze()` → `factor()` pipeline since amalgamation happens inside `analyze()`.

### scatter_original_entries_multi (existing function in numeric.rs)

**Current signature:**
```
fn scatter_original_entries_multi(
    frontal, matrix, perm_fwd, perm_inv, global_to_local, owned_ranges, scaling
)
```

**After change:** Gains an `amap` parameter (slice of the precomputed pairs for this supernode). When amap is available and the supernode has no delayed children that would shift positions, the function uses the amap for direct indexed scatter. Otherwise falls back to current logic.

### extend_add (existing function in numeric.rs)

**Current signature:**
```
fn extend_add(parent, child, global_to_local)
```

**After change:** Gains an optional `ea_row_map: Option<&[u32]>` parameter. When provided (zero-delay case), uses the precomputed row mapping instead of `global_to_local` lookups. When `None` (child has delays), falls back to current logic.

### extract_contribution (existing function in numeric.rs)

**Current:** Element-by-element column-major lower-triangle copy.

**After change:** Per-column `copy_from_slice` using `col_as_slice()`. No new parameters — the optimization is internal to the function body.

### extract_front_factors (existing function in numeric.rs)

**Current:** Element-by-element L11 and L21 copy.

**After change:** Per-column `copy_from_slice` for L21. Per-column slice copy for L11 (where pivot type allows contiguous segments). No new parameters.

### zero_frontal (existing method on FactorizationWorkspace)

**Current:** Nested loop with element-by-element indexed zeroing.

**After change:** Per-column `fill(0.0)` on column slices. No new parameters.

## Unchanged Entities

- **FrontalMatrix**: No changes to struct. Extraction functions are modified internally.
- **ContributionBlock**: No changes. Still owns its data via `Mat<f64>`.
- **FrontFactors**: No changes to stored data.
- **AptpKernelWorkspace**: No changes.
- **FactorizationWorkspace**: Only `zero_frontal` implementation changes.
- **SparseLDLT**: No changes to public API.
- **MixedDiagonal**: No changes.
- **AptpNumeric**: No changes to stored data. Assembly loop uses new maps.

## Entity Relationships

```
AptpSymbolic ---contains---> AssemblyMaps (amap + ea_map)
    |                              |
    | provides                     | provides
    v                              v
factor_single_supernode ----uses---> amap (per-supernode scatter)
    |                              |
    | calls                        | provides
    v                              v
scatter_original_entries_multi     extend_add (per-child row mapping)
    |                              |
    | scatter from                 | merge from
    v                              v
SparseColMat (original A)         ContributionBlock (child Schur complement)
    |                              |
    | into                         | into
    v                              v
FrontalMatrix (dense workspace)
```

## Memory Budget

For a matrix with `nnz(A)` nonzeros and `S` supernodes:

| Buffer | Size (bytes) | c-71 (nnz=860K) | c-big (nnz=2.3M) |
|--------|-------------|------------------|-------------------|
| amap_pairs | 8 * nnz(A) | 6.9 MB | 18.7 MB |
| amap_offsets | 8 * (S+1) | 50 KB | 252 KB |
| ea_map | 4 * sum(child_pattern_sizes) | ~2-5 MB | ~5-15 MB |
| ea_offsets | 8 * (total_children+1) | ~50 KB | ~250 KB |
| **Total** | | **~10-12 MB** | **~25-35 MB** |

This is a one-time cost paid during analysis, amortized across factorizations.
