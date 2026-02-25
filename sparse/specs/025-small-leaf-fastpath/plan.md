# Implementation Plan: Small Leaf Subtree Fast Path

**Branch**: `025-small-leaf-fastpath` | **Date**: 2026-02-24 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/025-small-leaf-fastpath/spec.md`

## Summary

Simplicial matrices (dixmaanl, mario001, bloweybq, etc.) are 1.5-3.3x slower than SPRAL because every supernode — even tiny ones with front sizes under 100 — pays the full cost of the general multifrontal machinery (frontal matrix allocation, zeroing, g2l mapping, contribution extraction). This feature classifies leaf subtrees where all supernodes have front_size < 256 and factors them via a streamlined code path that processes the entire subtree with a single small workspace, avoiding per-supernode dispatch overhead. The fast path produces identical `FrontFactors` and `ContributionBlock` output, so the solve path and parent extend-add assembly require no modification.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (dense LA, CSC), rayon 1.x (parallelism), serde/serde_json (diagnostic export)
**Storage**: N/A (in-memory numerical computation)
**Testing**: cargo test (unit + SuiteSparse integration), criterion (benchmarks), baseline_collection (SPRAL comparison)
**Target Platform**: Linux x86_64 (primary), cross-platform Rust
**Project Type**: Single Rust crate (library)
**Performance Goals**: Simplicial matrices within 1.5x SPRAL; no regression on non-simplicial matrices
**Constraints**: Backward error < 5e-11 for all 65 SuiteSparse matrices; median SPRAL ratio ≤ 1.0x
**Scale/Scope**: ~300 lines new code in numeric.rs, ~20 lines in solver.rs, ~200 lines tests

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Fast path produces same FrontFactors/ContributionBlock; validated by reconstruction tests + backward error on full SuiteSparse suite |
| II. Clean Room | PASS | SPRAL (BSD-3) consulted for algorithm architecture; academic papers cited. No HSL or GPL code. |
| III. TDD | PASS | Classification tests first (known tree structures), then fast-path correctness (hand-constructed + SuiteSparse), then performance |
| IV. Documentation | PASS | Algorithm references documented in spec. Implementation will cite SPRAL SmallLeafNumericSubtree.hxx (BSD-3) |
| V. Numerical Stability | PASS | Same APTP kernel (aptp_factor_in_place) with same threshold/fallback. Delayed pivots propagate identically. |
| VI. Structured Development | PASS | Phase 9.1f in ssids-plan.md. Follows phased plan. |
| VII. Code Quality | PASS | New types (SmallLeafSubtree) are minimal. No public API changes except FactorOptions.small_leaf_threshold. |

**Post-design re-check**: All gates still pass. The fast path is an internal optimization that doesn't change the solver's external behavior or API contract (except the new threshold option).

## Project Structure

### Documentation (this feature)

```text
specs/025-small-leaf-fastpath/
├── plan.md              # This file
├── spec.md              # Feature specification
├── research.md          # Phase 0 research findings
├── data-model.md        # Type changes and new types
├── quickstart.md        # Build/test/benchmark guide
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code (modifications)

```text
src/aptp/
├── numeric.rs           # PRIMARY: classification, fast-path function, level-set integration
│   ├── SupernodeInfo    # Add in_small_leaf: bool field
│   ├── SmallLeafSubtree # New struct (root, nodes, max_front_size, parent_of_root)
│   ├── classify_small_leaf_subtrees()     # New function
│   ├── factor_small_leaf_subtree()        # New function
│   └── factor_tree_levelset()             # Modified: pre-pass before level-set loop
└── solver.rs            # FactorOptions/SolverOptions: add small_leaf_threshold field
```

**Structure Decision**: All changes are within the existing `src/aptp/` module. No new files — the fast path is an internal optimization within the factorization loop. Classification and fast-path factorization are private functions in `numeric.rs`.

## Design

### Architecture Overview

```
AptpNumeric::factor()
  ├── build_supernode_info()           [existing]
  ├── amalgamate_supernodes()          [existing]
  ├── build_children_map()             [existing]
  ├── classify_small_leaf_subtrees()   [NEW — marks in_small_leaf, returns Vec<SmallLeafSubtree>]
  ├── build_assembly_maps()            [existing]
  └── factor_tree_levelset()           [MODIFIED]
        ├── Pre-pass: for each SmallLeafSubtree:
        │     └── factor_small_leaf_subtree()  [NEW]
        │           ├── Allocate SmallLeafWorkspace (threshold × threshold buffer)
        │           ├── For each node in postorder:
        │           │     ├── Zero relevant portion of buffer
        │           │     ├── Scatter original entries (using assembly maps or direct CSC)
        │           │     ├── Extend-add child contributions (within-subtree)
        │           │     ├── aptp_factor_in_place() on buffer submatrix
        │           │     ├── Extract FrontFactors (same as general path)
        │           │     └── Compute contribution GEMM → store for parent
        │           └── Return Vec<SupernodeResult> + root ContributionBlock
        ├── Adjust remaining_children counts (decrement for subtree root parents)
        ├── Build initial ready set (excluding in_small_leaf nodes)
        └── Main level-set loop                [existing, unchanged]
```

### Key Design Decisions

**D1: Pre-pass vs inline dispatch** (see research.md R3)
Pre-pass before level-set loop. Cleaner separation, no conditional in hot loop.

**D2: Workspace strategy** (see research.md R4)
Dedicated lightweight workspace per subtree. Fixed max_front = threshold (≤256). Buffer fits in L2 cache. Single allocation reused across all nodes in subtree.

**D3: Assembly within subtree**
Within a subtree, child contributions are consumed immediately by the parent node (no global contributions vector). Only the subtree root's contribution enters the global vector for the main level-set loop's extend-add.

**D4: Assembly maps reuse**
The existing `AssemblyMaps` are computed for all supernodes (including small-leaf ones) before the pre-pass. The fast path reuses the precomputed scatter maps (`amap_entries`) and extend-add maps (`ea_map`) where applicable (no-delay fast path). This avoids duplicating map computation logic.

**D5: Diagnostic support**
`PerSupernodeStats` are collected for each node in the fast path, identical to the general path. The `diagnostic` feature's timing instrumentation is applied within `factor_small_leaf_subtree()` using the same `#[cfg(feature = "diagnostic")]` guards.

### Interaction with Parallel Path

When `options.par != Par::Seq`, the level-set loop uses rayon for tree-level parallelism. The pre-pass factors small-leaf subtrees **before** entering the parallel section:
- Independent small-leaf subtrees could be parallelized via `par_iter` (each has its own workspace)
- But the target workload (simplicial matrices) is sequential anyway — the benefit is in cache locality and overhead reduction, not parallelism
- The main level-set loop's parallel path is unchanged; it simply sees fewer ready nodes initially

### Classification Algorithm Detail

```
classify_small_leaf_subtrees(supernodes, children_map, threshold) → Vec<SmallLeafSubtree>:

1. Compute front_size[s] = sum(owned_ranges[s].len()) + pattern[s].len() for all s
2. Bottom-up pass (postorder — s=0..n_supernodes-1):
     if front_size[s] >= threshold:
         in_small_leaf[s] = false
     else if children_map[s].is_empty():  // leaf
         in_small_leaf[s] = true
     else if children_map[s].iter().all(|&c| in_small_leaf[c]):
         in_small_leaf[s] = true
     else:
         in_small_leaf[s] = false
3. Identify subtree roots: in_small_leaf[s] = true AND (parent is None OR !in_small_leaf[parent])
4. For each root, collect descendant nodes via DFS into postorder Vec
5. Filter: discard subtrees with < 2 nodes
6. Compute max_front_size per subtree
7. Return Vec<SmallLeafSubtree>
```

### Fast Path Factorization Detail

```
factor_small_leaf_subtree(subtree, supernodes, matrix, perm, options, scaling, assembly_maps):

1. Allocate SmallLeafWorkspace:
     frontal_data = Mat::zeros(subtree.max_front_size, subtree.max_front_size)
     global_to_local = vec![NOT_IN_FRONT; n]
     contrib_buffer = Mat::new()  // lazy, recycled

2. Local contributions storage: Vec<Option<ContributionBlock>> indexed by subtree-local ID

3. For each node_id in subtree.nodes (postorder):
     sn = &supernodes[node_id]

     a. Collect delayed columns from children (within-subtree children only)
     b. Compute front structure: frontal_row_indices = [owned_cols, delayed, pattern]
     c. Zero frontal_data[0..m, 0..m]
     d. Build g2l mapping for m entries
     e. Scatter original entries using assembly_maps.amap_entries[node_id]
     f. Extend-add child contributions from local_contributions
     g. Call aptp_factor_in_place(frontal_data.submatrix_mut(0,0,m,m), k, options)
     h. Compute contribution GEMM (if not root or not fully eliminated)
     i. Extract FrontFactors (same extract_front_factors as general path)
     j. Extract ContributionBlock → store in local_contributions for parent
     k. Reset g2l entries
     l. Collect PerSupernodeStats

4. Return: Vec<(FrontFactors, Option<ContributionBlock>, PerSupernodeStats)>
   The last entry's ContributionBlock is the subtree root's contribution for the parent.
```

Note: Steps (c)-(k) are structurally identical to `factor_single_supernode()`. The optimization comes from:
- Single small workspace (cache-resident) vs per-node workspace from pool
- No workspace pool Mutex acquisition
- No workspace capacity recheck (max_front bounded)
- Sequential processing avoids rayon dispatch overhead
- Buffer stays warm across nodes (temporal locality)

## Algorithm References

See [spec.md](spec.md) Algorithm References section for full citations. Key references for implementation:

- **SPRAL `SymbolicSubtree.hxx:57-84`**: Classification algorithm (bottom-up, leaf-first)
- **SPRAL `SmallLeafNumericSubtree.hxx:187-446`**: Indefinite factorization path (assemble_pre/factor_node/assemble_post)
- **SPRAL `SymbolicSubtree.hxx:69-84`**: Subtree identification (walk from leaf to root)
- **Duff, Hogg & Lopez (2020)**: Algorithm 6.1 (find_subtree_partition)
- **Davis (2016) Section 11**: Multifrontal tree/node parallelism trade-offs

## Complexity Tracking

No constitution violations. No complexity tracking needed.
