# Research: Supernode Amalgamation

## Decision 1: Integration Point — Where to Apply Amalgamation

**Decision**: Post-process `Vec<SupernodeInfo>` after `build_supernode_info()` in `numeric.rs`, rather than modifying `AptpSymbolic` or faer.

**Rationale**:
- `build_supernode_info()` (numeric.rs:976-1029) is the single place where faer's `SymbolicCholeskyRaw::Supernodal` is consumed into our `SupernodeInfo` descriptors.
- The numeric factorization loop (`factor_tree_levelset`) only uses `Vec<SupernodeInfo>` and the children map — it never accesses faer's supernodal structure directly after build.
- The solve path operates on `FrontFactors` (col_indices, row_indices, L11, D11, L21) — completely opaque to supernode structure.
- No changes to faer required (no upstream PR, no API mismatch risk).

**Alternatives considered**:
- *Modify AptpSymbolic*: Would require owning supernodal metadata instead of delegating to faer. Higher coupling, more code churn, but would make amalgamation visible to symbolic analysis consumers. Rejected: the only consumer of supernode structure is numeric factorization.
- *Modify faer upstream*: Would require exposing relaxation parameters or adding a post-hoc merge API. Rejected: months of upstream review, and faer's consecutivity limitation is architectural, not a simple knob.
- *Build a separate `AmalgamatedSymbolic` type*: Rejected as over-engineering — `Vec<SupernodeInfo>` is already the right abstraction.

## Decision 2: Column Contiguity — No Renumbering Needed

**Decision**: Do not renumber columns after amalgamation. faer's fundamental supernodes already have contiguous column ranges, and merging parent-child pairs produces contiguous ranges when the child's `col_end == parent's col_begin`.

**Rationale**:
- SPRAL's `sperm` (supernode permutation) exists because SPRAL detects supernodes from scratch on the column-level elimination tree. rivrs starts from faer's contiguous fundamental supernodes.
- When merging child into parent where child is a direct predecessor in column order (`child.col_end == parent.col_begin`), the merged range `[child.col_begin, parent.col_end)` is contiguous.
- When merging supernodes that are NOT adjacent in column order (the nemin condition), the merged supernode's `col_begin..col_end` will span a range that includes columns belonging to OTHER supernodes. This is fine: `scatter_original_entries` loops only over the supernode's own columns, and the frontal matrix assembly uses `global_to_local` which handles arbitrary global indices.
- **Key insight**: The `frontal_rows` construction (numeric.rs:715-718) extends `col_begin..col_end` — for a non-contiguous merge this will include columns from other supernodes in the frontal matrix's fully-summed block. These extra columns have zero original entries but receive contributions from child extend-add, which is correct. The APTP kernel handles them like any other fully-summed column.

**Alternatives considered**:
- *Column renumbering pass*: Would produce truly contiguous merged supernodes but requires re-mapping ALL permutation arrays, row patterns, and the original matrix's column indices. Massive complexity for no correctness benefit. Rejected.

## Decision 3: Merge Predicate — SPRAL's do_merge()

**Decision**: Implement SPRAL's exact two-condition merge predicate from `core_analyse.f90:806-822`.

**Rationale**:
- Condition (a): `cc(par) == cc(node)-1 AND nelim(par) == 1` — this catches fundamental supernodes that faer's consecutivity check missed. These merges introduce zero fill-in.
- Condition (b): `nelim(par) < nemin AND nelim(node) < nemin` — this is the powerhouse for narrow-supernode matrices. Both nodes are small, so the fill-in cost of merging is bounded by `O(nemin^2)` per merge.
- SPRAL's `nemin_default = 32` is battle-tested across their full user base.
- No fill-in limit needed for nemin merges (SPRAL doesn't use one). The `ezero` tracking is for statistics only.

**Alternatives considered**:
- *Fill-ratio-limited merging*: Only merge when fill-in ratio is below threshold. Rejected: adds complexity and SPRAL's empirical evidence shows unconditional nemin merging works well.
- *faer-style relaxed merging with looser parameters*: Already experimentally shown to have no effect (ssids-plan.md Phase 9.1a analysis: c-71 went from 35,372 → 35,308 supernodes with maximally aggressive relaxation).

## Decision 4: Row Pattern Computation

**Decision**: Merged supernode's pattern = sorted union of constituent patterns, excluding columns that became fully-summed (i.e., columns in the merged supernode's own range).

**Rationale**:
- When child merges into parent, the parent's pattern gains any rows that were in the child's pattern but not in the parent's pattern.
- Rows corresponding to the child's own columns are now fully-summed (part of the merged supernode) and must be removed from the pattern.
- The fill entries (rows in the merged pattern that were not in either original pattern) are implicitly zero — the frontal matrix is initialized with `Mat::zeros` and only populated via scatter + extend-add.

## Decision 5: Tree Traversal After Amalgamation

**Decision**: Rebuild the children map and level-set scheduling from the amalgamated `Vec<SupernodeInfo>`. The existing `factor_tree_levelset` works unchanged.

**Rationale**:
- `build_children_map()` (numeric.rs) constructs children from `SupernodeInfo.parent` pointers.
- After amalgamation, parent pointers are updated: if child C had parent P, and P merged into grandparent G, then C's parent becomes G.
- The level-set scheduler (`factor_tree_levelset`) uses `remaining_children[s]` counts and the children map — both are rebuilt from the amalgamated structure.
- Postorder is maintained because merging parent-child pairs preserves the invariant that children come before parents.

## Decision 6: Where to Expose nemin Configuration

**Decision**: Add `nemin: usize` field to `AnalyzeOptions` (default 32) and propagate through `SolverOptions`. The amalgamation pass reads this value. Setting `nemin = 1` effectively disables amalgamation.

**Rationale**:
- Amalgamation is conceptually part of symbolic analysis (it transforms the supernode structure).
- `AnalyzeOptions` already holds `ordering: OrderingStrategy` — adding `nemin` is a natural extension.
- SPRAL exposes `options%nemin` at the same level (datatypes.f90:213).
- `nemin = 1` means "both nodes must have fewer than 1 eliminated column" which is never true, so no merges happen.

## Decision 7: Module Organization

**Decision**: New `src/aptp/amalgamation.rs` module containing the amalgamation pass. Called from `AptpNumeric::factor()` after `build_supernode_info()`.

**Rationale**:
- The amalgamation logic is ~600-800 LOC (per ssids-plan.md estimate) — large enough to warrant its own module.
- Clean separation: `symbolic.rs` handles faer interaction, `amalgamation.rs` handles post-hoc supernode merging, `numeric.rs` handles factorization.
- The function signature: `amalgamate(supernodes: Vec<SupernodeInfo>, nemin: usize) -> Vec<SupernodeInfo>`.
