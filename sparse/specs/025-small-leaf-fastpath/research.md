# Research: Small Leaf Subtree Fast Path

**Date**: 2026-02-24
**Feature**: 025-small-leaf-fastpath

## R1: Root Cause of Simplicial Matrix Slowdown

**Decision**: The dominant cost on simplicial matrices is per-supernode dispatch overhead, not algorithmic complexity.

**Rationale**: Phase 9.1c profiling of dixmaanl (max_front=68, 1,605 supernodes post-amalgamation) shows 23.3% "unaccounted" time — function call overhead, workspace setup, g2l construction/teardown, and assembly maps lookup. Each supernode processes a tiny (<68×68) frontal matrix where the actual BLAS work takes microseconds, but the per-node setup is fixed-cost. First supernode alone takes 6.6ms (18.3% of total 47ms factor time), likely first-use/setup overhead.

For comparison, c-71 (max_front=2475, 6.4K supernodes post-amalgamation) spends 73.4% of factor time in extend-add + contribution extraction — dominated by large O(n²) copies. The bottleneck is fundamentally different for small-front vs large-front matrices.

**Alternatives considered**:
- Reducing per-supernode overhead in the general path (e.g., cheaper g2l) — helps incrementally but doesn't eliminate the dispatch cost for thousands of tiny nodes
- Simplicial LDL^T via faer (column-by-column, no frontal matrices) — would bypass the overhead entirely but can't handle APTP pivoting (indefinite systems)

## R2: SPRAL's SmallLeafNumericSubtree Architecture

**Decision**: Adapt SPRAL's two-phase assembly pattern (assemble_pre → factor → assemble_post) but retain our `aptp_factor_in_place()` kernel with a reused small buffer.

**Rationale**: SPRAL's indefinite small-leaf path (`SmallLeafNumericSubtree.hxx:187-446`) uses three steps per node:
1. `assemble_pre`: Scatter original entries + delayed columns from children into lcol (direct scatter, no frontal matrix)
2. `factor_node`: Call `ldlt_tpp_factor()` in-place on lcol, compute contribution via GEMM
3. `assemble_post`: Scatter child contributions into the generated element (contribution block)

SPRAL scatters directly into the permanent factor storage (lcol) — no intermediate frontal buffer. Our architecture uses `aptp_factor_in_place(MatMut)` which expects a dense matrix buffer. Rather than rewriting the APTP kernel to scatter directly into factor storage (a deep refactoring), we use a single small reusable buffer (bounded by threshold² = 256² = 512KB) that stays cache-resident across all nodes in the subtree.

**Key SPRAL patterns to adopt**:
- Sequential processing within a subtree (no parallelism overhead)
- Single workspace reused across all nodes in the subtree
- Subtree processed as a unit before the main level-set loop
- `insmallleaf` flag on each symbolic node to skip during main loop
- Contribution block at subtree root flows to parent via normal extend-add

**Patterns NOT adopted**:
- Direct scatter into lcol (would require APTP kernel redesign)
- SPRAL's two-phase assembly split (assemble_pre/assemble_post) — our code already handles the fully-summed + generated element assembly in a single pass
- Flops-based threshold (using front-size-based per spec decision)

## R3: Fast Path Integration Point

**Decision**: Factor small-leaf subtrees as a pre-pass before the level-set loop in `factor_tree_levelset()`.

**Rationale**: The level-set loop (`numeric.rs:1478-1618`) starts by identifying leaf supernodes with `remaining_children[s] == 0`. A pre-pass can:
1. Identify all small-leaf subtrees
2. Factor them sequentially with a lightweight workspace
3. Store their results (FrontFactors + contribution blocks)
4. Decrement `remaining_children` for parents of subtree roots
5. Mark subtree nodes as already-factored

The main level-set loop then starts with an adjusted ready set that skips pre-factored nodes. This approach:
- Separates concerns cleanly (fast path vs general path)
- Avoids conditional dispatch inside the per-supernode hot loop
- Matches SPRAL's architecture (small-leaf subtrees processed before main task DAG)
- Doesn't require changes to `factor_single_supernode()`

**Alternatives considered**:
- Inline dispatch in level-set loop (conditional in hot path — less clean, mixes concerns)
- Recursive DFS within subtree (risk stack overflow on deep trees like dixmaanl)
- New iterative loop within subtree (adds complexity without benefit — subtrees are small enough for simple postorder iteration)

## R4: Workspace Strategy for Fast Path

**Decision**: Use a dedicated lightweight workspace per subtree, separate from the general FactorizationWorkspace pool.

**Rationale**: The general `FactorizationWorkspace` (numeric.rs:74-169) allocates `max_front × max_front` for `frontal_data`, which for c-71 means 2475² × 8 = 49MB. For small-leaf subtrees, the maximum front size is bounded by the threshold (256), so the workspace is at most 256² × 8 = 512KB — fits in L2 cache.

Key workspace components for the fast path:
- `frontal_data: Mat<f64>` — fixed at `max_front_in_subtree × max_front_in_subtree` (≤ 512KB)
- `frontal_row_indices: Vec<usize>` — capacity `max_front_in_subtree`
- `delayed_cols_buf: Vec<usize>` — capacity bounded by threshold
- `global_to_local: Vec<usize>` — length n (shared, same as general path)
- `contrib_buffer: Mat<f64>` — recycled across nodes within subtree

The workspace is allocated once per subtree and reused across all nodes. No Mutex contention (sequential processing). The g2l array is the only n-length allocation; it can be shared with the general workspace or allocated separately.

**Alternatives considered**:
- Reuse the general FactorizationWorkspace (over-allocates for small fronts, defeats cache locality)
- No workspace at all / SPRAL-style direct scatter into factor storage (requires APTP kernel redesign)
- Per-node allocation (what we're trying to avoid)

## R5: Contribution Block Compatibility

**Decision**: The fast path produces standard `ContributionBlock` structs identical to the general path.

**Rationale**: The `ContributionBlock` struct (`numeric.rs:479-486`) has three fields: `data: Mat<f64>`, `row_indices: Vec<usize>`, `num_delayed: usize`. The parent's extend-add (`extend_add()` or `extend_add_mapped()`) consumes these via the same interface regardless of how they were produced.

Within a small-leaf subtree, child→parent contributions are consumed inline (no need to store in the global contributions vector until the subtree root). At the subtree root, the final contribution block is stored in the global vector for the parent's extend-add in the main level-set loop.

**Key compatibility requirement**: The subtree root's contribution must have `row_indices` in the same global permuted index space, and `num_delayed` must correctly account for any columns that cascaded through the subtree without being eliminated.

## R6: Classification Algorithm

**Decision**: Bottom-up pass over post-amalgamation supernodes, marking `in_small_leaf: bool`.

**Rationale**: After amalgamation produces the final `Vec<SupernodeInfo>` and `children_map`, a single O(n_supernodes) pass identifies small-leaf subtrees:

```
For each supernode s in postorder:
    front_size(s) = owned_cols(s) + pattern(s).len()
    if front_size(s) >= threshold:
        in_small_leaf[s] = false
        continue
    if s is a leaf (no children):
        in_small_leaf[s] = true  // candidate (needs ≥2 nodes to qualify)
    else if all children have in_small_leaf = true:
        in_small_leaf[s] = true
    else:
        in_small_leaf[s] = false
```

Then, identify subtree roots: a node with `in_small_leaf = true` whose parent has `in_small_leaf = false` (or no parent). Remove single-node subtrees by checking if the root has no in_small_leaf children.

This is O(n_supernodes) and runs once after amalgamation. The `in_small_leaf` flag is stored on `SupernodeInfo` and reused across refactorizations.
