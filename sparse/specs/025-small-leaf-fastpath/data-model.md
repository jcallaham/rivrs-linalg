# Data Model: Small Leaf Subtree Fast Path

**Date**: 2026-02-24
**Feature**: 025-small-leaf-fastpath

## Modified Types

### SupernodeInfo (numeric.rs)

**Change**: Add `in_small_leaf` classification flag.

| Field | Type | Status | Description |
|-------|------|--------|-------------|
| col_begin | usize | Existing | First column index |
| col_end | usize | Existing | Past-end column index |
| pattern | Vec\<usize\> | Existing | Off-diagonal row indices |
| parent | Option\<usize\> | Existing | Parent supernode index |
| owned_ranges | Vec\<Range\<usize\>\> | Existing | Column ranges post-amalgamation |
| **in_small_leaf** | **bool** | **New** | True if this supernode belongs to a small-leaf subtree |

**Invariants**:
- If `in_small_leaf = true`, then `front_size(s) < small_leaf_threshold`
- If `in_small_leaf = true` and all children also have `in_small_leaf = true`, this node may be a subtree root or interior node
- A subtree root has `in_small_leaf = true` and either `parent = None` or `parent.in_small_leaf = false`
- Every small-leaf subtree contains ≥ 2 supernodes

### FactorOptions (solver.rs)

**Change**: Add small-leaf threshold configuration.

| Field | Type | Status | Default | Description |
|-------|------|--------|---------|-------------|
| threshold | f64 | Existing | 0.01 | APTP pivot threshold |
| fallback | AptpFallback | Existing | BunchKaufman | Fallback strategy |
| outer_block_size | usize | Existing | 256 | Two-level blocking |
| inner_block_size | usize | Existing | 32 | TPP switching point |
| par | Par | Existing | Par::Seq | Parallelism control |
| nemin | usize | Existing | 32 | Amalgamation threshold |
| **small_leaf_threshold** | **usize** | **New** | **256** | Front-size threshold for small-leaf subtrees. 0 = disabled. |

### SolverOptions (solver.rs)

**Change**: Mirror the new threshold.

| Field | Type | Status | Default | Description |
|-------|------|--------|---------|-------------|
| ... | ... | Existing | ... | ... |
| **small_leaf_threshold** | **usize** | **New** | **256** | Mirrors FactorOptions |

## New Types

### SmallLeafSubtree

Represents a classified small-leaf subtree for the fast path. Produced during classification, consumed during factorization.

| Field | Type | Description |
|-------|------|-------------|
| root | usize | Supernode ID of the subtree root (topmost node) |
| nodes | Vec\<usize\> | Supernode IDs in postorder (leaves first, root last) |
| max_front_size | usize | Maximum front size across all nodes in the subtree |
| parent_of_root | Option\<usize\> | Parent supernode (outside subtree) that receives the root's contribution |

**Invariants**:
- `nodes.len() >= 2`
- `nodes.last() == Some(&root)`
- All nodes have `in_small_leaf = true`
- `max_front_size < small_leaf_threshold`
- Nodes are in postorder: children appear before parents

### SmallLeafWorkspace

Lightweight workspace for processing a small-leaf subtree. Allocated once per subtree, reused across all nodes.

| Field | Type | Description |
|-------|------|-------------|
| frontal_data | Mat\<f64\> | Reusable frontal buffer (max_front × max_front for this subtree) |
| frontal_row_indices | Vec\<usize\> | Row indices for current node (capacity max_front) |
| delayed_cols_buf | Vec\<usize\> | Delayed columns from children |
| global_to_local | Vec\<usize\> | Global→local index map (length n, shared or dedicated) |
| contrib_buffer | Mat\<f64\> | Recycled contribution output buffer |

**Size bound**: `frontal_data` is at most `threshold × threshold` (256×256 = 512KB for default threshold). Fits in L2 cache.

## Unchanged Types (for reference)

These types are produced by the fast path in identical format to the general path:

- **FrontFactors**: L11, D11 (MixedDiagonal), L21, local_perm, col_indices, row_indices — per-supernode solve factors
- **ContributionBlock**: data (Mat\<f64\>), row_indices (Vec\<usize\>), num_delayed (usize) — Schur complement for parent assembly
- **PerSupernodeStats**: snode_id, front_size, num_eliminated, pivot counts, timing (diagnostic) — per-node diagnostics
- **FactorizationStats**: Aggregate statistics — accumulates across both fast and general paths

## Data Flow

```
classify_small_leaf_subtrees(supernodes, children_map, threshold)
  → Vec<SmallLeafSubtree>
  → sets in_small_leaf on SupernodeInfo

factor_small_leaf_subtree(subtree, ...)
  → Vec<(FrontFactors, Option<ContributionBlock>, PerSupernodeStats)>
  (one per node in subtree)

factor_tree_levelset(...)
  1. Pre-pass: factor all SmallLeafSubtrees → store results + root contributions
  2. Main loop: level-set factorization skipping in_small_leaf nodes
  3. Merge: combine results from both paths into final AptpNumeric
```
