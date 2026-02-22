# Data Model: Supernode Amalgamation

## Entities

### SupernodeInfo (existing, modified usage)

The existing `SupernodeInfo` struct (`numeric.rs:77-86`) is the primary data model. Amalgamation transforms a `Vec<SupernodeInfo>` by merging entries and does not introduce new types.

**Fields** (unchanged):
- `col_begin: usize` — First column (inclusive) in permuted ordering
- `col_end: usize` — Past-the-end column (exclusive)
- `pattern: Vec<usize>` — Off-diagonal row indices (global permuted), sorted
- `parent: Option<usize>` — Parent supernode index, or None for root

**Invariants** (must hold after amalgamation):
- `col_begin < col_end` for all supernodes
- `pattern` entries are sorted and unique
- No pattern entry falls within `[col_begin, col_end)` (those are fully-summed)
- `parent` is `None` for exactly one supernode (root) or zero (if virtual root)
- For all non-root supernodes: `parent.unwrap() > supernode_index` (postorder)
- `col_begin..col_end` ranges may overlap with other supernodes' ranges after nemin merging (the merged range includes columns from unrelated supernodes — this is handled by scatter's column-range loop)

### Amalgamation Merge Record (internal, transient)

Used during the amalgamation pass to track which original supernodes merged into which amalgamated supernode. Not persisted beyond the pass.

**Fields**:
- `original_index: usize` — Index in the pre-amalgamation Vec
- `merged_into: usize` — Index in the post-amalgamation Vec (or self if not merged)
- `nelim: usize` — Number of eliminated columns (accumulated during merge)

### AnalyzeOptions (existing, extended)

Extended with amalgamation threshold.

**New field**:
- `nemin: usize` — Minimum supernode size for amalgamation. Default: 32.

### SolverOptions (existing, extended)

Propagates nemin to the analyze phase.

**New field**:
- `nemin: usize` — Same semantics as AnalyzeOptions. Default: 32.

## State Transitions

```
faer SymbolicCholesky
    │
    ▼
build_supernode_info()  →  Vec<SupernodeInfo>  (fundamental supernodes)
    │
    ▼
amalgamate()            →  Vec<SupernodeInfo>  (merged supernodes, fewer entries)
    │
    ▼
build_children_map()    →  children: Vec<Vec<usize>>
    │
    ▼
factor_tree_levelset()  →  Vec<FrontFactors>  (per-supernode L, D, permutation)
    │
    ▼
aptp_solve()            →  solution vector
```

## Relationships

- Each amalgamated `SupernodeInfo` contains columns from 1 or more fundamental supernodes
- Parent-child relationships form a tree (assembly tree)
- Children map is derived from parent pointers: `children[parent] = [child1, child2, ...]`
- FrontFactors are 1:1 with amalgamated supernodes (not fundamental supernodes)
