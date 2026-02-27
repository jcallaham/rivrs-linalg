# Data Model: MC64 Matching & Scaling

**Feature**: 012-mc64-matching-scaling
**Date**: 2026-02-13

## Entities

### Mc64Result

The output of the MC64 matching and scaling algorithm.

**Attributes**:
- `matching: Perm<usize>` — The matching represented as a permutation vector. For row i, `matching.arrays().0[i]` gives the column σ(i) that row i is matched to. For symmetric matrices, decomposes into singletons (σ(i)=i) and 2-cycles (σ(i)=j, σ(j)=i).
- `scaling: Vec<f64>` — Symmetric scaling factors in linear domain. `scaling[i]` is the scaling factor s_i such that `A_scaled[i,j] = s_i * A[i,j] * s_j`. For the nonsingular matched block, the scaled matrix has unit diagonal and off-diagonals <= 1.
- `matched: usize` — Number of matched entries. Equals n for structurally nonsingular matrices. Less than n for structurally singular matrices.

**Relationships**:
- Matching feeds into future condensed ordering utility (not this phase)
- Scaling feeds into Phase 5 numeric factorization (`S A S` before factorize, unscale after solve)
- Matched count feeds into diagnostic reporting (structural rank information)

**Validation rules**:
- `matching` is a valid permutation of [0, n)
- `scaling.len() == n`
- All scaling factors are positive and finite
- `matched <= n`
- For matched indices: `|s_i * a_{i,σ(i)} * s_{σ(i)}| ≈ 1.0` (within floating-point tolerance)
- For all entries: `|s_i * a_{i,j} * s_j| <= 1.0` (within tolerance)

### Mc64Job (enum)

The optimization objective for the matching.

**Variants**:
- `MaximumProduct` — Maximize the product of diagonal entry magnitudes. Equivalent to minimizing the sum of `-log|a_ij|` costs. This is the default and primary objective for APTP.

**Notes**: `MaximumSum` variant may be added in a future phase. The enum is designed for extensibility.

### Internal: MatchingState

Working state during the Dijkstra augmenting path algorithm. Not exposed publicly.

**Attributes**:
- `row_match: Vec<Option<usize>>` — For each row i, the column it is currently matched to (None if unmatched)
- `col_match: Vec<Option<usize>>` — For each column j, the row it is currently matched to (None if unmatched)
- `u: Vec<f64>` — Row dual variables (log domain)
- `v: Vec<f64>` — Column dual variables (log domain)
- `d: Vec<f64>` — Dijkstra distances from current root column
- `parent: Vec<Option<usize>>` — Parent pointers for augmenting path reconstruction

**Lifecycle**:
1. Created during `mc64_matching()` call
2. Initialized via greedy heuristic (initial matching + dual variables)
3. Iteratively updated via Dijkstra augmentations
4. Consumed to produce `Mc64Result`
5. Not persisted or exposed

### Internal: CostGraph

The bipartite graph with logarithmic edge costs, derived from the input matrix.

**Attributes**:
- `col_ptr: Vec<usize>` — CSC column pointers for the full (symmetrized) matrix
- `row_idx: Vec<usize>` — CSC row indices
- `cost: Vec<f64>` — Edge costs: `c[i,j] = log(col_max_j) - log|a[i,j]|` where `col_max_j = max_k |a[k,j]|`
- `col_max: Vec<f64>` — Column maxima in log domain: `log(max_k |a[k,j]|)` for each column j

**Derivation**: Built once from the input `SparseColMat<usize, f64>` at the start of `mc64_matching()`. The upper-triangular input is expanded to full symmetric storage for bipartite graph construction.

**Validation**:
- All costs are non-negative (by construction: `c[i,j] >= 0`)
- Costs are finite for nonzero entries
- Zero entries are not stored (implicit infinite cost)

## Entity Relationships

```
SparseColMat<usize, f64> (input)
    │
    ├──► CostGraph (internal, built once)
    │        │
    │        ├──► MatchingState (internal, iteratively refined)
    │        │        │
    │        │        └──► Mc64Result (output)
    │        │                ├── matching: Perm<usize>
    │        │                ├── scaling: Vec<f64>
    │        │                └── matched: usize
    │        │
    │        └── [consumed after matching completes]
    │
    └── Mc64Job (input parameter)

Mc64Result
    ├── matching ──► [future] condensed_ordering() utility
    ├── scaling  ──► Phase 5: numeric factorization (S A S pipeline)
    └── matched  ──► diagnostic reporting (structural rank)
```
