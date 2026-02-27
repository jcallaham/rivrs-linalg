# API Contract: MC64 Matching & Scaling

**Feature**: 012-mc64-matching-scaling
**Date**: 2026-02-13

## Public API

### Module: `crate::aptp::matching`

Re-exported from `crate::aptp`: `mc64_matching`, `Mc64Result`, `Mc64Job`

---

### `mc64_matching`

Compute MC64 weighted bipartite matching and symmetric scaling for a sparse symmetric matrix.

**Signature** (notional):
```rust
pub fn mc64_matching(
    matrix: &SparseColMat<usize, f64>,
    job: Mc64Job,
) -> Result<Mc64Result, SparseError>
```

**Parameters**:
- `matrix`: Sparse symmetric matrix in upper-triangular CSC format (faer convention). Must be square. Numeric values are required (matching depends on entry magnitudes).
- `job`: Optimization objective. Currently only `Mc64Job::MaximumProduct`.

**Returns**:
- `Ok(Mc64Result)`: Matching permutation, scaling factors, and match count.
- `Err(SparseError::InvalidInput)`: Matrix is not square, has zero dimension, or has invalid structure.
- `Err(SparseError::AnalysisFailure)`: Internal algorithm failure (should not occur for valid input).

**Preconditions**:
- Matrix must be square (`nrows == ncols`)
- Matrix must be in upper-triangular CSC format (only upper triangle stored)
- Matrix entries must be finite (no NaN or Inf)

**Postconditions**:
- `result.matching` is a valid permutation of [0, n)
- `result.scaling` has length n, all entries positive and finite
- `result.matched <= n`
- For structurally nonsingular matrices: `result.matched == n`
- Scaling constraint: for all stored entries (i,j), `|s_i * a_ij * s_j| <= 1.0` (within tolerance)
- Diagonal constraint: for matched entries, `|s_i * a_{i,σ(i)} * s_{σ(i)}| ≈ 1.0`

**Complexity**: O(n(tau + n) log n) worst case, where n = dimension, tau = nnz. Practical cost ~O(0.2n(tau + n) log n) due to greedy initialization.

**Algorithm reference**: Duff & Koster (2001) Algorithm MPD; Duff & Pralet (2005) MC64SYM.

---

### `Mc64Result`

Result of MC64 matching-based preprocessing.

```rust
pub struct Mc64Result {
    /// Matching permutation: row i is matched to column σ(i).
    /// For symmetric matrices: decomposes into singletons (σ(i)=i)
    /// and 2-cycles (σ(i)=j, σ(j)=i).
    /// NOTE: This is NOT a fill-reducing ordering. Do not use directly
    /// with SymmetricOrdering::Custom.
    pub matching: Perm<usize>,

    /// Symmetric scaling factors (linear domain).
    /// A_scaled[i,j] = scaling[i] * A[i,j] * scaling[j]
    /// The scaled matrix has unit diagonal and off-diagonals <= 1
    /// for matched entries.
    pub scaling: Vec<f64>,

    /// Number of matched entries. Equals n for structurally nonsingular
    /// matrices. Less than n indicates structural singularity.
    pub matched: usize,
}
```

---

### `Mc64Job`

Optimization objective for the matching.

```rust
pub enum Mc64Job {
    /// Maximize the product of diagonal entry magnitudes.
    /// Equivalent to minimizing sum of -log|a_ij| costs.
    /// Default for APTP preprocessing.
    MaximumProduct,
}
```

---

## Internal Functions (not public, tested via unit tests)

### `build_cost_graph`

Build the bipartite cost graph from a symmetric sparse matrix.

```rust
fn build_cost_graph(matrix: &SparseColMat<usize, f64>) -> CostGraph
```

Expands upper-triangular CSC to full symmetric storage and computes logarithmic costs: `c[i,j] = log(col_max_j) - log|a[i,j]|`.

### `greedy_initial_matching`

Compute initial matching and dual variables using the greedy heuristic from Duff & Koster (2001).

```rust
fn greedy_initial_matching(graph: &CostGraph) -> MatchingState
```

Achieves ~80% matching cardinality on typical matrices.

### `dijkstra_augment`

Find shortest augmenting path from an unmatched column using Dijkstra on reduced costs.

```rust
fn dijkstra_augment(
    root_col: usize,
    graph: &CostGraph,
    state: &mut MatchingState,
) -> bool
```

Returns `true` if augmenting path found and matching updated, `false` if no path exists (contributes to structural singularity detection).

### `symmetrize_scaling`

Compute symmetric scaling factors from dual variables.

```rust
fn symmetrize_scaling(u: &[f64], v: &[f64]) -> Vec<f64>
```

Computes `D[i] = exp(-(u[i] + v[i]) / 2)` (MC64SYM from Duff & Pralet 2005).

### `duff_pralet_correction`

Apply scaling correction for unmatched indices in structurally singular matrices.

```rust
fn duff_pralet_correction(
    matrix: &SparseColMat<usize, f64>,
    scaling: &mut [f64],
    matched: &[bool],
)
```

For unmatched index i: `scaling[i] = 1.0 / max_k |a[i,k] * scaling[k]|` over matched k. Convention: `1/0 = 1.0`.

---

## Error Contracts

| Condition | Error | Context |
|-----------|-------|---------|
| Matrix not square | `NotSquare` | `dims: (m, n)` — uses dedicated variant with structured fields |
| Zero dimension | `InvalidInput` | "MC64 requires non-empty matrix" |
| NaN/Inf in values | `InvalidInput` | "MC64 requires finite matrix entries" |
| Internal failure | `AnalysisFailure` | "MC64 augmentation failed unexpectedly" |

Structural singularity is NOT an error — it produces `matched < n` in the result.

---

## Usage Examples

### Basic matching and scaling
```rust
use rivrs_sparse::aptp::{mc64_matching, Mc64Job};

let matrix = load_symmetric_matrix("example.mtx")?;
let result = mc64_matching(&matrix, Mc64Job::MaximumProduct)?;

println!("Matched {}/{} entries", result.matched, matrix.nrows());
println!("Scaling factors: {:?}", &result.scaling[..5]);
```

### Independent MC64 + METIS preprocessing
```rust
use rivrs_sparse::aptp::{mc64_matching, metis_ordering, AptpSymbolic, Mc64Job};
use faer::sparse::linalg::cholesky::SymmetricOrdering;

// MC64 for scaling (numerical quality)
let mc64 = mc64_matching(&matrix, Mc64Job::MaximumProduct)?;

// METIS for ordering (structural quality) — independent of MC64
let metis_perm = metis_ordering(matrix.symbolic())?;

// Symbolic analysis with METIS ordering
let symbolic = AptpSymbolic::analyze(
    matrix.symbolic(),
    SymmetricOrdering::Custom(metis_perm.as_ref()),
)?;

// Scaling factors stored for Phase 5 numeric factorization
// let numeric = factorize(&matrix, &symbolic, Some(&mc64.scaling))?;
```
