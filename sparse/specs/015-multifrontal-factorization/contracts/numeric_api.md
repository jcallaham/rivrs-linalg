# API Contract: AptpNumeric (Public)

**Module**: `src/aptp/numeric.rs`
**Visibility**: `pub`

## Entry Point

```
AptpNumeric::factor(
    symbolic: &AptpSymbolic,
    matrix: &SparseColMat<usize, f64>,
    options: &AptpOptions,
) -> Result<AptpNumeric, SparseError>
```

### Inputs

| Parameter | Type | Constraints |
|-----------|------|-------------|
| symbolic | &AptpSymbolic | Must be the result of analyzing a matrix with the same sparsity pattern as `matrix` |
| matrix | &SparseColMat\<usize, f64\> | Sparse symmetric matrix (lower triangle stored). Dimensions must match `symbolic.nrows()` |
| options | &AptpOptions | APTP configuration (threshold, fallback strategy). Shared across all supernodes |

### Output

| Field | Type | Description |
|-------|------|-------------|
| front_factors | Vec\<FrontFactors\> | Per-supernode factors, indexed by supernode ID |
| stats | FactorizationStats | Aggregate statistics |
| n | usize | Matrix dimension |

### Errors

| Error | Condition |
|-------|-----------|
| SparseError::DimensionMismatch | Matrix dimensions don't match symbolic analysis |
| SparseError::NumericalSingularity | Root supernode has unresolvable delayed columns |
| SparseError::AnalysisFailure | Symbolic analysis is inconsistent (internal error) |

### Postconditions

1. `front_factors.len() == n_supernodes` (one entry per supernode)
2. For each s: `front_factors[s].num_eliminated <= supernode_ncols(s) + delayed_from_children`
3. `stats.total_1x1_pivots + 2 * stats.total_2x2_pivots + root_delayed == n`
4. Reconstruction error `||P^T A P - L D L^T|| / ||A|| < 10^-12` (verified by tests, not enforced at runtime)

---

## Supporting Types

### FrontFactors

```
pub struct FrontFactors {
    l11: Mat<f64>,
    d11: MixedDiagonal,
    l21: Mat<f64>,
    local_perm: Vec<usize>,
    num_eliminated: usize,
    col_indices: Vec<usize>,
    row_indices: Vec<usize>,
}
```

Accessor methods: `l11()`, `d11()`, `l21()`, `local_perm()`, `num_eliminated()`, `col_indices()`, `row_indices()`

### FactorizationStats

```
pub struct FactorizationStats {
    pub total_1x1_pivots: usize,
    pub total_2x2_pivots: usize,
    pub total_delayed: usize,
    pub max_front_size: usize,
}
```

---

## Usage Example

```
let symbolic = AptpSymbolic::analyze(
    matrix.symbolic(),
    SymmetricOrdering::Amd,
)?;

let numeric = AptpNumeric::factor(
    &symbolic,
    &matrix,
    &AptpOptions::default(),
)?;

println!("Pivots: {} 1x1, {} 2x2, {} delayed",
    numeric.stats().total_1x1_pivots,
    numeric.stats().total_2x2_pivots,
    numeric.stats().total_delayed);
println!("Max front size: {}", numeric.stats().max_front_size);
```
