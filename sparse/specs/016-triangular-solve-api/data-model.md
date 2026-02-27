# Data Model: Triangular Solve & Solver API

**Feature**: 016-triangular-solve-api
**Date**: 2026-02-16

## New Types

### SparseLDLT (user-facing solver)

The high-level sparse symmetric indefinite solver struct.

| Field | Type | Description |
|-------|------|-------------|
| symbolic | AptpSymbolic | Reusable symbolic analysis (from Phase 3) |
| numeric | Option\<AptpNumeric\> | Numeric factors (None until factor() called) |
| scaling | Option\<Vec\<f64\>\> | MC64 scaling in elimination order (None if no scaling) |

**State transitions**:
- `analyze()` → symbolic populated, numeric = None → "Analyzed"
- `factor()` → numeric populated → "Factored" (can solve)
- `refactor()` → numeric replaced → "Factored" (updated)
- `solve()` requires "Factored" state, else error

**Invariants**:
- `symbolic` is always populated after construction
- `scaling.as_ref().map_or(true, |s| s.len() == symbolic.nrows())`
- `numeric.as_ref().map_or(true, |n| n.n() == symbolic.nrows())`

### AnalyzeOptions

Configuration for the symbolic analysis phase.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| ordering | OrderingStrategy | OrderingStrategy::Metis | Fill-reducing ordering to use |

### OrderingStrategy (enum)

| Variant | Description |
|---------|-------------|
| Amd | AMD ordering (faer built-in) |
| Metis | METIS ordering (via metis-sys) |
| MatchOrderMetis | MC64 matching + METIS (produces scaling factors) |
| UserSupplied(Perm\<usize\>) | User-provided ordering permutation |

### FactorOptions

Configuration for the numeric factorization phase.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| threshold | f64 | 0.01 | APTP pivot threshold (u parameter) |
| fallback | AptpFallback | AptpFallback::BunchKaufman | Fallback strategy for failed 1x1 pivots |

**Note**: This wraps the existing `AptpOptions` from Phase 5, adding no new fields.

### SolverOptions

Configuration for one-shot `solve_full()`. Combines analysis and factorization options.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| ordering | OrderingStrategy | OrderingStrategy::Metis | Fill-reducing ordering |
| threshold | f64 | 0.01 | APTP pivot threshold |
| fallback | AptpFallback | AptpFallback::BunchKaufman | Pivot fallback strategy |

## Modified Types

### MixedDiagonal (diagonal.rs)

**Change**: `solve_in_place()` method modified to handle zero pivots gracefully.

| Pivot condition | Current behavior | New behavior |
|----------------|-----------------|--------------|
| 1x1, d == 0.0 | debug_assert panic | x[col] = 0.0 |
| 2x2, det == 0.0 | debug_assert panic | x[col] = 0.0, x[partner] = 0.0 |
| Delayed column | unreachable! panic | unchanged (programming error) |

### AptpNumeric::factor() (numeric.rs)

**Change**: Add `scaling: Option<&[f64]>` parameter.

| Parameter | Type | Description |
|-----------|------|-------------|
| symbolic | &AptpSymbolic | Symbolic analysis |
| matrix | &SparseColMat\<usize, f64\> | Input matrix (lower triangle) |
| options | &AptpOptions | APTP configuration |
| scaling | Option\<&[f64]\> | **NEW**: Scaling factors in elimination order |

The scaling is applied in `scatter_original_entries` only. The APTP kernel, extend-add, and all other logic remain unchanged.

### SparseError (error.rs)

**New variant**:

| Variant | Fields | Description |
|---------|--------|-------------|
| SolveBeforeFactor | context: String | Solve called before factor() |

## Existing Types Used (no changes)

| Type | Module | Role in Phase 7 |
|------|--------|-----------------|
| AptpSymbolic | aptp/symbolic.rs | Provides perm(), supernode structure, n_supernodes() |
| AptpNumeric | aptp/numeric.rs | Provides front_factors(), stats(), n() |
| FrontFactors | aptp/numeric.rs | Per-supernode L11, D11, L21, col_indices, row_indices |
| FactorizationStats | aptp/numeric.rs | max_front_size (for workspace), zero_pivots |
| Inertia | aptp/inertia.rs | Returned by SparseLDLT::inertia() |
| AptpOptions | aptp/factor.rs | Wrapped by FactorOptions |
| AptpFallback | aptp/factor.rs | BunchKaufman / Delay enum |
| SparseError | error.rs | Error type for all operations |

## Relationships

```
SparseLDLT
  ├── AptpSymbolic (always present)
  │     ├── perm() → Option<PermRef>
  │     ├── supernode structure (for solve traversal)
  │     └── nrows() → matrix dimension
  ├── Option<AptpNumeric> (after factor)
  │     ├── front_factors() → &[FrontFactors]
  │     │     ├── l11() → &Mat<f64>
  │     │     ├── d11() → &MixedDiagonal
  │     │     ├── l21() → &Mat<f64>
  │     │     ├── col_indices() → &[usize]
  │     │     └── row_indices() → &[usize]
  │     └── stats() → &FactorizationStats
  └── Option<Vec<f64>> (scaling, after analyze with MC64)

aptp_solve(symbolic, numeric, rhs, stack)
  reads: AptpSymbolic (supernode traversal)
  reads: AptpNumeric::front_factors (per-supernode factors)
  mutates: rhs (in-place solve)
  uses: MemStack (temporary workspace)
```
