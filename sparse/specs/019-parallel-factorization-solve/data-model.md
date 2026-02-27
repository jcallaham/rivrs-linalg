# Data Model: Parallel Factorization & Solve (Phase 8.2)

**Date**: 2026-02-21

## Modified Entities

### FactorOptions (modified)

Existing struct in `src/aptp/solver.rs`. Add `par` field.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| threshold | f64 | 0.01 | APTP pivot threshold |
| fallback | AptpFallback | BunchKaufman | Fallback strategy |
| outer_block_size | usize | 256 | Two-level outer block size |
| inner_block_size | usize | 32 | Two-level inner block size |
| **par** | **Par** | **Par::Seq** | **Parallelism control (new)** |

### SolverOptions (modified)

Existing struct in `src/aptp/solver.rs`. Add `par` field.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| ordering | OrderingStrategy | MatchOrderMetis | Fill-reducing ordering |
| threshold | f64 | 0.01 | APTP pivot threshold |
| fallback | AptpFallback | BunchKaufman | Fallback strategy |
| **par** | **Par** | **Par::Seq** | **Parallelism control (new)** |

### AptpOptions (modified)

Internal struct in `src/aptp/factor.rs`. Add `par` field for threading through to BLAS-3 calls.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| threshold | f64 | 0.01 | Pivot threshold |
| small | f64 | 1e-20 | Singularity detection |
| fallback | AptpFallback | BunchKaufman | Fallback strategy |
| outer_block_size | usize | 256 | Outer block size |
| inner_block_size | usize | 32 | Inner block size |
| failed_pivot_method | FailedPivotMethod | Tpp | Failed pivot strategy |
| **par** | **Par** | **Par::Seq** | **BLAS parallelism (new)** |

## Internal Data Flow Changes

### Factorization (numeric.rs)

**Current** (sequential postorder loop):
```
global_to_local: Vec<usize>         ← shared, reset per supernode
contributions: Vec<Option<CB>>      ← index-based, consumed by parent
front_factors_vec: Vec<FrontFactors> ← sequential push
stats: FactorizationStats           ← incremental accumulation
```

**Parallel** (recursive rayon::scope):
```
contributions: Vec<Option<CB>>      ← protected by scope barrier (child→parent)
front_factors: Vec<Option<FF>>      ← pre-allocated, index-based write
per_sn_stats: Vec<Option<Stats>>    ← pre-allocated, index-based write
stats: FactorizationStats           ← computed from per_sn_stats after completion
global_to_local                     ← per-task allocation (cannot share across parallel tasks)
```

### Solve (solve.rs)

**Current** (sequential loops):
```
work: Vec<f64>     ← reused buffer (max_front_size)
work2: Vec<f64>    ← reused buffer (max_front_size)
```

**Parallel** (diagonal par_iter + tree-level scope for fwd/bwd):
```
# Diagonal: per-supernode allocation or thread-local pool
# Forward/backward: per-task work/work2 allocation within scope closures
```

## No New Public Types Required

All changes are additions to existing structs (`FactorOptions`, `SolverOptions`, `AptpOptions`). The `Par` type comes from faer — consistent with transparent composition principle.

## Constants

| Name | Value | Description |
|------|-------|-------------|
| INTRA_NODE_THRESHOLD | 256 | Front dimension below which intra-node BLAS uses Par::Seq |

This can be a constant in `numeric.rs` or a field on `FactorOptions` (const for now, field later if tuning needed).
