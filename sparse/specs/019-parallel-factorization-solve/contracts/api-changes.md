# API Contracts: Parallel Factorization & Solve (Phase 8.2)

**Date**: 2026-02-21

## Public API Changes

### 1. FactorOptions — add `par` field

```rust
// src/aptp/solver.rs
pub struct FactorOptions {
    pub threshold: f64,
    pub fallback: AptpFallback,
    pub outer_block_size: usize,
    pub inner_block_size: usize,
    pub par: Par,  // NEW: default Par::Seq
}
```

**Default**: `Par::Seq` (backward-compatible — existing code unchanged).

**Usage**:
```rust
let opts = FactorOptions {
    par: Par::rayon(4),  // Use 4 threads
    ..Default::default()
};
solver.factor(&matrix, &opts)?;
```

### 2. SolverOptions — add `par` field

```rust
// src/aptp/solver.rs
pub struct SolverOptions {
    pub ordering: OrderingStrategy,
    pub threshold: f64,
    pub fallback: AptpFallback,
    pub par: Par,  // NEW: default Par::Seq
}
```

### 3. SparseLDLT::solve_in_place — add `par` parameter

```rust
// Current:
pub fn solve_in_place(&self, rhs: &mut [f64]) -> Result<(), SparseError>

// New:
pub fn solve_in_place(&self, rhs: &mut [f64], par: Par) -> Result<(), SparseError>
```

**Note**: This is a breaking change to the solve signature. The `solve_full` convenience method should also accept `par` via `SolverOptions`.

### 4. aptp_solve — add `par` parameter

```rust
// Current:
pub fn aptp_solve(
    symbolic: &AptpSymbolic,
    numeric: &AptpNumeric,
    rhs: &mut [f64],
    stack: &mut MemStack,
) -> Result<(), SparseError>

// New:
pub fn aptp_solve(
    symbolic: &AptpSymbolic,
    numeric: &AptpNumeric,
    rhs: &mut [f64],
    stack: &mut MemStack,
    par: Par,
) -> Result<(), SparseError>
```

## Internal API Changes

### 5. AptpOptions — add `par` field

```rust
// src/aptp/factor.rs (internal)
pub struct AptpOptions {
    // ... existing fields ...
    pub par: Par,  // NEW: threaded to BLAS-3 calls
}
```

### 6. AptpNumeric::factor — threading par

The `par` from `FactorOptions` is:
1. Stored on `AptpOptions.par`
2. Used to decide tree-level parallelism (`Par::Seq` → sequential loop, `Par::Rayon` → rayon::scope)
3. Passed to `aptp_factor_in_place()` for fronts with dimension >= 256
4. Overridden to `Par::Seq` for fronts with dimension < 256

### 7. Per-supernode solve functions — add `par` parameter

```rust
// src/aptp/solve.rs (internal)
fn forward_solve_supernode(ff: &FrontFactors, rhs: &mut [f64], work: &mut [f64], work2: &mut [f64], par: Par)
fn diagonal_solve_supernode(ff: &FrontFactors, rhs: &mut [f64], work: &mut [f64])  // No par needed (scalar ops)
fn backward_solve_supernode(ff: &FrontFactors, rhs: &mut [f64], work: &mut [f64], work2: &mut [f64], par: Par)
```

## Backward Compatibility

- `FactorOptions::default()` has `par: Par::Seq` — existing callers unaffected
- `SolverOptions::default()` has `par: Par::Seq` — existing callers unaffected
- `SparseLDLT::solve_in_place` signature changes (breaking) — but library is pre-1.0
- `aptp_solve` signature changes (breaking) — internal API

## Cargo.toml Changes

```toml
[dependencies]
rayon = "1"  # NEW: tree-level parallelism (unconditional)
```
