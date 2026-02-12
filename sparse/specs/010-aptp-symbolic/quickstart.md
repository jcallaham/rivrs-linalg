# Quickstart: APTP Symbolic Analysis

**Feature**: 010-aptp-symbolic
**Date**: 2026-02-10

## What This Feature Does

AptpSymbolic performs the symbolic analysis step of the SSIDS sparse solver pipeline. Given a sparse symmetric matrix's sparsity pattern, it computes the elimination tree, predicts the factor structure, applies a fill-reducing ordering, and estimates APTP-specific pivot buffer requirements. This result is reusable across multiple numeric factorizations with the same sparsity pattern.

## Files to Create/Modify

### New Files

| File | Purpose |
|------|---------|
| `src/aptp/symbolic.rs` | AptpSymbolic struct, SymbolicStatistics, analyze() constructor |

### Modified Files

| File | Change |
|------|--------|
| `src/aptp/mod.rs` | Add `pub mod symbolic;` and re-export `AptpSymbolic`, `SymbolicStatistics` |
| `src/error.rs` | May need new error variant or context enrichment for analysis failures (evaluate during implementation) |
| `tests/` | New integration test file for symbolic analysis on test matrices |
| `benches/solver_benchmarks.rs` | Add symbolic analysis benchmark |

## Development Workflow (TDD)

Following constitution Principle III (Test-Driven Development):

### Step 1: Write failing tests
- Unit test: analyze 1×1 matrix → valid trivial result
- Unit test: analyze diagonal matrix → predicted_nnz == dimension
- Unit test: analyze hand-constructed arrow matrix → valid statistics
- Unit test: custom ordering passthrough
- Unit test: determinism (same input → same output)
- Unit test: invalid permutation → descriptive error
- Integration test: analyze all hand-constructed matrices with AMD
- Integration test: analyze SuiteSparse CI subset with AMD
- Property test: pivot buffer estimates non-negative for all matrices

### Step 2: Implement AptpSymbolic
- Define struct with inner SymbolicCholesky + etree + col_counts + pivot_buffer
- Implement `analyze()` constructor: validate → prefactorize → factorize_symbolic → buffer estimation
- Implement accessor methods (delegating to faer where possible)
- Implement SymbolicStatistics

### Step 3: Validate
- All tests pass
- `cargo clippy` clean
- `cargo fmt --check` clean
- `cargo doc` clean (all public items documented)
- Run on full test suite

## Usage Pattern (Notional)

```rust
use rivrs_sparse::aptp::AptpSymbolic;

// Load a sparse symmetric matrix
let matrix: SparseColMat<usize, f64> = /* ... */;

// Perform symbolic analysis with default AMD ordering
let symbolic = AptpSymbolic::analyze(
    matrix.symbolic(),
    SymmetricOrdering::Amd,
)?;

// Inspect results
let stats = symbolic.statistics();
println!("Dimension: {}", stats.dimension);
println!("Predicted NNZ: {}", stats.predicted_nnz);
println!("Supernodal: {}", stats.is_supernodal);

// Access permutation
if let Some(perm) = symbolic.perm() {
    let (fwd, inv) = perm.arrays();
    // ...
}

// Access elimination tree
let etree = symbolic.etree();  // parent pointer array

// Later: pass to numeric factorization (Phase 5-6)
// let numeric = AptpNumeric::factorize(&symbolic, &matrix)?;
```

## Key Design Decisions

1. **Transparent composition**: AptpSymbolic stores faer's SymbolicCholesky directly, delegates standard queries, adds only APTP-specific fields
2. **Two-call strategy**: `prefactorize_symbolic_cholesky` (for etree + col_counts) + `factorize_symbolic_cholesky` (for full symbolic result) — because faer's high-level API doesn't expose the etree publicly
3. **Dual-variant handling**: Must handle both simplicial and supernodal results from faer (pattern matching on `raw()`)
4. **10% buffer heuristic**: Initial estimate for delayed pivot workspace, to be refined empirically in Phase 5-6

## Dependencies

No new dependencies required. Uses existing:
- `faer 0.22` — symbolic Cholesky, SymmetricOrdering, PermRef, MemStack
- `rivrs_sparse::error::SparseError` — error handling
