# Quickstart: Assembly & Extraction Optimization (Phase 9.1c)

**Feature**: 022-assembly-extraction-opt
**Date**: 2026-02-22

## Overview

This feature optimizes the assembly and extraction phases of multifrontal LDL^T factorization. The optimizations are internal — the public API (`SparseLDLT::analyze/factor/solve`) is unchanged.

## Integration Scenarios

### Scenario 1: Normal usage (no API changes)

```rust
use rivrs_sparse::aptp::solver::SparseLDLT;

// Analyze (now precomputes scatter maps internally)
let solver = SparseLDLT::analyze_with_matrix(&matrix, options)?;

// Factor (now uses precomputed maps for assembly, bulk copies for extraction)
let solver = solver.factor(&matrix, factor_options)?;

// Solve (unchanged)
let x = solver.solve(&b)?;
```

Users see no API changes. The optimization is transparent:
- `analyze_with_matrix` takes slightly longer (computes scatter maps)
- `factor` runs faster (uses precomputed maps, bulk copies)
- `solve` is unchanged

### Scenario 2: Re-factorization with same sparsity pattern

```rust
// Analyze once
let analyzer = SparseLDLT::analyze_with_matrix(&matrix_v1, options)?;

// Factor multiple times with different values (scatter maps reused)
let solver_v1 = analyzer.factor(&matrix_v1, factor_options)?;
let solver_v2 = analyzer.factor(&matrix_v2, factor_options)?;  // maps reused
```

Scatter maps are computed once during analysis and reused across factorizations. This is the intended usage pattern for optimization solvers that solve many linear systems with the same sparsity pattern but different numerical values.

### Scenario 3: Profiling and diagnostics

```bash
# Profile a specific matrix to see assembly/extraction breakdown
cargo run --example profile_matrix --features diagnostic --release -- GHS_indef/c-71

# Expected output shows reduced assembly and extraction percentages:
# Assembly:       XXXX.XX ms  ( XX.X%)  # was 37.7%
# Kernel:         XXXX.XX ms  ( XX.X%)  # was 23.0%
# Extraction:     XXXX.XX ms  ( XX.X%)  # was 34.8%
```

### Scenario 4: Benchmark comparison

```bash
# Run full SPRAL comparison benchmark
cargo run --example spral_benchmark --release -- --threads 1 --rivrs

# Key metrics to check:
# c-71 ratio: should be < 3.0x (was 4.06x)
# c-big ratio: should be < 3.0x (was 4.19x)
# Median ratio: should be <= 1.01x
# All backward errors: < 5e-11
```

## Validation Checklist

After implementation, verify:

1. `cargo test` — all unit tests pass
2. `cargo test --features diagnostic` — diagnostic feature still compiles and works
3. `cargo test -- --ignored --test-threads=1` — all 65 SuiteSparse matrices pass
4. `cargo clippy --all-targets` — no warnings
5. `cargo clippy --all-targets --features diagnostic` — no warnings
6. `profile_matrix c-71` — assembly + extraction percentages decreased
7. `spral_benchmark` — c-71 and c-big ratios improved, no regressions > 5%
