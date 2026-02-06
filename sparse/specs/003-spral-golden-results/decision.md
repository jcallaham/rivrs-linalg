# Decision Record: Defer Phase 0.3 (SPRAL Golden Results)

**Feature Branch**: `003-spral-golden-results`
**Date**: 2026-02-06
**Decision**: Defer Phase 0.3. No SPRAL golden results infrastructure built now.

## Context

Phase 0.3 of the SSIDS development plan called for building SPRAL from source,
writing a C driver program, running SPRAL on all 82 test matrices, and capturing
comprehensive reference outputs (timing, inertia, pivot statistics, errors) as
JSON files.

During planning and API investigation, we assessed the cost-benefit of this work.

## Analysis

### What SPRAL's C API actually exposes

After investigation of `spral_ssids.h` and the SPRAL source:

- **Aggregate statistics only** from `inform` struct: num_factor, num_flops,
  num_delay, num_neg, num_two, matrix_rank, maxfront, maxsupernode, maxdepth
- **Block diagonal D** via `spral_ssids_enquire_indef()`: d[2*n] array with
  1x1 and 2x2 pivot entries
- **Pivot ordering** via `spral_ssids_enquire_indef()`: piv_order[n]
- **Solution vector** after solve (from which errors are computed)

**NOT exposed**: L factor, permutation P, elimination tree parent pointers,
supernode membership arrays, fill-in pattern. The actual factorization matrices
are internal to SPRAL's opaque `akeep`/`fkeep` handles.

### Cost

- Building SPRAL from source on arm64 (Fortran/C/meson, dependency chain)
- ~500 LOC of throwaway C driver code
- Shell orchestration scripts
- Dockerfile modifications
- Debugging potential arm64 build issues

### What golden results would provide

1. **Reference backward error** — but we can compute this independently
2. **Reference inertia for large matrices** — useful but not needed until Phases 2-8
3. **Difficulty characterization** — informative but not blocking
4. **Timing baseline** — arm64 vs SPRAL's x86 optimization makes this unreliable

### What we can validate WITHOUT golden results

The **strongest possible correctness test** is the **reconstruction test**:

> `||P^T A P - L D L^T|| / ||A|| < epsilon`

This proves mathematical correctness by definition — if our factorization
reconstructs A to machine precision, it is correct regardless of what SPRAL
would produce. This is strictly stronger than comparing against another solver.

Additional SPRAL-independent validation:
- **Backward error**: `||Ax - b|| / (||A|| ||x|| + ||b||)` — computed from our
  own solution, no reference needed
- **Hand-constructed matrices**: Phase 0.2 provides 15 matrices with analytically
  known LDL^T factorizations, inertia, and properties
- **Property-based tests**: Symmetry preservation, inertia consistency,
  permutation validity

## Decision

**Defer Phase 0.3** to a later point when we have a working Rust solver and
need SPRAL comparison for:
- Performance benchmarking (within 2x target)
- Inertia validation on large SuiteSparse matrices
- Cross-solver comparison for confidence on hard indefinite cases

## Impact on Development Plan

- Phase 0 exit criteria updated: remove "SPRAL reference results" requirement,
  replace with reconstruction test strategy
- Constitution updated: reconstruction tests are the primary correctness oracle;
  SPRAL comparison is a secondary validation layer added when the solver is ready
- Testing strategy: reconstruction (`P^T A P = L D L^T`) + backward error +
  hand-constructed matrices + property-based tests

## When to Revisit

Build SPRAL infrastructure when:
1. The Rust solver can factorize SuiteSparse matrices (Phases 2-8)
2. We need inertia validation on matrices too large for analytical verification
3. We want performance benchmarking against a production solver
