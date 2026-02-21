# API Contracts: Sequential Profiling & Optimization (Phase 8.1g)

**Date**: 2026-02-20

Phase 8.1g does not introduce new public API entry points. All changes extend existing internal data structures and add tooling. This document describes the contracts for extended types and the baseline collection tool interface.

## Extended Types (Internal)

### PerSupernodeStats — Timing Extension

The existing `PerSupernodeStats` struct gains three optional timing fields, conditionally compiled behind the `diagnostic` feature flag.

**Contract**:
- When `diagnostic` is NOT enabled: struct is unchanged from Phase 8.1f. No binary size or memory impact.
- When `diagnostic` IS enabled: three `Duration` fields are added. These are populated during `AptpNumeric::factor()` and accessible via `SparseLDLT::per_supernode_stats()`.
- Timing fields use `std::time::Duration` (not raw nanoseconds).
- Timing values represent wall-clock time for that supernode's phase, measured by `Instant::now()` differencing.

### FactorizationStats — Timing Extension

The existing `FactorizationStats` struct gains three aggregate timing fields, conditionally compiled behind `diagnostic`.

**Contract**:
- Aggregate timing fields are the sum of corresponding per-supernode timing fields.
- Accessible via `SparseLDLT::stats()` (existing API, no signature change).

## Baseline Collection Tool

### Interface: `examples/baseline_collection.rs`

**Invocation**:
```
cargo run --example baseline_collection --features diagnostic --release
cargo run --example baseline_collection --features diagnostic --release -- --ci-only
cargo run --example baseline_collection --features diagnostic --release -- --compare <previous_baseline.json>
```

**Output**:
- JSON file to `target/benchmarks/baselines/baseline-<timestamp>.json`
- Contains `BaselineSuite` structure (see data-model.md)
- Human-readable summary to stderr

**Behavioral contract**:
- Loads all available SuiteSparse matrices (65 full, 10 CI subset with `--ci-only`)
- Uses `MatchOrderMetis` ordering (the default)
- Generates a random RHS vector for each matrix, solves, computes backward error
- Records per-phase timing, per-supernode stats, and peak RSS
- When `--compare` is given, loads previous baseline and reports regressions/improvements

### Interface: `examples/workload_analysis.rs`

**Invocation**:
```
cargo run --example workload_analysis --features diagnostic --release
```

**Output**:
- Workload distribution report to stdout
- Per-matrix: front size histogram, time contribution by front size band, parallelism recommendation
- Summary: matrix classification into TreeLevel/IntraNode/Mixed categories
- Optional: Chrome Trace files to `target/profiles/` for representative matrices

**Behavioral contract**:
- Requires `diagnostic` feature for per-supernode timing
- Processes all available SuiteSparse matrices
- Computes `WorkloadProfile` for each matrix
- Produces the Phase 8.2 parallelism strategy recommendation

## Correctness Invariants (Unchanged)

All optimization changes must preserve:
- Reconstruction error: `||P^T A P - L D L^T|| / ||A|| < 1e-12` for all 65 SuiteSparse matrices
- Backward error: no regression from Phase 8.1f baselines
- Inertia counts: unchanged from Phase 8.1f
- Pivot statistics: unchanged (same algorithm, same decisions)

Allocation optimizations are strictly internal refactoring — no behavioral change at the API boundary.
