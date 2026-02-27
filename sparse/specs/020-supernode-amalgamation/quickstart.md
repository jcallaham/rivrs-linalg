# Quickstart: Supernode Amalgamation

## What This Feature Does

Adds a SPRAL-style supernode amalgamation pass that merges small supernodes after faer's symbolic analysis. This reduces assembly/extraction overhead for matrices with narrow supernodes (e.g., Schenk optimization matrices c-71, c-big).

## Key Files

| File | Role |
|------|------|
| `src/aptp/amalgamation.rs` | **NEW** — amalgamation pass (merge logic, pattern union) |
| `src/aptp/numeric.rs` | Modified — calls amalgamation after `build_supernode_info()` |
| `src/aptp/solver.rs` | Modified — adds `nemin` to `AnalyzeOptions` / `SolverOptions` |
| `src/aptp/mod.rs` | Modified — declares `amalgamation` module |

## How to Test

```bash
# Unit tests (includes new amalgamation tests)
cargo test

# Full SuiteSparse suite (includes c-71/c-big regression targets)
cargo test -- --ignored --test-threads=1

# Benchmark before/after
cargo run --example baseline_collection --features diagnostic --release -- --ci-only
```

## How to Configure

```rust
use rivrs_sparse::aptp::solver::{AnalyzeOptions, SparseLDLT};

// Default: nemin=32 (matches SPRAL)
let opts = AnalyzeOptions::default();

// More aggressive merging
let opts = AnalyzeOptions { nemin: 64, ..Default::default() };

// Disable amalgamation
let opts = AnalyzeOptions { nemin: 1, ..Default::default() };
```

## Algorithm Summary

For each parent supernode in the assembly tree, consider merging each child when:
1. **Structural match**: parent has 1 column and column count matches child's minus 1 (zero fill-in)
2. **Both small**: both parent and child have fewer than `nemin` eliminated columns

Reference: SPRAL `core_analyse.f90:806-822` (`do_merge` function).
