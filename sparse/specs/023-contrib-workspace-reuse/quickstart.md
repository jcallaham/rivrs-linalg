# Quickstart: Contribution Workspace Reuse

**Feature**: 023-contrib-workspace-reuse
**Date**: 2026-02-23

## Overview

This feature eliminates per-supernode heap allocation for contribution blocks in the multifrontal factorization loop. Phase 9.1c profiling identified contribution block allocation as consuming 73% of factorization time on large matrices (40.1% extraction + 33.3% extend-add). The optimization uses two complementary approaches:

1. **Contribution buffer pool**: Pre-allocated reusable dense matrices replace per-supernode `Mat::zeros` calls
2. **Direct extend-add from frontal workspace**: For the last child of each parent in DFS postorder, skip the contribution copy entirely and read the Schur complement in-place

## Key Files

| File | Changes |
|------|---------|
| `src/aptp/numeric.rs` | FactorizationWorkspace (new fields), contribution pool, DFS sequential path, direct extend-add, ContributionBlock pool lifecycle |
| `examples/profile_matrix.rs` | Updated sub-phase timing display for new code paths |
| `examples/baseline_collection.rs` | Updated baseline format for new timing fields |

## Build & Test

```bash
# Build (standard)
cargo build

# Build with diagnostic timing
cargo build --features diagnostic

# Unit tests
cargo test

# SuiteSparse CI subset (10 matrices)
cargo test -- --test-threads=1

# Full SuiteSparse suite (65 matrices, requires test-data/suitesparse/)
cargo test -- --ignored --test-threads=1

# Profile c-71 (requires diagnostic + suitesparse data)
cargo run --example profile_matrix --features diagnostic --release -- c-71

# Baseline comparison
cargo run --example baseline_collection --features diagnostic --release -- --compare target/benchmarks/baselines/prev.json
```

## Validation Checklist

1. `cargo test` — all unit tests pass
2. `cargo test -- --ignored --test-threads=1` — all 65 SuiteSparse matrices pass (backward error < 5e-11)
3. `cargo clippy --all-targets --features diagnostic` — no warnings
4. `cargo fmt --check` — formatting clean
5. Profile c-71: factor time ≤ 1.5× SPRAL, sys time < 10%
6. Baseline comparison: no matrix regresses > 5% vs Phase 9.1c

## Architecture Summary

### Sequential Path (DFS Postorder)

```
factor_subtree_dfs(s, active_buf, accum_buf):
  for each child c of s (except last):
    factor_subtree_dfs(c, active_buf, accum_buf)
    # c's contribution extracted into pool buffer
    # extend-add pool buffer into accum_buf (parent's frontal)
    # return pool buffer

  for last child c of s:
    factor_subtree_dfs(c, active_buf, accum_buf)
    # direct extend-add from active_buf into accum_buf (zero-copy)

  scatter original entries for s into accum_buf
  factor s in accum_buf
  swap(active_buf, accum_buf)  # s's data is now in "active"
```

### Parallel Path (Wave-based, unchanged traversal)

```
for each wave:
  par_iter over ready supernodes:
    factor_single_supernode(s, child_contribs, workspace)
    # workspace.contribution_pool used for extraction
    # ContributionBlock still owned (crosses thread boundaries)
    # pool buffers recycled within thread-local workspace
```
