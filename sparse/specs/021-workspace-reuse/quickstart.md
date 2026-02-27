# Quickstart: Workspace Reuse & Per-Supernode Allocation Optimization

**Feature**: 021-workspace-reuse
**Date**: 2026-02-22

## What This Feature Changes

This is an **internal optimization** — no changes to the public `SparseLDLT` API. Users get faster factorization without changing any code.

The optimization eliminates ~16-18 heap allocations per supernode in the multifrontal factorization loop by pre-allocating reusable workspace buffers sized to the maximum front dimension.

## How to Verify

### Run the existing test suite

```bash
cargo test                           # Unit tests (should all pass)
cargo test -- --ignored --test-threads=1   # Full SuiteSparse suite
```

### Run benchmarks to measure improvement

```bash
# Baseline collection (before and after)
cargo run --example baseline_collection --features diagnostic --release -- --ci-only

# Compare against previous baseline
cargo run --example baseline_collection --features diagnostic --release -- --compare target/benchmarks/baselines/prev.json
```

### Key matrices to watch

| Matrix | Pre-optimization | Target |
|--------|-----------------|--------|
| c-71 | 5.48x SPRAL | < 2x SPRAL |
| c-big | 6.29x SPRAL | < 2x SPRAL |
| bloweybq | ~2.0x SPRAL | < 1.5x SPRAL |
| dixmaanl | 2.51x SPRAL | < 1.5x SPRAL |

## Implementation Order

1. **FactorizationWorkspace struct** — Define the workspace type, constructor from AptpSymbolic
2. **Hoist frontal matrix allocation** — Modify `factor_single_supernode` to accept `&mut FactorizationWorkspace` and use workspace buffer instead of `Mat::zeros`
3. **Hoist APTP kernel temporaries** — Pre-allocate backup, l11_copy, ld workspace buffers; pass to `factor_inner` / `two_level_factor`
4. **Thread-local workspace for parallel path** — Add `Cell<FactorizationWorkspace>` thread-local alongside existing `G2L_BUF`
5. **Contribution block copy optimization** — Reduce or eliminate the element-by-element copy in `extract_contribution`
6. **Benchmark and validate** — Run full SuiteSparse suite, compare against baselines

## Files to Modify

| File | Changes |
|------|---------|
| `src/aptp/numeric.rs` | Add FactorizationWorkspace; modify `factor_single_supernode`, `factor_tree_levelset` |
| `src/aptp/factor.rs` | Add AptpKernelWorkspace; modify `factor_inner`, `two_level_factor`, `apply_and_check`, `update_trailing`, `update_cross_terms` |
| `src/aptp/solver.rs` | No public API changes; workspace allocated internally in `SparseLDLT::factor` |

## No New Dependencies

This feature uses only existing dependencies (faer, rayon). No new crates needed.
