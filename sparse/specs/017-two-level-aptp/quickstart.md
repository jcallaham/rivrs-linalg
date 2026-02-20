# Development Quickstart: Two-Level APTP Factorization

**Feature**: 017-two-level-aptp
**Branch**: `017-two-level-aptp`
**Date**: 2026-02-16

## Prerequisites

- Rust 1.87+ (edition 2024)
- Branch `017-two-level-aptp` checked out from `ssids`
- Test data available: `test-data/hand-constructed/` (git-tracked), `test-data/suitesparse-ci/` (git-tracked)

## Build & Test

```bash
cd /workspace/rivrs-linalg/sparse

# Build
cargo build

# Run unit tests (fast, no ignored tests)
cargo test

# Run full SuiteSparse tests (slow, requires test-data/suitesparse/)
cargo test -- --ignored --test-threads=1

# Run benchmarks
cargo bench

# Lint
cargo clippy
cargo fmt --check
```

## Key Files to Modify

| File | Change | Priority |
|------|--------|----------|
| `src/aptp/factor.rs` | Primary: complete pivoting, two-level kernel, backup, BLAS-3 ops | P0 |
| `src/aptp/factor.rs` (tests) | Add complete pivoting unit tests, two-level tests, block boundary tests | P0 |
| `src/aptp/solver.rs` | Add block size fields to `FactorOptions` | P1 |
| `src/aptp/numeric.rs` | Pass block size options through to kernel | P1 |
| `tests/multifrontal.rs` | Add block-boundary integration tests | P1 |
| `tests/solve.rs` | Verify no regression on end-to-end solve | P1 |
| `benches/solver_benchmarks.rs` | Add two-level vs single-level kernel benchmarks | P2 |

## Implementation Order

### Step 1: Complete Pivoting (Algorithm 4.1)
Write `complete_pivoting_factor` function and its unit tests. This is an independent function with no dependencies on the two-level infrastructure. Test against known matrices where complete pivoting order is analytically predictable.

### Step 2: factor_inner (Inner APTP with Complete Pivoting Leaves)
Refactor the existing `aptp_factor_in_place` main loop to call `complete_pivoting_factor` at the ib×ib diagonal level. At this point, the single-level kernel still works but uses complete pivoting for small diagonal blocks. All existing tests must pass.

### Step 3: BlockBackup
Implement the per-block backup and restore mechanism. Test with synthetic matrices that trigger backup/restore.

### Step 4: apply_and_check (TRSM + Threshold Check)
Implement the Apply phase: TRSM-based L21 computation with a posteriori threshold checking. Test in isolation with known L11/D11/A21 inputs.

### Step 5: update_trailing (GEMM)
Implement the Update phase: GEMM-based Schur complement. Test in isolation.

### Step 6: Two-Level Integration
Wire the outer block loop: backup → factor_inner → apply_and_check → restore_failed → update_trailing → update_delayed. Test on dense matrices of sizes 64, 256, 257, 512, 1024. Verify reconstruction < 1e-12.

### Step 7: Replace Single-Level
Replace the existing single-level main loop with the two-level dispatch. For p ≤ nb, the loop executes once (single block). Verify ALL existing tests pass.

### Step 8: Integration & Benchmarks
Run full SuiteSparse suite. Benchmark. Document crossover point.

## Key References

- **Algorithm 4.1** (complete pivoting): Duff et al. 2020, Section 4 → `/workspace/rivrs-linalg/references/ssids/duff2020.md`
- **Algorithm 3.1** (two-level APTP): Duff et al. 2020, Section 3 → same file
- **SPRAL kernel**: `/opt/references/spral/src/ssids/cpu/kernels/ldlt_app.cxx`
- **Current kernel**: `src/aptp/factor.rs` (lines 197-369: main loop, lines 460-654: helper functions)
- **Current tests**: `src/aptp/factor.rs` (lines 696-1629: embedded test module)

## Verification Commands

```bash
# After Step 1 (complete pivoting):
cargo test complete_pivoting

# After Step 2 (inner APTP refactor):
cargo test aptp_factor  # all existing factor tests

# After Step 6 (two-level integration):
cargo test two_level

# After Step 7 (replacement):
cargo test              # ALL tests must pass
cargo test -- --ignored --test-threads=1  # full SuiteSparse

# Benchmarks:
cargo bench -- kernel   # kernel-level benchmarks
cargo bench -- solver   # full pipeline benchmarks
```

## Debugging Tips

- Use `AptpOptions { outer_block_size: usize::MAX, ..Default::default() }` to force single-level behavior for comparison testing.
- Use `AptpStatistics` output (num_1x1, num_2x2, num_delayed) to compare pivot decisions between single-level and two-level.
- Use `AptpPivotRecord` per-column log to trace exactly where pivoting diverges.
- The profiling module (`src/profiling/`) can instrument Factor/Apply/Update phases to identify bottlenecks.
