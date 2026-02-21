# Quickstart: Sequential Profiling & Optimization (Phase 8.1g)

**Date**: 2026-02-20

## Prerequisites

- Rust 1.87+ (edition 2024)
- Full SuiteSparse test data extracted to `test-data/suitesparse/` (65 matrices)
- Branch: `018-sequential-profiling-optimization`

## Development Workflow

### Step 1: Verify baseline correctness (before any changes)

```bash
cd /workspace/rivrs-linalg/sparse

# Run all unit tests
cargo test

# Run full SuiteSparse validation (65 matrices)
cargo test -- --ignored --test-threads=1
```

All tests must pass before starting optimization work.

### Step 2: Add profiling instrumentation

Instrument `AptpNumeric::factor()` in `src/aptp/numeric.rs` with `ProfileSession`/`SectionGuard` calls, gated behind `cfg(feature = "diagnostic")`.

```bash
# Build with diagnostic feature
cargo build --features diagnostic

# Run tests with diagnostic feature to verify no breakage
cargo test --features diagnostic
```

### Step 3: Optimize allocation hotspots

Focus on `factor_inner` in `src/aptp/factor.rs`:
1. Hoist panel row permutation buffer (line 1664)
2. Hoist row permutation temp buffer (line 1692)
3. Hoist column order slice buffer (line 1704)

After each optimization:
```bash
# Verify correctness
cargo test
cargo test -- --ignored --test-threads=1

# Verify no clippy warnings
cargo clippy --all-targets
```

### Step 4: Collect baselines

```bash
# CI subset (quick, ~10 matrices)
cargo run --example baseline_collection --features diagnostic --release -- --ci-only

# Full suite (all 65 matrices)
cargo run --example baseline_collection --features diagnostic --release
```

### Step 5: Generate workload analysis

```bash
cargo run --example workload_analysis --features diagnostic --release
```

### Step 6: Produce performance report

The workload analysis output, combined with baseline data, forms the performance report documented in `docs/phase-8.1g-report.md`.

## Key Files

| File | Role |
|------|------|
| `src/aptp/numeric.rs` | Profiling instrumentation in factor loop |
| `src/aptp/factor.rs` | Allocation optimization in `factor_inner` |
| `src/profiling/session.rs` | Existing ProfileSession/SectionGuard API |
| `src/profiling/memory.rs` | Existing MemoryTracker API |
| `examples/baseline_collection.rs` | New: systematic baseline collection tool |
| `examples/workload_analysis.rs` | New: workload distribution analysis tool |
| `docs/phase-8.1g-report.md` | New: performance analysis report |

## Verification Checklist

- [ ] All 362+ unit tests pass
- [ ] All 65 SuiteSparse `--ignored` tests pass
- [ ] `cargo clippy --all-targets` clean
- [ ] `cargo clippy --all-targets --features diagnostic` clean
- [ ] Baseline JSON produced for full SuiteSparse suite
- [ ] Workload analysis report produced
- [ ] Performance report written with Phase 8.2 recommendation
