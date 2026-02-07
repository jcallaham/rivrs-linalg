# Quickstart: Continuous Integration Setup

**Feature**: 007-ci-setup
**Date**: 2026-02-07

## What This Feature Does

Enhances the existing GitHub Actions CI pipeline (`.github/workflows/ci.yml`) to add benchmark compilation verification for the sparse domain. The existing test, lint, and doc jobs already satisfy most spec requirements.

## Current State

The CI pipeline already provides:
- Test execution on MSRV (1.87) and stable Rust
- Clippy linting with `-D warnings`
- Formatting checks via `cargo fmt --check`
- Rustdoc build with `-D warnings`
- Build caching via Swatinem/rust-cache@v2
- Independent sparse/control domain jobs
- test-util feature activation via self-referencing dev-dependency

## What's Being Added

1. **bench-sparse job**: Compiles benchmark binary (`solver_benchmarks`) without execution, catching compilation regressions in the benchmarking infrastructure

## How to Verify

After implementation, the CI should show four sparse domain jobs on every PR:

```
test-sparse (stable)   ✓
test-sparse (1.87)     ✓
lint-sparse            ✓
doc-sparse             ✓
bench-sparse           ✓
```

### Local verification

```bash
# Replicate CI jobs locally:
cd sparse

# Test (matches test-sparse job)
cargo test --all-targets

# Lint (matches lint-sparse job)
cargo fmt --check
cargo clippy --all-targets -- -D warnings

# Doc (matches doc-sparse job)
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps

# Bench compilation (matches bench-sparse job)
cargo bench --no-run
```

## Files Modified

- `.github/workflows/ci.yml` — Add bench-sparse job

## Dependencies

- No new GitHub Actions dependencies
- Uses existing `dtolnay/rust-toolchain`, `Swatinem/rust-cache@v2`, `actions/checkout@v4`
