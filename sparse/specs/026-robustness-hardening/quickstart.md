# Quickstart: Robustness — Testing & Hardening

**Feature**: 026-robustness-hardening
**Date**: 2026-02-25

## Overview

This feature adds four categories of testing to the rivrs-sparse SSIDS solver:

1. **SPRAL test parity audit** — document mapping + test suite pruning
2. **Torture tests** — 500+ random instances with adversarial perturbations at the dense kernel level
3. **Property-based tests** — proptest-driven structural invariant verification
4. **Adversarial input tests** — edge cases and malformed inputs

## Prerequisites

- Rust 1.87+ (edition 2024)
- Working `cargo test` on current codebase (524 tests: 510 pass, 14 ignored)
- SPRAL source available at `/opt/references/spral/` (for audit)

## New Dependencies

```toml
# Cargo.toml [dev-dependencies]
proptest = "1.4"
```

## New Source Files

```
src/testing/
├── perturbations.rs   # cause_delays, make_singular, make_dblk_singular
└── strategies.rs      # proptest Strategy impls for matrix generation

docs/
└── spral-test-audit.md  # SPRAL test parity audit document
```

## Running New Tests

```bash
# Normal test suite (includes property tests, adversarial tests)
cargo test

# Torture tests (long-running, #[ignore])
cargo test -- --ignored torture --test-threads=1

# Property tests only
cargo test property

# Adversarial input tests only
cargo test adversarial
```

## Implementation Order

1. **Audit** (P1): Read SPRAL tests → write audit document → prune redundant tests
2. **Perturbation helpers** (P1): Implement `cause_delays`, `make_singular`, `make_dblk_singular`
3. **Torture tests** (P1): Wire perturbation helpers into 500+ instance test suites
4. **Proptest strategies** (P2): Create matrix generation strategies
5. **Property tests** (P2): Implement structural invariant properties
6. **Adversarial tests** (P2): Add edge case and malformed input tests
7. **Solver hardening** (P2): Fix any panics discovered by adversarial tests

## Validation

```bash
# Full validation after implementation
cargo test                                    # All non-ignored tests pass
cargo test -- --ignored --test-threads=1      # All ignored tests pass (SuiteSparse + torture)
cargo clippy --all-targets                    # No warnings
cargo clippy --all-targets --features diagnostic  # No warnings with diagnostic
```
