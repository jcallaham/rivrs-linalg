# Tasks: Module Restructure (027)

## Phase 1: File Moves & Module Creation

- [x] T1: Create `src/ordering/` directory and `mod.rs` with module declarations and re-exports
- [x] T2: Move `src/aptp/ordering.rs` → `src/ordering/metis.rs` and update internal imports
- [x] T3: Move `src/aptp/matching.rs` → `src/ordering/matching.rs` and update internal imports
- [x] T4: Move `src/aptp/perm.rs` → `src/ordering/perm.rs` and update internal imports

## Phase 2: Rename aptp → symmetric

- [x] T5: Rename `src/aptp/` → `src/symmetric/`, remove moved files, update `lib.rs` module declarations
- [x] T6: Update `src/symmetric/mod.rs` — remove ordering/matching/perm submodules, adjust re-exports
- [x] T7: Update `src/symmetric/solver.rs` — change ordering imports to `crate::ordering::`
- [x] T8: Update `src/symmetric/factor.rs` — change perm import to `crate::ordering::`
- [x] T9: Update any other internal `super::` or `crate::aptp::` references in symmetric/ files

## Phase 3: Update Crate-Internal Consumers

- [x] T10: Update `src/io/reference.rs` — `crate::aptp::Inertia` → `crate::symmetric::Inertia`
- [x] T11: Update `src/testing/mc64_validation.rs` — imports from `crate::aptp` → `crate::ordering`
- [x] T12: Update `src/lib.rs` — module declarations, any top-level re-exports

## Phase 4: Update Integration Tests

- [x] T13: Update `tests/aptp_data_structures.rs` [P]
- [x] T14: Update `tests/mc64_matching.rs` [P]
- [x] T15: Update `tests/mc64_profiling.rs` [P]
- [x] T16: Update `tests/match_order.rs` [P]
- [x] T17: Update `tests/metis_ordering.rs` [P]
- [x] T18: Update `tests/multifrontal.rs` [P]
- [x] T19: Update `tests/symbolic_analysis.rs` and `tests/symbolic_analysis_full.rs` [P]
- [x] T20: Update `tests/solve.rs` [P]
- [x] T21: Update `tests/hand_constructed.rs` [P]
- [x] T22: Update `tests/suitesparse_ci.rs` [P]
- [x] T23: Update `tests/property.rs` and `tests/adversarial.rs` [P]

## Phase 5: Update Examples, Benchmarks, Comparisons

- [x] T24: Update all examples (`basic_usage.rs`, `solve_timing.rs`, `profile_matrix.rs`, `parallel_scaling.rs`, `baseline_collection.rs`) [P]
- [x] T25: Update `benches/solver_benchmarks.rs` [P]
- [x] T26: Update `comparisons/src/spral_benchmark.rs` [P]

## Phase 6: Verify & Fix

- [x] T27: Run `cargo build` and fix any remaining import errors
- [x] T28: Run `cargo test` — all tests must pass (525 passed, 0 failed, 23 ignored)
- [x] T29: Run `cargo clippy --all-targets --features diagnostic` — zero warnings
- [x] T30: Run `cargo fmt` and verify clean formatting
- [x] T31: Run `cargo doc --no-deps` — docs build successfully

## Phase 7: Update Documentation

- [x] T32: Update `sparse/CLAUDE.md` Source Code Layout and module references
- [x] T33: Add restructure entry to `docs/ssids-log.md`
- [x] T34: Update any `aptp/` path references in other docs files
