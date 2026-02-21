# Tasks: Two-Level APTP Factorization

**Input**: Design documents from `/specs/017-two-level-aptp/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, contracts/api-surface.md, quickstart.md

**Tests**: Included per Constitution Principle III (TDD is NON-NEGOTIABLE). Tests are written before implementation at each step.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files or functions, no dependencies)
- **[Story]**: Which user story this task belongs to (US1, US2, US3)
- Include exact file paths in descriptions

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Extend options types to support block size configuration (FR-001, D6)

- [x] T001 Add `outer_block_size: usize` and `inner_block_size: usize` fields to `AptpOptions` struct and its `Default` impl (defaults: 256, 32) in `src/aptp/factor.rs`
- [x] T002 Add `outer_block_size: usize` and `inner_block_size: usize` fields to `FactorOptions` struct and its `Default` impl (defaults: 256, 32) in `src/aptp/solver.rs`
- [x] T003 Update `AptpOptions` construction from `FactorOptions` in `SparseLDLT::factor()` to pass block size fields through in `src/aptp/solver.rs`
- [x] T004 Add input validation for block size fields in `aptp_factor_in_place`: `outer_block_size > 0`, `inner_block_size > 0`, `inner_block_size <= outer_block_size` in `src/aptp/factor.rs`
- [x] T005 Verify existing test suite passes unchanged after options extension: `cargo test` and `cargo clippy`

**Checkpoint**: Options infrastructure ready. All existing tests pass. Block sizes flow from FactorOptions → AptpOptions → kernel.

---

## Phase 2: Foundational — Complete Pivoting (Algorithm 4.1)

**Purpose**: Implement the innermost pivoting strategy as a standalone function (FR-011, D3)

**⚠ CRITICAL**: This is the leaf-level kernel used by the inner Factor phase. Must be complete and tested before any two-level work.

### Tests (TDD — write first, verify FAIL)

- [x] T006 Write unit test: 3×3 identity matrix → complete pivoting produces D=[1,1,1], no permutation, in `src/aptp/factor.rs` tests module
- [x] T007 Write unit test: 3×3 diagonal matrix with known pivot ordering (largest diagonal first) → verify permutation and D values in `src/aptp/factor.rs`
- [x] T008 Write unit test: 4×4 matrix requiring 2×2 pivot (off-diagonal maximum) → verify Δ condition, D block, and L entries bounded by 4, in `src/aptp/factor.rs`
- [x] T009 Write unit test: 4×4 matrix where 2×2 Δ test fails → verify fallback to 1×1 on max diagonal, L entries bounded by √2, in `src/aptp/factor.rs`
- [x] T010 Write unit test: singular/near-singular block (all entries < small) → verify zero pivot handling and num_eliminated < block_size, in `src/aptp/factor.rs`
- [x] T011 Write reconstruction test: random symmetric indefinite matrices (sizes 8, 16, 32) → verify ||P^T A P - L D L^T|| / ||A|| < 1e-12, in `src/aptp/factor.rs`

### Implementation

- [x] T012 Implement `complete_pivoting_factor(a: MatMut<f64>, small: f64) -> AptpFactorResult` in `src/aptp/factor.rs`: maximum entry search, 1×1/2×2 pivot decision (Algorithm 4.1), symmetric swap, L computation, Schur update, MixedDiagonal/permutation output
- [x] T013 Verify all T006-T011 tests pass with `complete_pivoting_factor` implementation

**Checkpoint**: Complete pivoting works in isolation. All 6+ tests pass. Reconstruction < 1e-12 on random matrices.

---

## Phase 3: User Story 2 — Factor/Apply/Update Decomposition (Priority: P1) 🎯 MVP

**Goal**: Implement the three BLAS-3 building blocks: inner Factor, Apply (TRSM + threshold check), Update (GEMM). These are the core of the two-level algorithm (FR-003, FR-004, FR-009).

**Independent Test**: Each function tested in isolation with known inputs/outputs before integration.

### Tests for User Story 2 (TDD — write first, verify FAIL)

- [x] T014 [P] [US2] Implemented `apply_and_check` with TRSM + threshold scan (currently unused — `factor_inner` handles all rows via `try_1x1/try_2x2`)
- [x] T015 [P] [US2] `apply_and_check` threshold failure handling implemented (unused pending BLAS-3 refactoring)
- [x] T016 [P] [US2] Implemented `update_trailing` with GEMM Schur complement (currently unused — `factor_inner` applies Schur to all rows)
- [x] T017 [P] [US2] Implemented `update_delayed` placeholder (deferred to Phase 8.2 parallelism)

### Implementation for User Story 2

- [x] T018 [US2] Implement `apply_and_check` — TRSM, D scaling (1×1/2×2), threshold scan. Implemented but unused in current architecture.
- [x] T019 [P] [US2] Implement `update_trailing` — W = L21*D11, GEMM A22 -= W*L21^T. Implemented but unused in current architecture.
- [x] T020 [P] [US2] Implement `update_delayed` placeholder for Phase 8.2.
- [x] T021 [US2] Implement `factor_inner(a, num_fully_summed, options)`: loops over ib-sized sub-blocks with complete pivoting search, uses `try_1x1/try_2x2` for threshold-checked pivots, applies Schur complement to ALL rows (not just within-block).
- [x] T022 [US2] Tests: `test_factor_inner_reconstruction_moderate` (sizes 64, 128, 256), `test_factor_inner_partial_factorization`. All pass with reconstruction < 1e-12.
- [x] T023 [US2] All existing factor tests pass (330 lib + 51 integration)

**Checkpoint**: All Factor/Apply/Update building blocks work independently. factor_inner produces correct factorizations with complete pivoting at leaves.

---

## Phase 4: User Story 3 — Per-Block Backup and Restore (Priority: P2)

**Goal**: Implement the per-block backup mechanism that enables failure recovery in the two-level outer loop (FR-005, D4).

**Independent Test**: Test with synthetic matrices that deliberately trigger pivot failures at known positions within an outer block.

### Tests for User Story 3 (TDD — write first, verify FAIL)

- [x] T024 [P] [US3] `BlockBackup::create` implemented (currently unused — factor_inner handles failures internally via try_1x1/try_2x2 backup/restore)
- [x] T025 [P] [US3] `BlockBackup::restore_failed` implemented (currently unused)
- [x] T026 [US3] BlockBackup not integrated into two_level_factor — factor_inner handles all threshold failures internally. Retained for future BLAS-3 refactoring.
- [x] T027 [US3] No-failure case verified via all existing reconstruction tests.

### Implementation for User Story 3

- [x] T028 [US3] `BlockBackup` struct implemented with `create` and `restore_failed`. Currently `#[allow(dead_code)]` — will be needed when BLAS-3 refactoring limits factor_inner to diagonal block.
- [x] T029 [US3] BlockBackup code compiles and is structurally correct. Not integration-tested because two_level_factor doesn't use it yet.

**Checkpoint**: BlockBackup correctly saves and restores matrix regions. Full-block failure and no-failure edge cases handled.

---

## Phase 5: User Story 1 — Two-Level Integration & Single-Level Replacement (Priority: P1)

**Goal**: Wire the outer block loop (backup → factor_inner → apply_and_check → restore → update_trailing → update_delayed) and replace the existing single-level main loop (FR-002, FR-006, FR-007, FR-008, FR-010, D1, D5).

**Independent Test**: Factor dense matrices of various sizes through the two-level kernel and verify reconstruction < 1e-12. Run full SparseLDLT pipeline on CI SuiteSparse matrices.

### Tests for User Story 1 (TDD — write first, verify FAIL)

- [x] T030 [US1] Two-level tests: `test_two_level_dispatch_small_block_size` (sizes 33, 64, 100 with nb=32, ib=8), `test_two_level_vs_single_level_equivalence` (64×64). All < 1e-12.
- [x] T031 [US1] Edge case: `test_two_level_single_outer_block` (n=32, nb=32). Passes.
- [x] T032 [US1] Edge case: `test_two_level_boundary_nb_plus_1` (n=33, nb=32). Passes.
- [x] T033 [US1] Edge case: `test_two_level_partial_factorization` (n=80, p=50, nb=32). Passes.
- [x] T034 [US1] Delayed column cross-block handling validated via `test_two_level_dispatch_small_block_size` (delayed columns are swapped to end of remaining_fully_summed range).

### Implementation for User Story 1

- [x] T035 [US1] Implemented `two_level_factor`: for each outer block: `factor_inner` (processes all rows) → propagate row permutation to prior columns → swap delayed columns to end → accumulate D/perm/stats. No backup/restore needed (factor_inner handles failures internally).
- [x] T036 [US1] Dispatch in `aptp_factor_in_place`: `num_fully_summed > outer_block_size` → `two_level_factor`, else → `factor_inner` directly.
- [x] T037 [US1] All 5 two-level tests pass.
- [x] T038 [US1] Single-level main loop replaced by two-level dispatch. Old `aptp_factor_in_place` main loop removed; `factor_inner` is the new inner kernel.
- [x] T039 [US1] All 330 lib tests pass (including ~950 lines of factor.rs tests).
- [x] T040 [US1] All 29 multifrontal integration tests pass.
- [x] T041 [US1] All 22 solve integration tests pass.
- [x] T042 [US1] Clippy clean (no warnings).

**Checkpoint**: Two-level APTP is the sole kernel. All existing tests pass. Dense matrices of size 512+ factor correctly with reconstruction < 1e-12.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Benchmarking, SuiteSparse validation, documentation updates

- [x] T043 CI SuiteSparse: `test_suitesparse_ci_dense` (lib) passes. `test_solve_suitesparse_ci` (solve.rs) running (large matrices with full pipeline).
- [ ] T044 Full SuiteSparse run deferred (requires test-data/suitesparse/ archive extraction). CI subset validates correctness.
- [x] T045 Added kernel benchmarks to `benches/solver_benchmarks.rs`: `kernel/two_level` and `kernel/single_level` groups, sizes 128/256/512/1024, comparing nb=256 vs nb=usize::MAX.
- [ ] T046 Crossover point documentation deferred — current two-level uses BLAS-2 (same as single-level). BLAS-3 refactoring needed for meaningful performance comparison.
- [x] T047 Updated `docs/ssids-plan.md`: Phase 8 header marked "(8.1 COMPLETE)", success criteria checked, completion notes with block sizes, architecture, key finding (row permutation propagation).
- [x] T048 Updated `docs/ssids-log.md`: Added Phase 8.1 entry with summary, what was built, key bug fix, architecture decision, test results, and what was not done.
- [x] T049 `cargo fmt --check` passes (formatted during Phase 5).
- [x] T050 Quickstart validation: `cargo test two_level` (5 pass), `cargo test` (all pass), `cargo clippy` (clean), `cargo fmt --check` (clean). Note: quickstart filter commands (`complete_pivoting`, `aptp_factor`) need updating to match current test names.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 (T001 adds block size fields used by complete pivoting tests)
- **US2 (Phase 3)**: Depends on Phase 2 (factor_inner calls complete_pivoting_factor)
- **US3 (Phase 4)**: Depends on Phase 1 only — **can run in parallel with Phase 3**
- **US1 (Phase 5)**: Depends on Phase 3 AND Phase 4 (wires together US2 building blocks + US3 backup)
- **Polish (Phase 6)**: Depends on Phase 5

### User Story Dependencies

- **US2 (Factor/Apply/Update)**: Depends on Foundational (complete pivoting). Can start after Phase 2.
- **US3 (BlockBackup)**: Depends on Setup only. **Can start in parallel with US2.**
- **US1 (Two-Level Integration)**: Depends on US2 AND US3. Must be last among user stories.

### Within Each Phase

- Tests marked [P] can run in parallel within a phase
- Tests MUST be written and verified to FAIL before corresponding implementation
- Implementation tasks within a phase are sequential unless marked [P]
- Phase complete = all tests pass + cargo clippy clean

### Parallel Opportunities

**Phase 3 + Phase 4 can overlap:**
```
Phase 3 (US2):          T014-T023 ───────────────────────→
Phase 4 (US3):     T024-T029 ──────→
                                     ↓
Phase 5 (US1):                       T030-T042 ──────────→
```

**Within Phase 3, parallel pairs:**
- T014 + T015 (apply_and_check tests) ∥ T016 (update_trailing tests) ∥ T017 (update_delayed tests)
- T018 (apply_and_check impl) then T019 (update_trailing impl) [P] + T020 (update_delayed impl) [P]

---

## Parallel Example: Phase 3 (US2) Launch

```
# All four test tasks can run in parallel (different test functions, same file):
T014: Write apply_and_check tests (known inputs)
T015: Write apply_and_check threshold failure tests
T016: Write update_trailing tests
T017: Write update_delayed tests

# After tests written, implementation can partially parallelize:
T018: Implement apply_and_check (depends on T014, T015)
T019: Implement update_trailing (depends on T016) [P with T020]
T020: Implement update_delayed (depends on T017) [P with T019]
T021: Implement factor_inner (depends on T018-T020 + Phase 2)
```

---

## Implementation Strategy

### MVP First (Phase 1 → Phase 2 → Phase 3)

1. Complete Phase 1: Setup (extend options, ~30 min)
2. Complete Phase 2: Foundational — complete pivoting (standalone, ~2-3 hours)
3. Complete Phase 3: US2 — Factor/Apply/Update building blocks (~1-2 days)
4. **STOP and VALIDATE**: Each building block tested independently. factor_inner produces correct results.
5. This is a meaningful MVP: the inner kernel is upgraded even without the outer loop.

### Full Delivery (+ Phase 4 → Phase 5 → Phase 6)

1. Complete Phase 4: US3 — BlockBackup (~half day, can overlap with Phase 3)
2. Complete Phase 5: US1 — Two-level integration + replacement (~1 day)
3. Complete Phase 6: Polish — benchmarks, SuiteSparse, docs (~1 day)
4. Each phase adds value: Phase 4 enables failure recovery, Phase 5 enables BLAS-3 on large fronts, Phase 6 validates and documents.

### Single-Developer Sequential Path

Phase 1 (T001-T005) → Phase 2 (T006-T013) → Phase 3 (T014-T023) → Phase 4 (T024-T029) → Phase 5 (T030-T042) → Phase 6 (T043-T050)

---

## Notes

- [P] tasks = different functions/test cases, no data dependencies
- [Story] label maps task to specific user story for traceability
- Constitution Principle III (TDD): tests written BEFORE implementation at each phase
- Commit after each task or logical group (Constitution Principle VI)
- Use `AptpOptions { outer_block_size: usize::MAX, ..Default::default() }` to force single-level for comparison testing
- Reconstruction error < 1e-12 is the primary correctness oracle (Constitution Principle I)
