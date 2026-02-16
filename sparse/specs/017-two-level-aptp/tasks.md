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

- [ ] T006 Write unit test: 3×3 identity matrix → complete pivoting produces D=[1,1,1], no permutation, in `src/aptp/factor.rs` tests module
- [ ] T007 Write unit test: 3×3 diagonal matrix with known pivot ordering (largest diagonal first) → verify permutation and D values in `src/aptp/factor.rs`
- [ ] T008 Write unit test: 4×4 matrix requiring 2×2 pivot (off-diagonal maximum) → verify Δ condition, D block, and L entries bounded by 4, in `src/aptp/factor.rs`
- [ ] T009 Write unit test: 4×4 matrix where 2×2 Δ test fails → verify fallback to 1×1 on max diagonal, L entries bounded by √2, in `src/aptp/factor.rs`
- [ ] T010 Write unit test: singular/near-singular block (all entries < small) → verify zero pivot handling and num_eliminated < block_size, in `src/aptp/factor.rs`
- [ ] T011 Write reconstruction test: random symmetric indefinite matrices (sizes 8, 16, 32) → verify ||P^T A P - L D L^T|| / ||A|| < 1e-12, in `src/aptp/factor.rs`

### Implementation

- [ ] T012 Implement `complete_pivoting_factor(a: MatMut<f64>, small: f64) -> AptpFactorResult` in `src/aptp/factor.rs`: maximum entry search, 1×1/2×2 pivot decision (Algorithm 4.1), symmetric swap, L computation, Schur update, MixedDiagonal/permutation output
- [ ] T013 Verify all T006-T011 tests pass with `complete_pivoting_factor` implementation

**Checkpoint**: Complete pivoting works in isolation. All 6+ tests pass. Reconstruction < 1e-12 on random matrices.

---

## Phase 3: User Story 2 — Factor/Apply/Update Decomposition (Priority: P1) 🎯 MVP

**Goal**: Implement the three BLAS-3 building blocks: inner Factor, Apply (TRSM + threshold check), Update (GEMM). These are the core of the two-level algorithm (FR-003, FR-004, FR-009).

**Independent Test**: Each function tested in isolation with known inputs/outputs before integration.

### Tests for User Story 2 (TDD — write first, verify FAIL)

- [ ] T014 [P] [US2] Write unit tests for `apply_and_check`: (a) given a pre-factored 32×32 L11/D11 with all 1×1 pivots and a known 64×32 A21 panel, verify L21 = A21 * (L11 * D11)^{-T} matches expected values and threshold check returns correct nelim; (b) given L11/D11 containing mixed 1×1 and 2×2 pivot blocks in D11, verify D scaling correctly handles both pivot sizes (FR-009), in `src/aptp/factor.rs`
- [ ] T015 [P] [US2] Write unit tests for `apply_and_check` with threshold failure: given L21 entries that exceed 1/u, verify nelim is reduced to first failing column, in `src/aptp/factor.rs`
- [ ] T016 [P] [US2] Write unit tests for `update_trailing`: given a known L21 panel and D11, verify A22 -= L21 * D11 * L21^T matches direct computation, in `src/aptp/factor.rs`
- [ ] T017 [P] [US2] Write unit tests for `update_delayed`: given a factored block and a delayed column region, verify the rank-nelim update is applied correctly, in `src/aptp/factor.rs`

### Implementation for User Story 2

- [ ] T018 [US2] Implement `apply_and_check(a, col_start, block_nelim, block_cols, m, d, threshold) -> usize` in `src/aptp/factor.rs`: TRSM via faer's triangular solve, D scaling accounting for 1×1 and 2×2 pivots, column-by-column threshold scan
- [ ] T019 [P] [US2] Implement `update_trailing(a, col_start, nelim, m, d)` in `src/aptp/factor.rs`: compute W = L21 * D11 in workspace, then GEMM A22 -= W * L21^T via faer's matmul
- [ ] T020 [P] [US2] Implement `update_delayed(a, col_start, nelim, delayed_ranges, d)` in `src/aptp/factor.rs`: iterate over delayed regions, apply rank-nelim update using L and D from current block
- [ ] T021 [US2] Implement `factor_inner(a, num_fully_summed, options) -> Result<AptpFactorResult>` in `src/aptp/factor.rs`: loop over ib-sized sub-blocks; for each sub-block call `complete_pivoting_factor` on the ib×ib diagonal (processes all ib columns at once), then inner Apply (TRSM to update panel below diagonal within nb block) and inner Update (GEMM Schur complement within nb block)
- [ ] T022 [US2] Write tests for `factor_inner` on nb×nb matrices (128, 256) verifying reconstruction < 1e-12 and comparing statistics (num_1x1, num_2x2, num_delayed) against single-level baseline, in `src/aptp/factor.rs`
- [ ] T023 [US2] Verify T014-T017 tests pass and all existing factor tests still pass

**Checkpoint**: All Factor/Apply/Update building blocks work independently. factor_inner produces correct factorizations with complete pivoting at leaves.

---

## Phase 4: User Story 3 — Per-Block Backup and Restore (Priority: P2)

**Goal**: Implement the per-block backup mechanism that enables failure recovery in the two-level outer loop (FR-005, D4).

**Independent Test**: Test with synthetic matrices that deliberately trigger pivot failures at known positions within an outer block.

### Tests for User Story 3 (TDD — write first, verify FAIL)

- [ ] T024 [P] [US3] Write unit test for `BlockBackup::create`: back up a known 512×256 column of blocks, verify data is correctly copied, in `src/aptp/factor.rs`
- [ ] T025 [P] [US3] Write unit test for `BlockBackup::restore_failed`: given nelim=5 out of 32, verify columns 5-31 are restored from backup while columns 0-4 are untouched, in `src/aptp/factor.rs`
- [ ] T026 [US3] Write test for full-block failure: all columns in an outer block fail → entire block restored from backup, delayed columns list contains all block columns, in `src/aptp/factor.rs`
- [ ] T027 [US3] Write test for no-failure case: backup created but never restored, results satisfy reconstruction < 1e-12, in `src/aptp/factor.rs`

### Implementation for User Story 3

- [ ] T028 [US3] Implement `BlockBackup` struct with `create(a, col_start, block_cols, m) -> Self` and `restore_failed(a, col_start, nelim, block_cols, m)` in `src/aptp/factor.rs`
- [ ] T029 [US3] Verify T024-T027 tests pass

**Checkpoint**: BlockBackup correctly saves and restores matrix regions. Full-block failure and no-failure edge cases handled.

---

## Phase 5: User Story 1 — Two-Level Integration & Single-Level Replacement (Priority: P1)

**Goal**: Wire the outer block loop (backup → factor_inner → apply_and_check → restore → update_trailing → update_delayed) and replace the existing single-level main loop (FR-002, FR-006, FR-007, FR-008, FR-010, D1, D5).

**Independent Test**: Factor dense matrices of various sizes through the two-level kernel and verify reconstruction < 1e-12. Run full SparseLDLT pipeline on CI SuiteSparse matrices.

### Tests for User Story 1 (TDD — write first, verify FAIL)

- [ ] T030 [US1] Write two-level kernel tests for dense symmetric indefinite matrices at sizes 257, 512, 1024 verifying reconstruction < 1e-12 with default block sizes (nb=256, ib=32), in `src/aptp/factor.rs`
- [ ] T031 [US1] Write edge case test: frontal matrix dimension == nb (256) → single outer block, equivalent to factor_inner, in `src/aptp/factor.rs`
- [ ] T032 [US1] Write edge case test: frontal matrix dimension == nb+1 (257) → two blocks, second block has 1 column, in `src/aptp/factor.rs`
- [ ] T033 [US1] Write edge case test: frontal matrix with partial factorization (num_fully_summed < m) and dimension > nb, verifying contribution block (Schur complement) is correctly updated, in `src/aptp/factor.rs`
- [ ] T034 [US1] Write test with matrices that trigger delayed columns across outer block boundaries → verify UpdateNT/UpdateTN correctly updates delayed columns from earlier blocks, in `src/aptp/factor.rs`

### Implementation for User Story 1

- [ ] T035 [US1] Implement the two-level outer block loop in a new function (e.g. `two_level_factor`) in `src/aptp/factor.rs`: for each outer block: BlockBackup::create → factor_inner → apply_and_check → backup.restore_failed (if needed) → update_trailing → update_delayed; accumulate into global AptpFactorResult
- [ ] T036 [US1] Add dispatch logic to `aptp_factor_in_place`: if `num_fully_summed > options.outer_block_size` call `two_level_factor`, else call `factor_inner` directly, in `src/aptp/factor.rs`
- [ ] T037 [US1] Verify T030-T034 two-level tests pass
- [ ] T038 [US1] Replace existing single-level main loop in `aptp_factor_in_place` with the two-level dispatch as the sole code path in `src/aptp/factor.rs`
- [ ] T039 [US1] Verify ALL existing embedded tests in `src/aptp/factor.rs` pass (935 lines of tests)
- [ ] T040 [US1] Verify ALL multifrontal integration tests in `tests/multifrontal.rs` pass
- [ ] T041 [US1] Verify ALL end-to-end solve tests in `tests/solve.rs` pass
- [ ] T042 [US1] Run `cargo clippy` and fix any warnings introduced by new code

**Checkpoint**: Two-level APTP is the sole kernel. All existing tests pass. Dense matrices of size 512+ factor correctly with reconstruction < 1e-12.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Benchmarking, SuiteSparse validation, documentation updates

- [ ] T043 Run CI SuiteSparse test suite (9 matrices, MatchOrderMetis ordering) via `cargo test -- --ignored suitesparse_ci` and compare backward error to Phase 7 baseline (no regression target)
- [ ] T044 Run full SuiteSparse test suite (67 matrices, MatchOrderMetis) via `cargo test -- --ignored --test-threads=1` and document results in `docs/ssids-log.md`
- [ ] T045 [P] Add isolated kernel benchmarks (dense APTP on matrices 128, 256, 512, 1024, 2048) comparing two-level vs single-level (via `outer_block_size: usize::MAX` override) to `benches/solver_benchmarks.rs`
- [ ] T046 Document crossover point (front size where two-level outperforms single-level) based on benchmark results
- [ ] T047 Update `docs/ssids-plan.md` Phase 8.1 section: mark deliverables complete, record actual block sizes, record accuracy results
- [ ] T048 [P] Update `docs/ssids-log.md` with Phase 8.1 development log entry
- [ ] T049 Run `cargo fmt --check` and fix any formatting issues
- [ ] T050 Run quickstart.md validation: execute all verification commands listed in `specs/017-two-level-aptp/quickstart.md`

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
