# Tasks: Dense APTP Factorization Kernel

**Input**: Design documents from `/specs/014-dense-aptp-kernel/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/api.md

**Tests**: Required by Constitution Principle III (TDD — NON-NEGOTIABLE). Tests MUST be written and verified to FAIL before implementation begins.

**Organization**: Tasks grouped by user story. US1 = core 1x1 factorization with delay-only fallback (MVP). US2 = 2x2 Bunch-Kaufman fallback + configurable threshold. US3 = diagnostics (statistics, pivot log, inertia).

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup (Module Scaffolding)

**Purpose**: Create factor module with type definitions, wire into aptp module

- [X] T001 Define all public types (AptpOptions, AptpFallback, AptpFactorResult, AptpFactorization, AptpStatistics, AptpPivotRecord) with fields, Default impls, and rustdoc in `src/aptp/factor.rs`. Include module-level doc comment citing Duff, Hogg & Lopez (2020). All function bodies should be `todo!()` stubs.
- [X] T002 Add `pub mod factor;` to `src/aptp/mod.rs` and re-export public types: AptpOptions, AptpFallback, AptpFactorResult, AptpFactorization, AptpStatistics, AptpPivotRecord, aptp_factor_in_place, aptp_factor
- [X] T003 Verify `cargo build` succeeds with the new module (stub functions with `todo!()` are fine)

---

## Phase 2: Foundational (Test Infrastructure)

**Purpose**: Dense reconstruction test helper needed by ALL user stories. Existing NumericalValidator works with sparse matrices; the dense kernel needs a dense equivalent.

- [X] T004 Implement `reconstruct_dense_ldlt(l: &Mat<f64>, d: &MixedDiagonal, perm: &[usize]) -> Mat<f64>` helper function for tests. This computes `P^T L D L^T P` from the convenience API result. Place in `src/aptp/factor.rs` as `#[cfg(test)]` helper or in test module. Verify it produces identity-like results for trivial 2x2 diagonal input.
- [X] T005 Implement `dense_reconstruction_error(original: &Mat<f64>, l: &Mat<f64>, d: &MixedDiagonal, perm: &[usize]) -> f64` that computes `||A - P^T L D L^T P|| / ||A||`. This is the primary correctness oracle for all tests.

**Checkpoint**: Test infrastructure ready — user story implementation can begin

---

## Phase 3: User Story 1 — Core 1x1 Factorization with Delay Fallback (Priority: P1) MVP

**Goal**: Factor dense symmetric matrices using optimistic 1x1 pivots with a posteriori stability checking. Columns that fail the stability bound are delayed. No 2x2 pivots — delay-only fallback. This is the minimum viable APTP kernel.

**Independent Test**: Factor PD matrices (all 1x1, zero delays) and verify reconstruction error < 10^-12. Factor zero matrix (all delayed) and verify empty factorization returned.

### Tests for User Story 1

> **Constitution III**: Write tests FIRST, verify they FAIL before implementation

- [X] T006 [US1] Write test `test_1x1_trivial_diagonal` in `src/aptp/factor.rs` (tests module): 2x2 diagonal matrix [[4,0],[0,9]], verify both columns are 1x1 pivots, no delays, reconstruction error < 1e-12
- [X] T007 [US1] Write test `test_1x1_positive_definite_3x3` in `src/aptp/factor.rs`: 3x3 PD matrix (e.g., [[4,2,1],[2,5,3],[1,3,6]]), verify all 1x1 pivots, zero delays, reconstruction < 1e-12
- [X] T008 [US1] Write test `test_all_delayed_zero_matrix` in `src/aptp/factor.rs`: n x n zero matrix, verify num_eliminated == 0, all columns delayed, stats.num_delayed == n
- [X] T009 [US1] Write test `test_1x1_singularity_detection` in `src/aptp/factor.rs`: matrix with one near-zero diagonal entry (< `small` threshold), verify that column is delayed while others succeed
- [X] T010 [US1] Write test `test_stability_bound_enforced` in `src/aptp/factor.rs`: matrix where naive 1x1 pivot would produce |l_ij| > 1/threshold, verify column is delayed (not accepted with large L entries)
- [X] T011 [US1] Write test `test_1x1_matrix` in `src/aptp/factor.rs`: 1x1 matrix [[5.0]], verify single 1x1 pivot, reconstruction exact
- [X] T012 [US1] Write test `test_statistics_sum_invariant` in `src/aptp/factor.rs`: for any factorization result, verify num_1x1 + 2*num_2x2 + num_delayed == n
- [X] T013 [US1] Verify all US1 tests compile and FAIL (functions are still `todo!()` stubs)

### Implementation for User Story 1

- [X] T014 [US1] Implement input validation in `aptp_factor_in_place()` in `src/aptp/factor.rs`: check matrix is square, num_fully_summed <= nrows, threshold in valid range. Return `SparseError::InvalidInput` on failure.
- [X] T015 [US1] Implement `try_1x1_pivot()` private function in `src/aptp/factor.rs`: compute L column entries by dividing column k by diagonal d_kk, check singularity (|d_kk| < small → delay), check stability bound (max|l_ik| < 1/threshold), return Ok(d_value) or Err(max_l_entry)
- [X] T016 [US1] Implement `update_schur_1x1()` private function in `src/aptp/factor.rs`: symmetric rank-1 update of trailing submatrix using faer's `matmul_with_conj` or `rank_update`. Update both the remaining fully-summed columns AND the contribution block.
- [X] T017 [US1] Implement the main column loop in `aptp_factor_in_place()` in `src/aptp/factor.rs`: iterate columns 0..num_fully_summed, attempt 1x1 pivot via `try_1x1_pivot()`, on success record in MixedDiagonal and call `update_schur_1x1()`, on failure mark as Delayed. Track permutation as Vec\<usize\> (identity initially, delayed columns moved to end after loop).
- [X] T018 [US1] Implement post-loop finalization in `aptp_factor_in_place()`: collect delayed column indices, build final permutation (eliminated columns first, delayed last), construct AptpFactorResult with MixedDiagonal, perm, num_eliminated, delayed_cols, stats, pivot_log.
- [X] T019 [US1] Implement `extract_l()` private function in `src/aptp/factor.rs`: read L entries from the lower triangle of the mutated matrix respecting the column permutation and num_eliminated. Build a Mat\<f64\> with unit diagonal.
- [X] T020 [US1] Implement `aptp_factor()` convenience wrapper in `src/aptp/factor.rs`: copy input matrix, call `aptp_factor_in_place` with num_fully_summed = n, extract L via `extract_l()`, convert perm Vec to Perm\<usize\> via `perm_from_forward()`, return AptpFactorization.
- [X] T021 [US1] Run all US1 tests — verify T006-T012 all PASS. Fix any failures. Run `cargo clippy` and `cargo fmt --check`.

**Checkpoint**: Core APTP kernel works for PD matrices (all 1x1) and handles delayed columns. Reconstruction error < 10^-12 verified. MVP complete.

---

## Phase 4: User Story 2 — 2x2 Bunch-Kaufman Fallback & Configurable Threshold (Priority: P2)

**Goal**: When a 1x1 pivot fails stability check and BunchKaufman fallback is selected, attempt a 2x2 pivot with the best partner column. Support multiple threshold values for tuning.

**Independent Test**: Factor a 4x4 indefinite matrix known to require a 2x2 pivot. Verify reconstruction < 10^-12, correct TwoByTwo pivot type, and that BK fallback reduces delays compared to Delay-only fallback on the same matrix.

### Tests for User Story 2

- [X] T022 [US2] Write test `test_2x2_pivot_known_indefinite` in `src/aptp/factor.rs`: 4x4 symmetric indefinite matrix (e.g., [[1,0,0,2],[0,1,2,0],[0,2,-1,0],[2,0,0,-1]]) where 2x2 pivot is needed, verify at least one TwoByTwo pivot, reconstruction < 1e-12
- [X] T023 [US2] Write test `test_2x2_stability_test` in `src/aptp/factor.rs`: matrix where 2x2 block determinant condition (|det| >= 0.5*|a21|^2) should be checked. Construct case where 2x2 passes and case where it fails (both columns delayed).
- [X] T024 [US2] Write test `test_bk_vs_delay_fallback` in `src/aptp/factor.rs`: same indefinite matrix factored with AptpFallback::BunchKaufman and AptpFallback::Delay. Verify BK produces fewer delayed columns. Both must satisfy reconstruction tolerance.
- [X] T025 [US2] Write test `test_strict_threshold` in `src/aptp/factor.rs`: factor same matrix with threshold=0.01 and threshold=0.5. Verify stricter threshold produces more delays but both factorizations are correct.
- [X] T026 [US2] Write test `test_permutation_valid` in `src/aptp/factor.rs`: after factorization with 2x2 pivots (column swaps), verify perm is a valid bijection on 0..n and that P^T A P = L D L^T reconstruction holds.
- [X] T027 [US2] Verify all US2 tests compile and FAIL (2x2 path not yet implemented)

### Implementation for User Story 2

- [X] T028 [US2] Implement `select_2x2_partner()` private function in `src/aptp/factor.rs`: search uneliminated fully-summed columns for the one with largest |a_ik| (off-diagonal in column k). Return Some(partner_index) or None if no suitable partner.
- [X] T029 [US2] Implement `try_2x2_pivot()` private function in `src/aptp/factor.rs`: given columns k and partner, compute 2x2 block [[a_kk, a_pk], [a_pk, a_pp]], test determinant condition (|det| >= 0.5*|a_pk|^2 from Algorithm 4.1), compute L entries for both columns, check stability bound, return Ok(Block2x2) or Err(()).
- [X] T030 [US2] Implement `update_schur_2x2()` private function in `src/aptp/factor.rs`: precompute W = L_cols * D_22^{-1} (2 columns), then rank-2 update of trailing submatrix using faer's `matmul_with_conj`. Handle column swap (partner moved adjacent to k).
- [X] T031 [US2] Integrate BunchKaufman fallback into the main column loop in `aptp_factor_in_place()`: when 1x1 fails and options.fallback == BunchKaufman, call `select_2x2_partner()` then `try_2x2_pivot()`. On 2x2 success: swap partner column adjacent, record TwoByTwo in MixedDiagonal, call `update_schur_2x2()`, advance by 2. On 2x2 failure: delay column (and partner if attempted).
- [X] T032 [US2] Run all US1 + US2 tests — verify T006-T012 still PASS (no regression) and T022-T026 all PASS. Run `cargo clippy` and `cargo fmt --check`.

**Checkpoint**: Full APTP kernel with 1x1 and 2x2 pivots. Both BunchKaufman and Delay fallback strategies work. Configurable threshold verified.

---

## Phase 5: User Story 3 — Factorization Diagnostics (Priority: P3)

**Goal**: Report accurate summary statistics (pivot counts, max |L| entry) and per-column pivot log. Enable inertia computation from D factor.

**Independent Test**: Factor a PD matrix and verify stats report zero 2x2 pivots, zero delays. Factor an indefinite matrix and verify stats.num_1x1 + 2*stats.num_2x2 + stats.num_delayed == n. Verify max_l_entry matches independent scan.

### Tests for User Story 3

- [X] T033 [US3] Write test `test_pd_statistics` in `src/aptp/factor.rs`: PD matrix, verify stats.num_1x1 == n, num_2x2 == 0, num_delayed == 0, max_l_entry < 1/threshold
- [X] T034 [US3] Write test `test_max_l_entry_accuracy` in `src/aptp/factor.rs`: factor an indefinite matrix, independently scan all L entries to find actual max, verify stats.max_l_entry matches within floating-point tolerance
- [X] T035 [US3] Write test `test_pivot_log_completeness` in `src/aptp/factor.rs`: verify pivot_log has one entry per fully-summed column, each with correct pivot_type matching D factor, and was_fallback flag set correctly for BK-attempted columns
- [X] T036 [US3] Write test `test_inertia_from_d` in `src/aptp/factor.rs`: factor a known indefinite matrix (e.g., diagonal with values [3, -2, 1, -4, 0]), compute inertia via d.compute_inertia(), verify positive/negative/zero counts match expected
- [X] T037 [US3] Verify all US3 tests compile and FAIL (statistics not yet accumulated)

### Implementation for User Story 3

- [X] T038 [US3] Implement statistics accumulation in the main loop of `aptp_factor_in_place()` in `src/aptp/factor.rs`: track running max_l_entry across all pivot attempts, increment num_1x1/num_2x2/num_delayed counters, build AptpStatistics at end of loop
- [X] T039 [US3] Implement pivot log recording in the main loop: after each column's pivot decision, push an AptpPivotRecord with col (original index), pivot_type, max_l_entry for that column, and was_fallback flag
- [X] T040 [US3] Run all US1 + US2 + US3 tests — verify no regressions and T033-T036 all PASS. Run `cargo clippy` and `cargo fmt --check`.

**Checkpoint**: All three user stories complete. Full APTP kernel with diagnostics.

---

## Phase 6: Polish & Comprehensive Validation

**Purpose**: Stress testing, edge cases, integration with existing test infrastructure, documentation

- [X] T041 [P] Write test `test_random_pd_matrices` in `src/aptp/factor.rs`: generate 100+ random symmetric PD matrices (sizes 5, 10, 20, 50) using test-util generators, verify all produce zero delays, all 1x1, reconstruction < 1e-12 (SC-002, SC-004)
- [X] T042 [P] Write test `test_random_indefinite_matrices` in `src/aptp/factor.rs`: generate 100+ random symmetric indefinite matrices (sizes 5, 10, 20, 50), verify reconstruction < 1e-12, stability bound satisfied, statistics invariant holds (SC-002, SC-008)
- [X] T043 [P] Write test `test_edge_case_extreme_thresholds` in `src/aptp/factor.rs`: test threshold near 0.0 (accept almost everything) and near 1.0 (delay almost everything), verify valid factorizations in both cases
- [X] T044 [P] Write test `test_partial_factorization` in `src/aptp/factor.rs`: create a 10x10 dense matrix, call aptp_factor_in_place with num_fully_summed=6, verify only first 6 columns attempted, contribution block (rows/cols 6..10) updated correctly, reconstruction of eliminated portion is correct
- [X] T045 [P] Write test `test_both_fallback_strategies_valid` (SC-005) in `src/aptp/factor.rs`: factor 20+ indefinite matrices with both BunchKaufman and Delay fallback, verify both produce valid factorizations (reconstruction < 1e-12)
- [X] T046 Write integration test `test_hand_constructed_matrices` in `src/aptp/factor.rs` (#[ignore] test): load the 15 hand-constructed matrices from test-data/hand-constructed/ via the test infrastructure, convert to dense, factor each with aptp_factor(), verify reconstruction < 1e-12 (SC-001, SC-007)
- [X] T047 Write integration test `test_suitesparse_ci_dense` in `src/aptp/factor.rs` (#[ignore] test): load 10 SuiteSparse CI-subset matrices, convert to dense (if dimension <= 500), factor with aptp_factor(), verify reconstruction < 1e-12 (SC-003)
- [X] T048 Add comprehensive rustdoc to all public items in `src/aptp/factor.rs`: algorithm description citing Duff, Hogg & Lopez (2020), complexity analysis (O(n^3) time, O(n^2) space), SPRAL equivalent references, # Examples, # Errors, # Panics sections per Constitution IV
- [X] T049 Run full test suite: `cargo test` (unit + non-ignored), `cargo test -- --ignored --test-threads=1` (integration), `cargo clippy`, `cargo fmt --check`. Verify all pass.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 (T001-T003 complete)
- **US1 (Phase 3)**: Depends on Phase 2 (T004-T005 complete)
- **US2 (Phase 4)**: Depends on US1 complete (T021 verified) — 2x2 extends the 1x1 loop
- **US3 (Phase 5)**: Depends on US2 complete (T032 verified) — diagnostics wrap the full loop
- **Polish (Phase 6)**: Depends on US3 complete (T040 verified)

### User Story Dependencies

- **US1 (P1)**: Independent after Foundational — delivers working 1x1 kernel (MVP)
- **US2 (P2)**: Depends on US1 — adds 2x2 BK fallback to existing main loop
- **US3 (P3)**: Depends on US2 — adds statistics/logging to the complete main loop

Note: US2 and US3 depend on US1 because they modify the same function (`aptp_factor_in_place()`). The 2x2 fallback path integrates into the existing column loop, and diagnostics wrap the pivot decision logic.

### Within Each User Story

1. Tests MUST be written and verified to FAIL before implementation (Constitution III)
2. Internal helpers before main loop integration
3. In-place API before convenience wrapper
4. Verify all tests PASS before moving to next story

### Parallel Opportunities

- **Phase 1**: T001 and T002 are sequential (T002 depends on T001)
- **Phase 2**: T004 and T005 are sequential (T005 depends on T004)
- **Phase 3 tests**: T006-T012 can all be written in parallel [P] — they don't depend on each other
- **Phase 4 tests**: T022-T026 can all be written in parallel [P]
- **Phase 5 tests**: T033-T036 can all be written in parallel [P]
- **Phase 6**: T041-T045 can all be written in parallel [P] — different test functions in same file

---

## Parallel Example: User Story 1 Tests

```bash
# Write all US1 tests in parallel (different test functions, same file):
T006: test_1x1_trivial_diagonal
T007: test_1x1_positive_definite_3x3
T008: test_all_delayed_zero_matrix
T009: test_1x1_singularity_detection
T010: test_stability_bound_enforced
T011: test_1x1_matrix
T012: test_statistics_sum_invariant
```

## Parallel Example: Polish Phase

```bash
# Write all stress/edge case tests in parallel:
T041: test_random_pd_matrices
T042: test_random_indefinite_matrices
T043: test_edge_case_extreme_thresholds
T044: test_partial_factorization
T045: test_both_fallback_strategies_valid
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (module scaffolding) — T001-T003
2. Complete Phase 2: Foundational (reconstruction helper) — T004-T005
3. Complete Phase 3: User Story 1 (1x1 + delay) — T006-T021
4. **STOP and VALIDATE**: All PD matrices factor correctly, delayed columns work, reconstruction < 1e-12
5. This is a functional dense LDL^T kernel — can be used by Phase 6 for PD problems immediately

### Incremental Delivery

1. Setup + Foundational → Module compiles, test helpers ready
2. Add US1 → 1x1 factorization works → PD matrices solved correctly (MVP!)
3. Add US2 → 2x2 BK fallback → Indefinite matrices handled with fewer delays
4. Add US3 → Diagnostics → Statistics and inertia available for monitoring
5. Polish → Stress tests, integration, documentation → Phase 5 exit criteria met

### Success Criteria Mapping

| SC | Verified By |
|----|-------------|
| SC-001 (15 hand-constructed) | T046 |
| SC-002 (100+ random) | T041, T042 |
| SC-003 (SuiteSparse CI) | T047 |
| SC-004 (PD zero delays) | T041 |
| SC-005 (both fallbacks valid) | T045 |
| SC-006 (stats sum to n) | T012, T033 |
| SC-007 (inertia correct) | T036, T046 |
| SC-008 (max_l_entry accurate) | T034, T042 |

---

## Notes

- All implementation is in a single file: `src/aptp/factor.rs` (~500-800 lines expected)
- Only `src/aptp/mod.rs` is modified (re-exports) — no other existing files change
- `#[ignore]` tests (T046, T047) require test-data/ matrices — run with `--ignored` flag
- Constitution III TDD: every test task includes "verify FAIL before implementation"
- Commit after each phase checkpoint (T003, T005, T021, T032, T040, T049)
