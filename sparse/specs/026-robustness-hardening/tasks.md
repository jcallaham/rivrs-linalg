# Tasks: Robustness — Testing & Hardening

**Input**: Design documents from `specs/026-robustness-hardening/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, contracts/

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Add proptest dependency, create new module files, establish test infrastructure foundation

- [X] T001 Add `proptest = "1.4"` to `[dev-dependencies]` in `Cargo.toml`
- [X] T002 [P] Create `src/testing/perturbations.rs` with module skeleton (empty functions, imports, `#[cfg(feature = "test-util")]` gate) and register in `src/testing/mod.rs`
- [X] T003 [P] Create `src/testing/strategies.rs` with module skeleton (proptest imports, empty strategy functions) and register in `src/testing/mod.rs`
- [X] T004 Verify `cargo test` passes with new empty modules and proptest dependency (no regressions)

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Implement perturbation helpers that US1 audit gaps and US2 torture tests both depend on

**CRITICAL**: Perturbation helpers must be complete before torture tests (US2) can be written

- [X] T005 Implement `cause_delays(matrix: &mut Mat<f64>, rng: &mut impl Rng)` in `src/testing/perturbations.rs` — multiply n/8 random rows by 1000, n/8 random entries by 1000, first oversized row past first block when n > block_size. Include unit tests verifying symmetry preservation and that at least n/8 rows are scaled.
- [X] T006 Implement `make_singular(matrix: &mut Mat<f64>, col1: usize, col2: usize)` in `src/testing/perturbations.rs` — copy col1 scaled to col2, symmetrize. Include unit tests verifying rank deficiency (eigenvalue near zero) and symmetry preservation.
- [X] T007 Implement `make_dblk_singular(matrix: &mut Mat<f64>, block_row: usize, block_size: usize)` in `src/testing/perturbations.rs` — delegates to `make_singular` on first and last columns of specified diagonal block. Include unit tests verifying the targeted block becomes singular.
- [X] T008 Verify `cargo test` passes — all new perturbation helper unit tests green, existing tests unaffected

**Checkpoint**: Perturbation helpers ready — US1 audit and US2 torture tests can now proceed

---

## Phase 3: User Story 1 — SPRAL Test Parity Audit (Priority: P1) MVP

**Goal**: Produce a two-directional audit document mapping SPRAL kernel test scenarios to rivrs-sparse coverage, evaluate each existing test's independent value, prune redundant tests, and fill any coverage gaps with new tests.

**Independent Test**: Review the audit document and confirm (a) every SPRAL scenario mapped, (b) every removed test has documented rationale, (c) `cargo test` still passes.

### Implementation for User Story 1

- [X] T009 [US1] Catalog all SPRAL `ldlt_app.cxx` test scenarios (48 deterministic in categories A–H + 2 torture suites) in `docs/spral-test-audit.md` — for each scenario: test name, what it tests, matrix size, key assertions, perturbation used. Reference `research.md` Decision 1 and `/opt/references/spral/tests/ssids/kernels/ldlt_app.cxx`.
- [X] T010 [US1] Catalog all SPRAL `ldlt_tpp.cxx` test scenarios (23 deterministic in categories A–D + 2 torture suites) in `docs/spral-test-audit.md` — same format as T009. Reference `/opt/references/spral/tests/ssids/kernels/ldlt_tpp.cxx`.
- [X] T011 [US1] Map each SPRAL scenario to rivrs-sparse coverage in `docs/spral-test-audit.md` — for each SPRAL test: status (Covered / N/A / Gap), rivrs-sparse test name or N/A rationale. Use `research.md` SPRAL-to-rivrs Architecture Mapping table. Note: `app_aggressive` maps to N/A (not implemented in rivrs).
- [X] T012 [US1] Evaluate each existing rivrs-sparse test for independent value in `docs/spral-test-audit.md` — for each of the ~524 tests: status (Retain / Remove), brief rationale. Apply pruning criteria from `research.md` Decision 6. Tests covering meaningful scenarios not in SPRAL are explicitly retained.
- [X] T013 [US1] Fill identified gaps: for each Gap entry in the SPRAL mapping, add a new test to the appropriate test file (e.g., `src/aptp/factor.rs` for kernel-level gaps, `tests/` for integration gaps). Update audit document to reflect gaps filled.
- [X] T014 [US1] Remove or consolidate tests marked Remove in T012. For each removal: verify a neighboring test provides equivalent or broader coverage. Run `cargo test` after each batch of removals to ensure no regression.
- [X] T015 [US1] Run full validation: `cargo test` (all non-ignored pass), `cargo test -- --ignored --test-threads=1` (SuiteSparse 65/65 pass), `cargo clippy --all-targets`. Update audit document with final test count and summary.

**Checkpoint**: Audit document complete, test suite pruned, all tests pass. US1 independently verifiable.

---

## Phase 4: User Story 2 — SPRAL-Style Torture Testing (Priority: P1)

**Goal**: Add dense kernel-level torture tests with 500+ random instances per configuration, exercising the full pivot decision tree via adversarial perturbations.

**Independent Test**: Run `cargo test -- --ignored torture --test-threads=1` — all instances complete with zero panics and backward error within tolerance.

### Implementation for User Story 2

- [ ] T016 [US2] Implement `TortureTestConfig` struct in `src/testing/perturbations.rs` — fields: `num_instances`, `delay_probability` (0.70), `singular_probability` (0.20), `dblk_singular_probability` (0.10), `matrix_sizes`, `backward_error_threshold` (5e-11), `seed`. Include `Default` impl with SPRAL-matching values.
- [ ] T017 [US2] Implement APP torture test helper `ldlt_app_torture_test(config: &TortureTestConfig)` in `src/aptp/factor.rs` `#[cfg(test)]` module — generate random symmetric indefinite matrix, apply probabilistic perturbation (70/20/10 split), call `aptp_factor_in_place` with `num_fully_summed >= inner_block_size` (complete pivoting path), check backward error < threshold for non-singular, check inertia n_zero > 0 for singular, check q1+q2 = m. Loop for `num_instances` with reproducible seeded RNG.
- [ ] T018 [US2] Implement TPP torture test helper `ldlt_tpp_torture_test(config: &TortureTestConfig)` in `src/aptp/factor.rs` `#[cfg(test)]` module — same as T017 but call with `num_fully_summed < inner_block_size` (TPP path), perturbation mix: 70% delays + 20% singular (no dblk_singular for TPP, matching SPRAL). Additionally check L growth bound: max |L_ij| ≤ 1/threshold.
- [ ] T019 [US2] Add `#[test] #[ignore]` torture test entry points in `src/aptp/factor.rs` — APP sizes: (32,32), (64,64), (128,128), (128,48), (256,256); TPP sizes: (8,4), (33,21), (100,100), (100,50). 500 instances per configuration. Each with fixed seed for reproducibility. All kernel-level torture tests stay in `src/aptp/factor.rs` (not `tests/`) because they call `pub(crate)` functions (`aptp_factor_in_place`, `tpp_factor_as_primary`).
- [ ] T021 [US2] Run torture tests: `cargo test -- --ignored torture --test-threads=1`. Fix any failures (panics, backward error violations, inertia inconsistencies). Document any solver fixes required.
- [ ] T022 [US2] Verify no regressions: `cargo test` (non-ignored tests pass), `cargo clippy --all-targets`.

**Checkpoint**: Torture tests pass — 500+ instances per config, zero panics, backward error < 5e-11. US2 independently verifiable.

---

## Phase 5: User Story 3 — Property-Based Testing (Priority: P2)

**Goal**: Add proptest-driven property tests that verify structural invariants (backward error, inertia, permutation validity) across randomly generated matrices with automatic shrinking on failure.

**Independent Test**: Run `cargo test property` — all properties hold across 256+ generated instances (default proptest count).

### Implementation for User Story 3

- [ ] T023 [P] [US3] Implement `arb_symmetric_pd(size_range)` and `arb_symmetric_indefinite(size_range)` strategies in `src/testing/strategies.rs` — wrap existing `generate_random_symmetric` generator, return dense `Mat<f64>`, shrink toward smaller sizes. Include basic unit tests that generate a few matrices and verify symmetry + definiteness.
- [ ] T024 [P] [US3] Implement `arb_sparse_symmetric(size_range, density_range)` strategy in `src/testing/strategies.rs` — generate sparse symmetric CSC matrices suitable for `SparseLDLT`, shrink toward smaller sizes and lower density. Include unit tests.
- [ ] T025 [US3] Implement kernel-level property tests in `src/aptp/factor.rs` `#[cfg(test)]` module — `proptest!` block with properties: (1) reconstruction error < 1e-12 for PD matrices, (2) inertia sum = n for indefinite matrices, (3) permutation valid (each index once, fwd/inv are inverses), (4) no panics on perturbed matrices (using perturbation helpers from Phase 2).
- [ ] T026 [US3] Implement end-to-end property tests in `tests/property.rs` — `proptest!` block with properties: (1) backward error < 5e-11 for PD sparse matrices through `SparseLDLT::solve_full`, (2) backward error < 5e-11 or clean error for indefinite sparse matrices, (3) inertia consistency through full pipeline. Use `arb_sparse_symmetric` strategy.
- [ ] T027 [US3] Run property tests: `cargo test property`. Fix any failures (property violations trigger automatic shrinking — investigate minimal counterexamples). Tune proptest config if needed (case count, max shrink iters).
- [ ] T028 [US3] Verify no regressions: `cargo test` (all non-ignored pass), `cargo clippy --all-targets`.

**Checkpoint**: Property tests pass — structural invariants verified across 256+ generated matrices. US3 independently verifiable.

---

## Phase 6: User Story 4 — Adversarial & Edge-Case Input Testing (Priority: P2)

**Goal**: Add tests for malformed, extreme, and degenerate inputs ensuring zero panics and clean error returns. Harden solver where needed.

**Independent Test**: Run `cargo test adversarial` — all edge cases produce correct results or clean errors, zero panics.

### Implementation for User Story 4

- [ ] T029 [US4] Create `tests/adversarial.rs` with edge-case tests through `SparseLDLT` API — test cases: (1) 0×0 empty matrix → clean error, (2) 1×1 nonzero diagonal → correct solve, (3) 1×1 zero diagonal → error or zero-inertia factorization, (4) pure diagonal matrix → trivial correct factorization, (5) arrowhead matrix → correct factorization, (6) matrix requiring all 2×2 pivots (no stable 1×1 pivots) → correct factorization with all 2×2 pivots, (7) matrix with exact numerical cancellation during elimination → correct factorization or clean error.
- [ ] T030 [US4] Add extreme value tests in `tests/adversarial.rs` — test cases: (1) near-overflow entries (1e308) → error or correct, (2) near-underflow entries (1e-308) → error or correct, (3) matrix with NaN entries → clean error, (4) matrix with Inf entries → clean error. Each must not panic.
- [ ] T031 [US4] Add structural validity tests in `tests/adversarial.rs` — test cases: (1) non-symmetric sparsity pattern → clean error, (2) matrix with disconnected components → correct factorization, (3) matrix at power-of-2 boundaries (32, 64, 128, 256, 512) → correct factorization, (4) matrix where MC64 matching fails to find a perfect matching → clean error or fallback to unscaled ordering. Each must not panic.
- [ ] T032 [US4] Run adversarial tests, identify panics, and fix: add defensive guards in `src/aptp/solver.rs` (e.g., 0×0 check in `analyze_with_matrix`, NaN/Inf scan in `factor`) and `src/error.rs` (new error variants if needed). Each fix should be minimal — only guard the specific panic path.
- [ ] T033 [US4] Re-run adversarial tests after fixes: all tests pass (correct results or clean errors, zero panics). Verify no regressions: `cargo test` (all non-ignored pass).

**Checkpoint**: All adversarial inputs handled gracefully. US4 independently verifiable.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Final validation, documentation, cleanup

- [ ] T034 Run full validation suite: `cargo test` (non-ignored), `cargo test -- --ignored --test-threads=1` (SuiteSparse + torture), `cargo clippy --all-targets`, `cargo clippy --all-targets --features diagnostic`
- [ ] T035 [P] Update `docs/ssids-log.md` with Phase 9.2 entry: what was built, test counts before/after pruning, torture test results, property test results, adversarial findings and fixes
- [ ] T036 [P] Update `docs/ssids-plan.md` Phase 9.2 section with STATUS: COMPLETE and results summary
- [ ] T037 [P] Update `CLAUDE.md` Current Implementation Status with Phase 9.2 completion summary
- [ ] T038 Verify success criteria: SC-001 (audit 100% coverage), SC-002 (torture 500+ zero panics), SC-003 (proptest 1000+ no violations), SC-004 (adversarial zero panics), SC-005 (SuiteSparse 65/65 pass), SC-006 (test suite leaner)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on Setup (T001–T004)
- **US1 Audit (Phase 3)**: Depends on Foundational (T005–T008) for perturbation helpers (needed to fill gaps). Can start T009–T012 (audit/catalog) in parallel with Phase 2.
- **US2 Torture (Phase 4)**: Depends on Foundational (T005–T008) for perturbation helpers
- **US3 Property (Phase 5)**: Depends on Setup (T001–T004) for proptest dependency. Independent of US1/US2.
- **US4 Adversarial (Phase 6)**: Depends on Setup (T001–T004) only. Independent of US1/US2/US3.
- **Polish (Phase 7)**: Depends on all user stories complete

### User Story Dependencies

- **US1 (Audit)**: Can start cataloging (T009–T012) immediately; gap-filling (T013) needs perturbation helpers
- **US2 (Torture)**: Needs perturbation helpers (Phase 2). Independent of US1.
- **US3 (Property)**: Needs proptest dep (Phase 1). Independent of US1/US2.
- **US4 (Adversarial)**: Needs only project structure (Phase 1). Independent of all other stories.

### Parallel Opportunities

Within Phase 2: T005, T006, T007 can run in parallel (different functions, same file — separate development, merge)

Within Phase 3: T009 and T010 can run in parallel (different SPRAL files)

Within Phase 5: T023 and T024 can run in parallel (different strategy functions)

Across Phases: US3 (Phase 5) and US4 (Phase 6) are fully independent and can run in parallel with each other and with US2 (Phase 4)

---

## Implementation Strategy

### MVP First (US1 Audit + US2 Torture)

1. Complete Phase 1: Setup (T001–T004)
2. Complete Phase 2: Foundational perturbation helpers (T005–T008)
3. Complete Phase 3: US1 Audit (T009–T015)
4. Complete Phase 4: US2 Torture tests (T016–T022)
5. **STOP and VALIDATE**: Audit complete, torture tests pass, test suite pruned

### Incremental Delivery

1. Setup + Foundational → infrastructure ready
2. US1 Audit → audit document + pruned suite (MVP documentation deliverable)
3. US2 Torture → 500+ instance stress coverage (MVP test deliverable)
4. US3 Property → proptest structural invariants (complementary coverage)
5. US4 Adversarial → edge-case robustness (hardening deliverable)
6. Polish → docs, plan updates, final validation

---

## Notes

- Torture tests use `#[ignore]` — they won't slow down normal `cargo test`
- Property tests run with normal `cargo test` but are fast (256 cases default)
- Adversarial tests may uncover panics → T032 fixes are expected, not bugs
- Audit document is a reference artifact, not generated code — it lives in `docs/`
- All new test infrastructure code is behind `test-util` feature flag
- Commit after each task or logical group; never accumulate large unverified changes
