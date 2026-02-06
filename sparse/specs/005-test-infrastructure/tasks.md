# Tasks: Core Test Infrastructure

**Input**: Design documents from `/specs/005-test-infrastructure/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, contracts/

**Tests**: Included — constitution Principle III (TDD) requires test-first workflow for all components.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3, US4)
- Include exact file paths in descriptions

## User Story Mapping

| Story | Priority | Title | Key Files |
|-------|----------|-------|-----------|
| US1 | P1 | Validate Solver Component Correctness | `src/testing/harness.rs` |
| US2 | P1 | Numerical Validation with Configurable Tolerances | `src/testing/validator.rs` |
| US3 | P2 | Generate Random Test Matrices | `src/testing/generators.rs` |
| US4 | P2 | Unified Test Case Management | `src/testing/cases.rs` |

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Configure Cargo feature flag and create module skeleton

- [ ] T001 Add `test-util` feature flag and self-dev-dependency to `Cargo.toml`: add `[features] test-util = []` section and add `rivrs-sparse = { path = ".", features = ["test-util"] }` under `[dev-dependencies]`
- [ ] T002 Create `src/testing/mod.rs` with feature-gated module declarations for `cases`, `harness`, `validator`, `generators` submodules and public re-exports
- [ ] T003 Add feature-gated `pub mod testing` to `src/lib.rs` with `#[cfg(feature = "test-util")]` gate
- [ ] T004 Verify setup compiles: run `cargo test --no-run` and `cargo build` (production build must not include testing module)

**Checkpoint**: `src/testing/` module exists, feature gate works, production build unaffected

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core types that ALL user stories depend on — `TestResult`, `MetricResult`, `TestKind`, `SolverTestCase`, `TestMatrixProperties`

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

- [ ] T005 [P] Implement `TestKind` enum (`Analyze`, `Factor`, `Solve`, `Roundtrip`) with `Debug`, `Clone`, `PartialEq` derives in `src/testing/harness.rs`
- [ ] T006 [P] Implement `MetricResult` struct (`name: String`, `value: f64`, `threshold: f64`, `passed: bool`) with `Debug`, `Clone` derives and `Display` impl that formats as `"{name}: {value:.2e} (threshold: {threshold:.2e}) — {PASS|FAIL}"` in `src/testing/harness.rs`
- [ ] T007 [P] Implement `TestResult` struct (`passed: bool`, `test_kind: TestKind`, `matrix_name: String`, `metrics: Vec<MetricResult>`, `diagnostics: Vec<String>`) with `Debug`, `Clone` derives and a `Display` impl summarizing pass/fail with all diagnostics in `src/testing/harness.rs`
- [ ] T008 [P] Implement `TestMatrixProperties` struct (fields: `size`, `nnz`, `symmetric`, `positive_definite`, `indefinite`, `difficulty`, `structure`, `source`, `category`) with `Debug`, `Clone` derives in `src/testing/cases.rs`
- [ ] T009 [P] Implement `SolverTestCase` struct (`name: String`, `matrix: SparseColMat<usize, f64>`, `properties: TestMatrixProperties`, `reference: Option<ReferenceFactorization>`) with `Debug`, `Clone` derives in `src/testing/cases.rs`
- [ ] T010 Write unit tests for foundational types: construct `MetricResult` for pass/fail cases, verify `Display` output includes metric name and values; construct `TestResult` with mixed pass/fail metrics and verify `passed` field is false when any metric fails. Tests in `src/testing/harness.rs` `#[cfg(test)] mod tests`
- [ ] T011 Verify `cargo test` passes with all foundational types

**Checkpoint**: Foundation ready — all shared types compile and have basic tests

---

## Phase 3: User Story 4 — Unified Test Case Management (Priority: P2)

**Goal**: Load test matrices from the registry into `SolverTestCase` format with filtering by source, category, difficulty, CI-subset, and reference availability.

**Independent Test**: Load all hand-constructed test cases → verify 15 returned with correct properties. Filter by category → verify correct subset. (SC-005, SC-006)

**Why US4 before US1/US2**: The `SolverTestCase` type and `load_test_cases()` function provide the data that the validator (US2) and harness trait (US1) operate on. Building this first enables realistic testing of US2 and US1.

### Tests for User Story 4

- [ ] T012 [P] [US4] Write test `load_hand_constructed_returns_15` in `src/testing/cases.rs`: call `load_test_cases(&TestCaseFilter::hand_constructed())`, assert 15 cases returned, each with `reference.is_some()` and `properties.source == "hand-constructed"`
- [ ] T013 [P] [US4] Write test `load_ci_subset_returns_9` in `src/testing/cases.rs`: call `load_test_cases(&TestCaseFilter::ci_subset())`, assert 9 cases returned, each with `properties.source == "suitesparse"`
- [ ] T014 [P] [US4] Write test `filter_by_category` in `src/testing/cases.rs`: filter by `"hand-constructed"` category, verify all returned cases have matching category; filter by `"hard-indefinite"`, verify all returned are indefinite
- [ ] T015 [P] [US4] Write test `require_reference_filters_correctly` in `src/testing/cases.rs`: call with `require_reference(true)`, verify all returned cases have `reference.is_some()`
- [ ] T016 [P] [US4] Write test `load_performance_under_100ms` in `src/testing/cases.rs`: time a single matrix load via `load_test_cases`, assert elapsed < 100ms per matrix (FR-011)
- [ ] T017 [US4] Verify US4 tests FAIL (types exist but functions not yet implemented)

### Implementation for User Story 4

- [ ] T018 [US4] Implement `TestCaseFilter` struct and builder methods (`all()`, `hand_constructed()`, `ci_subset()`, `with_source()`, `with_category()`, `with_difficulty()`, `ci_only()`, `require_reference()`) in `src/testing/cases.rs`
- [ ] T019 [US4] Implement conversion from `registry::MatrixMetadata` + `registry::MatrixProperties` → `TestMatrixProperties` (helper function or `From` impl) in `src/testing/cases.rs`
- [ ] T020 [US4] Implement `load_test_cases(filter: &TestCaseFilter) -> Result<Vec<SolverTestCase>, SparseError>` in `src/testing/cases.rs`: load registry, apply filter predicates, call `registry::load_test_matrix()` for each match, convert to `SolverTestCase`, skip missing .mtx files gracefully
- [ ] T021 [US4] Verify all US4 tests PASS: `cargo test testing::cases`

**Checkpoint**: `load_test_cases()` returns correct, filtered `SolverTestCase` instances from the registry

---

## Phase 4: User Story 2 — Numerical Validation with Configurable Tolerances (Priority: P1)

**Goal**: Provide a `NumericalValidator` that wraps existing `validate.rs` functions with configurable tolerances and returns structured `MetricResult`/`TestResult` values with diagnostic context.

**Independent Test**: Validate hand-constructed matrices with default tolerances → all pass with recon < 1e-12. Perturb a factorization → validator reports failure with exact metric, value, and threshold. Custom tolerance → relaxed check passes. (SC-001, SC-002)

### Tests for User Story 2

- [ ] T022 [P] [US2] Write test `default_tolerances_match_constitution` in `src/testing/validator.rs`: create `NumericalValidator::new()`, assert `reconstruction_tol == 1e-12` and `backward_error_tol == 1e-10` (FR-004)
- [ ] T023 [P] [US2] Write test `check_reconstruction_passes_for_hand_constructed` in `src/testing/validator.rs`: load `arrow-5-pd` as `SolverTestCase`, call `check_reconstruction()`, assert `result.passed == true` and `result.value < 1e-12`
- [ ] T024 [P] [US2] Write test `check_reconstruction_fails_for_perturbed` in `src/testing/validator.rs`: load `arrow-5-pd`, perturb an L entry, call `check_reconstruction()`, assert `result.passed == false` and `result.name == "reconstruction_error"` (FR-005)
- [ ] T025 [P] [US2] Write test `custom_tolerance_relaxes_check` in `src/testing/validator.rs`: create validator with `with_reconstruction_tol(1e-6)`, validate a moderately perturbed factorization, assert it passes
- [ ] T026 [P] [US2] Write test `check_inertia_pass_and_fail` in `src/testing/validator.rs`: test matching inertia returns `passed == true`; mismatched inertia returns `passed == false` with diagnostic
- [ ] T027 [P] [US2] Write test `validate_all_15_hand_constructed` in `src/testing/validator.rs`: load all hand-constructed cases via `load_test_cases`, run `validate_factorization()` on each, assert all pass (SC-001)
- [ ] T028 [US2] Verify US2 tests FAIL (struct exists but methods not yet implemented)

### Implementation for User Story 2

- [ ] T029 [US2] Implement `NumericalValidator` struct with `new()`, `Default`, `with_reconstruction_tol()`, `with_backward_error_tol()` builder methods in `src/testing/validator.rs`
- [ ] T030 [US2] Implement `check_reconstruction(&self, matrix, reference) -> MetricResult` in `src/testing/validator.rs`: call `validate::reconstruction_error()`, compare against `self.reconstruction_tol`, populate `MetricResult` with name, value, threshold, passed
- [ ] T031 [US2] Implement `check_backward_error(&self, matrix, x, b) -> MetricResult` in `src/testing/validator.rs`: call `validate::backward_error()`, compare against `self.backward_error_tol`, populate `MetricResult`
- [ ] T032 [US2] Implement `check_inertia(&self, computed, expected) -> MetricResult` in `src/testing/validator.rs`: call `validate::check_inertia()`, return `MetricResult` with `name = "inertia"`, value = 0.0 or 1.0 (match/mismatch), threshold = 0.5
- [ ] T033 [US2] Implement `validate_factorization(&self, case: &SolverTestCase) -> TestResult` in `src/testing/validator.rs`: if reference available, run `check_reconstruction` + `check_inertia`; aggregate into `TestResult` with diagnostics including matrix name and dimensions (FR-005)
- [ ] T034 [US2] Verify all US2 tests PASS: `cargo test testing::validator`

**Checkpoint**: `NumericalValidator` validates factorizations with configurable tolerances and structured diagnostics

---

## Phase 5: User Story 1 — Validate Solver Component Correctness (Priority: P1)

**Goal**: Define the `SolverTest` trait and demonstrate it with a mock solver that uses identity permutation and the reference factorization directly.

**Independent Test**: Implement a trivial mock solver, run all four `SolverTest` methods on hand-constructed test cases, all should pass. (SC-004)

### Tests for User Story 1

- [ ] T035 [P] [US1] Write test `mock_solver_test_roundtrip_passes` in `src/testing/harness.rs`: implement `MockSolver` that uses identity permutation and reference factorization, run `test_roundtrip` on `arrow-5-pd`, assert `result.passed == true`
- [ ] T036 [P] [US1] Write test `mock_solver_test_analyze_passes` in `src/testing/harness.rs`: run `test_analyze` on `arrow-5-pd` with mock solver, assert valid permutation check passes
- [ ] T037 [P] [US1] Write test `mock_solver_test_factor_passes` in `src/testing/harness.rs`: run `test_factor` on `arrow-5-pd` with mock solver, assert reconstruction error metric passes
- [ ] T038 [P] [US1] Write test `mock_solver_test_solve_passes` in `src/testing/harness.rs`: run `test_solve` on `arrow-5-pd` with mock solver, assert backward error metric passes
- [ ] T039 [P] [US1] Write test `solver_test_all_hand_constructed` in `src/testing/harness.rs`: run mock solver `test_roundtrip` on all 15 hand-constructed cases, assert all pass (SC-004)
- [ ] T040 [US1] Verify US1 tests FAIL (trait defined but mock solver not yet implemented)

### Implementation for User Story 1

- [ ] T041 [US1] Define `SolverTest` trait with `test_analyze`, `test_factor`, `test_solve`, `test_roundtrip` method signatures in `src/testing/harness.rs` per contracts/testing-api.md
- [ ] T042 [US1] Implement `MockSolver` struct (test-only) that uses identity permutation and reference factorization from the test case, implementing all four `SolverTest` methods. `test_roundtrip` runs the full pipeline independently, using `NumericalValidator` for metrics. In `src/testing/harness.rs`
- [ ] T043 [US1] Verify all US1 tests PASS: `cargo test testing::harness`

**Checkpoint**: `SolverTest` trait is defined and validated with a mock implementation on all hand-constructed matrices

---

## Phase 6: User Story 3 — Generate Random Test Matrices (Priority: P2)

**Goal**: Provide random and structured sparse symmetric matrix generators with configurable properties (size, density, definiteness, pattern).

**Independent Test**: Generate PD matrix → verify symmetric, correct dimension, PD. Generate indefinite → verify mixed signs. Generate arrow/tridiagonal/banded → verify pattern. Performance < 1s for n ≤ 1000. (SC-003)

### Tests for User Story 3

- [ ] T044 [P] [US3] Write test `random_pd_matrix_properties` in `src/testing/generators.rs`: generate 100x100 PD matrix with ~500 nnz, assert symmetric (A == A^T via dense), correct dimension, PD (all diagonal entries > 0 after LDL^T), nnz within 20% of target
- [ ] T045 [P] [US3] Write test `random_indefinite_matrix_has_mixed_signs` in `src/testing/generators.rs`: generate 50x50 indefinite matrix, convert to dense, verify presence of both positive and negative diagonal entries (Gershgorin guarantee)
- [ ] T046 [P] [US3] Write test `generate_arrow_pattern` in `src/testing/generators.rs`: generate 20x20 arrow, verify first row/column are dense (n-1 off-diag entries), remaining rows have only diagonal
- [ ] T047 [P] [US3] Write test `generate_tridiagonal_pattern` in `src/testing/generators.rs`: generate 30x30 tridiagonal, verify each row has at most 3 nonzeros (diagonal + two neighbors)
- [ ] T048 [P] [US3] Write test `generate_banded_pattern` in `src/testing/generators.rs`: generate 40x40 banded with bandwidth 3, verify nonzeros only within band
- [ ] T049 [P] [US3] Write test `generation_performance_under_1s` in `src/testing/generators.rs`: generate 1000x1000 matrix, assert elapsed < 1 second (SC-003)
- [ ] T050 [P] [US3] Write test `infeasible_config_returns_error` in `src/testing/generators.rs`: request 1x1 indefinite matrix, assert error returned
- [ ] T051 [P] [US3] Write test `excessive_nnz_clamped` in `src/testing/generators.rs`: request 5x5 with 1000 nnz, verify actual nnz ≤ maximum possible (25 for full 5x5)
- [ ] T052 [US3] Verify US3 tests FAIL (function signatures exist but not implemented)

### Implementation for User Story 3

- [ ] T053 [US3] Implement `RandomMatrixConfig` struct in `src/testing/generators.rs`
- [ ] T054 [US3] Implement `generate_random_symmetric(config, rng)` in `src/testing/generators.rs`: generate random off-diagonal entries in lower triangle, symmetrize, apply diagonal dominance for PD or mixed signs for indefinite, clamp nnz if over max, build via `SparseColMat::try_new_from_triplets()`
- [ ] T055 [P] [US3] Implement `generate_arrow(size, positive_definite, rng)` in `src/testing/generators.rs`: create arrow sparsity pattern (dense first row/col + diagonal), apply diagonal dominance or mixed signs
- [ ] T056 [P] [US3] Implement `generate_tridiagonal(size, positive_definite, rng)` in `src/testing/generators.rs`: create tridiagonal pattern, apply diagonal dominance or mixed signs
- [ ] T057 [P] [US3] Implement `generate_banded(size, bandwidth, positive_definite, rng)` in `src/testing/generators.rs`: create banded pattern within specified bandwidth, apply diagonal dominance or mixed signs
- [ ] T058 [US3] Verify all US3 tests PASS: `cargo test testing::generators`

**Checkpoint**: All four generator functions produce valid matrices with correct structural properties

---

## Phase 7: Integration & Refactoring

**Purpose**: Refactor existing integration tests to use the new harness, demonstrating end-to-end value and eliminating duplicated validation logic.

- [ ] T059 Refactor `tests/hand_constructed.rs` to use `load_test_cases(TestCaseFilter::hand_constructed())` and `NumericalValidator` instead of manual registry calls and direct `validate::` function calls. Must produce identical pass/fail behavior to current test (SC-007, FR-012)
- [ ] T060 Refactor `tests/suitesparse_ci.rs` to use `load_test_cases(TestCaseFilter::ci_subset())` and `SolverTestCase` properties instead of manual registry loading and dimension checks (SC-007, FR-012)
- [ ] T061 Run `cargo test` and verify all existing tests still pass with identical behavior
- [ ] T062 Run `cargo test` with `--release` flag to verify optimized build works correctly
- [ ] T063 Run `cargo build` (without test-util feature) to verify production build excludes testing module

**Checkpoint**: Existing integration tests refactored; all tests pass; production build clean

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Final quality checks and documentation

- [ ] T064 [P] Add rustdoc comments to all public types and functions in `src/testing/mod.rs`, `harness.rs`, `validator.rs`, `cases.rs`, `generators.rs` — include `# Examples` sections per constitution Principle VII
- [ ] T065 [P] Run `cargo clippy -- -D warnings` and fix any warnings across all testing module files
- [ ] T066 Run `cargo test` full suite and verify total execution time increase is < 2 seconds over baseline (SC-006)
- [ ] T067 Verify quickstart.md examples compile by running them as doc tests or manual verification against the implemented API

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 — BLOCKS all user stories
- **US4 (Phase 3)**: Depends on Phase 2 — provides `SolverTestCase` and `load_test_cases` for US1/US2
- **US2 (Phase 4)**: Depends on Phase 2 + US4 `SolverTestCase` type — provides `NumericalValidator` for US1
- **US1 (Phase 5)**: Depends on US4 (test cases) + US2 (validator) — mock solver uses both
- **US3 (Phase 6)**: Depends on Phase 2 only — generators are independent of US1/US2/US4
- **Integration (Phase 7)**: Depends on US4 + US2 (uses `load_test_cases` and `NumericalValidator`)
- **Polish (Phase 8)**: Depends on all phases complete

### User Story Dependencies

```
Phase 1 (Setup)
    │
Phase 2 (Foundational types)
    │
    ├──────────────────────┐
    │                      │
Phase 3 (US4: Cases)    Phase 6 (US3: Generators) ← independent
    │
Phase 4 (US2: Validator)
    │
Phase 5 (US1: Harness trait)
    │
Phase 7 (Integration refactoring)
    │
Phase 8 (Polish)
```

### Within Each User Story

- Tests MUST be written and FAIL before implementation (constitution Principle III)
- Types/structs before functions that use them
- Core logic before edge case handling
- Story tests PASS before moving to next phase

### Parallel Opportunities

- **Phase 2**: T005, T006, T007, T008, T009 are all in different files/independent — run in parallel
- **Phase 3 (US4)**: T012-T016 (tests) can all run in parallel
- **Phase 4 (US2)**: T022-T027 (tests) can all run in parallel
- **Phase 5 (US1)**: T035-T039 (tests) can all run in parallel
- **Phase 6 (US3)**: T044-T051 (tests) can all run in parallel; T055-T057 (pattern generators) can run in parallel
- **Phase 6 (US3) runs in parallel with Phase 3-5** if US3 only needs Phase 2

---

## Parallel Example: Phase 2 (Foundational)

```
# All these tasks touch different sections and can run in parallel:
T005: TestKind enum in src/testing/harness.rs
T006: MetricResult struct in src/testing/harness.rs
T007: TestResult struct in src/testing/harness.rs
T008: TestMatrixProperties struct in src/testing/cases.rs
T009: SolverTestCase struct in src/testing/cases.rs

# Note: T005-T007 are in the same file but define independent types.
# In practice, write them sequentially within harness.rs, then T008-T009 in cases.rs in parallel.
```

## Parallel Example: US3 (Generators)

```
# Pattern generators are independent of each other:
T055: generate_arrow in src/testing/generators.rs
T056: generate_tridiagonal in src/testing/generators.rs
T057: generate_banded in src/testing/generators.rs

# These are in the same file but implement independent functions.
```

---

## Implementation Strategy

### MVP First (US4 + US2)

1. Complete Phase 1: Setup (T001-T004)
2. Complete Phase 2: Foundational types (T005-T011)
3. Complete Phase 3: US4 — `load_test_cases` with filtering (T012-T021)
4. Complete Phase 4: US2 — `NumericalValidator` (T022-T034)
5. **STOP and VALIDATE**: All 15 hand-constructed matrices load and validate with NumericalValidator
6. This is a useful MVP: future phases can already use `load_test_cases` + `NumericalValidator`

### Full Delivery

7. Complete Phase 5: US1 — `SolverTest` trait + mock (T035-T043)
8. Complete Phase 6: US3 — Generators (T044-T058) — can overlap with Phase 5
9. Complete Phase 7: Integration refactoring (T059-T063)
10. Complete Phase 8: Polish (T064-T067)

### Incremental Value

- After Phase 4: Can validate any factorization with configurable tolerances
- After Phase 5: Can test any `SolverTest` implementor on the full test suite
- After Phase 6: Can generate random matrices for property-based testing
- After Phase 7: Existing tests consolidated, codebase cleaner

---

## Notes

- Constitution Principle III (TDD) is NON-NEGOTIABLE: write tests first, verify they fail, then implement
- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Commit after each task or logical group (constitution Principle VI)
- Stop at any checkpoint to validate story independently
- Existing `validate.rs` functions remain unchanged — `NumericalValidator` wraps them
