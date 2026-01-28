# Tasks: Sylvester Equation Solver

**Feature**: 001-sylvester-solver
**Input**: Design documents from `/specs/001-sylvester-solver/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/sylvester.yaml

**Tests**: Not explicitly requested in specification. Implementation focuses on correctness verification through examples.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

## Path Conventions

Single Rust library project at repository root:
- Core implementation: `src/`
- Tests: `tests/`
- Examples: `examples/`
- Benchmarks: `benches/`

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project initialization and basic structure

- [ ] T001 Create Rust library crate structure in faer-rs/csrrs/ with Cargo.toml
- [ ] T002 Add faer dependency (>= 0.19) to Cargo.toml
- [ ] T003 [P] Add ndarray dependency (>= 0.16) to Cargo.toml for Python interop
- [ ] T004 [P] Configure Rust edition 2021 and minimum supported version 1.75+ in Cargo.toml
- [ ] T005 [P] Create src/lib.rs with module structure (error, sylvester, utils)
- [ ] T006 [P] Create basic README.md with project overview and clean room notice
- [ ] T007 [P] Setup examples/ directory for usage demonstrations

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core infrastructure that MUST be complete before ANY user story can be implemented

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

- [ ] T008 Define SylvesterError enum in src/error.rs with all error variants from data-model.md
- [ ] T009 Implement Display trait for SylvesterError in src/error.rs with informative messages
- [ ] T010 Implement Error trait for SylvesterError in src/error.rs
- [ ] T011 Define SylvesterSolution struct in src/sylvester/types.rs with generics over T: RealField
- [ ] T012 [P] Define Sign enum (Plus, Minus) in src/sylvester/types.rs
- [ ] T013 [P] Define Transpose enum (NoTrans, Trans) in src/sylvester/types.rs
- [ ] T014 [P] Define EquationType enum (Continuous, Discrete) in src/sylvester/types.rs
- [ ] T015 Create src/sylvester/validation.rs with validate_dimensions function
- [ ] T016 [P] Implement validate_finite function in src/sylvester/validation.rs
- [ ] T017 [P] Implement validate_quasi_triangular function in src/sylvester/validation.rs
- [ ] T018 Create src/sylvester/mod.rs and re-export public types and functions
- [ ] T019 Create tests/common/mod.rs with test matrix generation helpers
- [ ] T020 [P] Create tests/common/fixtures.rs with benchmark problem matrices from G&VL

**Checkpoint**: Foundation ready - user story implementation can now begin in parallel

---

## Phase 3: User Story 1 - Solve Standard Sylvester Equations (Priority: P1) 🎯 MVP

**Goal**: Implement continuous-time Sylvester solver (AX + XB = C) with Bartels-Stewart algorithm

**Independent Test**: Verify using Bartels-Stewart benchmark problems from G&VL Chapter 7, checking ||AX + XB - C|| < 10⁻¹²

### Implementation for User Story 1

**Step 1: Triangular Sylvester Solver (Core Algorithm)**

- [ ] T021 [US1] Implement 1×1 block solver in src/sylvester/triangular.rs with overflow prevention
- [ ] T022 [US1] Implement 2×2 block solver (lasy2-style) in src/sylvester/triangular.rs
- [ ] T023 [US1] Implement SCALE factor accumulation logic in src/sylvester/triangular.rs
- [ ] T024 [US1] Implement back-substitution loop for case NoTrans/NoTrans in src/sylvester/triangular.rs
- [ ] T025 [US1] Implement back-substitution loop for case Trans/NoTrans in src/sylvester/triangular.rs
- [ ] T026 [US1] Implement back-substitution loop for case Trans/Trans in src/sylvester/triangular.rs
- [ ] T027 [US1] Implement back-substitution loop for case NoTrans/Trans in src/sylvester/triangular.rs
- [ ] T028 [US1] Create solve_triangular_sylvester public function in src/sylvester/triangular.rs
- [ ] T029 [US1] Add near-singular detection (SMIN threshold checks) in src/sylvester/triangular.rs

**Step 2: High-Level Continuous Solver (Bartels-Stewart Algorithm)**

- [ ] T030 [US1] Implement compute_residual helper function in src/sylvester/utils.rs
- [ ] T031 [US1] Create solve_continuous skeleton in src/sylvester/continuous.rs
- [ ] T032 [US1] Add input validation calls in solve_continuous
- [ ] T033 [US1] Integrate faer Schur decomposition for matrix A in solve_continuous
- [ ] T034 [US1] Integrate faer Schur decomposition for matrix B in solve_continuous
- [ ] T035 [US1] Transform RHS: F = U₁ᵀ C U₂ using faer matmul in solve_continuous
- [ ] T036 [US1] Call solve_triangular_sylvester for TY + YS = F in solve_continuous
- [ ] T037 [US1] Back-transform solution: X = U₁ Y U₂ᵀ using faer matmul in solve_continuous
- [ ] T038 [US1] Compute residual norm ||AX + XB - C|| in solve_continuous
- [ ] T039 [US1] Construct and return SylvesterSolution struct in solve_continuous
- [ ] T040 [US1] Setup workspace allocation using dyn_stack in solve_continuous

**Step 3: Schur Form Solver (Advanced API)**

- [ ] T041 [US1] Implement solve_continuous_schur in src/sylvester/continuous.rs
- [ ] T042 [US1] Add quasi-triangular validation for schur_a and schur_b in solve_continuous_schur
- [ ] T043 [US1] Add orthogonality check for U and V matrices in solve_continuous_schur

**Step 4: Validation and Examples**

- [ ] T044 [US1] Create example basic_continuous.rs demonstrating solve_continuous
- [ ] T045 [US1] Create integration test test_continuous_small.rs with 2×2, 3×3 analytical cases
- [ ] T046 [US1] Create integration test test_continuous_gvl.rs with Golub & Van Loan examples
- [ ] T047 [US1] Create integration test test_continuous_benchmarks.rs comparing to SLICOT data
- [ ] T048 [US1] Add edge case test for dimension mismatches in tests/test_errors.rs
- [ ] T049 [US1] Add edge case test for NaN/Inf inputs in tests/test_errors.rs
- [ ] T050 [US1] Add edge case test for empty matrices (0×0) in tests/test_edge_cases.rs
- [ ] T051 [US1] Add rustdoc comments to all public functions with academic citations

**Checkpoint**: At this point, User Story 1 should be fully functional - continuous-time solver working with numerical validation

---

## Phase 4: User Story 2 - Handle Discrete-Time Sylvester Equations (Priority: P2)

**Goal**: Implement discrete-time Sylvester solver (AXB + X = C) using modified back-substitution

**Independent Test**: Verify using SLICOT SB04QD benchmark problems, checking ||AXB + X - C|| < 10⁻¹²

### Implementation for User Story 2

**Step 1: Discrete-Time Triangular Solver**

- [ ] T052 [US2] Implement discrete 1×1 block solver in src/sylvester/triangular_discrete.rs
- [ ] T053 [US2] Implement discrete 2×2 block solver in src/sylvester/triangular_discrete.rs
- [ ] T054 [US2] Implement discrete back-substitution loop (Y + SYT = F form) in src/sylvester/triangular_discrete.rs
- [ ] T055 [US2] Add SCALE factor logic for discrete case in src/sylvester/triangular_discrete.rs
- [ ] T056 [US2] Create solve_triangular_sylvester_discrete public function in src/sylvester/triangular_discrete.rs

**Step 2: High-Level Discrete Solver**

- [ ] T057 [US2] Create solve_discrete skeleton in src/sylvester/discrete.rs
- [ ] T058 [US2] Add Sign parameter handling (Plus vs Minus) in solve_discrete
- [ ] T059 [US2] Add input validation calls in solve_discrete
- [ ] T060 [US2] Integrate faer Schur decomposition for A and B in solve_discrete
- [ ] T061 [US2] Transform RHS: F = UᵀCV using faer matmul in solve_discrete
- [ ] T062 [US2] Call solve_triangular_sylvester_discrete for Y + SYT = F in solve_discrete
- [ ] T063 [US2] Back-transform solution: X = UYVᵀ in solve_discrete
- [ ] T064 [US2] Compute discrete residual norm ||AXB ± X - C|| in solve_discrete
- [ ] T065 [US2] Construct and return SylvesterSolution struct in solve_discrete

**Step 3: Validation and Examples**

- [ ] T066 [US2] Create example basic_discrete.rs demonstrating solve_discrete
- [ ] T067 [US2] Create integration test test_discrete_small.rs with small analytical cases
- [ ] T068 [US2] Create integration test test_discrete_slicot.rs comparing to SLICOT SB04QD data
- [ ] T069 [US2] Add test for Sign::Plus variant in tests/test_discrete_sign.rs
- [ ] T070 [US2] Add test for Sign::Minus variant in tests/test_discrete_sign.rs
- [ ] T071 [US2] Add rustdoc comments to discrete solver functions with academic citations

**Checkpoint**: At this point, User Stories 1 AND 2 should both work independently - continuous and discrete solvers complete

---

## Phase 5: User Story 3 - Detect and Report Singular Cases (Priority: P3)

**Goal**: Implement eigenvalue separation estimation and clear diagnostic error reporting

**Independent Test**: Construct matrices with known common eigenvalues, verify detection and informative error messages

### Implementation for User Story 3

**Step 1: Condition Estimation**

- [ ] T072 [P] [US3] Implement sep(A,B) estimation for continuous case in src/sylvester/condition.rs
- [ ] T073 [P] [US3] Implement sep(A,B) estimation for discrete case in src/sylvester/condition.rs
- [ ] T074 [US3] Create estimate_separation function in src/sylvester/condition.rs
- [ ] T075 [US3] Define separation threshold constants (e.g., 1e-14) in src/sylvester/condition.rs

**Step 2: Integration with Solvers**

- [ ] T076 [US3] Add eigenvalue separation check in solve_continuous before triangular solve
- [ ] T077 [US3] Add eigenvalue separation check in solve_discrete before triangular solve
- [ ] T078 [US3] Enhance near_singular flag setting based on separation estimates
- [ ] T079 [US3] Return CommonEigenvalues error when separation below threshold

**Step 3: Error Message Enhancement**

- [ ] T080 [US3] Add condition number to error messages in SylvesterError Display impl
- [ ] T081 [US3] Add suggested remediation (preconditioning, regularization) to error messages
- [ ] T082 [US3] Create helper function for eigenvalue diagnostics in src/sylvester/diagnostics.rs

**Step 4: Validation and Examples**

- [ ] T083 [US3] Create example singular_case.rs demonstrating error handling
- [ ] T084 [US3] Create test test_common_eigenvalues.rs with matrices sharing eigenvalues
- [ ] T085 [US3] Create test test_near_singular.rs with ill-conditioned matrices
- [ ] T086 [US3] Add test verifying informative error messages in tests/test_error_messages.rs
- [ ] T087 [US3] Add rustdoc comments documenting conditioning and error scenarios

**Checkpoint**: All user stories 1-3 should now work independently with robust error detection

---

## Phase 6: User Story 4 - Solve Large-Scale Problems Efficiently (Priority: P4)

**Goal**: Optimize performance for matrices n, m >= 500 using blocked Level-3 BLAS algorithm

**Independent Test**: Benchmark against SLICOT using 500×500 and 1000×1000 matrices, verify execution time and memory

### Implementation for User Story 4

**Step 1: Benchmark Infrastructure**

- [ ] T088 [P] [US4] Create benchmarks/sylvester_benchmarks.rs with criterion setup
- [ ] T089 [P] [US4] Add baseline benchmarks for current unblocked algorithm in benchmarks/
- [ ] T090 [P] [US4] Generate random test matrices for sizes 100, 200, 500, 1000 in benchmarks/
- [ ] T091 [P] [US4] Document baseline performance numbers in benchmarks/README.md

**Step 2: Blocked Triangular Solver (dtrsyl3-style)**

- [ ] T092 [US4] Implement blocked 1×1 solver with Level-3 BLAS in src/sylvester/triangular_blocked.rs
- [ ] T093 [US4] Implement blocked 2×2 solver with Level-3 BLAS in src/sylvester/triangular_blocked.rs
- [ ] T094 [US4] Implement blocked back-substitution with panel updates in src/sylvester/triangular_blocked.rs
- [ ] T095 [US4] Add block size selection heuristic (e.g., 64×64) in src/sylvester/triangular_blocked.rs
- [ ] T096 [US4] Optimize workspace allocation for blocked variant in src/sylvester/triangular_blocked.rs
- [ ] T097 [US4] Use faer matmul for off-diagonal block updates in src/sylvester/triangular_blocked.rs

**Step 3: Integration and Switchover Logic**

- [ ] T098 [US4] Add runtime switchover logic (n > threshold → use blocked) in solve_continuous
- [ ] T099 [US4] Add runtime switchover logic for discrete solver in solve_discrete
- [ ] T100 [US4] Tune switchover threshold via benchmarks (likely n ~ 100-200)
- [ ] T101 [US4] Add parallelism support via faer Par parameter in blocked solver

**Step 4: Performance Validation**

- [ ] T102 [US4] Add large-scale benchmarks (500×500, 1000×1000) in benchmarks/
- [ ] T103 [US4] Verify blocked algorithm gives identical results to unblocked in tests/
- [ ] T104 [US4] Profile memory usage using valgrind or similar tools
- [ ] T105 [US4] Compare performance to SLICOT and document results in benchmarks/PERFORMANCE.md
- [ ] T106 [US4] Verify execution time < 5s for 500×500 on standard hardware
- [ ] T107 [US4] Verify peak memory < 100MB beyond input for 1000×1000

**Checkpoint**: All user stories should now be complete with production-grade performance

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, cleanup, and final validation

- [ ] T108 [P] Create comprehensive API documentation with examples in docs/api.md
- [ ] T109 [P] Verify all public functions have rustdoc comments with citations
- [ ] T110 [P] Create migration guide from MATLAB/SLICOT in docs/migration.md
- [ ] T111 [P] Add crate-level documentation in src/lib.rs with overview and quickstart
- [ ] T112 Run cargo clippy and fix all warnings
- [ ] T113 Run cargo fmt to ensure consistent formatting
- [ ] T114 Verify code coverage reaches 95% target (SC-007) using cargo-tarpaulin
- [ ] T115 [P] Add CI/CD configuration (.github/workflows/) for automated testing
- [ ] T116 Verify all academic references cited in rustdoc match research.md
- [ ] T117 Run quickstart.md examples manually to validate documentation accuracy
- [ ] T118 Create CHANGELOG.md documenting initial release features
- [ ] T119 [P] Add LICENSE files (MIT/Apache-2.0 dual licensing)
- [ ] T120 Final integration test: Run all examples and verify outputs

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - can start immediately
- **Foundational (Phase 2)**: Depends on Setup completion - BLOCKS all user stories
- **User Stories (Phase 3-6)**: All depend on Foundational phase completion
  - User Story 1 (P1): Can start after Foundational - No dependencies on other stories
  - User Story 2 (P2): Can start after Foundational - Independent of US1 (separate algorithm)
  - User Story 3 (P3): Depends on US1 and US2 (adds error detection to existing solvers)
  - User Story 4 (P4): Depends on US1 and US2 (optimizes existing algorithms)
- **Polish (Phase 7)**: Depends on all desired user stories being complete

### User Story Dependencies

```
Foundational (Phase 2) [BLOCKS EVERYTHING]
    ↓
    ├─→ US1 (Continuous Solver) ──┐
    │                              ├─→ US3 (Error Detection)
    └─→ US2 (Discrete Solver) ────┘        ↓
                ↓                          US4 (Performance)
                └──────────────────────────┘
```

- **User Story 1 (P1)**: Can start immediately after Foundational
- **User Story 2 (P2)**: Can start immediately after Foundational (parallel with US1)
- **User Story 3 (P3)**: Requires US1 and US2 complete (adds diagnostics to both)
- **User Story 4 (P4)**: Requires US1 and US2 complete (optimizes both), can proceed in parallel with US3

### Within Each User Story

- Error types and validation → Core algorithm → High-level API → Tests → Documentation
- Triangular solver → High-level solver → Examples → Integration tests
- Models (types) → Services (algorithms) → Public API → Validation

### Parallel Opportunities

**Phase 1 (Setup)**:
- T003, T004, T006, T007 can all run in parallel (different files)

**Phase 2 (Foundational)**:
- T012, T013, T014 (type definitions) can run in parallel
- T016, T017 (validation functions) can run in parallel after T015
- T019, T020 (test infrastructure) can run in parallel

**Phase 3 (User Story 1)**:
- T044-T050 (examples and tests) can run in parallel after T043 (implementation complete)

**Phase 4 (User Story 2)**:
- T066-T071 (examples and tests) can run in parallel after T065 (implementation complete)
- **US2 can proceed entirely in parallel with US1** (separate files, no shared code)

**Phase 5 (User Story 3)**:
- T072, T073 (condition estimation) can run in parallel

**Phase 6 (User Story 4)**:
- T088-T091 (benchmarking setup) can run in parallel
- T102-T107 (performance validation) can run in parallel after T101

**Phase 7 (Polish)**:
- T108, T109, T110, T111, T115, T119 can all run in parallel (documentation tasks)

---

## Parallel Example: User Story 1

```bash
# After T020 completes, launch examples and tests in parallel:
Task T044: "Create example basic_continuous.rs"
Task T045: "Create integration test test_continuous_small.rs"
Task T046: "Create integration test test_continuous_gvl.rs"
Task T047: "Create integration test test_continuous_benchmarks.rs"
Task T048: "Add edge case test for dimension mismatches"
Task T049: "Add edge case test for NaN/Inf inputs"
Task T050: "Add edge case test for empty matrices"
Task T051: "Add rustdoc comments"

# All 8 tasks can proceed in parallel (different test files)
```

## Parallel Example: User Stories 1 & 2

```bash
# After Foundational (Phase 2) completes:

# Team Member A: Implements US1 (T021-T051)
Task T021-T051: User Story 1 (Continuous Solver)

# Team Member B: Implements US2 (T052-T071) in parallel
Task T052-T071: User Story 2 (Discrete Solver)

# These are completely independent implementations
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T007)
2. Complete Phase 2: Foundational (T008-T020) - CRITICAL BLOCKER
3. Complete Phase 3: User Story 1 (T021-T051)
4. **STOP and VALIDATE**:
   - Run examples/basic_continuous.rs
   - Verify all integration tests pass
   - Check residual norms < 10⁻¹²
   - Test with G&VL benchmark problems
5. **MVP COMPLETE**: Continuous-time solver ready for use

### Incremental Delivery

1. **Foundation** (T001-T020) → Project structure ready
2. **MVP: Continuous Solver** (T021-T051) → Can solve AX + XB = C
3. **Discrete Solver** (T052-T071) → Can solve AXB ± X = C
4. **Error Detection** (T072-T087) → Robust diagnostics
5. **Performance** (T088-T107) → Production-grade speed
6. **Polish** (T108-T120) → Documentation and CI/CD

Each increment adds value without breaking previous functionality.

### Parallel Team Strategy

With multiple developers:

1. **Week 1**: All developers work on Setup + Foundational together (T001-T020)
2. **Week 2-3**:
   - Developer A: US1 Continuous Solver (T021-T051)
   - Developer B: US2 Discrete Solver (T052-T071)
   - Both work independently in parallel
3. **Week 4**:
   - Developer A: US3 Error Detection (T072-T087)
   - Developer B: US4 Performance (T088-T107)
   - Can proceed in parallel after both US1 and US2 complete
4. **Week 5**: All developers work on Polish together (T108-T120)

---

## Notes

- **[P] marker**: Tasks marked [P] can run in parallel (different files, no conflicts)
- **[Story] label**: Maps task to specific user story from spec.md for traceability
- **File paths**: All tasks include exact file paths for implementation
- **Clean room compliance**: Never consult SLICOT source code (.f files), only documentation and test data
- **Academic citations**: Must cite Bartels-Stewart, G&VL, LAPACK sources in all rustdoc comments
- **Checkpoints**: Each user story phase ends with validation checkpoint
- **Testing**: Validation through examples and benchmark comparisons (no explicit TDD requested)
- **Commit strategy**: Commit after each task or logical group completion
- **Coverage target**: Aim for 95% code coverage per SC-007

---

## Summary

**Total Tasks**: 120
**Task Breakdown**:
- Setup (Phase 1): 7 tasks
- Foundational (Phase 2): 13 tasks (BLOCKS all user stories)
- User Story 1 (Phase 3): 31 tasks - Continuous solver
- User Story 2 (Phase 4): 20 tasks - Discrete solver
- User Story 3 (Phase 5): 16 tasks - Error detection
- User Story 4 (Phase 6): 20 tasks - Performance optimization
- Polish (Phase 7): 13 tasks

**Parallel Opportunities**: 35 tasks marked [P] can run in parallel
**MVP Scope**: Phases 1-3 (51 tasks) deliver continuous-time solver
**Independent Stories**: US1 and US2 can be developed in parallel after Foundational phase

**Critical Path**:
1. Foundational phase MUST complete first (13 tasks)
2. US1 + US2 can proceed in parallel (51 combined tasks)
3. US3 requires US1+US2 (16 tasks)
4. US4 requires US1+US2 (20 tasks), can overlap with US3
5. Polish requires all stories (13 tasks)

**Estimated Delivery**:
- MVP (Continuous solver): ~40-50% of total effort
- Full feature (all 4 stories): 100% of effort
- With 2 developers: ~60-70% time savings via parallelization
