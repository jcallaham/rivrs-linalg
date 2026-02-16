# Tasks: Triangular Solve & Solver API

**Input**: Design documents from `/specs/016-triangular-solve-api/`
**Prerequisites**: plan.md, spec.md, data-model.md, contracts/solve-api.md, contracts/internal-api.md, research.md, quickstart.md

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story. Within each story, **tests are written first** (RED phase) using function stubs with `todo!()` bodies, then implementation follows (GREEN phase), per Constitution Principle III (TDD).

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup

**Purpose**: Add error variants, module declarations, create function stubs for TDD workflow.

- [ ] T001 Add `SolveBeforeFactor` variant to `SparseError` in `src/error.rs`
- [ ] T002 [P] Add `pub(crate) mod solve;` and `pub(crate) mod solver;` declarations to `src/aptp/mod.rs`
- [ ] T003 [P] Create stub files `src/aptp/solve.rs` and `src/aptp/solver.rs` with module-level documentation and all function/type signatures from contracts (bodies as `todo!()`). Add re-exports to `src/aptp/mod.rs` for public types (`SparseLDLT`, `AnalyzeOptions`, `FactorOptions`, `SolverOptions`, `OrderingStrategy`). Stubs must compile so test tasks can reference them.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Modifications to existing code that MUST be complete before any solve implementation can work.

**CRITICAL**: No user story work can begin until this phase is complete.

- [ ] T004 Modify `MixedDiagonal::solve_in_place()` in `src/aptp/diagonal.rs` to handle zero pivots gracefully: set `x[col] = 0.0` for 1x1 zero pivots, set both components to 0.0 for 2x2 zero-determinant blocks (remove `debug_assert` panics per SPRAL `action=true` convention)
- [ ] T005 [P] Add `sparse_backward_error(a, x, b) -> f64` function to `src/validate.rs`. Since input matrices are lower-triangle CSC, compute the full symmetric product via two passes: `ax = A_lower*x + A_lower^T*x - diag(A)*x` using `sparse_dense_matmul` from `faer::sparse::linalg::matmul`. Compute `||A||_F` directly from sparse values with 2x multiplier for off-diagonal entries. See research.md R4 for details.
- [ ] T006 Add `scaling: Option<&[f64]>` parameter to `AptpNumeric::factor()` in `src/aptp/numeric.rs`. Apply scaling during `scatter_original_entries` as `scaled_val = scaling[perm_inv[orig_row]] * val * scaling[perm_inv[orig_col]]`. APTP kernel and extend-add remain unchanged.
- [ ] T007 Update all existing callers of `AptpNumeric::factor()` (in `tests/multifrontal.rs`, `tests/hand_constructed.rs`, and any other call sites) to pass `None` for the new `scaling` parameter

**Checkpoint**: Foundation ready — solve implementation can now begin.

---

## Phase 3: User Story 3 — Per-Supernode Triangular Solve Correctness (Priority: P1)

**Goal**: Implement the internal per-supernode forward/backward/diagonal solve functions and the top-level `aptp_solve()` that traverses the assembly tree.

**Independent Test**: For small hand-constructed matrices with known supernode structure, verify that each per-supernode forward solve step produces the analytically expected local vector and that the full forward → D → backward pipeline yields the correct solution in permuted coordinates.

### Tests for User Story 3 (write FIRST — RED)

> **TDD**: Write these tests against the function stubs from T003. Verify they compile but FAIL (`todo!()` panics) before implementing.

- [ ] T008 [US3] Add per-supernode unit tests in `tests/solve.rs`: test gather/scatter correctness on hand-constructed 2-supernode matrices with known L11, D11, L21, col_indices, row_indices. Test `forward_solve_supernode`, `diagonal_solve_supernode`, and `backward_solve_supernode` independently with analytically computed expected values.
- [ ] T009 [US3] Add `aptp_solve` integration test in `tests/solve.rs`: factor a hand-constructed matrix via `AptpNumeric::factor()`, then call `aptp_solve()` directly (in permuted coordinates), verify the solution matches `x_exact` after permutation.

### Implementation for User Story 3 (make tests pass — GREEN)

- [ ] T010 [US3] Implement `forward_solve_supernode(ff, rhs, work)` in `src/aptp/solve.rs`: gather via `col_indices`, solve `L11 * y = local` using `solve_unit_lower_triangular_in_place_with_conj`, write back, scatter via `L21 * local` using `matmul_with_conj` to `row_indices` entries (per contracts/internal-api.md)
- [ ] T011 [US3] Implement `diagonal_solve_supernode(ff, rhs, work)` in `src/aptp/solve.rs`: gather via `col_indices`, call `d11.solve_in_place(&mut local)`, write back
- [ ] T012 [US3] Implement `backward_solve_supernode(ff, rhs, work)` in `src/aptp/solve.rs`: gather local and tmp via `col_indices`/`row_indices`, update `local -= L21^T * tmp` using `matmul_with_conj`, solve `L11^T * z = local` via `solve_unit_upper_triangular_in_place_with_conj` with `l11.transpose()`, write back
- [ ] T013 [US3] Implement `aptp_solve(symbolic, numeric, rhs, stack)` in `src/aptp/solve.rs`: validate dimensions, iterate supernodes in postorder for forward solve, any order for D solve, reverse postorder for backward solve. Use `temp_mat_zeroed` from `MemStack` for per-supernode workspace.
- [ ] T014 [US3] Implement `aptp_solve_scratch(numeric, rhs_ncols)` in `src/aptp/solve.rs`: return `StackReq` based on `numeric.stats().max_front_size` using `temp_mat_scratch::<f64>(max_front_size, rhs_ncols)`

**Checkpoint**: Per-supernode solve is unit-tested and the internal solve pipeline works in permuted coordinates.

---

## Phase 4: User Story 1 — End-to-End Sparse Solve (Priority: P1)

**Goal**: Implement the `SparseLDLT` user-facing solver struct with analyze → factor → solve pipeline and one-shot convenience method. Deliver a working end-to-end solver validated by backward error.

**Independent Test**: Construct `b = A * x_exact`, solve via `SparseLDLT::solve_full()`, verify backward error < 10^-10 on all hand-constructed and SuiteSparse CI matrices.

**Dependencies**: Requires Phase 3 (US3) complete for `aptp_solve()`.

### Tests for User Story 1 (write FIRST — RED)

> **TDD**: Write these tests against the type/method stubs from T003. Verify they compile but FAIL before implementing.

- [ ] T015 [US1] Add end-to-end backward error tests in `tests/solve.rs`: test `SparseLDLT::solve_full()` on all 15 hand-constructed matrices, verify backward error < 10^-10 using `sparse_backward_error()`. Also verify `solver.inertia()` matches analytically known inertia for each test matrix (SC-004).
- [ ] T016 [US1] Add SuiteSparse CI subset backward error tests in `tests/solve.rs`: test `SparseLDLT` analyze → factor → solve on the 10 CI matrices (`test-data/suitesparse-ci/`), verify backward error < 10^-10
- [ ] T017 [US1] Add error handling and edge case tests in `tests/solve.rs`: (a) verify `SolveBeforeFactor` error when solve called before factor, (b) `DimensionMismatch` for wrong RHS length, (c) rank-deficient solve produces solution with zeroed components, (d) 0x0 matrix trivial solve, (e) matrix producing simplicial (single-column) supernodes, (f) API equivalence: for at least one matrix, assert `solve_full()` and `analyze→factor→solve` produce identical results (SC-009, SC-010).

### Implementation for User Story 1 (make tests pass — GREEN)

- [ ] T018 [US1] Implement `AnalyzeOptions`, `FactorOptions`, `SolverOptions`, `OrderingStrategy` types with `Default` impls in `src/aptp/solver.rs` (per contracts/solve-api.md)
- [ ] T019 [US1] Implement `SparseLDLT` struct (fields: `symbolic: AptpSymbolic`, `numeric: Option<AptpNumeric>`, `scaling: Option<Vec<f64>>`) and `SparseLDLT::analyze()` in `src/aptp/solver.rs`: compute ordering via `OrderingStrategy`, call `AptpSymbolic::new()`, return `SparseLDLT` with `numeric = None`
- [ ] T020 [US1] Implement `SparseLDLT::factor()` in `src/aptp/solver.rs`: call `AptpNumeric::factor(symbolic, matrix, options, scaling.as_deref())`, store result in `self.numeric`
- [ ] T021 [US1] Implement `SparseLDLT::solve_in_place()` in `src/aptp/solver.rs`: check factor state, permute RHS (`rhs_perm[perm_fwd[i]] = rhs[i]`), optionally scale (`rhs_perm[i] *= scaling[i]`), call `aptp_solve()`, optionally unscale, unpermute
- [ ] T022 [US1] Implement `SparseLDLT::solve()` (allocating wrapper), `solve_scratch()`, `solve_full()` (one-shot), `n()`, `inertia()`, `stats()` in `src/aptp/solver.rs`

**Checkpoint**: End-to-end solver works on hand-constructed and SuiteSparse CI matrices with backward error < 10^-10.

---

## Phase 5: User Story 2 — Three-Phase API with Reuse (Priority: P1)

**Goal**: Implement `refactor()` and validate that symbolic analysis is reusable across factorizations with different numeric values.

**Independent Test**: Analyze a matrix, factor+solve, then refactor with different values and solve again — both solutions meet backward error tolerance.

**Dependencies**: Requires Phase 4 (US1) complete for `SparseLDLT`.

### Tests for User Story 2 (write FIRST — RED)

- [ ] T023 [US2] Add refactoring test in `tests/solve.rs`: analyze a matrix, factor+solve, then refactor with different numeric values (same sparsity) and solve again, verify both solutions have backward error < 10^-10
- [ ] T024 [US2] Add multiple-RHS reuse test in `tests/solve.rs`: factor a matrix once, solve with 3 different RHS vectors reusing the same factorization, verify all solutions are correct

### Implementation for User Story 2 (make tests pass — GREEN)

- [ ] T025 [US2] Implement `SparseLDLT::refactor()` in `src/aptp/solver.rs` (delegates to `factor()` — provided for API clarity per contracts/solve-api.md)

**Checkpoint**: Three-phase API with reuse is validated. All P1 user stories are complete.

---

## Phase 6: User Story 4 — Scaling Integration (Priority: P2)

**Goal**: When MC64 matching/scaling is requested during analysis, correctly transform scaling factors to elimination order and apply them during factor and solve.

**Independent Test**: For matrices analyzed with `MatchOrderMetis`, verify backward error < 10^-10 against the original (unscaled) matrix.

**Dependencies**: Requires Phase 4 (US1) complete. Requires Phase 2 (T006) for scaling parameter in `factor()`.

### Tests for User Story 4 (write FIRST — RED)

- [ ] T026 [US4] Add MC64 scaling tests in `tests/solve.rs`: test `SparseLDLT` with `MatchOrderMetis` ordering on matrices where MC64 improves conditioning, verify backward error < 10^-10 against original matrix. Compare results with and without scaling.

### Implementation for User Story 4 (make tests pass — GREEN)

- [ ] T027 [US4] Implement scaling coordinate transform in `SparseLDLT::analyze()` in `src/aptp/solver.rs`: when `MatchOrderMetis` is selected, call `match_order_metis()`, permute scaling from original to elimination order (`elim_scaling[i] = orig_scaling[perm_fwd[i]]`), store as `self.scaling`. Verify that `solve_in_place()` (from T021) correctly applies `rhs_perm[i] *= scaling[i]` before forward solve and after backward solve (symmetric scaling in elimination order).

**Checkpoint**: Scaling integration works end-to-end. All P2 US4 is complete.

---

## Phase 7: User Story 5 — Workspace-Efficient Solve (Priority: P2)

**Goal**: Validate that solve uses stack-based workspace (no heap allocation in solve hot path) and that workspace queries return correct sizes.

**Independent Test**: Verify `solve_scratch()` returns valid `StackReq`, and solve succeeds with a `MemStack` of exactly that size.

**Dependencies**: Requires Phase 4 (US1) complete (solve_scratch and MemStack usage already built into aptp_solve).

### Tests for User Story 5

- [ ] T028 [US5] Add workspace correctness test in `tests/solve.rs`: verify `solve_scratch(1)` returns a `StackReq` that is sufficient for solve on all hand-constructed matrices. Verify solve succeeds with `MemBuffer::new(solve_scratch(1))`. Note: no-heap guarantee (SC-008) is verified by code inspection — `aptp_solve` uses only `MemStack`-based allocation.
- [ ] T029 [US5] Add workspace reuse test in `tests/solve.rs`: allocate `MemBuffer` once, create `MemStack`, call `solve_in_place()` twice with different RHS vectors reusing the same stack — verify both solutions correct.

**Checkpoint**: Workspace-efficient solve validated. All user stories complete.

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Full validation suite, documentation, and code quality.

- [ ] T030 [P] Add full SuiteSparse backward error test suite in `tests/solve_suitesparse.rs` (`#[ignore]`): test `SparseLDLT` on all 67 SuiteSparse matrices, verify >95% achieve backward error < 10^-10. Report per-matrix backward error and failure summary.
- [ ] T031 [P] Add rustdoc documentation to all public types and functions in `src/aptp/solver.rs` and `src/aptp/solve.rs`: algorithm references (Duff+2020, Liu 1992), SPRAL equivalents, complexity, usage examples
- [ ] T032 Run `cargo fmt --check` and `cargo clippy` on the full crate, fix any warnings
- [ ] T033 Run `cargo test` (all non-ignored tests pass) and `cargo test -- --ignored --test-threads=1` (full SuiteSparse suite)
- [ ] T034 Validate quickstart.md examples compile and run correctly (compile-test the code snippets from `specs/016-triangular-solve-api/quickstart.md`)
- [ ] T035 Update `docs/ssids-log.md` with Phase 7 completion entry and update `docs/ssids-plan.md` Phase 7 status

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 (T001 for error variant). T004, T005, T006 can run in parallel after T001.
- **US3 (Phase 3)**: Depends on Phase 2 (T004 for zero pivot handling, T006 for factor signature). Tests (T008-T009) written first against stubs from T003, implementation (T010-T014) follows.
- **US1 (Phase 4)**: Depends on Phase 3 (T013 for aptp_solve) and Phase 2 (T005 for sparse_backward_error). Tests (T015-T017) written first, implementation (T018-T022) follows.
- **US2 (Phase 5)**: Depends on Phase 4 (SparseLDLT implemented). Tests (T023-T024) first, implementation (T025) follows.
- **US4 (Phase 6)**: Depends on Phase 4 (SparseLDLT) and Phase 2 (T006 for scaling in factor). Test (T026) first, implementation (T027) follows.
- **US5 (Phase 7)**: Depends on Phase 4 (solve_scratch already implemented). Tests only (T028-T029).
- **Polish (Phase 8)**: Depends on all user story phases complete.

### User Story Dependencies

- **US3 (P1)**: First to implement — internal solve core. No dependency on other stories.
- **US1 (P1)**: Depends on US3 (uses aptp_solve). No dependency on US2/US4/US5.
- **US2 (P1)**: Depends on US1 (uses SparseLDLT). Independent of US4/US5.
- **US4 (P2)**: Depends on US1 (uses SparseLDLT). Independent of US2/US5.
- **US5 (P2)**: Depends on US1 (uses SparseLDLT). Independent of US2/US4.

### Within Each User Story (TDD Workflow)

1. **RED**: Write tests against function/type stubs (compile but fail via `todo!()`)
2. **GREEN**: Implement functions to make tests pass
3. **REFACTOR**: Clean up while tests remain green
4. Unit tests before integration tests within each phase

### Parallel Opportunities

- Phase 1: T002 and T003 can run in parallel
- Phase 2: T004, T005, T006 can run in parallel (different files)
- Phase 3: T010, T011, T012 can be developed in parallel (independent per-supernode helpers in same file)
- Phase 5/6/7: US2, US4, US5 can run in parallel after US1 is complete (different test concerns)
- Phase 8: T030, T031 can run in parallel

---

## Parallel Example: Phase 2 (Foundational)

```bash
# These three tasks modify different files and can run in parallel:
Task T004: "Modify MixedDiagonal::solve_in_place() in src/aptp/diagonal.rs"
Task T005: "Add sparse_backward_error() in src/validate.rs"
Task T006: "Add scaling parameter to factor() in src/aptp/numeric.rs"
```

## Parallel Example: After US1 Complete

```bash
# US2, US4, US5 can start in parallel:
Task T023: "Refactoring test in tests/solve.rs"       # US2 (RED)
Task T026: "MC64 scaling tests in tests/solve.rs"      # US4 (RED)
Task T028: "Workspace correctness test"                # US5
```

---

## Implementation Strategy

### MVP First (US3 + US1)

1. Complete Phase 1: Setup (T001-T003) — stubs enable TDD
2. Complete Phase 2: Foundational (T004-T007) — **CRITICAL: blocks all stories**
3. Complete Phase 3: US3 — Per-supernode solve (T008-T014, tests first)
4. Complete Phase 4: US1 — End-to-end solve (T015-T022, tests first)
5. **STOP and VALIDATE**: Run `cargo test` — all hand-constructed and SuiteSparse CI backward error tests pass

### Incremental Delivery

1. Setup + Foundational → Foundation ready
2. Add US3 → Internal solve validated → Core milestone
3. Add US1 → End-to-end solver works → **MVP! Working solver milestone**
4. Add US2 → Reuse API validated → Three-phase API complete
5. Add US4 → Scaling works → Hard indefinite problems supported
6. Add US5 → Workspace efficiency validated → Performance-ready
7. Polish → Full SuiteSparse validation, documentation → Phase 7 complete

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- **TDD enforced**: Each user story phase has a "Tests (RED)" section before "Implementation (GREEN)" section, per Constitution Principle III.
- US5 (workspace) is architecturally integral to aptp_solve — the MemStack design is built from T013/T014 onward, not bolted on later. Phase 7 tasks validate the workspace contract explicitly.
- US4 (scaling) depends on T006 (foundational scaling parameter) being complete but is otherwise independent of US2/US5.
- The `aptp_solve` function works in permuted coordinates; `SparseLDLT` wraps with permute/scale/unscale/unpermute.
- `col_indices` already encodes `local_perm` (absorbed during Phase 6 extract_front_factors), so no runtime permutation is needed in the per-supernode solve helpers.
- Test files: `tests/solve.rs` for unit + integration tests; `tests/solve_suitesparse.rs` for full SuiteSparse suite (#[ignore]).
