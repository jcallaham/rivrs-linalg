# Tasks: Match-Order Condensation Pipeline

**Input**: Design documents from `/specs/013-match-order-condensation/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, contracts/

**Tests**: Included — constitution principle III (TDD) requires test-first workflow for all solver components.

**Organization**: Tasks grouped by user story. Each story is independently testable after foundational phase.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup

**Purpose**: Create test scaffolding and result type

- [X] T001 Create integration test file `tests/match_order.rs` with imports for `rivrs_sparse::aptp`, `faer::sparse`, and common test helpers; add empty test functions as placeholders for US1/US2/US3 tests
- [X] T002 Define `MatchOrderResult` struct with fields (`ordering`, `scaling`, `matched`, `condensed_dim`, `singletons`, `two_cycles`) and rustdoc in `src/aptp/ordering.rs`; add re-exports (`MatchOrderResult`, `match_order_metis`) to `src/aptp/mod.rs`

---

## Phase 2: Foundational (Internal Building Blocks)

**Purpose**: Implement and test the three internal helper functions that all user stories depend on

**TDD discipline**: Each helper gets tests written first (red), then implementation (green).

### Cycle Decomposition

- [X] T003 Write unit tests for `split_matching_cycles()` in `src/aptp/ordering.rs` `#[cfg(test)]` module: (a) all-singletons matching (all matched, fwd[i]=i) → all partner=-1, condensed_dim=n; (b) pure 2-cycles [0↔1, 2↔3] (all matched) → partner[0]=1, partner[1]=0, condensed_dim=n/2; (c) 3-cycle [0→1→2→0] → split into 2-cycle {0,1} + singleton {2}; (d) 4-cycle [0→1→2→3→0] → split into 2-cycles {0,1} + {2,3}; (e) 5-cycle [0→1→2→3→4→0] → split into 2-cycles {0,1} + {2,3} + singleton {4}; (f) mixed: singletons + 2-cycles + unmatched (partner=-2, using is_matched=false); (g) trivial n=0 and n=1 cases. Verify tests FAIL (no implementation yet).
- [X] T004 Implement `CycleDecomposition` struct and `fn split_matching_cycles(matching_fwd: &[usize], is_matched: &[bool], n: usize) -> CycleDecomposition` in `src/aptp/ordering.rs`. The `is_matched` slice (from `Mc64Result.is_matched`) distinguishes true singletons (matched, fwd[i]==i) from unmatched indices (not matched, partner=-2). Algorithm: for each index, check is_matched first — unmatched indices get partner=-2 and are excluded from condensed graph; for matched indices, walk cycles via forward array, pair consecutive members into 2-cycles, odd-one-out becomes singleton. Build `old_to_new`/`new_to_old` mappings (matched only). Verify T003 tests PASS.

### Condensed Graph Construction

- [X] T005 Write unit tests for `build_condensed_adjacency()` in `src/aptp/ordering.rs` `#[cfg(test)]` module: (a) 4x4 arrow matrix with known 2-cycle → verify condensed dimension is 3 (2 singletons + 1 pair), edges are deduplicated, no self-loops; (b) diagonal matrix → condensed graph has zero edges; (c) fully-connected 4x4 with 2 pairs → verify all edges merged correctly. Verify tests FAIL.
- [X] T006 Implement `fn build_condensed_adjacency(matrix: SymbolicSparseColMatRef<'_, usize>, decomp: &CycleDecomposition) -> Result<(Vec<i32>, Vec<i32>), SparseError>` in `src/aptp/ordering.rs`. Returns `(xadj, adjncy)` in METIS CSR convention. Algorithm: iterate original columns, map row indices via `old_to_new`, use marker array for deduplication, skip self-loops and unmatched rows, build full symmetric CSR suitable for METIS. Verify T005 tests PASS.

### Ordering Expansion

- [X] T007 Write unit tests for `expand_ordering()` in `src/aptp/ordering.rs` `#[cfg(test)]` module: (a) 4-node condensed order [1,0] with 2 pairs → verify both pairs consecutive in output, valid permutation; (b) mixed singletons + pairs → verify singletons at correct positions, pairs consecutive; (c) with unmatched indices → verify unmatched appended at end. Note: `expand_ordering` takes `(condensed_order: &[i32], decomp: &CycleDecomposition, n: usize)` — inverse is computed internally. Verify tests FAIL.
- [X] T008 Implement `fn expand_ordering(condensed_order: &[i32], decomp: &CycleDecomposition, n: usize) -> Perm<usize>` in `src/aptp/ordering.rs`. Algorithm: build inverse of METIS output internally (`inv_order[order[i]] = i`), walk positions 0..condensed_dim, for each condensed node emit original index(es) via `new_to_old` + `partner`, matched pairs get consecutive positions, unmatched appended at end. Verify T007 tests PASS.

**Checkpoint**: All three internal building blocks implemented and unit-tested. Ready for orchestrator.

---

## Phase 3: User Story 1 — Combined Matching-Ordering with Pair Adjacency (P1)

**Goal**: End-to-end `match_order_metis()` pipeline that guarantees 2-cycle pairs are adjacent in the elimination order.

**Independent Test**: Run on any indefinite SuiteSparse matrix and verify every 2-cycle pair occupies consecutive positions in the output ordering.

### Tests (TDD red phase)

- [X] T009 [US1] Write pair adjacency integration tests in `tests/match_order.rs`: (a) hand-constructed 6x6 indefinite matrix with known 2-cycles → assert pairs consecutive in ordering; (b) hand-constructed PD matrix (all singletons) → assert valid permutation, no 2-cycles; (c) edge cases: n=0, n=1, diagonal matrix → assert trivial valid ordering
- [X] T010 [P] [US1] Write fill quality comparison test in `tests/match_order.rs`: for each SuiteSparse CI-subset matrix, run both `match_order_metis()` and `metis_ordering()`, compare `AptpSymbolic::analyze()` predicted nnz(L), assert `condensed_nnz <= unconstrained_nnz * 1.10` (SC-003 one-sided tolerance)
- [X] T011 [P] [US1] Write symbolic analysis integration test in `tests/match_order.rs`: for each SuiteSparse CI-subset matrix, verify `match_order_metis()` result produces valid `AptpSymbolic` via `SymmetricOrdering::Custom(result.ordering.as_ref())` — assert `predicted_nnz() > 0` and no errors (SC-007)

### Implementation

- [X] T012 [US1] Implement `pub fn match_order_metis(matrix: &SparseColMat<usize, f64>) -> Result<MatchOrderResult, SparseError>` in `src/aptp/ordering.rs`. Orchestrates: (1) trivial-case short-circuit for n<=1, (2) call `mc64_matching()`, (3) call `split_matching_cycles()` on matching forward array, (4) call `build_condensed_adjacency()` on matrix symbolic + decomposition, (5) call `METIS_NodeND` on condensed graph (reuse FFI pattern from `metis_ordering()`), (6) call `expand_ordering()`, (7) assemble `MatchOrderResult` with ordering, scaling, matched count, and diagnostics. Full rustdoc with algorithm references per contracts/api.md.
- [X] T013 [US1] Verify all US1 tests pass: run `cargo test --test match_order` and `cargo test ordering` (unit tests). Fix any failures. Then verify existing tests unaffected: `cargo test --test mc64_matching` and `cargo test --test metis_ordering` must still pass (SC-006).

**Checkpoint**: US1 complete — `match_order_metis()` produces valid orderings with guaranteed pair adjacency on hand-constructed and CI-subset matrices.

---

## Phase 4: User Story 2 — Structurally Singular Matrix Handling (P2)

**Goal**: Pipeline correctly handles partial matchings, placing unmatched indices at the end of the ordering.

**Independent Test**: Construct a structurally singular matrix, run pipeline, verify unmatched indices at positions [matched..n).

### Tests and Validation

- [X] T014 [US2] Write structurally singular integration tests in `tests/match_order.rs`: (a) hand-constructed 5x5 matrix with 1 structurally zero row → assert unmatched index at position n-1, matched indices at positions 0..n-1; (b) hand-constructed 6x6 with 2 unmatched → assert last 2 positions are unmatched; (c) verify scaling vector well-defined for all indices (positive, finite) including unmatched (Duff-Pralet correction)
- [X] T015 [US2] Write SuiteSparse singular matrix tests in `tests/match_order.rs` (use `#[ignore]` for full-set): for any SuiteSparse matrix where `mc64_matching().matched < n`, verify unmatched indices appear at ordering positions [matched..n) and output is valid permutation (SC-002)

**Checkpoint**: US2 complete — singular matrices handled gracefully with unmatched-at-end guarantee.

---

## Phase 5: User Story 3 — Reduced METIS Input Size for Performance (P3)

**Goal**: Demonstrate condensation reduces METIS input size and total pipeline overhead is bounded.

**Independent Test**: Time `match_order_metis()` vs separate `mc64_matching() + metis_ordering()` and verify overhead <= 1.5x.

### Benchmarks and Validation

- [X] T016 [US3] Add Criterion benchmark group `match_order` in `benches/solver_benchmarks.rs`: benchmark `match_order_metis()` on 3-5 CI-subset matrices (varying size/indefiniteness), report condensed_dim/n ratio, compare total time vs `mc64_matching() + metis_ordering()` separately (SC-005)
- [X] T017 [US3] Write condensation ratio validation test in `tests/match_order.rs`: for each CI-subset indefinite matrix, assert `result.condensed_dim < n` when `result.two_cycles > 0` (SC-004); log compression ratios for analysis
- [X] T018 [US3] Run full SuiteSparse validation (`cargo test --test match_order -- --ignored --test-threads=1`): verify pair adjacency (SC-001), valid permutation (FR-007), and symbolic analysis validity (SC-007) on all 65+ matrices

**Checkpoint**: US3 complete — performance characteristics validated, benchmarks available.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, plan updates, and final validation

- [X] T019 Add Phase 4.3 section to `docs/ssids-plan.md`: document match-order condensation as a completed deliverable with algorithm reference (SPRAL `match_order.f90`), success criteria results, and lessons learned
- [X] T020 Add Phase 4.3 entry to `docs/ssids-log.md`: summarize what was built (cycle splitting, condensed graph, METIS on condensed, expansion), key decisions (code in ordering.rs, marker-array dedup), and test results
- [X] T021 Run quickstart.md validation checklist: verify all 7 success criteria (SC-001 through SC-007) pass, document results; run `cargo fmt --check` and `cargo clippy` for code quality

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on T002 (MatchOrderResult struct) from Setup
- **US1 (Phase 3)**: Depends on ALL foundational tasks (T003-T008) — BLOCKS on internal helpers
- **US2 (Phase 4)**: Depends on T012 (match_order_metis implementation) — tests exercise the full pipeline
- **US3 (Phase 5)**: Depends on T012 (match_order_metis implementation) — benchmarks and validation
- **Polish (Phase 6)**: Depends on all user stories being complete

### User Story Dependencies

- **US1 (P1)**: Depends on Phase 2 completion. This is the MVP.
- **US2 (P2)**: Depends on US1 (T012). Can start as soon as `match_order_metis()` exists. Singular handling should work from the foundational implementation — these tests validate it.
- **US3 (P3)**: Depends on US1 (T012). Can run in parallel with US2 (different test files/concerns).

### Within Each Phase

- T003→T004, T005→T006, T007→T008: TDD red→green pairs (sequential within pair)
- T003, T005, T007: Test tasks are independent of each other but all in same file
- T009→T012→T013: Tests first, then implementation, then verification
- T010, T011: Can run in parallel with each other (different test concerns)
- T014, T015: Can run in parallel (different test scenarios)
- T016, T017: Can run in parallel (different files: benches/ vs tests/)

### Parallel Opportunities

Within Phase 2 (foundational), the three TDD pairs can proceed sequentially but the test-writing tasks could be authored simultaneously if working in separate functions within the test module.

Within Phase 3 (US1), T010 and T011 are marked [P] — they test different aspects (fill quality vs symbolic validity) and can be written in parallel.

US2 and US3 can proceed in parallel once US1 is complete.

---

## Parallel Example: User Story 1

```bash
# After Phase 2 complete, write US1 tests in parallel:
Task: "T010 [P] Fill quality comparison test in tests/match_order.rs"
Task: "T011 [P] Symbolic analysis integration test in tests/match_order.rs"

# Then implement orchestrator (depends on both test tasks being written):
Task: "T012 Implement match_order_metis() in src/aptp/ordering.rs"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T002)
2. Complete Phase 2: Foundational (T003-T008) — TDD cycles for internal helpers
3. Complete Phase 3: User Story 1 (T009-T013) — full pipeline with pair adjacency
4. **STOP and VALIDATE**: Run `cargo test --test match_order` — all pair adjacency and fill quality tests pass
5. This is a complete, shippable increment

### Incremental Delivery

1. Setup + Foundational → Internal helpers tested and working
2. Add US1 → Core pipeline validated on hand-constructed + CI-subset (MVP)
3. Add US2 → Singular matrix handling validated
4. Add US3 → Performance benchmarked, full SuiteSparse validated
5. Polish → Documentation updated, code quality verified

### Single Developer Strategy (this project)

Execute phases sequentially in priority order:
1. Phase 1 → Phase 2 → Phase 3 (US1) → validate MVP
2. Phase 4 (US2) → Phase 5 (US3) → validate comprehensive
3. Phase 6 (Polish) → commit and PR

---

## Notes

- [P] tasks = different files or independent test functions, no blocking dependencies
- [Story] label maps task to specific user story for traceability
- TDD pairs (T003→T004, T005→T006, T007→T008) must be sequential: write test, verify fail, implement, verify pass
- All code changes are in two files: `src/aptp/ordering.rs` (implementation) and `src/aptp/mod.rs` (re-exports)
- All tests in one integration test file: `tests/match_order.rs`
- Benchmarks in existing file: `benches/solver_benchmarks.rs`
- Commit after each TDD green phase and after each user story checkpoint
- Existing tests (`mc64_matching`, `metis_ordering`) must remain passing throughout (SC-006)
