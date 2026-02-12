# Tasks: METIS Nested Dissection Ordering

**Input**: Design documents from `/specs/011-metis-ordering/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, contracts/api.md

**Tests**: Required by constitution (Principle III: TDD is NON-NEGOTIABLE). Tests written first, verified to fail, then implementation follows.

**Organization**: Tasks grouped by user story. Each story is independently testable after completion.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Add METIS dependency and create module skeleton

- [ ] T001 Add `metis-sys` as a required dependency in Cargo.toml (vendored, no features needed)
- [ ] T002 [P] Add `pub mod ordering;` and `pub use ordering::metis_ordering;` to src/aptp/mod.rs
- [ ] T003 [P] Create src/aptp/ordering.rs with module-level doc comment (algorithm references: Karypis & Kumar 1998, George 1973), imports (`metis_sys`, `faer::perm::Perm`, `faer::sparse::SymbolicSparseColMatRef`, `crate::error::SparseError`), and stub function signatures for `metis_ordering` and `extract_adjacency` (return `todo!()`)
- [ ] T004 Verify `cargo check` succeeds with the new dependency and module structure (confirms vendored METIS compiles)

**Checkpoint**: Project compiles with METIS linked. Module skeleton in place.

---

## Phase 2: Foundational — Adjacency Extraction

**Purpose**: The CSC-to-CSR adjacency extraction is a blocking prerequisite for all user stories. It must be correct before METIS can be called.

**⚠️ CRITICAL**: `extract_adjacency` correctness determines the correctness of all METIS orderings downstream.

- [ ] T005 Write unit tests for `extract_adjacency` in src/aptp/ordering.rs: (a) 3x3 tridiagonal matrix — verify xadj/adjncy match expected CSR, no diagonal entries; (b) 5x5 arrow matrix — verify star-graph adjacency; (c) diagonal-only matrix — verify empty adjncy; (d) upper-triangle-only input — verify symmetric output; (e) verify adjacency invariants: xadj[0]==0, monotonic xadj, no self-loops, symmetry (for each (i,j) also (j,i) present)
- [ ] T006 Implement `extract_adjacency` in src/aptp/ordering.rs: accept `SymbolicSparseColMatRef`, symmetrize the structural pattern, exclude diagonal entries, produce `(Vec<i32>, Vec<i32>)` as (xadj, adjncy) in METIS CSR format with `i32` indices
- [ ] T007 Verify `cargo test` passes for adjacency extraction tests

**Checkpoint**: Adjacency extraction verified on hand-constructed matrices. Ready for METIS integration.

---

## Phase 3: User Story 1 — Compute METIS Ordering (Priority: P1) 🎯 MVP

**Goal**: A solver developer can call `metis_ordering(matrix.symbolic())` and receive a valid `Perm<usize>` that integrates with `AptpSymbolic::analyze()` via `SymmetricOrdering::Custom`.

**Independent Test**: Compute METIS ordering on any SuiteSparse CI-subset matrix, verify the result is a valid permutation, and pass it to `AptpSymbolic::analyze()` successfully.

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T008 [US1] Write unit tests for `metis_ordering` in src/aptp/ordering.rs: (a) 5x5 tridiagonal matrix — returns valid permutation (every index 0..4 appears exactly once in forward and inverse arrays); (b) hand-constructed 10x10 matrix — verify forward/inverse consistency (fwd[inv[i]] == i for all i); (c) dimension-0 matrix — returns empty permutation; (d) dimension-1 matrix — returns trivial permutation; (e) diagonal-only matrix (no off-diagonal) — returns identity permutation; (f) dimension exceeds i32::MAX — returns SparseError
- [ ] T009 [P] [US1] Write integration test in tests/metis_ordering.rs: for each CI-subset SuiteSparse matrix, compute `metis_ordering`, verify permutation validity (correct dimension, all indices present), then pass to `AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Custom(perm.as_ref()))` and verify symbolic analysis succeeds with non-zero predicted nnz

### Implementation for User Story 1

- [ ] T010 [US1] Implement `metis_ordering` in src/aptp/ordering.rs: (a) validate dimension fits in i32; (b) handle trivial cases (dim 0, 1, no off-diagonal → identity perm); (c) call `extract_adjacency`; (d) allocate output arrays, set METIS options (NUMBERING=0, rest defaults); (e) call `metis_sys::METIS_NodeND` in unsafe block; (f) check return code, map errors to SparseError; (g) convert i32 arrays to usize; (h) construct `Perm::new_checked(iperm_box, perm_box, n)` per research.md Decision 5 (METIS iperm = faer forward, METIS perm = faer inverse)
- [ ] T011 [US1] Add comprehensive rustdoc to `metis_ordering` in src/aptp/ordering.rs: module-level doc with algorithm references, function-level doc with `# Arguments`, `# Returns`, `# Errors`, `# Algorithm References` (Karypis & Kumar 1998, George 1973), `# Examples` section showing usage with `AptpSymbolic::analyze`
- [ ] T012 [US1] Run `cargo test` and verify all unit tests (T008) and integration tests (T009) pass

**Checkpoint**: `metis_ordering` works on all CI-subset matrices. Valid permutations produced. Integration with `AptpSymbolic::analyze()` verified. User Story 1 is independently functional.

---

## Phase 4: User Story 2 — Reproduce Published Fill Predictions (Priority: P2)

**Goal**: Predicted nnz(L) using METIS ordering matches Hogg et al. (2016) Table III values within 20%, and METIS produces less fill than AMD on >= 80% of test matrices.

**Independent Test**: Run symbolic analysis with METIS on CI-subset matrices, compare predicted nnz(L) against paper-reported values and against AMD.

### Tests for User Story 2

- [ ] T013 [US2] Add Hogg et al. (2016) Table III reference nnz(L) values as named constants in tests/metis_ordering.rs (for matrices overlapping with our SuiteSparse collection: bmwcra_1, crankseg_2, hood, ldoor, etc. — extract from `/workspace/rivrs-linalg/references/ssids/hogg2016.md`)
- [ ] T014 [US2] Write test `test_metis_nnz_matches_paper_values` in tests/metis_ordering.rs: for each Table III matrix in the CI-subset, compute METIS ordering, run `AptpSymbolic::analyze`, and assert `predicted_nnz()` is within 20% of the paper-reported value
- [ ] T015 [US2] Write test `test_metis_reduces_fill_vs_amd` in tests/metis_ordering.rs: for each CI-subset matrix, compute both AMD and METIS symbolic analysis, count how many matrices have METIS fill <= AMD fill, assert count >= 80% of total

### Implementation for User Story 2

> No new production code needed — US2 is purely validation of US1's output quality.

- [ ] T016 [US2] Run `cargo test` and verify T014 and T015 pass. If nnz(L) values are outside 20% tolerance, investigate: (a) check if METIS v5 vs v4 ordering difference is the cause; (b) adjust tolerance or document discrepancy in test comments

**Checkpoint**: METIS ordering quality validated against published benchmarks. Fill reduction vs AMD confirmed.

---

## Phase 5: User Story 3 — Full SuiteSparse Without AMD Limitations (Priority: P3)

**Goal**: All 67 SuiteSparse matrices complete symbolic analysis with METIS ordering, without the MAX_DIM_FOR_AMD guard, in under 2 minutes total.

**Independent Test**: Run `cargo test --test symbolic_analysis_full -- --ignored --test-threads=1` with METIS ordering and verify all matrices pass.

### Tests and Implementation for User Story 3

- [ ] T017 [US3] Add METIS ordering test function(s) to tests/symbolic_analysis_full.rs: parallel structure to existing AMD tests but using `metis_ordering` + `SymmetricOrdering::Custom` instead of `SymmetricOrdering::Amd`. No dimension cap. Include wall-clock timing assertion (total < 120 seconds).
- [ ] T018 [US3] Update existing AMD test functions in tests/symbolic_analysis_full.rs: make `MAX_DIM_FOR_AMD` guard apply only to AMD-specific tests (not to METIS tests). Add comments noting AMD limitation and referencing METIS as the preferred ordering.
- [ ] T019 [US3] Run `cargo test --test symbolic_analysis_full -- --ignored --test-threads=1` and verify all METIS ordering tests pass on the full 67-matrix collection

**Checkpoint**: Full SuiteSparse collection tested with METIS ordering. No dimension guard needed. All 67 matrices pass.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, code quality, and final validation

- [ ] T020 [P] Run `cargo fmt --check` and `cargo clippy` on all modified/new files, fix any warnings
- [ ] T021 [P] Update docs/ssids-log.md with Phase 4.1 development entry: what was built (METIS ordering via metis-sys), key decisions (required dep, vendored, idx_t conversion), test results (fill comparison stats)
- [ ] T022 [P] Update docs/ssids-plan.md: mark Phase 4.1 as COMPLETE, add any lessons learned
- [ ] T023 Run full test suite: `cargo test` (all unit + integration) and `cargo test -- --ignored --test-threads=1` (full SuiteSparse). Verify zero failures.
- [ ] T024 Validate quickstart.md example compiles and runs correctly (manual verification)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on Setup (T001-T004) — adjacency extraction blocks all METIS calls
- **US1 (Phase 3)**: Depends on Foundational (T005-T007) — core ordering implementation
- **US2 (Phase 4)**: Depends on US1 (T008-T012) — validation of ordering quality
- **US3 (Phase 5)**: Depends on US1 (T008-T012) — full-suite testing with METIS
- **Polish (Phase 6)**: Depends on US1-US3 completion

### User Story Dependencies

```
Phase 1 (Setup) → Phase 2 (Foundational: extract_adjacency)
                        ↓
                   Phase 3 (US1: metis_ordering)
                    ↙            ↘
      Phase 4 (US2: fill quality)  Phase 5 (US3: full suite)
                    ↘            ↙
                   Phase 6 (Polish)
```

- **US2 and US3 can proceed in parallel** after US1 is complete (different test files, no code dependencies)

### Within Each User Story

1. Tests written FIRST, verified to FAIL (constitution Principle III)
2. Implementation follows
3. Tests verified to PASS
4. Story complete before moving to next priority

### Parallel Opportunities

**Phase 1**: T002 and T003 can run in parallel (different files: mod.rs vs ordering.rs)
**Phase 3**: T008 and T009 can run in parallel (different files: src/ vs tests/)
**Phase 4-5**: US2 and US3 can proceed in parallel after US1 completes
**Phase 6**: T020, T021, T022 can run in parallel (different files)

---

## Parallel Example: User Story 1

```
# After Foundational phase completes, launch tests in parallel:
Agent 1: T008 — Unit tests for metis_ordering in src/aptp/ordering.rs
Agent 2: T009 — Integration tests in tests/metis_ordering.rs

# Then implement sequentially:
T010 — Implement metis_ordering in src/aptp/ordering.rs
T011 — Add rustdoc documentation in src/aptp/ordering.rs
T012 — Run tests and verify all pass
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T004) — ~30 min
2. Complete Phase 2: Foundational (T005-T007) — ~1 hour
3. Complete Phase 3: User Story 1 (T008-T012) — ~2-3 hours
4. **STOP and VALIDATE**: `metis_ordering` works on CI-subset matrices
5. Commit and verify CI passes

### Incremental Delivery

1. Setup + Foundational → METIS compiles and adjacency extraction works
2. User Story 1 → Core ordering function works → **MVP complete**
3. User Story 2 → Fill quality validated against literature → **Confidence in correctness**
4. User Story 3 → Full test suite unblocked → **Complete feature**
5. Polish → Documentation and code quality → **Ready for merge**

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Constitution requires TDD: tests MUST fail before implementation begins
- Commit after each phase checkpoint
- METIS `iperm` = faer forward array, METIS `perm` = faer inverse array (research.md Decision 5)
- `extract_adjacency` must handle upper-only, lower-only, and full symmetric storage
- Trivial cases (dim 0, 1, diagonal) handled in Rust without calling METIS FFI
