# Tasks: APTP Data Structures

**Input**: Design documents from `/specs/009-aptp-data-structures/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/aptp-api.md

**Tests**: TDD is REQUIRED by constitution (Principle III). Tests are written before implementation and must FAIL before the corresponding implementation task.

**Organization**: Tasks are grouped by user story. US1 (MixedDiagonal storage + solve) and US2 (pivot tracking + inertia) share types but are testable independently. US3 (perm_from_forward) is fully independent.

**References**: All implementations must include rustdoc citations per plan.md Traceability Requirement. Primary reference: Hogg, Duff & Lopez (2020) — `references/ssids/duff2020.md`. See plan.md References section for the complete citation table.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup (Module Structure)

**Purpose**: Create the `aptp/` module directory, hub file, and wire it into `lib.rs`.

- [x] T001 Create module directory `src/aptp/` and hub file `src/aptp/mod.rs` with submodule declarations (`pub mod pivot; pub mod diagonal; pub mod inertia; pub mod perm;`) and re-exports of all public types
- [x] T002 Add `pub mod aptp;` to `src/lib.rs` (after the existing `pub mod validate;` line)
- [x] T003 [P] Create empty submodule files: `src/aptp/pivot.rs`, `src/aptp/diagonal.rs`, `src/aptp/inertia.rs`, `src/aptp/perm.rs` — each with a module-level doc comment describing its purpose and citing the primary reference (Hogg, Duff & Lopez 2020)

**Checkpoint**: `cargo check` passes with the new (empty) module structure.

---

## Phase 2: Foundational (Inertia Relocation)

**Purpose**: Relocate `Inertia` from `io/reference.rs` to `aptp/inertia.rs` so it can be used by both reference loading and MixedDiagonal. This MUST complete before US1/US2 implementation.

- [x] T004 Move the `Inertia` struct definition (with all derives, impl block, and doc comments) from `src/io/reference.rs` to `src/aptp/inertia.rs`. Preserve the existing `Serialize`/`Deserialize` derives unchanged.
- [x] T005 Add `pub use crate::aptp::Inertia;` to `src/io/reference.rs` replacing the struct definition. Verify the re-export preserves the `LEntry` and `DBlock` types that remain in `io/reference.rs`.
- [x] T006 Run `cargo test` to confirm all existing tests pass with zero changes to `validate.rs`, `tests/hand_constructed.rs`, or any other file that imports `Inertia` via `crate::io::reference`.

**Checkpoint**: All 83+ existing unit tests and 2 integration tests pass. `Inertia` is importable from both `crate::aptp::Inertia` and `crate::io::reference::Inertia`.

---

## Phase 3: User Story 1 — Mixed 1x1/2x2 Block Diagonal (Priority: P1)

**Goal**: Implement PivotType, Block2x2, and MixedDiagonal with storage, query, and `solve_in_place`.

**Independent Test**: Construct known mixed diagonals, solve D x = b, verify ||D x - b|| / ||b|| < 10^-14.

### Tests for US1 (TDD — write FIRST, verify they FAIL)

- [x] T007 [P] [US1] Write unit tests for `PivotType` in `src/aptp/pivot.rs` (`#[cfg(test)] mod tests`): verify `OneByOne`, `TwoByTwo { partner }`, and `Delayed` variants; test `Copy`, `Clone`, `PartialEq`, `Eq` derives; test Debug formatting
- [x] T008 [P] [US1] Write unit tests for `Block2x2` in `src/aptp/pivot.rs`: verify field storage (`first_col`, `a`, `b`, `c`); test `determinant()` and `trace()` computed properties with known values; test `Copy`, `Clone`, `PartialEq` derives
- [x] T009 [US1] Write unit tests for `MixedDiagonal` construction and query in `src/aptp/diagonal.rs` (`#[cfg(test)] mod tests`): test `new(n)` creates all-Delayed; test `set_1x1` and `set_2x2` mark correct pivot types; test `get_pivot_type`, `get_1x1`, `get_2x2` return correct values; test `num_delayed`, `num_1x1`, `num_2x2_pairs` counts; test `dimension()`
- [x] T010 [US1] Write unit tests for `MixedDiagonal::solve_in_place` in `src/aptp/diagonal.rs`: test all-1x1 diagonal solve (trivial: x[i] = b[i] / d[i]); test all-2x2 diagonal solve (Cramer's rule); test mixed 1x1+2x2 solve; verify relative error < 10^-14 for each case; test dimension-0 solve is a no-op
- [x] T011 [US1] Write edge case tests for `MixedDiagonal` in `src/aptp/diagonal.rs`: test dimension 0; test dimension 1 (single 1x1); test dimension 2 (single 2x2 block); test all-2x2 with even n; test that `solve_in_place` debug-asserts on delayed columns (use `#[should_panic]` in debug mode); test that `set_2x2` at column n-1 with odd n triggers debug-assert (`#[should_panic]`); test MixedDiagonal at n=10,000 with random mixed pivot pattern and verify solve round-trip accuracy (SC-001)

### Implementation for US1

- [x] T012 [P] [US1] Implement `PivotType` enum in `src/aptp/pivot.rs` per contracts/aptp-api.md: three variants (`OneByOne`, `TwoByTwo { partner: usize }`, `Delayed`), derives (`Debug, Clone, Copy, PartialEq, Eq`), rustdoc citing Hogg, Duff & Lopez (2020) for APTP pivot classification
- [x] T013 [P] [US1] Implement `Block2x2` struct in `src/aptp/pivot.rs` per contracts/aptp-api.md: four fields (`first_col`, `a`, `b`, `c`), `determinant()` and `trace()` methods, derives (`Debug, Clone, Copy, PartialEq`), rustdoc citing Bunch & Kaufman (1977) for 2x2 symmetric block structure
- [x] T014 [US1] Implement `MixedDiagonal` struct and construction/query methods in `src/aptp/diagonal.rs` per contracts/aptp-api.md and data-model.md: private fields (`pivot_map`, `diag_1x1`, `blocks_2x2`, `n`); `new(n)`, `set_1x1`, `set_2x2`; all query methods (`dimension`, `get_pivot_type`, `get_1x1`, `get_2x2`, `num_delayed`, `num_1x1`, `num_2x2_pairs`); all debug-asserts per contracts/aptp-api.md debug-assert table; rustdoc citing Hogg, Duff & Lopez (2020) for mixed diagonal D storage in APTP
- [x] T015 [US1] Implement `MixedDiagonal::solve_in_place` in `src/aptp/diagonal.rs` per research.md R3: iterate `pivot_map`; for 1x1: `x[col] /= diag_1x1[col]`; for 2x2: Cramer's rule (`det = a*c - b*b`, `x1 = (c*r1 - b*r2)/det`, `x2 = (a*r2 - b*r1)/det`); skip partner columns (already processed); debug-assert no Delayed columns and no singular pivots (`det != 0` for 2x2, `value != 0` for 1x1); rustdoc citing Cramer's rule for 2x2 systems and Hogg, Duff & Lopez (2020) Section 3 for D-solve in APTP context
- [x] T016 [US1] Run `cargo test` for US1 — verify all T007-T011 tests pass, `cargo clippy` clean, `cargo fmt --check` clean

**Checkpoint**: MixedDiagonal fully functional with construction, query, and solve. All US1 acceptance scenarios verified.

---

## Phase 4: User Story 2 — Pivot Tracking & Inertia Computation (Priority: P1)

**Goal**: Add `compute_inertia` to MixedDiagonal and verify inertia matches known values from hand-constructed reference factorizations.

**Independent Test**: Construct MixedDiagonal from reference factorization data, compute inertia, compare against known Inertia values.

### Tests for US2 (TDD — write FIRST, verify they FAIL)

- [x] T017 [US2] Write unit tests for `MixedDiagonal::compute_inertia` in `src/aptp/diagonal.rs`: test all-positive-1x1 → inertia (n, 0, 0); test mixed-sign-1x1 → correct positive/negative counts; test 2x2 block with det < 0 → one positive, one negative; test 2x2 block with det > 0, trace > 0 → both positive; test 2x2 block with det > 0, trace < 0 → both negative; test mixed 1x1+2x2 combined inertia
- [x] T018 [US2] Write integration test in `tests/aptp_data_structures.rs` that loads all 15 hand-constructed reference factorizations (via `io::registry`), constructs a MixedDiagonal from each reference's `d_blocks` and `permutation`, computes inertia via `compute_inertia()`, and asserts it matches the reference `inertia` field (SC-003). Note: includes a local test helper to convert reference `DBlock` entries into MixedDiagonal `set_1x1`/`set_2x2` calls, since `DBlock` (IO deserialization) and `MixedDiagonal` (live factorization) are intentionally separate types (FR-012)

### Implementation for US2

- [x] T019 [US2] Implement `MixedDiagonal::compute_inertia` in `src/aptp/diagonal.rs` per research.md R2: iterate pivot_map; for 1x1 pivots classify sign of diagonal value; for 2x2 blocks use trace/determinant table (det > 0 + trace > 0 → +2 positive, det > 0 + trace < 0 → +2 negative, det < 0 → +1/+1, det = 0 cases → zero eigenvalue handling); debug-assert no Delayed columns; return `Inertia { positive, negative, zero }`; rustdoc citing eigenvalue sign classification from standard linear algebra and the role of inertia in APTP per Hogg, Duff & Lopez (2020) Section 2
- [x] T020 [US2] Run `cargo test` for US2 — verify all T017-T018 tests pass including integration test against 15 hand-constructed matrices

**Checkpoint**: Inertia computation verified against all 15 hand-constructed reference factorizations. US2 acceptance scenarios verified.

---

## Phase 5: User Story 3 — Permutation Construction from Ordering Output (Priority: P2)

**Goal**: Implement `perm_from_forward` utility that bridges ordering algorithm output to faer's `Perm<usize>`.

**Independent Test**: Construct permutations from forward arrays, verify round-trip identity, test error handling.

### Tests for US3 (TDD — write FIRST, verify they FAIL)

- [x] T021 [P] [US3] Write unit tests for `perm_from_forward` in `src/aptp/perm.rs` (`#[cfg(test)] mod tests`): test identity permutation [0,1,2,...]; test non-trivial permutation with round-trip (apply perm then inverse = original); test composition of two permutations; test empty array (dimension 0); test single element [0]; test error on duplicate indices; test error on out-of-bounds indices; test that returned `Perm` arrays satisfy forward/inverse relationship; test perm_from_forward at n=10,000 and verify round-trip identity (SC-004)
- [x] T022 [P] [US3] Write integration test in `tests/aptp_data_structures.rs` that loads hand-constructed reference factorizations, passes their `permutation` field through `perm_from_forward`, and verifies the resulting `Perm<usize>` has correct forward/inverse arrays (SC-004)

### Implementation for US3

- [x] T023 [US3] Implement `perm_from_forward` in `src/aptp/perm.rs` per data-model.md and research.md R1: validate input via `validate::validate_permutation`; compute inverse array (`inv[fwd[i]] = i`); convert to `Box<[usize]>`; construct via `Perm::new_checked(fwd_box, inv_box, n)`; return `Result<Perm<usize>, SparseError>`; rustdoc documenting the function's purpose as a bridge between ordering algorithms and faer's permutation type
- [x] T024 [US3] Run `cargo test` for US3 — verify all T021-T022 tests pass

**Checkpoint**: `perm_from_forward` functional with validation and error handling. US3 acceptance scenarios verified.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Final validation, documentation quality, and CI readiness.

- [x] T025 Run full `cargo test --all-targets` to verify all tests pass (unit + integration)
- [x] T026 Run `cargo clippy --all-targets -- -D warnings` and fix any warnings
- [x] T027 Run `cargo fmt --check` and fix any formatting issues
- [x] T028 Run `cargo doc --no-deps` with `RUSTDOCFLAGS="-D warnings"` and fix any doc warnings; verify all public items have rustdoc with academic citations per plan.md Traceability Requirement
- [x] T029 Review quickstart.md examples against actual API — update if any signatures changed during implementation
- [x] T030 Update `dev/ssids-log.md` with Phase 2 completion entry documenting what was built, design decisions made, and reference to this feature spec

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 — Inertia relocation BLOCKS US1/US2
- **US1 (Phase 3)**: Depends on Phase 2 — PivotType, Block2x2, MixedDiagonal storage + solve
- **US2 (Phase 4)**: Depends on Phase 3 (T014 specifically) — adds `compute_inertia` to MixedDiagonal
- **US3 (Phase 5)**: Depends on Phase 1 only — fully independent of US1/US2, can run in parallel with Phase 3-4
- **Polish (Phase 6)**: Depends on all user stories being complete

### User Story Dependencies

- **US1 (P1)**: Requires Phase 2. Implements PivotType, Block2x2, MixedDiagonal core.
- **US2 (P1)**: Requires US1 (T014 — MixedDiagonal struct must exist to add `compute_inertia`).
- **US3 (P2)**: Requires Phase 1 only. Independent of US1/US2 — can run fully in parallel.

### Within Each User Story

- Tests MUST be written first and FAIL before implementation (TDD per constitution)
- Type definitions (PivotType, Block2x2) before struct using them (MixedDiagonal)
- Construction/query before operations (solve, inertia)
- Unit tests before integration tests

### Parallel Opportunities

```
Phase 1 (T001-T003): T003 parallel with T001/T002
Phase 2 (T004-T006): Sequential (relocation must be verified)
Phase 3 (T007-T016): T007, T008 parallel (different test groups in same file)
                      T012, T013 parallel (PivotType and Block2x2 in same file, independent)
Phase 4 (T017-T020): Sequential (tests → implementation → verify)
Phase 5 (T021-T024): T021, T022 parallel (unit vs integration tests)
                      ** Phase 5 can run entirely in parallel with Phases 3-4 **
```

---

## Implementation Strategy

### MVP First (US1 Only)

1. Complete Phase 1: Module structure
2. Complete Phase 2: Inertia relocation
3. Complete Phase 3: US1 (MixedDiagonal with solve)
4. **STOP and VALIDATE**: `cargo test`, `cargo clippy`, manual verification of solve accuracy
5. MixedDiagonal is independently useful for downstream Phase 5 (numeric factorization) development

### Incremental Delivery

1. Phase 1 + Phase 2 → Foundation ready
2. Phase 3 (US1) → MixedDiagonal with solve → Validate (MVP)
3. Phase 4 (US2) → Add inertia computation → Validate against 15 reference matrices
4. Phase 5 (US3) → Add perm_from_forward → Validate (can run in parallel with Phase 3-4)
5. Phase 6 → Polish, docs, CI → Feature complete

### Optimal Parallel Strategy

```
                    Phase 1 (Setup)
                         |
                    Phase 2 (Inertia relocation)
                    /                    \
        Phase 3 (US1)              Phase 5 (US3)    ← parallel
              |                          |
        Phase 4 (US2)              (complete)
              |                    /
        Phase 6 (Polish) ←--------
```

---

## Notes

- [P] tasks = different files, no dependencies on incomplete tasks
- [Story] label maps task to specific user story for traceability
- TDD is mandatory per constitution — test tasks exist for each implementation task
- All rustdoc must cite academic references per plan.md Traceability Requirement
- Commit after each task or logical group (Phase completion at minimum)
- Debug-asserts are used for internal invariants; `Result` returns for user-facing errors
- The existing `DBlock` in `io/reference.rs` is NOT modified (FR-012)
