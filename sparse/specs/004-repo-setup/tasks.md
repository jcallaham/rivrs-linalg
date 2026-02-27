# Tasks: Repository Setup for Solver Development

**Input**: Design documents from `/specs/004-repo-setup/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, contracts/

**Tests**: Included — the project constitution (Principle III) mandates TDD for all solver components.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Dependencies, error types, and module wiring that all user stories depend on

- [x] T001 Add `serde` and `serde_json` dependencies to `Cargo.toml` (serde = { version = "1", features = ["derive"] }, serde_json = "1")
- [x] T002 Extend `SparseError` in `src/error.rs` with IO variants: `IoError { source: String, path: String }` and `ParseError { reason: String, path: String, line: Option<usize> }`. Add `From<std::io::Error>` and `From<serde_json::Error>` impls. Existing variants unchanged.
- [x] T003 Create module structure: `src/io.rs` (pub mod mtx, reference, registry), `src/io/` directory with empty `mtx.rs`, `reference.rs`, `registry.rs`. Create `src/validate.rs` (empty). Update `src/lib.rs` to add `pub mod io; pub mod validate;`

**Checkpoint**: `cargo build` compiles with new module skeleton and dependencies

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Data types shared across multiple user stories — must be complete before story implementation

**CRITICAL**: No user story work can begin until this phase is complete

- [x] T004 [P] Define `Inertia` struct (positive, negative, zero: usize) with `#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]` in `src/io/reference.rs`. Include `Inertia::dimension()` method returning sum. Export from `src/io.rs`.
- [x] T005 [P] Define `LEntry` struct (row, col: usize, value: f64) with serde derives in `src/io/reference.rs`
- [x] T006 [P] Define `DBlock` enum with `OneByOne { value: f64 }` and `TwoByTwo { values: [[f64; 2]; 2] }` variants in `src/io/reference.rs`. Implement custom serde deserializer that reads the `"size"` field to determine variant (size=1 → scalar values, size=2 → 2x2 matrix values). See `test-data/hand-constructed/stress-delayed-pivots.json` for 2x2 example.
- [x] T007 [P] Define `ReferenceFactorization` struct (matrix_name: String, permutation: Vec<usize>, l_entries: Vec<LEntry>, d_blocks: Vec<DBlock>, inertia: Inertia, notes: String) with serde derives in `src/io/reference.rs`
- [x] T008 [P] Define `MatrixProperties` struct (symmetric, positive_definite, indefinite: bool, difficulty: String, structure: Option<String>) and `MatrixMetadata` struct (name, source, category, path: String, size, nnz: usize, in_repo: bool, properties: MatrixProperties, factorization_path: Option<String>, plus other optional fields using `#[serde(default)]`) with serde derives in `src/io/registry.rs`. Schema must match `test-data/metadata.json` structure.
- [x] T009 [P] Define `TestMatrix` struct (metadata: MatrixMetadata, matrix: SparseColMat<usize, f64>, reference: Option<ReferenceFactorization>) in `src/io/registry.rs`

**Checkpoint**: `cargo build` compiles. All shared data types are defined and serializable.

---

## Phase 3: User Story 1 — Load Test Matrices in Rust Tests (Priority: P1) MVP

**Goal**: Developers can load .mtx and .json files into faer types and verify matrix properties in tests

**Independent Test**: `cargo test --test hand_constructed` passes, loading all 15 matrices with correct dimensions

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [x] T010 [P] [US1] Write unit tests for MTX parser in `src/io/mtx.rs` (as `#[cfg(test)] mod tests`): test loading `arrow-5-pd.mtx` produces 5x5 matrix with 9 stored entries (lower tri), test that malformed input returns `SparseError::ParseError`, test that unsupported format header returns error. Use `CARGO_MANIFEST_DIR` to locate test-data/.
- [x] T011 [P] [US1] Write unit tests for JSON reference loader in `src/io/reference.rs` (as `#[cfg(test)] mod tests`): test loading `arrow-5-pd.json` produces correct inertia (5,0,0) and 10 l_entries, test loading `stress-delayed-pivots.json` produces 5 `TwoByTwo` d_blocks, test invalid JSON returns error.
- [x] T012 [P] [US1] Write unit tests for registry in `src/io/registry.rs` (as `#[cfg(test)] mod tests`): test `load_registry()` returns 82 entries, test `load_test_matrix("arrow-5-pd")` returns Some with correct dimensions, test `load_test_matrix("nonexistent")` returns error, test `load_test_matrix` for a gitignored matrix returns Ok(None) when file absent.
- [x] T013 [US1] Write integration test `tests/hand_constructed.rs`: iterate all 15 hand-constructed matrices from registry (filter by source == "hand-constructed"), load each, assert dimensions and nnz match metadata, assert reference factorization loads for each. Verify all 15 load in < 1 second (use `std::time::Instant`).

### Implementation for User Story 1

- [x] T014 [US1] Implement `load_mtx(path: &Path) -> Result<SparseColMat<usize, f64>, SparseError>` in `src/io/mtx.rs`: parse `%%MatrixMarket matrix coordinate real symmetric` header, skip `%` comment lines, parse `nrows ncols nnz` size line, parse 1-indexed triplets converting to 0-indexed, symmetrize off-diagonal entries (add (j,i,v) for each (i,j,v) where i != j), construct via `SparseColMat::try_new_from_triplets(nrows, ncols, &triplets)`. Return descriptive `ParseError` on any failure with file path and line number.
- [x] T015 [US1] Implement `load_reference(path: &Path) -> Result<ReferenceFactorization, SparseError>` in `src/io/reference.rs`: read file to string, deserialize with serde_json, validate l_entries have col < row, validate permutation is a valid permutation of 0..n (contains each index exactly once). Return `ParseError` with file path on failure.
- [x] T016 [US1] Implement registry functions in `src/io/registry.rs`: `load_registry()` reads `test-data/metadata.json` relative to `CARGO_MANIFEST_DIR`, deserializes into `Vec<MatrixMetadata>`. `load_test_matrix(name)` finds entry by name, resolves .mtx path relative to test-data dir, returns Ok(None) if .mtx file doesn't exist, loads .mtx via `load_mtx`, loads .json via `load_reference` if `factorization_path` is Some, constructs `TestMatrix`.
- [x] T017 [US1] Run `cargo test` — verify all US1 unit tests (T010-T012) and integration test (T013) pass. Fix any issues. Run `cargo clippy --all-targets -- -D warnings` and `cargo fmt --check`.

**Checkpoint**: All 15 hand-constructed matrices load with correct dimensions. `cargo test --test hand_constructed` passes. This is the MVP — all subsequent stories build on this.

---

## Phase 4: User Story 2 — CI Pipeline Validates Sparse Code (Priority: P2)

**Goal**: CI automatically runs test/lint/fmt/doc checks for the sparse domain on every PR

**Independent Test**: Push a commit, verify sparse CI jobs appear and pass on GitHub Actions

### Implementation for User Story 2

- [x] T018 [US2] Add `test-sparse` job to `.github/workflows/ci.yml`: mirror `test-control` structure, run `cd sparse && cargo test --all-targets`, test matrix with `[stable, "1.87"]` for MSRV validation. Ensure checkout@v4 and rust-cache are used.
- [x] T019 [P] [US2] Add `lint-sparse` job to `.github/workflows/ci.yml`: mirror `lint-control` structure, run `cd sparse && cargo fmt --check` and `cd sparse && cargo clippy --all-targets -- -D warnings`.
- [x] T020 [P] [US2] Add `doc-sparse` job to `.github/workflows/ci.yml`: mirror `doc-control` structure, run `cd sparse && cargo doc --no-deps` with `RUSTDOCFLAGS: -D warnings`.

**Checkpoint**: CI YAML is valid. All sparse jobs mirror the existing control/ pattern. Locally verify: `cargo test`, `cargo clippy --all-targets -- -D warnings`, `cargo fmt --check`, `cargo doc --no-deps` all pass.

---

## Phase 5: User Story 3 — Numerical Validation Utilities (Priority: P3)

**Goal**: Working reconstruction error, backward error, and inertia comparison functions tested against known factorizations

**Independent Test**: `cargo test --test hand_constructed` passes with reconstruction error < 10^-12 for all matrices

### Tests for User Story 3

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [x] T021 [P] [US3] Write unit tests for `reconstruction_error` in `src/validate.rs` (as `#[cfg(test)] mod tests`): test with `arrow-5-pd` (load matrix + reference, compute error, assert < 1e-12), test with `stress-delayed-pivots` (2x2 d_blocks, assert < 1e-12), test with intentionally perturbed L entries (assert error > 0.01).
- [x] T022 [P] [US3] Write unit tests for `backward_error` in `src/validate.rs`: construct a small known system Ax=b (e.g., use arrow-5-pd matrix, set b = A*x_exact for a known x_exact, verify backward_error(A, x_exact, b) < 1e-14), test with perturbed x that backward error is measurably larger.
- [x] T023 [P] [US3] Write unit tests for `check_inertia` in `src/validate.rs`: test matching inertia returns true, test mismatched inertia returns false.

### Implementation for User Story 3

- [x] T024 [US3] Implement `reconstruction_error(a, reference)` in `src/validate.rs`: (1) get n from a.nrows(), (2) build dense L as n×n identity then fill strict lower triangle from reference.l_entries, (3) build dense D as n×n zero then fill diagonal blocks from reference.d_blocks (1x1 → single diagonal entry, 2x2 → 2x2 block on diagonal), (4) convert sparse A to dense via `.to_dense()`, (5) apply permutation P to get P^T*A*P by reindexing rows and cols of A_dense using reference.permutation, (6) compute LDLt = L * D * L.transpose() via faer dense multiply, (7) compute diff = pap - ldlt, (8) return diff.norm_l2() / a_dense.norm_l2(). Handle edge case: if ||A|| == 0, return 0.0.
- [x] T025 [US3] Implement `backward_error(a, x, b)` in `src/validate.rs`: (1) compute r = A*x - b using faer's sparse_dense_matmul for A*x then dense subtraction, (2) compute norms: ||r||, ||A||_F (via a.to_dense().norm_l2()), ||x|| (as Col norm), ||b|| (as Col norm), (3) return ||r|| / (||A||*||x|| + ||b||). Handle edge case: if denominator == 0, return 0.0.
- [x] T026 [US3] Implement `check_inertia(computed, expected)` in `src/validate.rs`: field-wise equality comparison of positive, negative, and zero counts.
- [x] T027 [US3] Update `tests/hand_constructed.rs` to add reconstruction error and inertia validation: for each hand-constructed matrix with a reference factorization, compute `reconstruction_error` and assert < 1e-12, compute `check_inertia` comparing reference inertia against itself (sanity check — real solver inertia testing comes in later phases).
- [x] T028 [US3] Run `cargo test` — verify all US3 unit tests (T021-T023) and updated integration test (T027) pass. SC-008 validated: reconstruction error < 10^-12 for all hand-constructed matrices.

**Checkpoint**: All validation functions work end-to-end. `reconstruction_error` returns < 10^-12 for all 15 hand-constructed matrices with known factorizations.

---

## Phase 6: User Story 4 — Benchmark Scaffold (Priority: P4)

**Goal**: `cargo bench` runs at least one benchmark loading a test matrix

**Independent Test**: `cargo bench` completes without error and produces timing output

### Implementation for User Story 4

- [x] T029 [US4] Create `benches/matrix_loading.rs` with criterion benchmark group: benchmark `load_mtx` for `test-data/hand-constructed/arrow-10-indef.mtx`, benchmark `load_mtx` + `reconstruction_error` for same matrix (requires loading reference too). Add `[[bench]] name = "matrix_loading" harness = false` to `Cargo.toml` if not already present.
- [x] T030 [US4] Run `cargo bench` — verify benchmarks execute and produce timing output. Verify no clippy warnings.

**Checkpoint**: `cargo bench` produces timing output for matrix loading.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Final validation, documentation, and cleanup

- [x] T031 Write integration test `tests/suitesparse_ci.rs`: load all 10 CI-subset matrices from registry (filter by path containing "suitesparse-ci"), verify dimensions and nnz match metadata, verify matrix is symmetric (row indices within bounds, all columns have valid entries). This validates the MTX parser against real-world SuiteSparse format.
- [x] T032 [P] Ensure all public functions in `src/io/mtx.rs`, `src/io/reference.rs`, `src/io/registry.rs`, and `src/validate.rs` have rustdoc comments with `# Errors` sections. Run `cargo doc --no-deps` with `RUSTDOCFLAGS="-D warnings"` and fix any warnings.
- [x] T033 [P] Run full validation suite: `cargo test --all-targets`, `cargo clippy --all-targets -- -D warnings`, `cargo fmt --check`, `cargo doc --no-deps`. All must pass with zero warnings. This validates SC-001 through SC-004.
- [x] T034 Update `dev/ssids-log.md` with Phase 0.4 completion entry: what was built (IO modules, validation utilities, CI pipeline, benchmark scaffold), decisions made (custom MTX parser, faer-based validation, three-tier test data loading).
- [x] T035 Update `dev/ssids-plan.md` to mark Phase 0.4 success criteria as checked and note completion status.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 (T001-T003) — BLOCKS all user stories
- **US1 (Phase 3)**: Depends on Phase 2 — BLOCKS US2 (CI needs passing tests), US3 (validation needs loaders), US4 (benchmarks need loaders)
- **US2 (Phase 4)**: Depends on US1 completion (CI must have tests to run)
- **US3 (Phase 5)**: Depends on US1 completion (validation functions use loaded matrices)
- **US4 (Phase 6)**: Depends on US1 and US3 (benchmarks load matrices and run validation)
- **Polish (Phase 7)**: Depends on all user stories complete

### User Story Dependencies

- **US1 (P1)**: Foundation only → MVP milestone
- **US2 (P2)**: US1 must be complete (CI runs the tests US1 created)
- **US3 (P3)**: US1 must be complete (validation uses IO modules from US1)
- **US4 (P4)**: US1 + US3 must be complete (benchmarks use both loading and validation)

Note: US2 (CI) and US3 (validation) can run in parallel after US1 completes.

### Within Each User Story

- Tests MUST be written and FAIL before implementation (Constitution Principle III)
- Data types before loaders
- Loaders before integration tests
- Core implementation before integration
- Run full check (`test + clippy + fmt`) before marking story complete

### Parallel Opportunities

Within Phase 2: T004-T009 are all [P] (different structs in different files or same file independent sections)
Within US1 tests: T010-T012 are all [P] (different test files/modules)
Within US3 tests: T021-T023 are all [P] (different test functions)
After US1: US2 and US3 can proceed in parallel
Within US2: T019 and T020 are [P] (independent CI job definitions)

---

## Parallel Example: Phase 2 (Foundational)

```
# Launch all struct definitions in parallel:
Task T004: "Define Inertia struct in src/io/reference.rs"
Task T005: "Define LEntry struct in src/io/reference.rs"
Task T006: "Define DBlock enum in src/io/reference.rs"
Task T007: "Define ReferenceFactorization struct in src/io/reference.rs"
Task T008: "Define MatrixMetadata/Properties in src/io/registry.rs"
Task T009: "Define TestMatrix struct in src/io/registry.rs"
```

## Parallel Example: After US1 Completion

```
# US2 and US3 can run in parallel:
Developer A: Phase 4 (US2 — CI pipeline)
Developer B: Phase 5 (US3 — Validation utilities)
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T003)
2. Complete Phase 2: Foundational (T004-T009)
3. Complete Phase 3: US1 — Load Test Matrices (T010-T017)
4. **STOP and VALIDATE**: `cargo test --test hand_constructed` passes, all 15 matrices load correctly
5. This is a deliverable milestone — all subsequent work builds on working IO

### Incremental Delivery

1. Setup + Foundational → Skeleton compiles
2. US1 → Matrix loading works → MVP!
3. US2 + US3 (parallel) → CI protection + validation utilities
4. US4 → Benchmarks ready
5. Polish → Docs, final checks

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Constitution Principle III requires TDD: write tests first, verify they fail, then implement
- Commit after each task or logical group (Constitution Principle VI)
- All matrices use `coordinate real symmetric` Matrix Market format
- The `d_blocks` JSON field has polymorphic deserialization (size=1 vs size=2) — see T006
- `CARGO_MANIFEST_DIR` env var is used for locating test-data/ in both tests and benchmarks
