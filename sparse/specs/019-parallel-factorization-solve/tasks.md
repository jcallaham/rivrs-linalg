# Tasks: Parallel Factorization & Solve (Phase 8.2)

**Input**: Design documents from `/specs/019-parallel-factorization-solve/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/api-changes.md

**Tests**: Tests are included per the constitution's TDD requirement (Principle III). New parallel correctness tests validate each user story independently.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing. US1 (intra-node parallelism) depends on foundational API changes. US2 (tree-level) depends on US1 for the API surface but adds its own factorization loop refactoring. US3 (parallel solve) depends on US1's API surface. US4 (benchmarking) depends on US1+US2+US3.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3, US4)
- Include exact file paths in descriptions

---

## Phase 1: Setup (Rayon Dependency)

**Purpose**: Add Rayon dependency and verify the build compiles cleanly.

- [X] T001 Add `rayon = "1"` to `Cargo.toml` dependencies section and verify `cargo build` succeeds

---

## Phase 2: Foundational (API Surface — Par Threading)

**Purpose**: Add `par: Par` to all option structs and propagate through the call chain. No behavioral change — `Par::Seq` default everywhere. All 358 existing tests must pass unchanged.

**CRITICAL**: No user story work can begin until this phase is complete.

- [X] T002 Add `par: Par` field (default `Par::Seq`) to `FactorOptions` and `SolverOptions` in `src/aptp/solver.rs`, update `Default` impls
- [X] T003 Add `par: Par` field (default `Par::Seq`) to `AptpOptions` in `src/aptp/factor.rs`, update `Default` impl
- [X] T004 Thread `par` from `FactorOptions` into `AptpOptions` in `SparseLDLT::factor()` and `SparseLDLT::factor_with_options()` in `src/aptp/solver.rs`
- [X] T005 Add `par: Par` parameter to `aptp_solve()` and `SparseLDLT::solve_in_place()` in `src/aptp/solve.rs` and `src/aptp/solver.rs`; thread `par` to per-supernode solve functions `forward_solve_supernode()`, `backward_solve_supernode()` in `src/aptp/solve.rs`
- [X] T006 Update `SparseLDLT::solve_full()` in `src/aptp/solver.rs` to propagate `par` from `SolverOptions` to `solve_in_place()`
- [X] T007 Update all call sites across tests, examples, and benchmarks to use new signatures: `src/aptp/solve.rs` tests, `src/aptp/numeric.rs` tests, `src/aptp/solver.rs` tests, all `examples/*.rs`, `benches/solver_benchmarks.rs`
- [X] T008 Define `INTRA_NODE_THRESHOLD` constant (value: 256) in `src/aptp/numeric.rs`
- [X] T009 Run `cargo test` (all 358 tests pass), `cargo clippy --all-targets`, `cargo clippy --all-targets --features diagnostic`, and `cargo bench --no-run`

**Checkpoint**: API surface complete. All existing tests pass with Par::Seq defaults. No behavioral change.

---

## Phase 3: User Story 1 — Parallel Factorization of Large Frontal Matrices (Priority: P1) MVP

**Goal**: Replace hardcoded `Par::Seq` in the dense APTP kernel's BLAS-3 operations (TRSM, GEMM) with `options.par`, gated by the 256 front-size threshold. Large IntraNode-class fronts benefit from parallel dense BLAS.

**Independent Test**: Factor `sparsine` (max_front=11,125) with `Par::rayon(4)` and verify backward error < 5e-11 and wall-clock speedup > 1x.

### Tests for User Story 1

- [X] T010 [US1] Add parallel correctness test: factor + solve CI subset matrices with `Par::rayon(4)`, assert backward error < 5e-11 for each, in `src/aptp/solver.rs` tests (new `#[test] fn test_parallel_factor_ci_subset`)
- [X] T011 [US1] Add determinism test: factor same matrix twice with same `Par::rayon(4)`, assert bitwise-identical `FrontFactors` (pivot sequences, L, D entries), in `src/aptp/numeric.rs` tests (new `#[test] fn test_parallel_factor_determinism`)
- [X] T012 [US1] Add small-matrix overhead test: time factorization of small matrices (n < 1000) with `Par::Seq` vs `Par::rayon(4)`, assert < 5% slowdown, in `src/aptp/solver.rs` tests (new `#[test] fn test_parallel_overhead_small_matrices`)

### Implementation for User Story 1

- [X] T013 [US1] Add `par` parameter to `apply_and_check()` in `src/aptp/factor.rs` and replace hardcoded `Par::Seq` at TRSM call (line ~1278) with `par`
- [X] T014 [US1] Add `par` parameter to `update_trailing()` in `src/aptp/factor.rs` and replace hardcoded `Par::Seq` at GEMM call (line ~1426) with `par`
- [X] T015 [US1] Thread `par` from `AptpOptions` through `factor_inner()` → `apply_and_check()` + `update_trailing()` in `src/aptp/factor.rs`
- [X] T016 [US1] Thread `par` through `two_level_factor()` → `factor_inner()` in `src/aptp/factor.rs`
- [X] T017 [US1] In `AptpNumeric::factor()` in `src/aptp/numeric.rs`, compute `effective_par` per supernode: `Par::Seq` if front dimension < `INTRA_NODE_THRESHOLD` (256), else `options.par`. Create a per-supernode `AptpOptions` copy with the effective par and pass to `aptp_factor_in_place()`
- [X] T018 [US1] In per-supernode solve functions in `src/aptp/solve.rs`, replace the 4 hardcoded `Par::Seq` calls (lines ~147, ~170, ~250, ~261) with the `par` parameter, applying the same front-size threshold (use `Par::Seq` if `ne < INTRA_NODE_THRESHOLD`)
- [X] T019 [US1] Run `cargo test`, all SuiteSparse `--ignored` tests with `--test-threads=1`, and verify all pass. Run `cargo clippy --all-targets --features diagnostic`

**Checkpoint**: Intra-node BLAS parallelism active for large fronts. 63% of SuiteSparse matrices benefit. SC-001, SC-002, SC-004, SC-005, SC-006, SC-007 can be partially validated.

---

## Phase 4: User Story 2 — Tree-Level Parallel Factorization (Priority: P2)

**Goal**: Refactor the sequential postorder factorization loop into a recursive function with `rayon::scope` at tree forks. Independent subtrees process in parallel. Combined with US1's intra-node parallelism for Mixed-class matrices.

**Independent Test**: Factor Mixed-class matrices (e.g., `ncvxqp3`, `cfd2`) and TreeLevel-class matrices (e.g., `bloweybq`) with `Par::rayon(4)` and verify backward error < 5e-11 and measurable speedup on Mixed-class.

**Depends on**: Phase 3 (US1) for Par API surface and intra-node kernel.

### Tests for User Story 2

- [X] T020 [US2] Add regression test: verify sequential recursive traversal produces identical results to original postorder loop, in `src/aptp/numeric.rs` tests (new `#[test] fn test_recursive_factor_matches_sequential`)
- [X] T021 [US2] Add tree-level parallel correctness test: factor all 65 SuiteSparse matrices with `Par::rayon(4)`, assert backward error < 5e-11 for each, in `src/aptp/solver.rs` tests (new `#[ignore] #[test] fn test_parallel_tree_level_all_suitesparse`)
- [X] T022 [US2] Add tree-level determinism test: factor a Mixed-class matrix twice with `Par::rayon(4)`, assert bitwise-identical factors, in `src/aptp/numeric.rs` tests (new `#[test] fn test_parallel_tree_determinism`)

### Implementation for User Story 2

- [X] T023 [US2] Extract the per-supernode factorization loop body in `AptpNumeric::factor()` in `src/aptp/numeric.rs` into a standalone `fn factor_supernode(s, supernodes, children, contributions, matrix, perm_fwd, perm_inv, options, scaling, global_to_local) -> (FrontFactors, Option<ContributionBlock>, PerSupernodeStats)` helper function
- [X] T024 [US2] Replace sequential `front_factors_vec: Vec<FrontFactors>` with pre-allocated `Vec<Option<FrontFactors>>` (indexed by supernode) in `src/aptp/numeric.rs`, and same for `per_sn_stats`
- [X] T025 [US2] Refactor `AptpNumeric::factor()` to use a recursive `fn factor_subtree(root, ...)` function that: (a) recursively processes children, (b) assembles and factors the root supernode, (c) stores results at index `root` in pre-allocated vectors. When `par == Par::Seq`, recurse sequentially. In `src/aptp/numeric.rs`
- [X] T026 [US2] Replace shared `global_to_local: Vec<usize>` with per-task allocation in `factor_supernode()` in `src/aptp/numeric.rs` — each parallel supernode allocates its own `vec![NOT_IN_FRONT; n]`, builds its mapping, and drops it after extraction
- [X] T027 [US2] Add `rayon::scope` dispatch in `factor_subtree()` when `par` is `Par::Rayon`: spawn each child via `s.spawn(move |_| factor_subtree(child, ...))`, implicit barrier ensures all children complete before parent assembly. In `src/aptp/numeric.rs`
- [X] T028 [US2] Compute aggregate `FactorizationStats` from collected `Vec<Option<PerSupernodeStats>>` after the tree traversal completes, replacing the incremental accumulation pattern. In `src/aptp/numeric.rs`
- [X] T029 [US2] Handle `contributions` via return-value pattern: `factor_supernode()` returns `(FrontFactors, Option<ContributionBlock>, PerSupernodeStats)` as owned values. `factor_subtree()` collects children's returned `ContributionBlock`s after the `rayon::scope` barrier (children return via scoped closures), then passes them to the parent supernode's assembly. No shared mutable `Vec<Option<ContributionBlock>>` is needed — all data flows through return values. No `unsafe` code. In `src/aptp/numeric.rs`
- [X] T030 [US2] Run `cargo test`, all SuiteSparse `--ignored` tests with `--test-threads=1`, and `cargo clippy --all-targets --features diagnostic`

**FR-008 coverage**: All parallel code uses safe Rust (no `unsafe`, no `UnsafeCell`). Rust's type system statically prevents data races at compile time. The return-value pattern for contributions avoids shared mutable state entirely. No dedicated ThreadSanitizer task is needed — successful compilation is sufficient proof of FR-008 compliance.

**Checkpoint**: Tree-level parallelism active. All 65 SuiteSparse matrices pass. Mixed/TreeLevel classes show measurable speedup. SC-003, SC-006 can be validated.

---

## Phase 5: User Story 3 — Parallel Triangular Solve (Priority: P3)

**Goal**: Parallelize the three-phase solve (forward, diagonal, backward). Diagonal solve is embarrassingly parallel. Forward/backward use tree-level scheduling with per-task workspace.

**Independent Test**: Factor then solve with `Par::rayon(4)`, verify backward error < 5e-11 and solve time reduction on matrices with wide trees.

**Depends on**: Phase 3 (US1) for Par API in solve functions. Phase 4 (US2) for tree structure (children map) needed by tree-level solve.

### Tests for User Story 3

- [X] T031 [US3] Add parallel solve correctness test: solve CI subset with `Par::rayon(4)`, assert backward error < 5e-11, in `src/aptp/solve.rs` tests (new `#[test] fn test_parallel_solve_correctness`)
- [X] T032 [US3] Add parallel solve determinism test: solve same matrix twice with `Par::rayon(4)`, assert bitwise-identical solution vectors, in `src/aptp/solve.rs` tests (new `#[test] fn test_parallel_solve_determinism`)

### Implementation for User Story 3

- [X] T033 [US3] Implement parallel diagonal solve in `aptp_solve()` in `src/aptp/solve.rs`: when `par` is `Par::Rayon`, use Rayon parallel iteration over supernodes for the diagonal solve phase (each supernode reads/writes only its own `col_indices`, no overlaps). Allocate per-task `work` buffer within each parallel iteration.
- [X] T034 [US3] Expose children map from `AptpNumeric` (or pass `AptpSymbolic` + children to `aptp_solve`) so tree-level solve can use the same tree structure as factorization. In `src/aptp/numeric.rs` and `src/aptp/solve.rs`
- [X] T035 [US3] Implement parallel forward solve in `aptp_solve()` in `src/aptp/solve.rs`: recursive tree traversal with `rayon::scope` — children are processed in parallel (scope barrier), then parent performs its scatter. Per-task `work`/`work2` allocation within scope closures.
- [X] T036 [US3] Implement parallel backward solve in `aptp_solve()` in `src/aptp/solve.rs`: reverse recursive tree traversal with `rayon::scope` — parent processes first, then children in parallel (scope barrier). Per-task workspace allocation.
- [X] T037 [US3] When `par == Par::Seq`, keep existing sequential solve loops (no rayon overhead). In `src/aptp/solve.rs`
- [X] T038 [US3] Run `cargo test`, all SuiteSparse `--ignored` tests with `--test-threads=1`, and `cargo clippy --all-targets --features diagnostic`

**Checkpoint**: All three solve phases parallelized. Solve produces identical backward error. SC-004 (solve determinism), SC-005, SC-006 can be validated.

---

## Phase 6: User Story 4 — Performance Benchmarking & Scaling Report (Priority: P4)

**Goal**: Build a benchmark tool that produces a structured parallel scaling report comparing performance across thread counts and matrix classes.

**Independent Test**: Run the benchmark tool on CI subset and verify the report contains per-matrix speedup, efficiency, and workload classification.

**Depends on**: Phases 3-5 (US1+US2+US3 for parallel factor + solve to benchmark).

### Implementation for User Story 4

- [X] T039 [US4] Create `examples/parallel_scaling.rs` with `required-features = ["diagnostic"]`: accepts `--ci-only` / `--all` flags and `--threads 1,2,4,8` parameter. Loads matrices from registry, runs factorization and solve at each thread count, collects timing.
- [X] T040 [US4] Implement speedup and efficiency computation in `examples/parallel_scaling.rs`: speedup = T_1 / T_n, efficiency = speedup / n. Per-matrix and per-class (IntraNode/Mixed/TreeLevel) aggregation.
- [X] T041 [US4] Implement structured JSON output in `examples/parallel_scaling.rs` to `target/benchmarks/parallel/scaling-<timestamp>.json`: per-matrix records with name, n, nnz, class, and per-thread-count timing (factor_ms, solve_ms, backward_error, speedup, efficiency).
- [X] T042 [US4] Implement human-readable summary table in `examples/parallel_scaling.rs`: per-class average speedup at each thread count, matrices meeting SC-001/SC-002/SC-003 targets, any regressions.
- [X] T043 [US4] Add `[[example]]` entry for `parallel_scaling` with `required-features = ["diagnostic"]` in `Cargo.toml`
- [X] T044 [US4] Run the benchmark tool on CI subset with `--ci-only --threads 1,2,4,8` and verify report is produced with expected fields. Document results in `dev/phase-8.2-report.md`.

**Checkpoint**: Parallel scaling report produced. SC-008 validated. Performance targets (SC-001 through SC-003) evaluated.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, cleanup, and final validation across all user stories.

- [X] T045 Update `dev/ssids-log.md` with Phase 8.2 development log entry
- [X] T046 Update `dev/ssids-plan.md` with Phase 8.2 completion status
- [X] T047 Update `CLAUDE.md` (sparse) with Phase 8.2 status, Rayon dependency, parallelism API notes
- [X] T048 [P] Add rustdoc documentation for all new/modified public APIs: `FactorOptions.par`, `SolverOptions.par`, `SparseLDLT::solve_in_place` `par` parameter, with `# Examples` sections showing parallel usage
- [X] T049 [P] Add academic attribution comments to parallel factorization and solve code: cite Duff 2020 Sections 5-6, Hogg 2016 Section 2.5/3.5, and SPRAL `NumericSubtree.hxx` (BSD-3)
- [X] T050 Run full validation: `cargo test`, `cargo test --features diagnostic`, `cargo test -- --ignored --test-threads=1`, `cargo clippy --all-targets --features diagnostic -- -D warnings`, `cargo fmt --check`, `cargo doc --no-deps`

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately
- **Phase 2 (Foundational)**: Depends on Phase 1 — BLOCKS all user stories
- **Phase 3 (US1 — Intra-Node)**: Depends on Phase 2
- **Phase 4 (US2 — Tree-Level)**: Depends on Phase 3 (needs intra-node Par threading in kernel)
- **Phase 5 (US3 — Parallel Solve)**: Depends on Phase 3 (needs Par in solve API). Can start in parallel with Phase 4 if tree structure is exposed separately.
- **Phase 6 (US4 — Benchmarking)**: Depends on Phases 3+4+5 (needs all parallelism active)
- **Phase 7 (Polish)**: Depends on all user stories being complete

### User Story Dependencies

```
US1 (Intra-Node BLAS) ──→ US2 (Tree-Level Factor) ──→ US4 (Benchmarking)
                     └──→ US3 (Parallel Solve)   ──→ US4 (Benchmarking)
```

- **US1 (P1)**: Can start after Foundational (Phase 2). No story dependencies.
- **US2 (P2)**: Depends on US1 for Par API surface and intra-node kernel.
- **US3 (P3)**: Depends on US1 for Par API in solve. Can partially overlap with US2.
- **US4 (P4)**: Depends on US1+US2+US3 (benchmarks all parallelism strategies).

### Within Each User Story

- Tests written first (TDD per constitution)
- Internal helpers before integration
- Correctness verification before moving to next story

### Parallel Opportunities

- T002 and T003 can run in parallel (different files: solver.rs vs factor.rs)
- T013 and T014 can run in parallel (different functions in factor.rs)
- T031 and T032 (test writing) can run in parallel
- T039-T043 (benchmark tool) tasks are sequential within the example file
- T045-T049 (polish) tasks marked [P] can run in parallel

---

## Parallel Example: User Story 1

```
# After Phase 2 foundational tasks complete:

# These can run in parallel (different functions in factor.rs):
Task T013: "Add par to apply_and_check() in src/aptp/factor.rs"
Task T014: "Add par to update_trailing() in src/aptp/factor.rs"

# Then sequentially:
Task T015: "Thread par through factor_inner()"
Task T016: "Thread par through two_level_factor()"
Task T017: "Compute effective_par per supernode in numeric.rs"
Task T018: "Replace Par::Seq in solve.rs per-supernode functions"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001)
2. Complete Phase 2: Foundational API (T002-T009)
3. Complete Phase 3: User Story 1 — Intra-Node BLAS (T010-T019)
4. **STOP and VALIDATE**: Factor `sparsine` with 4 threads, verify speedup and correctness
5. This alone delivers value for 63% of SuiteSparse matrices

### Incremental Delivery

1. Setup + Foundational → API ready, all tests pass with Par::Seq
2. Add US1 (Intra-Node) → Test → **Immediate speedup for large fronts** (MVP!)
3. Add US2 (Tree-Level) → Test → **Broader coverage for Mixed/TreeLevel matrices**
4. Add US3 (Parallel Solve) → Test → **Solve phase benefits**
5. Add US4 (Benchmarking) → **Quantified scaling report validates all targets**
6. Polish → **Documentation, attribution, final validation**

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Each user story can be validated independently via backward error checks
- Constitution Principle III (TDD): tests written before implementation in each phase
- Commit after each task or logical group
- SuiteSparse `--ignored` tests MUST use `--test-threads=1` to avoid memory pressure
- The `INTRA_NODE_THRESHOLD` (256) is a constant for now; can be promoted to a configurable field later if tuning is needed
