# Tasks: Continuous Integration Setup

**Input**: Design documents from `specs/007-ci-setup/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, quickstart.md

**Tests**: Not explicitly requested in the feature specification. Test tasks are omitted.

**Organization**: Tasks are grouped by user story. Due to the gap analysis finding that 8 of 10 requirements are already satisfied, most user stories require only verification, not implementation. The single implementation task (bench-sparse job) maps to User Story 5.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

## Phase 1: Setup

**Purpose**: No setup needed — the CI workflow file and all GitHub Actions dependencies already exist.

*(No tasks — `.github/workflows/ci.yml` and all referenced actions are already in place)*

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: No foundational work needed — the existing CI configuration provides the foundation (triggers, environment variables, caching, independent domain jobs) on which the new job builds.

*(No tasks — existing workflow structure is the foundation)*

---

## Phase 3: User Story 5 - Performance Regression Awareness (Priority: P3) — The Only Implementation Work

**Goal**: Add a `bench-sparse` job to `.github/workflows/ci.yml` that compiles the benchmark binary without executing it, preventing the benchmarking infrastructure from bit-rotting.

**Independent Test**: After merging, open a PR that modifies benchmark code and confirm the bench-sparse job appears in the CI checks and passes. Alternatively, verify locally with `cd sparse && cargo bench --no-run`.

**Why this is the only implementation phase**: The gap analysis (plan.md, research.md) determined that User Stories 1-4 (automated testing, code quality, documentation, feature-gated coverage) are already fully satisfied by the existing CI configuration:

| User Story | Spec Requirement | Existing CI Job | Status |
|------------|-----------------|-----------------|--------|
| US1 (P1) — Test validation | FR-001, FR-002 | test-sparse (stable + 1.87 matrix) | Already satisfied |
| US2 (P1) — Code quality | FR-003, FR-004 | lint-sparse (fmt + clippy) | Already satisfied |
| US3 (P2) — Documentation | FR-005 | doc-sparse (rustdoc -D warnings) | Already satisfied |
| US4 (P2) — Feature-gated tests | FR-006 | test-sparse (self-ref dev-dep activates test-util) | Already satisfied |
| US5 (P3) — Benchmarks | FR-007 | **Missing** | **Needs implementation** |

Cross-cutting requirements FR-008 (caching), FR-009 (clear output), FR-010 (independent jobs) are also already satisfied.

### Implementation for User Story 5

- [x] T001 [US5] Add `bench-sparse` job to `.github/workflows/ci.yml` that runs `cargo bench --no-run` in the sparse directory, using stable toolchain and Swatinem/rust-cache@v2, following the existing job naming and structure conventions (see doc-sparse job as template)

---

## Phase 4: Polish & Cross-Cutting Concerns

**Purpose**: Verification and documentation updates

- [x] T002 Verify all existing CI jobs still pass by running local equivalents: `cd sparse && cargo test --all-targets && cargo fmt --check && cargo clippy --all-targets -- -D warnings && RUSTDOCFLAGS="-D warnings" cargo doc --no-deps && cargo bench --no-run`
- [x] T003 Update `dev/ssids-plan.md` to mark Phase 1.3 success criteria as checked (matching the pattern used for Phases 0.1-0.4, 1.1-1.2)
- [x] T004 Add Phase 1.3 entry to `dev/ssids-log.md` documenting what was built, key decisions, and issues encountered

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: N/A — no tasks
- **Phase 2 (Foundational)**: N/A — no tasks
- **Phase 3 (US5 Implementation)**: No dependencies — can start immediately
- **Phase 4 (Polish)**: T002 depends on T001 completion. T003 and T004 can run in parallel with T002.

### User Story Dependencies

- **US1-US4**: Already satisfied — no implementation work required
- **US5**: Independent, no dependencies on other stories

### Within Phase 3

Single task (T001) — no internal parallelism needed.

### Parallel Opportunities

- T003 and T004 can run in parallel with each other (different files)
- T002 must run after T001 (verification depends on implementation)

---

## Implementation Strategy

### MVP First

1. Complete T001 — add bench-sparse job (the only code change)
2. Complete T002 — verify nothing broke
3. **STOP and VALIDATE**: Push branch, open PR, confirm all 5 sparse CI jobs pass
4. Complete T003-T004 — documentation updates

### Scope Note

This is a minimal-scope feature. The gap analysis revealed that most of the Phase 1.3 plan requirements were already implemented during Phases 0.4 and 1.1-1.2. The remaining work is a single additive YAML change (~15 lines) plus documentation updates.

---

## Notes

- The bench-sparse job should follow the exact structure of the existing doc-sparse job (checkout → toolchain → cache → run) for consistency
- Use `stable` toolchain only for bench-sparse (not MSRV matrix) — benchmarks are a development tool, MSRV compilation adds no value
- The `cargo bench --no-run` command compiles all bench targets including criterion harness without executing
- Do not add path filtering, retry logic, or multi-OS support — these are explicitly deferred per research.md decisions R4
