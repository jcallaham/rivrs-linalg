# Implementation Plan: Core Test Infrastructure

**Branch**: `005-test-infrastructure` | **Date**: 2026-02-06 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/005-test-infrastructure/spec.md`

## Summary

Build a reusable test harness for validating SSIDS solver components across all development phases (2-11). The harness provides: a `SolverTest` trait with per-phase validation methods, a `NumericalValidator` with configurable tolerances wrapping existing `validate.rs` functions, unified `SolverTestCase` loading from the matrix registry with filtering, and random/structured matrix generators for property-based testing. All test infrastructure is gated behind a `test-util` Cargo feature flag.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (sparse/dense LA), serde + serde_json (JSON parsing), rand + rand_distr (random generation, dev-dependency)
**Storage**: Filesystem (test-data/ directory with .mtx and .json files, metadata.json registry)
**Testing**: cargo test + criterion benchmarks
**Target Platform**: Linux (Docker dev container), cross-platform Rust
**Project Type**: Single Rust library crate with `test-util` feature flag
**Performance Goals**: Matrix loading < 100ms per matrix; random generation < 1s for n ≤ 1000; harness overhead < 2s total
**Constraints**: No new runtime dependencies (test infra is dev-only); constitution accuracy standards (reconstruction < 1e-12, backward error < 1e-10)
**Scale/Scope**: 82 test matrices (15 hand-constructed, 67 SuiteSparse); ~6 new source files; ~1500 LOC estimated

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Harness encodes constitution accuracy standards as defaults (recon < 1e-12, berr < 1e-10). All validation uses mathematically rigorous checks. |
| II. Clean Room | PASS | Test infrastructure has no algorithm implementation. References only permissive sources (faer MIT, SPRAL BSD-3 for test patterns). |
| III. TDD | PASS | Test infrastructure itself will be test-driven: write tests for harness types before implementing them. Harness enables TDD for all future phases. |
| IV. Documentation | PASS | All public types and functions will have rustdoc with examples. Research.md documents design decisions. |
| V. Numerical Stability | PASS | Validator uses existing numerically stable validation functions. Generators use diagonal dominance (proven PD guarantee). |
| VI. Structured Development | PASS | This is Phase 1.1 per ssids-plan.md. Phase 0 is complete. |
| VII. Code Quality | PASS | Uses faer types, Result-based error handling, builder patterns. Feature-gated to keep production API clean. |

**Post-design re-check**: All gates still pass. No constitution violations.

## Project Structure

### Documentation (this feature)

```text
specs/005-test-infrastructure/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0 research decisions
├── data-model.md        # Entity definitions
├── quickstart.md        # Usage examples
├── contracts/
│   └── testing-api.md   # API contract (Rust type signatures)
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # (created by /speckit.tasks)
```

### Source Code (repository root)

```text
src/
├── lib.rs               # Add: pub mod testing (feature-gated)
├── error.rs             # Existing (unchanged)
├── io.rs                # Existing (unchanged)
├── io/
│   ├── mtx.rs           # Existing (unchanged)
│   ├── reference.rs     # Existing (unchanged)
│   └── registry.rs      # Existing (unchanged)
├── validate.rs          # Existing (unchanged, functions remain public)
└── testing/             # NEW: test infrastructure module
    ├── mod.rs           # Module declarations and re-exports
    ├── cases.rs         # SolverTestCase, TestMatrixProperties, TestCaseFilter, load_test_cases
    ├── harness.rs       # SolverTest trait, TestResult, MetricResult, TestKind
    ├── validator.rs     # NumericalValidator (wraps validate.rs)
    └── generators.rs    # Random matrix generators, pattern generators

tests/
├── common/
│   └── mod.rs           # Existing (may be simplified or emptied)
├── hand_constructed.rs  # REFACTORED: use testing::load_test_cases + NumericalValidator
└── suitesparse_ci.rs    # REFACTORED: use testing::load_test_cases + TestCaseFilter

benches/
└── matrix_loading.rs    # Existing (unchanged for now)

Cargo.toml               # MODIFIED: add [features] test-util = []
```

**Structure Decision**: Single Rust library crate. New `testing/` module under `src/` with four submodules. Feature-gated behind `test-util` to keep production builds clean. Integration tests refactored to use new harness.

## Implementation Order

### Layer 1: Foundation Types (no dependencies on each other)
1. `testing/mod.rs` — module declarations, re-exports
2. `testing/harness.rs` — `TestResult`, `MetricResult`, `TestKind` structs; `SolverTest` trait definition
3. `testing/cases.rs` — `SolverTestCase`, `TestMatrixProperties`, `TestCaseFilter` structs

### Layer 2: Core Logic (depends on Layer 1)
4. `testing/validator.rs` — `NumericalValidator` (depends on harness types + existing validate.rs)
5. `testing/cases.rs` — `load_test_cases()` function (depends on SolverTestCase + registry)

### Layer 3: Generators (independent of Layers 1-2 for types, uses SparseColMat)
6. `testing/generators.rs` — random symmetric, arrow, tridiagonal, banded generators

### Layer 4: Integration
7. Cargo.toml — add `test-util` feature, self-dev-dependency
8. `src/lib.rs` — add feature-gated `pub mod testing`
9. Refactor `tests/hand_constructed.rs` to use harness
10. Refactor `tests/suitesparse_ci.rs` to use harness
11. Add unit tests for all testing module components

## Key Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Module visibility | Cargo `test-util` feature flag | Accessible to both unit and integration tests; excluded from production builds |
| validate.rs relationship | Wrap, don't replace | Standalone functions remain public; NumericalValidator adds configurable tolerances on top |
| test_roundtrip semantics | Independent end-to-end | Validates full pipeline as one operation; does not compose individual phase tests |
| PD matrix generation | Diagonal dominance | Guaranteed PD by Gershgorin theorem; no eigenvalue computation needed |
| Indefinite generation | Mixed diagonal signs | Controllable inertia by construction; simpler than random-and-hope |
| Module structure | 4 submodules under testing/ | Clean separation of concerns; each file < 400 LOC |

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| faer API changes in future versions | Low | Medium | Pin faer 0.22 in Cargo.toml; API contract documents current signatures |
| SolverTest trait design doesn't fit Phase 2 solver | Medium | Low | Trait is minimal (4 methods); can be extended with default methods later |
| Random generator produces degenerate matrices | Low | Low | Diagonal dominance guarantees PD; test generators themselves |
| Integration test refactoring breaks existing behavior | Low | Medium | Run existing tests before AND after refactoring; SC-007 requires identical behavior |
