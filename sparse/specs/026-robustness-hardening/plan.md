# Implementation Plan: Robustness — Testing & Hardening

**Branch**: `026-robustness-hardening` | **Date**: 2026-02-25 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/026-robustness-hardening/spec.md`

## Summary

Add robustness testing to the rivrs-sparse SSIDS solver: a two-directional test audit against SPRAL's kernel test suites plus independent value assessment of existing tests; SPRAL-style torture tests (500+ random instances with adversarial perturbations at the dense kernel level); property-based tests via proptest (structural invariants across generated matrices); and adversarial/edge-case input tests (zero panics on malformed inputs). One new dev-dependency (proptest 1.4), two new source modules (perturbations, strategies), one audit document, and solver hardening fixes as needed.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22, rayon 1.x, proptest 1.4 (new dev-dep)
**Storage**: N/A (in-memory numerical computation)
**Testing**: cargo test, proptest, criterion (benchmarks)
**Target Platform**: Linux (primary), macOS (CI)
**Project Type**: Single Rust library crate
**Performance Goals**: N/A (testing feature — no runtime performance impact)
**Constraints**: Torture tests must complete in < 10 minutes; property tests in < 60 seconds for CI
**Scale/Scope**: ~700 lines new code, ~2 new source modules, 1 audit document, test suite pruning

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | **PASS** | This feature directly strengthens correctness validation |
| II. Clean Room | **PASS** | Consulting SPRAL tests (BSD-3) for audit; perturbation helpers implemented from scratch |
| III. Test-Driven Development | **PASS** | This feature *is* the TDD enhancement phase |
| IV. Algorithm Documentation | **PASS** | Audit document provides traceability; academic references cited in spec |
| V. Numerical Stability | **PASS** | Torture tests specifically exercise stability edge cases |
| VI. Structured Development | **PASS** | Phase 9.2 in the phased plan; Phase 9.1 exit criteria met |
| VII. Code Quality | **PASS** | New code follows existing patterns; proptest is MIT-licensed |

**Post-Phase 1 Re-check**: No violations. All new code is test infrastructure behind `test-util` feature flag. No production code changes except defensive guards for adversarial inputs (FR-007).

## Project Structure

### Documentation (this feature)

```text
specs/026-robustness-hardening/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # SPRAL test inventory, integration decisions
├── data-model.md        # Test entities and relationships
├── quickstart.md        # Implementation quickstart
├── contracts/
│   └── test-apis.md     # Perturbation, strategy, and test contracts
└── tasks.md             # Phase 2 output (/speckit.tasks)
```

### Source Code (new and modified files)

```text
src/testing/
├── mod.rs                  # MODIFIED: add perturbations, strategies modules
├── generators.rs           # EXISTING: random matrix generators (foundation)
├── perturbations.rs        # NEW: cause_delays, make_singular, make_dblk_singular
└── strategies.rs           # NEW: proptest Strategy impls for matrices

src/aptp/
├── factor.rs               # MODIFIED: add torture_tests and property_tests cfg(test) modules
└── solver.rs               # MODIFIED: defensive guards for adversarial inputs (0×0, NaN, etc.)

src/error.rs                # POSSIBLY MODIFIED: new error variants if needed for adversarial inputs

tests/
├── property.rs             # NEW: end-to-end property tests via SparseLDLT
└── adversarial.rs          # NEW: edge case and malformed input tests

docs/
└── spral-test-audit.md     # NEW: two-directional audit document
```

**Structure Decision**: Follows existing project layout. New test modules in `src/testing/` behind `test-util` feature flag. New integration tests in `tests/`. Audit document in `docs/`.

## Complexity Tracking

No constitution violations — this section is intentionally empty.
