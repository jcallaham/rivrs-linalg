# Implementation Plan: Benchmarking Framework

**Branch**: `006-benchmarking-framework` | **Date**: 2026-02-07 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/006-benchmarking-framework/spec.md`

## Summary

Build a Criterion.rs-based benchmarking framework for the SSIDS sparse solver that measures individual solver phases (analyze, factor, solve) and end-to-end pipelines across the test matrix collection. The framework provides a `Benchmarkable` trait for solver integration, reuses the existing `TestCaseFilter` for matrix selection, records process-wide peak RSS, exports results to CSV/JSON, detects performance regressions against saved baselines, and generates Markdown summary tables. External solver comparison (SPRAL) is deferred.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22, criterion 0.5, serde/serde_json (existing)
**Storage**: JSON files for baselines (`target/benchmarks/baselines/`), CSV for exports
**Testing**: `cargo test` (unit tests for data structures and reporting), `cargo bench` (integration via Criterion)
**Target Platform**: Linux (peak RSS via `/proc/self/status`; benchmark harness is cross-platform)
**Project Type**: Single Rust library with benchmark binaries
**Performance Goals**: CI matrix subset completes in <5 minutes; framework overhead negligible vs solver computation
**Constraints**: Zero overhead on production builds (dev-dependency or feature-gated); no new required dependencies for the library crate
**Scale/Scope**: 15 hand-constructed + 10 CI-subset + 67 full SuiteSparse matrices; 4 benchmark phases; 1 solver (initially)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### Pre-Research Gate

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | **PASS** | Benchmarking is infrastructure — does not affect numerical correctness. Framework measures performance of solver phases without modifying them. |
| II. Clean Room | **PASS** | No algorithm implementation involved. Framework is tooling, not solver code. No restricted sources needed. |
| III. TDD | **PASS** | Framework will be developed test-first: data structures tested via unit tests, Criterion integration verified with mock solver. |
| IV. Documentation | **PASS** | All benchmark modules will have rustdoc with usage examples. Academic attribution N/A (tooling, not algorithms). |
| V. Numerical Stability | **N/A** | Framework does not perform numerical computation. |
| VI. Structured Development | **PASS** | This is Phase 1.2 per `docs/ssids-plan.md`. Phase 0.5 (test infrastructure) is complete. |
| VII. Code Quality | **PASS** | Will follow Rust idioms, use `Result` for fallible ops, feature-gate non-production code. |

### Post-Design Re-Check

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | **PASS** | No change — framework is measurement infrastructure only. |
| III. TDD | **PASS** | Data model types testable via serde round-trip. Regression detection testable with synthetic data. Criterion integration testable with mock solver. |
| VI. Structured Development | **PASS** | Design reuses Phase 0.5 infrastructure (`TestCaseFilter`, `SolverTestCase`). No premature solver code. |
| VII. Code Quality | **PASS** | `Benchmarkable` trait is generic; `Option` returns handle unimplemented phases gracefully. Feature-gated behind `test-util`. |

No violations. Complexity Tracking section not needed.

## Project Structure

### Documentation (this feature)

```text
specs/006-benchmarking-framework/
├── spec.md
├── plan.md              # This file
├── research.md          # Phase 0: Criterion patterns, RSS measurement, trait design
├── data-model.md        # Phase 1: Entity definitions and relationships
├── quickstart.md        # Phase 1: Usage guide
├── contracts/
│   ├── benchmarkable-trait.md   # Benchmarkable trait contract
│   └── benchmark-harness.md     # Harness function contracts
└── checklists/
    └── requirements.md  # Spec quality checklist
```

### Source Code (repository root)

```text
src/
├── lib.rs                          # Add benchmarking module export (behind test-util)
├── benchmarking/
│   ├── mod.rs                      # Module re-exports
│   ├── config.rs                   # BenchmarkConfig, BenchmarkPhase
│   ├── traits.rs                   # Benchmarkable trait definition
│   ├── results.rs                  # BenchmarkResult, BenchmarkSuiteResult, SkippedBenchmark
│   ├── baseline.rs                 # Baseline save/load, regression detection
│   ├── report.rs                   # Markdown table generation, CSV/JSON export
│   └── rss.rs                      # Peak RSS measurement (Linux /proc/self/status)
├── testing/                        # Existing — reused for TestCaseFilter, SolverTestCase
│   ├── mod.rs
│   ├── cases.rs
│   ├── harness.rs
│   ├── validator.rs
│   └── generators.rs

benches/
├── matrix_loading.rs               # Existing — retained as-is or refactored
├── solver_benchmarks.rs            # NEW — Criterion benchmark binary using Benchmarkable trait
```

**Structure Decision**: New `src/benchmarking/` module within the existing library, feature-gated behind `test-util` (same as `testing/`). Benchmark binary `benches/solver_benchmarks.rs` uses both `benchmarking` and `testing` modules. The existing `benches/matrix_loading.rs` is retained for backward compatibility and can be refactored later.
