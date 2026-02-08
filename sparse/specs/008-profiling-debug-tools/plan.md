# Implementation Plan: Profiling and Debug Tools

**Branch**: `008-profiling-debug-tools` | **Date**: 2026-02-08 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/008-profiling-debug-tools/spec.md`

## Summary

Build profiling, memory tracking, and debug visualization tools to support SSIDS solver development through Phases 2-8. The profiler provides hierarchical timing with Chrome Trace export, the memory tracker provides per-section RSS snapshots, and visualization tools render sparsity patterns and elimination trees as text. All tools are feature-gated behind `test-util`, compile to zero-cost no-ops when disabled, and are thread-safe from the start.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (sparse matrix types), serde/serde_json (Chrome Trace JSON export), std only for timing/threading (no new external deps)
**Storage**: N/A (in-memory only; Chrome Trace exported as JSON file)
**Testing**: cargo test (unit tests with hand-constructed matrices and mock workloads)
**Target Platform**: Linux primary (RSS from /proc/self/status); graceful degradation on other platforms
**Project Type**: Single Rust library crate (existing `rivrs-sparse`)
**Performance Goals**: <1% overhead when profiling active; zero overhead when feature disabled
**Constraints**: No new external dependencies; all gated behind existing `test-util` feature flag
**Scale/Scope**: Development tooling for single-developer workflow; will scale to parallel workloads in Phase 8

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Dev tooling — no solver correctness impact. Tools help verify correctness of future phases. |
| II. Clean Room | PASS | No algorithm implementation from restricted sources. Pure infrastructure code. |
| III. TDD | PASS | All components will have tests written before implementation. Profiler tested with mock workloads, visualizations tested against known matrices. |
| IV. Documentation | PASS | All public APIs will have rustdoc with examples and academic attribution where applicable. |
| V. Numerical Stability | N/A | No numerical computation — measurement and visualization only. |
| VI. Structured Development | PASS | Phase 1.4 of the plan, building on completed Phases 0-1.3. |
| VII. Code Quality | PASS | Follows established patterns (feature gating, module structure, error handling). No new dependencies. |

**Post-Phase 1 re-check**: All gates still pass. Design introduces no new dependencies, no solver code, no restricted-source patterns. Thread-safety design uses only std library primitives.

## Project Structure

### Documentation (this feature)

```text
specs/008-profiling-debug-tools/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0 research findings
├── data-model.md        # Entity definitions
├── quickstart.md        # Usage guide
├── contracts/
│   ├── profiling-api.md # Profiler API contract
│   ├── memory-api.md    # Memory tracker API contract
│   └── debug-viz-api.md # Visualization API contract
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code (repository root)

```text
src/
├── lib.rs                    # Add profiling + debug module declarations
├── profiling/
│   ├── mod.rs                # Module root, re-exports
│   ├── session.rs            # ProfileSession, SectionGuard, FinishedSession
│   ├── section.rs            # ProfileSection, ProfileEvent (data types)
│   ├── report.rs             # Summary report generation + Chrome Trace export
│   └── memory.rs             # MemoryTracker, MemoryReport, MemorySnapshot
├── debug/
│   ├── mod.rs                # Module root, re-exports
│   ├── sparsity.rs           # SparsityDisplay with downsampling
│   └── etree.rs              # ETreeDisplay, EliminationTreeStats
└── benchmarking/
    └── rss.rs                # (existing) RSS reading — shared with profiling/memory.rs
```

**Structure Decision**: Two new peer modules (`profiling/`, `debug/`) alongside existing `benchmarking/` and `testing/`. The existing `benchmarking/rss.rs` RSS reader is shared by importing from `benchmarking` module — no code move needed since both modules are behind `test-util`. A dedicated shared utility module can be factored out later if the dependency direction becomes awkward.

## Complexity Tracking

No constitution violations to justify. All design choices align with existing patterns:
- Feature gating: reuses `test-util` (established in Phase 0.5)
- Module structure: peer modules (same pattern as `benchmarking/`, `testing/`)
- No new external dependencies
- Thread-safety via std primitives only (`thread_local!`, `RefCell`, `Instant`)
