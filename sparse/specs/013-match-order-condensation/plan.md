# Implementation Plan: Match-Order Condensation Pipeline

**Branch**: `013-match-order-condensation` | **Date**: 2026-02-13 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/013-match-order-condensation/spec.md`

## Summary

Implement SPRAL-style condensed matching-ordering pipeline that combines MC64 matching
with METIS ordering via cycle condensation. The pipeline guarantees matched 2-cycle pairs
occupy consecutive positions in the elimination order — critical for APTP's 2x2 pivot
detection. Algorithm follows SPRAL's `match_order.f90:mo_split` (BSD-3): decompose
matching into singletons/2-cycles, build condensed graph (~n/2 dimension), run METIS
on condensed graph, expand ordering back to original indices.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (sparse matrix types, permutations), metis-sys 0.3.x (vendored METIS 5.x)
**Storage**: N/A (in-memory graph algorithms)
**Testing**: cargo test (unit + integration), criterion (benchmarks)
**Target Platform**: Linux (x86_64), CI via GitHub Actions
**Project Type**: Single Rust library crate
**Performance Goals**: Condensation overhead negligible vs METIS cost; total pipeline <= 1.5x (MC64 + METIS independently)
**Constraints**: No heap allocation beyond graph construction; existing MC64/METIS APIs unchanged
**Scale/Scope**: Matrices up to ~500K dimension (SuiteSparse collection range)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Evidence |
|-----------|--------|----------|
| I. Correctness First | PASS | SC-001 (100% pair adjacency), SC-002 (unmatched at end), SC-007 (valid symbolic analysis). Comprehensive validation strategy. |
| II. Clean Room | PASS | Algorithm from SPRAL `match_order.f90` (BSD-3). Academic references: Hogg & Scott (2013) HSL_MC80 origin. No HSL source consulted. |
| III. TDD | PASS | Tests designed before implementation: pair adjacency, cycle decomposition, fill quality, round-trip expansion. Hand-constructed + SuiteSparse. |
| IV. Documentation | PASS | Academic attribution to Duff & Koster (2001), Duff & Pralet (2005), SPRAL `match_order.f90`. Rustdoc planned for all public items. |
| V. Numerical Stability | PASS | Feature is structural (graph/permutation), not numerical. Scaling vector passed through unchanged. No floating-point sensitivity. |
| VI. Structured Development | PASS | Phase 4.3 deliverable, builds on completed Phases 4.1 (METIS) and 4.2 (MC64). Exit criteria defined in spec. |
| VII. Code Quality | PASS | Follows established patterns: `Perm<usize>` composition, `SparseError` returns, `Result` propagation. No new dependencies. |

All gates pass. No violations to justify.

## Project Structure

### Documentation (this feature)

```text
specs/013-match-order-condensation/
├── plan.md              # This file
├── research.md          # Phase 0 output
├── data-model.md        # Phase 1 output
├── quickstart.md        # Phase 1 output
├── contracts/           # Phase 1 output
│   └── api.md           # Rust API contract
└── tasks.md             # Phase 2 output (/speckit.tasks)
```

### Source Code (repository root)

```text
src/
└── aptp/
    ├── mod.rs               # Add re-exports for new public items
    ├── matching.rs          # Existing MC64 (unchanged)
    ├── ordering.rs          # Existing METIS (unchanged) + new match_order functions
    └── ...                  # Existing modules (unchanged)

tests/
├── match_order.rs           # New: integration tests for condensed pipeline
├── mc64_matching.rs         # Existing (unchanged, SC-006)
├── metis_ordering.rs        # Existing (unchanged, SC-006)
└── ...

benches/
└── solver_benchmarks.rs     # Add match_order benchmark group
```

**Structure Decision**: New code goes into existing `src/aptp/ordering.rs` alongside
`metis_ordering()`, since match-order condensation is an ordering pipeline variant.
The cycle decomposition and condensed graph logic are internal (not `pub`) helper
functions within `ordering.rs`. Only the top-level `match_order_metis()` and
`MatchOrderResult` are public. This follows the established pattern of keeping
related algorithms in the same module (matching algorithms in `matching.rs`,
ordering algorithms in `ordering.rs`).
