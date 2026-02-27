# Implementation Plan: APTP Data Structures

**Branch**: `009-aptp-data-structures` | **Date**: 2026-02-08 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/009-aptp-data-structures/spec.md`

## Summary

Implement the data structures unique to indefinite APTP factorization (Phase 2 of ssids-plan.md): `PivotType` enum for pivot classification, `Block2x2` for 2x2 symmetric diagonal blocks, `MixedDiagonal` for the mixed D factor with solve and inertia computation, `perm_from_forward` utility for constructing faer `Perm<usize>` from ordering output, and relocation of `Inertia` to its proper domain module.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (sparse/dense LA), serde + serde_json (serialization for Inertia)
**Storage**: N/A (in-memory data structures only)
**Testing**: cargo test (unit tests in module, integration tests in tests/)
**Target Platform**: Linux (development), cross-platform (library)
**Project Type**: Single Rust library crate
**Performance Goals**: Not performance-critical at this phase — data structures support correctness validation. `solve_in_place` will be called once per solve, not in an inner loop.
**Constraints**: Numerical accuracy: `solve_in_place` relative error < 10^-14. No new external dependencies.
**Scale/Scope**: Dimensions up to ~1.6M (largest SuiteSparse matrix). 5 new source files, ~500-700 lines of production code, ~400-600 lines of tests.

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Evidence |
|-----------|--------|----------|
| I. Correctness First | PASS | Primary validation via reconstruction tests (SC-002 < 10^-14), inertia cross-check against 15 hand-constructed references (SC-003), debug-asserts for invariants |
| II. Clean Room | PASS | Data structures derived from academic literature (Hogg/Duff/Lopez 2020, Bunch-Kaufman). SPRAL (BSD-3) consulted for design patterns. No HSL source consulted. |
| III. TDD | PASS | Tests will be written before implementation per constitution. Unit tests for each type, solve round-trip tests, edge case tests. |
| IV. Documentation | PASS | All public types and functions will have rustdoc with academic references, algorithm descriptions, and examples. |
| V. Numerical Stability | PASS | 2x2 solve uses Cramer's rule (exact for 2x2). Debug-asserts catch singular pivots. Inertia uses trace/determinant (no sqrt needed). |
| VI. Structured Development | PASS | Phase 2 follows completed Phases 0-1. Exit criteria defined in ssids-plan.md. |
| VII. Code Quality | PASS | Follows Rust idioms, transparent composition with faer, `Result` for fallible operations, debug-asserts for invariants. |

**Post-design re-check**: All gates still pass. No new dependencies, no complex patterns, no architecture violations.

## Project Structure

### Documentation (this feature)

```text
specs/009-aptp-data-structures/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0: faer API research, 2x2 eigenvalue math
├── data-model.md        # Phase 1: entity definitions and relationships
├── quickstart.md        # Phase 1: usage examples and build instructions
├── contracts/
│   └── aptp-api.md      # Phase 1: complete public API contract
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code (repository root)

```text
src/
├── lib.rs                  # Modified: add `pub mod aptp;`
├── aptp/
│   ├── mod.rs              # New: module hub, re-exports
│   ├── pivot.rs            # New: PivotType, Block2x2
│   ├── diagonal.rs         # New: MixedDiagonal
│   ├── inertia.rs          # New: Inertia (relocated from io/reference.rs)
│   └── perm.rs             # New: perm_from_forward
├── io/
│   └── reference.rs        # Modified: add `pub use crate::aptp::Inertia;`
└── ... (all other files unchanged)

tests/
└── aptp_data_structures.rs # New: integration tests for APTP types
```

**Structure Decision**: New `aptp/` module directory within existing `src/`, following the same pattern as `io/`, `testing/`, `profiling/`, `debug/`. This module is NOT feature-gated behind `test-util` — it is production code used by the solver.

## Complexity Tracking

No constitution violations. No complexity justifications needed.

## Design Decisions

### D1: Module organization — flat submodules under `aptp/`

Each concern gets its own file (`pivot.rs`, `diagonal.rs`, `inertia.rs`, `perm.rs`) with a hub `mod.rs` that re-exports public items. This matches the existing codebase pattern (`io/mtx.rs`, `io/reference.rs`, etc.).

### D2: MixedDiagonal internal storage — parallel arrays

`pivot_map: Vec<PivotType>` (one per column) + `diag_1x1: Vec<f64>` (length n, unused entries at 2x2/delayed cols) + `blocks_2x2: Vec<Block2x2>` (one per pair). This is cache-friendly for the solve inner loop (iterates pivot_map sequentially, accesses diag_1x1 or blocks_2x2 based on type).

Alternative rejected: single `Vec<DiagonalEntry>` enum — would require matching on every access and has worse data locality due to mixed-size entries.

### D3: Inertia relocation strategy — move + re-export

Move `Inertia` struct definition to `aptp/inertia.rs`. In `io/reference.rs`, replace the struct definition with `pub use crate::aptp::Inertia;`. All existing code (`validate.rs`, tests) that imports `crate::io::reference::Inertia` continues to work without changes. The `Serialize`/`Deserialize` derives stay on Inertia (needed for JSON reference data).

### D4: perm_from_forward signature — `Vec<usize>` input, `Result` return

Takes ownership of the forward array (avoids a clone when converting to `Box<[usize]>`). Returns `Result<Perm<usize>, SparseError>` to handle invalid input. Validates using the existing `validate_permutation` function before computing the inverse.

### D5: 2x2 inertia classification — trace/determinant

Uses `det = a*c - b*b` and `trace = a + c` to classify eigenvalue signs without computing actual eigenvalues. This is exact in IEEE 754 arithmetic for well-conditioned blocks (no sqrt, no iteration). See research.md R2 for the complete decision table.

## References

### Traceability Requirement

All algorithm implementations in this feature MUST include rustdoc citations tracing to specific academic references or permissively-licensed source code. This is required by Constitution principles II (Clean Room) and IV (Documentation & Attribution). Each public function's doc comment should cite the relevant paper(s) from the list below by author/year shorthand (e.g., "Hogg, Duff & Lopez 2020, Section 3.2").

### Academic Sources

| Citation | Paper | Local Reference | Relevance to This Feature |
|----------|-------|-----------------|--------------------------|
| Hogg, Duff & Lopez (2020) | "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting", SIAM J. Sci. Comput. 42(4) | `references/ssids/duff2020.md` | **Primary**: APTP algorithm definition, pivot classification (1x1/2x2/delayed), mixed diagonal D storage, delayed-column bookkeeping |
| Hogg, Ovtchinnikov & Scott (2016) | "A Sparse Symmetric Indefinite Direct Solver for GPU Architectures", ACM TOMS | `references/ssids/hogg2016.md` | SSIDS architecture, pivot handling in supernodal context |
| Duff & Pralet (2005) | "Strategies for scaling and pivoting for sparse symmetric indefinite problems" | `references/ssids/duff2005.md` | Pivoting strategies, scaling + pivoting interaction |
| Schenk & Gärtner (2006) | "On fast factorization pivoting methods for sparse symmetric indefinite systems" | `references/ssids/schenk2006.md` | Alternative pivoting approaches (context for APTP design choices) |
| Bunch & Kaufman (1977) | "Some Stable Methods for Calculating Inertia and Solving Symmetric Linear Systems", Math. Comp. | *(not in local refs — standard textbook reference)* | 2x2 pivot strategy, Bunch-Kaufman pivot selection |

### Permissive Reference Code

- **SPRAL (BSD-3)** — `/opt/references/spral/`, pivot handling patterns in `src/ssids/cpu/kernels/ldlt_app.hxx`
- **faer 0.22 (MIT)** — `/opt/references/faer-rs/faer/src/perm/permown.rs` for `Perm::new_checked` API; `linalg/cholesky/bunch_kaufman/` for 2x2 block patterns
