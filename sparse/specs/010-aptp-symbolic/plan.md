# Implementation Plan: APTP Symbolic Analysis

**Branch**: `010-aptp-symbolic` | **Date**: 2026-02-10 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/010-aptp-symbolic/spec.md`

## Summary

Build `AptpSymbolic` as a thin composition over faer's `SymbolicCholesky<usize>`, supplemented with elimination tree parent pointers and APTP-specific delayed-pivot buffer estimates. This is Phase 3 of the SSIDS plan — the "analyze" step of the three-phase analyze→factorize→solve pipeline. The implementation delegates all standard symbolic Cholesky computation to faer and adds only what faer doesn't provide: direct etree access, per-column counts, and heuristic pivot buffer sizing.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (symbolic Cholesky, AMD ordering, MemStack), serde/serde_json (existing, not new)
**Storage**: N/A (in-memory data structures only)
**Testing**: cargo test (unit + integration), criterion 0.5 (benchmarks)
**Target Platform**: Linux (CI), cross-platform (pure Rust)
**Project Type**: Single Rust library crate (existing `rivrs-sparse`)
**Performance Goals**: Symbolic analysis < 5% of numeric factorization time (baseline measured now, validated in Phase 5-6)
**Constraints**: No new dependencies. No heap allocation beyond the result struct itself. Must handle matrices up to ~100K dimension (SuiteSparse test suite).
**Scale/Scope**: ~300-500 lines of new code (one module file + tests)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Validation via full test suite (25 matrices), determinism checks, property tests. No numerical computation in this phase (symbolic only). |
| II. Clean Room | PASS | Implementation uses faer (MIT) public API + academic papers. All references documented in spec. No HSL/GPL code consulted. |
| III. TDD | PASS | Test-first workflow planned: write failing tests → implement → validate. Test categories: unit (edge cases), integration (full suite), property (invariants). |
| IV. Academic Attribution | PASS | Algorithm references in spec: Liu (1990) etree, Gilbert et al. (1992, 1994) symbolic factorization, Hogg et al. (2016) SSIDS analyze phase. |
| V. Numerical Stability | PASS (N/A) | Symbolic phase has no floating-point computation. Buffer estimates are integer heuristics. |
| VI. Structured Development | PASS | Phase 2 complete (exit criteria met). Phase 3 builds on Phase 2 types. |
| VII. Code Quality | PASS | Idiomatic Rust, Result-based errors, rustdoc for all public items, faer types at boundary. |

**Post-Phase 1 re-check**: All gates still PASS. The two-call strategy (prefactorize + factorize_symbolic) is the minimum needed to obtain etree access while using faer's high-level API. No over-engineering detected.

## Project Structure

### Documentation (this feature)

```text
specs/010-aptp-symbolic/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0: faer API research, design decisions
├── data-model.md        # Phase 1: entity definitions, type relationships
├── quickstart.md        # Phase 1: development workflow, usage patterns
├── contracts/
│   └── aptp-symbolic-api.md  # Phase 1: public API contract
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code (repository root)

```text
src/
├── aptp/
│   ├── mod.rs              # MODIFY: add `pub mod symbolic;` + re-exports
│   ├── symbolic.rs         # NEW: AptpSymbolic, SymbolicStatistics
│   ├── pivot.rs            # existing (Phase 2)
│   ├── diagonal.rs         # existing (Phase 2)
│   ├── inertia.rs          # existing (Phase 2)
│   └── perm.rs             # existing (Phase 2)
├── error.rs                # REVIEW: may need context enrichment
├── lib.rs                  # existing (already exports aptp module)
└── ...

tests/
├── symbolic_analysis.rs    # NEW: integration tests for AptpSymbolic
├── hand_constructed.rs     # existing
├── suitesparse_ci.rs       # existing
└── ...

benches/
└── solver_benchmarks.rs    # MODIFY: add symbolic analysis benchmark group
```

**Structure Decision**: Single new file (`src/aptp/symbolic.rs`) added to the existing `aptp` module. AptpSymbolic is part of the APTP algorithm's analyze phase, so it belongs with the other APTP types. Integration tests go in a new `tests/symbolic_analysis.rs` file.

## Key Design Decisions

### D1: Two-Call Strategy for Etree Access

**Problem**: faer's high-level `factorize_symbolic_cholesky` does not expose the elimination tree through its public API. The etree is computed internally but stored in `pub(crate)` fields.

**Solution**: Call two faer functions:
1. `prefactorize_symbolic_cholesky` — gives us `EliminationTreeRef` (parent pointers) and column counts directly
2. `factorize_symbolic_cholesky` — gives us the full `SymbolicCholesky` (permutation, L structure, simplicial/supernodal)

**Cost**: The prefactorize step is O(nnz · α(n)) — almost linear and negligible compared to the full symbolic factorization. The etree is computed on the **permuted** matrix structure (after ordering), so we need to apply the permutation before calling prefactorize — OR call prefactorize on the original matrix and accept that the etree corresponds to the unpermuted structure.

**Important subtlety**: `factorize_symbolic_cholesky` internally computes the etree on the permuted matrix. If we call `prefactorize_symbolic_cholesky` on the original matrix, we get the **unpermuted** etree — which is different from what faer uses internally. For the APTP buffer heuristic and Phase 6 assembly tree, we need the **permuted** etree.

**Resolution**: Apply the permutation to the matrix's symbolic structure first (reindex using the permutation), then call `prefactorize_symbolic_cholesky` on the permuted structure. Alternatively, if faer doesn't provide a convenient way to permute symbolic structure, we can compute it ourselves from the CSC column pointers and the permutation arrays.

### D2: Dual-Variant Handling (Simplicial vs Supernodal)

**Problem**: `SymbolicCholeskyRaw` is an enum — the result is always one variant or the other, decided by faer based on density heuristics.

**Solution**: AptpSymbolic handles both variants transparently:
- `pivot_buffer` length adapts: per-supernode (supernodal) or per-column (simplicial)
- Accessors like `n_supernodes()` return `Option<usize>` (None for simplicial)
- `is_supernodal()` reports which variant was selected
- Pattern match in `analyze()` to compute appropriate buffer estimates

### D3: Pivot Buffer Heuristic

**Formula**: For each unit (supernode or column), estimate:
```
buffer = ceil(fraction * column_count)
```
where `fraction = 0.10` (10%) initially, and `column_count` is the number of nonzeros in that unit's factor columns.

**Tuning**: The fraction is a constant in the implementation, easily adjustable. It will be validated empirically when numeric factorization exists (Phase 5-6).

### D4: Error Propagation

Input validation happens before faer calls:
- Non-square → `SparseError::NotSquare`
- Custom perm wrong dimension → `SparseError::DimensionMismatch`
- faer failure → `SparseError::AnalysisFailure` (wrapping faer's error message)

## Algorithm References

| Phase 3 Concept | Reference | Location |
|------------------|-----------|----------|
| Elimination tree | Liu (1990), "The role of elimination trees in sparse factorization" | foundational theory |
| Elimination tree algorithm | Gilbert et al. (1992), Section 3.3.4 | `references/ssids/gilbert1992.md` |
| Column count prediction | Gilbert, Ng, Peyton (1994) | `references/ssids/gilbert1994.md` |
| Assembly tree / supernodes | Liu (1992), "The Multifrontal Method" | `references/ssids/liu1992.md` |
| SSIDS analyze phase | Hogg et al. (2016), Section 2.2 | `references/ssids/hogg2016.md` |
| Multifrontal assembly tree | Duff, Reid (1984), Section 2 | `references/ssids/duff1984.md` |
| Ordering impact | George (1973) | `references/ssids/george1973.md` |
| faer API | `factorize_symbolic_cholesky`, `prefactorize_symbolic_cholesky` | faer 0.22 source |

## Complexity Tracking

> No constitution violations to justify. All gates pass cleanly.
