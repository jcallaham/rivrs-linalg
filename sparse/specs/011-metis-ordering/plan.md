# Implementation Plan: METIS Nested Dissection Ordering

**Branch**: `011-metis-ordering` | **Date**: 2026-02-12 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/011-metis-ordering/spec.md`

## Summary

Integrate METIS nested dissection ordering into the sparse solver via `metis-sys` FFI bindings. The feature adds a single public function (`metis_ordering`) that extracts the graph adjacency structure from a symmetric sparse matrix, calls `METIS_NodeND`, and returns a `Perm<usize>` suitable for `AptpSymbolic::analyze()` via `SymmetricOrdering::Custom`. No changes to the existing symbolic analysis API. METIS is a required dependency with vendored C source (no system library needed).

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (existing), metis-sys 0.3.x (new — vendored METIS 5.x C source)
**Storage**: N/A (in-memory graph algorithms)
**Testing**: `cargo test` (unit + integration), `cargo test -- --ignored` (full SuiteSparse suite)
**Target Platform**: Linux (Docker development environment), general x86_64
**Project Type**: Single Rust library crate
**Performance Goals**: Full SuiteSparse symbolic analysis (67 matrices) with METIS ordering < 2 minutes total
**Constraints**: `metis_sys::idx_t` = `i32` — matrix dimension must fit in `i32` (not a practical limitation for our test suite, max ~150K)
**Scale/Scope**: 1 new source file (~200-300 lines), 1 new test file, modifications to 3 existing files (Cargo.toml, aptp/mod.rs, tests/symbolic_analysis_full.rs)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Validated via permutation validity checks, fill comparison against Hogg et al. (2016) Table III values, and integration with existing reconstruction test infrastructure |
| II. Clean Room | PASS | METIS source is Apache-2.0. `metis-sys` is dual MIT/Apache-2.0. Algorithm reference: Karypis & Kumar (1998). No restricted source code involved |
| III. TDD | PASS | Tests written before implementation: permutation validity, fill reduction vs AMD, paper-reported nnz(L) comparison |
| IV. Documentation | PASS | Algorithm references cited in spec. Rustdoc with academic attribution planned for all public functions |
| V. Numerical Stability | N/A | Ordering is a graph algorithm (combinatorial), not a numerical computation. No floating-point arithmetic involved |
| VI. Structured Development | PASS | Phase 4.1 follows completed Phase 3. Exit criteria defined (SC-001 through SC-006) |
| VII. Code Quality | PASS | Safe Rust wrapper around unsafe FFI. `Result<T, SparseError>` for all errors. Follows existing `aptp` module patterns |

**Post-design re-check**: No changes — design adds one module following existing patterns. No new architectural complexity.

## Project Structure

### Documentation (this feature)

```text
specs/011-metis-ordering/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0: crate selection, FFI details, permutation semantics
├── data-model.md        # Phase 1: adjacency structure, type conversions
├── quickstart.md        # Phase 1: usage examples
├── contracts/
│   └── api.md           # Phase 1: public API contract
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code (repository root)

```text
src/
├── aptp/
│   ├── mod.rs              # MODIFIED: add `pub mod ordering;` and re-export
│   ├── ordering.rs         # NEW: metis_ordering(), extract_adjacency()
│   ├── symbolic.rs         # UNCHANGED
│   ├── diagonal.rs         # UNCHANGED
│   ├── inertia.rs          # UNCHANGED
│   ├── perm.rs             # UNCHANGED
│   └── pivot.rs            # UNCHANGED
├── error.rs                # POSSIBLY MODIFIED: add METIS error variant if needed
├── lib.rs                  # UNCHANGED
└── ...

tests/
├── metis_ordering.rs           # NEW: METIS ordering unit + integration tests
├── symbolic_analysis_full.rs   # MODIFIED: add METIS ordering tests, update MAX_DIM guard
└── ...

Cargo.toml                     # MODIFIED: add metis-sys dependency
```

**Structure Decision**: Single new file (`src/aptp/ordering.rs`) following the existing `aptp` module pattern. The ordering function is a natural addition to the APTP pipeline — it produces the input that `AptpSymbolic::analyze()` consumes. One new test file plus modifications to the existing full-suite tests.

## Key Design Decisions

### 1. Module placement: `aptp::ordering` (not a top-level module)

The ordering function produces input for `AptpSymbolic::analyze()`. Placing it in the `aptp` module keeps related functionality together and follows the existing pattern where all APTP pipeline components live under `src/aptp/`.

### 2. Input type: `SymbolicSparseColMatRef` (not `SparseColMat`)

The ordering only needs the sparsity pattern, not numerical values. Accepting `SymbolicSparseColMatRef` makes this explicit and avoids requiring the caller to have a fully numeric matrix. This matches how `AptpSymbolic::analyze()` accepts its matrix argument.

### 3. Direct `Perm::new_checked` construction (not `perm_from_forward`)

METIS provides both forward (`iperm`) and inverse (`perm`) arrays. Using `Perm::new_checked` directly avoids redundant inverse computation that `perm_from_forward()` would perform.

### 4. Adjacency extraction handles any triangle storage

The `extract_adjacency` helper must handle matrices stored as:
- Upper triangle only
- Lower triangle only
- Full symmetric structure

It symmetrizes the structural pattern and excludes diagonal entries. This ensures robustness regardless of how the input matrix was constructed or loaded.

### 5. Trivial cases handled before FFI call

Matrices with dimension 0, 1, or no off-diagonal entries are handled directly in Rust (return identity permutation) without calling METIS. This avoids edge case behavior in the C library.

## Complexity Tracking

No constitution violations to justify. The design adds minimal complexity: one function, one helper, one dependency.
