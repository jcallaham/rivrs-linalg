# Feature Specification: Module Restructure

**Feature Branch**: `027-module-restructure`
**Created**: 2026-02-26
**Status**: Draft
**Input**: User description: "Restructure module hierarchy per roadmap items 1-3: extract ordering/ from aptp/, rename aptp/ to symmetric/, update all imports, tests, examples, and documentation. Optimize for future reusability, modularity, and composability."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Library Consumer Uses Ordering Independently (Priority: P1)

A developer building a custom sparse solver (e.g., iterative solver with ILU preconditioner) needs METIS ordering and MC64 matching/scaling without pulling in any symmetric-indefinite solver internals. They import ordering utilities from a top-level `ordering` module.

**Why this priority**: Ordering and matching are the highest-value reusable assets in the crate. Making them independently accessible is the single biggest modularity win and the primary motivation for the restructure.

**Independent Test**: Can be tested by importing `rivrs_sparse::ordering::{metis_ordering, mc64_matching, match_order_metis}` and verifying they work without any `symmetric` module involvement.

**Acceptance Scenarios**:

1. **Given** a sparse matrix in CSC format, **When** a developer calls `rivrs_sparse::ordering::metis_ordering()`, **Then** they receive a valid fill-reducing permutation without needing any APTP or symmetric solver types.
2. **Given** a sparse matrix, **When** a developer calls `rivrs_sparse::ordering::mc64_matching()`, **Then** they receive matching and scaling results usable with any solver.
3. **Given** the previous `rivrs_sparse::aptp::metis_ordering` import path, **When** a developer updates to the new version, **Then** the function is available at `rivrs_sparse::ordering::metis_ordering` with the same signature and behavior.

---

### User Story 2 - Library Consumer Uses Symmetric Solver (Priority: P1)

A developer using the symmetric indefinite solver (`SparseLDLT`) continues to access all solver functionality through a clearly named module path (`symmetric`) that describes the problem domain rather than the implementation detail (`aptp`).

**Why this priority**: The solver is the primary product of this crate. Renaming must not break the API or degrade discoverability.

**Independent Test**: Can be tested by importing `rivrs_sparse::symmetric::{SparseLDLT, AnalyzeOptions, FactorOptions, OrderingStrategy}` and running the existing end-to-end solve tests.

**Acceptance Scenarios**:

1. **Given** the existing `SparseLDLT` API, **When** the module is renamed from `aptp` to `symmetric`, **Then** all public types and functions remain accessible via `rivrs_sparse::symmetric::*` with identical signatures and behavior.
2. **Given** APTP-specific internals (PivotType, MixedDiagonal, dense kernel), **When** a developer explores the `symmetric` module, **Then** these are available under `rivrs_sparse::symmetric::*` (flat re-export, not requiring knowledge of internal submodule structure).
3. **Given** the existing test suite (546+ tests), **When** all import paths are updated, **Then** every test continues to pass with identical results.

---

### User Story 3 - Crate Maintainer Adds a Second Solver (Priority: P2)

A future developer adding a new solver (e.g., unsymmetric LU, iterative CG) can reuse the `ordering` module directly and can see clear boundaries between general infrastructure (`ordering/`, `io/`, `validate`, `profiling/`) and symmetric-specific code (`symmetric/`).

**Why this priority**: This is the forward-looking architectural goal — the restructure should make the codebase legible and composable for future extensions. It's P2 because it's validated by code organization, not by a runtime test.

**Independent Test**: Can be validated by inspecting the module hierarchy and confirming that no symmetric-specific types leak into general modules.

**Acceptance Scenarios**:

1. **Given** the restructured crate, **When** a developer reads `src/lib.rs`, **Then** the module hierarchy clearly separates general infrastructure from solver-specific code.
2. **Given** the `ordering` module, **When** inspecting its public API, **Then** it has no dependencies on `symmetric` (no `PivotType`, `MixedDiagonal`, `AptpSymbolic`, etc.).
3. **Given** the `symmetric` module, **When** it needs ordering functionality, **Then** it imports from `crate::ordering`, not from internal submodules.

---

### Edge Cases

- What happens to `pub use` re-exports that currently exist at `aptp::*`? They move to `symmetric::*` with identical names.
- What happens to `io/reference.rs` which re-exports `Inertia` from `aptp`? It updates to import from `symmetric`.
- What happens to `testing/mc64_validation.rs` which imports `aptp::matching::count_cycles`? It imports from `ordering::matching::count_cycles`.
- What happens to `amalgamation.rs` (internal to numeric.rs)? It stays inside `symmetric/` as an internal helper — it's inherently tied to the symmetric multifrontal solver.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: Crate MUST extract `ordering.rs`, `matching.rs`, and `perm.rs` from `aptp/` into a new top-level `ordering/` module.
- **FR-002**: Crate MUST rename `aptp/` directory to `symmetric/` and update the module declaration in `lib.rs`.
- **FR-003**: Crate MUST update all `use crate::aptp::*` imports throughout the crate to reference either `crate::ordering::*` or `crate::symmetric::*` as appropriate.
- **FR-004**: Crate MUST update all integration tests (`tests/*.rs`) to use new import paths (`rivrs_sparse::ordering::*` or `rivrs_sparse::symmetric::*`).
- **FR-005**: Crate MUST update all examples (`examples/*.rs`) to use new import paths.
- **FR-006**: Crate MUST update all benchmarks (`benches/*.rs`) to use new import paths.
- **FR-007**: Crate MUST update the comparison driver (`comparisons/src/spral_benchmark.rs`) to use new import paths.
- **FR-008**: Crate MUST preserve the public re-export surface: every type/function previously available at `rivrs_sparse::aptp::X` MUST be available at either `rivrs_sparse::ordering::X` or `rivrs_sparse::symmetric::X`.
- **FR-009**: Crate MUST NOT extract `multifrontal/` as a separate module — there is no second consumer, and premature extraction risks designing wrong abstractions.
- **FR-010**: The `ordering` module MUST have zero dependencies on the `symmetric` module (clean dependency direction: symmetric depends on ordering, not vice versa).
- **FR-011**: Crate MUST update `CLAUDE.md` source code layout documentation to reflect the new structure.
- **FR-012**: Crate MUST update `docs/ssids-log.md` to document the restructure.
- **FR-013**: All existing tests MUST pass after the restructure with bit-identical results.
- **FR-014**: `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` MUST pass with no warnings.
- **FR-015**: `cargo fmt --check` MUST pass.
- **FR-016**: `cargo doc --no-deps` MUST succeed.

### Proposed New Module Structure

```
sparse/src/
├── lib.rs                          # Updated module declarations
├── error.rs                        # General (unchanged)
├── validate.rs                     # General (unchanged)
│
├── io/                             # General (unchanged internally)
│   ├── mod.rs
│   ├── mtx.rs
│   ├── registry.rs
│   └── reference.rs                # Update: import Inertia from symmetric
│
├── ordering/                       # NEW: extracted from aptp/
│   ├── mod.rs                      # Re-exports: metis_ordering, mc64_matching, etc.
│   ├── metis.rs                    # Was aptp/ordering.rs
│   ├── matching.rs                 # Was aptp/matching.rs (moved as-is)
│   └── perm.rs                     # Was aptp/perm.rs (moved as-is)
│
├── symmetric/                      # RENAMED from aptp/
│   ├── mod.rs                      # Re-exports (same items as old aptp/mod.rs minus ordering)
│   ├── pivot.rs                    # Unchanged
│   ├── diagonal.rs                 # Unchanged
│   ├── inertia.rs                  # Unchanged
│   ├── symbolic.rs                 # Unchanged
│   ├── factor.rs                   # Update: use crate::ordering::perm
│   ├── numeric.rs                  # Unchanged (no ordering deps)
│   ├── solve.rs                    # Unchanged
│   ├── solver.rs                   # Update: use crate::ordering::{match_order_metis, metis_ordering}
│   └── amalgamation.rs             # Unchanged (internal)
│
├── profiling/                      # General (unchanged)
├── debug/                          # General (unchanged)
├── testing/                        # General (update mc64_validation imports)
└── benchmarking/                   # General (unchanged)
```

### Key Design Decisions

1. **Rename `aptp/ordering.rs` to `ordering/metis.rs`**: The file contains `metis_ordering()` and `match_order_metis()` — its content is METIS-specific, not generic ordering. The filename `metis.rs` is more descriptive than `ordering.rs` inside an `ordering/` module.

2. **Do NOT create `multifrontal/`**: Per roadmap Section 4's "extract on second consumer" principle. The workspace, assembly maps, tree traversal, and amalgamation code are all tightly coupled to the symmetric factorization. Extracting them now would require designing generic interfaces without a concrete second use case.

3. **Flat re-exports in `symmetric/mod.rs`**: Users should not need to know about the internal `pivot`, `diagonal`, `factor` submodule structure. All public types are re-exported at the `symmetric` level, matching the current `aptp/mod.rs` pattern.

4. **`ordering` depends on nothing solver-specific**: The `ordering` module only depends on `crate::error` (for `SolverError`). The `symmetric` module depends on `crate::ordering` for `match_order_metis`, `metis_ordering`, and `perm_from_forward`.

### Key Entities

- **`ordering` module**: Top-level module providing fill-reducing orderings (METIS), weighted bipartite matching (MC64), and permutation utilities. Independent of any specific solver.
- **`symmetric` module**: Top-level module providing the symmetric indefinite direct solver (SparseLDLT) with APTP pivoting. Depends on `ordering` for preprocessing.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: All existing tests pass after restructure (`cargo test` and `cargo test --features diagnostic`).
- **SC-002**: All 65 SuiteSparse integration tests pass (`cargo test -- --ignored --test-threads=1`).
- **SC-003**: Zero `use crate::aptp::` or `rivrs_sparse::aptp::` references remain in the codebase after restructure.
- **SC-004**: The `ordering` module has zero import dependencies on the `symmetric` module (verified by grep).
- **SC-005**: `cargo clippy --all-targets --features diagnostic` produces zero warnings.
- **SC-006**: All examples compile and run correctly (`cargo build --examples`).
- **SC-007**: `cargo doc --no-deps` succeeds with zero warnings.
- **SC-008**: `CLAUDE.md` source code layout section accurately reflects the new module structure.
