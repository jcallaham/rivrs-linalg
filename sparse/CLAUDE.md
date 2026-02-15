# CLAUDE.md - Sparse Linear Algebra Solvers

This file provides guidance to Claude Code when working on sparse linear algebra solvers in this directory.

## Project Overview

This directory contains sparse linear algebra solver implementations for rivrs-linalg, focused on a Sparse Symmetric Indefinite Direct Solver (SSIDS) inspired by SPRAL. The algorithms are implemented using modern Rust with `faer` for high-performance linear algebra operations.

**Parent Project**: rivrs-linalg - Numerical Linear Algebra for Rivrs
**Domain**: Sparse direct solvers (SSIDS, LDL^T factorization, APTP pivoting)
**Current Status**: Phase 5 complete — Dense APTP factorization kernel. Ready for Phase 6 (Multifrontal Numeric Factorization).

### Development docs

The development plan for the SSIDS solver lives in docs/ssids-plan.md and is a **living document**.
It should not be updated with arbitrary progress (will be tracked in specs/), but should remain an updated reflection of the plan as it evolves (for instance, if some decision is made to change the plan, this should be included in the plan doc).

The docs/ folder also includes a ssids-log.md that should be used as a sort of changelog for development updates, including what was built, changed, fixed and why, organized by development phases corresponding to the plan document.

## Licensing Strategy - Clean Room Implementation

**This is a clean room implementation to maintain licensing flexibility.**

- **SPRAL (BSD-3)**: May be freely consulted for implementation patterns
- **LAPACK (BSD-3)**: May be freely consulted for dense kernels
- **HSL documentation**: Consult for algorithm descriptions only
- **Primary sources**: Academic papers, textbooks
- **Always document** academic sources used for each implementation
- This approach allows rivrs-sparse to be released under Apache-2.0

## Reference Materials

Available in the parent `references/` directory:

- **faer-rs/**: Pure Rust linear algebra library - primary dependency for matrix operations and sparse infrastructure
- **lapack/**: Reference LAPACK implementation (BSD-licensed) - consult freely for dense kernels
- **ssids/**: Academic literature converted to Markdown format
- **spral/**: SSIDS development plan and reference materials

Key faer infrastructure to leverage:
- CSC (Compressed Sparse Column) storage format
- Elimination tree computation
- AMD (Approximate Minimum Degree) ordering
- Permutation utilities
- Workspace management patterns

## Technology Stack

**Core Linear Algebra:**
- Use `faer` for dense linear algebra and sparse infrastructure (CSC, elimination trees, ordering)
- Reference LAPACK for dense kernel patterns
- Reference SPRAL for sparse solver patterns

**Future Python Bindings:**
- PyO3 for Rust-Python integration (not yet implemented)
- `ndarray` for zero-copy data sharing with NumPy

## Algorithm Architecture

Following the SSIDS plan phases:

- **Symbolic analysis**: Ordering (AMD/metis), elimination tree, symbolic factorization
- **Numeric factorization**: LDL^T with A Posteriori Threshold Pivoting (APTP)
- **Triangular solve**: Forward/backward substitution using factored form
- **Supernodal optimization** (future): Blocked operations on dense frontal matrices

### Key Algorithm: APTP (A Posteriori Threshold Pivoting)

Unlike traditional threshold pivoting (which decides pivots before elimination), APTP:
1. Performs elimination optimistically (assuming 1x1 pivots)
2. Checks stability after the fact
3. Falls back to 2x2 (Bunch-Kaufman) pivots or delays columns when needed

This allows better data locality and parallelism while maintaining numerical stability.

## faer Integration: Transparent Composition

This library is a **specialized extension of faer**, not a competing implementation.
The guiding principle is **transparent composition**: use the highest-level faer API
that gives us what we need, expose faer types at our boundary, and only add new types
for concepts faer doesn't have.

**Concrete guidelines:**

- **Accept faer types as inputs**: `SparseColMat<usize, f64>`, `PermRef`, `SymmetricOrdering`, `SymbolicSparseColMatRef`
- **Return faer types where they fit**: permutations as `Perm<usize>` (not custom wrappers), sparse matrices as `SparseColMat`
- **Compose, don't wrap**: our structs contain faer results as fields with accessor methods that delegate, rather than copying data into our own arrays (e.g., `AptpSymbolic` stores `SymbolicCholesky<usize>` and delegates `perm()`, `predicted_nnz()`)
- **Only add types for genuinely new concepts**: `PivotType`, `MixedDiagonal`, `AptpSymbolic` — things faer doesn't model
- **Prefer high-level faer APIs over low-level primitives**: use `factorize_symbolic_cholesky` over manual `prefactorize_symbolic_cholesky` + `factorize_simplicial_symbolic_cholesky`, unless we need fine-grained control

**Anti-patterns to avoid:**

- Defining a custom `Permutation` struct when `faer::perm::Perm<usize>` works
- Adding type aliases like `AptpMatrix<T>` for `SparseColMat<usize, f64>` — faer types are already well-named
- Re-implementing CSC storage, elimination tree computation, or AMD ordering
- Hiding faer types behind opaque wrappers that prevent interop with the faer ecosystem

## Code Architecture

**Numerical Considerations:**
- Prioritize numerical stability for indefinite systems
- Use structured decompositions rather than direct inversion
- Implement overflow prevention and condition monitoring
- Support both positive definite (fast path) and indefinite (APTP) factorization

**API Design:**
- Three-phase API: analyze → factorize → solve
- Symbolic analysis result is reusable across multiple factorizations with same sparsity pattern
- Use builder patterns for solver configuration (pivot threshold, ordering strategy, etc.)
- Leverage Rust's type system for safety (e.g., factorized matrix types)

**Error Handling:**
- Use Result types for all fallible operations
- Distinguish structural singularity (analysis) from numerical singularity (factorization)
- Provide informative diagnostics (pivot delays, factorization statistics)

## Development Workflow

When implementing a new component:

1. Review SPRAL source (BSD-3) to understand the algorithm structure
2. Identify relevant academic papers and textbooks
3. Review LAPACK source for any dense kernel patterns needed
4. Design idiomatic Rust API leveraging faer's sparse infrastructure
5. Implement using academic references and permissively-licensed code
6. Document academic references and SPRAL routines consulted
7. Add comprehensive tests (hand-constructed matrices, SuiteSparse collection)
8. Benchmark against SPRAL and other solvers
9. Create Python bindings (future)

## Rust guidelines

- Sort imports by: std, external, workspace, crate, super
- Limit scope to minimum necessary (private over `pub(crate)` over `pub`, etc.)
- No `.unwrap()` in non-test code
- Use `cargo fmt --check`, `cargo clippy`, and `cargo test` throughout development to check code quality and test status
- Tests that run the full SuiteSparse test set use `#[ignore]` and should be run with `cargo test -- --ignored --test-threads=1` before considering an implementation complete
- Pay special attention to Rust edition-specific reserved keywords (e.g., `gen` is reserved in Rust 2024)
- When using faer or other Rust crates, always verify API types by reading the actual source/docs before using them. Do not assume tuple-based constructors — check for dedicated types like `Triplet`.

Preferences (can be violated if needed)
- Immutability over mutation
- Iterators over manual loops
- Enums over dynamic dispatch (e.g. `dyn` trait)
- Prefer `Result` propagation over `.expect()`
- Prefer `thiserror` over custom `Result`/`Error` traits
- Importing with `use` over inline imports (exception: module imports are okay if several functions or types are used)
- Try to avoid allocation in performance-critical code

## Git Commit Practices

- **Commit frequently**: After completing logical units (phase completion, tests passing, new module)
- **Commit after verification**: Always commit when tests pass or milestones are verified
- **Descriptive messages**: Use clear commit messages
- **Never accumulate**: Don't accumulate large amounts of uncommitted work

## Testing Strategy

- Unit tests with hand-constructed matrices (known factorizations)
- Integration tests with SuiteSparse Matrix Collection
- Property-based tests for numerical stability
- Comparison tests against SPRAL results
- Benchmark suite for performance regression testing
- Python integration tests (future)

### Validation strategy

Primary correctness validation (established in Phase 0.3 decision):

1. **Reconstruction tests**: `||P^T A P - L D L^T|| / ||A|| < 10^-12` (primary oracle)
2. **Backward error**: `||Ax - b|| / (||A|| ||x|| + ||b||) < 10^-10` (solve pipeline)
3. **Hand-constructed matrices**: 15 matrices with analytically known factorizations
4. **Property-based tests**: Inertia, symmetry preservation, permutation validity
5. **SPRAL comparison**: Deferred to Phases 2-8 (performance benchmarking, large-matrix inertia)

### Test infrastructure

- `SolverTest` trait with `MockSolver` for pre-solver validation
- `NumericalValidator` with configurable tolerances
- `TestCaseFilter` for composable test case selection
- Random matrix generators (PD and indefinite) behind `test-util` feature

### Testing Discipline for Implementation Phases

Five principles guiding test design for Phases 3+ (implementation-heavy):

1. **Validation proportional to novelty**: Code that delegates to faer gets sanity checks and is included in integration tests. Code written from scratch (permutation remapping, assembly tree derivation, APTP kernel) gets mathematical proof-level tests with exact expected values.

2. **Refactor for testability**: Extract non-trivial logic into standalone functions with minimal inputs. Example: `permute_symbolic_upper_triangle` was extracted from `compute_permuted_etree_and_col_counts` so permutation correctness can be tested independently of faer's symbolic Cholesky.

3. **Regression tests before refactoring**: Before optimizing or restructuring code, encode its current correct behavior with exact expected values. These regression tests catch silent breakage during future changes.

4. **Cross-validation with faer**: When both our code and faer compute overlapping quantities, assert consistency. Example: `col_counts.sum() == predicted_nnz` verifies our permuted structure matches faer's internal computation.

5. **Property-based testing**: For implementation-dependent outputs (supernode boundaries, assembly tree shape), test structural properties rather than exact values. Example: assembly tree parent pointers must satisfy postorder (parent > child), total children + roots = total supernodes.


## Documentation Standards

- Document the mathematical operation performed
- **Cite academic references** used during implementation
- Include references to sparse direct solver literature
- Provide usage examples
- Document numerical stability characteristics and failure modes
- Cross-reference SPRAL routine names for users migrating from Fortran
- Add proper attribution as dictated by open-source community standards

## Current Implementation Status

**Completed:**
- ✅ Phase 0.1: Literature review and reference library
- ✅ Phase 0.2: Test matrix collection (82 matrices)
- ✅ Phase 0.3: Deferred (reconstruction tests adopted as primary oracle)
- ✅ Phase 0.4: Repository setup (IO modules, validation, CI, benchmarks)
- ✅ Phase 1.1: Test infrastructure (harness, validator, generators, test-util feature)

## Active Technologies
- Rust 1.87+ (edition 2024) + faer 0.22 (sparse/dense LA)
- serde + serde_json (JSON parsing for metadata and reference factorizations)
- rand + rand_distr (random matrix generation, optional via `test-util` feature)
- approx (test comparisons, dev dependency)
- criterion (benchmarking, dev dependency)
- Filesystem (test-data/ directory with .mtx and .json files, metadata.json registry)
- Rust 1.87+ (edition 2024) + faer 0.22, criterion 0.5, serde/serde_json (existing) (006-benchmarking-framework)
- JSON files for baselines (`target/benchmarks/baselines/`), CSV for exports (006-benchmarking-framework)
- YAML (GitHub Actions workflow) + Rust 1.87+ (edition 2024) + GitHub Actions (`actions/checkout@v4`, `dtolnay/rust-toolchain@master`, `Swatinem/rust-cache@v2`) (007-ci-setup)
- Rust 1.87+ (edition 2024) + faer 0.22 (sparse matrix types), serde/serde_json (Chrome Trace JSON export), std only for timing/threading (no new external deps) (008-profiling-debug-tools)
- Rust 1.87+ (edition 2024) + faer 0.22 (sparse/dense LA), serde + serde_json (serialization for Inertia) (009-aptp-data-structures)
- N/A (in-memory data structures only) (009-aptp-data-structures)
- Rust 1.87+ (edition 2024) + faer 0.22 (symbolic Cholesky, AMD ordering, MemStack), serde/serde_json (existing, not new) (010-aptp-symbolic)
- Rust 1.87+ (edition 2024) + faer 0.22 (existing), metis-sys 0.3.x (new — vendored METIS 5.x C source) (011-metis-ordering)
- N/A (in-memory graph algorithms) (011-metis-ordering)
- Rust 1.87+ (edition 2024) + faer 0.22 (sparse matrix types, permutations), std::collections::BinaryHeap (Dijkstra priority queue) (012-mc64-matching-scaling)
- Rust 1.87+ (edition 2024) + faer 0.22 (sparse matrix types, permutations), metis-sys 0.3.x (vendored METIS 5.x) (013-match-order-condensation)
- Rust 1.87+ (edition 2024) + faer 0.22 (dense matrix types, matmul, triangular solve), Phase 2 types (MixedDiagonal, PivotType, Block2x2, Inertia) (014-dense-aptp-kernel)
- N/A (in-memory dense matrices only) (014-dense-aptp-kernel)

## Recent Changes
- 010-aptp-symbolic: Added Rust 1.87+ (edition 2024) + faer 0.22 (symbolic Cholesky, AMD ordering, MemStack), serde/serde_json (existing, not new)
- 009-aptp-data-structures: Added Rust 1.87+ (edition 2024) + faer 0.22 (sparse/dense LA), serde + serde_json (serialization for Inertia)
- 009-aptp-data-structures: Added [if applicable, e.g., PostgreSQL, CoreData, files or N/A]
