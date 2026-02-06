# CLAUDE.md - Sparse Linear Algebra Solvers

This file provides guidance to Claude Code when working on sparse linear algebra solvers in this directory.

## Project Overview

This directory contains sparse linear algebra solver implementations for rivrs-linalg, focused on a Sparse Symmetric Indefinite Direct Solver (SSIDS) inspired by SPRAL. The algorithms are implemented using modern Rust with `faer` for high-performance linear algebra operations.

**Parent Project**: rivrs-linalg - Numerical Linear Algebra for Rivrs
**Domain**: Sparse direct solvers (SSIDS, LDL^T factorization, APTP pivoting)
**Current Status**: Phase 0 - Scaffolding and literature review.  Phase 0.1 complete

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

**Planned:**
- 📋 Phase 1: Symbolic analysis (ordering, elimination tree)
- 📋 Phase 2-8: Simplicial APTP solver
- 📋 Phase 9: Supernodal optimization
- 📋 Phase 10-11: Polish and release

## Active Technologies
- Rust 1.87+ with faer (>= 0.22) for linear algebra and sparse infrastructure
- approx for test comparisons (dev dependency)
- criterion for benchmarking (dev dependency)
- Rust 1.87+ (edition 2024) + faer 0.22 (sparse/dense LA), serde + serde_json (JSON parsing) (004-repo-setup)
- Filesystem (test-data/ directory with .mtx and .json files) (004-repo-setup)

## Testing Strategy

Primary correctness validation (established in Phase 0.3 decision):
1. **Reconstruction tests**: `||P^T A P - L D L^T|| / ||A|| < 10^-12` (primary oracle)
2. **Backward error**: `||Ax - b|| / (||A|| ||x|| + ||b||) < 10^-10` (solve pipeline)
3. **Hand-constructed matrices**: 15 matrices with analytically known factorizations
4. **Property-based tests**: Inertia, symmetry preservation, permutation validity
5. **SPRAL comparison**: Deferred to Phases 2-8 (performance benchmarking, large-matrix inertia)

## Recent Changes
- 004-repo-setup: Added Rust 1.87+ (edition 2024) + faer 0.22 (sparse/dense LA), serde + serde_json (JSON parsing)
- 003-spral-golden-results: Phase 0.3 deferred. SPRAL golden results infrastructure not built; reconstruction tests adopted as primary correctness oracle (stronger than cross-solver comparison). Constitution updated to v1.1.0.
- 002-test-matrix-collection: 82 test matrices (15 hand-constructed + 67 SuiteSparse). Three-tier storage: hand-constructed in git, 10-matrix CI subset in `suitesparse-ci/` (plain git), full collection in gitignored `suitesparse/` (extracted from `references/ssids/suitesparse.tar.gz` at container build). No Git LFS.
