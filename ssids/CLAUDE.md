# CLAUDE.md - Sparse Symmetric Indefinite Direct Solvers

This file provides guidance to Claude Code when working on sparse solver algorithms in this directory.

## Project Overview

This directory contains sparse symmetric indefinite direct solver implementations for rivrs-linalg, based on algorithms from SPRAL (Sparse Parallel Robust Algorithms Library). The focus is on multifrontal LDLT factorization for solving symmetric indefinite linear systems Ax = b.

**Parent Project**: rivrs-linalg - Numerical Linear Algebra for Rivrs
**Domain**: Sparse symmetric indefinite direct solvers (similar to HSL MA27/MA57/MA97)
**Current Status**: Phase 1 - Early scaffolding, no algorithms implemented yet

## Licensing Strategy - Clean Room Implementation

**This is a clean room implementation to maintain licensing flexibility.**

HSL library source code is proprietary/restrictive. To avoid licensing issues:

- **SPRAL (BSD-3-Clause)**: May be freely consulted - this is our primary reference
- **HSL documentation**: May be consulted for algorithm descriptions and purpose
- **HSL source code**: NEVER read or copy HSL library source code
- **Primary sources**: Academic papers (Duff & Reid, Hogg et al.), SPRAL source, textbooks
- **Always document** academic sources used for each implementation
- This approach allows rivrs-linalg to be released under Apache-2.0

## Reference Materials

To be added to parent `references/` directory:

- **SPRAL**: Clone from https://github.com/ralgr/spral (BSD-3-Clause licensed)
  - Consult freely for implementation patterns
  - Primary reference for multifrontal LDLT

Available in parent `references/`:
- **faer-rs/**: Pure Rust linear algebra for dense operations within frontal matrices
- **lapack/**: Reference LAPACK for dense BLAS operations

## Technology Stack

**Sparse Matrix Representation:**
- Use `sprs` for sparse matrix data structures (CSC/CSR format)
- Consider custom sparse formats optimized for multifrontal methods

**Dense Linear Algebra:**
- Use `faer` for dense operations within frontal matrices
- Reference LAPACK for dense factorizations

**Graph Algorithms:**
- Implement or use existing crates for:
  - Minimum degree ordering (AMD, approximate minimum degree)
  - Nested dissection
  - Elimination tree construction
  - Assembly tree (supernodal analysis)

**Future Parallelism:**
- Consider `rayon` for task parallelism in multifrontal DAG
- Thread-safe data structures for concurrent factorization

## Algorithm Roadmap

### Phase 1: Basic Infrastructure (Current)
- Sparse matrix representation and validation
- Symbolic analysis framework
- Elimination tree construction

### Phase 2: Core Algorithm
- AMD ordering implementation
- Numeric LDLT factorization (single-threaded)
- Forward/backward substitution
- Basic pivoting strategies (threshold partial pivoting)

### Phase 3: Advanced Features
- Supernodal techniques
- Block operations for cache efficiency
- Memory-efficient out-of-core solvers
- Iterative refinement

### Phase 4: Parallelism
- Task-parallel multifrontal factorization
- DAG scheduling
- NUMA-aware memory management

## Code Architecture

**Key Components:**

1. **Symbolic Analysis**
   - Matrix ordering (AMD, METIS, nested dissection)
   - Elimination tree construction
   - Symbolic factorization (predict fill-in)
   - Supernodal analysis (future)

2. **Numeric Factorization**
   - Multifrontal LDLT with pivoting
   - Assembly of frontal matrices
   - Partial factorization and update matrices
   - Extend-add operations

3. **Solve Phase**
   - Forward substitution (L y = P b)
   - Diagonal solve (D z = y)
   - Backward substitution (L^T P^T x = z)

4. **Utilities**
   - Matrix validation (symmetry check)
   - Condition estimation
   - Residual computation
   - Iterative refinement (future)

**Numerical Considerations:**
- Threshold partial pivoting for numerical stability (Duff & Reid)
- Block pivoting for 2×2 pivots (indefinite systems)
- Delayed pivoting to improve performance
- Static pivoting strategies (future)

**API Design:**
- Builder pattern for solver configuration
- Separate symbolic and numeric phases (analyze-once, factorize-many)
- In-place vs. out-of-place solve options
- Iterative refinement as optional post-processing

## Development Workflow

When implementing sparse solver algorithms:

1. Review SPRAL source code freely (BSD-3-Clause)
2. Review HSL documentation for algorithm descriptions (NOT source code)
3. Identify relevant academic papers:
   - Duff & Reid (1983) - multifrontal method
   - Hogg, Reid & Scott (2010) - DAG-based parallelism
   - Amestoy, Davis & Duff (2004) - approximate minimum degree
4. Design Rust API that is idiomatic and leverages modern features
5. Implement using SPRAL and academic references
6. Document which references informed the implementation
7. Add comprehensive tests with known solutions
8. Benchmark against SPRAL and other libraries (SuiteSparse, Pardiso)
9. Profile and optimize critical paths

## Git Commit Practices

- **Commit frequently**: After completing logical units
- **Commit after verification**: Always commit when tests pass
- **Descriptive messages**: Clear commit messages describing what was accomplished
- **Never accumulate**: Don't accumulate large amounts of uncommitted work

## Testing Strategy

- Unit tests for individual components (ordering, tree construction, etc.)
- Integration tests with known sparse systems from:
  - Matrix Market collection
  - SuiteSparse Matrix Collection
  - SPRAL test suite
- Property-based tests (e.g., A x = b satisfaction)
- Benchmark suite comparing to SPRAL and SuiteSparse
- Numerical accuracy validation (forward error, backward error)

## Documentation Standards

- Document the mathematical operation performed
- **Cite academic references** used during implementation (critical for clean room)
- Reference SPRAL routines consulted
- Include performance characteristics (complexity, memory usage)
- Document pivoting strategies and when they apply
- Provide usage examples for common workflows
- Cross-reference HSL routine names (MA27, MA57, MA97) for users familiar with HSL

## Current Implementation Status

**Completed:**
- ✅ Project scaffold and build configuration
- ✅ Error type definitions

**In Progress:**
- (Nothing yet - awaiting first implementation)

**Planned:**
- 📋 Sparse matrix validation and structure checking
- 📋 Elimination tree construction
- 📋 AMD ordering implementation
- 📋 Symbolic factorization
- 📋 Numeric LDLT factorization
- 📋 Solve phase (forward/backward substitution)

## Key Academic References

### Multifrontal Method
- **Duff, I.S. & Reid, J.K. (1983)**. "The Multifrontal Solution of Indefinite
  Sparse Symmetric Linear Systems". *ACM Transactions on Mathematical Software*,
  9(3):302-325. DOI: 10.1145/356044.356047

### Parallel Multifrontal
- **Hogg, J.D., Reid, J.K. & Scott, J.A. (2010)**. "Design of a Multicore Sparse
  Cholesky Factorization Using DAGs". *SIAM Journal on Scientific Computing*,
  32(6):3627-3649. DOI: 10.1137/090757216

### Ordering Algorithms
- **Amestoy, P.R., Davis, T.A. & Duff, I.S. (2004)**. "Algorithm 837: AMD, an
  Approximate Minimum Degree Ordering Algorithm". *ACM Transactions on Mathematical
  Software*, 30(3):381-388. DOI: 10.1145/1024074.1024081

### Pivoting Strategies
- **Ashcraft, C., Grimes, R.G. & Lewis, J.G. (1998)**. "Accurate Symmetric
  Indefinite Linear Equation Solvers". *SIAM Journal on Matrix Analysis and
  Applications*, 20(2):513-561. DOI: 10.1137/S0895479896296921

## Active Technologies
- Rust 1.87+ with sprs for sparse matrices
- faer (>= 0.22) for dense linear algebra within frontal matrices
- criterion for benchmarking
