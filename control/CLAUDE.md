# CLAUDE.md - Control Systems Algorithms

This file provides guidance to Claude Code when working on control systems algorithms in this directory.

## Project Overview

This directory contains control systems algorithm implementations for rivrs-linalg, focused on routines similar to those in SLICOT (Subroutine Library in Systems and Control Theory). The algorithms are implemented using modern Rust with `faer` for high-performance linear algebra operations.

**Parent Project**: rivrs-linalg - Numerical Linear Algebra for Rivrs
**Domain**: Control systems algorithms (Sylvester equations, Lyapunov equations, Riccati solvers, etc.)
**Current Status**: Phase 1 - Isolated development

## Licensing Strategy - Clean Room Implementation

**This is a clean room implementation to maintain licensing flexibility.**

SLICOT source code (both GPL-licensed slicot/ and BSD-3-licensed SLICOT-Reference/) was consulted during implementation:

- **SLICOT-Reference (BSD-3)**: May be freely consulted for implementation patterns
- **SLICOT GPL**: Documentation and test cases only - NEVER read source code (slicot/src/*.f files)
- **Primary sources**: Academic papers, control theory textbooks, BLAS/LAPACK source (BSD-licensed)
- **Always document** academic sources used for each implementation
- This approach allows rivrs-linalg to be released under Apache-2.0

## Reference Materials

Available in the parent `references/` directory:

- **faer-rs/**: Pure Rust linear algebra library - primary dependency for matrix operations
- **lapack/**: Reference LAPACK implementation (BSD-licensed) - consult freely
- **SLICOT-Reference/**: BSD-3-licensed SLICOT implementation - consult freely
- **slicot/**: GPL-licensed SLICOT (documentation and test cases only, NOT source code)

## Technology Stack

**Core Linear Algebra:**
- Use `faer` for high-performance dense linear algebra (decompositions, eigenvalues, matrix ops)
- Use `nalgebra` for Schur decomposition and other specialized algorithms
- Reference LAPACK and SLICOT-Reference for implementation patterns

**Future Python Bindings:**
- PyO3 for Rust-Python integration (not yet implemented)
- `ndarray` for zero-copy data sharing with NumPy

## Algorithm Categories

Following SLICOT's organization:

- **Sylvester equations** (current): AX + XB = C (continuous), AXB + X = C (discrete)
- **Lyapunov equations** (planned): Special cases of Sylvester equations
- **Riccati equations** (planned): Algebraic Riccati equation solvers
- **State-space analysis** (future): Canonical forms, controllability/observability
- **Model reduction** (future): Balancing and truncation methods

## Code Architecture

**Numerical Considerations:**
- Prioritize numerical stability - consult LAPACK and academic literature
- Use structured decompositions (Schur, QR, SVD) rather than direct inversion
- Implement workspace allocation patterns from faer
- Consider condition number estimation and backward error analysis

**API Design:**
- Provide both low-level routines and high-level composable APIs
- Use builder patterns for operations with many parameters
- Leverage Rust's type system for dimensional correctness
- Design for both batch and single-system workflows

**Error Handling:**
- Use Result types for numerically fallible operations
- Provide informative error messages
- Consider error recovery strategies

## Development Workflow

When implementing a new algorithm:

1. Review SLICOT-Reference/doc/ and slicot/doc/ to understand the algorithm's purpose
2. Identify relevant academic papers and textbooks
3. Review LAPACK source and SLICOT-Reference source for implementation patterns
4. **DO NOT examine GPL SLICOT source code** (slicot/src/*.f)
5. Design idiomatic Rust API leveraging faer
6. Implement using academic references and permissively-licensed code
7. Document academic references and LAPACK/SLICOT-Reference routines consulted
8. Add comprehensive tests using SLICOT test cases (slicot/examples/)
9. Benchmark against SLICOT
10. Create Python bindings (future)

## Git Commit Practices

- **Commit frequently**: After completing logical units (phase completion, tests passing, new module)
- **Commit after verification**: Always commit when tests pass or milestones are verified
- **Descriptive messages**: Use clear commit messages
- **Never accumulate**: Don't accumulate large amounts of uncommitted work

## Testing Strategy

- Unit tests for routines with known analytical solutions
- Integration tests comparing outputs to SLICOT examples
- Property-based tests for numerical stability
- Benchmark suite vs SLICOT and MATLAB Control Systems Toolbox
- Python integration tests (future)

## Documentation Standards

- Document the mathematical operation performed
- **Cite academic references** used during implementation (critical for clean room)
- Include references to standard control theory texts
- Provide usage examples
- Document numerical stability characteristics and failure modes
- Cross-reference SLICOT routine names for users migrating from Fortran/MATLAB

## Current Implementation Status

**Completed:**
- ✅ Continuous Sylvester solver (AX + XB = C) - Bartels-Stewart algorithm
- ✅ Discrete Sylvester solver (AXB + X = C) - Modified Bartels-Stewart
- ✅ Blocked Level-3 BLAS variant for continuous case
- ✅ Overflow prevention via scaling
- ✅ Eigenvalue separation estimation

**In Progress:**
- 🚧 Lyapunov equation solvers (special case of Sylvester)

**Planned:**
- 📋 Algebraic Riccati equation solvers
- 📋 Discrete-time Riccati solvers
- 📋 Generalized Sylvester and Lyapunov equations

## Active Technologies
- Rust 1.87+ with faer (>= 0.22) for linear algebra
- nalgebra (>= 0.34) for Schur decomposition
- approx for test comparisons (dev dependency)
- criterion for benchmarking (dev dependency)
