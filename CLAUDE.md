# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

rivrs-linalg (Numerical Linear Algebra for Rivrs) is a scientific computing library providing numerical linear algebra implementations for the Rivrs symbolic-numeric framework. Currently focused on control systems algorithms similar to those in SLICOT (Subroutine Library in Systems and Control Theory), with plans to expand to sparse solvers and other numerical methods. The project uses modern Rust with `ndarray` and `faer` for linear algebra operations and will provide Python bindings.

## Licensing Strategy - Clean Room Implementation

**This is a clean room implementation to maintain licensing flexibility.** SLICOT is GPL-licensed, which would force any derivative work to also be GPL. To avoid this:

- **NEVER read or copy SLICOT source code** (slicot/src/*.f files) during implementation
- Implement algorithms from: academic papers, control theory textbooks, BLAS/LAPACK source (BSD-licensed), and published standards
- Use SLICOT only for: documentation of what algorithms do, test cases with expected outputs, and numerical accuracy criteria
- Always document the academic sources used for each implementation
- This approach allows rivrs-linalg to be released under a permissive license

## Reference Materials

The repository contains three reference implementations in subdirectories:

- **faer-rs/**: Pure Rust linear algebra library (https://faer.veganb.tw) - primary dependency for matrix operations
- **lapack/**: Reference LAPACK implementation (BSD-licensed) - consult freely for understanding numerical algorithms and implementation patterns
- **slicot/**: Reference SLICOT Fortran implementation organized by categories (AB=state-space analysis, BB=benchmark problems, FB=frequency weighted balancing, IB=identification and balancing, MA=mathematical routines, MB=model reduction, SB=synthesis benchmarks, TB=transformation and block algorithms, etc.)

These directories should NOT be modified.

**CRITICAL - Clean Room Implementation:**
To maintain licensing flexibility, SLICOT source code (.f files in slicot/src/) must NEVER be consulted during implementation. SLICOT is GPL-licensed and examining the source during development could create copyright contamination. Instead:
- Use SLICOT only for non-copyrightable information: documentation (slicot/doc/), test cases (slicot/examples/), and expected numerical pass criteria
- Implement algorithms from BLAS/LAPACK source (permissive BSD license), academic papers, and control theory textbooks
- Document which academic references were used for each implementation

## Technology Stack

**Core Linear Algebra:**
- Use `faer` for high-performance dense linear algebra operations (decompositions, eigenvalues, matrix operations)
- Use `ndarray` for array structure and Python interoperability via PyO3
- Reference faer-rs workspace structure: note the separation of faer-traits, faer (core), faer-macros, and faer-ffi

**Python Bindings:**
- PyO3 for Rust-Python integration
- Follow ndarray's pattern for zero-copy data sharing between Rust and NumPy

## Code Architecture Principles

**Algorithm Categories (following SLICOT structure):**
- State-space analysis and transformations (canonical forms, controllability/observability)
- System interconnections (cascade, feedback, parallel)
- Model reduction and balancing
- Synthesis routines (Riccati equations, H-infinity control, LQG/LQR)
- Frequency domain analysis
- Time/frequency domain conversions
- Descriptor systems

**Numerical Considerations:**
- Prioritize numerical stability - consult LAPACK implementations and academic literature for proven algorithms
- Use structured decompositions (Schur, QR, SVD) rather than direct inversion where possible
- Implement workspace allocation patterns from faer for stack-based temporary storage
- Follow faer's approach to SIMD and performance optimization
- Consider condition number estimation and backward error analysis

**API Design:**
- Provide both low-level routines (à la SLICOT individual functions) and high-level composable APIs
- Use builder patterns for complex operations with many parameters
- Leverage Rust's type system to enforce physical correctness (e.g., state-space system dimensions)
- Design for both batch operations and single-system workflows

**Error Handling:**
- Use Result types for operations that can fail numerically (singular matrices, non-convergence)
- Provide informative error messages that guide users to solutions
- Consider error recovery strategies (e.g., regularization hints)

## Git Commit Practices

- **Commit frequently**: Create git commits after completing each logical unit of work (e.g., after finishing a phase, after getting tests passing, after adding a new module).
- **Commit after verification**: Always commit when tests pass or a milestone is verified working.
- **Descriptive messages**: Use clear commit messages that describe what was accomplished.
- **Never skip commits**: Do not accumulate large amounts of uncommitted work.

## Development Workflow

When implementing a new control systems algorithm:
1. Review slicot/doc/ HTML documentation to understand what the algorithm does and its purpose
2. Identify relevant academic papers, textbooks, or standards that describe the algorithm (e.g., IEEE papers, Laub, Skogestad & Postlethwaite, Golub & Van Loan, etc.)
3. Review LAPACK source code (lapack/) if the algorithm builds on standard linear algebra operations
4. **DO NOT examine SLICOT source code (slicot/src/*.f files)** - this would contaminate the clean room implementation
5. Design a Rust API that is idiomatic and leverages faer's capabilities
6. Implement using faer's matrix types and decompositions based on academic references and LAPACK patterns
7. Document which academic references and LAPACK routines informed the implementation
8. Add comprehensive tests comparing outputs against SLICOT test cases (slicot/examples/) - test data and pass criteria are non-copyrightable
9. Benchmark against SLICOT to verify performance is competitive
10. Create Python bindings that expose the functionality naturally to SciPy/NumPy users

## Rust Idioms for Scientific Computing

- Prefer generic implementations over f32/f64 duplicated code (use faer's trait abstractions)
- Use const generics for fixed-size systems where appropriate
- Leverage zero-cost abstractions (views, strided access) from faer and ndarray
- Implement Clone/Copy judiciously - large matrices should be borrowed
- Use in-place operations to minimize allocations
- Consider no_std compatibility for embedded control applications (see faer-no-std-test for patterns)

## Testing Strategy

- Unit tests for individual routines with known analytical solutions
- Integration tests comparing outputs to SLICOT examples
- Property-based tests for numerical stability properties
- Benchmark suite comparing performance to SLICOT and MATLAB Control Systems Toolbox
- Python integration tests ensuring bindings work correctly with NumPy arrays

## Performance Optimization

Reference faer's optimization patterns:
- Check faer/Cargo.toml for profile settings (note opt-level=3 for dev profile)
- Use faer's parallelism support for large-scale operations
- Profile before optimizing - use cargo-flamegraph or similar tools
- Consider adding faer-traits implementations for custom types
- Use faer's SIMD capabilities through generic abstractions

## Documentation Standards

- Every public function should document the mathematical operation it performs
- **Cite academic references used during implementation** (papers, textbooks, LAPACK routines) - this is critical for clean room documentation
- Include references to standard control theory texts (e.g., Skogestad & Postlethwaite, Zhou et al., Doyle et al.)
- Provide usage examples for common control systems workflows
- Document numerical stability characteristics and failure modes
- Cross-reference equivalent SLICOT routine names for users migrating from Fortran/MATLAB, but note that rivrs-linalg implementations are based on independent academic sources

## Active Technologies
- Rust 1.75+ (MSRV to be determined during setup) + faer (>= 0.19) for linear algebra, ndarray (>= 0.16) for array structures, approx for test comparisons (001-sylvester-solver)
- N/A (in-memory numerical computation library) (001-sylvester-solver)

## Recent Changes
- 001-sylvester-solver: Added Rust 1.75+ (MSRV to be determined during setup) + faer (>= 0.19) for linear algebra, ndarray (>= 0.16) for array structures, approx for test comparisons
