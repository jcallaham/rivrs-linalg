<!--
Sync Impact Report - Constitution v1.0.0
========================================

Version Change: Initial creation → v1.0.0

Principles Established:
- I. Clean Room Implementation (Licensing Integrity)
- II. Algorithm Documentation & Academic Attribution
- III. Test-Driven Validation (NON-NEGOTIABLE)
- IV. Numerical Correctness & Stability
- V. Performance & Scalability
- VI. Code Quality & Rust Best Practices

Templates Requiring Updates:
✅ plan-template.md - Constitution Check section aligns with new principles
✅ spec-template.md - Requirements section compatible with validation standards
✅ tasks-template.md - Task organization supports TDD and validation workflow
⚠ No command files exist yet in .specify/templates/commands/

Follow-up TODOs:
- None - all placeholders filled

Ratification: 2026-01-27 (initial adoption)
-->

# CSRRS Constitution

## Core Principles

### I. Clean Room Implementation (Licensing Integrity)

**NON-NEGOTIABLE**: To maintain licensing flexibility and avoid GPL virality, CSRRS MUST be implemented as a clean room project.

**Rules**:
- NEVER read, copy, or reference SLICOT Fortran source code (slicot/src/*.f files) during algorithm implementation
- NEVER examine GPL-licensed code during the development of equivalent algorithms
- ALWAYS implement algorithms from: academic papers, control theory textbooks, BLAS/LAPACK source (BSD-licensed), SLICOT-Reference (BSD-3 licensed), SciPy source (BSD-licensed), and published standards
- PERMITTED uses of SLICOT: documentation review (slicot/doc/), test case construction (slicot/examples/), and numerical accuracy criteria
- FREELY access SLICOT-Reference subdirectory - it is BSD-3 licensed and can be consulted without restriction
- Document ALL academic sources used for each implementation (papers, textbooks, standards, permissively-licensed code)

**Rationale**: GPL's copyleft provisions would force CSRRS to adopt GPL licensing, limiting adoption in commercial and permissively-licensed projects. Clean room implementation allows MIT/Apache-2.0 dual licensing while still benefiting from SLICOT's test suite and documentation.

**Enforcement**: Code reviews MUST verify that academic references are cited and no GPL source code was consulted. Pull requests without proper attribution will be rejected.

### II. Algorithm Documentation & Academic Attribution

Every algorithm implementation MUST include comprehensive documentation with academic attribution.

**Required Documentation**:
- Mathematical description of the algorithm's operation
- Academic sources consulted: specific papers (author, title, journal, year), textbook sections (author, title, edition, pages), or permissively-licensed reference implementations (project name, file path, license, commit hash)
- For direct ports from LAPACK, SciPy, or SLICOT-Reference: explicit attribution with original source file path and license
- Numerical stability characteristics and known failure modes
- Input/output specifications with dimension requirements
- Cross-reference to equivalent SLICOT routine names (for user migration from Fortran/MATLAB)

**Rationale**: Academic attribution establishes provenance for clean room implementation, aids future maintainers in understanding algorithm choices, and respects intellectual property of original researchers. It also helps users assess algorithm quality and make informed choices.

**Example**:
```rust
/// Computes the Schur decomposition of a real matrix using the QR algorithm.
///
/// # Algorithm
/// Implements the Francis double-shift QR algorithm as described in:
/// - Golub & Van Loan, "Matrix Computations" (4th Ed), Algorithm 7.5.2, pp. 378-382
/// - Based on LAPACK's DGEES routine (lapack/SRC/dgees.f, BSD-3-Clause)
///
/// # References
/// - G. H. Golub and C. F. Van Loan (2013). Matrix Computations, 4th ed. Johns Hopkins University Press.
/// - LAPACK Working Note 41: "A Parallel Algorithm for the Symmetric Eigenvalue Problem" (Anderson et al., 1992)
///
/// # SLICOT Equivalent
/// Comparable to SLICOT routine MB03QD (Schur form computation), but implemented independently from academic sources.
```

### III. Test-Driven Validation (NON-NEGOTIABLE)

All algorithms MUST be validated against known test cases before being considered complete.

**Validation Requirements**:
- MUST include unit tests with analytical solutions where available
- MUST include integration tests comparing outputs to reference implementations (SLICOT, MATLAB Control Systems Toolbox, SciPy, or published test cases)
- MUST verify numerical accuracy meets established criteria (e.g., relative error < 10⁻¹⁰ for double precision)
- MUST test edge cases: singular matrices, ill-conditioned systems, boundary conditions, empty inputs
- For direct ports from LAPACK/SciPy: MUST implement validation tests confirming equivalence

**Test Organization**:
- `tests/unit/` - Individual routine tests with analytical solutions
- `tests/integration/` - Cross-routine tests and reference implementation comparisons
- `tests/validation/` - SLICOT test case comparisons (use test data, not source)
- `benches/` - Performance benchmarks against SLICOT and MATLAB

**Rationale**: Scientific computing requires absolute correctness. Test-driven validation ensures algorithms match theoretical expectations and published results. Reference comparisons catch implementation errors that might not be obvious from code review alone.

**Red-Green-Refactor Cycle**:
1. Write tests that encode expected behavior from academic references or test cases
2. Verify tests FAIL with clear diagnostic output
3. Implement algorithm based on academic sources
4. Verify tests PASS
5. Refactor for clarity and performance while maintaining passing tests

### IV. Numerical Correctness & Stability

Implementations MUST prioritize numerical stability over naive implementations.

**Standards**:
- Use structured decompositions (Schur, QR, QZ, SVD) rather than direct matrix inversion
- Implement workspace allocation patterns from faer for stack-based temporary storage
- Avoid forming normal equations (A^T A) when A is ill-conditioned
- Consider condition number estimation and backward error analysis
- Use pivoting strategies to enhance stability (partial pivoting for LU, column pivoting for QR)
- For iterative algorithms: implement convergence checks and iteration limits
- Document numerical properties: complexity, stability guarantees, failure modes

**Rationale**: Control systems routines often deal with ill-conditioned matrices from discretization or high system order. Naive implementations can produce numerically meaningless results. Following established numerical linear algebra practices ensures robustness.

**Validation**: Tests MUST include ill-conditioned test cases and verify that algorithms either:
- Produce accurate results within tolerance, OR
- Fail gracefully with informative error messages suggesting remediation (regularization, preconditioning, etc.)

### V. Performance & Scalability

Implementations MUST be competitive with reference implementations while maintaining correctness.

**Performance Requirements**:
- Benchmark critical routines against SLICOT and MATLAB Control Systems Toolbox
- Target performance: within 2x of optimized Fortran for equivalent algorithms
- Use faer's SIMD capabilities through generic abstractions
- Leverage zero-cost abstractions (views, strided access) from faer and ndarray
- Implement parallelism for large-scale operations using faer's parallel support
- Profile before optimizing - use cargo-flamegraph or perf

**Scalability Considerations**:
- Design APIs for batch operations (multiple systems, parameter sweeps)
- Support both stack-based (small systems) and heap-based (large systems) allocation
- Consider memory layout for cache efficiency
- Document computational complexity (time and space)

**Rationale**: Rust's zero-cost abstractions enable performance comparable to Fortran while providing memory safety. Users migrating from MATLAB/SLICOT expect performance, not just correctness. However, premature optimization undermines maintainability - benchmark first.

**Benchmarking Standard**: Each major algorithm should include a benchmark comparing:
- CSRRS Rust implementation
- SLICOT via Python bindings (if available)
- MATLAB Control Systems Toolbox (where applicable)
- SciPy equivalent (where applicable)

### VI. Code Quality & Rust Best Practices

Code MUST follow Rust idioms and scientific computing best practices.

**Rust Standards**:
- Prefer generic implementations over f32/f64 duplicated code (use faer's trait abstractions)
- Use const generics for fixed-size systems where appropriate
- Implement Clone/Copy judiciously - large matrices should be borrowed, not copied
- Use in-place operations to minimize allocations (provide both in-place and allocating variants)
- Leverage type system to enforce correctness (e.g., state-space system dimensions, symmetric matrices)
- Consider no_std compatibility for embedded control applications
- Follow Rust API guidelines (https://rust-lang.github.io/api-guidelines/)

**Error Handling**:
- Use `Result<T, E>` for operations that can fail numerically (singular matrices, non-convergence, dimension mismatches)
- Provide informative error messages that guide users to solutions
- Consider error recovery strategies and document them (e.g., "Consider regularization if matrix is singular")
- Never panic in library code except for clear programmer errors (dimension mismatches caught at compile time when possible)

**Documentation**:
- Every public function MUST have rustdoc comments
- Include `# Examples` section for non-trivial functions
- Include `# Panics` section if function can panic
- Include `# Errors` section describing failure modes
- Add `# Safety` section for unsafe code with invariants

**Code Organization**:
- Separate concerns: core algorithms (`src/`), Python bindings (`python/`), benchmarks (`benches/`), validation (`tests/validation/`)
- Follow faer-rs workspace pattern: traits (`csrrs-traits`), core (`csrrs`), bindings (`csrrs-py`)
- Use modules to organize by algorithm category (state_space, frequency_domain, synthesis, reduction, etc.)

## Technical Standards

### Dependency Management

**Primary Dependencies**:
- `faer` (>= 0.19) - High-performance dense linear algebra
- `ndarray` (>= 0.16) - Array structure and Python interoperability
- `pyo3` (>= 0.22) - Python bindings
- `approx` - Floating-point comparison in tests
- `criterion` - Benchmarking

**Dependency Rules**:
- Minimize dependencies to reduce supply chain risk
- Prefer pure Rust implementations over FFI bindings
- All dependencies MUST have permissive licenses (MIT, Apache-2.0, BSD-2/3-Clause)
- Document rationale for each major dependency

### Python Bindings Strategy

Follow ndarray's pattern for zero-copy data sharing:
- Accept NumPy arrays as inputs via PyO3
- Return NumPy arrays as outputs
- Use ndarray as intermediate representation
- Provide Pythonic API with keyword arguments and default values
- Include type stubs (`.pyi` files) for IDE support
- Write docstrings following NumPy documentation standard

### Version Compatibility

- Support current stable Rust + 2 prior minor versions
- Support Python 3.9+ (following NumPy's support policy)
- Document MSRV (Minimum Supported Rust Version) in README and CI

## Development Workflow

### Implementation Process

When implementing a new control systems algorithm:

1. **Research Phase**:
   - Review slicot/doc/ HTML documentation to understand algorithm purpose and I/O
   - Identify relevant academic papers, textbooks, or standards describing the algorithm
   - Review LAPACK source (lapack/) or SciPy source if algorithm builds on standard linear algebra
   - Review SLICOT-Reference (BSD-3 licensed) freely as reference implementation
   - **DO NOT examine slicot/src/*.f files or other GPL-licensed source code**

2. **Design Phase**:
   - Design Rust API that is idiomatic and leverages faer's capabilities
   - Choose generic abstractions vs. concrete types
   - Decide on in-place vs. allocating variants
   - Document expected numerical properties

3. **Test Phase (BEFORE Implementation)**:
   - Extract test cases from SLICOT examples, academic papers, or create analytical test cases
   - Write unit tests and integration tests
   - Verify tests FAIL with clear diagnostics
   - Get user/reviewer approval on test coverage

4. **Implementation Phase**:
   - Implement using faer's matrix types and decompositions based on academic references
   - Document which academic references informed the implementation
   - For direct ports: include explicit attribution with file path and license
   - Write comprehensive rustdoc comments

5. **Validation Phase**:
   - Verify all tests PASS
   - Add additional edge case tests
   - Run benchmarks and compare to reference implementations
   - Document numerical accuracy achieved (e.g., "Achieves relative error < 10⁻¹² for well-conditioned matrices")

6. **Integration Phase**:
   - Create Python bindings if algorithm is user-facing
   - Write Python-level tests and documentation
   - Update examples and tutorials
   - Update CHANGELOG.md

### Code Review Requirements

Every PR MUST be reviewed for:

- [ ] **Clean Room Compliance**: No GPL source code consulted; academic references cited
- [ ] **Testing**: All tests pass; edge cases covered; numerical accuracy verified
- [ ] **Documentation**: Rustdoc complete; academic attribution present; examples included
- [ ] **Numerical Stability**: Algorithm uses stable methods; ill-conditioned cases handled
- [ ] **Performance**: Benchmarks show competitive performance (if applicable)
- [ ] **API Quality**: Follows Rust API guidelines; ergonomic for users

### Continuous Integration

CI MUST run:
- All unit, integration, and validation tests
- Benchmarks (informational, not blocking)
- `cargo clippy` with `-D warnings`
- `cargo fmt -- --check`
- `cargo doc --no-deps` (verify documentation builds)
- License compliance check (no GPL dependencies)

## Governance

### Amendment Process

This constitution can be amended through:

1. Proposal in GitHub issue or RFC document
2. Discussion with maintainers and community
3. PR updating constitution with rationale
4. Approval from 2+ maintainers
5. Update `LAST_AMENDED_DATE` and increment `CONSTITUTION_VERSION`

### Version Semantics

- **MAJOR**: Backward incompatible governance changes (e.g., removing a core principle, changing licensing strategy)
- **MINOR**: New principle added or materially expanded guidance (e.g., adding security requirements)
- **PATCH**: Clarifications, wording improvements, typo fixes

### Complexity Justification

When a PR violates simplicity principles (e.g., adds unnecessary abstraction, duplicates functionality), it MUST include justification:

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|-------------------------------------|
| [Specific complexity] | [Technical necessity] | [Why simpler approach insufficient] |

### Compliance Review

- All PRs MUST verify compliance with constitution principles
- Quarterly audits of codebase for compliance drift
- New contributors MUST read constitution before first contribution
- Violations require documented remediation plan

### Licensing

CSRRS is dual-licensed under MIT/Apache-2.0. All contributions MUST be compatible with this licensing model. Use of GPL or LGPL dependencies is PROHIBITED.

**Version**: 1.0.0 | **Ratified**: 2026-01-27 | **Last Amended**: 2026-01-27
