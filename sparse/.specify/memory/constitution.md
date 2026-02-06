<!--
Sync Impact Report - Constitution v1.1.0
========================================

Version Change: v1.0.0 → v1.1.0

Changes in v1.1.0:
- Principle I (Correctness First): Updated validation language from
  "SPRAL reference results" as prerequisite to reconstruction tests
  as primary oracle with SPRAL comparison as secondary layer
- Principle III (TDD): Reorganized test categories to establish
  reconstruction tests (P^T A P = L D L^T) as the primary correctness
  oracle. SPRAL comparison tests moved to secondary category, added
  when infrastructure is available. Added reconstruction error standard
  (< 10^-12). Backward error formula updated to scaled form.

Rationale: Phase 0.3 analysis revealed that SPRAL's C API does not
expose full factorization matrices (L, P), only aggregate statistics.
Reconstruction tests are mathematically stronger than cross-solver
comparison and require no external infrastructure.

Ratification: 2026-02-05 (initial adoption)
Amendment: 2026-02-06 (v1.1.0 - reconstruction test primacy)
-->

# rivrs-sparse Constitution

## Core Principles

### I. Correctness First (NON-NEGOTIABLE)

**This project prioritizes correctness above all other concerns.**
rivrs-sparse is a foundational building block for higher-level numerical
tools (optimization solvers, differential equation integrators, control
system design). An incorrect sparse solver that ships quickly will
silently poison every downstream consumer. A correct solver that ships
slowly is infinitely more valuable.

**Rules**:
- Every algorithm MUST produce mathematically correct results, verified
  against analytical solutions, reconstruction tests, and published
  test cases before being considered complete
- No implementation is "done" until it passes comprehensive validation:
  reconstruction tests (P^T A P = L D L^T), backward error checks,
  hand-constructed matrices with known factorizations, and SuiteSparse
  Matrix Collection test cases. SPRAL comparison is a secondary
  validation layer added when available, not a prerequisite.
- When correctness conflicts with performance, correctness wins.
  Performance optimization MUST NOT compromise numerical accuracy
- When correctness conflicts with schedule, correctness wins.
  There is no deadline that justifies shipping a solver that gives
  wrong answers
- Bugs in numerical output MUST be treated as P0 severity regardless
  of how rarely they occur or how small the error appears
- Every code path that affects numerical results MUST be covered by
  tests that verify correctness, not just absence of crashes

**Rationale**: A sparse solver used in interior point methods or
differential equation solvers will be called millions of times in
optimization loops. A subtle numerical error (wrong pivot decision,
incorrect inertia count, silent loss of precision) can cause upstream
solvers to converge to wrong solutions, diverge, or produce physically
meaningless results — with no indication that the sparse solver is the
root cause. The cost of debugging such failures vastly exceeds the cost
of rigorous upfront validation.

### II. Clean Room Implementation (Licensing Integrity)

**NON-NEGOTIABLE**: To maintain licensing flexibility and avoid
proprietary/GPL contamination, rivrs-sparse MUST be implemented as a
clean room project.

**Rules**:
- NEVER read, copy, or reference HSL source code (proprietary/
  restrictive license) during algorithm implementation
- NEVER examine GPL-licensed code during development of equivalent
  algorithms
- ALWAYS implement algorithms from: academic papers, textbooks, SPRAL
  source (BSD-3 licensed), LAPACK source (BSD-3 licensed), and faer
  source (MIT licensed)
- FREELY consult SPRAL source code — it is BSD-3 licensed and the
  primary reference implementation for SSIDS
- FREELY consult HSL *documentation* for algorithm descriptions only;
  NEVER consult HSL *source code*
- Document ALL academic sources used for each implementation (papers,
  textbooks, permissively-licensed code references)

**Rationale**: GPL copyleft or HSL's proprietary terms would force
rivrs-sparse to adopt restrictive licensing, preventing adoption in
commercial and permissively-licensed projects. Clean room implementation
from academic papers and BSD-licensed references allows Apache-2.0
licensing while producing equivalent algorithms.

**Enforcement**: Code reviews MUST verify that academic references are
cited and no restricted source code was consulted. Pull requests without
proper attribution will be rejected.

### III. Test-Driven Development (NON-NEGOTIABLE)

All algorithms MUST follow a strict test-first workflow. Tests encode
the specification; implementation follows.

**Red-Green-Refactor Cycle**:
1. Write tests that encode expected behavior from academic references,
   known factorizations, or SPRAL reference results
2. Verify tests FAIL with clear diagnostic output
3. Implement algorithm based on academic sources and permissive
   references
4. Verify tests PASS
5. Refactor for clarity and performance while maintaining passing tests

**Test Categories (all REQUIRED for solver components)**:
- **Reconstruction tests** (primary oracle): Verify that the
  factorization reconstructs the original matrix:
  `||P^T A P - L D L^T|| / ||A|| < epsilon`. This is the strongest
  possible correctness test — it proves mathematical correctness by
  definition, independent of any reference solver.
- **Unit tests**: Small hand-constructed matrices with analytically
  known results (e.g., 5x5 arrow matrix with known LDL^T factorization,
  inertia, and pivot structure)
- **Backward error tests**: Verify solution quality via
  `||Ax - b|| / (||A|| ||x|| + ||b||)` — computed independently,
  no reference solver needed
- **Property-based tests**: Verify structural invariants
  (correct inertia, symmetry preservation, permutation validity)
- **Edge case tests**: Singular matrices, zero pivots, maximally
  delayed pivots, ill-conditioned systems, empty inputs, 1x1 matrices
- **SPRAL comparison tests** (secondary, when available): Compare
  outputs against SPRAL (inertia, backward error, performance).
  Added when SPRAL build infrastructure is in place (deferred from
  Phase 0.3 to Phases 2-8). Not a prerequisite for correctness.
- **Regression tests**: Any bug fix MUST include a test that reproduces
  the bug before the fix and passes after

**Numerical Accuracy Standards**:
- Reconstruction error (||P^T A P - L D L^T|| / ||A||) MUST be
  < 10^-12 for double-precision problems
- Backward error (||Ax - b|| / (||A|| ||x|| + ||b||)) MUST be
  < 10^-10 for well-conditioned double-precision problems
- Inertia counts MUST be correct (verified analytically for
  hand-constructed matrices; cross-checked against SPRAL when available)
- Pivot decisions (1x1 vs 2x2) need not match SPRAL exactly, but
  the resulting factorization MUST satisfy accuracy requirements

**Rationale**: Scientific computing demands absolute correctness. TDD
ensures that expected behavior is encoded before implementation begins,
preventing "implementation-guided testing" where tests merely confirm
what was written rather than what should be true. For a foundational
library, this discipline catches errors at the cheapest possible point
in the development cycle.

### IV. Algorithm Documentation & Academic Attribution

Every algorithm implementation MUST include comprehensive documentation
tracing it to its academic origins.

**Required Documentation**:
- Mathematical description of the algorithm's operation
- Academic sources consulted: specific papers (author, title, journal,
  year), textbook sections (author, title, edition, pages), or
  permissively-licensed reference implementations (project name, file
  path, license, commit hash)
- For patterns adapted from SPRAL: explicit attribution with original
  source file path and BSD-3 license acknowledgment
- Numerical stability characteristics and known failure modes
- Input/output specifications with dimension and sparsity requirements
- Complexity analysis (time and space, in terms of matrix dimensions
  and nonzeros)
- Cross-reference to equivalent SPRAL function names for users
  migrating from Fortran

**Example**:
```rust
/// Performs A Posteriori Threshold Pivoting on a dense frontal matrix.
///
/// # Algorithm
/// Implements the APTP strategy as described in:
/// - Hogg, Duff, Lopez (2020), "A New Sparse LDL^T Solver Using
///   A Posteriori Threshold Pivoting", SIAM J. Sci. Comput. 42(4)
/// - Dense kernel pattern adapted from SPRAL ssids/cpu/kernels/ldlt_app.hxx
///   (BSD-3-Clause, commit abc1234)
///
/// # SPRAL Equivalent
/// Corresponds to the inner loop of spral_ssids_factor_cpu.
```

**Rationale**: Academic attribution establishes provenance for clean
room implementation, aids future maintainers in understanding algorithm
choices, and enables users to consult the original literature for
deeper understanding of behavior and limitations.

### V. Numerical Stability & Robustness

Implementations MUST prioritize numerical stability. Sparse indefinite
systems are inherently challenging; naive implementations will silently
produce garbage on real-world problems.

**Standards**:
- Use structured LDL^T decomposition with threshold pivoting — NEVER
  attempt to invert sparse matrices directly
- Implement APTP correctly: optimistic 1x1 pivots, a posteriori
  stability checks, fallback to 2x2 Bunch-Kaufman pivots, column
  delay when neither option is stable
- Monitor and report pivot statistics (number of 1x1 pivots, 2x2
  pivots, delayed columns) as factorization diagnostics
- Compute and report inertia (positive, negative, zero eigenvalue
  counts) — this is critical for interior point method consumers
- Implement overflow prevention for large entries and near-zero pivots
- Support both positive definite fast path (no pivoting needed) and
  full indefinite path (APTP) with runtime detection

**Failure Handling**:
- Structural singularity (detected at analysis phase) MUST return a
  descriptive error, not panic
- Numerical singularity (zero or near-zero pivots exceeding delay
  capacity) MUST return an error with diagnostic information
- Ill-conditioning MUST be detectable via factorization statistics,
  enabling callers to decide on remediation

**Rationale**: Sparse indefinite systems arise in interior point
optimization and constrained mechanics. These problems inherently have
saddle-point structure, meaning indefiniteness is the normal case, not
an edge case. The solver MUST handle this robustly or it has no
practical value.

### VI. Structured Development Discipline

Development MUST follow the phased plan in `docs/ssids-plan.md`. This
project is not on a deadline, and rushing phases creates compounding
technical debt in a domain where correctness is paramount.

**Rules**:
- Each phase MUST meet its documented exit criteria before proceeding
  to the next phase
- Phase 0 (literature review and test data) MUST be complete before
  any solver code is written — understanding the algorithm thoroughly
  prevents costly rework
- New phases MUST NOT be started while the current phase has failing
  tests or unresolved correctness issues
- Shortcuts that sacrifice correctness for velocity are PROHIBITED.
  If a phase is taking longer than expected, the response is to
  understand why, not to skip validation steps
- Progress is tracked in `docs/ssids-log.md`; the plan document
  (`docs/ssids-plan.md`) is updated when decisions change the approach
- Commit frequently after logical units; never accumulate large amounts
  of unverified work

**Rationale**: Sparse direct solvers are complex multi-component systems
where later phases depend critically on earlier ones being correct. A
bug in symbolic analysis corrupts every factorization. A bug in dense
APTP kernels corrupts every solve. Phased development with exit criteria
ensures each foundation layer is solid before building on it.

### VII. Code Quality & Rust Best Practices

Code MUST follow Rust idioms and scientific computing best practices.

**Rust Standards**:
- Use `faer` types and patterns as the foundation (CSC storage,
  elimination trees, workspace management, SIMD abstractions)
- Prefer generic implementations via faer's trait abstractions over
  duplicated f32/f64 code
- Use `Result<T, E>` for all fallible operations (singular matrices,
  non-convergence, dimension mismatches)
- Never panic in library code except for clear programmer errors
  (debug assertions for internal invariants are acceptable)
- Provide informative error messages that guide users to solutions
- Use the type system to enforce correctness where possible (e.g.,
  distinct types for symbolic vs. numeric factorization results)
- Follow Rust API guidelines for public API design

**API Design**:
- Three-phase API: `analyze → factorize → solve`
- Symbolic analysis result is reusable across multiple factorizations
  with the same sparsity pattern
- Builder pattern for solver configuration (pivot threshold, ordering
  strategy, parallelism settings)
- In-place operations to minimize allocations; provide both in-place
  and allocating variants for key operations

**Error Handling**:
- Distinguish structural singularity (analysis phase) from numerical
  singularity (factorization phase) via distinct error types
- Include diagnostic context in errors (matrix dimensions, pivot
  statistics, phase where failure occurred)
- Document all error conditions in rustdoc `# Errors` sections

**Documentation**:
- Every public function MUST have rustdoc comments
- Include `# Examples` section for non-trivial functions
- Include `# Panics` section if function can panic
- Include `# Safety` section for any unsafe code with invariants

## Technical Standards

### Dependency Management

**Primary Dependencies**:
- `faer` (>= 0.22) — High-performance linear algebra and sparse
  infrastructure (CSC, elimination trees, AMD ordering)
- `approx` — Floating-point comparison in tests (dev dependency)
- `criterion` — Benchmarking (dev dependency)

**Future Dependencies** (add when needed, not before):
- `rayon` — Parallel factorization (when supernodal phase begins)
- `pyo3` — Python bindings (when core solver is validated)

**Dependency Rules**:
- Minimize dependencies to reduce supply chain risk
- Prefer pure Rust implementations over FFI bindings
- All dependencies MUST have permissive licenses (MIT, Apache-2.0,
  BSD-2/3-Clause)
- Document rationale for each major dependency addition

### Performance Standards

Performance optimization is welcome but MUST NOT compromise correctness.

**Guidelines**:
- Benchmark against SPRAL reference results before and after
  optimization
- Target: within 2x of SPRAL sequential performance for equivalent
  algorithms (aspirational, not blocking)
- Profile before optimizing — use cargo-flamegraph or perf
- Document computational complexity for all algorithms
- Performance regressions > 10% require investigation

**Benchmarking**:
- Use criterion.rs for reproducible benchmarks
- Benchmark suite covers: symbolic analysis, dense APTP kernel,
  full solve pipeline
- Track performance over commits via CI

## Development Workflow

### Implementation Process

When implementing a new solver component:

1. **Literature Phase**: Review relevant academic papers, SPRAL source
   (BSD-3), and faer source. Document algorithm understanding.
   DO NOT examine HSL source code.

2. **Test Phase (BEFORE Implementation)**: Write tests encoding
   expected behavior from references, known solutions, and SPRAL
   golden results. Verify tests FAIL.

3. **Implementation Phase**: Implement from academic references and
   permissive code. Document all sources consulted. Write comprehensive
   rustdoc comments.

4. **Validation Phase**: Verify all tests PASS. Add edge case and
   stress tests. Run benchmarks. Document numerical accuracy achieved.

5. **Review Phase**: Verify clean room compliance, test coverage,
   documentation completeness, and numerical stability.

### Code Review Checklist

Every PR MUST be reviewed for:

- [ ] **Correctness**: All tests pass; numerical accuracy meets
  standards; no silent failure modes
- [ ] **Clean Room Compliance**: No restricted source code consulted;
  academic references cited
- [ ] **Testing**: TDD workflow followed; edge cases covered; reference
  comparisons included
- [ ] **Documentation**: Rustdoc complete; academic attribution present;
  complexity documented
- [ ] **Numerical Stability**: Algorithm uses stable methods;
  ill-conditioned cases handled or diagnosed
- [ ] **API Quality**: Follows Rust API guidelines; error types are
  informative; types enforce invariants

## Governance

### Amendment Process

This constitution can be amended through:

1. Proposal documenting the change and rationale
2. Discussion with maintainers
3. PR updating constitution with rationale
4. Update `LAST_AMENDED_DATE` and increment `CONSTITUTION_VERSION`

### Version Semantics

- **MAJOR**: Backward incompatible governance changes (e.g., removing a
  core principle, changing licensing strategy)
- **MINOR**: New principle added or materially expanded guidance
- **PATCH**: Clarifications, wording improvements, typo fixes

### Compliance Review

- All PRs MUST verify compliance with constitution principles
- New contributors MUST read this constitution before first contribution
- Violations of NON-NEGOTIABLE principles (I, II, III) require
  immediate remediation before merge

### Licensing

rivrs-sparse is licensed under Apache-2.0. All contributions MUST be
compatible with this licensing model. Use of GPL, LGPL, or proprietary
dependencies is PROHIBITED.

**Version**: 1.1.0 | **Ratified**: 2026-02-05 | **Last Amended**: 2026-02-06
