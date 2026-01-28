# Implementation Plan: Sylvester Equation Solver

**Branch**: `001-sylvester-solver` | **Date**: 2026-01-28 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/001-sylvester-solver/spec.md`

## Summary

Implement a numerically stable solver for continuous-time Sylvester equations (AX + XB = C) and discrete-time variants (AXB + X = C) using the Bartels-Stewart algorithm based on Schur decomposition. The solver will be implemented in Rust using the faer linear algebra library, following clean room implementation practices to maintain MIT/Apache-2.0 licensing. Implementation will reference Golub & Van Loan Chapter 7.6.3, LAPACK's BSD-licensed dtrsyl routines, and academic literature, with validation against SLICOT test cases.

## Technical Context

**Language/Version**: Rust 1.75+ (MSRV to be determined during setup)
**Primary Dependencies**: faer (>= 0.19) for linear algebra, ndarray (>= 0.16) for array structures, approx for test comparisons
**Storage**: N/A (in-memory numerical computation library)
**Testing**: cargo test, criterion for benchmarks, comparison against SLICOT validation data
**Target Platform**: Cross-platform (Linux, macOS, Windows) with IEEE 754 floating-point arithmetic
**Project Type**: Single library project (scientific computing)
**Performance Goals**:
- 100×100 matrices in <100ms
- 500×500 matrices in <5s
- Within 2× of SLICOT Fortran performance
**Constraints**:
- Residual norms <10⁻¹² for well-conditioned double precision problems
- Memory usage <100MB overhead for 1000×1000 matrices
- Clean room implementation (no GPL code consultation)
**Scale/Scope**: Core library implementation, ~2000-3000 LOC for solver + tests

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### Clean Room Implementation (Principle I)

**Status**: ✅ PASS

**Compliance**:
- ✅ Specification explicitly prohibits consulting SLICOT Fortran source (slicot/src/*.f)
- ✅ Reference materials identified: G&VL Chapter 7.6.3, LAPACK BSD-licensed code, academic papers
- ✅ SLICOT usage limited to: documentation review, test case extraction, benchmark comparisons
- ✅ Dependencies section explicitly lists permissible sources

**Academic Sources Identified**:
- Golub & Van Loan "Matrix Computations" Chapter 7.6.3 (Algorithm 7.6.2)
- Bartels & Stewart (1972) "Solution of the Matrix Equation AX + XB = C" (CACM 15(9):820-826)
- LAPACK dtrsyl.f, strsyl.f (BSD-3-Clause)
- TOMS Algorithm 432 (validation only)

**Validation Approach**: Compare numerical outputs against SLICOT test cases without reading source code.

### Algorithm Documentation & Academic Attribution (Principle II)

**Status**: ✅ PASS

**Plan**:
- All public functions will include rustdoc with algorithm descriptions
- Citations to G&VL sections, LAPACK routines, and academic papers
- Cross-reference to SLICOT routine names (e.g., SB04MD, SB04ND) for user migration
- Document numerical properties: stability, complexity, failure modes

### Test-Driven Validation (Principle III)

**Status**: ✅ PASS

**Testing Strategy**:
- Unit tests: Small analytical test cases (2×2, 3×3 matrices with known solutions)
- Integration tests: G&VL textbook examples, LAPACK test suite cases
- Validation tests: SLICOT benchmark problems (extract test data from slicot/examples/)
- Edge cases: Singular matrices, ill-conditioned systems, dimension mismatches, NaN/Inf inputs
- Target coverage: 95% code coverage per SC-007

**Test Organization**:
```
tests/
├── unit/           # Small analytical cases
├── integration/    # Cross-routine tests, G&VL examples
├── validation/     # SLICOT benchmark comparisons
└── benches/        # Performance benchmarks
```

### Numerical Correctness & Stability (Principle IV)

**Status**: ✅ PASS

**Stability Measures**:
- Use Schur decomposition (structured, numerically stable) rather than eigendecomposition
- Solve triangular Sylvester systems via back-substitution (stable for well-separated eigenvalues)
- Implement condition number estimation based on eigenvalue separation (sep function)
- Avoid explicit matrix inversion
- Provide residual norm computation for solution verification

**Reference**: G&VL Section 7.6.4 discusses sensitivity and conditioning of Sylvester equations.

### Performance & Scalability (Principle V)

**Status**: ✅ PASS (with deferred optimization)

**Initial Implementation**:
- P1 (Priority 1): Correct implementation of Bartels-Stewart algorithm
- P4 (Priority 4): Performance optimization for large-scale problems

**Benchmark Suite**:
- Compare against SLICOT (via Python bindings if available)
- Compare against SciPy's solve_sylvester
- Matrix sizes: 10×10, 50×50, 100×100, 250×250, 500×500, 1000×1000

**Optimization Strategy**:
- Start with straightforward implementation using faer's Schur decomposition
- Profile and optimize bottlenecks if performance targets not met
- Consider Level-3 BLAS blocked algorithms (dtrsyl3 style) if needed for large matrices

### Code Quality & Rust Best Practices (Principle VI)

**Status**: ✅ PASS

**Rust Standards**:
- Generic over f32/f64 using faer's trait abstractions
- Result<T, E> for error handling (singular matrices, dimension mismatches, convergence failures)
- Comprehensive rustdoc comments with Examples, Errors, and Panics sections
- Follow Rust API guidelines
- In-place and allocating variants where appropriate

**Error Types**:
```rust
pub enum SylvesterError {
    DimensionMismatch { expected: (usize, usize), got: (usize, usize) },
    SingularMatrix { matrix: String },
    CommonEigenvalues { separation: f64 },
    InvalidInput { reason: String },
    ConvergenceFailure { iterations: usize },
}
```

### Constitution Check Summary

**Overall Status**: ✅ ALL GATES PASS

No violations of constitution principles. Implementation can proceed to Phase 0 (Research).

## Project Structure

### Documentation (this feature)

```text
specs/001-sylvester-solver/
├── spec.md              # Feature specification (complete)
├── plan.md              # This file (implementation plan)
├── research.md          # Phase 0: Algorithm research and design decisions
├── data-model.md        # Phase 1: Data structures and types
├── quickstart.md        # Phase 1: Usage examples and getting started
├── contracts/           # Phase 1: API contracts and public interfaces
│   └── sylvester.yaml   # Public API specification
├── checklists/          # Quality validation checklists
│   └── requirements.md  # Specification quality checklist (complete)
└── tasks.md             # Phase 2: Implementation tasks (created by /speckit.tasks)
```

### Source Code (repository root)

```text
src/
├── lib.rs              # Public API exports and module declarations
├── error.rs            # Error types (SylvesterError, Result type aliases)
├── sylvester/          # Sylvester equation solvers
│   ├── mod.rs          # Module exports
│   ├── continuous.rs   # Continuous-time solver (AX + XB = C)
│   ├── discrete.rs     # Discrete-time solver (AXB + X = C)
│   ├── triangular.rs   # Triangular Sylvester solver (core algorithm)
│   └── condition.rs    # Condition number estimation (sep function)
├── schur/              # Schur decomposition utilities (if not using faer's directly)
│   └── mod.rs
└── utils/              # Shared utilities
    ├── validation.rs   # Input validation, dimension checks
    └── residual.rs     # Residual norm computation

tests/
├── unit/
│   ├── test_triangular.rs      # Small triangular system tests
│   ├── test_continuous.rs      # Continuous-time analytical tests
│   ├── test_discrete.rs        # Discrete-time analytical tests
│   └── test_validation.rs      # Input validation tests
├── integration/
│   ├── test_gvl_examples.rs    # Golub & Van Loan textbook examples
│   └── test_lapack_suite.rs    # LAPACK test cases (if accessible)
└── validation/
    ├── test_slicot_continuous.rs   # SLICOT benchmark: continuous-time
    ├── test_slicot_discrete.rs     # SLICOT benchmark: discrete-time
    └── data/                        # Test matrices extracted from SLICOT examples

benches/
├── sylvester_continuous.rs     # Benchmark continuous-time solver
├── sylvester_discrete.rs       # Benchmark discrete-time solver
└── comparison.rs               # Comparison vs SciPy (if integration possible)

docs/
├── gvl-ch7.md                  # Golub & Van Loan Chapter 7 reference (already present)
└── 432.f                       # TOMS Algorithm 432 (validation reference, already present)

lapack/                         # LAPACK reference source (already present)
slicot/                         # SLICOT reference (already present, GPL - do not read src/)
```

**Structure Decision**: Single library project structure. This is a scientific computing library providing numerical routines, not an application with frontend/backend. The structure follows the constitution's guidance for code organization with clear separation of concerns: core algorithms (src/), validation tests (tests/validation/), and benchmarks (benches/).

## Complexity Tracking

No constitution violations requiring justification. All complexity is inherent to the numerical algorithm domain.

---

## Phase 0: Research & Algorithm Design

**Goal**: Resolve all technical unknowns and document algorithm design decisions.

### Research Tasks

1. **Bartels-Stewart Algorithm Deep Dive**
   - Study G&VL Algorithm 7.6.2 in detail
   - Understand the three main steps:
     1. Schur decomposition of A and B
     2. Transform C using Schur vectors
     3. Solve triangular Sylvester system via back-substitution
   - Document algorithmic steps in pseudocode
   - Note: G&VL Section 7.6.3 describes FZ - ZG = C; need to adapt for AX + XB = C

2. **LAPACK dtrsyl Implementation Review**
   - Read LAPACK dtrsyl.f header documentation
   - Understand parameter conventions (TRANA, TRANB, ISGN)
   - Note workspace requirements
   - Understand SCALE parameter for overflow prevention
   - Review dtrsyl3.f for Level-3 BLAS blocked variant (deferred to P4 optimization)

3. **faer Library Capabilities**
   - Confirm faer provides real Schur decomposition (QR algorithm)
   - Understand faer's matrix types: Mat<f64>, MatRef, MatMut
   - Check if faer provides triangular system solvers (trsm equivalent)
   - Understand workspace allocation patterns
   - Review faer documentation for generic abstractions over f32/f64

4. **Discrete-Time Variant**
   - Research transformation from AXB + X = C to standard form
   - G&VL and LAPACK may not directly cover discrete-time; check SciPy documentation
   - Possible approach: Transform to AXB^T - (-X)B^T = C or similar
   - Alternative: Implement as separate routine with modified back-substitution

5. **Condition Number Estimation**
   - Review G&VL Section 7.6.4 on sep(A, B) function
   - Understand sep(A, B) = min ||AX - XB||_F / ||X||_F
   - Research estimation techniques (don't need exact computation, estimate sufficient)
   - LAPACK dtrsyl provides SCALE but not sep estimate; may need separate routine

6. **Test Case Generation**
   - Extract test matrices from slicot/examples/ (identify which examples test Sylvester solvers)
   - Find G&VL examples or create small analytical cases
   - Generate random well-conditioned and ill-conditioned matrices
   - Create edge cases: singular, near-singular, dimension mismatches

### Design Decisions to Document in research.md

- **API Design**: Function signatures for public API
- **Generic vs Concrete**: Use faer's trait abstractions or concrete f64/f32 types?
- **In-place vs Allocating**: Provide both? Start with allocating for simplicity?
- **Workspace Management**: Pre-allocated workspace or internal allocation?
- **Error Handling**: Which conditions are Result::Err vs panic?
- **Discrete-Time Implementation**: Separate function or parameter flag?

**Output**: `research.md` with all decisions documented

---

## Phase 1: Data Model & API Contracts

**Prerequisites**: research.md complete

### Data Model (data-model.md)

Extract from spec and design decisions:

**Core Types**:
- `SylvesterSolver` - Main solver struct (may be stateless, just functions)
- Input matrices: A (n×n), B (m×m), C (n×m)
- Output matrix: X (n×m)
- `SylvesterError` - Error enumeration
- `SylvesterResult<T>` - Result type alias

**Constraints**:
- Matrix dimensions must satisfy: A.nrows() == A.ncols(), B.nrows() == B.ncols(), C.nrows() == A.nrows(), C.ncols() == B.ncols()
- Matrices must not contain NaN or Inf
- Eigenvalues of (A, -B) should not overlap (checked via sep estimate or Schur form inspection)

**Validation Rules**:
- Dimension compatibility checks on input
- NaN/Inf detection
- Eigenvalue separation warning if sep < threshold

### API Contracts (contracts/sylvester.yaml)

Generate OpenAPI-style specification for public Rust API:

```yaml
# Conceptual contract (actual implementation is Rust, not REST API)
solve_sylvester_continuous:
  inputs:
    a: Matrix<T> (n×n, square)
    b: Matrix<T> (m×m, square)
    c: Matrix<T> (n×m, rectangular)
  outputs:
    x: Result<Matrix<T>, SylvesterError>
  constraints:
    - a.is_square()
    - b.is_square()
    - c.nrows() == a.nrows()
    - c.ncols() == b.ncols()
  errors:
    - DimensionMismatch
    - InvalidInput (NaN/Inf)
    - CommonEigenvalues (sep too small)
  postconditions:
    - ||AX + XB - C||_F < ε

solve_sylvester_discrete:
  # Similar structure for AXB + X = C
  ...

compute_residual:
  inputs:
    a, b, c, x: Matrices
    equation_type: Continuous | Discrete
  outputs:
    residual_norm: f64
  description: "Computes ||AX + XB - C|| or ||AXB + X - C||"
```

### Quickstart Guide (quickstart.md)

Minimal working example:

```rust
use csrrs::sylvester::solve_continuous;
use faer::prelude::*;

// Define matrices A (2×2), B (2×2), C (2×2)
let a = mat![[1.0, 0.5], [0.0, -1.0]];
let b = mat![[2.0, 0.0], [0.0, 3.0]];
let c = mat![[1.0, 2.0], [3.0, 4.0]];

// Solve AX + XB = C
let x = solve_continuous(&a, &b, &c)?;

// Verify solution
let residual = compute_residual(&a, &b, &c, &x, EquationType::Continuous);
assert!(residual < 1e-10);
```

### Agent Context Update

Run `.specify/scripts/bash/update-agent-context.sh claude` to update Claude-specific context with:
- faer (linear algebra library)
- ndarray (array structures)
- Sylvester equation solver domain knowledge

**Output**: `data-model.md`, `contracts/sylvester.yaml`, `quickstart.md`, updated agent context

---

## Phase 2: Task Generation

**Prerequisites**: Phase 0 and Phase 1 complete

This phase is handled by the `/speckit.tasks` command, which will generate `tasks.md` with dependency-ordered implementation tasks based on the design from Phase 1.

Expected task categories:
1. Project setup (Cargo.toml, module structure)
2. Error type definitions
3. Input validation utilities
4. Triangular Sylvester solver (core algorithm)
5. Schur decomposition integration (using faer)
6. Continuous-time solver (high-level API)
7. Discrete-time solver (high-level API)
8. Residual computation utilities
9. Test suite implementation (unit, integration, validation)
10. Documentation and examples
11. Benchmark suite
12. Python bindings (deferred to future feature)

---

## Performance Optimization Strategy

### Two-Phase Approach: Correctness First, Then Speed

**Philosophy**: Implement correct, numerically stable solver first using unblocked (dtrsyl-style) algorithm. Optimize for performance second using blocked (dtrsyl3-style) algorithm after correctness is verified.

### Phase 1: Unblocked Algorithm (dtrsyl-style) - Priority 1-3

**Target**: Numerical correctness and stability
**Expected Performance**: Within 10× of SLICOT (acceptable for n < 500)

**Implementation**:
- Sequential block-by-block processing (1×1 and 2×2 blocks)
- Single-level SCALE factor logic (follows LAPACK dtrsyl pattern)
- Use faer's optimized primitives: Schur decomposition, matmul, triangular solvers
- Leverage faer's automatic SIMD for RealField operations

**Why Unblocked First**:
1. **Simpler algorithm**: Direct mapping to Bartels-Stewart (1972) and G&VL Algorithm 7.6.2
2. **Clean room clarity**: Fewer papers to reference, clearer provenance
3. **Easier verification**: Sequential logic easier to test and debug
4. **Sufficient for target domain**: Most control systems use n < 500

### Phase 2: Blocked Algorithm (dtrsyl3-style) - Priority 4

**Target**: Performance competitive with SLICOT for all matrix sizes
**Expected Performance**: Within 2× of SLICOT for n ≥ 500

**Implementation**:
- Block partitioning with dynamic block size (typically 16-64)
- Detect 2×2 diagonal blocks and adjust partition boundaries
- Two-level scaling: local SCALOC + global SCALE + consistency pass
- Level-3 BLAS updates using faer's matmul (cache-efficient blocked operations)
- Parallelize independent operations via faer's Par parameter

**Academic References**:
- Quintana-Orti & Van De Geijn (2003): "Formal derivation of algorithms: The triangular Sylvester equation," ACM TOMS 29:218-243
- Schwarz & Kjelgaard Mikkelsen (2020): "Robust Task-Parallel Solution of the Triangular Sylvester Equation," LNCS 12043:82-92
- LAPACK Working Note 107 and dtrsyl3.f implementation (BSD-3-Clause)

**Performance Gains**:
- n = 100: ~2× faster (blocked BLAS beats scalar operations)
- n = 500: ~10× faster (cache efficiency dominates)
- n = 1000: ~50× faster (O(n³) with better constants)

### Alignment with faer Performance Patterns

**Where We Get Performance Automatically**:
1. ✅ **Generic SIMD**: `T: RealField` enables automatic SIMD for f32/f64
2. ✅ **Optimized Schur**: faer's QR algorithm already parallel and blocked
3. ✅ **Optimized matmul**: faer's matrix multiplication cache-efficient
4. ✅ **Parallel Schur**: Use `Par::Rayon` for large decompositions

**Where We Must Optimize**:
- ⚠️ **Triangular Sylvester solve**: Custom algorithm we implement
  - Phase 1 (unblocked): Sequential, but uses faer's matmul for updates
  - Phase 2 (blocked): Use DGEMM-style blocked updates via faer's matmul

**faer Internal Code Study**:
- faer has private `lasy2` (2×2 Sylvester solver) in `real_schur.rs`
- faer has private `solve_sylvester_single_block` in `qz_real/mod.rs` (generalized form)
- Both are internal utilities, not public API
- We can study their patterns for numerical stability techniques

### API Design for Performance

**Workspace Control** (matches faer pattern):
```rust
// Simple API (internal allocation)
pub fn solve_continuous<T: RealField>(
    a: MatRef<T>, b: MatRef<T>, c: MatRef<T>
) -> Result<Mat<T>, SylvesterError>

// Advanced API (user workspace, zero-allocation)
pub fn solve_continuous_with_work<T: RealField>(
    a: MatRef<T>, b: MatRef<T>, c: MatRef<T>,
    stack: &mut MemStack,
) -> Result<Mat<T>, SylvesterError>

// Workspace query (like faer)
pub fn solve_continuous_req<T: Entity>(n: usize, m: usize) -> StackReq
```

**Algorithm Selection** (internal, transparent to user):
```rust
// Future: Automatically choose unblocked vs blocked based on size
fn select_algorithm(n: usize, m: usize) -> Algorithm {
    if n < 100 && m < 100 {
        Algorithm::Unblocked  // Simplicity over speed
    } else {
        Algorithm::Blocked    // Speed critical
    }
}
```

### Benchmark Strategy

**Matrix Sizes**: 10, 50, 100, 250, 500, 750, 1000
**Comparisons**: SLICOT (via Python bindings), SciPy, MATLAB (if available)
**Metrics**: Time, memory, numerical accuracy (residual norms)

**Acceptance Criteria**:
- Phase 1 (unblocked): Pass all tests, achieve 10⁻¹² accuracy, within 10× of SLICOT
- Phase 2 (blocked): Within 2× of SLICOT for all sizes, maintain accuracy

---

## Re-evaluation: Constitution Check (Post-Design)

**Status**: ✅ ALL GATES PASS (no changes from initial check)

The design maintains clean room compliance, follows test-driven validation principles, prioritizes numerical stability, and adheres to Rust best practices. Staged optimization approach (correctness first, speed second) aligns with constitution principles. No violations introduced during planning.

---

## Next Steps

1. ✅ Specification complete (`spec.md`)
2. ✅ Implementation plan complete (`plan.md`) - **YOU ARE HERE**
3. ⏭️ Run `/speckit.tasks` to generate dependency-ordered task list (`tasks.md`)
4. ⏭️ Begin implementation following tasks.md
5. ⏭️ Iterate with `/speckit.clarify` if ambiguities arise during implementation

---

**Plan Status**: Complete and ready for task generation
**Blocking Issues**: None
**Next Command**: `/speckit.tasks`
