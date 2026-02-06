# Research: Core Test Infrastructure

**Branch**: `005-test-infrastructure` | **Date**: 2026-02-06

## R1: Module Visibility — Cargo Feature Flag Pattern

**Decision**: Use a `test-util` Cargo feature flag to gate test infrastructure in `src/`.

**Rationale**: `#[cfg(test)]` modules are only visible within the same crate (unit tests), not to integration tests in `tests/`. A feature flag allows both unit and integration tests to access harness types. The pattern `[dev-dependencies] rivrs-sparse = { path = ".", features = ["test-util"] }` ensures the test code is only compiled during `cargo test`, not in production builds.

**Alternatives considered**:
- `tests/common/mod.rs`: Already exists but can only be used by integration tests, not unit tests in `src/`. Also, types defined there can't be used in `src/` modules.
- Always-public module: Would pollute the production API with test-only types.

**Implementation**: Add `test-util = []` to `[features]` in Cargo.toml. Gate the new `testing` module with `#[cfg(feature = "test-util")]`.

## R2: faer Sparse Matrix Construction API

**Decision**: Use `SparseColMat::try_new_from_triplets()` for random matrix generation.

**Rationale**: This is faer's primary constructor for building sparse matrices from (row, col, value) triplets. It handles deduplication and CSC construction internally. The `Triplet::new(row, col, val)` helper creates entries.

**Key API surface**:
- `SparseColMat::try_new_from_triplets(nrows, ncols, &[Triplet<usize, usize, f64>]) -> Result<Self, CreationError>`
- faer stores full symmetric matrices (both triangles), not just lower/upper
- Symmetrization must be done manually: for each off-diagonal (i,j,v), also add (j,i,v)

**Alternatives considered**:
- Manual CSC construction from col_ptr/row_idx/values arrays: More error-prone, no benefit for test matrices.

## R3: Positive Definite Matrix Generation

**Decision**: Use diagonal dominance to guarantee positive definiteness.

**Rationale**: A real symmetric matrix with diagonal entries strictly greater than the sum of absolute values of off-diagonal entries in the same row is positive definite (Gershgorin circle theorem). This is cheap to compute and doesn't require eigenvalue decomposition.

**Algorithm**:
1. Generate random off-diagonal entries in lower triangle
2. Mirror to upper triangle for symmetry
3. Compute row sums of absolute off-diagonal values
4. Set diagonal = row_sum + margin (margin > 0)

**Alternatives considered**:
- A^T A construction: Always PD but produces denser matrices, harder to control sparsity pattern.
- Cholesky of a known L: Gives exact LDL^T but doesn't test ordering.

## R4: Indefinite Matrix Generation

**Decision**: Generate symmetric matrices with explicitly negative diagonal entries mixed with positive ones.

**Rationale**: Setting some diagonal entries negative while keeping others positive (with small off-diagonal entries) guarantees both positive and negative eigenvalues by Gershgorin. This is simpler and more controllable than random generation with hoped-for indefiniteness.

**Algorithm**:
1. Generate sparse symmetric structure
2. For first n/2 pivots: set positive diagonal dominance
3. For remaining pivots: set negative diagonal entries
4. Result has mixed inertia by construction

## R5: NumericalValidator Wrapping Pattern

**Decision**: Keep existing `validate::*` functions public. `NumericalValidator` wraps them with configurable thresholds and structured result reporting.

**Rationale**: The standalone functions (`reconstruction_error`, `backward_error`, `check_inertia`) are pure mathematical computations that callers may want directly. The `NumericalValidator` adds tolerance checking, pass/fail logic, and diagnostic context on top.

**Pattern**:
```
validate::reconstruction_error()  →  returns f64 (raw metric)
NumericalValidator::check_reconstruction()  →  returns ValidationResult { passed, metric, threshold, context }
```

## R6: SolverTest Trait Design

**Decision**: Trait with four methods returning `TestResult`. The `test_roundtrip` method is independent end-to-end, not a composition of the other three.

**Rationale**: Individual phase tests (`test_analyze`, `test_factor`, `test_solve`) isolate bugs to specific components. The `test_roundtrip` test validates the full pipeline as a single operation and checks the final mathematical result. Keeping them independent avoids coupling phase tests to the roundtrip.

**Design insight from faer**: faer's test suite uses reconstruction tests (A = LDL^T) as the primary oracle and random-RHS solve tests for backward error. Our `SolverTest` trait encodes these same patterns in a reusable interface.

## R7: Test Case Loading and Filtering

**Decision**: Build on existing `registry::load_registry()` and `registry::load_test_matrix()` with a new `load_test_cases(filter)` function.

**Rationale**: The registry already handles metadata parsing, path resolution, CI-subset fallback, and reference factorization loading. The new function adds filtering by (source, category, difficulty) and wraps results in `SolverTestCase`.

**Current pattern to supersede**: Both `hand_constructed.rs` and `suitesparse_ci.rs` manually call `load_registry()`, filter, then loop with `load_test_matrix()`. This is ~30 lines of boilerplate per test file.

## R8: Pattern Matrix Generators

**Decision**: Implement arrow, tridiagonal, and banded matrix generators as standalone functions returning `SparseColMat<usize, f64>`.

**Rationale**: These three patterns cover the most common structures in the hand-constructed test set and represent fundamentally different sparsity characteristics:
- **Arrow**: Dense first row/column, sparse remainder (tests ordering sensitivity)
- **Tridiagonal**: Minimal fill-in, easy to analyze (tests baseline correctness)
- **Banded**: Tunable fill-in via bandwidth parameter (tests scaling)

Each returns a matrix with known structural properties that tests can assert against.
