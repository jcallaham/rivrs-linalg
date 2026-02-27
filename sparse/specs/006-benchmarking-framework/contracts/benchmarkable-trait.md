# Contract: Benchmarkable Trait

**Purpose**: Defines the interface any solver must implement to participate in benchmarking.

## Trait: `Benchmarkable`

The core abstraction allowing the benchmark harness to time solver operations without knowing solver internals.

### Methods

#### `bench_analyze`
- **Input**: Reference to a sparse symmetric matrix (CSC format)
- **Output**: `Option<Box<dyn Any>>` — opaque analysis result, or `None` if not implemented
- **Semantics**: Performs symbolic analysis (ordering, elimination tree, symbolic factorization). The returned value is passed to `bench_factor` if the full pipeline is being benchmarked.

#### `bench_factor`
- **Input**: Reference to a sparse symmetric matrix, optional analysis result from `bench_analyze`
- **Output**: `Option<Box<dyn Any>>` — opaque factorization result, or `None` if not implemented
- **Semantics**: Performs numeric factorization (LDL^T with APTP pivoting). If analysis result is `None`, the method should perform analysis internally.

#### `bench_solve`
- **Input**: Reference to a sparse symmetric matrix, factorization result from `bench_factor`, right-hand side vector
- **Output**: `Option<Vec<f64>>` — solution vector, or `None` if not implemented
- **Semantics**: Performs triangular solve using the factored form.

#### `bench_roundtrip`
- **Input**: Reference to a sparse symmetric matrix, right-hand side vector
- **Output**: `Option<Vec<f64>>` — solution vector, or `None` if not implemented
- **Semantics**: Performs complete analyze → factor → solve pipeline. Default implementation chains the three methods above.

### Design Rationale

- **`Option` return**: Signals whether a phase is implemented. `None` causes the benchmark harness to skip with a diagnostic message (FR-009).
- **`Box<dyn Any>` for intermediate state**: Avoids coupling the benchmark trait to specific solver types. The harness only needs to pass opaque state between phases; it never inspects it.
- **Separate from `SolverTest`**: `SolverTest` validates correctness (returns `TestResult` with metrics). `Benchmarkable` measures performance (returns computation results for Criterion to time). Different concerns, different traits.

### Relationship to `SolverTest`

| Aspect | `SolverTest` | `Benchmarkable` |
|--------|-------------|-----------------|
| Purpose | Correctness validation | Performance measurement |
| Output | `TestResult` with metrics | Opaque computation result |
| Timing | Not measured | Measured by Criterion |
| Used by | Test harness (`#[test]`) | Benchmark harness (`cargo bench`) |
| Input | `SolverTestCase` | Matrix + optional prior-phase output |

A solver implementation will typically implement both traits.
