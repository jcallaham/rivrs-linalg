# Sparse Symmetric Indefinite Solver for Rust
## Project Specification and Development Phases

---

## Project Overview

**Objective:** Port SPRAL's SSIDS (Sparse Symmetric Indefinite Direct Solver) to Rust, leveraging faer for dense linear algebra operations and modern Rust ecosystem tooling.

**Target Applications:** Interior point optimization methods, with consideration for differential equation solvers.

**Key Technical Goals:**
- Implement A Posteriori Threshold Pivoting (APTP) algorithm
- Achieve numerical correctness on hard indefinite problems
- Match or exceed SPRAL's performance (sequential and shared-memory parallel)
- Provide safe, idiomatic Rust API
- Build on faer ecosystem and design patterns

**Architecture Strategy:**
- **Leverage faer infrastructure** (~70% reuse): CSC storage, elimination trees, AMD ordering, permutation utilities, workspace management, supernodal symbolic analysis
- **Build APTP-specific components** (~30% new): 2x2 pivot logic, mixed diagonal storage, multifrontal APTP factorization kernel
- **Multifrontal from the start**: APTP is inherently a blocked algorithm operating on dense frontal matrices — a simplicial (column-by-column) implementation would be throwaway code. faer already provides simplicial LDL^T as a fallback for positive-definite or mildly indefinite problems. Our value-add is the APTP kernel, which requires frontal matrices. Phases 2-7 build the multifrontal solver directly; Phase 8 adds parallelism and performance optimization.
- **Clean room implementation**: All code derived from BSD-licensed references (LAPACK, SLICOT-Reference, SPRAL) and academic papers

**Development Timeline:**
- Phases 0-1: Foundation, infrastructure, and tooling
- Phases 2-5: APTP data structures, symbolic analysis, ordering, dense kernel
- Phases 6-7: Multifrontal factorization, solve & end-to-end integration
- Phase 8: Parallelization and performance optimization
- Phase 9: Polish and release

---

## Phase 0: Foundation & Literature Review (**COMPLETE**)

### Objectives
Establish comprehensive understanding of algorithms, gather reference implementations, and compile authoritative test data before writing any solver code.

### Deliverables

#### 0.1: Literature Review & Reference Library (**COMPLETE**)
**Task:** Build a complete technical reference library

**Documents to compile:**
1. **Core APTP Papers:**
   - Hogg, Duff, Lopez (2020) - "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting" (SIAM J. Sci. Comput.)
   - Hogg, Ovtchinnikov, Scott (2016) - "A sparse symmetric indefinite direct solver for GPU architectures" (ACM TOMS)
   - RAL Technical Reports on SSIDS version 2.0

2. **Multifrontal Method Foundations:**
   - Duff & Reid (1983) - "The multifrontal solution of indefinite sparse symmetric linear equations"
   - Liu (1992) - "The Multifrontal Method for Sparse Matrix Solution: Theory and Practice" (SIAM Review)
   - Davis survey paper on direct methods

3. **Pivoting Strategies:**
   - Schenk & Gärtner (2006) - "On fast factorization pivoting methods for sparse symmetric indefinite systems"
   - Duff & Pralet (2005) - "Strategies for scaling and pivoting for sparse symmetric indefinite problems"

4. **Ordering & Analysis:**
   - George & Liu - Nested dissection papers
   - Gilbert, Ng, Peyton - Symbolic factorization methods
   - Duff & Koster - Matching-based orderings (MC64)

5. **faer Documentation:**
   - faer-core design documentation
   - faer-sparse implementation details
   - JOSS paper on faer
   - Review existing LBLT/Bunch-Kaufman implementation

**Deliverable Format:**
```
/docs/references/
  ├── papers/
  │   ├── aptp/
  │   ├── multifrontal/
  │   ├── pivoting/
  │   └── ordering/
  ├── algorithms/
  │   ├── APTP-ALGORITHM.md
  │   ├── MULTIFRONTAL-METHOD.md
  │   ├── SYMBOLIC-ANALYSIS.md
  │   └── ORDERING-STRATEGIES.md
  └── notes/
      ├── SPRAL-CODE-REVIEW.md
      └── FAER-INTEGRATION-NOTES.md
```

**Success Criteria:**
- [x] All key papers obtained and organized
- [x] Algorithm pseudocode extracted and documented
- [x] SPRAL source code reviewed with annotations
- [x] faer integration points identified and documented
- [x] Team member can explain APTP algorithm from memory

#### 0.2: Test Matrix Collection Assembly
**Task:** Build comprehensive test suite from all available sources

**Sources:**

1. **SuiteSparse Matrix Collection:**
   - Query for symmetric indefinite matrices
   - Download "easy indefinite" set (~30 matrices)
   - Download "hard indefinite" set (~15 matrices)
   - Include positive definite set for validation (~20 matrices)

2. **SPRAL Repository:**
   - Extract test matrices from SPRAL tests
   - Document test cases from examples/
   - Capture any synthetic test generators

3. **Academic Paper Test Sets:**
   - Matrices specifically mentioned in APTP papers
   - Test cases from HSL MA57/MA86/MA97 papers (if available)
   - Interior point method test problems

4. **Hand-Constructed Matrices:**
   - Small matrices with known factorizations (5×5 to 20×20)
   - Matrices with specific properties:
     - Arrow matrices (known supernode structure)
     - Block diagonal (independent subsystems)
     - Tridiagonal (simple elimination tree)
   - Matrices designed to stress specific components:
     - Maximum delayed pivots
     - Worst-case fill-in
     - Ill-conditioned problems

**Deliverable Format:**
```
/test-data/
  ├── metadata.json          # Index of all matrices
  ├── hand-constructed/
  │   ├── arrow-10.mtx
  │   ├── tridiagonal-20.mtx
  │   └── ...
  ├── suitesparse/
  │   ├── easy-indefinite/
  │   ├── hard-indefinite/
  │   └── positive-definite/
  ├── spral-tests/
  └── interior-point/
```

**Metadata Schema:**
```json
{
  "name": "GHS_indef/aug3dcqp",
  "source": "SuiteSparse",
  "size": 35543,
  "nnz": 77829,
  "properties": {
    "symmetric": true,
    "indefinite": true,
    "difficulty": "hard",
    "expected_delays": "high"
  },
  "reference_results": {
    "spral_residual": 1.2e-14,
    "spral_inertia": [12000, 23543, 0],
    "spral_time_analyze": 0.15,
    "spral_time_factor": 2.34
  }
}
```

**Success Criteria:**
- [x] Minimum 70 test matrices collected (82 total: 15 hand-constructed + 67 SuiteSparse)
- [x] All matrices in standard Matrix Market format
- [x] Complete metadata for each matrix (metadata.json with properties, paper references)
- [x] Range from 1×1 to 1.6M×1.6M problems (4+ orders of magnitude)
- [x] Known "killer" cases identified (stokes128, ncvxqp3, c-big, c-71, plus hand-constructed stress tests)

#### 0.3: SPRAL Golden Results Generation (**DEFERRED**)
**Task:** ~~Run SPRAL on all test matrices to create reference outputs~~

**Status:** Deferred. See `specs/003-spral-golden-results/decision.md` for full rationale.

**Summary of Decision:**
After investigating SPRAL's C API, we determined that the golden results
infrastructure (building SPRAL from source, writing a C driver, running on
82 matrices) provides limited incremental value over what we can validate
independently:

- SPRAL does not expose the L factor, permutation P, or elimination tree
  structure — only aggregate statistics and the block diagonal D
- The **reconstruction test** (`||P^T A P - L D L^T|| / ||A|| < epsilon`)
  is a strictly stronger correctness oracle than comparing against SPRAL's output
- Backward error (`||Ax - b|| / (||A|| ||x|| + ||b||)`) is computed from
  our own solution, no reference needed
- Phase 0.2's hand-constructed matrices already have analytically known
  factorizations and inertia
- Building SPRAL on arm64 carries non-trivial infrastructure risk for
  limited return at this stage

**Deferred to:** When the Rust solver can factorize SuiteSparse matrices
(Phases 2-7), for performance benchmarking and inertia validation on large
matrices where analytical verification is impractical.

**Original Success Criteria (deferred, not abandoned):**
- [ ] SPRAL successfully factors all test matrices
- [ ] Complete results captured for each matrix
- [ ] Results are reproducible (run twice, verify identical)
- [ ] Summary statistics generated (median times, failure modes)
- [ ] Challenging cases documented

#### 0.4: Initial Repository Setup (**COMPLETE**)
**Task:** Establish project structure and development environment

**Actual Structure** (evolved from plan to match monorepo layout):
```
sparse/
├── Cargo.toml              # faer 0.22, serde, serde_json
├── src/
│   ├── lib.rs              # pub mod io, validate, error
│   ├── error.rs            # SparseError with IoError, ParseError
│   ├── io.rs               # pub mod mtx, reference, registry
│   ├── io/mtx.rs           # Matrix Market parser
│   ├── io/reference.rs     # JSON factorization loader
│   ├── io/registry.rs      # Test matrix catalog (metadata.json)
│   └── validate.rs         # reconstruction_error, backward_error, check_inertia
├── tests/
│   ├── hand_constructed.rs # 15-matrix integration test with reconstruction validation
│   └── suitesparse_ci.rs  # 10-matrix CI-subset integration test
├── benches/
│   └── matrix_loading.rs  # Criterion benchmarks
├── test-data/              # From Phase 0.2
└── .github/workflows/ci.yml  # test-sparse, lint-sparse, doc-sparse jobs
```

**Dependencies:**
```toml
[dependencies]
faer = "0.22"
serde = { version = "1", features = ["derive"] }
serde_json = "1"

[dev-dependencies]
criterion = "0.5"
approx = "0.5"
rand = "0.8"
rand_distr = "0.4"
```

**Success Criteria:**
- [x] Repository structure established (IO modules, validation, error types)
- [x] CI pipeline configured (test + clippy + fmt + doc for sparse domain)
- [x] All test matrices accessible from tests (15 hand-constructed, 10 CI-subset)
- [x] Documentation framework in place (rustdoc with -D warnings)
- [x] Reconstruction error < 10^-12 for all 15 hand-constructed matrices (SC-008)

### Phase 0 Exit Criteria

**Required Outcomes:**
1. Complete understanding of APTP algorithm documented
2. 70+ test matrices collected with metadata and difficulty classification
3. Hand-constructed matrices have analytically known factorizations (L, D, P, inertia)
4. Correctness validation strategy defined: reconstruction tests as primary oracle,
   backward error as secondary, SPRAL comparison deferred to Phases 2-8
5. Development environment configured
6. Team can articulate:
   - How APTP differs from TPP
   - Why multifrontal method is used
   - What makes indefinite problems "hard"

**Correctness Validation Hierarchy** (established in Phase 0.3 decision):
1. **Reconstruction test** (primary): `||P^T A P - L D L^T|| / ||A|| < epsilon`
   — Proves mathematical correctness by definition
2. **Backward error** (secondary): `||Ax - b|| / (||A|| ||x|| + ||b||) < 10^-10`
   — Validates the full solve pipeline
3. **Analytical verification**: Inertia and factorization match known results
   for hand-constructed matrices
4. **SPRAL comparison** (deferred): Performance benchmarking and inertia
   cross-validation on large SuiteSparse matrices — added when Rust solver
   can process those matrices

**Validation Questions:**
- Can we explain the APTP algorithm to someone unfamiliar with it?
- Do we have test cases that will expose errors in each component?
- Is our validation strategy (reconstruction + backward error) sufficient
  to catch all categories of solver bugs?

**Time Investment:** This phase should not be rushed. Better to spend extra time here than discover gaps later.

---

## Phase 1: Testing & Benchmarking Infrastructure (**COMPLETE**)

### Objectives
Build production-quality testing and benchmarking framework before implementing any solver components. This infrastructure will be used throughout all subsequent phases.

### Deliverables

#### 1.1: Core Test Infrastructure (**COMPLETE**)
**Task:** Create reusable testing framework

**Test Harness Design:**

```rust
// sparse-ldlt-test/src/lib.rs

/// Represents a complete test case for the solver
pub struct SolverTestCase {
    pub name: String,
    pub matrix: SparseMatrix,
    pub properties: MatrixProperties,
    pub reference: Option<ReferenceResults>,
}

pub struct MatrixProperties {
    pub size: usize,
    pub nnz: usize,
    pub symmetric: bool,
    pub definite: Option<bool>,
    pub condition_number: Option<f64>,
}

pub struct ReferenceResults {
    pub source: ResultSource,  // SPRAL, Exact, etc.
    pub inertia: Inertia,
    pub solution: Option<Vec<f64>>,
    pub residual_norm: f64,
}

/// Core test operations
pub trait SolverTest {
    fn test_analyze(&self, case: &SolverTestCase) -> TestResult;
    fn test_factor(&self, case: &SolverTestCase) -> TestResult;
    fn test_solve(&self, case: &SolverTestCase) -> TestResult;
    fn test_roundtrip(&self, case: &SolverTestCase) -> TestResult;
}

/// Numerical validation
pub struct NumericalValidator {
    rtol: f64,
    atol: f64,
}

impl NumericalValidator {
    pub fn check_residual(&self, a: &Matrix, x: &[f64], b: &[f64])
        -> ValidationResult;

    pub fn check_inertia(&self, computed: Inertia, expected: Inertia)
        -> ValidationResult;

    pub fn check_forward_error(&self, x: &[f64], x_ref: &[f64])
        -> ValidationResult;
}
```

**Test Utilities:**

```rust
// Test matrix generators
pub fn load_test_matrix(name: &str) -> SolverTestCase;
pub fn generate_random_matrix(n: usize, nnz: usize,
                              props: MatrixProperties) -> SparseMatrix;

// Comparison tools
pub fn compare_with_spral(result: &SolverResult,
                         reference: &SPRALReference) -> Comparison;

// Property-based testing helpers
pub fn arbitrary_sparse_symmetric() -> impl Strategy<Value = SparseMatrix>;
```

**Success Criteria:**
- [ ] Can load any test matrix in <100ms
- [ ] Numerical validation with configurable tolerances
- [ ] Clear, actionable error messages
- [ ] Property-based test generators for random matrices
- [ ] Integration with Rust test framework

#### 1.2: Benchmarking Framework (**COMPLETE**)
**Task:** Create consistent performance measurement infrastructure

**Benchmark Suite Design:**

> **Note:** The notional API below was the original design sketch. The actual
> implementation uses a `Benchmarkable` trait + Criterion integration.
> `comparison_solvers` and `parallel_configs` are deferred until the solver
> supports parallelism and comparison benchmarks become meaningful (Phases 2-8).

```rust
// sparse-ldlt-bench/src/lib.rs (notional — see src/benchmarking/ for actual API)

pub struct BenchmarkConfig {
    pub matrices: Vec<String>,
    pub operations: Vec<Operation>,
    pub comparison_solvers: Vec<Solver>,      // Deferred to Phases 2-8
    pub parallel_configs: Vec<ParallelConfig>, // Deferred to parallel implementation
}

pub enum Operation {
    Analyze,
    Factor,
    Solve,
    FullRoundtrip,
}

pub enum Solver {
    SparseLDLT,
    SPRAL,
    // Future: MUMPS, PARDISO, etc.
}

pub struct BenchmarkResult {
    pub matrix: String,
    pub operation: Operation,
    pub solver: Solver,
    pub config: ParallelConfig,
    pub time: Duration,
    pub memory_peak: usize,
    pub metadata: HashMap<String, Value>,
}

/// Run comprehensive benchmark suite
pub fn benchmark_suite(config: BenchmarkConfig) -> BenchmarkResults;

/// Generate comparison report
pub fn generate_report(results: &BenchmarkResults) -> Report;
```

**Benchmarks to Implement:**

1. **Component Benchmarks:**
   - Symbolic analysis time vs matrix size
   - Dense APTP performance vs block size
   - Tree traversal overhead
   - Memory allocation patterns

2. **End-to-End Benchmarks:**
   - Full solve time on test suite
   - Scaling with number of threads
   - Memory usage vs matrix size

3. **Comparison Benchmarks** *(deferred to Phases 2-8)*:
   - SparseLDLT vs SPRAL on same matrices
   - Sequential vs parallel performance
   - Different ordering strategies

**Criterion.rs Integration:**

```rust
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

fn bench_analyze(c: &mut Criterion) {
    let mut group = c.benchmark_group("analyze");

    for matrix_name in ["small", "medium", "large"] {
        let matrix = load_test_matrix(matrix_name);

        group.bench_with_input(
            BenchmarkId::new("SparseLDLT", matrix_name),
            &matrix,
            |b, m| b.iter(|| analyze(m))
        );

        group.bench_with_input(
            BenchmarkId::new("SPRAL", matrix_name),
            &matrix,
            |b, m| b.iter(|| spral_analyze(m))
        );
    }

    group.finish();
}
```

**Output Format:**
- Criterion's HTML reports for detailed analysis
- CSV exports for tracking over time
- Markdown tables for documentation
- Automated regression detection

**Success Criteria:**
- [ ] Can benchmark any solver component
- [ ] Direct comparison with SPRAL possible
- [ ] Results reproducible across runs
- [ ] Performance regression detection configured
- [ ] Reports suitable for publication

#### 1.3: Continuous Integration Setup (**COMPLETE**)
**Task:** Automate testing and validation

**Actual Implementation** (evolved from plan — see `specs/007-ci-setup/` for full spec):

The existing CI pipeline (`.github/workflows/ci.yml`, established in Phase 0.4) already
provided test, lint, and documentation jobs. Phase 1.3 gap analysis found that 8 of 10
spec requirements were already satisfied. The remaining work was adding benchmark
compilation verification.

**CI Jobs (sparse domain):**
- `test-sparse` — `cargo test --all-targets` on MSRV (1.87) + stable, `fail-fast: false`
- `lint-sparse` — `cargo fmt --check` + `cargo clippy --all-targets -- -D warnings`
- `doc-sparse` — `cargo doc --no-deps` with `RUSTDOCFLAGS: -D warnings`
- `bench-sparse` — `cargo bench --no-run` (compile benchmarks, don't execute)

All jobs use `actions/checkout@v4`, `dtolnay/rust-toolchain`, `Swatinem/rust-cache@v2`.
Feature-gated `test-util` code is tested via self-referencing dev-dependency in Cargo.toml.

**Deferred from original plan:**
- Multi-OS testing (macOS) — deferred until platform-specific code paths exist
- SPRAL comparison in CI — deferred to Phases 2-8 per Phase 0.3 decision
- Code coverage tracking — deferred (no solver code to measure coverage of)
- Python-based reporting scripts — replaced by native cargo/clippy/rustdoc output

**Success Criteria:**
- [x] Full test suite runs in <10 minutes (83 unit + 2 integration tests, ~20s)
- [x] Benchmark infrastructure validated (bench-sparse compiles criterion harness)
- [x] Clippy warnings treated as errors in CI
- [x] Clear failure reports with reproducible commands (native cargo output)

#### 1.4: Profiling and Debug Tools ✅
**Status:** Complete (branch `008-profiling-debug-tools`)

Implemented as `src/profiling/` and `src/debug/` modules behind `test-util` feature.

**What was built (vs original plan):**
- `ProfileSession` with RAII `SectionGuard` (replaces `ProfileRecorder` closure API — guard pattern is more ergonomic and handles panics)
- Chrome Trace export (implemented); flamegraph export (deferred — Chrome Trace viewers provide equivalent functionality)
- `MemoryTracker` with RSS snapshots (replaces allocation-tracking `MemoryTracker` — RSS is simpler and sufficient for phase-level memory analysis)
- `SparsityDisplay` text-based visualization (replaces `DebugVisualizer::visualize_sparsity` image output)
- `ETreeDisplay` text tree + statistics (replaces GraphViz DOT — text output is self-contained, no external dependency)
- Factorization animation deferred to later phases when factorization exists

**Success Criteria:**
- [x] Can profile any component in isolation (ProfileSession + SectionGuard)
- [x] Memory usage traceable to specific operations (MemoryTracker snapshots)
- [x] Visualization tools for debugging (SparsityDisplay + ETreeDisplay)
- [ ] Integration with `perf`, `valgrind`, etc. (deferred — orthogonal to library tooling)

### Phase 1 Exit Criteria

**Required Outcomes:**
1. ~~Complete test infrastructure can validate any solver component~~ — done (1.1)
2. ~~Benchmarking framework ready for use~~ — done (1.2)
3. ~~CI pipeline validates every commit~~ — done (1.3)
4. Debugging and profiling tools — in progress (1.4)

**Actual outcomes (Phases 1.1–1.3):**
- `SolverTest` trait with `MockSolver` for pre-solver validation
- `NumericalValidator` with configurable tolerances
- `TestCaseFilter` for composable test case selection
- Random matrix generators (PD and indefinite) behind `test-util` feature
- `Benchmarkable` trait + Criterion integration with baseline tracking
- CI: test (MSRV + stable), lint (fmt + clippy), doc, bench-compile

**Checkpoint:** Infrastructure validated with MockSolver and matrix loading benchmarks.

---

## Phase 2: APTP Data Structures (**COMPLETE**)

### Objectives
Define the data structures unique to indefinite APTP factorization: mixed 1×1/2×2 diagonal
storage, pivot tracking, and delayed-column bookkeeping. These structures are consumed by
the numeric factorization (Phase 5) and solve (Phase 7) phases.

### What was already completed elsewhere

The original plan included three subphases (2.1–2.3). Two were substantially completed
during earlier phases and are recorded here for traceability:

- **Original 2.1 (faer-sparse adoption + Matrix Market I/O)** — **Absorbed into Phase 0.4.**
  `io/mtx.rs` parses Matrix Market files into `SparseColMat<usize, f64>`.
  `io/registry.rs` provides a metadata-driven catalog backed by `metadata.json`.
  All 82 test matrices are loadable. Direct faer types are used throughout (see
  Permutation Decision below for rationale on not adding type aliases).

- **Original 2.3 (test infrastructure integration)** — **Absorbed into Phase 1.1.**
  `testing/cases.rs` defines `SolverTestCase` with a `SparseColMat<usize, f64>` field.
  `testing/harness.rs` provides the `SolverTest` trait, `TestResult`, `MetricResult`.
  `TestCaseFilter` supports composable loading from the registry.
  Test infrastructure already works end-to-end with faer types.

### Design Decisions

#### Permutation: Use faer's `Perm<usize>` directly (no custom wrapper)

APTP tracks three conceptual permutations:

| Permutation | Source | Phase | Nature |
|---|---|---|---|
| P_mc64 (matching) | MC64 algorithm | 4.2 | Standard index permutation |
| P_ord (fill-reducing) | AMD / METIS | 4.1 | Standard index permutation |
| P_piv (pivots) | APTP factorization | 5 | Sequence of 1×1/2×2 pivot decisions |

The first two (P_mc64, P_ord) are composed before factorization and applied as a
symmetric permutation `P^T A P`. faer's `Perm<usize>` supports all required operations:

- **Both arrays stored**: forward + inverse (zero-copy `inverse()` swaps references)
- **Composition**: `p1 * p2` via operator overloading
- **Symmetric permutation**: `perm.inverse() * A * perm.as_ref()` on `SparseColMat`
- **Vector permutation**: `perm * col_vector`
- **Checked construction**: `new_checked(forward, inverse, dim)` validates correctness

The third (P_piv) is **not a permutation array** — it is a sequence of pivot decisions
(1×1, 2×2 Bunch-Kaufman, or delayed). This is handled by `PivotType` and `MixedDiagonal`,
not by a permutation type.

**One gap:** faer requires both forward and inverse arrays at construction. External
ordering algorithms typically produce only the forward array. A small helper function
bridges this:

```rust
/// Construct a `Perm<usize>` from a forward permutation array,
/// computing the inverse automatically.
pub fn perm_from_forward(fwd: Vec<usize>) -> Perm<usize> {
    let n = fwd.len();
    let mut inv = vec![0usize; n];
    for (i, &j) in fwd.iter().enumerate() {
        inv[j] = i;
    }
    Perm::new_checked(fwd.into(), inv.into(), n)
}
```

**Consequence for Phase 4:** The `Permutation` struct shown in Phase 4's notional API
should be replaced with `faer::perm::Perm<usize>` and `perm_from_forward()`. The
`compose()` method becomes `p1 * p2`.

#### Type aliases: Not adding

The codebase already uses `SparseColMat<usize, f64>` and `Perm<usize>` directly.
Adding aliases like `AptpMatrix<T>` would be a layer of indirection without clear
benefit — faer types are already well-named and the codebase is small enough that
direct usage is readable. This can be revisited if the API surface grows significantly.

### Deliverables

**Task:** Define data structures unique to indefinite APTP factorization

These structures are the foundation for the numeric factorization kernel (Phase 5)
and the triangular solve (Phase 7).

**Key structures:**

- `PivotType` — enum tracking whether each pivot is 1×1, 2×2, or delayed
- `Block2x2<T>` — storage for a 2×2 symmetric block `[a, b; b, c]`
- `MixedDiagonal<T>` — the D factor in LDL^T, supporting mixed 1×1 and 2×2 blocks,
  with a solve method (`D x = b`)
- `DelayedColumn` / pivot delay tracking — bookkeeping for columns that fail the
  APTP stability check and must be passed to an ancestor node
- `perm_from_forward()` — helper for constructing `Perm<usize>` from ordering output

**Notional API** (to be refined during speccing):

```rust
/// Tracks whether each pivot is 1×1 or 2×2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PivotType {
    /// Standard 1×1 pivot.
    OneByOne,
    /// 2×2 Bunch-Kaufman pivot. The value identifies the paired column.
    TwoByTwo { partner: usize },
    /// Column was delayed (failed stability check).
    Delayed,
}

/// A 2×2 symmetric block [a, b; b, c].
#[derive(Debug, Clone)]
pub struct Block2x2<T> {
    pub first_col: usize,
    pub a: T,
    pub b: T,
    pub c: T,
}

/// Storage for the D factor in P^T A P = L D L^T where D has mixed
/// 1×1 and 2×2 diagonal blocks.
pub struct MixedDiagonal<T> {
    /// Per-column pivot classification.
    pivot_map: Vec<PivotType>,
    /// Diagonal values for 1×1 pivots (indexed by column).
    diag_1x1: Vec<T>,
    /// 2×2 blocks (one entry per pair; the lower-indexed column owns it).
    blocks_2x2: Vec<Block2x2<T>>,
    /// Number of columns (matrix dimension).
    n: usize,
}

impl<T: faer::RealField> MixedDiagonal<T> {
    pub fn new(n: usize) -> Self;
    pub fn set_1x1(&mut self, col: usize, value: T);
    pub fn set_2x2(&mut self, first_col: usize, block: Block2x2<T>);
    pub fn get_pivot_type(&self, col: usize) -> PivotType;

    /// Solve D x = b in-place where D has mixed 1×1/2×2 blocks.
    pub fn solve_in_place(&self, x: &mut [T]);

    /// Number of delayed pivots.
    pub fn num_delayed(&self) -> usize;
}
```

**Testing strategy:**
- Unit tests with hand-constructed mixed diagonals (known solves)
- Round-trip: construct D, solve D x = b, verify x
- Edge cases: all 1×1, all 2×2, single delayed column, empty matrix
- Property: `MixedDiagonal::solve_in_place` inverts `MixedDiagonal` multiplication

**Success Criteria:**
- [x] MixedDiagonal correctly stores mixed 1×1/2×2 blocks
- [x] `solve_in_place` produces correct results for mixed block patterns
- [x] PivotType tracks 1×1, 2×2, and delayed states
- [x] Block2x2 solve is numerically correct (2×2 symmetric system)
- [x] Efficient access patterns (no unnecessary allocation in solve)
- [x] `perm_from_forward()` helper tested (round-trip with faer `Perm` operations)

**Time Estimate:** 2–3 days

### Phase 2 Exit Criteria

**Required Outcomes:**
1. APTP-specific structures (`MixedDiagonal`, `PivotType`, `Block2x2`) implemented and tested
2. `perm_from_forward()` helper available for Phase 4 ordering integration
3. Design decisions documented (Permutation strategy, type alias decision — see above)

**Validation Questions:**
- Do APTP-specific structures handle all pivot patterns correctly?
- Does `MixedDiagonal::solve_in_place` match analytical solutions on hand-constructed cases?
- Are the structures efficient enough for the factorization inner loop?

**Checkpoint:** Verify MixedDiagonal solve correctness on hand-constructed examples with
mixed 1×1/2×2 patterns. Demonstrate PivotType tracking across a sequence of set operations.

---

## Phase 3: Symbolic Analysis (**COMPLETE**)

### Objectives
Build symbolic analysis using faer's elimination tree construction.

### Design Decisions

#### Transparent composition with faer (decided)

faer's `factorize_symbolic_cholesky` bundles ordering + elimination tree + sparsity
pattern computation into a single call, returning `SymbolicCholesky<usize>`. The
symbolic phase for indefinite LDL^T is **identical** to SPD Cholesky — pivoting is
purely a numeric-phase concern.

Rather than using faer's low-level primitives (`prefactorize_symbolic_cholesky`,
`factorize_simplicial_symbolic_cholesky`) and reconstructing the pipeline ourselves,
`AptpSymbolic` composes faer's high-level result with APTP-specific metadata:

- `AptpSymbolic` stores faer's `SymbolicCholesky<usize>` as an inner field
- Accessors delegate to faer (permutation, L structure, predicted nnz)
- APTP-specific fields are added alongside (pivot buffer sizing)
- Ordering is an **input parameter** using faer's `SymmetricOrdering` enum

This follows the project-wide **transparent composition** principle (see CLAUDE.md).

**Consequence for Phase 3/4 boundary:** Ordering is an input to symbolic analysis, not
a separate downstream phase. Phase 3 uses `SymmetricOrdering::Amd` as default. Phase 4
becomes "compute alternative orderings (METIS, MC64) that produce a `Perm<usize>`
passed as `SymmetricOrdering::Custom(perm)`."

#### SimplicialLFactor superseded by faer (decided)

The original plan proposed a custom `SimplicialLFactor` struct (CSC for L). faer already
provides `SymbolicSimplicialCholesky<I>` with the same data (col_ptr, row_idx) plus
methods (`len_val()`, `col_ptr()`, `row_idx()`, `factor()`). No need to duplicate this —
it's accessible through `SymbolicCholesky::raw()`.

#### Supernodal symbolic for multifrontal path (decided)

`SymbolicCholesky<usize>` computes both simplicial and supernodal symbolic structure.
The multifrontal factorization (Phase 6) needs access to the supernodal structure:
supernode column ranges, assembly tree, row structure per supernode, and postorder
traversal. `AptpSymbolic` should expose these through accessor methods that delegate
to faer's `SymbolicSupernodalCholesky<I>` (accessible via `SymbolicCholesky::raw()`).

This absorbs old Phase 8.1 (Supernode Detection) — faer handles supernode detection
internally. If APTP-specific supernode adjustments are needed (e.g., due to delayed
pivots crossing supernode boundaries), they can be handled at the numeric level.

### Deliverables

**Task:** Build `AptpSymbolic` as a thin composition over faer's symbolic pipeline,
adding APTP-specific delayed-pivot buffer estimation

This consolidates the original 3.1 (elimination tree), 3.2 (sparsity pattern), and
3.3 (API integration) into a single deliverable, since faer handles all three internally
and they share a single output type.

**Algorithm References:**
- Liu (1990) — "The role of elimination trees in sparse factorization"
- faer: `factorize_symbolic_cholesky`, `SymbolicCholesky<I>`, `SymbolicSimplicialCholesky<I>`

**Notional API** (to be refined during speccing):

```rust
use faer::sparse::linalg::cholesky::{
    factorize_symbolic_cholesky,
    SymbolicCholesky,
    CholeskySymbolicParams,
};
use faer::sparse::linalg::SymmetricOrdering;

/// APTP symbolic analysis result.
///
/// Composes faer's `SymbolicCholesky` (etree, L structure, permutation)
/// with APTP-specific metadata (delayed-pivot buffer estimates).
pub struct AptpSymbolic {
    /// faer's symbolic result (contains perm, L structure, etree).
    inner: SymbolicCholesky<usize>,

    /// APTP-specific: estimated extra space per column for delayed pivots.
    pivot_buffer: Vec<usize>,
}

impl AptpSymbolic {
    /// Run symbolic analysis with the given ordering strategy.
    ///
    /// Uses faer's `factorize_symbolic_cholesky` internally.
    /// Default ordering is AMD (`SymmetricOrdering::Amd`).
    pub fn analyze(
        matrix: SymbolicSparseColMatRef<'_, usize>,
        ordering: SymmetricOrdering<'_, usize>,
    ) -> Result<Self, SparseError> {
        let inner = factorize_symbolic_cholesky(
            matrix,
            Side::Upper,
            ordering,
            CholeskySymbolicParams::default(),
        )?;

        // Heuristic: 10% buffer per column for delayed pivots
        let pivot_buffer = /* computed from inner's col structure */;

        Ok(Self { inner, pivot_buffer })
    }

    // Delegate to faer
    pub fn perm(&self) -> Option<PermRef<'_, usize>> { self.inner.perm() }
    pub fn predicted_nnz(&self) -> usize { self.inner.len_val() }

    /// Statistics for debugging and diagnostics.
    pub fn statistics(&self) -> SymbolicStatistics { /* ... */ }
}

pub struct SymbolicStatistics {
    pub dimension: usize,
    pub predicted_nnz: usize,
    pub average_col_count: f64,
}
```

**Testing:**

```rust
#[test]
fn test_symbolic_analysis_all_matrices() {
    for case in all_test_cases() {
        let symbolic = AptpSymbolic::analyze(
            case.matrix().symbolic(),
            SymmetricOrdering::Amd,
        ).unwrap();

        let stats = symbolic.statistics();
        assert_eq!(stats.dimension, case.matrix().nrows());
        assert!(stats.predicted_nnz > 0);
        assert!(stats.average_col_count >= 1.0);
    }
}

#[test]
fn test_reproducible() {
    let matrix = load_test_matrix("arrow-10-pd").unwrap();

    let sym1 = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd).unwrap();
    let sym2 = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd).unwrap();

    assert_eq!(sym1.predicted_nnz(), sym2.predicted_nnz());
}

#[test]
fn test_custom_ordering_accepted() {
    // Verify that a user-supplied permutation is propagated through
    let matrix = load_test_matrix("arrow-5-pd").unwrap();
    let identity_perm = /* identity PermRef */;

    let symbolic = AptpSymbolic::analyze(
        matrix.symbolic(),
        SymmetricOrdering::Custom(identity_perm),
    ).unwrap();

    assert!(symbolic.perm().is_some());
}
```

**Success Criteria:**
- [x] `AptpSymbolic::analyze` works on all test matrices with AMD ordering
- [x] Custom ordering (`SymmetricOrdering::Custom`) accepted and propagated
- [x] Predicted NNZ matches faer's Cholesky prediction (since symbolic phase is identical)
- [x] Reproducible (deterministic)
- [x] Statistics reported correctly
- [x] Analysis time < 5% of factor time

**Time Estimate:** 3–5 days

### Phase 3 Exit Criteria

**Required Outcomes:**
1. `AptpSymbolic` wraps faer's `SymbolicCholesky` and adds APTP pivot buffer estimates
2. Symbolic analysis works on all test matrices (hand-constructed + SuiteSparse CI subset)
3. Ordering is accepted as input parameter (AMD default, custom supported)
4. Design decisions documented (transparent composition, SimplicialLFactor superseded)

**Validation Questions:**
- Does `AptpSymbolic` correctly delegate to faer for etree/structure/permutation?
- Are pivot buffer estimates reasonable (non-negative, proportional to column counts)?
- Does the API cleanly support Phase 4's alternative orderings via `SymmetricOrdering::Custom`?

**Checkpoint:** Run symbolic analysis on entire test suite with AMD ordering. Verify
predicted NNZ is reasonable. Demonstrate custom ordering passthrough.

**Time Estimate:** 3–5 days

#### Lessons Learned: AMD Ordering Quality on Full SuiteSparse Collection

Post-completion testing on the full 65-matrix SuiteSparse collection revealed that
AMD produces catastrophically poor orderings for several benchmark matrices. Comparing
our AMD-based predicted nnz(L) against paper-reported values (which used METIS):

| Matrix | Dim | AMD nnz(L) | METIS nnz(L) | Ratio |
|--------|-----|-----------|-------------|-------|
| nd3k | 9K | 22.8M | 12.9M | 1.8× |
| Si10H16 | 17K | 87.8M | 30.6M | 2.9× |
| nd6k | 18K | 72.8M | 39.8M | 1.8× |
| Si5H12 | 20K | 125.2M | 44.1M | 2.8× |
| sparsine | 50K | **1,037M** | ~50-100M (est.) | **~10-20×** |

Paper-reported nnz(L) values from Hogg, Ovtchinnikov, Scott (2016) Table III. All
papers use METIS (SPRAL's default ordering), not AMD.

**Key finding**: SPRAL uses METIS by default (`options%ordering = 1`), not AMD.
AMD is a simpler O(nnz) heuristic; METIS uses graph partitioning (nested dissection)
which produces dramatically better orderings for matrices with geometric structure
(FEM meshes, quantum chemistry). For sparsine, AMD produces 10-20× more fill than
METIS, making even symbolic analysis take 12+ seconds and the eventual numeric
factorization infeasible.

**Impact on plan**: METIS integration elevated from "consider for Phase 8" to
Phase 4.1 (before MC64). Without METIS, the solver cannot practically handle many
of the benchmark matrices from the APTP papers. This does not change faer
infrastructure reuse — METIS produces a `Perm<usize>` that plugs into
`SymmetricOrdering::Custom` unchanged. It is purely an input quality improvement.

---

## Phase 4: Ordering & Preprocessing (**COMPLETE**)

### Objectives
Implement fill-reducing ordering (METIS) and matching-based scaling (MC64) for
indefinite systems. Phase 3 testing revealed AMD produces catastrophically poor
orderings on many benchmark matrices — METIS is essential for practical performance.
MC64 provides complementary numerical preprocessing (diagonal dominance improvement)
for the APTP factorization kernel.

### What was already completed / absorbed

The original plan had four sub-deliverables (4.1–4.4). After the Phase 3 transparent
composition decision, most are absorbed:

- **4.1 (Ordering framework)**: Absorbed by Phase 3. faer's `SymmetricOrdering` enum
  already provides `Amd`, `Identity`, and `Custom(PermRef)`. No framework to build.
- **4.3 (Combined ordering)**: Thin — just `mc64_perm * amd_perm` via faer's `Perm`
  multiplication, then `SymmetricOrdering::Custom(combined.as_ref())`. Documented as
  a usage pattern, not a separate deliverable.
- **4.4 (Integration with analysis)**: Completely absorbed. `AptpSymbolic::analyze()`
  already accepts ordering as input and stores the permutation. No `AnalysisWithOrdering`
  struct needed.

### Design Decisions

#### Use `Perm<usize>` not custom `Permutation` (decided, Phase 2)

The original plan defined a custom `Permutation` struct with forward/inverse arrays,
`identity()`, `from_vec()`, `compose()`, etc. Per the Phase 2 decision, all permutation
handling uses faer's `Perm<usize>` directly. The custom struct is deleted. See Phase 2
design decisions for details.

#### METIS ordering: elevated to Phase 4.1

Originally deferred to Phase 8, METIS was elevated to Phase 4.1 after Phase 3
testing showed AMD produces 2-20× more fill than METIS on benchmark matrices.
SPRAL uses METIS by default — without it, the solver cannot practically handle
many of the APTP paper benchmark matrices.

**Integration path:** METIS produces a fill-reducing permutation via nested
dissection. The Rust `metis` crate provides FFI bindings. The output is a
`Perm<usize>` passed to `AptpSymbolic::analyze()` as
`SymmetricOrdering::Custom(perm.as_ref())`. No API changes needed — this is
purely an input quality improvement.

**Validation:** The existing full SuiteSparse test suite (`symbolic_analysis_full.rs`)
can directly validate METIS by comparing predicted nnz(L) against paper-reported
values. Structural property tests (etree, supernode structure, assembly tree) are
ordering-independent and apply unchanged.

#### Scaling flows to numeric phase, not symbolic (decided)

MC64 produces two outputs with different consumers:
- **Permutation** (`Perm<usize>`) → feeds into Phase 3 as `SymmetricOrdering::Custom`
- **Scaling** (`Vec<f64>`) → feeds into Phase 5 numeric factorization

For symmetric matrices, scaling is Â_ij = s_i · A_ij · s_j (same vector for rows and
columns). Scaling does not affect sparsity structure, so Phase 3 (symbolic) never sees it.
The numeric pipeline applies scaling before factorization and reverses it after solve:

1. Scale: Â = S A S
2. Factorize: Â = P^T L D L^T P
3. Solve: b̂ = S b, solve Â x̂ = b̂, then x = S x̂

MC64 returns `(Perm<usize>, Vec<f64>)` as a pure function. The exact API for threading
scaling into numeric factorization is deferred to Phase 5 design.

#### Structural singularity handling in MC64 (decided)

When the input matrix is structurally singular (`matched < n`), MC64 follows the
Duff & Pralet (2005, Section 4.2.3) algorithm — this is a warning, not an error:

1. **Initial matching**: MC64 returns a maximum cardinality matching of size `matched < n`.
2. **Submatrix extraction** (Property 4.2): The indices in the matching identify a
   structurally nonsingular submatrix. The restriction of A to matched-row × matched-row
   indices is guaranteed to be structurally nonsingular.
3. **Re-matching**: MC64 is re-run on the nonsingular submatrix to get a complete weighted
   matching and dual variables (u, v) for scaling.
4. **Scaling correction**: Matched indices get standard MC64SYM scaling
   (`s_i = √(exp(u_i) · exp(v_i))`). Unmatched indices get the Duff-Pralet correction:
   `s_i = 1 / max_k |a_ik · s_k|` over matched k, with convention 1/0 = 1.
5. **Permutation ordering**: Unmatched rows/columns are placed at the end of the
   permutation. SPRAL's `match_order.f90` does this explicitly before passing to METIS.

**Downstream consequences for Phases 5–8:**

The unmatched rows at the end of the permutation will appear as structurally zero
diagonal entries in the reordered matrix. The numeric factorization (Phase 5/6) will
encounter these as zero pivots. This is expected behavior:

- APTP's pivot delay mechanism naturally handles this — zero pivots are delayed, and
  the delayed pivot count in `AptpNumeric` will reflect the structural rank deficiency.
- The solve phase (Phase 7) must handle the rank-deficient case: either the system is
  consistent (solution exists but is not unique) or inconsistent (no solution). The
  `Inertia` from Phase 2 reports the number of zero eigenvalues, which equals `n - matched`.
- SPRAL reports this via `inform%flag = SSIDS_WARNING_ANAL_SINGULAR` — a warning that
  propagates to the caller but does not abort the solve.

No special handling is needed in Phase 4.2 beyond returning the `matched` count and
placing unmatched indices last. The downstream phases handle rank deficiency through
their existing pivot delay and inertia reporting mechanisms.

### Deliverables

#### 4.1: METIS Nested Dissection Ordering — COMPLETE

**Task:** Integrate METIS graph partitioning for fill-reducing ordering
**Branch:** `011-metis-ordering`

**Algorithm Reference:**
- Karypis & Kumar (1998) — "A Fast and High Quality Multilevel Scheme for Partitioning
  Irregular Graphs", SIAM J. Sci. Comput.
- George (1973) — Nested dissection theory

**Approach:**
METIS computes a fill-reducing ordering via multilevel graph partitioning (nested
dissection). For sparse symmetric matrices with geometric structure (FEM, quantum
chemistry, optimization), METIS typically produces 2-10× less fill than AMD.

The `metis` Rust crate provides safe FFI bindings to libmetis. The integration is
thin: extract the adjacency structure from `SparseColMat`, call METIS, wrap the
result as `Perm<usize>`.

**Notional API** (to be refined during speccing):

```rust
use faer::perm::Perm;
use faer::sparse::SparseColMat;

/// Compute a METIS nested dissection ordering for a symmetric sparse matrix.
///
/// Returns a fill-reducing permutation suitable for use with
/// `SymmetricOrdering::Custom(perm.as_ref())`.
pub fn metis_ordering(
    matrix: &SparseColMat<usize, f64>,
) -> Result<Perm<usize>, SparseError>;
```

**Testing:**

```rust
#[test]
fn test_metis_reduces_fill_vs_amd() {
    // On SuiteSparse matrices, METIS should produce less fill than AMD
    // for the majority of cases (>= 80%)
    for case in suitesparse_ci_subset() {
        let sym_amd = AptpSymbolic::analyze(case.symbolic(), SymmetricOrdering::Amd)?;
        let metis_perm = metis_ordering(&case)?;
        let sym_metis = AptpSymbolic::analyze(
            case.symbolic(),
            SymmetricOrdering::Custom(metis_perm.as_ref()),
        )?;
        // METIS should produce less or equal fill
    }
}

#[test]
fn test_metis_nnz_matches_paper_values() {
    // Compare predicted nnz(L) against values reported in Hogg et al. (2016)
    // Table III, which used METIS ordering
}
```

**Success Criteria:**
- [x] METIS ordering produces valid permutations on all test matrices
- [x] Predicted nnz(L) within tolerance of paper-reported values (Hogg et al. 2016 Table III)
- [x] METIS reduces fill vs AMD on >= 80% of SuiteSparse matrices (89% achieved)
- [x] Full SuiteSparse symbolic analysis completes in < 5 minutes total with METIS (146s achieved)
- [x] Integrates via `SymmetricOrdering::Custom` (no API changes)

**Lessons Learned:**
- METIS `perm` uses MATLAB convention `A' = A(perm, perm)`, meaning `perm[new] = old`
  (forward permutation). The inverse `iperm[old] = new`. Always verify FFI semantics.
- `metis-sys` vendors METIS 5.x C source (no system install). Compile adds ~30s.
- Only 2 of 35 Table III matrices are in the CI-subset. Wider tolerance (5×) needed
  because our symbolic prediction (Cholesky pattern) overestimates for indefinite matrices.

**Time Estimate:** 3-5 days

#### 4.2: MC64 Matching & Scaling — COMPLETE

**Task:** Implement weighted bipartite matching for indefinite matrix preprocessing

**Algorithm References:**
- Duff & Koster (2001) — "On algorithms for permuting large entries to the diagonal of a
  sparse matrix", SIAM J. Matrix Anal. Appl. — core MC64 weighted matching algorithm
  (shortest augmenting paths with Dijkstra on reduced costs).
  (`/workspace/rivrs-linalg/references/ssids/duff2001.md`)
- Duff & Pralet (2005) — "Strategies for scaling and pivoting for sparse symmetric
  indefinite problems" — symmetric adaptation (MC64SYM), scaling symmetrization,
  structural singularity handling (Property 4.2), Duff-Pralet scaling correction.
  (`/workspace/rivrs-linalg/references/ssids/duff2005.md`)

**Approach:**
MC64 computes a weighted bipartite matching that permutes large entries onto the diagonal,
improving diagonal dominance for indefinite systems. For APTP, this reduces the number of
delayed pivots. The algorithm also produces symmetric scaling factors. For symmetric
matrices, the MC64SYM approach (Duff & Pralet 2005) applies the unsymmetric MC64
algorithm and symmetrizes the scaling via geometric mean of row/column dual variables.
Structurally singular matrices are handled via partial matching with the Duff-Pralet
correction (see "Structural singularity handling" decision above).

faer has no matching or scaling functionality — this is genuinely new code.

**Notional API** (to be refined during speccing):

```rust
use faer::perm::Perm;
use faer::sparse::SparseColMat;

/// Result of MC64 matching-based preprocessing.
pub struct Mc64Result {
    /// Matching permutation (maximizes diagonal product).
    pub perm: Perm<usize>,
    /// Symmetric scaling factors: Â_ij = s_i * A_ij * s_j.
    pub scaling: Vec<f64>,
    /// Number of matched diagonal entries.
    pub matched: usize,
}

pub enum Mc64Job {
    /// Maximize product of diagonal entries (default for APTP).
    MaximumProduct,
    /// Maximize sum of diagonal entries.
    MaximumSum,
}

/// Compute MC64 matching and scaling for a symmetric sparse matrix.
pub fn mc64_matching(
    matrix: &SparseColMat<usize, f64>,
    job: Mc64Job,
) -> Result<Mc64Result, SparseError>;
```

**Usage with Phase 3 (ordering) and Phase 5 (scaling):**

```rust
// MC64 produces permutation + scaling
let mc64 = mc64_matching(&matrix, Mc64Job::MaximumProduct)?;

// Permutation feeds into symbolic analysis via SymmetricOrdering::Custom
let symbolic = AptpSymbolic::analyze(
    matrix.symbolic(),
    SymmetricOrdering::Custom(mc64.perm.as_ref()),
)?;

// Scaling feeds into numeric factorization (Phase 5 — API TBD)
// let numeric = factorize(&matrix, &symbolic, Some(&mc64.scaling))?;
```

**Combined ordering (MC64 + AMD):**

```rust
// When both matching and fill-reduction are desired:
// 1. Compute MC64 matching permutation
let mc64 = mc64_matching(&matrix, Mc64Job::MaximumProduct)?;
// 2. Run symbolic analysis with AMD on the matched pattern
//    (faer applies AMD internally when using SymmetricOrdering::Amd)
// 3. The effective ordering is the composition of both permutations
//    This is handled internally by factorize_symbolic_cholesky when
//    the matrix is pre-permuted, or can be composed explicitly:
let combined = mc64.perm * amd_perm;  // faer Perm composition
```

**Testing:**

```rust
#[test]
fn test_mc64_produces_valid_matching() {
    let matrix = create_badly_scaled_matrix();
    let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();

    // Should find full matching
    assert_eq!(result.matched, matrix.nrows());
    // Scaling factors should be positive
    assert!(result.scaling.iter().all(|&s| s > 0.0));
}

#[test]
fn test_mc64_improves_diagonal_dominance() {
    for case in test_cases_indefinite() {
        let result = mc64_matching(&case.matrix, Mc64Job::MaximumProduct).unwrap();

        // Apply permutation and scaling
        let scaled = apply_mc64(&case.matrix, &result);

        // Diagonal entries should be larger relative to off-diagonal
        let dominance_before = diagonal_dominance_metric(&case.matrix);
        let dominance_after = diagonal_dominance_metric(&scaled);
        assert!(dominance_after >= dominance_before);
    }
}

#[test]
fn test_mc64_perm_works_with_symbolic_analysis() {
    let matrix = load_test_matrix("arrow-10-indef").unwrap();
    let mc64 = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();

    // Should integrate cleanly with Phase 3 API
    let symbolic = AptpSymbolic::analyze(
        matrix.symbolic(),
        SymmetricOrdering::Custom(mc64.perm.as_ref()),
    ).unwrap();

    assert!(symbolic.predicted_nnz() > 0);
}
```

**Success Criteria:**
- [ ] MC64 produces valid matching (full matching on test suite)
- [ ] Improves diagonal dominance on indefinite matrices
- [ ] Scaling factors are positive and well-conditioned
- [ ] Permutation integrates with Phase 3 via `SymmetricOrdering::Custom`
- [ ] Returns `(Perm<usize>, Vec<f64>)` — no custom `Permutation` type

**Time Estimate:** 5–7 days

#### 4.3: Match-Order Condensation Pipeline — COMPLETE

**Task:** Combine MC64 matching with METIS ordering via cycle condensation
**Branch:** `013-match-order-condensation`

**Algorithm References:**
- Hogg & Scott (2013), HSL_MC80 — cycle condensation approach
- SPRAL `match_order.f90` (BSD-3) — reference implementation (`mo_split`, `mo_order`)

**Approach:**
SPRAL's `ordering=2` mode combines MC64 matching with METIS ordering, guaranteeing
that matched 2-cycle pairs are adjacent in the elimination order. This is critical
for APTP's 2x2 pivot detection efficiency. The pipeline:
1. MC64 matching → singletons + 2-cycles + longer cycles
2. Cycle splitting → decompose longer cycles into 2-cycles + singletons (SPRAL `mo_split`)
3. Condensed graph → fuse 2-cycle pairs into super-nodes (~n/2 dimension)
4. METIS ordering on condensed graph
5. Expand → map back to original indices, 2-cycle pairs get consecutive positions

**Implementation:**
- `match_order_metis()` in `src/aptp/ordering.rs` — public orchestrator
- `CycleDecomposition` — internal struct with partner/old_to_new/new_to_old mappings
- `split_matching_cycles()` — cycle decomposition using SPRAL's `mo_split` algorithm
- `build_condensed_adjacency()` — marker-array deduplication for condensed graph
- `expand_ordering()` — maps condensed METIS output back to full-size Perm<usize>
- `MatchOrderResult` — public result struct with ordering, scaling, and diagnostics

**Results:**
- 9/9 CI-subset matrices: valid pair adjacency (SC-001) and symbolic analysis (SC-007)
- Condensation ratios: 0.56-1.0 depending on 2-cycle count
- Fill quality: most matrices within 3% of unconstrained METIS; some matrices with heavy
  condensation (40%+ reduction) show 2-2.5x fill regression — expected trade-off for
  pair adjacency guarantee. SPRAL uses this mode by default for indefinite systems.
- Existing MC64 and METIS tests unaffected (SC-006)

**Lessons Learned:**
- Matrices with many 2-cycles (e.g., cvxqp3: 7475 pairs in 17500 nodes) show significant
  fill regression because condensation fundamentally changes the graph's separator structure.
  This is the inherent trade-off of the match-order approach. SPRAL accepts this because
  pair adjacency improves factorization quality (fewer delayed pivots) even if symbolic
  fill is higher.
- The `is_matched` field on `Mc64Result` is essential — `build_singular_permutation()`
  assigns unmatched rows to free columns, making `fwd[i]==i` unreliable for distinguishing
  singletons from unmatched indices.

### Phase 4 Exit Criteria

**Required Outcomes:**
1. METIS ordering produces fill predictions within 20% of paper-reported values
2. Full SuiteSparse symbolic analysis completes in reasonable time with METIS
3. MC64 matching produces valid matchings on all indefinite test matrices
4. Scaling improves diagonal dominance measurably
5. Both orderings integrate via `SymmetricOrdering::Custom`
6. Combined METIS + MC64 ordering demonstrated
7. Scaling vector stored and documented for Phase 5 consumption

**Validation Questions:**
- Does METIS ordering match the fill predictions reported in the APTP papers?
- Does MC64 find full matching on hard indefinite problems?
- Does scaling improve diagonal dominance (reduce off-diagonal/diagonal ratio)?
- Does the MC64 permutation reduce delayed pivots compared to METIS alone? (May need Phase 5
  to fully validate — can check indirectly via diagonal dominance metrics.)

**Checkpoint:** Run METIS ordering on full SuiteSparse collection, verify fill predictions
match paper values. Run MC64 on indefinite test matrices, verify full matching, positive
scaling factors, and improved diagonal dominance. Demonstrate integration with Phase 3
symbolic analysis.

---

## Phase 5: Dense APTP Factorization Kernel (**COMPLETE**)

### Objectives
Implement A Posteriori Threshold Pivoting for dense symmetric indefinite matrices. This is the core numerical kernel that Phase 6 (sparse factorization) will call on each frontal matrix.

### What was already completed / absorbed

- **5.4 (Integration with faer)**: Per transparent composition, we use faer's dense
  matmul and triangular solve from the start — not as a separate integration step.
  There is no "naive" version followed by a "faer" version.
- **5.5 (Validation on hard indefinite)**: Validation testing is exit criteria for
  this feature, using existing test infrastructure (NumericalValidator, TestCaseFilter,
  SuiteSparse collection). Not a separate deliverable.
- **5.3 (Two-level APTP)**: Deferred to Phase 8 (performance optimization). Two-level
  blocking is a cache performance optimization for large frontal matrices. In the
  simplicial solver (Phases 5–7), fronts are small (column-level). Two-level APTP
  becomes relevant when supernodal fronts reach hundreds of rows.

### Design Decisions

#### Use Phase 2 types for D storage (decided)

The plan originally stored D as `Mat<f64>` (dense matrix). Phase 2 already defines
`MixedDiagonal<T>` specifically for mixed 1×1/2×2 block diagonal storage with a
`solve_in_place` method. The APTP factorization result should use `MixedDiagonal`,
not `Mat<f64>`.

Phase 5 also defines a `Pivot` enum with data (`OneByOne { index, value }`,
`TwoByTwo { indices, block }`). This is complementary to Phase 2's `PivotType`
(classification tag): `PivotType` goes in `MixedDiagonal`'s pivot map, while the
APTP kernel returns a richer pivot log for diagnostics.

#### Use faer dense BLAS from day one (decided)

Per transparent composition, the Schur complement update uses `faer::linalg::matmul`
and triangular solves use faer's TRSM. The parallelism parameter is `Par` (not the
old `Parallelism` enum). There is no separate "faer integration" deliverable.

#### FallbackStrategy::Delayed semantics (decided)

In the dense kernel, "delay" means "mark this column as uneliminated and return it
to the caller." The kernel does not know about the elimination tree or parent nodes.
The caller (Phase 6 sparse factorization) decides what to do with delayed columns
(pass to parent node in the multifrontal assembly).

#### faer's Bunch-Kaufman as reference for fallback path (noted)

faer's `linalg::cholesky::bunch_kaufman` module provides a full dense indefinite
factorization with multiple pivoting strategies (Partial, PartialDiag, Rook,
RookDiag, Full). This is a valuable reference for the 2×2 pivot selection logic
in APTP's fallback path. The APTP kernel may delegate to faer's BK implementation
for the fallback, or implement a simplified version — to be decided during speccing.

### Deliverables

**Task:** Implement the APTP algorithm for dense symmetric indefinite matrices,
including fallback pivoting strategies

**Algorithm Reference:**
- Hogg, Duff, Lopez (2020) - Section 3 "A posteriori threshold pivoting" (Algorithm 1)
- Hogg, Duff, Lopez (2020) - Section 4.1 "Fallback strategies"
- faer `linalg::cholesky::bunch_kaufman` module — reference for 2×2 pivot selection

**Approach:**

The APTP algorithm is fundamentally different from both faer's dynamic regularization
(which corrects pivots without recomputing) and standard Bunch-Kaufman (which decides
pivot size before elimination). APTP:

1. Factors a block optimistically, assuming all 1×1 pivots
2. Checks stability a posteriori: `|l_ij| < 1/threshold` for all entries
3. If any entry fails: falls back to 2×2 Bunch-Kaufman pivot, or marks column as delayed
4. Delayed columns are returned to the caller (the sparse factorization in Phase 6)

This is genuinely new work — faer has no APTP implementation.

**Notional API** (to be refined during speccing):

```rust
use faer::{Mat, MatRef, MatMut, Par};

/// Result of dense APTP factorization on a frontal matrix.
pub struct AptpFactorization {
    /// Unit lower triangular factor (dense).
    l: Mat<f64>,
    /// Block diagonal factor — mixed 1×1 and 2×2 blocks (from Phase 2).
    d: MixedDiagonal<f64>,
    /// Columns that failed the stability check and could not be eliminated.
    /// The caller (sparse factorization) decides what to do with these.
    delayed_cols: Vec<usize>,
    /// Per-step pivot record for diagnostics.
    pivot_log: Vec<AptpPivotRecord>,
    /// Summary statistics.
    stats: AptpStatistics,
}

/// Diagnostic record for a single pivot step.
pub struct AptpPivotRecord {
    pub col: usize,
    pub pivot_type: PivotType,  // From Phase 2
    pub max_l_entry: f64,       // Worst stability metric
    pub was_fallback: bool,     // True if 2×2 fallback was needed
}

pub struct AptpStatistics {
    pub num_1x1: usize,
    pub num_2x2: usize,
    pub num_delayed: usize,
    pub max_l_entry: f64,
}

/// Configuration for the APTP kernel.
pub struct AptpOptions {
    /// Stability threshold: entries must satisfy |l_ij| < 1/threshold.
    /// Typical value: 0.01 (allowing growth factor up to 100).
    pub threshold: f64,
    /// What to do when a 1×1 pivot fails the stability check.
    pub fallback: AptpFallback,
}

pub enum AptpFallback {
    /// Try 2×2 Bunch-Kaufman pivot; if that also fails, delay the column.
    BunchKaufman,
    /// Use complete pivoting on the failed block.
    CompletePivoting,
    /// Immediately delay the column (return to caller).
    Delay,
}

/// Factor a dense symmetric matrix using APTP.
///
/// Uses faer's matmul for Schur complement updates and TRSM for
/// triangular solves internally.
pub fn aptp_factor(
    a: MatRef<f64>,
    options: &AptpOptions,
) -> Result<AptpFactorization, SparseError>;
```

**Usage example:**

```rust
let a = load_dense_frontal_matrix();
let options = AptpOptions {
    threshold: 0.01,
    fallback: AptpFallback::BunchKaufman,
};

let result = aptp_factor(a.as_ref(), &options)?;

// Check stability guarantee
assert!(result.stats.max_l_entry < 1.0 / options.threshold);

// Delayed columns are returned to the caller
if !result.delayed_cols.is_empty() {
    // In sparse context: pass these to the parent node
    println!("{} columns delayed", result.delayed_cols.len());
}

// Solve D x = b using Phase 2's MixedDiagonal
let mut x = b.clone();
result.d.solve_in_place(&mut x);
```

**Testing strategy:**

- **Positive definite**: should factor with no delays, all 1×1 pivots
- **Indefinite**: should factor correctly, some 2×2 pivots expected
- **Stability bounds**: `|l_ij| < 1/threshold` verified on 100+ random matrices
- **Reconstruction**: `||A - LDL^T|| / ||A|| < 10^-12` (using existing NumericalValidator)
- **Fallback comparison**: all AptpFallback variants produce valid factorizations
- **Hard indefinite**: matrices from SuiteSparse collection (via TestCaseFilter)
- **Pathological cases**: near-singular, highly indefinite, clustered eigenvalues

**Success Criteria:**
- [ ] Factors positive definite matrices with no delays
- [ ] Handles indefinite matrices correctly (2×2 pivots used as needed)
- [ ] Stability bound `|l_ij| < 1/threshold` enforced for all entries
- [ ] Reconstruction error < 10^-12 on all test matrices
- [ ] All fallback strategies produce valid factorizations
- [ ] Delayed columns correctly identified and returned
- [ ] Uses faer matmul/TRSM for Schur complement (no naive implementation)
- [ ] Successfully factors hard indefinite SuiteSparse matrices

**Time Estimate:** 2–3 weeks

### Phase 5 Exit Criteria

**Required Outcomes:**
1. Dense APTP kernel implemented with fallback strategies
2. Uses Phase 2's `MixedDiagonal` and `PivotType` for D storage
3. Uses faer dense BLAS (matmul, TRSM) for all heavy computation
4. Hard indefinite problems factor correctly with acceptable residuals
5. Stability bounds enforced and verified

**Validation Questions:**
- Does APTP maintain stability (`|l_ij| < 1/threshold`) on all test cases?
- Are delayed columns correctly reported to the caller?
- Does reconstruction error meet tolerance on hard indefinite matrices?

**Checkpoint: HIGH RISK VALIDATION**

This is the critical checkpoint. APTP is the foundation of the entire solver.
Run comprehensive validation using the existing test infrastructure:

- All hand-constructed matrices (15) via `TestCaseFilter::hand_constructed()`
- All hard indefinite SuiteSparse matrices via `TestCaseFilter` difficulty filter
- 100+ random symmetric indefinite matrices via generators (Phase 0.5)
- Pathological cases: near-singular, highly indefinite, clustered eigenvalues

**If this checkpoint fails:** Investigate before proceeding. Options:
1. Tune threshold parameter
2. Improve fallback strategy (try CompletePivoting vs BunchKaufman)
3. Consult SPRAL's dense kernel implementation (BSD-3) for reference
4. Check for numerical issues in Schur complement update

**If it passes:** Working core factorization kernel. The remaining phases
are assembly (Phase 6) and solve + end-to-end integration (Phase 7).

#### Phase 5 Post-Completion Review (added after Phase 5 SPRAL comparison)

The following Phase 5 design decisions were made for correctness and simplicity.
Later phases should revisit them if benchmarks indicate a need:

1. **Per-column backup/restore** (vs SPRAL's per-block backup): Phase 5 uses
   simple Vec backup for each column attempt, which is the degenerate case of
   SPRAL's per-block backup with block size = 1. When Phase 8.1 adds two-level
   blocking, the backup strategy must change to per-block — see Phase 8.1 notes.

2. **MixedDiagonal with PivotType enum** (vs SPRAL's flat d[2n] with Inf sentinel):
   Our approach trades a small amount of memory (the pivot_map vec) for type safety.
   SPRAL's Inf sentinel allows O(1) pivot-type detection during solve without a
   separate type array. If Phase 8.2 benchmarks show solve is a bottleneck (>25%
   of total time), revisit D storage format — see Phase 9.1 notes.

3. **Column-by-column Schur complement** (BLAS-2, not BLAS-3): Single-level APTP
   uses rank-1/rank-2 updates. BLAS-3 blocking gives ~2-5x speedup on fronts >128
   rows. This is the primary motivation for Phase 8.1.

4. **No TPP (Traditional Threshold Pivoting) fallback**: SPRAL can re-factorize
   heavily-delayed columns using TPP. We return delayed columns to the parent node
   instead. If Phase 8.2 benchmarks show excessive delays propagating up the tree
   on real matrices, consider adding TPP — see Phase 9.1 notes.

---

## Phase 6: Multifrontal Numeric Factorization (**COMPLETE**)

### Objectives
Implement the multifrontal factorization loop: assemble frontal matrices from
original matrix entries and child contributions, factor each front using Phase 5's
dense APTP kernel, and propagate contribution blocks up the assembly tree.

### Design Decisions

#### Multifrontal is the primary factorization path (decided)

APTP is inherently a blocked algorithm — it factors a block of columns optimistically,
then checks stability a posteriori. A column-by-column "simplicial APTP" would be a
different algorithm from Hogg et al. (2020). The multifrontal approach directly matches
the paper's algorithm: each supernodal front is a dense block that Phase 5's
`aptp_factor()` operates on.

faer provides `factorize_simplicial_numeric_ldlt` with dynamic regularization for
users who need a simpler indefinite factorization. Our value-add is the APTP algorithm,
which requires blocked/frontal computation.

#### Supernode detection via faer (decided)

faer's `SymbolicCholesky<usize>` (wrapped by Phase 3's `AptpSymbolic`) already computes
supernodal structure via `SymbolicSupernodalCholesky<I>`. This includes supernode column
ranges, assembly tree, row structure per supernode, and postorder traversal. We delegate
to faer rather than reimplementing supernode detection.

Old Phase 8.1 (Supernode Detection) is absorbed here.

#### Frontal matrices use Phase 2 types (decided)

The D factor within each front uses `MixedDiagonal<f64>` (Phase 2), not `Mat<f64>`
or `Vec<f64>`. Pivot decisions are tracked as `PivotType` (Phase 2).

#### Delayed pivot propagation (noted)

When `aptp_factor()` returns delayed columns from a front, they are added to the
parent front's row structure as additional fully-summed columns. This is the standard
multifrontal delayed pivot mechanism (Duff & Reid 1983). The exact propagation
mechanism will be refined during speccing.

### Deliverables

**Task:** Implement the multifrontal factorization loop that processes each
supernode's frontal matrix in assembly-tree postorder

**Algorithm Reference:**
- Duff & Reid (1983) - "The multifrontal solution of indefinite sparse symmetric
  linear equations"
- Liu (1992) - "The Multifrontal Method for Sparse Matrix Solution: Theory and
  Practice"
- Hogg, Duff, Lopez (2020) - Sections 2-3 (APTP within multifrontal framework)

**Approach:**

For each supernode in assembly-tree postorder:
1. **Create** frontal matrix (F11/F21/F22 partitioning from supernode structure)
2. **Scatter** original matrix entries into frontal matrix
3. **Extend-add** contribution blocks from child supernodes
4. **Factor** F11 using Phase 5's `aptp_factor()` → L11, D11, delayed columns
5. **Solve** for L21: L21 = F21 × L11^{-T} × D11^{-1} (using faer TRSM)
6. **Compute** Schur complement: contribution = F22 - L21 × D11 × L21^T (using faer matmul)
7. **Store** L11, D11, L21 as this supernode's factors
8. **Propagate** contribution block and delayed columns to parent

**Notional API** (to be refined during speccing):

```rust
/// Frontal matrix for a supernode.
///
/// Structured as:
///   [F11  |  ---]
///   [F21  |  F22]
///
/// F11 (fully_summed × fully_summed): factored with APTP.
/// F21 (remaining × fully_summed): solved via triangular solve.
/// F22 (remaining × remaining): Schur complement → contribution to parent.
pub struct FrontalMatrix {
    /// Dense storage.
    data: Mat<f64>,
    /// Number of fully summed columns (supernode size + any delayed from children).
    fully_summed: usize,
    /// Global row indices (maps local row → global column index).
    row_indices: Vec<usize>,
}

impl FrontalMatrix {
    pub fn f11(&self) -> MatRef<f64>;
    pub fn f21(&self) -> MatRef<f64>;
    pub fn f22(&self) -> MatRef<f64>;
    pub fn contribution(&self) -> MatRef<f64> { self.f22() }
}

/// Per-supernode factorization result.
pub struct FrontFactors {
    /// Dense L11 (unit lower triangular).
    l11: Mat<f64>,
    /// D11 (mixed 1×1/2×2 from Phase 2).
    d11: MixedDiagonal<f64>,
    /// Dense L21 (subdiagonal block).
    l21: Mat<f64>,
    /// Global indices of delayed columns (passed to parent).
    delayed_cols: Vec<usize>,
}

/// Complete numeric factorization result.
pub struct AptpNumeric {
    /// Per-supernode factors, indexed by supernode ID.
    front_factors: Vec<FrontFactors>,
    /// Factorization statistics.
    stats: FactorizationStats,
}

pub struct FactorizationStats {
    pub total_1x1_pivots: usize,
    pub total_2x2_pivots: usize,
    pub total_delayed: usize,
    pub max_front_size: usize,
}

impl AptpNumeric {
    /// Factor a sparse symmetric matrix using multifrontal APTP.
    ///
    /// Requires symbolic analysis from Phase 3.
    pub fn factor(
        symbolic: &AptpSymbolic,
        matrix: &SparseColMat<usize, f64>,
        options: &AptpOptions,
    ) -> Result<Self, SparseError>;
}
```

**Usage example:**

```rust
// Symbolic analysis (Phase 3)
let symbolic = AptpSymbolic::analyze(
    matrix.symbolic(),
    SymmetricOrdering::Amd,
)?;

// Numeric factorization (Phase 6)
let numeric = AptpNumeric::factor(
    &symbolic,
    &matrix,
    &AptpOptions { threshold: 0.01, fallback: AptpFallback::BunchKaufman },
)?;

println!("Stats: {} delays, max front size {}",
    numeric.stats.total_delayed,
    numeric.stats.max_front_size);
```

**Testing strategy:**

- **Small matrices**: compare factorization against dense APTP (convert to dense,
  factor, verify same reconstruction error)
- **Assembly correctness**: verify assembled fronts match expected dense blocks
  for hand-constructed elimination trees
- **Contribution propagation**: verify Schur complement is correctly passed to parent
- **Delayed pivots**: construct matrices where delays occur, verify they propagate
  correctly and the factorization still produces correct results
- **All test matrices**: reconstruction error < 10^-12 using NumericalValidator

**Success Criteria:**
- [x] Multifrontal factorization produces correct factors on all hand-constructed matrices
- [x] Reconstruction error < 10^-12 on all test matrices (where dense reconstruction feasible)
- [x] Assembly correctly scatters original entries + child contributions
- [x] Delayed pivots propagate to parent and are eventually eliminated
- [x] Uses Phase 5's `aptp_factor()` for intra-front factorization
- [x] Uses Phase 2's `MixedDiagonal` for D storage
- [x] Matches dense APTP factorization on small problems (converted to dense)
- [x] Factors all CI-subset SuiteSparse matrices correctly

### Phase 6 Exit Criteria (**MET**)

**Required Outcomes:**
1. ~~Multifrontal numeric factorization working end-to-end~~ ✓
2. ~~All hand-constructed and easy indefinite matrices factor correctly~~ ✓
3. ~~Reconstruction error meets tolerance on all test cases~~ ✓
4. ~~Delayed pivot mechanism working~~ ✓

**Validation Answers:**
- Factorization matches dense APTP on small test matrices (test_dense_equivalence, test_single_supernode_matches_dense)
- Contribution blocks computed correctly (test_multi_level_contribution_flow, test_extract_contribution_structure)
- Delayed pivots eventually get eliminated (test_delayed_pivot_propagation)

**Note:** Full SuiteSparse backward error validation deferred to Phase 7
(requires triangular solve). Dense reconstruction is O(n²) and impractical
for large matrices. SPRAL uses solve-based backward error exclusively.

---

## Phase 7: Triangular Solve & Solver API (**COMPLETE**)

### Objectives
Implement forward/backward substitution through the multifrontal factor structure
and assemble the complete solver pipeline into a user-facing API. This is the
"working solver" milestone: analyze → factor → solve, validated end-to-end against
the full SuiteSparse test suite.

### Key Findings: Accuracy & Ordering

Post-implementation benchmarking revealed that **ordering strategy is the dominant
factor in accuracy** for our single-level APTP kernel:

| Ordering | bratu3d delays | backward error | factor time |
|----------|---------------:|---------------:|------------:|
| METIS | 53,841 | 5.59e-3 | 7.0s |
| MatchOrderMetis | 1 | 1.00e-9 | 1.6s |

**Root cause**: Matrices like bratu3d (3D Bratu problem, saddle-point structure)
cause massive pivot delays with plain METIS ordering. MC64 matching+scaling
(MatchOrderMetis) pre-pairs difficult pivots, eliminating virtually all delays.

**Default ordering changed to `MatchOrderMetis`** — this matches SPRAL's
recommendation for indefinite problems (`ordering=2` mode). All new tests and
benchmarks should use this default.

**Random matrix validation**: 30+ random symmetric indefinite matrices (n=50..1000)
achieve backward error ~1e-17 with both METIS and MatchOrderMetis, confirming the
dense APTP kernel is correct. Accuracy degradation is ordering- and front-size
dependent.

**CI SuiteSparse results** (8 matrices, MatchOrderMetis, excluding sparsine):

| Matrix | n | backward error | status (5e-11) |
|--------|--:|---------------:|:--------------:|
| t2dal | 4,257 | 2.71e-18 | PASS |
| bloweybq | 10,001 | 2.99e-11 | PASS |
| ncvxqp1 | 12,111 | 4.69e-14 | PASS |
| cfd2 | 123,440 | 1.03e-18 | PASS |
| bratu3d | 27,792 | 1.00e-9 | WARN |
| cvxqp3 | 17,500 | 1.51e-7 | WARN |
| stokes128 | 49,666 | 6.08e-6 | FAIL |
| ncvxqp3 | 75,000 | 1.39e-7 | WARN |

Matrices with zero delayed pivots achieve machine precision. The remaining gap
to SPRAL's 5e-11 threshold correlates with delay count and is expected to improve
with two-level APTP (Phase 8.1), which reduces accumulated rounding error per front.

**Testing guidance**: Use `MatchOrderMetis` for all production and integration
testing. `OrderingStrategy::Metis` remains available for comparison but should not
be used as default in tests or examples.

### Rationale for merging old Phases 7 & 8

The triangular solve cannot be meaningfully validated without end-to-end tests
(backward error requires the full factor → solve pipeline), and the `SparseLDLT`
wrapper is a thin facade over existing types (`AptpSymbolic`, `AptpNumeric`). The
API design decisions (MemStack workspace, solver options) are straightforward
applications of faer's existing conventions. There is no natural stopping point
where "solve works but the API doesn't exist yet" is a useful deliverable.

### What was absorbed from other phases

- **Old Phase 7.2** (User-Facing API): originally split into a separate phase
- **Old Phase 7.3** (Integration Tests): originally split into a separate phase

### Design Decisions

#### SparseLDLT follows faer's solver wrapper pattern (decided)

faer's `Llt<I, T>` in `sparse/solvers.rs` wraps `SymbolicLlt<I>` + `numeric: Vec<T>`
and provides solve methods. Our `SparseLDLT` is the APTP equivalent, wrapping
`AptpSymbolic` + `Option<AptpNumeric>` + optional scaling factors. The `Option`
allows the analyze-then-factor-later pattern. Users familiar with faer's solver API
will find ours immediately recognizable.

#### Stack-based workspace allocation for solve (decided)

faer uses `MemStack` for all temporary workspace during solve — no heap allocations
in the hot path. All scratch space is queried upfront via `solve_in_place_scratch()`
and allocated from a reusable `MemBuffer`. Our **solve** methods should accept
`&mut MemStack` and expose `*_scratch() -> StackReq` methods following faer's
convention. The one-shot `solve_full` convenience method can allocate internally.

**Scoping note**: MemStack applies to the **solve** hot path in this phase.
Factorization currently allocates `FrontalMatrix` (via `Mat::zeros`) per supernode
in Phase 6's assembly loop — this is inherent to the multifrontal algorithm and
making it MemStack-based would require arena allocation, which is deferred to
Phase 9.1. The `factor()` method on `SparseLDLT` does **not** take `&mut MemStack`
in this phase.

#### Scaling integration at the SparseLDLT level (decided)

Phase 4's `match_order_metis()` returns `MatchOrderResult { ordering, scaling, ... }`.
The scaling factors need to persist between analyze and solve. Neither `AptpSymbolic`
(which wraps faer's `SymbolicCholesky`) nor `AptpNumeric` knows about MC64 scaling.

Scaling is applied/unapplied at the `SparseLDLT` level:
- `SparseLDLT` stores `Option<Vec<f64>>` scaling factors alongside the symbolic
  analysis (populated when `AnalyzeOptions` requests MC64 ordering)
- `factor()` applies scaling to matrix entries during assembly (SPRAL's approach:
  scale entries as they're scattered into frontal matrices, via
  `scaled_value = s[i] * a[i][j] * s[j]`)
- `solve()` applies `S * rhs` before forward solve and `S * solution` after
  backward solve
- Phase 6's `AptpNumeric::factor()` remains scaling-unaware — `SparseLDLT::factor()`
  pre-scales and delegates

#### Solve is a free function, not an AptpNumeric method (decided)

Phase 6's `AptpNumeric` is a data-only type: `{ front_factors, stats, n }`. The
solve requires both `AptpSymbolic` (supernode structure, permutation, assembly tree)
and `AptpNumeric` (per-supernode L, D factors). Rather than putting the solve on
`AptpNumeric` (which would need `&AptpSymbolic` passed in) or duplicating it on
`SparseLDLT`, the core solve is a free function:

```rust
pub(crate) fn aptp_solve(
    symbolic: &AptpSymbolic,
    numeric: &AptpNumeric,
    rhs: ColMut<f64>,
    stack: &mut MemStack,
) -> Result<(), SparseError>;
```

`SparseLDLT::solve_in_place()` adds scaling/unscaling and permutation around this.

#### Refactor is trivially factor (noted)

Given Phase 6's API, `refactor()` with the same sparsity pattern is identical to
calling `factor()` again — the symbolic analysis is reusable by construction, and
there is no internal state to reset. `refactor()` exists for API clarity and to
mirror SPRAL's interface, but delegates directly to `factor()`.

### Deliverables

**Task 1:** Implement the per-supernode triangular solve

**⚠️ Planning note — per-supernode index mapping is the hard part:**

The high-level solve algorithm (permute → forward → D solve → backward → unpermute)
is straightforward, but the per-supernode mechanics involve significant index-mapping
complexity that must be carefully designed. Each `FrontFactors` stores:
- `l11` (ne × ne unit lower triangular)
- `d11` (`MixedDiagonal` with `num_delayed == 0`)
- `l21` (r × ne subdiagonal block)
- `local_perm` (APTP pivot permutation within the front)
- `col_indices` (global column positions of eliminated columns)
- `row_indices` (global row positions for L21 rows)

The forward solve for supernode s must:
1. Gather entries from the global RHS vector using `col_indices` into a local vector
2. Apply `local_perm` to reorder the local vector to match APTP's pivot order
3. Solve `L11 * y_local = rhs_local` via dense TRSV (faer)
4. Scatter updates to non-fully-summed rows: for each row i in `row_indices`,
   subtract `L21[i, :] * y_local` from `rhs[row_indices[i]]`

The D solve applies `d11.solve_in_place()` to each supernode's eliminated entries
(per-supernode, not a single global D solve).

The backward solve reverses: gather, apply L21^T updates, solve L11^T, scatter.

**This index-mapping logic is the most error-prone part of the solve and must be
designed with extreme care.** During speccing:
- Study SPRAL's `ldlt_app_solve_fwd` and `ldlt_app_solve_bwd` in
  `spral/src/ssids/cpu/kernels/ldlt_app.cxx` (BSD-3) for reference
- Study SPRAL's `solve()` in `spral/src/ssids/cpu/subtree.hxx` for the
  assembly-tree traversal and per-supernode dispatch
- Study faer's `factorize_supernodal_numeric_intranode_lblt` solve path for
  how faer handles per-supernode gather/scatter in its own LBLT solver
- Write exhaustive unit tests for the index-mapping functions **before**
  integrating them into the full solve
- Test with hand-constructed matrices where the expected per-supernode local
  vectors can be computed analytically

**High-level approach:**

Given P^T A P = L D L^T from the multifrontal factorization, solve Ax = b via:
1. Permute: b̂ = P b (using `AptpSymbolic::perm()`)
2. Scale (if MC64): b̂ = S b̂ (using stored scaling factors)
3. Forward solve: L y = b̂ (postorder traversal, per-supernode gather/TRSV/scatter)
4. Diagonal solve: D z = y (per-supernode `d11.solve_in_place()`)
5. Backward solve: L^T w = z (reverse postorder, per-supernode gather/TRSV/scatter)
6. Unscale: w = S w
7. Unpermute: x = P^T w

**Task 2:** Create the user-facing `SparseLDLT` API and validate end-to-end

**Notional API** (to be refined during speccing):

```rust
/// High-level sparse symmetric indefinite solver.
///
/// Implements the three-phase API: analyze → factor → solve.
/// The symbolic analysis is reusable across factorizations with
/// the same sparsity pattern.
pub struct SparseLDLT {
    symbolic: AptpSymbolic,
    numeric: Option<AptpNumeric>,
    scaling: Option<Vec<f64>>,  // MC64 scaling factors (if ordering used MC64)
}

impl SparseLDLT {
    /// Symbolic analysis phase.
    pub fn analyze(
        matrix: SymbolicSparseColMatRef<'_, usize>,
        options: &AnalyzeOptions,
    ) -> Result<Self, SparseError>;

    /// Numeric factorization phase.
    /// Heap-allocates per-supernode frontal matrices internally.
    /// MemStack-based factorization workspace is deferred to Phase 9.
    pub fn factor(
        &mut self,
        matrix: &SparseColMat<usize, f64>,
        options: &FactorOptions,
    ) -> Result<(), SparseError>;

    /// Solve phase (allocating).
    pub fn solve(&self, rhs: ColRef<f64>, stack: &mut MemStack) -> Result<Col<f64>, SparseError>;

    /// Solve phase (in-place).
    pub fn solve_in_place(&self, rhs: ColMut<f64>, stack: &mut MemStack) -> Result<(), SparseError>;

    /// Workspace requirement for solve.
    pub fn solve_scratch(&self, rhs_ncols: usize) -> StackReq;

    /// One-shot: analyze + factor + solve (allocates workspace internally).
    pub fn solve_full(
        matrix: &SparseColMat<usize, f64>,
        rhs: ColRef<f64>,
        options: &SolverOptions,
    ) -> Result<Col<f64>, SparseError>;

    /// Refactor with same sparsity pattern (reuses symbolic analysis).
    /// Equivalent to calling `factor()` again — provided for API clarity.
    pub fn refactor(
        &mut self,
        matrix: &SparseColMat<usize, f64>,
        options: &FactorOptions,
    ) -> Result<(), SparseError>;

    /// Get inertia from the factorization.
    pub fn inertia(&self) -> Option<Inertia>;

    /// Get factorization statistics.
    pub fn stats(&self) -> Option<&FactorizationStats>;
}

pub struct SolverOptions {
    pub ordering: SymmetricOrdering<'static, usize>,
    pub threshold: f64,
    pub fallback: AptpFallback,
}
```

**Testing strategy:**

- **Per-supernode index mapping**: unit tests for gather/scatter/local_perm
  operations with hand-computed expected values **before** integrating into full solve
- **Known solutions**: construct b = A x_exact, solve, verify x ≈ x_exact
- **Backward error**: ||Ax - b|| / (||A|| ||x|| + ||b||) < 10^-10
  using existing NumericalValidator
- **Backward error as supplementary oracle for factorization** *(added after Phase 5
  SPRAL comparison)*: Phase 5 validates factorization via reconstruction error
  (||A - P^T L D L^T P|| / ||A||). SPRAL instead validates via backward error
  through the full pipeline (factor → solve → check residual), which is a more
  end-to-end test. This phase should add backward error checks to all existing
  Phase 5 test matrices (hand-constructed, SuiteSparse CI, random) as a
  supplementary oracle that validates the full factor+solve pipeline, not just
  factorization in isolation.
- **Scaling round-trip**: verify scaling/unscaling is exact inverse on all matrices
  that use MC64 ordering
- **Multiple RHS**: solve with several right-hand sides reusing same factorization
- **All test matrices**: backward error on hand-constructed + SuiteSparse
- **API ergonomics**: one-shot solve, analyze-factor-solve, refactor
- **Error handling**: solve before factor, pattern mismatch on refactor
- **Documentation examples**: compile and run
- **SuiteSparse**: all easy/hard indefinite matrices, backward error, factor success
- **Random matrices** (via Phase 0.5 generators): backward error
- **Comprehensive report**: using existing test infrastructure

**Success Criteria:**
- [ ] Backward error < 10^-10 on all test matrices
- [ ] Backward error oracle applied to all Phase 5 test matrices (factorization
  validated end-to-end, not just via reconstruction)
- [ ] Per-supernode index mapping correct (gather/scatter/local_perm unit tested)
- [ ] Permutation and scaling correctly applied/unapplied
- [ ] Forward and backward traversals use correct supernode ordering
- [ ] D solve per-supernode via `FrontFactors::d11().solve_in_place()`
- [ ] Multiple RHS handled efficiently (reusing factorization)
- [ ] Three-phase `SparseLDLT` API working end-to-end
- [ ] MemStack-based workspace for solve (no heap allocation in solve hot path)
- [ ] One-shot convenience method works (allocates internally)
- [ ] Refactoring with same sparsity pattern works
- [ ] Error messages clear and actionable
- [ ] >95% of test matrices solve successfully
- [ ] Median backward error < 10^-9
- [ ] Inertia matches reference when available

**Time Estimate:** 3–4 weeks

### Phase 7 Exit Criteria

**Required Outcomes:**
1. Complete solve pipeline working: analyze → factor → solve
2. Backward error meets tolerance on all test matrices
3. Permutation and optional scaling correctly integrated
4. Per-supernode forward/backward solve with correct index mapping
5. Complete solver API working end-to-end via `SparseLDLT`
6. API follows faer's MemStack convention (solve only; factor allocates internally)
7. Integration tests comprehensive and passing

**Validation Questions:**
- Does solve produce correct solutions on all test matrices?
- Is backward error acceptable on hard indefinite problems?
- Are per-supernode gather/scatter/local_perm operations unit tested?
- Does refactoring with same sparsity pattern work correctly?
- Can users solve problems with a simple API call?
- Are error messages clear and actionable?

**Checkpoint:** Run full test suite through public API. Generate comprehensive
test report. This is the "working solver" milestone.

---

## Phase 8: Performance Optimization (8.1a-f COMPLETE)

### Objectives
Optimize the working solver for performance: two-level blocking for large fronts,
shared-memory parallelism for independent subtree processing, and benchmarking
against reference implementations.

### Motivation
The solver from Phases 2–7 is correct but sequential and uses single-level APTP
blocking. This phase adds:
- Cache-efficient blocking for large frontal matrices (two-level APTP)
- Shared-memory parallelism for independent subtree processing
- Performance benchmarking and comparison

**NOTE**: METIS integration was originally deferred to this phase but has been
elevated to Phase 4.1 after Phase 3 testing revealed AMD produces catastrophically
poor orderings on many benchmark matrices (2-20× more fill than METIS). See Phase 3
"Lessons Learned" and Phase 4.1 for details.

### Design Decisions

#### Reuse faer's `Par` for parallelism control (decided)

faer controls parallelism via the `Par` enum (`Par::Seq` for sequential,
`Par::rayon(nthreads)` for parallel). Our factor and solve methods should accept
a `Par` parameter rather than inventing a new parallelism mechanism. This is
consistent with transparent composition: users set parallelism the same way they
do for faer operations. Internally, we use Rayon for assembly-tree level-set
scheduling (new infrastructure), but the user-facing control surface is faer's `Par`.

### Deliverables

#### 8.1: Two-Level APTP & BLAS-3 Factorization (Deferred from Phase 5)
**Task:** Implement nested blocking for the dense APTP kernel and refactor the
multifrontal assembly loop for BLAS-3 performance

**Algorithm Reference:**
- Hogg, Duff, Lopez (2020) - Section 4.2 "Two-level APTP"

**Motivation:**
The single-level APTP kernel from Phase 5 processes columns one at a time within
each block. For small frontal matrices this is sufficient, but large fronts (hundreds
of rows) benefit from cache-efficient nested blocking. Two-level APTP applies the
APTP algorithm recursively: an outer loop processes large blocks, and within each
block an inner loop processes smaller sub-blocks, maximizing use of faer's blocked
matmul.

**BLAS-3 scope note:** Phase 6 passes the entire frontal matrix (F11 + F21 + F22)
to `aptp_factor_in_place()`, relying on the Phase 5 kernel's implicit column-by-column
Schur complement propagation to update all trailing rows. This avoids a separate
TRSM for L21 and GEMM for the Schur complement, but uses BLAS-2 rank-1/rank-2
updates instead of BLAS-3 blocked operations. Phase 8.1 should also evaluate
refactoring the assembly loop to do explicit TRSM (`L21 = F21 * L11^{-T}`) and
GEMM (`F22 -= L21 * D * L21^T`) per outer block, which is the natural decomposition
for BLAS-3 and the prerequisite for parallelizing within a supernode.

**Notional API:**

```rust
pub struct TwoLevelAptpOptions {
    /// Outer block size (e.g., 256). Processes this many columns at a time.
    pub outer_block_size: usize,
    /// Inner block size (e.g., 32). Sub-blocks within each outer block.
    pub inner_block_size: usize,
    /// Stability threshold (same semantics as single-level).
    pub threshold: f64,
    /// Fallback strategy (same as single-level).
    pub fallback: AptpFallback,
}

/// Two-level APTP: same interface as single-level, different internal blocking.
pub fn two_level_aptp_factor(
    a: MatRef<f64>,
    options: &TwoLevelAptpOptions,
) -> Result<AptpFactorization, SparseError>;
```

**Testing:**
- Correctness: two-level must match single-level results (same delays, same
  reconstruction error) on random matrices and hand-constructed cases
- Accuracy: CI SuiteSparse suite (MatchOrderMetis ordering) — target all 8 matrices
  below SPRAL's 5e-11 backward error threshold (Phase 7 baseline: 4/8 pass)
- Performance: benchmark single-level vs two-level on fronts of size 64, 128, 256,
  512, 1024, 2048 (informed by profiling: 41 matrices have max_front 1001–5000)
- Full SuiteSparse evaluation (67 matrices, MatchOrderMetis) after two-level is
  working — sparsine (max_front 11,125) is the ultimate stress test
- Identify crossover point where two-level blocking wins

**Revisiting `matmul` on APTP kernel:** The inner dense APTP kernel uses an O(n^2) naive
operation rather than dispatching to `faer::linalg::matmul`. Phase 7 profiling shows
that single-level factorization of fronts with max_front > 1000 dominates total time
(e.g., ncvxqp3: 53s for max_front 2447). The BLAS-3 refactoring in this phase
(explicit TRSM/GEMM per outer block) is the primary fix; the inner kernel's BLAS-2
operations matter only for the inner block size (~32), where dispatch overhead may
dominate any benefit from vectorized matmul.

**Implementation notes from Phase 5 SPRAL comparison:**

The following topics must be addressed during Phase 8.1 research/speccing:

1. **Backup strategy must change from per-column to per-block**: Phase 5's
   single-level APTP backs up one column at a time (simple `Vec<f64>` in
   `try_1x1_pivot` and `try_2x2_pivot`). Two-level APTP processes outer blocks
   of `nb` columns at a time. SPRAL's `APP_BLOCK` mode backs up one block column
   before optimistic elimination, then restores on failure. SPRAL also has
   `APP_AGGRESSIVE` (back up entire matrix, eliminate all, restore on global
   failure) — decide which mode(s) to implement. See SPRAL's `CopyBackup<T>`
   and `PoolBackup<T>` classes in `ldlt_app.cxx`.

2. **Factor/Apply/Update decomposition**: The Duff, Hogg & Lopez (2020) paper
   Section 4.2 describes two-level APTP as three phases per outer block:
   Factor (inner APTP on the diagonal block), Apply (update the panel below
   the diagonal block using the newly factored portion), Update (rank-nb
   Schur complement on the trailing submatrix). This decomposition is what
   enables BLAS-3 performance and is also the natural parallelism boundary
   (Apply and Update can be parallelized). Phase 8.1 research should map
   this to our existing internal function structure.

3. **Front size profiling results (Phase 7 data — skip criteria resolved)**:
   Full SuiteSparse profiling (80 matrices, METIS ordering) shows two-level APTP
   is **strongly justified**:

   | Category | Count | Description |
   |----------|------:|-------------|
   | Simplicial | 24 | No supernodes (hand-constructed, tridiagonal, etc.) |
   | max_front ≤ 256 | 1 | Single-level fine |
   | max_front 257–1000 | 10 | Borderline — two-level helps modestly |
   | max_front 1001–5000 | **41** | Needs two-level for performance |
   | max_front > 5000 | 4 | Critical for two-level |

   **51 of 56 supernodal matrices** (91%) have max_front > 256. The top 10 by max
   front size: sparsine (11,125), H2O (9,258), nd12k (7,387), Si5H12 (5,452),
   nd6k (4,369), Si10H16 (4,239), apache2 (4,132), c-big (4,064), offshore (3,411),
   shipsec5 (3,249).

   However, most fronts in any given matrix are small: **median front size is
   typically 4–70**, meaning the vast majority of supernodes need no blocking.
   The performance win is concentrated in a long tail of large fronts at the
   top of the elimination tree.

   **Block size recommendations based on profiling**:
   - **Outer block size 256**: Captures 91% of matrices that benefit from blocking.
     Fronts ≤ 256 stay single-level (no overhead).
   - **Inner block size 32**: Fits L1 cache for sub-block operations. 32×32 blocks
     give good BLAS-3 utilization without excessive loop overhead.
   - **Threshold for two-level dispatch**: Apply two-level only when front_size > 256.
     Below this, single-level kernel is adequate (confirmed by Phase 7 benchmarks
     showing excellent accuracy on small fronts).

   **Accuracy note**: Phase 7 benchmarks show backward error correlates with delay
   count, not front size per se. Two-level APTP should improve accuracy on the
   WARN/FAIL matrices (bratu3d 1e-9, cvxqp3 1e-7, stokes128 6e-6) by reducing
   accumulated rounding error within large fronts. Target: all CI matrices below
   SPRAL's 5e-11 threshold with MatchOrderMetis ordering.

**Success Criteria:**
- [x] Two-level APTP produces same factorization quality as single-level
- [ ] Performance improvement on large frontal matrices (n > ~128) — **deferred**: current implementation uses complete pivoting + threshold checking within factor_inner, which processes ALL rows (not just diagonal block). BLAS-3 Apply/Update separation deferred to future refactoring.
- [ ] Cache-friendly memory access patterns verified via profiling — **deferred**: requires BLAS-3 refactoring
- [x] Per-block backup/restore implemented (not per-column) — BlockBackup struct implemented; not integrated because factor_inner handles failures internally via try_1x1/try_2x2

**Phase 8.1 Completion Notes (2026-02-17):**
- **Branch**: `017-two-level-aptp`
- **Block sizes**: outer_block_size=256 (default), inner_block_size=32 (default)
- **Dispatch**: `aptp_factor_in_place` dispatches to `two_level_factor` when `num_fully_summed > outer_block_size`, else to `factor_inner` directly
- **Architecture**: `factor_inner` uses complete pivoting (Algorithm 4.1, Duff et al. 2020) within `factor_block_diagonal`. `two_level_factor` calls `factor_inner` on submatrix views and propagates row permutations to already-factored columns. `factor_inner` processes the entire tile as one block (matching SPRAL's `block_ldlt`), then calls `apply_and_check` (TRSM) for panel rows and `update_trailing` (GEMM) for the trailing matrix.
- **Key finding**: Submatrix views in the outer loop require explicit row permutation propagation — `swap_symmetric` within a submatrix view doesn't reach previously-factored columns. Without this fix, reconstruction errors are O(1).
- **Backward error fix**: Initial implementation broke the tile into `inner_block_size=32` sub-blocks, restricting pivot search to 32 rows. This caused severe backward error regressions on hard indefinite matrices (bratu3d: 8.32e-4 vs 1e-9 single-level). Fix: process the entire tile as one block (`block_size = end_pos - k`), matching SPRAL's `block_ldlt` which searches all tile rows. See `docs/two-level-backward-error-investigation.md` for full details.
- **BLAS-3 status**: `apply_and_check` (TRSM) and `update_trailing` (GEMM) are called in the hot path but use manual loops instead of faer's optimized BLAS-3 kernels. Replacing with `faer::linalg::triangular_solve` and `faer::linalg::matmul::matmul` is straightforward and planned as immediate follow-up.
- **Accuracy**: All 336 lib tests pass. Reconstruction error < 1e-12 on all test matrices. Two-level produces equivalent backward error to single-level on all tested matrices (bratu3d 1e-9, stokes128 1.38e-9, bloweybq 4.27e-11).
- **Single-level removal**: Old single-level main loop in `aptp_factor_in_place` removed. `factor_inner` is the new inner kernel (subsumes old single-level behavior for small fronts).

**Inner blocking within tiles (future consideration):**
SPRAL uses `INNER_BLOCK_SIZE=32` within each tile for minor cache locality, but the
search scope still covers all tile rows. Re-adding inner blocking to our `factor_inner`
would require solving the backup/restore problem for cross-block swaps (when a pivot
candidate is found beyond the current sub-block). The investigation in
`docs/two-level-backward-error-investigation.md` documents the approaches tried and
why they failed. The performance benefit is small compared to BLAS-3 dispatch for
TRSM/GEMM. Consider revisiting after Phase 8.2 profiling identifies bottlenecks.

**Time Estimate:** 1 week

**Post-8.1a sub-phases (8.1b-f):** After the initial two-level APTP implementation
(8.1a), a series of debugging, validation, and pipeline-hardening sub-phases were
carried out on branch `017-two-level-aptp`. These are documented in `docs/ssids-log.md`:

- **8.1b**: BLAS-3 refactoring of `factor_inner` (factor_block_diagonal → apply_and_check
  → update_trailing) + comprehensive accuracy audit. Critical `extract_front_factors`
  bug fix (L21 extraction from wrong row range → 10^12-14× backward error improvement).
- **8.1c**: MC64 Dijkstra heap bug fix — replaced Rust `BinaryHeap` with lazy deletion
  with SPRAL's exact indexed binary heap (`heap_update_inline`/`heap_delete_inline`).
  Eliminated dual feasibility violations in the Hungarian matching algorithm.
- **8.1d**: Fix `col_order` tracking in `two_level_factor` — delay swap must happen
  before `block_perm` update to maintain correct column-to-original mapping.
- **8.1e**: TPP (Threshold Partial Pivoting) fallback for APTP failed columns —
  when block-scoped complete pivoting cannot find acceptable pivots, TPP retries with
  exhaustive serial search.
- **8.1f**: TPP as primary factorization for small fronts, matching SPRAL's dispatch
  (`ldlt_tpp_factor` for ncol < 32, `block_ldlt` for full aligned blocks). Result:
  **65/65 SuiteSparse matrices pass** with MatchOrderMetis. d_pretok backward error
  improved from 2.50e-6 to 7.21e-19.

#### 8.1g: Sequential Profiling & Optimization

**Task:** Profile the complete solver pipeline, establish performance baselines,
and optimize sequential bottlenecks before adding parallelism.

**Motivation:**
Phase 8.1a-f produced a correct, fully-validated sequential solver (65/65 SuiteSparse
pass). Before adding parallel infrastructure, we need to:
1. Understand where time is spent (which supernodes, which operations)
2. Identify and fix allocation hotspots that would cause contention under threads
3. Establish clean sequential baselines for measuring parallel speedup
4. Make informed decisions about which parallelism strategy matters most

Without this data, parallelism decisions are guesswork. Profiling reveals whether
tree-level parallelism (many independent supernodes) or intra-node BLAS-3 parallelism
(few large supernodes) is more important for a given matrix class.

**Deliverables:**

1. **Instrumented profiling of factorization loop**
   - Add `ProfileSession` instrumentation to per-supernode factor loop
   - Measure: time per supernode, frontal matrix allocation, assembly, APTP kernel,
     contribution extraction
   - Generate Chrome Trace output for representative matrices (small, medium, large)
   - Identify top-10 supernodes by time for each CI matrix

2. **Allocation audit and optimization**
   - Per-supernode `Mat::zeros(m, m)` frontal matrix: evaluate reusable scratch buffer
   - Per-row `Vec` allocations in `factor_inner` row permutation loops: hoist outside loop
   - `BlockBackup::create` dense copies: evaluate diagonal-only backup
   - Per-column backup `Vec<f64>` in `try_1x1_pivot`/`try_2x2_pivot`: pre-allocate and reuse
   - Measure allocation count and RSS before/after fixes

3. **Sequential performance baselines**
   - Factor time vs matrix dimension (full SuiteSparse suite)
   - Solve time vs matrix dimension
   - Peak RSS per matrix
   - Time breakdown: ordering / symbolic / numeric / solve
   - Comparison with faer's built-in simplicial LDL^T (where applicable)

4. **Performance report**
   - Workload distribution across supernodes (histogram of per-supernode time)
   - Front size distribution vs time contribution
   - Allocation pressure analysis (count, total bytes, per-supernode)
   - Recommendations for Phase 8.2 parallelism strategy based on profiling data

**Testing:**
- All optimizations must preserve correctness (65/65 SuiteSparse, reconstruction < 1e-12)
- No regressions in backward error on any CI matrix
- Memory usage must not increase (should decrease)

**Success Criteria:**
- [ ] Profiling instrumentation in factorization hot path
- [ ] Chrome Trace profiles generated for CI matrix suite
- [ ] Allocation hotspots identified and top candidates fixed
- [ ] Sequential performance baselines documented
- [ ] Workload distribution analysis complete (informs 8.2 parallelism strategy)
- [ ] No correctness regressions

**Time Estimate:** 1-2 weeks

#### 8.2: Parallel Factorization & Solve
**Task:** Add shared-memory parallelism to factorization and solve, informed by
Phase 8.1g profiling data

**Algorithm Reference:**
- Duff & Reid (1983) — assembly tree parallelism
- Liu (1992) — level-set scheduling for multifrontal methods

**Approach — Parallel Factorization:**
Independent subtrees in the assembly tree can be factored in parallel.
Use level-set parallelism (all supernodes at the same tree level are independent)
or work-stealing via Rayon's `par_iter`. Accept faer's `Par` for user-facing
parallelism control; use Rayon internally for tree-level scheduling.

**Approach — Parallel Solve:**
Level-set parallelism for independent supernodes at the same tree level.
Limited parallelism due to data dependencies, but effective for wide trees.
Use faer's dense TRSM parallelism within large supernodes via `Par`.

**Notional API:**

```rust
impl AptpNumeric {
    /// Factor with parallelism control.
    /// Par::Seq gives sequential; Par::rayon(n) uses n threads.
    pub fn factor(
        symbolic: &AptpSymbolic,
        matrix: &SparseColMat<usize, f64>,
        options: &AptpOptions,
        par: Par,
        stack: &mut MemStack,
    ) -> Result<Self, SparseError>;
}

impl SparseLDLT {
    /// Factor and solve accept Par for parallelism control.
    pub fn factor(
        &mut self,
        matrix: &SparseColMat<usize, f64>,
        options: &FactorOptions,
        par: Par,
        stack: &mut MemStack,
    ) -> Result<(), SparseError>;

    pub fn solve_in_place(
        &self,
        rhs: ColMut<f64>,
        par: Par,
        stack: &mut MemStack,
    ) -> Result<(), SparseError>;
}
```

**NOTE**: The `Par` parameter will also be added to Phase 6–7 APIs retroactively
(with `Par::Seq` default until Phase 8 adds the parallel implementation). The
exact threading boundary — which methods accept `Par` — will be refined during
speccing.

**Benchmarking deliverable:** Parallel scaling report (building on 8.1g baselines):
- Parallel scaling curves (1, 2, 4, 8 threads) for factor and solve
- Speedup vs sequential baselines from Phase 8.1g
- Load balancing analysis (workload imbalance across threads)
- Comparison with faer's built-in solvers and SPRAL (if available)

**Testing:**
- Parallel results must be deterministic
- Parallel results must match sequential results exactly (bitwise)
- Speedup observed on large problems

**Success Criteria:**
- [ ] Deterministic parallel execution
- [ ] Results identical to sequential (bitwise)
- [ ] Speedup on factorization (≥3× on 4 cores for large problems)
- [ ] Speedup on solve for wide assembly trees
- [ ] Load balancing effective for unbalanced trees
- [ ] ≥75% parallel efficiency on 4 cores
- [ ] No correctness regressions

**Time Estimate:** 2-3 weeks

### Phase 8 Exit Criteria

**Required Outcomes:**
1. Two-level APTP correct and validated on full SuiteSparse (8.1a-f: DONE)
2. Sequential performance profiled, baselined, and optimized (8.1g)
3. Parallel factorization and solve working and deterministic (8.2)
4. No correctness regressions from optimization or parallelism
5. All CI SuiteSparse matrices (MatchOrderMetis) achieve backward error < 5e-11

**Validation Questions:**
- Where does sequential time go? (8.1g profiling → informs 8.2 strategy)
- Is parallel speedup significant? (8.2)
- Is performance competitive with SPRAL? (8.1g sequential + 8.2 parallel)

**Checkpoint:** Run full SuiteSparse benchmark suite (67+ matrices, MatchOrderMetis).
Generate performance report with sequential baselines (8.1g) and parallel scaling
plots (8.2). Compare Phase 7 baseline (4/8 CI pass) vs Phase 8 results.

**Time Estimate:** 4-5 weeks (8.1g: 1-2 weeks, 8.2: 2-3 weeks)

---

## Phase 9: Polish & Release

### Objectives
Harden the solver for production use (memory optimization, robustness testing) and
prepare for public release (documentation, examples, packaging).

### What was absorbed or superseded from original Phases 10–11

- **Old 10.1 (Performance Profiling)**: Profiling tools already built in Phase 1.4.
  Benchmarking already covered by Phase 8.2. Targeted bottleneck fixes are part of
  9.1 below.
- **Old 10.2 (Memory Optimization)**: Workspace allocation addressed by Phase 7's
  MemStack pattern. Compact factor storage addressed by Phase 2's MixedDiagonal and
  Phase 6's FrontFactors. Arena allocation for frontal matrices remains (9.1 below).
- **Old 10.3 (API Refinement)**: Superseded by Phase 7's SparseLDLT API design and
  Phase 8's Par integration. Builder pattern is minor sugar, deferred to 9.2 if needed.
- **Old 11.3 (Benchmarks)**: Fully superseded by Phase 8.2's performance report.
- **Old 11.4 (CI/CD)**: Already done in Phase 1.3.

### Deliverables

#### 9.1: Solver Hardening
**Task:** Improve robustness and memory efficiency of the working solver

**Memory optimization:**
- Arena allocation for frontal matrices during factorization — allocate from a
  pre-sized memory pool rather than individual heap allocations per front. This
  complements Phase 7's MemStack (which handles solve-phase workspace) by
  optimizing the factorization-phase memory pattern.
- Peak memory tracking and validation against symbolic predictions
- Memory usage within 50% of prediction from AptpSymbolic

**SPRAL testing parity:**
- Perform a comprehensive analysis of SPRAL's test suite
- Ensure that all SPRAL tests map conceptually to one of our unit tests if feasible
- Review our test suite for unnecessary or duplicate tests, or tests that were used as part of a TDD protocol but are no longer independently necessary

**Property-based testing (proptest):**
- Generate random symmetric matrices (varying size, density, definiteness)
- Verify backward error < tolerance for all solvable systems
- Verify no panics on singular or near-singular matrices

**Fuzzing:**
- Malformed inputs (invalid sparsity patterns, non-symmetric matrices)
- Extreme values (near-overflow, near-underflow, exact zeros)
- Singular and structurally singular matrices
- Goal: no panics, only clean error returns

**Iterative refinement:**
- Implement one or two steps of iterative refinement after the initial solve:
  1. Compute residual `r = b - A*x` (using sparse matvec)
  2. Solve `A*z = r` using the existing factorization (forward/diagonal/backward)
  3. Update `x += z`
  4. Optionally repeat until `||r|| / (||A||*||x|| + ||b||) < tol`
- SPRAL achieves backward error < 5e-11 on bratu3d; our current solver achieves
  1e-9. The ~20× gap is very likely explained by iterative refinement recovering
  1-2 digits of accuracy per step.
- API: add `SolveOptions { max_refinement_steps: usize, refinement_tol: f64 }`
  to control refinement. Default: 2 steps (matching SPRAL). Zero steps gives
  the raw solve for users who want maximum speed.
- The factorization is already computed, so each refinement step costs only
  O(nnz) for the sparse matvec + O(n) for the triangular solves — negligible
  compared to the O(n*nnz) factorization cost.
- Decision criterion: implement if any CI matrices remain above SPRAL's 5e-11
  threshold after BLAS-3 dispatch is complete (Phase 8.1 follow-up).

**Targeted performance fixes:**
- Use Phase 1.4 profiling tools and Phase 8.2 benchmark results to identify
  top 3 bottlenecks
- Optimize allocation patterns, cache utilization, unnecessary copies
- Verify no performance regressions

**Items from Phase 5 SPRAL comparison** *(added after Phase 5 completion)*:

*TPP (Traditional Threshold Pivoting) fallback*:
SPRAL has a TPP fallback (`ldlt_tpp_factor` in `ldlt_app.cxx`) — when APTP delays
too many columns in a frontal matrix, SPRAL can re-factorize the delayed portion
using traditional threshold pivoting rather than propagating all delays to the
parent node. Our Phase 5 kernel returns delayed columns to the caller (Phase 6
multifrontal assembly), which passes them up the elimination tree. If Phase 8.2
benchmarks reveal excessive delay propagation (e.g., >20% of columns delayed on
>10% of frontal matrices across the benchmark suite), consider adding a TPP
fallback kernel. This would be a new function alongside `aptp_factor_in_place`
that uses traditional pivoting on the delayed submatrix before returning results.
Decision criterion: add TPP if delay propagation measurably increases fill
(>15% more nonzeros in parent fronts) or worsens backward error (>10x worse)
compared to a hypothetical TPP-enabled run.

*D storage format optimization*:
Phase 5 uses `MixedDiagonal` with a `PivotType` enum and parallel arrays
(pivot_map, diag, off_diag). SPRAL uses a flat `d[2n]` array with an `Inf`
sentinel to mark 2x2 block boundaries, enabling O(1) pivot-type detection during
solve without a separate type array. If Phase 8.2 benchmarks show the triangular
solve phase is >25% of total solve time, profile whether D access is a bottleneck.
If so, consider adding a compact `pack_d()` method that produces a flat array for
the solve phase while keeping `MixedDiagonal` for the factorization phase.

*Torture testing (SPRAL-style stress tests)*:
SPRAL has comprehensive torture testing for the APTP kernel that goes beyond our
current random matrix stress tests. Add SPRAL-style perturbation functions:
- `cause_delays()`: randomly multiply n/8 rows by 1000 to force pivot delays
- `make_singular()`: scale one column and copy to another to create rank deficiency
- `make_dblk_singular()`: make a specific diagonal block singular
- Probabilistic test generation: ~70% chance of delays, ~20% singularity, ~10%
  singular diagonal blocks (following SPRAL's distribution)
- Generate 500+ random instances per configuration
Reference: `spral/tests/ssids/kernels/ldlt_app.cxx` lines 78-134 (perturbation
helpers), lines 423-451 (`ldlt_torture_test` with probabilistic generation),
and line 593 (500 tests of 128x128 matrices). Also `ldlt_tpp.cxx` lines 343-369
for the TPP torture test variant. These tests caught subtle numerical issues in
SPRAL that deterministic tests missed.

**MC64 scaling quality revisit:**
- Phase 4.2's MC64 matching produces degraded row_max quality (rows where
  `max_j |s_i·a_ij·s_j| < 0.75`) on 5/47 fully-matched SuiteSparse matrices
  (TSOPF_FS_b39_c7, d_pretok, nd12k, ship_003, thread). This is an inherent
  limitation of the MC64SYM symmetric averaging formula (Duff & Pralet 2005),
  not an implementation bug — SPRAL uses the same formula.
- With end-to-end solve benchmarks available from Phase 8.2, measure whether
  the quality degradation causes measurably more delayed pivots or worse
  backward error on these matrices compared to identity scaling.
- If impact is significant: investigate maintaining column duals during Dijkstra
  (as SPRAL does) to produce more singleton-dominated matchings with shorter
  cycles and less averaging error. See `docs/mc64-scaling-notes.md` for the
  full analysis.

**Success Criteria:**
- [ ] Arena allocation reduces peak memory on large problems
- [ ] Peak memory within 50% of symbolic prediction
- [ ] Property-based tests pass on random matrices (size 5–500)
- [ ] No panics on any invalid input (clean error returns)
- [ ] >90% code coverage on core solver code (`aptp/`, `error.rs`, `validate.rs`)
- [ ] Top bottlenecks identified and addressed

**Time Estimate:** 2 weeks

#### 9.2: Release Preparation
**Task:** Documentation, examples, and packaging for public release

**Documentation:**
- README with quick start guide
- User guide covering the three-phase API, options, and common patterns
- API reference (rustdoc with examples on key types and methods)
- Performance guide (when to use APTP vs faer's built-in solvers, tuning options)
- Migration guide for users coming from SPRAL, MUMPS, or HSL solvers

**Examples** (standalone `examples/` directory):
- `basic_solve.rs` — one-shot solve with `SparseLDLT::solve_full`
- `three_phase.rs` — analyze/factor/solve with workspace reuse
- `multiple_rhs.rs` — solving with several right-hand sides
- `refactorization.rs` — reusing symbolic analysis for updated values
- `interior_point.rs` — integration with an interior point optimization loop

**Packaging:**
- Version 0.1.0 on crates.io
- Issue templates and contributing guide
- Verify CI passes on Linux/macOS (Windows stretch goal)

**Success Criteria:**
- [ ] README and user guide complete
- [ ] All examples compile and run correctly
- [ ] Rustdoc examples on public types and methods
- [ ] Published on crates.io
- [ ] Contributing guide and issue templates in place

**Time Estimate:** 2 weeks

### Phase 9 Exit Criteria

**Required Outcomes:**
1. Solver robust against malformed inputs (no panics)
2. Memory usage predictable and efficient
3. Documentation and examples complete
4. Published on crates.io

**Validation Questions:**
- Is the library ready for public use?
- Can a new user solve a problem from the README alone?
- Are there any remaining robustness issues?

**Checkpoint:** Run full test suite including property-based tests. Verify all
examples compile and run. Publish to crates.io.

**Time Estimate:** 4 weeks

---

## Success Metrics (Overall Project)

### Numerical Quality
- [ ] >95% of test matrices solve successfully
- [ ] Median residual < 1e-9
- [ ] Inertia computation accurate
- [ ] Stable on hard indefinite problems

### Performance
- [ ] Sequential performance within 2× of SPRAL
- [ ] Parallel speedup: 3× on 4 cores
- [ ] Memory usage within 50% of prediction

### Code Quality
- [ ] >90% test coverage
- [ ] No clippy warnings
- [ ] Comprehensive documentation
- [ ] Examples for common use cases

### Ecosystem Integration
- [ ] Works with faer types
- [ ] Published on crates.io
- [ ] CI passing on Linux/macOS/Windows
- [ ] Ready for community contributions

---

This specification provides a clear roadmap from foundation to release, with well-defined checkpoints and success criteria at each phase.
