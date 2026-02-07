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
- **Leverage faer infrastructure** (~70% reuse): CSC storage, elimination trees, AMD ordering, permutation utilities, workspace management
- **Build APTP-specific components** (~30% new): 2x2 pivot logic, mixed diagonal storage, indefinite factorization kernel
- **Simplicial first, supernodal later**: Start with column-by-column factorization (Phases 2-8), add supernodal optimization in Phase 9
- **Clean room implementation**: All code derived from BSD-licensed references (LAPACK, SLICOT-Reference, SPRAL) and academic papers

**Development Timeline:**
- Phases 0-8: Simplicial APTP solver (~12-16 weeks with faer reuse)
- Phase 9: Supernodal optimization (~4-6 weeks)
- Phases 10-11: Polish and release (~4-6 weeks)
- **Total**: ~20-28 weeks to production-ready solver

---

## Phase 0: Foundation & Literature Review

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
(Phases 2-8), for performance benchmarking and inertia validation on large
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

## Phase 1: Testing & Benchmarking Infrastructure

### Objectives
Build production-quality testing and benchmarking framework before implementing any solver components. This infrastructure will be used throughout all subsequent phases.

### Deliverables

#### 1.1: Core Test Infrastructure
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

#### 1.2: Benchmarking Framework
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

#### 1.4: Profiling and Debug Tools
**Task:** Build tools for performance analysis and debugging

**Profiling Helpers:**

```rust
pub struct ProfileRecorder {
    events: Vec<ProfileEvent>,
}

impl ProfileRecorder {
    pub fn record_section<F, R>(&mut self, name: &str, f: F) -> R
    where F: FnOnce() -> R;

    pub fn record_memory_snapshot(&mut self, label: &str);

    pub fn export_chrome_trace(&self) -> String;
    pub fn export_flamegraph(&self) -> String;
}

// Usage in implementation:
fn factor_node(&mut self, node: NodeId, profiler: &mut ProfileRecorder) {
    profiler.record_section("dense_factor", || {
        self.aptp_factor(node)
    });

    profiler.record_section("update_ancestors", || {
        self.propagate_updates(node)
    });
}
```

**Debug Visualization:**

```rust
pub struct DebugVisualizer;

impl DebugVisualizer {
    /// Generate GraphViz DOT file of elimination tree
    pub fn visualize_elimination_tree(tree: &EliminationTree) -> String;

    /// Generate sparsity pattern visualization
    pub fn visualize_sparsity(matrix: &SparseMatrix) -> Image;

    /// Animate factorization process
    pub fn animate_factorization(steps: &[FactorStep]) -> Animation;
}
```

**Memory Debugging:**

```rust
#[cfg(feature = "memory-debug")]
pub struct MemoryTracker {
    allocations: HashMap<String, AllocationInfo>,
}

impl MemoryTracker {
    pub fn track_allocation(&mut self, label: &str, bytes: usize);
    pub fn generate_report(&self) -> MemoryReport;
    pub fn detect_leaks(&self) -> Vec<Leak>;
}
```

**Success Criteria:**
- [ ] Can profile any component in isolation
- [ ] Memory usage traceable to specific operations
- [ ] Visualization tools for debugging
- [ ] Integration with `perf`, `valgrind`, etc.

### Phase 1 Exit Criteria

**Required Outcomes:**
1. Complete test infrastructure can validate any solver component
2. Benchmarking framework ready for use
3. Can measure performance regression automatically
4. Debugging and profiling tools available
5. CI pipeline validates every commit

**Validation Questions:**
- Can we detect a 1% performance regression automatically?
- Can we identify which component caused a test failure?
- Can we profile a single matrix through the entire solve?
- Is the test suite fast enough to run frequently?

**Checkpoint:** Run infrastructure on trivial "pass-through" solver to verify it works.

---

## Phase 2: Sparse Matrix Infrastructure (Leverage faer)

### Objectives
Adopt faer-sparse types with APTP-specific extensions. Minimal reimplementation; focus on what's unique to indefinite APTP.

### Deliverables

#### 2.1: faer-sparse Adoption with Type Aliases
**Task:** Use faer's sparse matrix infrastructure directly

**Approach:**
Leverage faer's mature sparse matrix implementation instead of building from scratch:
- Use `faer::sparse::SparseColMat<I, T>` as primary storage
- Use faer's Matrix Market I/O
- Add only APTP-specific storage structures

**Implementation:**

```rust
// Type aliases for clarity and consistency
pub type AptpMatrix<T = f64> = faer::sparse::SparseColMat<usize, T>;
pub type AptpMatrixRef<'a, T = f64> = faer::sparse::SparseColMatRef<'a, usize, T>;
pub type AptpMatrixMut<'a, T = f64> = faer::sparse::SparseColMatMut<'a, usize, T>;

// Re-export commonly used faer functionality
pub use faer::sparse::SymbolicSparseColMat;

// Matrix Market I/O via faer (or matrixmarket-rs crate)
pub mod io {
    use super::*;

    pub fn read_matrix_market<P: AsRef<Path>>(
        path: P
    ) -> Result<AptpMatrix<f64>> {
        // Use faer or matrixmarket-rs
        todo!()
    }

    pub fn write_matrix_market<P: AsRef<Path>>(
        path: P,
        matrix: AptpMatrixRef<'_, f64>
    ) -> Result<()> {
        todo!()
    }
}
```

**Testing:**
```rust
#[test]
fn test_load_all_test_matrices() {
    for entry in glob("test-data/**/*.mtx").unwrap() {
        let path = entry.unwrap();
        let result = io::read_matrix_market(&path);
        assert!(result.is_ok(), "Failed to load {:?}", path);
    }
}

#[test]
fn test_faer_integration() {
    let mat = io::read_matrix_market("test-data/arrow-10.mtx").unwrap();

    // Verify it's a valid faer matrix
    assert_eq!(mat.nrows(), 10);
    assert_eq!(mat.ncols(), 10);
    assert!(mat.nrows() * mat.ncols() >= mat.compute_nnz());
}
```

**Success Criteria:**
- [ ] Can load all 70+ test matrices using faer
- [ ] Type aliases provide clean API
- [ ] Matrix Market I/O working
- [ ] Basic faer operations accessible

**Time Estimate:** 1-2 days

#### 2.2: APTP-Specific Storage Structures
**Task:** Define data structures unique to indefinite APTP factorization

**Implementation:**

```rust
/// Tracks whether each pivot is 1x1 or 2x2
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PivotType {
    One,        // 1x1 pivot
    Two(usize), // 2x2 pivot, points to block index
}

/// Storage for mixed 1x1 and 2x2 diagonal blocks
pub struct MixedDiagonal<T> {
    /// 1x1 diagonal entries (sparse - only for columns with 1x1 pivots)
    one_by_one: Vec<(usize, T)>,  // (column_index, value)

    /// 2x2 diagonal blocks (stored as [a, b; b, c])
    two_by_two: Vec<Block2x2<T>>,

    /// Map from column index to pivot type
    pivot_map: Vec<PivotType>,
}

/// A 2x2 symmetric block [a, b; b, c]
#[derive(Debug, Clone)]
pub struct Block2x2<T> {
    pub first_col: usize,  // Starting column
    pub a: T,              // (0,0) entry
    pub b: T,              // (0,1) and (1,0) entry
    pub c: T,              // (1,1) entry
}

impl<T: ComplexField> MixedDiagonal<T> {
    pub fn new(n: usize) -> Self;

    pub fn set_1x1(&mut self, col: usize, value: T);
    pub fn set_2x2(&mut self, first_col: usize, block: Block2x2<T>);

    pub fn get_pivot_type(&self, col: usize) -> PivotType;

    /// Solve D*x = b where D is mixed 1x1/2x2
    pub fn solve(&self, b: &mut [T]);
}

/// Permutation wrapper with inversion
pub struct Permutation {
    perm: Vec<usize>,      // forward: old -> new
    inv_perm: Vec<usize>,  // inverse: new -> old
}

impl Permutation {
    pub fn identity(n: usize) -> Self;
    pub fn from_vec(perm: Vec<usize>) -> Result<Self>;

    pub fn apply<T: Clone>(&self, x: &[T]) -> Vec<T>;
    pub fn apply_inverse<T: Clone>(&self, x: &[T]) -> Vec<T>;

    pub fn is_valid(&self) -> bool;
}
```

**Testing:**
```rust
#[test]
fn test_mixed_diagonal() {
    let mut d = MixedDiagonal::<f64>::new(4);

    d.set_1x1(0, 2.0);
    d.set_2x2(1, Block2x2 { first_col: 1, a: 3.0, b: 1.0, c: 4.0 });
    d.set_1x1(3, 5.0);

    assert_eq!(d.get_pivot_type(0), PivotType::One);
    assert_eq!(d.get_pivot_type(1), PivotType::Two(0));
    assert_eq!(d.get_pivot_type(2), PivotType::Two(0));
    assert_eq!(d.get_pivot_type(3), PivotType::One);
}

#[test]
fn test_mixed_diagonal_solve() {
    let mut d = MixedDiagonal::new(3);
    d.set_1x1(0, 2.0);
    d.set_2x2(1, Block2x2 { first_col: 1, a: 4.0, b: 1.0, c: 3.0 });

    let mut x = vec![4.0, 7.0, 5.0];
    d.solve(&mut x);

    // Check D * x = original_b
    // For 1x1: x[0] = 4.0 / 2.0 = 2.0
    // For 2x2: solve [4,1;1,3]*[x1;x2] = [7;5]
    assert_close(x[0], 2.0, 1e-14);
}
```

**Success Criteria:**
- [ ] MixedDiagonal correctly stores mixed 1x1/2x2 blocks
- [ ] Solve with mixed D works correctly
- [ ] Permutation wrapper validated
- [ ] Efficient access patterns

**Time Estimate:** 2-3 days

#### 2.3: Integration with Test Infrastructure
**Task:** Connect faer matrices to test framework

**Implementation:**

```rust
impl SolverTestCase {
    pub fn from_matrix_market<P: AsRef<Path>>(path: P)
        -> Result<Self> {
        let matrix = crate::io::read_matrix_market(path)?;
        // Extract metadata, create test case
        todo!()
    }

    pub fn matrix(&self) -> AptpMatrixRef<'_, f64>;
}

// Test suite accessor using faer matrices
pub fn all_test_cases() -> Vec<SolverTestCase> {
    // Load all matrices from test-data/ using faer
}
```

**Success Criteria:**
- [ ] Test infrastructure works with faer types
- [ ] All test matrices loadable
- [ ] Easy to add new matrices

**Time Estimate:** 1 day

### Phase 2 Exit Criteria

**Required Outcomes:**
1. faer-sparse types integrated and usable
2. APTP-specific structures (MixedDiagonal, pivot tracking) implemented and tested
3. All test matrices loadable via faer
4. Test infrastructure using faer types

**Validation Questions:**
- Can we load all test matrices using faer?
- Do APTP-specific structures work correctly?
- Is the API clean and ergonomic?

**Checkpoint:** Load entire test suite using faer. Verify MixedDiagonal solve correctness on hand-constructed examples.

**Time Estimate:** 3-5 days (reduced from 2-3 weeks by leveraging faer)

---

## Phase 3: Symbolic Analysis (Leverage faer, Defer Supernodes)

### Objectives
Build symbolic analysis using faer's elimination tree construction. Focus on simplicial APTP initially; defer supernodal optimization to Phase 9.

### Architecture Decision
**Start with simplicial (column-by-column) APTP factorization:**
- Simpler data structures
- Easier 2x2 pivot handling
- Faster path to working solver
- Supernodal optimization becomes Phase 9

### Deliverables

#### 3.1: Elimination Tree (Reuse faer)
**Task:** Use faer's elimination tree construction with APTP adaptations

**Approach:**
The elimination tree structure for indefinite LDL^T is **identical** to SPD Cholesky - pivoting is a numeric-phase concern. Leverage faer's mature implementation.

**Algorithm References:**
- Liu (1990) - "The role of elimination trees in sparse factorization"
- faer implementation: `faer::sparse::linalg::cholesky::prefactorize_symbolic_cholesky`

**Implementation:**

```rust
use faer::sparse::linalg::cholesky::{
    prefactorize_symbolic_cholesky,
    EliminationTreeRef,
};

/// APTP symbolic analysis (wraps faer's etree + APTP-specific metadata)
pub struct AptpSymbolic<I = usize> {
    /// Elimination tree from faer (reused directly)
    etree: Vec<I::Signed>,  // parent[i] = parent of i, or -1 if root

    /// Column counts (from faer)
    col_counts: Vec<I>,

    /// APTP-specific: extra space for delayed pivots
    pivot_buffer_size: Vec<I>,

    /// Total predicted NNZ (may be approximate for indefinite)
    predicted_nnz: usize,
}

impl AptpSymbolic {
    pub fn analyze(
        matrix: SymbolicSparseColMatRef<'_, I>
    ) -> Result<Self> {
        // Use faer's prefactorization (computes etree + col_counts)
        let (etree, col_counts) = prefactorize_symbolic_cholesky(
            matrix,
            PodStack::new(&mut GlobalPodBuffer::new(
                StackReq::new::<I>(matrix.ncols())
            )),
        )?;

        // APTP-specific: estimate buffer for delayed pivots
        let pivot_buffer_size = col_counts
            .iter()
            .map(|&c| (c as f64 * 0.1) as I)  // 10% buffer
            .collect();

        let predicted_nnz = col_counts.iter().sum();

        Ok(Self {
            etree: etree.into_inner().to_vec(),
            col_counts: col_counts.to_vec(),
            pivot_buffer_size,
            predicted_nnz,
        })
    }

    pub fn parent(&self, node: usize) -> Option<usize> {
        let p = self.etree[node];
        if p < 0 {
            None
        } else {
            Some(p as usize)
        }
    }

    pub fn predicted_nnz(&self) -> usize {
        self.predicted_nnz
    }
}
```

**Testing:**

```rust
#[test]
fn test_elimination_tree_via_faer() {
    let matrix = io::read_matrix_market("test-data/arrow-10.mtx").unwrap();
    let symbolic = AptpSymbolic::analyze(matrix.symbolic()).unwrap();

    // Tree should be valid
    assert!(symbolic.predicted_nnz() > matrix.compute_nnz());
}

#[test]
fn test_against_spral_reference() {
    for case in test_cases_by_difficulty(Difficulty::Easy) {
        let symbolic = AptpSymbolic::analyze(
            case.matrix().symbolic()
        ).unwrap();

        if let Some(ref_results) = case.reference {
            // Predicted NNZ should be close to SPRAL's
            let error = (symbolic.predicted_nnz() as f64
                        - ref_results.predicted_nnz as f64).abs()
                       / ref_results.predicted_nnz as f64;
            assert!(error < 0.1, "NNZ prediction off by {:.1}%", error * 100.0);
        }
    }
}

#[test]
fn test_reproducible() {
    let matrix = io::read_matrix_market("test-data/symmetric-10x10.mtx").unwrap();

    let sym1 = AptpSymbolic::analyze(matrix.symbolic()).unwrap();
    let sym2 = AptpSymbolic::analyze(matrix.symbolic()).unwrap();

    assert_eq!(sym1.etree, sym2.etree);
    assert_eq!(sym1.col_counts, sym2.col_counts);
}
```

**Success Criteria:**
- [ ] Can construct etree via faer for all test matrices
- [ ] Predicted NNZ within 10% of SPRAL (indefinite may differ slightly)
- [ ] Reproducible (deterministic)
- [ ] Analysis time < 5% of factor time

**Time Estimate:** 2-3 days

#### 3.2: Sparsity Pattern Prediction (Simplicial)
**Task:** Predict L sparsity pattern for simplicial factorization

**Approach:**
For simplicial (column-by-column) factorization, we need to know which entries of L will be nonzero. This is determined by the elimination tree and reachability analysis.

**Implementation:**

```rust
impl AptpSymbolic {
    /// Allocate storage for L factor
    pub fn allocate_factor_storage(&self) -> SimplicialLFactor {
        // L has same structure as lower triangle of symbolic Cholesky factor
        let total_nnz = self.predicted_nnz;

        SimplicialLFactor {
            col_ptr: vec![0; self.etree.len() + 1],
            row_idx: Vec::with_capacity(total_nnz),
            values: Vec::with_capacity(total_nnz),
        }
    }
}

/// Storage for simplicial L factor (CSC format, lower triangular)
pub struct SimplicialLFactor {
    pub col_ptr: Vec<usize>,
    pub row_idx: Vec<usize>,
    pub values: Vec<f64>,
}

impl SimplicialLFactor {
    pub fn ncols(&self) -> usize {
        self.col_ptr.len() - 1
    }

    pub fn col(&self, j: usize) -> &[usize] {
        let start = self.col_ptr[j];
        let end = self.col_ptr[j + 1];
        &self.row_idx[start..end]
    }

    pub fn col_values(&self, j: usize) -> &[f64] {
        let start = self.col_ptr[j];
        let end = self.col_ptr[j + 1];
        &self.values[start..end]
    }
}
```

**Testing:**

```rust
#[test]
fn test_factor_storage_allocation() {
    let matrix = io::read_matrix_market("test-data/arrow-10.mtx").unwrap();
    let symbolic = AptpSymbolic::analyze(matrix.symbolic()).unwrap();

    let l_factor = symbolic.allocate_factor_storage();

    assert_eq!(l_factor.ncols(), 10);
    assert!(l_factor.values.capacity() >= matrix.compute_nnz());
}
```

**Success Criteria:**
- [ ] Storage correctly allocated
- [ ] Size predictions reasonable
- [ ] Ready for numeric factorization

**Time Estimate:** 2-3 days

#### 3.3: Complete Symbolic Analysis Interface
**Task:** Integrate symbolic components into clean API

**Implementation:**

```rust
/// Options for symbolic analysis
pub struct AnalysisOptions {
    // Future: could add relaxation params when supernodes added
}

impl Default for AnalysisOptions {
    fn default() -> Self {
        Self {}
    }
}

impl AptpSymbolic {
    /// Convenience method for full analysis pipeline
    pub fn analyze_with_options(
        matrix: SymbolicSparseColMatRef<'_, usize>,
        _options: AnalysisOptions,
    ) -> Result<Self> {
        // For now, options unused (no supernodes yet)
        Self::analyze(matrix)
    }

    /// Statistics for debugging
    pub fn statistics(&self) -> SymbolicStatistics {
        SymbolicStatistics {
            dimension: self.etree.len(),
            predicted_nnz: self.predicted_nnz,
            average_col_count: self.predicted_nnz as f64 / self.etree.len() as f64,
        }
    }
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
fn test_complete_analysis() {
    for case in all_test_cases() {
        let symbolic = AptpSymbolic::analyze(
            case.matrix().symbolic()
        ).unwrap();

        let stats = symbolic.statistics();

        // Sanity checks
        assert_eq!(stats.dimension, case.matrix().nrows());
        assert!(stats.predicted_nnz > 0);
        assert!(stats.average_col_count >= 1.0);

        // Compare with reference if available
        if let Some(ref_results) = case.reference {
            let error = (stats.predicted_nnz as f64
                        - ref_results.predicted_nnz as f64).abs()
                       / ref_results.predicted_nnz as f64;
            // Allow 10% error for indefinite (delayed pivots unpredictable)
            assert!(error < 0.1);
        }
    }
}

#[test]
fn test_analysis_reproducible() {
    let case = load_test_case("medium-matrix");

    let sym1 = AptpSymbolic::analyze(case.matrix().symbolic()).unwrap();
    let sym2 = AptpSymbolic::analyze(case.matrix().symbolic()).unwrap();

    assert_eq!(sym1.etree, sym2.etree);
    assert_eq!(sym1.predicted_nnz, sym2.predicted_nnz);
}
```

**Success Criteria:**
- [ ] Complete analysis works on all test matrices
- [ ] Results within 10% of SPRAL (indefinite harder to predict)
- [ ] Reproducible (deterministic)
- [ ] Clean, documented API

**Time Estimate:** 1-2 days

### Phase 3 Exit Criteria

**Required Outcomes:**
1. Symbolic analysis using faer's etree complete and tested
2. Sparsity pattern prediction working for simplicial factorization
3. Integration tests pass on full test suite
4. Performance acceptable (< 5% of total solve time)
5. Predictions within 10% of SPRAL (indefinite is harder to predict exactly)

**Validation Questions:**
- Do we get the same elimination tree as SPRAL?
- Are NNZ predictions reasonable?
- Is factor storage correctly allocated?

**Checkpoint:** Run symbolic analysis on entire test suite. Verify etree matches SPRAL. Check that predicted NNZ is within 10% of SPRAL reference results.

**Time Estimate:** 5-7 days (reduced from 3-4 weeks by leveraging faer and deferring supernodes)

---

## Phase 4: Ordering & Scaling (Partial faer Reuse)

### Objectives
Integrate ordering algorithms: reuse faer's AMD, add METIS wrapper, implement MC64 matching-based ordering for indefinite systems.

### Deliverables

#### 4.1: Ordering Framework (Reuse AMD from faer)
**Task:** Set up ordering interface and leverage faer's AMD implementation

**Approach:**
faer provides AMD; we add METIS and matching-based orderings.

**Implementation:**

```rust
use faer::sparse::linalg::amd;
use metis::Graph;

pub enum OrderingMethod {
    Natural,          // No reordering
    Amd,             // AMD via faer
    Metis,           // Nested dissection (new)
    Mc64Amd,         // MC64 + AMD for indefinite (new)
    Custom(Vec<usize>),
}

pub fn compute_ordering(
    matrix: SymbolicSparseColMatRef<'_, usize>,
    method: OrderingMethod,
) -> Result<Permutation> {
    match method {
        OrderingMethod::Natural => {
            Ok(Permutation::identity(matrix.nrows()))
        }
        OrderingMethod::Amd => {
            // Use faer's AMD directly
            let perm_data = amd::order(matrix, Default::default())?;
            Ok(Permutation::from_faer(perm_data))
        }
        OrderingMethod::Metis => {
            // New implementation
            metis_nested_dissection(matrix)
        }
        OrderingMethod::Mc64Amd => {
            // New: matching + AMD
            mc64_then_amd(matrix)
        }
        OrderingMethod::Custom(perm) => {
            Permutation::from_vec(perm)
        }
    }
}

pub struct Permutation {
    perm: Vec<usize>,      // forward: old -> new
    inv_perm: Vec<usize>,  // inverse: new -> old
}

impl Permutation {
    pub fn identity(n: usize) -> Self {
        let perm: Vec<usize> = (0..n).collect();
        Self {
            inv_perm: perm.clone(),
            perm,
        }
    }

    pub fn from_vec(perm: Vec<usize>) -> Result<Self> {
        let n = perm.len();
        let mut inv_perm = vec![0; n];
        for (new_pos, &old_pos) in perm.iter().enumerate() {
            inv_perm[old_pos] = new_pos;
        }
        Ok(Self { perm, inv_perm })
    }

    pub fn from_faer(perm: faer::perm::PermRef<'_, usize>) -> Self {
        Self::from_vec(perm.arrays().0.to_vec()).unwrap()
    }

    pub fn apply<T: Clone>(&self, x: &[T]) -> Vec<T> {
        self.perm.iter().map(|&i| x[i].clone()).collect()
    }

    pub fn is_valid(&self) -> bool {
        let n = self.perm.len();
        let mut seen = vec![false; n];
        for &p in &self.perm {
            if p >= n || seen[p] {
                return false;
            }
            seen[p] = true;
        }
        true
    }
}

// METIS-specific implementation (new code)
fn metis_nested_dissection(
    matrix: SymbolicSparseColMatRef<'_, usize>
) -> Result<Permutation> {
    // Convert to METIS graph format
    let graph = matrix_to_metis_graph(matrix)?;

    // Call METIS for nested dissection ordering
    let perm = graph.part_nd()?;

    Ok(Permutation::from_vec(perm)?)
}
```

**Testing:**

```rust
#[test]
fn test_metis_ordering() {
    let pattern = load_pattern("laplacian-2d-100");

    let perm = compute_ordering(&pattern, OrderingMethod::METIS)
        .unwrap();

    // Permutation should be valid
    assert!(perm.is_valid());
    assert_eq!(perm.len(), pattern.n);

    // Should reduce fill-in
    let natural_nnz = predict_fill(&pattern, &Permutation::identity(pattern.n));
    let metis_nnz = predict_fill(&pattern, &perm);

    assert!(metis_nnz <= natural_nnz);
}

#[test]
fn test_ordering_quality() {
    for case in test_cases_by_difficulty(Difficulty::Hard) {
        let pattern = case.matrix.pattern();

        let natural = Permutation::identity(pattern.n);
        let metis = compute_ordering(pattern, OrderingMethod::METIS)
            .unwrap();

        let natural_analysis = SymbolicAnalysis::analyze_with_ordering(
            pattern, &natural
        );
        let metis_analysis = SymbolicAnalysis::analyze_with_ordering(
            pattern, &metis
        );

        // METIS should not increase fill-in
        assert!(metis_analysis.predicted_nnz
               <= natural_analysis.predicted_nnz * 1.1);
    }
}
```

**Success Criteria:**
- [ ] METIS ordering reduces fill-in on test matrices
- [ ] Integration with metis-rs crate works
- [ ] Ordering quality comparable to SPRAL

#### 4.2: MC64 Matching Algorithm
**Task:** Implement maximum matching for indefinite matrices

**Algorithm Reference:**
- Duff & Koster (2001) - "On algorithms for permuting large entries to the diagonal"

**Implementation:**

```rust
pub struct MatchingResult {
    pub permutation: Permutation,
    pub scaling: Vec<f64>,
    pub matched: usize,  // Number of matched entries
}

pub fn mc64_matching(
    matrix: &SparseColMat<f64>,
    job: MC64Job
) -> Result<MatchingResult>;

pub enum MC64Job {
    MaximumProduct,      // Maximize product of diagonal
    MaximumSum,          // Maximize sum (for scaling)
    MaximumMin,          // Maximize minimum diagonal
}

// Hungarian algorithm implementation
fn hungarian_algorithm(
    cost_matrix: &[Vec<f64>]
) -> Vec<usize>;

// Auction algorithm (faster alternative)
pub fn auction_scaling(
    matrix: &SparseColMat<f64>,
    params: AuctionParams
) -> Result<MatchingResult>;

pub struct AuctionParams {
    pub epsilon: f64,
    pub max_iterations: usize,
}
```

**Testing:**

```rust
#[test]
fn test_mc64_basic() {
    // Matrix with small diagonal elements
    let matrix = create_badly_scaled_matrix();

    let result = mc64_matching(&matrix, MC64Job::MaximumProduct)
        .unwrap();

    // Should find full matching
    assert_eq!(result.matched, matrix.nrows());

    // Diagonal should be larger after permutation + scaling
    let scaled = matrix
        .permute_symmetric(&result.permutation)
        .scale_symmetric(&result.scaling);

    let diag_before = matrix.diagonal_max();
    let diag_after = scaled.diagonal_max();

    assert!(diag_after >= diag_before);
}

#[test]
fn test_auction_vs_mc64() {
    for case in test_cases_indefinite() {
        let mc64_result = mc64_matching(
            &case.matrix,
            MC64Job::MaximumProduct
        ).unwrap();

        let auction_result = auction_scaling(
            &case.matrix,
            AuctionParams::default()
        ).unwrap();

        // Both should find full matching
        assert_eq!(mc64_result.matched, case.matrix.nrows());
        assert_eq!(auction_result.matched, case.matrix.nrows());

        // Quality should be similar (within 10%)
        let mc64_quality = measure_scaling_quality(
            &case.matrix, &mc64_result
        );
        let auction_quality = measure_scaling_quality(
            &case.matrix, &auction_result
        );

        assert!((mc64_quality - auction_quality).abs() / mc64_quality < 0.1);
    }
}
```

**Success Criteria:**
- [ ] MC64 produces valid matching
- [ ] Improves diagonal dominance
- [ ] Auction algorithm is faster with similar quality

#### 4.3: Combined Ordering Strategy
**Task:** Implement matching-based ordering for indefinite systems

**Implementation:**

```rust
pub fn matching_based_ordering(
    matrix: &SparseColMat<f64>,
    params: MatchingOrderingParams
) -> Result<(Permutation, Vec<f64>)> {
    // 1. Compute matching to identify large off-diagonal entries
    let matching = mc64_matching(matrix, MC64Job::MaximumProduct)?;

    // 2. Apply matching permutation
    let matched_matrix = matrix.permute_symmetric(&matching.permutation);

    // 3. Constrained METIS ordering that keeps matched entries near diagonal
    let ordering = constrained_metis_ordering(
        &matched_matrix,
        &matching,
        params
    )?;

    // 4. Combine permutations
    let combined_perm = matching.permutation.compose(&ordering);

    Ok((combined_perm, matching.scaling))
}

pub struct MatchingOrderingParams {
    pub matching_algorithm: MatchingAlgorithm,
    pub ordering_method: OrderingMethod,
    pub scaling_method: ScalingMethod,
}

pub enum MatchingAlgorithm {
    MC64,
    Auction,
}

pub enum ScalingMethod {
    None,
    MC64,
    Auction,
    EquilibrationIterative,
}
```

**Testing:**

```rust
#[test]
fn test_matching_ordering_on_hard_indefinite() {
    for case in test_cases_by_difficulty(Difficulty::Hard) {
        let (perm, scaling) = matching_based_ordering(
            &case.matrix,
            MatchingOrderingParams::default()
        ).unwrap();

        // Analyze with this ordering
        let reordered = case.matrix
            .permute_symmetric(&perm)
            .scale_symmetric(&scaling);

        let analysis = SymbolicAnalysis::analyze(
            reordered.pattern(),
            AnalysisOptions::default()
        );

        // Should have reasonable fill-in
        if let Some(ref_results) = case.reference {
            // Should be within 20% of SPRAL's fill-in
            let ratio = analysis.predicted_nnz as f64
                       / ref_results.predicted_nnz as f64;
            assert!(ratio < 1.2);
        }
    }
}
```

**Success Criteria:**
- [ ] Matching-based ordering works on hard indefinite problems
- [ ] Fill-in comparable to SPRAL
- [ ] Numerical stability improved

#### 4.4: Integration with Analysis
**Task:** Connect ordering/scaling to symbolic analysis

**Implementation:**

```rust
pub struct AnalysisWithOrdering {
    pub analysis: SymbolicAnalysis,
    pub permutation: Permutation,
    pub scaling: Option<Vec<f64>>,
}

impl AnalysisWithOrdering {
    pub fn analyze_full(
        matrix: &SparseColMat<f64>,
        options: AnalysisOptions
    ) -> Result<Self> {
        // 1. Compute ordering and scaling
        let (perm, scaling) = match options.ordering_method {
            OrderingMethod::METIS => {
                let p = compute_ordering(matrix.pattern(), OrderingMethod::METIS)?;
                (p, None)
            }
            OrderingMethod::MatchingAMD => {
                let (p, s) = matching_based_ordering(matrix, options.matching_params)?;
                (p, Some(s))
            }
            // ... other methods
        };

        // 2. Apply permutation to pattern
        let reordered_pattern = matrix.pattern()
            .permute_symmetric(&perm);

        // 3. Perform symbolic analysis
        let analysis = SymbolicAnalysis::analyze(
            &reordered_pattern,
            options
        );

        Ok(Self {
            analysis,
            permutation: perm,
            scaling,
        })
    }
}
```

**Testing:**

```rust
#[test]
fn test_full_analysis_pipeline() {
    for case in all_test_cases() {
        let result = AnalysisWithOrdering::analyze_full(
            &case.matrix,
            AnalysisOptions {
                ordering_method: OrderingMethod::MatchingAMD,
                ..Default::default()
            }
        ).unwrap();

        // Should produce valid analysis
        assert!(result.analysis.predicted_nnz > 0);

        // Compare with reference
        if let Some(ref_results) = case.reference {
            compare_analysis(&result.analysis, &ref_results);
        }
    }
}
```

**Success Criteria:**
- [ ] Full analysis pipeline works end-to-end
- [ ] Ordering choices affect fill-in as expected
- [ ] Results match SPRAL's analysis

### Phase 4 Exit Criteria

**Required Outcomes:**
1. METIS ordering integrated and tested
2. MC64 matching produces quality results
3. Matching-based ordering works on hard indefinite matrices
4. Full analysis pipeline validated on test suite
5. Ordering quality comparable to SPRAL

**Validation Questions:**
- Does ordering reduce fill-in significantly?
- Does matching improve stability on hard indefinite problems?
- Are we making good ordering choices automatically?

**Checkpoint:** Run full analysis (with ordering/scaling) on entire test suite. Compare fill-in predictions with SPRAL. Should be within 20% on hard indefinite problems.

---

## Phase 5: Dense APTP Factorization

### Objectives
Implement A Posteriori Threshold Pivoting for dense symmetric indefinite matrices. This is the core numerical kernel.

### Deliverables

#### 5.1: Basic APTP Implementation (Single-Level)
**Task:** Implement simplest version of APTP

**Algorithm Reference:**
- Hogg, Duff, Lopez (2020) - Section 3 "A posteriori threshold pivoting"
- Algorithm 1 in the paper

**Implementation:**

```rust
pub struct APTPFactorization {
    pub l: Mat<f64>,        // Unit lower triangular
    pub d: Mat<f64>,        // Block diagonal (1×1 and 2×2 blocks)
    pub pivot_sequence: Vec<Pivot>,
    pub num_delays: usize,
}

pub enum Pivot {
    OneByOne { index: usize, value: f64 },
    TwoByTwo { indices: [usize; 2], block: [[f64; 2]; 2] },
}

pub fn aptp_factor(
    a: MatRef<f64>,
    threshold: f64,
    options: APTPOptions
) -> Result<APTPFactorization>;

pub struct APTPOptions {
    pub block_size: usize,
    pub threshold: f64,
    pub fallback_strategy: FallbackStrategy,
}

pub enum FallbackStrategy {
    CompleteBlocking,   // Use complete pivoting on failed block
    TPP,                // Fall back to threshold partial pivoting
    Delayed,            // Delay to parent node
}

// Core APTP algorithm
fn aptp_block(
    a: MatMut<f64>,
    threshold: f64
) -> BlockResult {
    // 1. Attempt to factor block assuming diagonal pivots
    // 2. Check stability: |l_ij| < 1/threshold
    // 3. If any entry fails, mark as failed pivot
    // 4. Return factorization + list of failed pivots
}
```

**Testing:**

```rust
#[test]
fn test_aptp_simple() {
    // Small positive definite matrix (should succeed with no delays)
    let a = Mat::from_fn(4, 4, |i, j| {
        if i == j { 4.0 } else if i.abs_diff(j) == 1 { 1.0 } else { 0.0 }
    });

    let result = aptp_factor(a.as_ref(), 0.01, APTPOptions::default())
        .unwrap();

    // Should have no delays
    assert_eq!(result.num_delays, 0);

    // Check factorization: A = LDL'
    let reconstructed = reconstruct_ldlt(&result);
    assert_matrix_close(&a, &reconstructed, 1e-12);
}

#[test]
fn test_aptp_indefinite() {
    // Indefinite matrix (some negative eigenvalues)
    let a = create_indefinite_matrix(10);

    let result = aptp_factor(a.as_ref(), 0.01, APTPOptions::default())
        .unwrap();

    // Should successfully factor
    let reconstructed = reconstruct_ldlt(&result);
    assert_matrix_close(&a, &reconstructed, 1e-10);
}

#[test]
fn test_stability_bounds() {
    for _ in 0..100 {
        let a = random_symmetric_matrix(20);
        let threshold = 0.01;

        let result = aptp_factor(a.as_ref(), threshold, APTPOptions::default())
            .unwrap();

        // Check stability: all |l_ij| < 1/threshold
        let max_l_entry = result.l.as_ref()
            .iter()
            .map(|&x| x.abs())
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();

        assert!(max_l_entry < 1.0 / threshold + 1e-10);
    }
}
```

**Success Criteria:**
- [ ] Factors positive definite matrices with no delays
- [ ] Handles indefinite matrices correctly
- [ ] Stability bounds enforced
- [ ] Reconstructed matrix matches input (within tolerance)

#### 5.2: Fallback Pivoting Strategies
**Task:** Implement strategies for failed pivots

**Implementation:**

```rust
// Complete pivoting for small failed blocks
fn complete_pivoting_fallback(
    block: MatMut<f64>
) -> FallbackResult {
    let (i, j, max_val) = find_absolute_maximum(block.as_ref());

    if i == j {
        // Diagonal pivot
        pivot_1x1(block, i)
    } else {
        // Off-diagonal pivot -> use 2×2
        pivot_2x2(block, i, j)
    }
}

// TPP fallback
fn tpp_fallback(
    block: MatMut<f64>,
    threshold: f64
) -> FallbackResult {
    // Standard threshold partial pivoting
    // Similar to MA57 algorithm
}

// Delayed pivoting
fn delay_pivots(
    block: MatRef<f64>,
    failed_indices: &[usize]
) -> DelayResult {
    // Mark these variables for elimination at parent node
    // Return contribution block
}
```

**Testing:**

```rust
#[test]
fn test_complete_pivoting_fallback() {
    // Create block designed to fail APTP
    let block = create_aptp_failure_case();

    let result = complete_pivoting_fallback(block.as_mut());

    // Should still factor correctly
    assert!(result.is_stable());
}

#[test]
fn test_fallback_strategy_comparison() {
    for case in challenging_matrices() {
        let cp_result = aptp_factor(
            case.as_ref(),
            0.01,
            APTPOptions {
                fallback_strategy: FallbackStrategy::CompleteBlocking,
                ..Default::default()
            }
        ).unwrap();

        let tpp_result = aptp_factor(
            case.as_ref(),
            0.01,
            APTPOptions {
                fallback_strategy: FallbackStrategy::TPP,
                ..Default::default()
            }
        ).unwrap();

        // Both should produce valid factorizations
        assert!(cp_result.is_valid());
        assert!(tpp_result.is_valid());

        // Complete pivoting may have fewer delays
        assert!(cp_result.num_delays <= tpp_result.num_delays);
    }
}
```

**Success Criteria:**
- [ ] All fallback strategies produce valid factorizations
- [ ] Complete pivoting handles worst-case blocks
- [ ] TPP fallback matches traditional behavior

#### 5.3: Two-Level APTP (Recursive)
**Task:** Implement nested APTP structure from SPRAL

**Algorithm Reference:**
- Hogg, Duff, Lopez (2020) - Section 4.2 "Two-level APTP"

**Implementation:**

```rust
pub struct TwoLevelAPTPOptions {
    pub outer_block_size: usize,     // e.g., 256
    pub inner_block_size: usize,     // e.g., 32
    pub outer_threshold: f64,
    pub inner_threshold: f64,
}

pub fn two_level_aptp_factor(
    a: MatRef<f64>,
    options: TwoLevelAPTPOptions
) -> Result<APTPFactorization> {
    let n = a.nrows();
    let nb_outer = options.outer_block_size;
    let nb_inner = options.inner_block_size;

    // Outer loop: process nb_outer columns at a time
    for outer_block in (0..n).step_by(nb_outer) {
        // Recursive APTP on this outer block
        let outer_result = aptp_block_recursive(
            a.submatrix(/* ... */),
            nb_inner,
            options.outer_threshold
        )?;

        // Handle any failed pivots
        if !outer_result.failed_pivots.is_empty() {
            handle_failed_outer_pivots(/* ... */)?;
        }

        // Update trailing submatrix
        update_schur_complement(/* ... */);
    }

    Ok(/* ... */)
}

fn aptp_block_recursive(
    block: MatMut<f64>,
    inner_block_size: usize,
    threshold: f64
) -> BlockResult {
    // Inner loop: process inner_block_size columns at a time
    // Use complete pivoting for truly small blocks
    // ...
}
```

**Testing:**

```rust
#[test]
fn test_two_level_vs_single_level() {
    let a = random_symmetric_matrix(256);

    let single_result = aptp_factor(
        a.as_ref(),
        0.01,
        APTPOptions::default()
    ).unwrap();

    let two_level_result = two_level_aptp_factor(
        a.as_ref(),
        TwoLevelAPTPOptions {
            outer_block_size: 64,
            inner_block_size: 16,
            outer_threshold: 0.01,
            inner_threshold: 0.01,
        }
    ).unwrap();

    // Both should produce equivalent results
    assert_eq!(single_result.num_delays, two_level_result.num_delays);

    // Factorizations should be equally accurate
    let recon_single = reconstruct_ldlt(&single_result);
    let recon_two_level = reconstruct_ldlt(&two_level_result);

    assert_matrix_close(&recon_single, &recon_two_level, 1e-12);
}
```

**Success Criteria:**
- [ ] Two-level APTP produces same results as single-level
- [ ] Performance improves for large blocks
- [ ] Memory access patterns improved

#### 5.4: Integration with faer
**Task:** Use faer for dense operations

**Implementation:**

```rust
use faer::{MatMut, MatRef, Parallelism};
use faer::linalg::{matmul, triangular_solve};

// Use faer's GEMM for Schur complement updates
fn update_schur_with_faer(
    schur: MatMut<f64>,
    l21: MatRef<f64>,
    d: MatRef<f64>,
    parallelism: Parallelism
) {
    // schur -= L21 * D * L21'
    let temp = matmul::matmul(
        Mat::zeros(l21.nrows(), l21.ncols()),
        l21,
        d,
        None,
        1.0,
        parallelism
    );

    matmul::matmul(
        schur,
        temp.as_ref(),
        l21.transpose(),
        Some(1.0),
        -1.0,
        parallelism
    );
}

// Use faer's TRSM for triangular solves
fn solve_with_l_faer(
    l: MatRef<f64>,
    rhs: MatMut<f64>,
    parallelism: Parallelism
) {
    triangular_solve::solve_lower_triangular(
        l,
        rhs,
        parallelism
    );
}
```

**Testing:**

```rust
#[test]
fn test_faer_integration() {
    let a = random_symmetric_matrix(100);

    // Factor using APTP with faer operations
    let result = aptp_factor_with_faer(
        a.as_ref(),
        0.01,
        Parallelism::Rayon(4)
    ).unwrap();

    // Should produce identical results to non-faer version
    let result_no_faer = aptp_factor(
        a.as_ref(),
        0.01,
        APTPOptions::default()
    ).unwrap();

    assert_matrix_close(
        &reconstruct_ldlt(&result),
        &reconstruct_ldlt(&result_no_faer),
        1e-14
    );
}

#[test]
fn test_faer_performance() {
    let sizes = [100, 500, 1000];

    for n in sizes {
        let a = random_symmetric_matrix(n);

        let start = Instant::now();
        let _ = aptp_factor_with_faer(
            a.as_ref(),
            0.01,
            Parallelism::Rayon(4)
        );
        let faer_time = start.elapsed();

        let start = Instant::now();
        let _ = aptp_factor_naive(a.as_ref(), 0.01);
        let naive_time = start.elapsed();

        println!("n={}: faer={:.2}ms, naive={:.2}ms, speedup={:.2}x",
                 n, faer_time.as_millis(), naive_time.as_millis(),
                 naive_time.as_secs_f64() / faer_time.as_secs_f64());
    }
}
```

**Success Criteria:**
- [ ] faer integration produces correct results
- [ ] Performance improved over naive implementation
- [ ] Parallelism through faer works correctly

#### 5.5: Validation on Hard Indefinite Problems
**Task:** Stress test APTP on challenging matrices

**Critical Test:**

```rust
#[test]
fn test_hard_indefinite_suite() {
    let hard_cases = test_cases_by_difficulty(Difficulty::Hard);

    for case in hard_cases {
        println!("Testing: {}", case.name);

        // Extract a large frontal matrix from the problem
        // (This simulates what would happen in multifrontal factorization)
        let front = extract_representative_front(&case.matrix, 500);

        let result = aptp_factor(
            front.as_ref(),
            0.01,
            APTPOptions::default()
        );

        match result {
            Ok(factorization) => {
                // Check residual
                let residual = compute_residual(&front, &factorization);

                println!("  SUCCESS: residual={:.2e}, delays={}",
                         residual, factorization.num_delays);

                assert!(residual < 1e-6,
                        "Poor residual on {}: {:.2e}",
                        case.name, residual);
            }
            Err(e) => {
                panic!("FAILED to factor {}: {:?}", case.name, e);
            }
        }
    }
}

#[test]
fn test_against_spral_dense() {
    // For matrices where we have SPRAL dense factorizations
    for case in get_spral_dense_references() {
        let result = aptp_factor(
            case.matrix.as_ref(),
            case.threshold,
            APTPOptions::default()
        ).unwrap();

        // Number of delays should match SPRAL (or be better)
        assert!(result.num_delays <= case.spral_delays);

        // Residual should be similar
        let our_residual = compute_residual(&case.matrix, &result);
        assert!(our_residual <= case.spral_residual * 10.0);
    }
}
```

**Success Criteria:**
- [ ] Successfully factors all hard indefinite matrices
- [ ] Residuals acceptable (< 1e-6)
- [ ] Number of delays reasonable
- [ ] Performance within 2× of SPRAL on dense problems

### Phase 5 Exit Criteria

**Required Outcomes:**
1. APTP algorithm implemented and validated
2. All fallback strategies working
3. Two-level structure improves performance
4. faer integration provides speedup
5. Hard indefinite problems factor correctly

**Validation Questions:**
- Does APTP maintain stability on all test cases?
- Are delayed pivots handled correctly?
- Is performance competitive with SPRAL's dense factorization?

**Checkpoint: HIGH RISK VALIDATION**

This is the critical checkpoint. Run comprehensive tests:

```rust
#[test]
#[ignore]  // Run manually with: cargo test critical_checkpoint -- --ignored
fn critical_checkpoint_aptp_validation() {
    println!("\n========================================");
    println!("CRITICAL CHECKPOINT: APTP VALIDATION");
    println!("========================================\n");

    let mut report = ValidationReport::new();

    // Test 1: All hard indefinite matrices
    for case in test_cases_by_difficulty(Difficulty::Hard) {
        let result = test_aptp_on_case(&case);
        report.add_result(&case.name, result);
    }

    // Test 2: Random matrices
    for i in 0..100 {
        let matrix = random_symmetric_indefinite(100);
        let result = test_aptp_stability(&matrix);
        report.add_result(&format!("random_{}", i), result);
    }

    // Test 3: Pathological cases
    let pathological = vec![
        create_near_singular_matrix(),
        create_highly_indefinite_matrix(),
        create_clustered_eigenvalue_matrix(),
    ];
    for (i, matrix) in pathological.iter().enumerate() {
        let result = test_aptp_on_matrix(matrix);
        report.add_result(&format!("pathological_{}", i), result);
    }

    // Generate report
    report.print();

    // Pass/fail criteria
    assert!(report.success_rate() > 0.95,
            "APTP failed on too many cases: {:.1}%",
            100.0 - report.success_rate() * 100.0);

    println!("\n✓ CRITICAL CHECKPOINT PASSED");
}
```

**If this checkpoint fails:** Investigate deeply before proceeding. APTP is the foundation of the entire solver. Options:
1. Tune thresholds and parameters
2. Improve fallback strategies
3. Consider alternative pivoting strategies
4. Consult with numerical analysis experts

**If it passes:** You have a working core factorization kernel. The rest is "just" assembly and orchestration.

---

## Phase 6: Multifrontal Assembly

### Objectives
Implement the multifrontal method: assemble frontal matrices, perform dense factorizations, propagate updates up the assembly tree.

### Deliverables

#### 6.1: Frontal Matrix Structure
**Task:** Define frontal matrix data structure

**Implementation:**

```rust
pub struct FrontalMatrix {
    pub node_id: usize,
    pub supernode: Supernode,

    // Frontal matrix structure: [F11  F12]
    //                           [F21  F22]
    // where F11 is the fully summed part
    pub nrows: usize,
    pub ncols: usize,
    pub fully_summed_cols: usize,

    // Dense storage
    pub data: Mat<f64>,

    // Row indices (mapping to global indices)
    pub row_indices: Vec<usize>,
}

impl FrontalMatrix {
    pub fn new(supernode: Supernode, row_indices: Vec<usize>) -> Self;

    pub fn f11(&self) -> MatRef<f64>;
    pub fn f11_mut(&mut self) -> MatMut<f64>;

    pub fn f21(&self) -> MatRef<f64>;
    pub fn f21_mut(&mut self) -> MatMut<f64>;

    pub fn f22(&self) -> MatRef<f64>;
    pub fn f22_mut(&mut self) -> MatMut<f64>;

    pub fn contribution_block(&self) -> MatRef<f64> {
        self.f22()
    }
}

pub struct FactorStorage {
    // Store L and D factors for each supernode
    factors: Vec<SupernodeFactors>,
}

pub struct SupernodeFactors {
    pub l: Mat<f64>,        // Lower triangular
    pub d: Vec<f64>,        // Diagonal (1×1 and 2×2 blocks)
    pub pivot_sequence: Vec<Pivot>,
}
```

**Testing:**

```rust
#[test]
fn test_frontal_matrix_structure() {
    let supernode = Supernode {
        first_col: 0,
        last_col: 3,
        row_structure: vec![0, 1, 2, 3, 5, 7],
    };

    let front = FrontalMatrix::new(
        supernode,
        vec![0, 1, 2, 3, 5, 7]
    );

    assert_eq!(front.nrows, 6);
    assert_eq!(front.ncols, 4);
    assert_eq!(front.fully_summed_cols, 4);

    // Check submatrix views
    assert_eq!(front.f11().nrows(), 4);
    assert_eq!(front.f11().ncols(), 4);
}
```

**Success Criteria:**
- [ ] Frontal matrix structure correct
- [ ] Submatrix views work properly
- [ ] Memory layout efficient

#### 6.2: Assembly Process
**Task:** Assemble contributions from children into parent front

**Implementation:**

```rust
pub fn assemble_front(
    front: &mut FrontalMatrix,
    child_contributions: &[ContributionBlock],
    original_matrix: &SparseColMat<f64>
) {
    // 1. Initialize front with original matrix entries
    initialize_from_sparse(front, original_matrix);

    // 2. Add contributions from child nodes
    for contrib in child_contributions {
        add_contribution(front, contrib);
    }
}

pub struct ContributionBlock {
    pub row_indices: Vec<usize>,
    pub data: Mat<f64>,
}

fn initialize_from_sparse(
    front: &mut FrontalMatrix,
    matrix: &SparseColMat<f64>
) {
    // Extract entries of A corresponding to this front
    for col in front.supernode.columns() {
        for (row, val) in matrix.column(col) {
            if let Some(local_row) = front.global_to_local_row(row) {
                front.data[(local_row, col)] = val;
            }
        }
    }
}

fn add_contribution(
    front: &mut FrontalMatrix,
    contrib: &ContributionBlock
) {
    // Map contribution's row indices to front's row indices
    // Add contribution block to appropriate part of front

    for (i, &global_i) in contrib.row_indices.iter().enumerate() {
        if let Some(local_i) = front.global_to_local_row(global_i) {
            for (j, &global_j) in contrib.row_indices.iter().enumerate() {
                if let Some(local_j) = front.global_to_local_row(global_j) {
                    front.data[(local_i, local_j)] += contrib.data[(i, j)];
                }
            }
        }
    }
}
```

**Testing:**

```rust
#[test]
fn test_assembly_simple() {
    // Create simple 2-level assembly tree
    let child_front = create_child_front();
    let child_factors = factor_front(&child_front);
    let child_contrib = compute_contribution(&child_factors);

    let mut parent_front = create_parent_front();
    assemble_front(
        &mut parent_front,
        &[child_contrib],
        &original_matrix
    );

    // Parent should contain original entries + child contribution
    verify_assembly(&parent_front, &original_matrix, &child_contrib);
}

#[test]
fn test_assembly_commutativity() {
    // Assembly order shouldn't matter
    let contribs = vec![contrib1, contrib2, contrib3];

    let mut front1 = create_front();
    for c in &contribs {
        add_contribution(&mut front1, c);
    }

    let mut front2 = create_front();
    for c in contribs.iter().rev() {
        add_contribution(&mut front2, c);
    }

    assert_matrix_close(&front1.data, &front2.data, 1e-14);
}
```

**Success Criteria:**
- [ ] Assembly process correct
- [ ] Contributions added properly
- [ ] Order independence verified

#### 6.3: Front Factorization
**Task:** Factor each frontal matrix using APTP

**Implementation:**

```rust
pub fn factor_front(
    front: &FrontalMatrix,
    options: APTPOptions
) -> Result<FrontFactorization> {
    // 1. Factor F11 (fully summed part)
    let f11_factors = aptp_factor(
        front.f11(),
        options.threshold,
        options
    )?;

    // 2. Solve for F21: F21 = F21 * inv(L11)
    let mut f21 = front.f21().to_owned();
    solve_with_l(f11_factors.l.as_ref(), f21.as_mut());

    // 3. Compute contribution block: F22 - F21 * D11 * F21'
    let contrib = compute_schur_complement(
        front.f22(),
        f21.as_ref(),
        f11_factors.d.as_ref()
    );

    Ok(FrontFactorization {
        l11: f11_factors.l,
        d11: f11_factors.d,
        l21: f21,
        contribution: contrib,
    })
}

pub struct FrontFactorization {
    pub l11: Mat<f64>,
    pub d11: Mat<f64>,
    pub l21: Mat<f64>,
    pub contribution: ContributionBlock,
}
```

**Testing:**

```rust
#[test]
fn test_front_factorization() {
    let front = create_test_front(10, 10);

    let factors = factor_front(&front, APTPOptions::default())
        .unwrap();

    // Verify: F11 = L11 * D11 * L11'
    let f11_reconstructed = reconstruct_ldlt_from_front(&factors);
    assert_matrix_close(&f11_reconstructed, front.f11(), 1e-12);

    // Verify contribution block
    let expected_contrib = compute_expected_contribution(&front, &factors);
    assert_matrix_close(
        &factors.contribution.data,
        &expected_contrib,
        1e-12
    );
}
```

**Success Criteria:**
- [ ] Front factorization correct
- [ ] Contribution blocks computed properly
- [ ] Numerical accuracy maintained

#### 6.4: Tree Traversal
**Task:** Implement bottom-up tree traversal

**Implementation:**

```rust
pub fn factor_tree(
    assembly_tree: &AssemblyTree,
    matrix: &SparseColMat<f64>,
    options: FactorOptions
) -> Result<TreeFactorization> {
    let mut factors = FactorStorage::new();

    // Postorder traversal (children before parents)
    for node_id in assembly_tree.postorder() {
        let node = assembly_tree.node(node_id);

        // 1. Create frontal matrix for this node
        let mut front = FrontalMatrix::new(
            node.supernode.clone(),
            node.row_indices.clone()
        );

        // 2. Gather contributions from children
        let child_contribs: Vec<_> = node.children
            .iter()
            .map(|&child_id| factors.get_contribution(child_id))
            .collect();

        // 3. Assemble front
        assemble_front(&mut front, &child_contribs, matrix);

        // 4. Factor front
        let front_factors = factor_front(&front, options.aptp_options)?;

        // 5. Store factors
        factors.store(node_id, front_factors);
    }

    Ok(TreeFactorization {
        storage: factors,
        tree: assembly_tree.clone(),
    })
}
```

**Testing:**

```rust
#[test]
fn test_tree_factorization_small() {
    let matrix = load_test_matrix("laplacian-5x5");
    let analysis = AnalysisWithOrdering::analyze_full(
        &matrix,
        AnalysisOptions::default()
    ).unwrap();

    let factorization = factor_tree(
        &analysis.analysis.assembly_tree,
        &matrix,
        FactorOptions::default()
    ).unwrap();

    // Should have factors for every node
    assert_eq!(
        factorization.storage.len(),
        analysis.analysis.assembly_tree.num_nodes()
    );
}

#[test]
fn test_against_dense_factorization() {
    // For small matrices, compare with dense APTP
    let matrix = create_test_matrix(20, 20);

    let dense_factors = aptp_factor(
        matrix.to_dense().as_ref(),
        0.01,
        APTPOptions::default()
    ).unwrap();

    let analysis = analyze(&matrix);
    let tree_factors = factor_tree(
        &analysis.assembly_tree,
        &matrix,
        FactorOptions::default()
    ).unwrap();

    // Should produce equivalent factorizations
    let dense_reconstructed = reconstruct_from_dense(&dense_factors);
    let tree_reconstructed = reconstruct_from_tree(&tree_factors);

    assert_matrix_close(&dense_reconstructed, &tree_reconstructed, 1e-10);
}
```

**Success Criteria:**
- [ ] Tree traversal in correct order
- [ ] All nodes factored successfully
- [ ] Equivalent to dense factorization (for small problems)

### Phase 6 Exit Criteria

**Required Outcomes:**
1. Multifrontal assembly implemented
2. Tree traversal working correctly
3. Full factorization on small matrices
4. Numerical correctness validated

**Validation Questions:**
- Does factorization work on entire test suite?
- Are contribution blocks computed correctly?
- Does tree factorization match dense factorization?

**Checkpoint:** Factor all "easy indefinite" test matrices. Verify residuals < 1e-10.

---

## Phase 7: Solve Phase & Integration

### Objectives
Implement forward/backward substitution, integrate all components, create user-facing API.

### Deliverables

#### 7.1: Triangular Solve
**Task:** Implement solve with L and D

**Implementation:**

```rust
pub fn solve_ldlt(
    factorization: &TreeFactorization,
    rhs: &[f64]
) -> Vec<f64> {
    let mut x = rhs.to_vec();

    // 1. Permute RHS
    apply_permutation(&mut x, &factorization.permutation);

    // 2. Scale if needed
    if let Some(scaling) = &factorization.scaling {
        apply_scaling(&mut x, scaling);
    }

    // 3. Forward solve: L * y = Pb
    forward_solve_tree(&factorization, &mut x);

    // 4. Diagonal solve: D * z = y
    diagonal_solve(&factorization, &mut x);

    // 5. Backward solve: L' * x = z
    backward_solve_tree(&factorization, &mut x);

    // 6. Unscale
    if let Some(scaling) = &factorization.scaling {
        apply_scaling(&mut x, scaling);  // S is its own inverse
    }

    // 7. Unpermute
    apply_inverse_permutation(&mut x, &factorization.permutation);

    x
}

fn forward_solve_tree(
    factorization: &TreeFactorization,
    x: &mut [f64]
) {
    // Postorder traversal
    for node_id in factorization.tree.postorder() {
        let node_factors = factorization.storage.get(node_id);

        // Solve with this node's L factor
        forward_solve_supernode(node_factors, x);
    }
}

fn backward_solve_tree(
    factorization: &TreeFactorization,
    x: &mut [f64]
) {
    // Reverse postorder traversal
    for node_id in factorization.tree.postorder().iter().rev() {
        let node_factors = factorization.storage.get(node_id);

        // Solve with this node's L' factor
        backward_solve_supernode(node_factors, x);
    }
}
```

**Testing:**

```rust
#[test]
fn test_solve_correctness() {
    for case in all_test_cases() {
        let analysis = analyze_full(&case.matrix);
        let factorization = factor_tree(&analysis, &case.matrix).unwrap();

        // Create random RHS
        let b = random_vector(case.matrix.nrows());

        // Solve
        let x = solve_ldlt(&factorization, &b);

        // Check residual: ||Ax - b|| / ||b||
        let ax = case.matrix.mul_vec(&x);
        let residual = compute_residual(&ax, &b);

        println!("{}: residual = {:.2e}", case.name, residual);

        assert!(residual < 1e-10,
                "Poor solve accuracy on {}: {:.2e}",
                case.name, residual);
    }
}

#[test]
fn test_solve_multiple_rhs() {
    let matrix = load_test_matrix("medium-matrix");
    let analysis = analyze_full(&matrix);
    let factorization = factor_tree(&analysis, &matrix).unwrap();

    let nrhs = 10;
    let b_mat = random_matrix(matrix.nrows(), nrhs);

    // Solve each RHS
    let mut x_mat = Mat::zeros(matrix.nrows(), nrhs);
    for j in 0..nrhs {
        let b_j = b_mat.col(j).try_as_slice().unwrap();
        let x_j = solve_ldlt(&factorization, b_j);
        x_mat.col_mut(j).copy_from_slice(&x_j);
    }

    // Check all solutions
    for j in 0..nrhs {
        let b_j = b_mat.col(j).try_as_slice().unwrap();
        let x_j = x_mat.col(j).try_as_slice().unwrap();
        let ax = matrix.mul_vec(x_j);
        let residual = compute_residual(&ax, b_j);

        assert!(residual < 1e-10);
    }
}
```

**Success Criteria:**
- [ ] Solve produces correct solutions
- [ ] Residuals < 1e-10 on all test matrices
- [ ] Multiple RHS handled efficiently

#### 7.2: User-Facing API
**Task:** Create ergonomic public API

**Implementation:**

```rust
// High-level solver interface
pub struct SparseLDLT {
    analysis: AnalysisWithOrdering,
    factorization: Option<TreeFactorization>,
}

impl SparseLDLT {
    /// Create solver from sparse matrix
    pub fn new(matrix: SparseColMat<f64>) -> Result<Self> {
        Self::with_options(matrix, SolverOptions::default())
    }

    pub fn with_options(
        matrix: SparseColMat<f64>,
        options: SolverOptions
    ) -> Result<Self> {
        // Analyze
        let analysis = AnalysisWithOrdering::analyze_full(
            &matrix,
            options.analysis_options
        )?;

        Ok(Self {
            analysis,
            factorization: None,
        })
    }

    /// Perform numerical factorization
    pub fn factor(&mut self, matrix: &SparseColMat<f64>) -> Result<()> {
        let factorization = factor_tree(
            &self.analysis.analysis.assembly_tree,
            matrix,
            FactorOptions::default()
        )?;

        self.factorization = Some(factorization);
        Ok(())
    }

    /// Solve Ax = b
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>> {
        let factorization = self.factorization.as_ref()
            .ok_or(Error::NotFactored)?;

        Ok(solve_ldlt(factorization, b))
    }

    /// One-shot solve: analyze + factor + solve
    pub fn solve_matrix(
        matrix: SparseColMat<f64>,
        b: &[f64]
    ) -> Result<Vec<f64>> {
        let mut solver = Self::new(matrix.clone())?;
        solver.factor(&matrix)?;
        solver.solve(b)
    }

    /// Get inertia (num positive, num negative, num zero eigenvalues)
    pub fn inertia(&self) -> Option<(usize, usize, usize)> {
        self.factorization.as_ref().map(|f| f.inertia())
    }

    /// Refactor with same sparsity pattern
    pub fn refactor(&mut self, matrix: &SparseColMat<f64>) -> Result<()> {
        // Check pattern matches
        if !same_pattern(&self.analysis.analysis.assembly_tree, matrix) {
            return Err(Error::PatternMismatch);
        }

        self.factor(matrix)
    }
}

pub struct SolverOptions {
    pub ordering: OrderingMethod,
    pub scaling: ScalingMethod,
    pub threshold: f64,
    pub analysis_options: AnalysisOptions,
    pub factor_options: FactorOptions,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            ordering: OrderingMethod::MatchingAMD,
            scaling: ScalingMethod::MC64,
            threshold: 0.01,
            analysis_options: AnalysisOptions::default(),
            factor_options: FactorOptions::default(),
        }
    }
}
```

**Testing:**

```rust
#[test]
fn test_api_basic() {
    let matrix = load_test_matrix("small-indefinite");
    let b = vec![1.0; matrix.nrows()];

    // One-shot solve
    let x = SparseLDLT::solve_matrix(matrix, &b).unwrap();

    // Verify
    let residual = compute_residual_sparse(&matrix, &x, &b);
    assert!(residual < 1e-10);
}

#[test]
fn test_api_refactor() {
    let matrix1 = load_test_matrix("truss-1");
    let matrix2 = load_test_matrix("truss-2");  // Same pattern

    let mut solver = SparseLDLT::new(matrix1.clone()).unwrap();
    solver.factor(&matrix1).unwrap();

    let b = vec![1.0; matrix1.nrows()];
    let x1 = solver.solve(&b).unwrap();

    // Refactor with same pattern
    solver.refactor(&matrix2).unwrap();
    let x2 = solver.solve(&b).unwrap();

    // Both should be correct solutions
    assert!(compute_residual_sparse(&matrix1, &x1, &b) < 1e-10);
    assert!(compute_residual_sparse(&matrix2, &x2, &b) < 1e-10);
}

#[test]
fn test_api_error_handling() {
    let matrix = load_test_matrix("small-matrix");

    let solver = SparseLDLT::new(matrix.clone()).unwrap();

    // Should error if try to solve before factoring
    let b = vec![1.0; matrix.nrows()];
    assert!(solver.solve(&b).is_err());
}
```

**Success Criteria:**
- [ ] API is ergonomic and safe
- [ ] Error handling is clear
- [ ] Documentation is comprehensive
- [ ] Examples work correctly

#### 7.3: Integration Tests
**Task:** End-to-end testing on full test suite

**Implementation:**

```rust
#[test]
fn test_full_suite_sequential() {
    let mut report = TestReport::new();

    for case in all_test_cases() {
        println!("Testing: {}", case.name);

        let result = test_full_solve(&case);
        report.add(case.name.clone(), result);
    }

    report.print();

    // Success criteria
    assert!(report.success_rate() > 0.95);
    assert!(report.median_residual() < 1e-9);
}

fn test_full_solve(case: &SolverTestCase) -> TestResult {
    // 1. Solve
    let b = random_vector(case.matrix.nrows());

    let solver_result = SparseLDLT::solve_matrix(
        case.matrix.clone(),
        &b
    );

    let x = match solver_result {
        Ok(x) => x,
        Err(e) => return TestResult::Failed(format!("{:?}", e)),
    };

    // 2. Check residual
    let residual = compute_residual_sparse(&case.matrix, &x, &b);

    // 3. Check inertia if available
    let inertia_ok = if let Some(ref_inertia) = case.reference.as_ref()
        .and_then(|r| r.inertia)
    {
        let solver = SparseLDLT::new(case.matrix.clone()).unwrap();
        let computed_inertia = solver.inertia().unwrap();
        computed_inertia == ref_inertia
    } else {
        true
    };

    // 4. Compare with SPRAL if available
    let spral_comparison = if let Some(ref_results) = case.reference {
        compare_with_spral_results(&case, &ref_results)
    } else {
        None
    };

    TestResult::Success {
        residual,
        inertia_ok,
        spral_comparison,
    }
}
```

**Success Criteria:**
- [ ] >95% of test matrices solve successfully
- [ ] Median residual < 1e-9
- [ ] Inertia matches reference when available
- [ ] Performance within 2× of SPRAL (sequential)

### Phase 7 Exit Criteria

**Required Outcomes:**
1. Complete solver working end-to-end
2. All test matrices solve correctly
3. API is usable and documented
4. Integration tests pass

**Validation Questions:**
- Can users solve their problems with simple API?
- Are error messages clear and actionable?
- Is performance acceptable?

**Checkpoint:** Run full test suite and generate comprehensive report. Should match or exceed SPRAL's numerical quality.

---

## Phase 8: Parallelization

### Objectives
Add parallel execution using Rayon, following faer's parallelism patterns.

### Deliverables

#### 8.1: Parallel Tree Traversal
**Task:** Parallelize bottom-up tree factorization

**Implementation:**

```rust
pub fn factor_tree_parallel(
    assembly_tree: &AssemblyTree,
    matrix: &SparseColMat<f64>,
    options: FactorOptions,
    parallelism: Parallelism
) -> Result<TreeFactorization> {
    let mut factors = Arc::new(Mutex::new(FactorStorage::new()));

    match parallelism {
        Parallelism::None => {
            factor_tree_sequential(assembly_tree, matrix, options)
        }
        Parallelism::Rayon(n_threads) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n_threads)
                .build()
                .unwrap()
                .install(|| {
                    factor_tree_level_by_level(
                        assembly_tree,
                        matrix,
                        options,
                        factors
                    )
                })
        }
    }
}

fn factor_tree_level_by_level(
    assembly_tree: &AssemblyTree,
    matrix: &SparseColMat<f64>,
    options: FactorOptions,
    factors: Arc<Mutex<FactorStorage>>
) -> Result<TreeFactorization> {
    // Process tree level-by-level for parallelism
    let levels = assembly_tree.level_order();

    for level_nodes in levels {
        // All nodes at same level can be factored in parallel
        level_nodes.par_iter().try_for_each(|&node_id| {
            let node = assembly_tree.node(node_id);

            // Create and assemble front
            let mut front = create_and_assemble_front(
                node,
                matrix,
                &factors.lock().unwrap()
            );

            // Factor
            let front_factors = factor_front(&front, options.aptp_options)?;

            // Store
            factors.lock().unwrap().store(node_id, front_factors);

            Ok::<(), Error>(())
        })?;
    }

    Ok(TreeFactorization {
        storage: Arc::try_unwrap(factors).unwrap().into_inner().unwrap(),
        tree: assembly_tree.clone(),
    })
}
```

**Testing:**

```rust
#[test]
fn test_parallel_determinism() {
    let matrix = load_test_matrix("medium-matrix");
    let analysis = analyze_full(&matrix);

    // Run 10 times in parallel
    let results: Vec<_> = (0..10)
        .into_par_iter()
        .map(|_| {
            factor_tree_parallel(
                &analysis.assembly_tree,
                &matrix,
                FactorOptions::default(),
                Parallelism::Rayon(4)
            ).unwrap()
        })
        .collect();

    // All results should be identical
    for (i, result) in results.iter().enumerate() {
        assert_eq!(result.inertia(), results[0].inertia(),
                   "Run {} differs from run 0", i);

        // Solve and check residuals match
        let b = vec![1.0; matrix.nrows()];
        let x_i = solve_ldlt(result, &b);
        let x_0 = solve_ldlt(&results[0], &b);

        assert_vector_close(&x_i, &x_0, 1e-14);
    }
}

#[test]
fn test_parallel_vs_sequential() {
    for case in test_cases_medium() {
        let analysis = analyze_full(&case.matrix);

        let seq = factor_tree(
            &analysis.assembly_tree,
            &case.matrix,
            FactorOptions::default()
        ).unwrap();

        let par = factor_tree_parallel(
            &analysis.assembly_tree,
            &case.matrix,
            FactorOptions::default(),
            Parallelism::Rayon(4)
        ).unwrap();

        // Should produce identical results
        assert_eq!(seq.inertia(), par.inertia());

        // Test solve
        let b = vec![1.0; case.matrix.nrows()];
        let x_seq = solve_ldlt(&seq, &b);
        let x_par = solve_ldlt(&par, &b);

        assert_vector_close(&x_seq, &x_par, 1e-13);
    }
}
```

**Success Criteria:**
- [ ] Parallel execution is deterministic
- [ ] Results identical to sequential
- [ ] Speedup observed on large problems

#### 8.2: Load Balancing
**Task:** Implement work-stealing for unbalanced trees

**Implementation:**

```rust
pub struct LoadBalancingStrategy {
    pub method: BalancingMethod,
    pub min_work_per_thread: f64,  // flops
}

pub enum BalancingMethod {
    LevelByLevel,      // Simple: all nodes at same level
    WorkStealing,      // Dynamic: steal from other threads
    StaticPartition,   // Map subtrees to threads
}

fn static_partition_tree(
    assembly_tree: &AssemblyTree,
    n_threads: usize
) -> Vec<Vec<usize>> {
    // Partition tree into subtrees of roughly equal work
    let total_flops = assembly_tree.total_flops();
    let target_flops_per_thread = total_flops / n_threads as f64;

    // Greedy assignment of subtrees to threads
    let mut thread_assignments = vec![Vec::new(); n_threads];
    let mut thread_loads = vec![0.0; n_threads];

    for node_id in assembly_tree.postorder() {
        let node_flops = assembly_tree.node(node_id).flops;

        // Assign to thread with least work
        let min_thread = thread_loads
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;

        thread_assignments[min_thread].push(node_id);
        thread_loads[min_thread] += node_flops;
    }

    thread_assignments
}
```

**Testing:**

```rust
#[test]
fn test_load_balancing() {
    let unbalanced_matrix = create_unbalanced_tree_matrix();
    let analysis = analyze_full(&unbalanced_matrix);

    // Time with different strategies
    let strategies = vec![
        LoadBalancingStrategy {
            method: BalancingMethod::LevelByLevel,
            min_work_per_thread: 1e6,
        },
        LoadBalancingStrategy {
            method: BalancingMethod::StaticPartition,
            min_work_per_thread: 1e6,
        },
    ];

    for strategy in strategies {
        let start = Instant::now();
        let _ = factor_tree_parallel_balanced(
            &analysis.assembly_tree,
            &unbalanced_matrix,
            FactorOptions::default(),
            Parallelism::Rayon(4),
            strategy
        ).unwrap();
        let time = start.elapsed();

        println!("{:?}: {:.2}ms", strategy.method, time.as_millis());
    }
}
```

**Success Criteria:**
- [ ] Load balancing improves performance on unbalanced trees
- [ ] Thread utilization >80%
- [ ] Overhead acceptable

#### 8.3: Parallel Solve
**Task:** Parallelize forward/backward substitution

**Implementation:**

```rust
pub fn solve_ldlt_parallel(
    factorization: &TreeFactorization,
    rhs: &[f64],
    parallelism: Parallelism
) -> Vec<f64> {
    match parallelism {
        Parallelism::None => solve_ldlt(factorization, rhs),
        Parallelism::Rayon(n_threads) => {
            solve_ldlt_parallel_impl(factorization, rhs, n_threads)
        }
    }
}

fn solve_ldlt_parallel_impl(
    factorization: &TreeFactorization,
    rhs: &[f64],
    n_threads: usize
) -> Vec<f64> {
    // Forward/backward solve has limited parallelism
    // due to dependencies
    //
    // Strategy: parallelize within supernodes for large fronts
    // use level-set parallelism for small fronts

    let mut x = rhs.to_vec();

    // Permute/scale
    apply_permutation(&mut x, &factorization.permutation);
    if let Some(scaling) = &factorization.scaling {
        apply_scaling(&mut x, scaling);
    }

    // Forward solve with level-set parallelism
    for level in factorization.tree.level_order() {
        level.par_iter().for_each(|&node_id| {
            // Independent operations at same level
            forward_solve_node(factorization, node_id, &mut x);
        });
    }

    // Similar for backward solve...

    x
}
```

**Testing:**

```rust
#[test]
fn test_parallel_solve() {
    for case in test_cases_large() {
        let analysis = analyze_full(&case.matrix);
        let factorization = factor_tree_parallel(
            &analysis.assembly_tree,
            &case.matrix,
            FactorOptions::default(),
            Parallelism::Rayon(4)
        ).unwrap();

        let b = random_vector(case.matrix.nrows());

        let x_seq = solve_ldlt(&factorization, &b);
        let x_par = solve_ldlt_parallel(
            &factorization,
            &b,
            Parallelism::Rayon(4)
        );

        // Should produce same results
        assert_vector_close(&x_seq, &x_par, 1e-13);
    }
}
```

**Success Criteria:**
- [ ] Parallel solve correct
- [ ] Speedup on large problems
- [ ] Deterministic results

#### 8.4: Parallel Benchmarking
**Task:** Measure parallel performance

**Implementation:**

```rust
#[bench]
fn bench_parallel_scaling(b: &mut Bencher) {
    let matrix = load_test_matrix("large-3d-fem");
    let analysis = analyze_full(&matrix);

    let thread_counts = [1, 2, 4, 8];

    for &n_threads in &thread_counts {
        b.iter_custom(|iters| {
            let mut total_time = Duration::ZERO;

            for _ in 0..iters {
                let start = Instant::now();
                let _ = factor_tree_parallel(
                    &analysis.assembly_tree,
                    &matrix,
                    FactorOptions::default(),
                    Parallelism::Rayon(n_threads)
                );
                total_time += start.elapsed();
            }

            total_time
        });

        // Compute speedup
        let t1 = benchmark_time_one_thread();
        let tn = benchmark_time_n_threads(n_threads);
        let speedup = t1 / tn;

        println!("Threads: {}, Speedup: {:.2}x", n_threads, speedup);
    }
}

#[test]
fn test_parallel_scaling_report() {
    let mut report = ParallelScalingReport::new();

    for case in test_cases_large() {
        let scaling_data = measure_parallel_scaling(&case, &[1, 2, 4, 8]);
        report.add(case.name, scaling_data);
    }

    report.print();
    report.save_csv("parallel_scaling.csv");
}
```

**Success Criteria:**
- [ ] 3× speedup on 4 cores for large problems
- [ ] Efficient: >75% parallel efficiency
- [ ] No performance regressions

### Phase 8 Exit Criteria

**Required Outcomes:**
1. Parallel factorization working
2. Deterministic results
3. Measurable speedup
4. Load balancing effective

**Validation Questions:**
- Is parallelism adding value?
- Are there race conditions?
- Is the overhead acceptable?

**Checkpoint:** Run parallel scalability study on large test matrices. Generate scaling plots.

---

## Phase 9: Supernodal APTP Optimization

### Objectives
Implement supernodal (block-wise) APTP factorization for improved performance on large matrices. This phase adds the multifrontal infrastructure deferred from Phase 3.

### Motivation
Simplicial (column-by-column) APTP works but is inefficient for large dense blocks:
- Poor cache locality
- Limited use of BLAS-3 operations
- Lower parallelism potential

Supernodal factorization processes dense blocks together, achieving 2-5x speedup on large problems.

### Deliverables

#### 9.1: Supernode Detection (Adapt faer)
**Task:** Identify supernodes and build assembly tree

**Approach:**
Adapt faer's supernodal infrastructure for APTP with 2x2 pivots:
- Start with faer's `SymbolicSupernodalCholesky` as reference
- Modify to accommodate potential supernode splitting due to delayed pivots
- Use relaxed supernodes for better performance

**Implementation:**

```rust
use faer::sparse::linalg::cholesky::SymbolicSupernodalParams;

pub struct SupernodalAptpSymbolic {
    /// Elimination tree (from Phase 3)
    etree: Vec<isize>,

    /// Supernode structure
    supernode_begin: Vec<usize>,  // supernode_begin[s] = first col of supernode s
    supernode_end: Vec<usize>,    // supernode_end[s] = last col + 1

    /// Assembly tree (parent-child relationships)
    supernode_parent: Vec<Option<usize>>,

    /// Row structure for each supernode
    supernode_row_indices: Vec<usize>,
    supernode_row_ptrs: Vec<usize>,
}

impl SupernodalAptpSymbolic {
    pub fn from_simplicial(
        simplicial: &AptpSymbolic,
        params: SupernodalParams,
    ) -> Self {
        // Detect fundamental supernodes
        let fundamental = detect_fundamental_supernodes(
            &simplicial.etree,
            &simplicial.col_counts,
        );

        // Apply relaxation
        let relaxed = relax_supernodes(fundamental, params);

        // Build assembly tree
        Self::build_assembly_tree(relaxed, &simplicial.etree)
    }
}

pub struct SupernodalParams {
    pub relax_factor: f64,     // How much extra fill to tolerate
    pub min_size: usize,       // Minimum supernode size
}
```

**Testing:**

```rust
#[test]
fn test_supernode_detection() {
    for case in test_cases() {
        let simplicial = AptpSymbolic::analyze(case.matrix().symbolic()).unwrap();

        let supernodal = SupernodalAptpSymbolic::from_simplicial(
            &simplicial,
            SupernodalParams::default(),
        );

        // Should have fewer "nodes" than columns
        assert!(supernodal.num_supernodes() < simplicial.dimension());
    }
}
```

**Success Criteria:**
- [ ] Supernode detection working on all test matrices
- [ ] Matches faer's supernode counts on SPD matrices
- [ ] Assembly tree correctly constructed

**Time Estimate:** 1 week

#### 9.2: Supernodal APTP Factorization Kernel
**Task:** Implement dense APTP on frontal matrices

**Approach:**
Process supernodes as dense blocks, applying APTP within each supernode:
- Dense Bunch-Kaufman-style pivoting within supernode
- Handle supernode splitting when 2x2 pivots cross boundaries
- Update propagation to parent supernodes

**Implementation:**

```rust
pub struct SupernodalAptpNumeric {
    symbolic: SupernodalAptpSymbolic,
    supernode_factors: Vec<SupernodeFactor>,
    d_blocks: MixedDiagonal<f64>,
}

pub struct SupernodeFactor {
    supernode_id: usize,
    /// Dense L block (column-major)
    l_dense: Vec<f64>,
    nrows: usize,
    ncols: usize,
}

impl SupernodalAptpNumeric {
    pub fn factor(
        symbolic: SupernodalAptpSymbolic,
        matrix: AptpMatrixRef<'_, f64>,
        options: AptpOptions,
    ) -> Result<Self> {
        // Process supernodes in postorder
        for snode_id in symbolic.postorder() {
            // Assemble frontal matrix
            let mut frontal = assemble_frontal_matrix(snode_id, &symbolic, &matrix);

            // Apply APTP to dense frontal
            let (l_block, d_block, pivots) = factor_frontal_aptp(
                &mut frontal,
                options.threshold,
            )?;

            // Store factors
            // Update parent supernodes
        }

        Ok(Self { /* ... */ })
    }
}

/// Factor dense frontal matrix with APTP
fn factor_frontal_aptp(
    frontal: &mut DenseFrontalMatrix,
    threshold: f64,
) -> Result<(Vec<f64>, Vec<Block2x2<f64>>, Vec<PivotType>)> {
    // Dense APTP algorithm (from Phase 5) on frontal matrix
    todo!()
}
```

**Success Criteria:**
- [ ] Supernodal factorization produces correct factors
- [ ] Residuals match simplicial version
- [ ] Performance improvement on large matrices

**Time Estimate:** 2-3 weeks

#### 9.3: Supernodal Solve
**Task:** Implement triangular solves with supernodal structure

**Implementation:**

```rust
impl SupernodalAptpNumeric {
    pub fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        let mut x = rhs.to_vec();

        // Forward solve: L * y = P * rhs
        for snode_id in self.symbolic.postorder() {
            solve_supernode_forward(snode_id, &self.supernode_factors, &mut x);
        }

        // D solve
        self.d_blocks.solve(&mut x);

        // Backward solve: L^T * x = y
        for snode_id in self.symbolic.postorder().rev() {
            solve_supernode_backward(snode_id, &self.supernode_factors, &mut x);
        }

        x
    }
}
```

**Success Criteria:**
- [ ] Solve produces correct results (matches simplicial)
- [ ] Faster than simplicial solve

**Time Estimate:** 1 week

#### 9.4: Performance Comparison
**Task:** Benchmark simplicial vs. supernodal

**Deliverable:** Performance report showing:
- Factor time: simplicial vs. supernodal
- Solve time: simplicial vs. supernodal
- Memory usage comparison
- Identification of problem sizes where supernodal wins

**Success Criteria:**
- [ ] Supernodal 2-5x faster on large problems (n > 10,000)
- [ ] Simplicial still correct and available as fallback

**Time Estimate:** 3-4 days

### Phase 9 Exit Criteria

**Required Outcomes:**
1. Supernodal APTP working and tested
2. Performance improvement demonstrated
3. Both simplicial and supernodal available as solver options
4. Clear guidance on when to use each

**Validation Questions:**
- Does supernodal match simplicial accuracy?
- Is the speedup significant?
- Are there any correctness issues?

**Checkpoint:** Run full benchmark suite comparing simplicial vs. supernodal. Document performance characteristics.

**Time Estimate:** 4-6 weeks

---

## Phase 10: Optimization & Polish

### Objectives
Performance optimization, memory efficiency, API refinement, comprehensive documentation.

### Deliverables

#### 10.1: Performance Profiling
**Task:** Identify and optimize bottlenecks

**Process:**

```rust
// Use built-in profiling tools
#[test]
#[ignore]
fn profile_full_solve() {
    let cases = test_cases_by_difficulty(Difficulty::Hard);

    for case in cases {
        println!("\nProfiling: {}", case.name);

        let profiler = ProfileRecorder::new();

        // Analyze
        profiler.start_section("analyze");
        let analysis = analyze_full(&case.matrix);
        profiler.end_section();

        // Factor
        profiler.start_section("factor");
        let factorization = factor_tree(
            &analysis.assembly_tree,
            &case.matrix,
            FactorOptions::default()
        ).unwrap();
        profiler.end_section();

        // Solve
        profiler.start_section("solve");
        let b = random_vector(case.matrix.nrows());
        let _ = solve_ldlt(&factorization, &b);
        profiler.end_section();

        // Report
        profiler.print_report();
        profiler.export_flamegraph(&format!("{}.svg", case.name));
    }
}
```

**Optimization Targets:**
1. Memory allocation patterns
2. Cache utilization
3. SIMD utilization (via faer)
4. Minimize copies
5. Thread synchronization overhead

**Success Criteria:**
- [ ] Identified top 3 bottlenecks
- [ ] Optimization opportunities documented
- [ ] Performance improvement plan created

#### 10.2: Memory Optimization
**Task:** Reduce memory footprint and improve locality

**Strategies:**

```rust
// Use arena allocators for frontal matrices
pub struct FrontalArena {
    memory_pool: Vec<f64>,
    allocations: Vec<AllocationInfo>,
}

impl FrontalArena {
    pub fn allocate_front(&mut self, nrows: usize, ncols: usize)
        -> FrontalMatrix;

    pub fn free_front(&mut self, front: FrontalMatrix);

    pub fn peak_usage(&self) -> usize;
}

// Memory-efficient factor storage
pub struct CompactFactorStorage {
    // Store only lower triangle
    // Pack 1×1 and 2×2 blocks efficiently
}
```

**Testing:**

```rust
#[test]
fn test_memory_usage() {
    for case in test_cases() {
        let before = get_memory_usage();

        let analysis = analyze_full(&case.matrix);
        let factorization = factor_tree(
            &analysis.assembly_tree,
            &case.matrix,
            FactorOptions::default()
        ).unwrap();

        let peak = get_peak_memory_usage();
        let after = get_memory_usage();

        println!("{}: peak={:.1}MB, final={:.1}MB, predicted={:.1}MB",
                 case.name,
                 peak / 1e6,
                 after / 1e6,
                 analysis.analysis.predicted_memory / 1e6);

        // Peak should be close to prediction
        let ratio = peak as f64 / analysis.analysis.predicted_memory as f64;
        assert!(ratio < 1.5, "Memory usage higher than predicted");
    }
}
```

**Success Criteria:**
- [ ] Memory usage predictable
- [ ] Peak memory within 50% of prediction
- [ ] No memory leaks

#### 10.3: API Refinement
**Task:** Polish public API based on feedback

**Documentation:**

```rust
/// Sparse symmetric indefinite linear system solver.
///
/// Uses LDL^T factorization with A Posteriori Threshold Pivoting (APTP)
/// for numerical stability on indefinite systems.
///
/// # Examples
///
/// Basic usage:
/// ```
/// use sparse_ldlt::SparseLDLT;
///
/// let matrix = /* load your matrix */;
/// let b = vec![1.0; matrix.nrows()];
///
/// let x = SparseLDLT::solve_matrix(matrix, &b)?;
/// ```
///
/// Reusing factorization:
/// ```
/// let mut solver = SparseLDLT::new(matrix.clone())?;
/// solver.factor(&matrix)?;
///
/// let x1 = solver.solve(&b1)?;
/// let x2 = solver.solve(&b2)?;  // Reuse factorization
/// ```
pub struct SparseLDLT { /* ... */ }

// Add builder pattern for options
impl SparseLDLT {
    pub fn builder() -> SparseLDLTBuilder {
        SparseLDLTBuilder::new()
    }
}

pub struct SparseLDLTBuilder {
    ordering: OrderingMethod,
    scaling: ScalingMethod,
    threshold: f64,
    parallelism: Parallelism,
}

impl SparseLDLTBuilder {
    pub fn ordering(mut self, method: OrderingMethod) -> Self {
        self.ordering = method;
        self
    }

    pub fn threshold(mut self, threshold: f64) -> Self {
        self.threshold = threshold;
        self
    }

    pub fn parallel(mut self, n_threads: usize) -> Self {
        self.parallelism = Parallelism::Rayon(n_threads);
        self
    }

    pub fn build(self, matrix: SparseColMat<f64>) -> Result<SparseLDLT> {
        // ...
    }
}
```

**Examples:**

```rust
// examples/basic_solve.rs
//! Basic sparse LDL^T solve example

use sparse_ldlt::*;

fn main() -> Result<()> {
    // Load matrix
    let matrix = io::read_matrix_market("test.mtx")?;

    // Create RHS
    let b = vec![1.0; matrix.nrows()];

    // Solve
    let x = SparseLDLT::solve_matrix(matrix, &b)?;

    println!("Solution norm: {:.6}", x.iter().map(|&v| v*v).sum::<f64>().sqrt());

    Ok(())
}
```

**Success Criteria:**
- [ ] API is intuitive
- [ ] Documentation comprehensive
- [ ] Examples cover common use cases
- [ ] Error messages actionable

#### 10.4: Comprehensive Testing
**Task:** Expand test coverage

**Additional Tests:**

```rust
// Property-based testing
#[cfg(test)]
mod property_tests {
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_solve_always_correct(
            matrix in arbitrary_symmetric_matrix(10, 20)
        ) {
            if let Ok(solver) = SparseLDLT::new(matrix.clone()) {
                let b = random_vector(matrix.nrows());
                let x = solver.solve(&b).unwrap();

                let residual = compute_residual_sparse(&matrix, &x, &b);
                prop_assert!(residual < 1e-6);
            }
        }
    }
}

// Fuzzing
#[test]
fn test_fuzz_inputs() {
    // Test with malformed inputs
    // - Invalid matrix patterns
    // - Singular matrices
    // - Extreme values
    // Should not panic
}

// Regression tests
#[test]
fn test_known_issues() {
    // Test cases that previously failed
    // Ensure they now work
}
```

**Success Criteria:**
- [ ] >90% code coverage
- [ ] Property tests pass
- [ ] No panics on invalid inputs
- [ ] All known issues resolved

### Phase 9 Exit Criteria

**Required Outcomes:**
1. Performance optimized
2. Memory usage efficient
3. API polished and documented
4. Comprehensive test suite

**Validation Questions:**
- Is performance competitive with SPRAL?
- Is the library ready for public use?
- Are there any remaining bugs?

**Final Checkpoint:**
- Run full benchmark suite
- Compare with SPRAL on all metrics
- Generate final report
- Create publication-ready plots

---

## Phase 11: Release Preparation

### Objectives
Prepare for public release: documentation, examples, tutorials, packaging.

### Deliverables

#### 11.1: Documentation
- [ ] README with quick start
- [ ] User guide
- [ ] API reference (rustdoc)
- [ ] Performance guide
- [ ] Migration guide (from other solvers)

#### 11.2: Examples
- [ ] Basic solve
- [ ] Interior point method integration
- [ ] Iterative refinement
- [ ] Multiple RHS
- [ ] Refactorization

#### 11.3: Benchmarks
- [ ] Comparison with SPRAL
- [ ] Performance plots
- [ ] Scaling studies

#### 11.4: Release
- [ ] Version 0.1.0 on crates.io
- [ ] CI/CD setup
- [ ] Issue templates
- [ ] Contributing guide

### Phase 11 Exit Criteria

**Required Outcomes:**
1. Published on crates.io
2. Documentation complete
3. Examples working
4. Community ready

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
