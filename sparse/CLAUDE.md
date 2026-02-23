# CLAUDE.md - Sparse Linear Algebra Solvers

This file provides guidance to Claude Code when working on sparse linear algebra solvers in this directory.

## Project Overview

This directory contains sparse linear algebra solver implementations for rivrs-linalg, focused on a Sparse Symmetric Indefinite Direct Solver (SSIDS) inspired by SPRAL. The algorithms are implemented using modern Rust with `faer` for high-performance linear algebra operations.

**Parent Project**: rivrs-linalg - Numerical Linear Algebra for Rivrs
**Domain**: Sparse direct solvers (SSIDS, LDL^T factorization, APTP pivoting)
**Current Status**: Phase 9.1c complete — assembly/extraction optimization + sub-phase profiling. c-71: 2.48× SPRAL (was 4.06×). Identified contribution block allocation as dominant bottleneck (40.1% + 33.3% of factor time). 65/65 SuiteSparse matrices passing.

### Development docs

- `docs/ssids-plan.md` — **Living document** describing the overall development plan. Update when the plan itself changes (new decisions, scope adjustments), not for arbitrary progress tracking.
- `docs/ssids-log.md` — Changelog for development updates, organized by phase. Records what was built, changed, fixed, and why.
- `docs/phase-8.1g-report.md` — Performance analysis report with workload distribution and Phase 8.2 parallelism recommendations.
- `docs/phase9/phase-9.1c-profiling-report.md` — Sub-phase profiling + perf stat analysis identifying contribution block allocation as dominant bottleneck.

## Licensing Strategy - Clean Room Implementation

**This is a clean room implementation to maintain licensing flexibility.**

- **SPRAL (BSD-3)**: May be freely consulted for implementation patterns
- **LAPACK (BSD-3)**: May be freely consulted for dense kernels
- **HSL documentation**: Consult for algorithm descriptions only
- **Primary sources**: Academic papers, textbooks
- **Always document** academic sources used for each implementation
- This approach allows rivrs-sparse to be released under Apache-2.0

## Reference Materials

Available in the parent `references/` directory (symlink to `/opt/references`, not git-tracked):

- **faer-rs/**: Pure Rust linear algebra library — primary dependency
- **lapack/**: Reference LAPACK implementation (BSD-licensed)
- **ssids/**: Academic literature converted to Markdown format
- **spral/**: SPRAL source code (BSD-3) — primary reference implementation

## Technology Stack

- **Rust 1.87+** (edition 2024)
- **faer 0.22** — dense and sparse linear algebra (CSC, elimination trees, AMD ordering, permutations, workspace management, matmul, triangular solve)
- **metis-sys 0.3** — METIS graph partitioning (vendored C source, used for nested dissection ordering)
- **serde + serde_json** — JSON serialization for matrix metadata, reference factorizations, Chrome Trace profiling export
- **rand + rand_distr** — random matrix generation (optional, behind `test-util` feature)
- **criterion 0.5** — benchmarking (dev dependency)
- **approx** — floating-point test comparisons (dev dependency)

### Cargo Feature Flags

| Feature | Purpose | Activates |
|---------|---------|-----------|
| `test-util` | Test infrastructure: random matrix generators, debug displays, profiling, benchmarking modules | `rand`, `rand_distr` |
| `diagnostic` | Per-supernode timing instrumentation in factorization. Zero overhead when disabled. | Timing fields on `PerSupernodeStats`/`FactorizationStats`, `ProfileSession` on `AptpNumeric` |

The `test-util` feature is automatically enabled in tests/benchmarks via the self-referential dev-dependency in `Cargo.toml`.

## Source Code Layout

```
src/
├── lib.rs              # Crate root: SolverPhase enum, module declarations
├── error.rs            # Error types (SolverError, thiserror)
├── validate.rs         # sparse_backward_error() — correctness validation
├── io.rs               # IO module root
├── io/                 # Matrix I/O
│   ├── mtx.rs          # MatrixMarket (.mtx) reader (mirrors for full symmetric CSC)
│   └── registry.rs     # Test matrix registry (metadata.json, CI-path fallback)
├── aptp/               # Core APTP solver
│   ├── mod.rs          # Public re-exports
│   ├── pivot.rs        # PivotType, Block2x2
│   ├── diagonal.rs     # MixedDiagonal (D factor: mixed 1x1/2x2, solve, inertia)
│   ├── inertia.rs      # Inertia (eigenvalue sign counts)
│   ├── perm.rs         # perm_from_forward()
│   ├── symbolic.rs     # AptpSymbolic (wraps faer's SymbolicCholesky + pivot buffer estimates)
│   ├── ordering.rs     # metis_ordering(), match_order_metis() (METIS + MC64 pipeline)
│   ├── matching.rs     # mc64_matching() — weighted bipartite matching & scaling
│   ├── factor.rs       # Dense APTP kernel: aptp_factor_in_place(), factor_inner (BLAS-3)
│   ├── numeric.rs      # AptpNumeric::factor() — multifrontal factorization loop
│   ├── solve.rs        # aptp_solve() — per-supernode forward/diagonal/backward solve
│   └── solver.rs       # SparseLDLT — user-facing API (analyze/factor/solve)
├── profiling/          # Performance profiling (behind test-util or diagnostic)
│   ├── session.rs      # ProfileSession, SectionGuard (RAII), FinishedSession
│   ├── section.rs      # Section timing data structures
│   ├── memory.rs       # MemoryTracker (RSS snapshots)
│   └── report.rs       # Chrome Trace JSON export
├── debug/              # Debug visualization (behind test-util)
│   ├── sparsity.rs     # SparsityDisplay (density chars, downsampling)
│   └── etree.rs        # ETreeDisplay (tree + stats)
├── benchmarking/       # Benchmark infrastructure (behind test-util or diagnostic)
│   ├── config.rs       # BenchmarkConfig
│   ├── baseline.rs     # Baseline management
│   ├── results.rs      # BenchmarkResult
│   ├── report.rs       # Report generation
│   ├── rss.rs          # read_peak_rss_kb()
│   └── traits.rs       # BenchmarkMatrix trait
└── testing/            # Test infrastructure (behind test-util)
    ├── harness.rs      # SolverTest trait, MockSolver
    ├── validator.rs    # NumericalValidator
    ├── cases.rs        # TestCaseFilter
    └── generators.rs   # Random matrix generators (PD, indefinite)
```

### Key Examples

```
examples/
├── profile_matrix.rs       # Single-matrix profiling: per-supernode timing + Chrome Trace (requires diagnostic)
├── baseline_collection.rs  # Structured JSON baseline collection (requires diagnostic)
├── workload_analysis.rs    # Workload distribution + parallelism classification (requires diagnostic)
├── export_frontal.rs       # Chrome Trace export for frontal matrices (requires diagnostic)
├── solve_timing.rs         # End-to-end solve timing
├── spral_comparison.rs     # Compare against SPRAL reference results
├── accuracy_benchmark.rs   # Backward error across SuiteSparse suite
└── front_sizes.rs          # Front size distribution analysis
```

## Commands Reference

### Building

```bash
cargo build                          # Standard build
cargo build --features diagnostic    # Build with per-supernode timing instrumentation
cargo build --release                # Optimized build (needed for large matrices)
```

Note: `dev` and `test` profiles already use `opt-level = 3` because faer's dense LA is unusably slow at opt-level 0.

### Testing

```bash
# Unit tests (358 tests, ~40s)
cargo test

# Unit tests with diagnostic feature (same tests, verifies cfg builds)
cargo test --features diagnostic

# SuiteSparse CI subset tests (10 small matrices, ~19MB) run as regular tests.
# Full SuiteSparse collection tests (65 matrices) are #[ignore]:
cargo test -- --ignored --test-threads=1

# Run a specific test by name
cargo test test_name

# Run tests in a specific module
cargo test aptp::factor
```

**Important**: Full SuiteSparse `--ignored` tests MUST use `--test-threads=1` to avoid memory pressure from concurrent large-matrix factorizations.

### Linting

```bash
cargo fmt --check                              # Check formatting
cargo clippy --all-targets                     # Lint without diagnostic
cargo clippy --all-targets --features diagnostic  # Lint with diagnostic (catches cfg issues)
```

### Benchmarking

```bash
# Criterion benchmarks (compile check only — useful for CI)
cargo bench --no-run

# Run Criterion benchmarks
cargo bench

# Baseline collection — structured JSON with per-phase timing, per-supernode stats, backward error
cargo run --example baseline_collection --features diagnostic --release -- --ci-only
cargo run --example baseline_collection --features diagnostic --release

# Compare against previous baseline
cargo run --example baseline_collection --features diagnostic --release -- --compare target/benchmarks/baselines/prev.json
```

Baseline JSON output goes to `target/benchmarks/baselines/baseline-<timestamp>.json`.

### Profiling & Analysis

```bash
# Profile a single matrix — phase timing, factor breakdown, top supernodes by time
cargo run --example profile_matrix --features diagnostic --release -- <matrix-name>
cargo run --example profile_matrix --features diagnostic --release -- --path <file.mtx>

# Export Chrome Trace for a single matrix (open in https://ui.perfetto.dev)
cargo run --example profile_matrix --features diagnostic --release -- <matrix-name> --trace /tmp/trace.json

# List available matrices
cargo run --example profile_matrix --features diagnostic --release -- --list

# Workload analysis — classifies matrices by parallelism strategy (TreeLevel/IntraNode/Mixed)
cargo run --example workload_analysis --features diagnostic --release

# Export assembled frontal matrix for SPRAL comparison (open in chrome://tracing or Perfetto)
cargo run --example export_frontal --features diagnostic --release

# Front size distribution
cargo run --example front_sizes --release

# Solve timing
cargo run --example solve_timing --release

# Per-supernode statistics
cargo run --example supernode_stats --release
```

### CI (GitHub Actions)

The CI workflow (`.github/workflows/ci.yml` at repo root) runs on push to `main` and PRs to `main`/`ssids`:

- `cargo test --all-targets` on stable + MSRV (1.87)
- `cargo fmt --check` + `cargo clippy --all-targets -- -D warnings`
- `cargo doc --no-deps` with `-D warnings`
- `cargo bench --no-run` (compile check)

CI does NOT run `--ignored` tests (SuiteSparse matrices not available in CI environment).

## Algorithm Architecture

The solver implements a multifrontal LDL^T factorization with A Posteriori Threshold Pivoting (APTP). The pipeline is:

1. **Analyze** (`SparseLDLT::analyze_with_matrix`): Ordering (MC64 matching + METIS nested dissection), symbolic factorization (elimination tree, supernodal structure, column counts)
2. **Factor** (`SparseLDLT::factor` → `AptpNumeric::factor`): Multifrontal assembly + dense APTP kernel per supernode (BLAS-3 pipeline: factor_block_diagonal → TRSM → threshold check → GEMM)
3. **Solve** (`SparseLDLT::solve` → `aptp_solve`): Per-supernode forward substitution, diagonal solve (MixedDiagonal), backward substitution

### Key Algorithm: APTP (A Posteriori Threshold Pivoting)

Unlike traditional threshold pivoting (which decides pivots before elimination), APTP:
1. Performs elimination optimistically (assuming 1x1 pivots)
2. Checks stability after the fact
3. Falls back to 2x2 (Bunch-Kaufman) pivots or delays columns when needed

The two-level architecture uses TPP (Threshold Partial Pivoting) as the primary strategy for small blocks (< 32 columns) and complete pivoting with BLAS-3 blocking for larger fronts.

### Key Types

- `SparseLDLT` — user-facing solver (analyze/factor/solve), handles MC64 scaling
- `AptpSymbolic` — wraps faer's `SymbolicCholesky` + APTP-specific pivot buffer estimates
- `AptpNumeric` — multifrontal factorization result (per-supernode L, D, permutations)
- `MixedDiagonal` — the D factor with mixed 1x1/2x2 blocks (`solve_in_place`, `compute_inertia`)
- `OrderingStrategy` — `Amd`, `Metis`, `MatchOrderMetis` (default)
- `FactorizationStats` — aggregate pivot counts, delays, max front size, timing (with `diagnostic`)

## faer Integration: Transparent Composition

This library is a **specialized extension of faer**, not a competing implementation.
The guiding principle is **transparent composition**: use the highest-level faer API
that gives us what we need, expose faer types at our boundary, and only add new types
for concepts faer doesn't have.

**Concrete guidelines:**

- **Accept faer types as inputs**: `SparseColMat<usize, f64>`, `PermRef`, `SymmetricOrdering`, `SymbolicSparseColMatRef`
- **Return faer types where they fit**: permutations as `Perm<usize>` (not custom wrappers), sparse matrices as `SparseColMat`
- **Compose, don't wrap**: our structs contain faer results as fields with accessor methods that delegate, rather than copying data into our own arrays (e.g., `AptpSymbolic` stores `SymbolicCholesky<usize>` and delegates `perm()`, `predicted_nnz()`)
- **Only add types for genuinely new concepts**: `PivotType`, `MixedDiagonal`, `AptpSymbolic` — things faer doesn't model
- **Prefer high-level faer APIs over low-level primitives**: use `factorize_symbolic_cholesky` over manual `prefactorize_symbolic_cholesky` + `factorize_simplicial_symbolic_cholesky`, unless we need fine-grained control

**Anti-patterns to avoid:**

- Defining a custom `Permutation` struct when `faer::perm::Perm<usize>` works
- Adding type aliases like `AptpMatrix<T>` for `SparseColMat<usize, f64>` — faer types are already well-named
- Re-implementing CSC storage, elimination tree computation, or AMD ordering
- Hiding faer types behind opaque wrappers that prevent interop with the faer ecosystem

## Code Architecture

**Numerical Considerations:**
- Prioritize numerical stability for indefinite systems
- Use structured decompositions rather than direct inversion
- Implement overflow prevention and condition monitoring
- Support both positive definite (fast path) and indefinite (APTP) factorization

**API Design:**
- Three-phase API: analyze → factorize → solve
- Symbolic analysis result is reusable across multiple factorizations with same sparsity pattern
- Configuration via `AnalyzeOptions`, `FactorOptions` structs
- `MemStack` for solve workspace (factor allocates internally)

**Error Handling:**
- Use Result types for all fallible operations
- Distinguish structural singularity (analysis) from numerical singularity (factorization)
- Provide informative diagnostics (pivot delays, factorization statistics)

**CSC Storage Convention:**
- All matrices store FULL symmetric CSC (both upper and lower triangles)
- `.mtx` reader mirrors entries to produce full symmetric storage
- `sparse_backward_error` and matvec helpers must NOT mirror entries — regular matvec is correct
- `scatter_original_entries` in `numeric.rs` has upper-triangle skip logic to avoid double-counting during frontal matrix assembly

## Development Workflow

When implementing a new component:

1. Review SPRAL source (BSD-3) to understand the algorithm structure
2. Identify relevant academic papers and textbooks
3. Review LAPACK source for any dense kernel patterns needed
4. Design idiomatic Rust API leveraging faer's sparse infrastructure
5. Implement using academic references and permissively-licensed code
6. Document academic references and SPRAL routines consulted
7. Add comprehensive tests (hand-constructed matrices, SuiteSparse collection)
8. Benchmark against SPRAL and other solvers

## Rust Guidelines

- Sort imports by: std, external, workspace, crate, super
- Limit scope to minimum necessary (private over `pub(crate)` over `pub`, etc.)
- No `.unwrap()` in non-test code
- Pay special attention to Rust edition-specific reserved keywords (e.g., `gen` is reserved in Rust 2024 — use `rng.r#gen::<f64>()`)
- When using faer or other Rust crates, always verify API types by reading the actual source/docs before using them. Do not assume tuple-based constructors — check for dedicated types like `Triplet`.

Preferences (can be violated if needed):
- Immutability over mutation
- Iterators over manual loops
- Enums over dynamic dispatch (e.g. `dyn` trait)
- Prefer `Result` propagation over `.expect()`
- Prefer `thiserror` over custom `Result`/`Error` traits
- Importing with `use` over inline imports (exception: module imports are okay if several functions or types are used)
- Try to avoid allocation in performance-critical code

## Git Commit Practices

- **Commit frequently**: After completing logical units (phase completion, tests passing, new module)
- **Commit after verification**: Always commit when tests pass or milestones are verified
- **Descriptive messages**: Use clear commit messages
- **Never accumulate**: Don't accumulate large amounts of uncommitted work

## Testing Strategy

### Validation strategy

Primary correctness validation (established in Phase 0.3 decision):

1. **Reconstruction tests**: `||P^T A P - L D L^T|| / ||A|| < 10^-12` (primary oracle)
2. **Backward error**: `||Ax - b|| / (||A|| ||x|| + ||b||) < 5e-11` (SPRAL's threshold)
3. **Hand-constructed matrices**: 15 matrices with analytically known factorizations
4. **Property-based tests**: Inertia, symmetry preservation, permutation validity
5. **SPRAL comparison**: Deferred (performance benchmarking, large-matrix inertia)

### Ordering for tests and benchmarks

**Use `MatchOrderMetis` (the default) for all production and integration testing.**
Phase 7 benchmarking showed that plain METIS causes massive pivot delays on hard
indefinite matrices (bratu3d: 53K delays, backward error 5e-3), while MatchOrderMetis
(MC64 matching+scaling) eliminates the delays (1 delay, backward error 1e-9).
`OrderingStrategy::Amd` should never be used in production; it exists only for
unit tests of the symbolic analysis and factorization kernel on small matrices.

### Test matrix storage (three-tier, no Git LFS)

| Tier | Path | Contents | Size | Git status |
|------|------|----------|------|------------|
| Hand-constructed | `test-data/hand-constructed/` | 15 matrices + .json factorizations | ~144KB | Tracked |
| CI subset | `test-data/suitesparse-ci/` | 10 small SuiteSparse matrices (easy+hard indefinite) | ~19MB | Tracked |
| Full suite | `test-data/suitesparse/` | 67 SuiteSparse matrices | ~500MB | Gitignored |

`test-data/metadata.json` is the matrix registry with `ci_subset: true` flags. The `registry.rs` module has CI-path fallback logic.

### Test infrastructure

- `SolverTest` trait with `MockSolver` for pre-solver validation
- `NumericalValidator` with configurable tolerances
- `TestCaseFilter` for composable test case selection
- Random matrix generators (PD and indefinite) behind `test-util` feature

### Testing discipline

1. **Validation proportional to novelty**: Code that delegates to faer gets sanity checks. Code written from scratch gets mathematical proof-level tests with exact expected values.
2. **Refactor for testability**: Extract non-trivial logic into standalone functions with minimal inputs.
3. **Regression tests before refactoring**: Encode current correct behavior with exact expected values before optimizing.
4. **Cross-validation with faer**: Assert consistency when both our code and faer compute overlapping quantities.
5. **Property-based testing**: For implementation-dependent outputs, test structural properties rather than exact values.

## Documentation Standards

- Document the mathematical operation performed
- **Cite academic references** used during implementation
- Include references to sparse direct solver literature
- Provide usage examples
- Document numerical stability characteristics and failure modes
- Cross-reference SPRAL routine names for users migrating from Fortran
- Add proper attribution as dictated by open-source community standards

## Current Implementation Status

**Completed:**
- Phase 0: Foundation — literature review, test matrix collection (82 matrices), repository setup
- Phase 1: Infrastructure — test harness, benchmarking framework (Criterion), CI (GitHub Actions), profiling & debug tools
- Phase 2: APTP data structures — PivotType, Block2x2, MixedDiagonal, Inertia, perm_from_forward
- Phase 3: Symbolic analysis — AptpSymbolic wrapping faer's SymbolicCholesky + pivot buffer estimates
- Phase 4: Ordering & preprocessing — METIS ordering, MC64 matching & scaling, match-order condensation pipeline
- Phase 5: Dense APTP kernel — aptp_factor_in_place with swap-delayed-to-end architecture
- Phase 6: Multifrontal numeric factorization — AptpNumeric, FrontalMatrix, extend-add, scatter
- Phase 7: Triangular solve & solver API — aptp_solve, SparseLDLT (analyze/factor/solve), OrderingStrategy
- Phase 8.1: Two-level APTP + BLAS-3 refactoring — TPP for small fronts, complete pivoting for large fronts, critical extract_front_factors bug fix. 65/65 SuiteSparse matrices pass.
- Phase 8.1g: Sequential profiling & optimization — per-supernode timing instrumentation, allocation hotspot fixes, baseline collection tool, workload analysis (41/65 IntraNode, 19 Mixed, 5 TreeLevel)
- Phase 8.2: Parallel factorization & solve — rayon + faer Par, intra-node BLAS parallelism (TRSM/GEMM) for large fronts, tree-level par_iter for independent subtrees, parallel diagonal solve. All safe Rust (no unsafe). parallel_scaling.rs benchmark tool.
- Phase 9.1a: Supernode amalgamation — SPRAL-style post-faer merge pass (do_merge predicate, nemin=32 default). c-71: 35K→6.4K supernodes. 65/65 SuiteSparse pass. Non-contiguous merges via owned_ranges on SupernodeInfo. Configurable via AnalyzeOptions.nemin.
- Phase 9.1b: Workspace reuse — Two-tier pre-allocated workspace (FactorizationWorkspace in numeric.rs, AptpKernelWorkspace in factor.rs). FrontalMatrix changed from owned to borrowed types. Thread-local workspace via Cell for parallel path. CI subset: median 1.16× factor speedup, 24% RSS reduction, bit-exact backward errors.
- Phase 9.1c: Assembly & extraction optimization — Precomputed scatter maps (AssemblyMaps), bulk column-slice copies, fill(0.0) zeroing, sub-phase timing instrumentation. c-71: 4.06×→2.48× SPRAL. Sub-phase profiling identified contribution block allocation as dominant bottleneck (extract_contrib 40.1% + extend-add 33.3% of factor time). perf stat: 644B dTLB misses, 3.1s/32% sys time from mmap churn.

**Next:**
- Phase 9.1d: Contribution workspace reuse — eliminate per-supernode Mat allocation for contribution blocks
- Phase 9.1e: Small leaf subtree fast path
- Phase 9.2: Release preparation (docs, examples, crates.io)

## Recent Changes
- 023-contrib-workspace-reuse: Added Rust 1.87+ (edition 2024) + faer 0.22, rayon 1.x, serde/serde_json (diagnostic export)
- 022-assembly-extraction-opt: Added Rust 1.87+ (edition 2024) + faer 0.22 (existing — no new deps)
- 021-workspace-reuse: Added N/A (internal optimization; no new technologies)

## Active Technologies
- Rust 1.87+ (edition 2024) + faer 0.22, rayon 1.x, serde/serde_json (diagnostic export) (023-contrib-workspace-reuse)
- N/A (in-memory buffers; no persistent storage changes) (023-contrib-workspace-reuse)
