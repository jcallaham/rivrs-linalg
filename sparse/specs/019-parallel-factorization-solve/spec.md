# Feature Specification: Parallel Factorization & Solve (Phase 8.2)

**Feature Branch**: `019-parallel-factorization-solve`
**Created**: 2026-02-21
**Status**: Draft
**Input**: User description: "Implement Phase 8.2 — parallel factorization and solve for the SSIDS sparse symmetric indefinite direct solver"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Parallel Factorization of Large Frontal Matrices (Priority: P1)

A library user has a large sparse symmetric indefinite system (e.g., FEM/CFD problems with n > 10,000) where factorization is dominated by a handful of large dense fronts. The user passes a parallelism parameter to the factorization API and observes a significant speedup on the dense BLAS-3 operations (TRSM, GEMM) within those large fronts, without any change to their existing code beyond adding the parallelism parameter.

**Why this priority**: Phase 8.1g profiling shows 41/65 (63%) of SuiteSparse test matrices are "IntraNode" — the top 10% of supernodes by front size account for >80% of total factorization time. Parallelizing the dense operations within these large fronts is the highest-impact optimization. Matrices like `sparsine` (max_front=11,125, 112s sequential factor time) would benefit enormously.

**Independent Test**: Can be fully tested by factoring any large SuiteSparse matrix (e.g., `sparsine`, `nd12k`) with 1, 2, and 4 threads and measuring wall-clock speedup. Delivers immediate performance value for the majority of real-world matrices.

**Acceptance Scenarios**:

1. **Given** a large sparse indefinite matrix with at least one frontal matrix of dimension > 500, **When** the user calls `SparseLDLT::factor` with parallel execution enabled (e.g., 4 threads), **Then** the factorization completes faster than the sequential baseline and produces numerically identical results.
2. **Given** the same matrix and parallelism configuration, **When** the factorization is run multiple times, **Then** the results are deterministic (identical pivot sequences, identical factors).
3. **Given** a small matrix where all fronts are < 100, **When** the user enables parallelism, **Then** the factorization still completes correctly with negligible overhead compared to sequential.

---

### User Story 2 - Tree-Level Parallel Factorization (Priority: P2)

A library user has a matrix where factorization work is distributed across many medium-sized fronts rather than concentrated in a few large ones (the "Mixed" or "TreeLevel" workload pattern). The solver processes independent subtrees of the assembly tree in parallel, achieving speedup even when individual fronts are too small for intra-node BLAS parallelism to be effective.

**Why this priority**: 19/65 (29%) of SuiteSparse matrices are classified as "Mixed" — moderate workload concentration (30-80% in top 10% of fronts). These matrices have many medium-sized fronts (100-500) rather than a few huge ones. Both tree-level and intra-node parallelism contribute to speedup. 5 additional matrices (8%) are "TreeLevel" where this is the only effective parallelism strategy. Tree-level parallelism is more complex to implement (dependency tracking, contribution synchronization) but necessary for comprehensive coverage.

**Independent Test**: Can be tested by factoring Mixed-class matrices (e.g., `G3_circuit`, `ldoor`, `msdoor`) and TreeLevel-class matrices (e.g., `bloweybq`, `dixmaanl`) with multiple threads and measuring speedup. Delivers value for matrices where intra-node parallelism alone is insufficient.

**Acceptance Scenarios**:

1. **Given** a matrix with a wide assembly tree (many independent subtrees), **When** the user enables parallel factorization, **Then** independent subtrees are processed concurrently and wall-clock time decreases relative to sequential execution.
2. **Given** a matrix where tree-level parallelism is the dominant strategy (TreeLevel class), **When** 4 threads are used, **Then** measurable speedup is observed despite small individual front sizes.
3. **Given** any matrix, **When** both tree-level and intra-node parallelism are active, **Then** they compose correctly without data races or synchronization errors, and results remain deterministic and identical to sequential.

---

### User Story 3 - Parallel Triangular Solve (Priority: P3)

A library user performs repeated solves after a single factorization (e.g., in an iterative refinement or interior-point optimization loop). The solve phase benefits from parallelism in two ways: (a) diagonal solve is fully parallelized across supernodes, and (b) forward/backward substitution exploits both tree-level independence and intra-node dense parallelism for large supernodes.

**Why this priority**: Solve time is typically 1-5% of factor time (per Phase 8.1g baselines), so the absolute impact is lower than factorization parallelism. However, for applications with many sequential solves (e.g., iterative refinement, multiple right-hand sides), cumulative solve time becomes significant. The diagonal solve is embarrassingly parallel; forward/backward solve parallelism is limited by tree dependencies but still beneficial for wide trees.

**Independent Test**: Can be tested by timing solve with 1 vs 4 threads on matrices with different tree structures. The diagonal solve parallelism alone should show improvement; forward/backward solve improvement depends on tree width.

**Acceptance Scenarios**:

1. **Given** a factorized matrix, **When** the user calls `SparseLDLT::solve_in_place` with parallel execution enabled, **Then** the solve produces results identical to sequential solve.
2. **Given** a matrix with a wide assembly tree, **When** parallel solve is enabled, **Then** wall-clock solve time decreases relative to sequential for the forward and backward substitution phases.
3. **Given** a matrix with large supernodes, **When** parallel solve is enabled, **Then** the per-supernode dense triangular solve and matrix-vector operations within forward/backward substitution exploit parallel BLAS.

---

### User Story 4 - Performance Benchmarking & Scaling Report (Priority: P4)

A library developer or user evaluates the parallel solver's performance characteristics across a range of matrix types and thread counts. A benchmarking tool produces a scaling report that compares parallel performance against the Phase 8.1g sequential baselines, characterizes load balancing, and identifies the effective parallelism strategy for different matrix classes.

**Why this priority**: Benchmarking is essential for validating that parallelism delivers the expected speedup and for identifying remaining bottlenecks. It builds on the existing `baseline_collection.rs` and `workload_analysis.rs` tools from Phase 8.1g. The report informs whether the >=3x on 4 cores target is met and guides future optimization.

**Independent Test**: Can be tested by running the benchmark tool on the CI subset and full SuiteSparse suite and verifying that the report contains the expected data fields (speedup, efficiency, thread scaling).

**Acceptance Scenarios**:

1. **Given** the parallel solver and sequential baselines from Phase 8.1g, **When** the benchmark tool is run with 1, 2, 4, and 8 threads, **Then** it produces a structured report with per-matrix speedup, parallel efficiency, and thread scaling curves.
2. **Given** the benchmark results, **When** results are compared across matrix classes (IntraNode, Mixed, TreeLevel), **Then** the report identifies which parallelism strategy (intra-node, tree-level, or combined) is most effective for each class.

---

### Edge Cases

- What happens when the thread count is 1? The solver MUST behave identically to the current sequential implementation with no measurable overhead.
- What happens when the matrix is very small (n < 100)? Parallel overhead MUST NOT cause slowdown; the solver should effectively fall back to sequential behavior for small problems.
- What happens when a supernode has delayed columns that propagate to a parent being processed on a different thread? Contribution blocks and delayed column data MUST be correctly synchronized across threads.
- What happens when the assembly tree is a single chain (no independent subtrees)? Tree-level parallelism provides no benefit, but the solver MUST still function correctly. Intra-node parallelism should still apply.
- What happens during the solve phase when the RHS vector is modified by concurrent supernodes? The solve MUST prevent data races on the shared RHS array through appropriate scheduling or synchronization.
- What happens if the user requests more threads than available cores? The solver MUST handle this gracefully (the thread pool should manage oversubscription without correctness issues).

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The solver MUST accept a parallelism configuration parameter via faer's `Par` enum added as a field on `FactorOptions` and `SolverOptions` (and propagated to internal routines). `Par::Seq` (the default) provides sequential execution; `Par::Rayon(n)` enables parallel execution with `n` threads. This maintains transparent composition with faer's parallelism model.
- **FR-002**: The dense APTP kernel MUST support parallel execution of BLAS-3 operations (triangular solve, matrix multiplication) within individual frontal matrices, using the parallelism configuration to control thread count. Intra-node parallelism MUST be skipped (use `Par::Seq`) for fronts with dimension < 256 to avoid thread dispatch overhead. This threshold aligns with the existing `outer_block_size` default.
- **FR-003**: The multifrontal factorization loop MUST support processing independent subtrees of the assembly tree in parallel, with correct synchronization of contribution blocks between child and parent supernodes. Tree-level scheduling uses a bottom-up recursive `rayon::scope` strategy: children of each supernode are spawned as parallel Rayon tasks within a scope, the scope barrier ensures all children complete (contribution blocks ready) before the parent processes itself. Dependency correctness is enforced structurally by the scope barrier (no explicit synchronization primitives needed). This pattern handles N-ary trees naturally (unlike `rayon::join` which takes exactly 2 closures).
- **FR-004**: The triangular solve phase MUST support parallel execution. At minimum, the diagonal solve (which is fully independent across supernodes) MUST be parallelized. Forward and backward substitution SHOULD exploit tree-level parallelism where dependencies allow.
- **FR-005**: Parallel factorization and solve MUST produce deterministic results at a fixed thread count — identical outputs given identical inputs and the same `Par` configuration, regardless of thread scheduling order.
- **FR-006**: Parallel results MUST satisfy the same backward error threshold as sequential results (< 5e-11). Bitwise identity across different thread counts is NOT required, because faer's parallel BLAS may use different floating-point summation order than sequential. Pivot sequences may differ at ULP-level tie-breaks.
- **FR-007**: A benchmarking tool MUST produce a parallel scaling report comparing performance across thread counts (1, 2, 4, 8) against Phase 8.1g sequential baselines, covering factorization and solve phases separately.
- **FR-008**: The solver MUST NOT introduce data races, deadlocks, or undefined behavior under any parallelism configuration or input matrix.
- **FR-009**: Parallel overhead for small problems (where parallelism provides no benefit) MUST be negligible — no measurable slowdown compared to sequential execution on matrices with n < 1,000.
- **FR-010**: The existing sequential API MUST remain backward-compatible. Users who do not opt into parallelism MUST experience no behavioral change.

### Key Entities

- **Parallelism Configuration**: faer's `Par` enum, added as a `par` field on `FactorOptions` and `SolverOptions`. `Par::Seq` (default) for sequential, `Par::Rayon(n)` for parallel with `n` threads. Propagated to both intra-node BLAS operations and tree-level scheduling.
- **Assembly Tree Schedule**: Bottom-up recursive traversal using `rayon::scope` at tree forks. Children of each supernode execute as parallel Rayon tasks within a scope; the scope barrier ensures all children complete before the parent proceeds. Rayon's work-stealing handles load balancing. Handles N-ary trees naturally. Matches SPRAL's OpenMP task model (Duff 2020, Figure 3) translated to Rayon's fork-join.
- **Contribution Synchronization**: The mechanism ensuring that a parent supernode's frontal matrix assembly is complete (all children's contribution blocks have been received) before the parent begins factorization.
- **Parallel Scaling Report**: A structured output (extending Phase 8.1g baselines) containing per-matrix speedup, parallel efficiency, load balance metrics, and thread scaling curves for both factorization and solve.

## Assumptions

- The Rayon work-stealing thread pool is the appropriate parallelism runtime for tree-level scheduling in Rust. This is consistent with the broader Rust ecosystem and avoids introducing custom thread pool management. Rayon is an unconditional dependency (not feature-gated) to avoid conditional compilation complexity.
- faer's `Par` enum provides sufficient control for intra-node BLAS-3 parallelism and composes correctly with an outer Rayon-based tree-level scheduler.
- The existing profiling infrastructure (TLS-based `ProfileSession`, `flush_thread_events`) is already thread-safe and requires no structural changes for parallel execution.
- Bitwise determinism at a fixed thread count is achievable because: (a) the assembly tree schedule is fixed (not work-stealing-order-dependent for numerics), (b) intra-node BLAS operations in faer produce deterministic results for a given thread count, and (c) pivot decisions depend only on matrix values, not thread scheduling. Bitwise identity across different thread counts is NOT assumed — faer's parallel BLAS may change floating-point summation order.
- The Docker development environment has sufficient cores (at least 4) for meaningful parallel scaling tests.

## Algorithm References

The following academic papers inform the parallelism strategy. Markdown versions are available in the reference library:

- **Duff, Hogg, & Lopez (2020)** — "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting" (`/workspace/rivrs-linalg/references/ssids/duff2020.md`). Central reference for APTP algorithm. Section 5 describes the multifrontal implementation using OpenMP 4 tasks. Section 6 describes the high-level parallel strategy with NUMA-aware tree partitioning (Algorithm 6.1). Figure 3 shows task-based pseudocode with dependency clauses. Figure 4 shows nested parallel regions.
- **Hogg, Ovtchinnikov, & Scott (2016)** — "A Sparse Symmetric Indefinite Direct Solver for GPU Architectures" (`/workspace/rivrs-linalg/references/ssids/hogg2016.md`). Describes tree parallelism vs node parallelism tradeoffs (Section 2.5). Section 3.1 covers sparse assembly synchronization. Section 3.5 covers parallel solve with tree and node parallelism.
- **Duff & Pralet (2005)** — "Strategies for Scaling and Pivoting for Sparse Symmetric Indefinite Problems" (`/workspace/rivrs-linalg/references/ssids/duff2005.md`). Section 2 provides overview of multifrontal approach. Demonstrates how MC64 + compressed-graph ordering reduces delayed pivots, which is critical because delays cause extra work and tree imbalance.
- **Davis, Rajamanickam, & Sid-Lakhdar (2016)** — "A Survey of Direct Methods for Sparse Linear Systems" (`/workspace/rivrs-linalg/references/ssids/davis2016.md`). Section 11 covers multifrontal methods comprehensively. Section 12.2 covers parallel triangular solve.
- **Duff & Reid (1983)** — Assembly tree parallelism (cited in ssids-plan.md).
- **Liu (1992)** — Level-set scheduling for multifrontal methods (cited in ssids-plan.md).

## Phase 8.1g Profiling Context

The following findings from Phase 8.1g (see `dev/phase-8.1g-report.md`) directly inform the parallelism strategy:

**Workload Distribution** (65 SuiteSparse matrices, MatchOrderMetis ordering):
- **IntraNode** (41 matrices, 63%): Top 10% of supernodes by front size account for >80% of factorization time. Max front sizes range from 577 to 11,125. These benefit most from parallel dense BLAS within large fronts.
- **Mixed** (19 matrices, 29%): Moderate concentration (30-80% in top 10%). Many medium-size fronts (100-500). Both tree-level and intra-node parallelism contribute.
- **TreeLevel** (5 matrices, 8%): Many small supernodes (max_front 10-12), very even workload. Tree-level parallelism is the only effective strategy, but these are already fast.

**Recommended Phase 8.2 Order** (from report Section 7):
1. Intra-node BLAS-3 parallelism first (low complexity, high impact on 63% of matrices)
2. Tree-level scheduling second (higher complexity, benefits remaining 29%)
3. Benchmark both strategies across representative matrices from each class

**Sequential Baselines** (CI subset, representative):
- `sparsine`: 112s factor (max_front=11,125) — primary intra-node target
- `ncvxqp3`: 4.4s factor — moderate fronts, mixed parallelism
- `bloweybq`: 17.5ms factor (max_front=10-12) — tree-level only
- All backward errors well below SPRAL's 5e-11 threshold

## Clarifications

### Session 2026-02-21

- Q: What shape should the parallelism configuration API take? → A: Use faer's `Par` enum directly as a field on `FactorOptions` and `SolverOptions` (transparent composition with faer).
- Q: Should parallel results be bitwise-identical across thread counts (FR-006)? → A: No. Relax to backward-error equivalence across thread counts; bitwise identity required only at fixed thread count (FR-005).
- Q: What front size threshold should gate intra-node BLAS parallelism? → A: 256 (matches outer_block_size default; below this, use Par::Seq regardless of user setting).
- Q: What tree-level scheduling strategy for parallel factorization? → A: Bottom-up recursive `rayon::scope` at tree forks (children as parallel tasks within scope, scope barrier ensures completion, then parent processes). Handles N-ary trees naturally.
- Q: Should Rayon be feature-gated or unconditional? → A: Unconditional (always-on dependency, no cfg gates).

## Scope & Dependencies

**In scope**:
- Adding parallelism configuration to factor and solve APIs
- Intra-node BLAS-3 parallelism in the dense APTP kernel
- Tree-level parallelism for independent subtrees in the multifrontal factorization loop
- Parallel diagonal solve (embarrassingly parallel)
- Forward/backward solve parallelism (tree-level and/or intra-node)
- Parallel scaling benchmark tool and report
- Thread-safety of all shared data structures in the factorization and solve paths

**Out of scope**:
- GPU acceleration (future work, if ever)
- NUMA-aware memory allocation or thread pinning (Phase 9 or beyond; Duff 2020 Section 6 describes this but it adds significant complexity)
- Distributed-memory parallelism (MPI)
- Changes to the ordering or symbolic analysis phases (already sequential and relatively fast)
- Arena memory allocation for frontal matrices (Phase 9.1)
- Iterative refinement (Phase 9.1)

**Dependencies**:
- Phase 8.1g sequential baselines (COMPLETE — provides comparison targets)
- Phase 8.1 two-level APTP kernel (COMPLETE — provides the dense kernel to parallelize)
- Phase 7 triangular solve and solver API (COMPLETE — provides the solve path to parallelize)
- Rayon crate for tree-level work-stealing thread pool (unconditional dependency, not feature-gated)
- faer's `Par` enum for intra-node BLAS parallelism

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Factorization of large IntraNode-class matrices (max_front > 500) achieves at least 3x speedup on 4 cores compared to sequential baselines from Phase 8.1g.
- **SC-002**: Parallel factorization achieves at least 75% parallel efficiency on 4 cores for IntraNode-class matrices (i.e., speedup >= 3.0 on 4 threads).
- **SC-003**: Mixed-class matrices show measurable speedup (>1.5x on 4 cores) from combined tree-level and intra-node parallelism.
- **SC-004**: Parallel results are deterministic at a fixed thread count — running the same factorization and solve twice with the same `Par` configuration produces bitwise-identical outputs.
- **SC-005**: Parallel results satisfy the same backward error threshold (< 5e-11) as sequential results. Bitwise identity across thread counts is not required; numerical equivalence (same order-of-magnitude backward error) suffices.
- **SC-006**: All 65 SuiteSparse matrices pass correctness checks (backward error < 5e-11) under parallel execution, with no regressions from sequential.
- **SC-007**: Parallel overhead on small matrices (n < 1,000) is negligible — no more than 5% slowdown compared to sequential.
- **SC-008**: A parallel scaling report is produced covering the full SuiteSparse suite with 1, 2, 4, and 8 thread configurations, including per-matrix speedup, parallel efficiency, and workload class analysis.
