# Feature Specification: Sequential Profiling & Optimization (Phase 8.1g)

**Feature Branch**: `018-sequential-profiling-optimization`
**Created**: 2026-02-20
**Status**: Draft
**Input**: User description: "Phase 8.1g: Sequential Profiling and Optimization — Profile the complete solver pipeline, establish performance baselines, and optimize sequential bottlenecks before adding parallelism"

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Instrument Factorization Hot Path (Priority: P1)

A solver developer needs to understand *where time is spent* during multifrontal factorization so that optimization and parallelism efforts target the right bottlenecks. Today, only phase-level timing (analyze / factor / solve) is available. Per-supernode timing — breaking down assembly, APTP kernel execution, and contribution extraction — does not exist.

The developer instruments the factorization loop with the existing `ProfileSession` / `SectionGuard` infrastructure and generates Chrome Trace output for representative matrices. The trace shows a timeline of per-supernode work, revealing whether a few large fronts dominate or many small fronts contribute equally.

**Why this priority**: Without per-supernode profiling data, all parallelism decisions in Phase 8.2 are guesswork. This is the foundational deliverable that every other story depends on.

**Independent Test**: Can be fully tested by factoring a CI-subset matrix with profiling enabled and verifying that the Chrome Trace output contains per-supernode timing sections with correct nesting hierarchy.

**Acceptance Scenarios**:

1. **Given** the existing `ProfileSession` and `SectionGuard` infrastructure, **When** a developer enables profiling and factors any SuiteSparse CI matrix, **Then** the resulting Chrome Trace JSON contains per-supernode timing entries for at least: frontal matrix assembly, APTP kernel, and contribution extraction.
2. **Given** a Chrome Trace output file, **When** opened in a trace viewer, **Then** the developer can identify which supernodes consume the most time and see the hierarchical breakdown of work within each supernode.
3. **Given** profiling instrumentation in the factorization loop, **When** all 65 SuiteSparse matrices are factored, **Then** no correctness regression occurs (reconstruction < 1e-12, backward error unchanged from Phase 8.1f baselines).

---

### User Story 2 — Reduce Allocation Pressure in Factorization (Priority: P2)

A solver developer identifies that the factorization loop performs many short-lived heap allocations per supernode and per inner block — temporary vectors for row permutation, backup matrices for pivot restoration, and workspace matrices for BLAS-3 operations. These allocations cause unnecessary memory-management overhead in the sequential solver and would cause contention under threads.

The developer audits and optimizes the top allocation hotspots so that frequently-allocated temporary buffers are hoisted out of tight loops and reused across iterations.

**Why this priority**: Allocation reduction both improves sequential performance and is a prerequisite for contention-free parallel factorization. However, it depends on profiling data (Story 1) to confirm which allocations actually matter.

**Independent Test**: Can be tested by running the full SuiteSparse suite before and after optimization, comparing allocation counts and peak RSS, and verifying no correctness regressions.

**Acceptance Scenarios**:

1. **Given** the current factorization code with per-row and per-block temporary allocations, **When** the top allocation hotspots are optimized, **Then** the number of heap allocations during factorization decreases measurably (validated by profiling or allocation counting).
2. **Given** allocation optimizations applied, **When** the full 65-matrix SuiteSparse suite is factored, **Then** peak RSS does not increase and all correctness tests pass (reconstruction < 1e-12, backward error < 5e-11 for MatchOrderMetis matrices).
3. **Given** a large matrix with many supernodes (e.g., bloweybq, bratu3d), **When** factored before and after optimization, **Then** factorization time does not regress and allocation pressure is reduced.

---

### User Story 3 — Establish Sequential Performance Baselines (Priority: P3)

A solver developer needs documented sequential performance baselines — factor time, solve time, peak RSS, and phase breakdown — for the full SuiteSparse suite. These baselines serve as the reference point for measuring parallel speedup in Phase 8.2 and for detecting future performance regressions.

The developer runs a systematic baseline collection across all 65 SuiteSparse matrices and records timing and memory data in a structured, reproducible format.

**Why this priority**: Baselines are the measuring stick for Phase 8.2's parallelism work. Without them, there is no way to quantify speedup or detect regressions. However, they can be collected after instrumentation (Story 1) and optimization (Story 2) are complete, so the baselines reflect the optimized sequential solver.

**Independent Test**: Can be tested by running the baseline collection tool and verifying that structured output is produced for all 65 matrices with all required metrics (timing per phase, peak memory).

**Acceptance Scenarios**:

1. **Given** the instrumented and optimized solver, **When** the baseline collection is run on the full SuiteSparse suite, **Then** structured output is produced containing: ordering time, symbolic analysis time, numeric factorization time, solve time, and peak RSS for each matrix.
2. **Given** baseline data for all matrices, **When** the developer examines the results, **Then** workload distribution is visible: which supernodes/fronts consume the most time, what fraction of total time is spent in the top N% of fronts.
3. **Given** baseline data, **When** the developer reviews the results, **Then** there is enough information to make an informed recommendation about whether tree-level parallelism or intra-node BLAS-3 parallelism should be prioritized in Phase 8.2.

---

### User Story 4 — Produce Performance Analysis Report (Priority: P4)

A solver developer produces a summary report from profiling and baseline data that characterizes the solver's workload distribution. The report includes front size histograms, per-supernode time contribution analysis, allocation pressure statistics, and a recommendation for Phase 8.2's parallelism strategy.

**Why this priority**: The report synthesizes all other deliverables into actionable guidance for Phase 8.2. It is the final output but depends on all prior work being complete.

**Independent Test**: Can be tested by verifying that the report is generated and contains all required sections with non-trivial content derived from actual profiling data.

**Acceptance Scenarios**:

1. **Given** profiling data and baselines for the SuiteSparse suite, **When** the performance report is generated, **Then** it contains: front size distribution, per-supernode time contribution histogram, allocation pressure analysis, and a parallelism strategy recommendation.
2. **Given** the performance report, **When** reviewed by the developer, **Then** it clearly answers: "For matrix class X, is tree-level parallelism or intra-node parallelism more beneficial?" for at least three representative matrix classes (small/sparse, medium, large/dense fronts).

---

### Edge Cases

- What happens when a matrix has only one supernode (single dense front)? Profiling should still produce meaningful output showing the single front's breakdown.
- How does profiling overhead affect timing accuracy on very small matrices (sub-millisecond factorization)? Instrumentation overhead should be documented and accounted for.
- What happens when a matrix has zero delayed pivots? The profiling and baseline data should correctly handle the case where no columns are delayed.
- How are matrices with extremely skewed front size distributions handled (e.g., one huge root front and thousands of tiny leaf fronts)? The analysis should capture this distribution shape.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide per-supernode timing instrumentation within the multifrontal factorization loop, measuring at minimum: frontal matrix assembly, dense APTP kernel execution, and contribution block extraction.
- **FR-002**: System MUST generate Chrome Trace JSON output from instrumented factorization runs, compatible with the existing `ProfileSession` / `SectionGuard` API and viewable in standard trace viewers.
- **FR-003**: System MUST identify and reduce the top allocation hotspots in the factorization inner loop — specifically, per-iteration temporary vector allocations, per-block backup matrix allocations, and per-block workspace matrix allocations — by hoisting allocations out of tight loops and reusing buffers across iterations.
- **FR-004**: System MUST collect sequential performance baselines for the full 65-matrix SuiteSparse suite, recording per-matrix: ordering time, symbolic analysis time, numeric factorization time, solve time, and peak RSS.
- **FR-005**: System MUST produce a workload distribution analysis showing: front size distribution across supernodes, time contribution per front size band, and the fraction of total factorization time attributable to the largest N% of fronts.
- **FR-006**: System MUST produce a parallelism strategy recommendation based on profiling data, characterizing whether tree-level parallelism (many independent supernodes) or intra-node BLAS-3 parallelism (few large supernodes) is more impactful for different matrix classes.
- **FR-007**: All optimizations MUST preserve solver correctness — no regressions on any of the 65 SuiteSparse matrices (reconstruction < 1e-12 where applicable, backward error unchanged from Phase 8.1f).
- **FR-008**: Profiling instrumentation MUST NOT be present in production builds — it should be gated behind the existing feature-flag mechanism.
- **FR-009**: System MUST record per-supernode timing data in a programmatically accessible form (extending the existing per-supernode statistics), not only via Chrome Trace files.

### Key Entities

- **ProfileSession**: Existing thread-safe hierarchical timing recorder with Chrome Trace export. Phase 8.1g extends its use to per-supernode granularity within the factorization loop.
- **PerSupernodeStats**: Existing per-front statistics (front size, pivots, delays, max L entry). Extended with per-supernode timing fields for assembly, kernel, and extraction.
- **FactorizationStats**: Existing aggregate statistics. Extended with total timing breakdown across phases.
- **Performance Baseline**: New structured record per matrix capturing timing per phase, peak RSS, and workload distribution metrics.
- **Performance Report**: New summary document synthesizing baselines and profiling data into parallelism strategy guidance for Phase 8.2.

## Assumptions

- The existing `ProfileSession` / `SectionGuard` / Chrome Trace infrastructure is sufficient for per-supernode instrumentation and does not need architectural changes.
- The existing `MemoryTracker` (RSS snapshots) is sufficient for peak memory measurement on the development platform.
- The full 65-matrix SuiteSparse test set is available in the development environment (extracted from archive into `test-data/suitesparse/`).
- "Top allocation hotspots" refers to the most frequently executed allocations in the inner factorization loop, not all allocations in the codebase. The profiling data from Story 1 will confirm which specific allocations to target.
- Performance baselines are collected on a single representative machine and are relative comparisons (not absolute benchmarks). Phase 8.2 will re-baseline on the same machine for fair parallel vs. sequential comparison.
- The `diagnostic` feature flag (already in Cargo.toml) is the appropriate gate for profiling instrumentation in the factorization hot path. This flag is semantically correct for profiling and avoids adding overhead to regular test builds.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Per-supernode profiling data is collected for all 65 SuiteSparse matrices, and the developer can identify which supernodes consume the most factorization time for any given matrix.
- **SC-002**: Allocation hotspot optimizations measurably reduce the number of heap allocations during factorization (validated by before/after comparison on representative matrices).
- **SC-003**: Peak RSS does not increase after allocation optimizations (and ideally decreases) for any SuiteSparse matrix.
- **SC-004**: Sequential performance baselines are recorded for all 65 SuiteSparse matrices with per-phase timing breakdown (ordering, symbolic, numeric, solve) and peak RSS.
- **SC-005**: A workload distribution analysis is produced that shows what fraction of total factorization time is attributable to the top 10% of fronts (by size) for each matrix — this directly answers whether large-front parallelism matters.
- **SC-006**: No correctness regressions: all 65 SuiteSparse matrices pass with reconstruction < 1e-12 and backward error matching or improving on Phase 8.1f results.
- **SC-007**: The performance report contains a clear, data-driven recommendation for Phase 8.2's parallelism strategy (tree-level vs. intra-node) with supporting evidence from at least three matrix classes.
