# Feature Specification: Profiling and Debug Tools

**Feature Branch**: `008-profiling-debug-tools`
**Created**: 2026-02-08
**Status**: Draft
**Input**: User description: "Implement Phase 1.4 in ssids-plan.md — profiling and debug tools for SSIDS solver development"

## Clarifications

### Session 2026-02-08

- Q: Should profiling/debug tools reuse the existing `test-util` feature flag or get a dedicated flag? → A: Reuse `test-util` for all profiling/debug tools, consistent with benchmarking and testing modules. Split out later if needed.
- Q: Should profiling/debug tools integrate with or extend the existing benchmarking module? → A: Separate peer modules with no dependency on `benchmarking/`. Profiling is a general-purpose dev tool; benchmarking is Criterion-specific. Shared utilities (e.g., RSS reading) can be factored into a common location.
- Q: How should the profiler handle repeated sections with the same label (e.g., called in a loop)? → A: Aggregate — same-named sections accumulate total time, call count, min/max/mean per call. Individual invocations preserved in trace export.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Profile Solver Component Execution Times (Priority: P1)

A solver developer is implementing a new phase (e.g., symbolic analysis, numeric factorization) and wants to understand where time is being spent within that component. They wrap sections of their code with profiling instrumentation, run the solver on a test matrix, and receive a structured breakdown of wall-clock time by labeled section. This helps them identify bottlenecks without leaving the Rust development workflow.

**Why this priority**: Time profiling is the most fundamental development tool for a performance-oriented numerical library. Every subsequent solver phase (Phases 2-8) will need this capability to guide optimization decisions. Without it, developers must rely on external tools for every performance question, which disrupts the development cycle.

**Independent Test**: Can be fully tested by instrumenting the existing `MockBenchmarkable` or any standalone function with profiling sections and verifying that timing data is collected and reported correctly.

**Acceptance Scenarios**:

1. **Given** a developer has instrumented a function with profiling sections, **When** the function executes, **Then** the profiler records the wall-clock duration for each labeled section.
2. **Given** profiling data has been collected for a function with nested sections, **When** the developer inspects the results, **Then** they see a hierarchical breakdown showing each section's time and its proportion of the parent section.
3. **Given** profiling data has been collected, **When** the developer exports the results, **Then** a trace format is produced that can be opened in a standard trace viewer for visual analysis.
4. **Given** profiling instrumentation is present in code, **When** profiling is not active (feature disabled), **Then** the instrumentation compiles to no-ops with zero runtime cost.

---

### User Story 2 - Track Memory Usage During Solver Operations (Priority: P2)

A solver developer wants to understand how much memory each part of the solver allocates and whether peak memory usage is acceptable for their target matrix sizes. They annotate solver sections with memory tracking markers, run the solver, and receive a report showing cumulative allocations per labeled section and the high-water mark for the overall process.

**Why this priority**: Sparse solvers are memory-intensive. Understanding allocation patterns early prevents late-stage redesigns. The existing `read_peak_rss_kb` function provides only a single global snapshot; per-section tracking fills a critical gap.

**Independent Test**: Can be fully tested by running memory-tracked operations on hand-constructed matrices and verifying that reported allocation sizes are consistent with expected data structure sizes.

**Acceptance Scenarios**:

1. **Given** a developer has annotated code sections with memory tracking, **When** the sections execute, **Then** the tracker records a memory snapshot (process RSS) at each annotation point.
2. **Given** memory tracking data has been collected, **When** the developer generates a report, **Then** the report shows the peak RSS observed, the delta between each annotated point, and the labels associated with each measurement.
3. **Given** the solver runs on a matrix of known size, **When** memory tracking is active, **Then** the reported memory deltas are plausible relative to the expected data structure sizes (within 2x of theoretical minimum).

---

### User Story 3 - Visualize Sparse Matrix Structure (Priority: P3)

A solver developer or researcher wants to inspect the sparsity pattern of a matrix to verify correctness of reordering, understand fill-in patterns, or debug unexpected behavior. They pass a sparse matrix to a visualization function and receive a textual or image representation of the non-zero structure.

**Why this priority**: Visual inspection of sparsity patterns is invaluable for debugging symbolic analysis (Phase 2-3) and verifying that ordering algorithms produce expected structures. However, the solver can be developed without this — the testing infrastructure already validates correctness numerically.

**Independent Test**: Can be fully tested by generating visualizations for hand-constructed matrices with known sparsity patterns and verifying the output matches expected structure.

**Acceptance Scenarios**:

1. **Given** a sparse matrix loaded from the test suite, **When** the developer generates a sparsity pattern visualization, **Then** a text-based representation is produced showing the non-zero locations in a grid format suitable for terminal display.
2. **Given** a sparse matrix of size larger than the terminal display area, **When** the developer generates a visualization, **Then** the output is downsampled or summarized to fit a configurable maximum dimension while preserving the overall structural pattern.
3. **Given** two matrices (e.g., before and after reordering), **When** the developer generates visualizations for both, **Then** the outputs are directly comparable (same scale, same format) to support visual diffing.

---

### User Story 4 - Visualize Elimination Tree Structure (Priority: P4)

A solver developer working on symbolic analysis (Phase 3) wants to inspect the elimination tree to verify its structure, identify deep chains that could limit parallelism, or debug tree construction. They pass an elimination tree to a visualization function and receive a structured text representation.

**Why this priority**: Elimination tree visualization directly supports Phase 3 (Symbolic Analysis) debugging. It is less urgent than profiling and memory tracking because tree correctness can also be validated through numerical tests.

**Independent Test**: Can be fully tested by constructing elimination trees for small hand-constructed matrices and verifying the output matches the expected parent-child relationships.

**Acceptance Scenarios**:

1. **Given** an elimination tree for a small matrix (n < 20), **When** the developer generates a visualization, **Then** a text-based tree representation is produced showing parent-child relationships.
2. **Given** an elimination tree for a larger matrix, **When** the developer generates a summary, **Then** statistics are produced including tree depth, number of leaves, branching factor distribution, and subtree sizes.

---

### Edge Cases

- What happens when profiling sections are nested more than 10 levels deep? The profiler should support arbitrary nesting without stack overflow or significant overhead.
- What happens when a profiled section panics? Partial timing data collected before the panic should still be accessible for debugging.
- What happens when memory tracking is requested on a platform that doesn't expose RSS (non-Linux)? The system should gracefully degrade, returning `None` or a clear "not supported" indication rather than failing.
- What happens when visualizing an empty matrix (0x0 or all zeros)? The visualization should produce a meaningful "empty" representation rather than crashing.
- What happens when the profiler is used in a multi-threaded context? The profiler MUST support concurrent recording from multiple threads from the start, avoiding a redesign when parallelism is introduced in Phase 8. This requires thread-safe data structures or thread-local storage with post-execution merging.

## Requirements *(mandatory)*

### Functional Requirements

**Profiling (P1):**

- **FR-001**: The system MUST provide a mechanism to record wall-clock duration of labeled code sections with sub-microsecond precision.
- **FR-002**: The system MUST support hierarchical (nested) profiling sections, recording each section's time relative to its parent.
- **FR-003**: The system MUST export profiling data in Chrome Trace Event format (JSON) for viewing in trace visualization tools.
- **FR-004**: The system MUST compile profiling instrumentation to zero-cost no-ops when the profiling feature is not enabled.
- **FR-005**: The system MUST provide a summary report showing each section's total time, call count, min/max/mean per-call duration, and percentage of parent section time.
- **FR-005a**: When a labeled section is entered multiple times, the system MUST aggregate invocations (total time, call count, min/max/mean) in summary reports. Individual invocations MUST be preserved in the trace export.

**Memory Tracking (P2):**

- **FR-006**: The system MUST record process-level memory snapshots (RSS) at developer-annotated points during execution.
- **FR-007**: The system MUST compute and report the peak RSS observed across all snapshots in a tracking session.
- **FR-008**: The system MUST compute memory deltas between consecutive annotated points, attributing changes to labeled operations.
- **FR-009**: The system MUST generate a textual report of memory usage by section, suitable for terminal display.
- **FR-010**: The system MUST gracefully handle platforms where RSS measurement is unavailable, returning an explicit "unavailable" indication rather than failing.

**Sparsity Visualization (P3):**

- **FR-011**: The system MUST produce a text-based sparsity pattern representation for any sparse matrix, showing non-zero positions in a grid format.
- **FR-012**: The system MUST support configurable output dimensions, downsampling large matrices to fit a specified maximum width and height.
- **FR-013**: The system MUST display matrix dimensions and non-zero count alongside the pattern.

**Elimination Tree Visualization (P4):**

- **FR-014**: The system MUST produce a text-based tree representation showing parent-child relationships for elimination trees.
- **FR-015**: The system MUST compute and display summary statistics for an elimination tree: depth, leaf count, and branching factor distribution.

**Cross-Cutting:**

- **FR-016**: All profiling and debug tools MUST be gated behind the `test-util` feature flag so they are excluded from production builds, consistent with the existing benchmarking and testing modules.
- **FR-017**: All profiling and debug tools MUST be usable independently — enabling one tool does not require enabling others.
- **FR-018**: The profiler MUST support concurrent recording from multiple threads without data races or lost measurements, so that it is ready for parallel solver phases without redesign.

### Key Entities

- **Profile Session**: A collection of timed sections from a single execution, with hierarchical parent-child relationships. Attributes: root section, start time, sections list.
- **Profile Section**: A labeled measurement within a profile session. When the same label is entered multiple times (e.g., in a loop), invocations are aggregated. Attributes: name, total duration, call count, min/max/mean per-call duration, parent section, child sections.
- **Memory Snapshot**: A point-in-time measurement of process memory. Attributes: label, timestamp, RSS value, delta from previous snapshot.
- **Memory Report**: An aggregated view of memory snapshots from a tracking session. Attributes: snapshots list, peak RSS, total tracked sections.
- **Sparsity Pattern**: A visual representation of non-zero structure. Attributes: source matrix dimensions, display dimensions, character grid.

## Assumptions

- **Thread-safe profiling from the start**: Although the current solver architecture is simplicial (column-by-column), the profiler will support concurrent recording from multiple threads from day one. This avoids a disruptive redesign when parallelism is introduced in Phase 8.
- **Text-based visualization is sufficient**: Graphical (image/SVG) output is not needed at this stage. Text-based output works in terminals, CI logs, and documentation. Image export can be added later if needed.
- **Separate from benchmarking module**: Profiling and debug tools live in their own peer modules alongside `benchmarking/` and `testing/`, with no coupling between them. Shared utilities (e.g., RSS reading) may be factored into a common location accessible to both.
- **RSS is the primary memory metric**: Tracking individual allocations via a global allocator wrapper is a more invasive approach. Process RSS snapshots provide useful coarse-grained data with minimal complexity.
- **Chrome Trace Event format for profiling export**: This is a widely supported format viewable in Chrome's `chrome://tracing`, Perfetto, and other tools. It provides hierarchical visualization with no additional dependencies.
- **Feature gating reuses `test-util`**: All profiling and debug tools are gated behind the existing `test-util` feature flag, consistent with the benchmarking and testing modules. No dedicated feature flag is introduced at this stage. A separate flag can be split out later if profiling-in-production becomes a real need.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A developer can profile any function's internal sections and receive a timing breakdown within 5 minutes of adding instrumentation, with no external tools required.
- **SC-002**: Profiling instrumentation adds less than 1% overhead to measured operations when active, and zero overhead when the feature is disabled.
- **SC-003**: Memory tracking reports correctly identify the peak-RSS section to within 10% accuracy on at least 90% of test runs (accounting for OS measurement granularity).
- **SC-004**: Sparsity pattern visualizations for all 15 hand-constructed test matrices produce correct representations that match their known structures.
- **SC-005**: Exported Chrome Trace files load successfully in at least one standard trace viewer and display the correct section hierarchy and timing data.
- **SC-006**: All profiling and debug tools are fully excluded from production builds (verified by checking that the release binary size is unchanged).
