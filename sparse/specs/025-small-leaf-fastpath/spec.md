# Feature Specification: Small Leaf Subtree Fast Path

**Feature Branch**: `025-small-leaf-fastpath`
**Created**: 2026-02-24
**Status**: Draft
**Input**: User description: "Phase 9.1f: Small leaf subtree fast path for simplicial matrices"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Faster Factorization on Simplicial Matrices (Priority: P1)

A developer solving a sparse symmetric indefinite system whose sparsity structure produces many small supernodes (e.g., 3D finite-difference grids, graph Laplacians) expects factorization performance competitive with SPRAL. Currently these "simplicial" matrices are 1.5-3.3x slower than SPRAL because every supernode — even tiny ones with front sizes under 100 — pays the full cost of frontal matrix allocation, zeroing, global-to-local mapping, contribution block extraction, and extend-add assembly. The solver should detect when a subtree of the elimination tree consists entirely of small supernodes and factor it using a streamlined code path that avoids the per-supernode overhead of the general multifrontal machinery.

**Why this priority**: This is the core optimization. Simplicial matrices (dixmaanl, mario001, bloweybq, linverse, spmsrtls, rail_79841) are the worst-performing category in the benchmark suite relative to SPRAL. The full frontal matrix machinery is disproportionately expensive when supernodes are tiny, and profiling shows significant unaccounted overhead (23.3% on dixmaanl) attributable to per-supernode setup cost on thousands of small fronts.

**Independent Test**: Can be tested by factorizing known simplicial SuiteSparse matrices (dixmaanl, bloweybq, mario001) and verifying that (a) backward error remains below the established 5e-11 threshold and (b) factorization time decreases relative to the current implementation.

**Acceptance Scenarios**:

1. **Given** a sparse symmetric indefinite matrix whose elimination tree produces leaf subtrees with all front sizes below the classification threshold, **When** the solver factorizes the matrix, **Then** it uses the small-leaf fast path for those subtrees and produces a numerically correct result (backward error < 5e-11), equivalent to the general path within floating-point tolerance.
2. **Given** a matrix with a mix of small leaf subtrees and large supernodes (e.g., a 3D FEM matrix), **When** the solver factorizes, **Then** only the qualifying subtrees use the fast path; large supernodes use the existing general path. The overall factorization result is correct (backward error < 5e-11).
3. **Given** the simplicial benchmark matrices (dixmaanl, bloweybq, mario001, linverse, spmsrtls), **When** factorized with the small-leaf fast path enabled, **Then** factorization time is at most 1.5x SPRAL's factorization time (currently 2-3.3x).

---

### User Story 2 - Subtree Classification (Priority: P2)

At the start of factorization, after amalgamation produces the final supernodal structure, the solver should classify which subtrees of the elimination tree qualify as "small leaf subtrees" so that the factorization can dispatch to the appropriate code path. This classification uses only symbolic information (front sizes derived from the post-amalgamation supernode structure) and runs as a cheap O(n_supernodes) pass. Because amalgamation runs per-factorization, the classification also runs per-factorization; caching it at the symbolic analysis level is a possible future optimization but not required for the initial implementation.

**Why this priority**: Subtree classification is a prerequisite for the fast path dispatch. The classification pass is lightweight (O(n_supernodes), no allocations beyond the output vector) and its cost is negligible relative to the factorization itself.

**Independent Test**: Can be tested by running symbolic analysis on known matrices and verifying the classification output: which supernodes are marked as belonging to small leaf subtrees, and that the classification is consistent with the front size threshold.

**Acceptance Scenarios**:

1. **Given** a matrix where all supernodes have front sizes below the threshold, **When** classification completes, **Then** the entire tree (or all leaf-to-root paths that qualify) is classified as small-leaf subtrees.
2. **Given** a matrix with a deep elimination tree where leaves are small but the root subtree has large fronts, **When** classification completes, **Then** only the leaf portions below the first large supernode are classified as small-leaf subtrees.
3. **Given** a matrix that is refactorized with different numeric values but the same sparsity, **When** the solver runs classification again, **Then** it produces the same subtree assignments (classification is deterministic for a given supernodal structure). Note: caching classification across refactorizations is a future optimization; the initial implementation re-runs the O(n_supernodes) pass each time.

---

### User Story 3 - No Regression on Non-Simplicial Matrices (Priority: P2)

For matrices where most supernodes are large (bulk FEM, large supernodal structures), the small-leaf fast path should not introduce any performance regression. These matrices either have no qualifying subtrees (so the fast path is never invoked) or have only a few trivial leaf subtrees that contribute negligible total work.

**Why this priority**: The existing benchmark suite includes 65 SuiteSparse matrices spanning diverse sparsity structures. Any optimization must not degrade performance on the matrices that are already performing well.

**Independent Test**: Can be tested by running the full SuiteSparse benchmark suite and comparing factorization times before and after the change. No matrix should show a statistically significant regression (beyond measurement noise of ~5%).

**Acceptance Scenarios**:

1. **Given** a large supernodal matrix (e.g., c-71, c-big, G3_circuit), **When** factorized with the small-leaf fast path enabled, **Then** factorization time is within 5% of the time without the fast path (noise tolerance).
2. **Given** the full 65-matrix SuiteSparse suite, **When** factorized, **Then** all 65 matrices produce backward error below 5e-11 (no correctness regression) and the median SPRAL ratio does not increase.

---

### Edge Cases

- What happens when a subtree has exactly one supernode at the threshold boundary (front size equals the threshold)? The classification should use a strict less-than comparison — front size < threshold qualifies, front size >= threshold does not.
- What happens when delayed pivots from a small-leaf supernode cascade into a parent outside the small-leaf subtree? The fast path must produce a contribution block in the same format as the general path so the parent's extend-add works unchanged.
- What happens when all children of a supernode qualify as small-leaf but the supernode itself does not (front size above threshold)? Only the children use the fast path. The parent assembles their contributions via the normal extend-add.
- What happens when a small-leaf subtree contains supernodes with zero original matrix entries (pure fill-in from amalgamation)? The fast path must still handle supernodes with no original entries correctly.
- How does the fast path interact with MC64 scaling? Scaling factors must be applied identically to the general path (multiply original entries by `scaling[row] * scaling[col]` during assembly).
- What happens when the small-leaf threshold is set to 0 (effectively disabling the fast path)? The solver should fall back to the general path for all supernodes with no behavioral change.

## Clarifications

### Session 2026-02-24

- Q: Should single-supernode subtrees qualify for the small-leaf fast path? → A: Exclude single-node subtrees (match SPRAL). Require at least 2 supernodes per subtree.
- Q: Does the solve path (aptp_solve) need modification for small-leaf supernodes? → A: No. The fast path produces identical FrontFactors per supernode; aptp_solve works unchanged.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The solver MUST classify subtrees of the post-amalgamation assembly tree as "small leaf" when every supernode in the subtree has a front size (number of fully-summed columns plus number of pattern rows) below a configurable threshold, and the subtree contains at least 2 supernodes. Single-supernode subtrees are excluded.
- **FR-002**: Subtree classification MUST be performed using symbolic information (post-amalgamation supernode structure). The classification runs as a cheap O(n_supernodes) pass at the start of factorization, after amalgamation. Caching the classification across refactorizations with the same sparsity pattern is a future optimization, not required initially.
- **FR-003**: The solver MUST provide a streamlined factorization path for small-leaf subtrees that reduces per-supernode overhead by amortizing workspace allocation across the entire subtree (single cache-resident buffer reused for all nodes), eliminating workspace pool Mutex contention, and processing nodes sequentially without parallel dispatch overhead.
- **FR-004**: The small-leaf fast path MUST produce numerically correct results — backward error within the established 5e-11 threshold, consistent pivot decisions (same APTP threshold and fallback strategy as the general path).
- **FR-005**: The small-leaf fast path MUST produce contribution blocks at the subtree's root that are compatible with the general path's extend-add assembly, so that parent supernodes outside the small-leaf subtree can assemble them without modification.
- **FR-006**: The solver MUST correctly handle the boundary between small-leaf subtrees and the general path — the contribution block from the topmost supernode in a small-leaf subtree flows into the parent supernode's normal extend-add assembly.
- **FR-007**: The small-leaf classification threshold MUST be a front-size-based parameter: a subtree qualifies when every supernode has front_size < threshold. The default threshold MUST be 256 (matching the existing intra-node parallelism threshold). The threshold MUST be configurable via solver options, and setting it to 0 disables the fast path.
- **FR-008**: The solver MUST preserve support for the `diagnostic` feature — per-supernode timing statistics must still be collected for supernodes processed by the small-leaf fast path, enabling profiling and performance analysis.
- **FR-009**: The fast path MUST correctly apply MC64 scaling factors during assembly of original matrix entries, identically to the general path.
- **FR-010**: The fast path MUST handle pivot delays correctly — if a column cannot be eliminated at a small-leaf supernode, the delayed column must propagate to the parent supernode's fully-summed set, whether the parent is in the small-leaf subtree or the general path.
- **FR-011**: The fast path MUST work correctly with both the sequential and parallel factorization paths. In the parallel case, independent small-leaf subtrees may be processed concurrently.

### Key Entities

- **Small Leaf Subtree**: A maximal contiguous subtree of the post-amalgamation assembly tree where every supernode has a front size below the classification threshold and the subtree contains at least 2 supernodes. Single-supernode subtrees are excluded (matching SPRAL's behavior) because a lone small supernode has minimal extend-add overhead and doesn't benefit from the streamlined path. The subtree is "leaf" in the sense that it has no children outside itself in the partitioned tree — it sits at the bottom of the assembly tree. Multiple independent small-leaf subtrees may exist in a single matrix.
- **Classification Threshold**: A configurable front-size parameter (default 256). A supernode qualifies when its front size (fully-summed columns + pattern rows) is strictly less than the threshold. Applied per-supernode to the post-amalgamation supernode structure. SPRAL uses a flops-based threshold (4M flops per subtree) which bounds total work; the simpler front-size approach was chosen because our target matrices are fully simplicial (entire tree qualifies) and a flops cap can be added later as a parallelism tuning refinement if needed.
- **Contribution Block**: The Schur complement output from a supernode (or from the root of a small-leaf subtree), containing un-eliminated rows and delayed columns. Must be produced in a format compatible with the general path's extend-add assembly regardless of which code path produced it.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Simplicial benchmark matrices (dixmaanl, bloweybq, mario001, linverse, spmsrtls) achieve factorization time at most 1.5x SPRAL's factorization time (currently 2-3.3x).
- **SC-002**: All 65 SuiteSparse matrices produce backward error below 5e-11 (no correctness regression from the fast path).
- **SC-003**: Non-simplicial matrices (c-71, c-big, G3_circuit, and other large supernodal matrices) show no statistically significant factorization time regression (within 5% measurement noise).
- **SC-004**: The full SuiteSparse suite median SPRAL ratio does not increase (currently 0.98x; must remain at or below 1.0x).

## Scope Boundary

- **In scope**: Subtree classification, streamlined factorization path for small-leaf subtrees, contribution block compatibility with general path.
- **Out of scope**: Solve-path changes. The fast path produces the same `FrontFactors` per supernode as the general path, so `aptp_solve` requires no modification.

## Assumptions

- The current amalgamation pass (Phase 9.1a, nemin=32) has already reduced trivial 1-column supernodes, so the small-leaf fast path operates on post-amalgamation supernodes with front sizes typically in the range 32-256.
- The primary source of overhead on simplicial matrices is per-supernode setup cost (frontal matrix allocation, zeroing, g2l mapping) rather than the dense APTP kernel itself, as shown by Phase 9.1c profiling (23.3% unaccounted time on dixmaanl, first supernode alone 18.3% of total).
- The existing workspace reuse infrastructure (FactorizationWorkspace, Phase 9.1b) provides the foundation — the optimization is about bypassing the full frontal-matrix round-trip for qualifying supernodes, not about better workspace management.
- SPRAL's SmallLeafNumericSubtree approach (sequential processing, direct scatter assembly, GEMM-based contribution) is the primary algorithmic reference. The implementation should follow the same high-level strategy while using idiomatic Rust and faer APIs.
- Phase 9.1e (direct GEMM into contribution buffer) may or may not be complete before this feature. The spec is written to be independent of 9.1e's status — the fast path is a separate code path that does not depend on optimizations to the general path.

## Dependencies

- **Phase 9.1a (supernode amalgamation)**: Required — the fast path operates on post-amalgamation supernodes.
- **Phase 9.1b (workspace reuse)**: Required — provides FactorizationWorkspace infrastructure.
- **Phase 9.1e (direct GEMM into contribution)**: Independent — may improve the general path but does not affect the fast path design.

## Algorithm References

The following academic and reference materials inform this feature:

### Academic Papers

- **Duff, Hogg, & Lopez (2020)**: "A New Sparse Symmetric Indefinite Direct Solver" — describes SSIDS's leaf/root subtree partitioning strategy (Algorithm 6.1: `find_subtree_partition`), NUMA-aware scheduling. Reference markdown: `/workspace/rivrs-linalg/references/ssids/duff2020.md` (lines 291-320).
- **Hogg, Duff, & Lopez (2016)**: "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting" — discusses small frontal matrices at leaf nodes and two-level parallelism (tree parallelism near leaves, node parallelism near root). Reference markdown: `/workspace/rivrs-linalg/references/ssids/hogg2016.md` (lines 52, 80-95, 160+).
- **Davis (2016)**: "Direct Methods for Sparse Linear Systems" — Section 11 covers multifrontal method, amalgamation, and tree/node parallelism trade-offs for leaf vs root portions of the elimination tree. Reference markdown: `/workspace/rivrs-linalg/references/ssids/davis2016.md` (lines 1696-1800+).
- **Liu (1992)**: "The Multifrontal Method for Sparse Matrix Solution: Theory and Practice" — foundational reference for multifrontal assembly tree structure, contribution block flow, and tree parallelism. Reference markdown: `/workspace/rivrs-linalg/references/ssids/liu1992.md`.

### SPRAL Source References (BSD-3, freely consultable)

- `SmallLeafSymbolicSubtree.hxx:18-109` — threshold-based leaf subtree identification, symbolic setup, flops estimation bottom-up
- `SmallLeafNumericSubtree.hxx:38-220` (PD path), `187-446` (indefinite path) — sequential factorization without full frontal matrices, direct scatter assembly, GEMM-based contribution computation (`calcLD` + `host_gemm`)
- `NumericSubtree.hxx:53-89, 96-155` — threshold selection (`small_subtree_threshold = 4*10^6 flops`), subtree classification, task scheduling that skips already-handled small-leaf nodes
- `SymbolicNode.hxx` — `insmallleaf` flag on each symbolic node
- `datatypes.f90:95` — default threshold constant (`small_subtree_threshold`)
- `SymbolicSubtree.hxx:57-84` — bottom-up flops estimation with parent penalty to prevent overly broad subtrees
