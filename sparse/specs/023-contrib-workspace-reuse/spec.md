# Feature Specification: Contribution Workspace Reuse

**Feature Branch**: `023-contrib-workspace-reuse`
**Created**: 2026-02-23
**Status**: Draft
**Input**: User description: "Phase 9.1d: Contribution workspace reuse — eliminate per-supernode allocation for contribution blocks (both Approaches 1 & 2)"

## User Scenarios & Testing

### User Story 1 — Eliminate contribution block allocation churn (Priority: P1)

A solver user factorizes a large sparse symmetric indefinite matrix (e.g., optimization matrices with deep assembly trees and large fronts like c-71 or c-big) and expects factorization to complete without spending a dominant fraction of time on memory allocation and deallocation rather than numerical computation.

**Why this priority**: Phase 9.1c profiling identified contribution block allocation as the single largest remaining bottleneck — `extract_contribution` accounts for 40.1% and `extend_add` for 33.3% of factorization time on c-71. Diagnostic testing (R8) confirmed the dominant cost is **first-touch page faults** during faer's `Mat::zeros` zero-write pass on freshly mapped memory (~15-25% of factor time), compounded by the column-by-column data copy (~10-15%). Buffer reuse eliminates page faults by keeping pages physically resident (empirically measured 1.7-2.3x speedup on extraction). Direct extend-add eliminates both the zero-write and the data copy entirely for last-child supernodes. Combined, these optimizations target a ~30-50% reduction in factor time.

**Independent Test**: Can be validated by factorizing c-71 and measuring (a) factor time relative to SPRAL and (b) system time percentage. Success is visible without any other optimization.

**Acceptance Scenarios**:

1. **Given** a matrix with many large supernodes (c-71: 6,350 supernodes post-amalgamation, max front ~2,475), **When** factorizing on the sequential path, **Then** the solver does not allocate a new `Mat<f64>` for every supernode's contribution block — instead reusing pre-allocated workspace or reading directly from the frontal workspace.
2. **Given** the same matrix, **When** measuring system-level overhead, **Then** kernel (`sys`) time drops from ~32% of wall time to under 10%, confirming that first-touch page fault overhead from per-supernode allocation has been eliminated.
3. **Given** the full 65-matrix SuiteSparse benchmark suite, **When** factorizing each matrix, **Then** backward error remains within established tolerances (< 5e-11 for MatchOrderMetis) and no matrix regresses in correctness.

---

### User Story 2 — Direct extend-add from frontal workspace (Priority: P1)

When a supernode has been factored and its Schur complement (trailing submatrix) is still resident in the frontal workspace, the solver reads that data directly for extend-add into the parent's frontal matrix rather than first copying it into an intermediate contribution block.

**Why this priority**: The contribution block copy is redundant in the sequential path — the data already exists in the frontal workspace. Skipping the copy eliminates the 40.1% extract_contribution cost entirely and reduces the extend-add cost by removing an intermediate data movement step. This is the highest-impact single optimization identified by Phase 9.1c profiling.

**Independent Test**: Can be validated by running c-71 with only the sequential factorization path and comparing sub-phase timing: `extract_contrib` time should drop to near zero for supernodes whose contributions are read directly from the frontal workspace.

**Acceptance Scenarios**:

1. **Given** a supernode that has just been factored (sequential path), **When** the parent's frontal matrix is already allocated in a separate buffer, **Then** the extend-add reads directly from the trailing `(m - ne) x (m - ne)` submatrix of the child's frontal workspace without performing a separate extraction/copy step.
2. **Given** a supernode with delayed columns (APTP pivoting decided some columns could not be eliminated), **When** performing direct extend-add, **Then** the row index mapping correctly accounts for the APTP permutation and delayed column positions so that each Schur complement entry lands in the correct position in the parent's frontal matrix.
3. **Given** a parent supernode with multiple children, **When** not all children can use direct extend-add (because only the most recent child's data is still in the frontal workspace), **Then** earlier children's contributions are handled via the pre-allocated contribution workspace (US1) and only the last child processed before the parent uses direct extend-add.

---

### User Story 3 — Parallel path contribution workspace (Priority: P2)

When running with parallel factorization (rayon tree-level scheduling), a per-thread contribution workspace is available so that each rayon worker avoids per-supernode allocation for contribution blocks within its own subtree processing.

**Why this priority**: The parallel path already uses thread-local `FactorizationWorkspace` via `Cell` semantics. Adding a contribution workspace field to this structure extends the same allocation-reuse pattern to contribution blocks. The parallel path is the production-recommended mode for large matrices. However, direct extend-add (US2) is more complex in the parallel case because contributions must cross thread boundaries — so owned `ContributionBlock` remains necessary when child and parent are on different threads. The per-thread contribution workspace still eliminates the per-supernode allocation for the intra-thread case.

**Independent Test**: Can be validated by running c-71 with parallel factorization and verifying factor time improvement and correctness (65/65 SuiteSparse passing).

**Acceptance Scenarios**:

1. **Given** a parallel factorization where a thread processes multiple supernodes in a subtree, **When** a supernode's contribution does not need to cross a thread boundary, **Then** the contribution is written into the thread-local contribution workspace instead of allocating a new `Mat<f64>`.
2. **Given** a parallel factorization where a child's contribution must be transferred to a parent processed on a different thread, **When** the contribution crosses a thread boundary, **Then** an owned `ContributionBlock` is allocated (or the contribution workspace is swapped out) to enable safe transfer.
3. **Given** the full benchmark suite under parallel factorization, **When** comparing backward error against sequential results, **Then** backward errors are bit-exact (identical to sequential path results for the same matrix and ordering).

---

### Edge Cases

- What happens when a supernode is fully eliminated (ne == m, no contribution)? The solver must skip contribution extraction entirely, as it does today — no workspace write is needed.
- What happens when delayed columns change the contribution row indices compared to the symbolic prediction? The row index mapping must be computed from the actual APTP result permutation, not from precomputed assembly maps. Direct extend-add must handle this dynamic mapping.
- What happens when the parent's frontal matrix is larger than the child's? The parent's buffer is sized independently (from its own front size estimate). The extend-add scatter only writes to positions in the parent that correspond to the child's contribution rows.
- What happens when a single-child chain exists (linear chain in the assembly tree)? This is the ideal case for direct extend-add: the child's frontal workspace remains valid until the parent starts, and only one child's contribution needs to be merged. No contribution workspace is needed.
- What happens when workspace capacity estimates are exceeded due to delayed columns? The workspace must be resizable (as `FactorizationWorkspace` already handles via `ensure_capacity`). A contribution that exceeds the pre-allocated workspace dimension must trigger a grow operation.

## Requirements

### Functional Requirements

- **FR-001**: The solver MUST provide a pre-allocated contribution workspace (a dense buffer and row index buffer) as part of the per-thread factorization workspace, sized to accommodate the largest expected contribution block dimension.
- **FR-002**: The `extract_contribution` operation MUST be able to write into the pre-allocated contribution workspace instead of allocating a new `Mat<f64>` for each supernode.
- **FR-003**: On the sequential factorization path, the solver MUST support reading the Schur complement directly from the child's frontal workspace (the trailing `(m - ne) x (m - ne)` submatrix) for extend-add into the parent's frontal matrix, without an intermediate copy.
- **FR-004**: Direct extend-add from the frontal workspace MUST correctly handle delayed columns by computing the contribution row index mapping from the APTP factorization result (permutation and number of eliminated columns).
- **FR-005**: The factorization loop MUST ensure that the parent's frontal matrix buffer is available (allocated and zeroed in the appropriate region) before the child's direct extend-add writes into it. The solver MUST use a dual-buffer strategy: two frontal `Mat<f64>` buffers in the workspace, alternated between parent and child tree levels. The child's frontal matrix occupies one buffer while the parent's frontal matrix occupies the other, enabling direct extend-add from the child's buffer into the parent's buffer without overwriting.
- **FR-006**: For parent supernodes with multiple children, the solver MUST handle the case where earlier children's contributions have been overwritten in the frontal workspace by later siblings. Earlier children MUST use the contribution workspace (FR-001/FR-002) while the last-processed child before the parent MAY use direct extend-add (FR-003).
- **FR-007**: The parallel factorization path MUST provide per-thread contribution workspaces via the existing thread-local mechanism, extending `FactorizationWorkspace` with contribution buffer fields.
- **FR-008**: When a contribution must cross a thread boundary in the parallel path, the solver MUST produce an owned contribution block that can be safely transferred between threads.
- **FR-009**: All existing correctness guarantees MUST be preserved: backward error < 5e-11 for all 65 SuiteSparse matrices, bit-exact reconstruction test results, and correct inertia counts.
- **FR-010**: The contribution workspace MUST be resizable to handle cases where delayed columns cause the actual contribution dimension to exceed the symbolic estimate.

### Key Entities

- **FactorizationWorkspace**: The existing pre-allocated reusable buffer structure (frontal data, row indices, delayed columns buffer, global-to-local mapping). Extended with contribution data and contribution row indices fields.
- **ContributionBlock**: The existing Schur complement type (data + row_indices + num_delayed). Remains as the owned representation for the parallel cross-thread transfer case, but is no longer the primary storage for sequential-path contributions.
- **FrontalMatrix**: The existing borrowed view into the frontal workspace. The trailing submatrix of a factored frontal matrix is the in-place contribution — the new "direct extend-add" reads from this view.
- **AssemblyMaps**: The existing precomputed scatter index structure. May need extension to support direct extend-add row mappings (mapping from frontal workspace positions to parent frontal positions, accounting for the `ne` offset).

## Success Criteria

### Measurable Outcomes

- **SC-001**: c-71 factorization time is within 1.5x of SPRAL factor time (currently 2.48x after Phase 9.1c; target improvement of ~40-70% of current factor time).
- **SC-002**: System (`sys`) time during factorization of c-71 drops from ~32% of wall time to under 10%, confirming elimination of first-touch page fault overhead from per-supernode contribution block allocation.
- **SC-003**: All 65 SuiteSparse benchmark matrices pass with backward error < 5e-11 under MatchOrderMetis ordering (no correctness regressions).
- **SC-004**: No matrix in the benchmark suite regresses in factor time by more than 5% compared to Phase 9.1c baseline.
- **SC-005**: The `extract_contrib` sub-phase timing (diagnostic feature) drops to near zero for supernodes handled via direct extend-add on the sequential path.
- **SC-006**: dTLB-load-misses for c-71 drop measurably compared to Phase 9.1c baseline (644 billion), if measurable in the deployment environment. (R8 diagnostic showed TLB shootdowns are a minor contributor; this metric is secondary to SC-001/SC-002.)

## Algorithm References

The following academic references and their locations in the references directory inform this feature:

| Reference | Key Contribution | File |
|-----------|-----------------|------|
| Liu (1992), "The Multifrontal Method for Sparse Matrix Solution" | Sections 4-6: frontal matrix assembly, extend-add operator, stack-based workspace management for contribution blocks, supernodal multifrontal Algorithm 6.1 | `/workspace/rivrs-linalg/references/ssids/liu1992.md` |
| Duff, Hogg & Lopez (2020), "A New Sparse LDL^T Solver Using APTP" | Section 5: workspace reuse strategy, memory management (stack allocator for factors, buddy allocator for contribution blocks) | `/workspace/rivrs-linalg/references/ssids/duff2020.md` |
| Hogg, Duff & Scott (2016), "A DAG-based parallel solver" | Section 2.3: factorize phase with contribution blocks stored on temporary stack for parent; Section 3.1: sparse assembly using mapping arrays; Fig. 3: assembly of contribution blocks into parent frontal matrix | `/workspace/rivrs-linalg/references/ssids/hogg2016.md` |
| Duff & Reid (1983/1984), "The multifrontal solution of indefinite sparse symmetric linear equations" | Foundational multifrontal method: assembly tree traversal, contribution block propagation to parent | `/workspace/rivrs-linalg/references/ssids/duff1984.md` |

### SPRAL Source References (BSD-3, may be freely consulted)

- `BlockPool.hxx:16-82` — fixed-size block pool for workspace reuse
- `BuddyAllocator.hxx:18-148` — buddy allocator for variable-size workspace (contribution blocks)
- `NumericSubtree.hxx:75-81` — per-thread workspace initialization
- `factor_cpu.cxx:73-168` — factor_subtree with workspace reuse
- `NumericSubtree.hxx:267-340` — assembly path (node_ilist, aval, amap)

## Assumptions

- The sequential factorization path processes supernodes in strict postorder (bottom-up). This means that when processing supernode `s`, all children of `s` have already been processed and their frontal workspaces are no longer needed by any other supernode except `s`'s parent. This invariant is fundamental to the direct extend-add approach.
- The contribution dimension `(m - ne)` is bounded by the maximum front size minus one. The pre-allocated contribution workspace can be sized using the same `estimate_max_front_size` function used for the frontal workspace.
- The existing `AssemblyMaps` precomputation (Phase 9.1c) provides child-to-parent row mappings that can be extended or adapted for the direct extend-add case, avoiding the need for a completely new index mapping scheme.
- The `diagnostic` feature flag will be used to instrument the new code paths with sub-phase timing, continuing the pattern established in Phase 9.1c.
- For the parallel path, the existing `Cell<FactorizationWorkspace>` pattern (take/set move semantics) will be extended to include the contribution workspace fields. No new synchronization primitives are needed.
- Huge pages (`madvise(MADV_HUGEPAGE)`) are deferred to a separate optimization — this feature focuses on eliminating allocations, not on improving TLB behavior of the remaining allocations.
