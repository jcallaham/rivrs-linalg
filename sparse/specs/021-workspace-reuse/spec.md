# Feature Specification: Workspace Reuse & Per-Supernode Allocation Optimization

**Feature Branch**: `021-workspace-reuse`
**Created**: 2026-02-22
**Status**: Draft
**Input**: User description: "Phase 9.1b — Workspace reuse and per-supernode allocation overhead optimization"

## Context & Motivation

Phase 9.1a (supernode amalgamation) successfully reduced supernode counts — c-71 went from 35,372 to 6,350 supernodes (fewer than SPRAL's 7,697). However, c-71 remains 5.5x slower than SPRAL and c-big is 6.3x slower. Post-amalgamation profiling identified per-supernode allocation overhead as the dominant remaining cost: approximately 16-18 heap allocations per supernode through the full frontal matrix machinery (frontal matrix, MixedDiagonal, L11/D11/L21 extraction, contribution blocks).

This overhead affects all matrix categories:
- **Supernodal (c-71, c-big)**: 5.5-6.3x slower despite equal/fewer supernodes than SPRAL
- **Simplicial (bloweybq, dixmaanl, linverse)**: 2.0-2.5x slower; amalgamation helped (was 3.3-7.8x) but remaining supernodes still traverse full frontal machinery
- **Bulk FEM (most matrices)**: 1.0-1.2x range; small per-supernode overhead is dwarfed by dense kernel time

SPRAL addresses this through stack allocation for factors, a buddy system allocator for contribution blocks and thread workspaces, and small-leaf subtree specialization that avoids allocation within L2-cache-resident subtrees.

### Algorithm References

- **Duff, Hogg, Lopez (2020)** — "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting", Section on memory management describing stack allocation for factors and buddy system for workspaces. `/workspace/rivrs-linalg/references/ssids/duff2020.md` (lines 246-289)
- **Liu (1992)** — "The Multifrontal Method for Sparse Matrix Solution", Section 7 on stack storage management and Section 8 on out-of-core/limited-memory variants. `/workspace/rivrs-linalg/references/ssids/liu1992.md` (lines 718-794, esp. 814)
- **Davis (2016)** — survey on direct methods, sequential memory layout and contribution stack. `/workspace/rivrs-linalg/references/ssids/davis2016.md`

### SPRAL Source References

- `SmallLeafSymbolicSubtree.hxx:18-109` — leaf subtree symbolic setup (L2-cache-resident)
- `SmallLeafNumericSubtree.hxx:38-220` — sequential factorization without full frontal matrices
- `BlockPool.hxx:16-82` — fixed-size block pool for workspace reuse
- `BuddyAllocator.hxx:18-148` — buddy allocator for variable-size workspace
- `NumericSubtree.hxx:75-81` — per-thread workspace initialization

### Post-Amalgamation Baseline (Sequential, Single Thread)

| Metric | Value |
|--------|-------|
| Median ratio vs SPRAL | ~1.13x |
| c-71 | 5.48x (was 29.75x pre-amalgamation) |
| c-big | 6.29x (was 11.11x pre-amalgamation) |
| Matrices beating SPRAL | 10 of 65 (ratio < 1.0x) |
| Worst simplicial | dixmaanl 2.51x, rail_79841 2.24x, mario001 2.13x |
| G3_circuit (1.6M) | 1.51x |

### Current Per-Supernode Allocation Inventory

During `factor_single_supernode()` each supernode incurs:

| Component | Allocations | Dominant Size |
|-----------|-------------|---------------|
| Frontal assembly | 1 Mat + 2-3 Vec | O(m^2) |
| APTP kernel (col_order, pivot_log, d, delayed_cols) | ~7 Vec | O(m) |
| Factor extraction (L11, D11, L21, row_indices) | 2 Mat + 4 Vec | O(ne^2) + O(r x ne) |
| Contribution extraction (data, row_indices) | 1 Mat + 1 Vec | O((m-ne)^2) |

Total: ~4 dense matrices + ~14 vectors per supernode, all heap-allocated and freed.

The solve path already uses workspace reuse: two work buffers allocated once and reused across all supernodes.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Workspace Reuse for Frontal Matrix Factorization (Priority: P1)

A solver user factorizes a sparse symmetric indefinite matrix. The factorization loop reuses pre-allocated workspace buffers across supernodes rather than allocating and freeing dense matrices and vectors for each supernode independently. This eliminates the ~16-18 heap allocations per supernode that currently dominate factorization cost for matrices with many supernodes.

**Why this priority**: This is the single highest-leverage optimization. It affects every matrix (not just outliers) and is identified as the likely main reason c-71/c-big remain 5-6x slower than SPRAL despite comparable supernode counts. SPRAL pre-allocates workspace once and reuses it; rivrs currently allocates fresh on every supernode.

**Independent Test**: Can be tested by running the full SuiteSparse benchmark suite and comparing factorization times against pre-optimization baselines. The existing backward error and reconstruction correctness tests validate that workspace reuse does not introduce numerical errors.

**Acceptance Scenarios**:

1. **Given** the solver factorizes any matrix from the 65-matrix SuiteSparse suite, **When** factorization completes, **Then** backward error remains below 5e-11 (identical to pre-optimization behavior)
2. **Given** the solver factorizes c-71 (6,350 supernodes post-amalgamation), **When** factorization completes, **Then** total factorization time is measurably reduced vs the pre-optimization baseline
3. **Given** the solver factorizes a simplicial-dominated matrix (bloweybq, dixmaanl), **When** factorization completes, **Then** factorization time is measurably reduced vs the pre-optimization baseline
4. **Given** the solver factorizes the same matrix twice (refactorization), **When** both factorizations complete, **Then** both produce results within the established tolerance thresholds (workspace reuse is transparent)

---

### User Story 2 - Contribution Block Copy Optimization (Priority: P2)

Instead of extracting the Schur complement via an element-by-element copy into a separately allocated dense matrix after each supernode's factorization, the solver uses bulk column-wise copy operations (e.g., faer's `copy_from` or `submatrix().to_owned()`) to reduce the overhead of contribution block extraction. The contribution block must still own its data (views into the frontal workspace are infeasible because the workspace is reused for the next supernode before the parent processes the contribution).

**Why this priority**: Contribution block extraction is the second-largest allocation category per supernode (after the frontal matrix itself). Eliminating the copy reduces both allocation overhead and memory bandwidth consumption. However, this requires careful lifetime management since the frontal workspace is reused for the next supernode.

**Independent Test**: Can be tested by factorizing matrices with large contribution blocks (FEM matrices with high fill-in) and verifying that numerical results are identical while elapsed time is reduced. The existing `extract_contribution` test cases cover correctness.

**Acceptance Scenarios**:

1. **Given** a matrix with large Schur complements (e.g., sparsine with max_front 11,125), **When** factorization completes, **Then** backward error is identical to pre-optimization baseline
2. **Given** any matrix in the SuiteSparse suite, **When** factorization completes, **Then** peak memory usage does not increase relative to the pre-optimization baseline
3. **Given** a matrix with many children per parent supernode, **When** the extend-add assembly accumulates contributions from multiple children, **Then** the numerical result matches the current separate-allocation approach

---

### Deferred: Simplicial Fast Path

A simplicial fast path (column-by-column LDL^T on CSC storage, bypassing frontal matrix machinery for single-column supernodes) was considered as P3 but is **deferred to post-release**. Rationale: amalgamation already merges most single-column supernodes, US1 workspace reuse addresses remaining per-supernode overhead, and adding a separate factorization code path increases maintenance and test surface for diminishing returns. Phase 9.1c profiling will determine if this is still needed.

---

### Edge Cases

- What happens when the maximum front size is 1 (purely diagonal matrix)? Workspace pre-allocation should handle degenerate dimensions gracefully.
- What happens when a supernode has zero delayed columns? The workspace should still work correctly with the "fast" case of no delays.
- What happens when all children of a parent supernode have large contribution blocks that overlap? In-place contribution handling must accumulate multiple overlapping contributions correctly.
- What happens during parallel factorization? Each rayon worker has its own thread-local workspace (Cell-based move semantics, matching the existing g2l pattern). Workspace must be sized to the max front in the worker's assigned subtree.
- What happens when workspace reuse is combined with amalgamated supernodes that have non-contiguous `owned_ranges`? The `scatter_original_entries_multi` path must work correctly with reused buffers.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The factorization loop MUST pre-allocate reusable workspace buffers (frontal matrix, APTP kernel scratch, extraction buffers) sized to the maximum front dimension known from symbolic analysis, and reuse them across all supernodes in a sequential factorization pass
- **FR-002**: Workspace reuse MUST produce numerically equivalent results to the current per-supernode allocation approach, defined as: backward error remains below 5e-11 and reconstruction error remains below 1e-12 on all test matrices. Bit-for-bit identity is not required; floating-point operation reordering (e.g., from changed contribution block accumulation order) is acceptable provided tolerances are met
- **FR-003**: Workspace reuse MUST be compatible with the existing parallel factorization architecture (rayon tree-level + faer intra-node BLAS). Each rayon worker MUST have its own thread-local workspace, following the existing `Cell`-based move semantics pattern used for `g2l` buffers in tree-level factorization
- **FR-004**: The contribution block handling MUST avoid unnecessary dense matrix copies where the trailing submatrix can be retained or transferred without allocation
- **FR-005**: The solver MUST NOT increase peak memory usage relative to the current implementation on any matrix in the test suite
- **FR-006**: ~~Deferred~~ — Simplicial fast path removed from scope (see Deferred section under User Scenarios). Workspace reuse (FR-001) addresses simplicial overhead; separate code path deferred to post-release
- **FR-007**: All 65 SuiteSparse matrices MUST continue to pass the strict backward error threshold (< 5e-11) after optimization
- **FR-008**: The 15 hand-constructed matrices MUST continue to pass reconstruction tests (||P^T A P - L D L^T|| / ||A|| < 1e-12)
- **FR-009**: Workspace buffers MUST be properly zeroed or overwritten before each supernode's use to prevent data leakage between supernodes
- **FR-010**: The `diagnostic` feature timing instrumentation MUST continue to work correctly with workspace reuse, reporting accurate per-supernode timing

### Key Entities

- **Factorization Workspace**: A pre-allocated set of dense matrix buffers and vector scratch space, sized to the maximum front dimensions from symbolic analysis, reused across supernodes during factorization
- **Contribution Block**: The Schur complement (trailing submatrix) produced by partial elimination of a frontal matrix, passed from child to parent supernodes during multifrontal assembly. Currently a separately allocated `Mat<f64>` + `Vec<usize>`; targeted for in-place handling
- **Simplicial Supernode**: *(deferred — see Deferred section)* A supernode with exactly one column; fast-path handling removed from this feature's scope

## Clarifications

### Session 2026-02-22

- Q: Should numerical equivalence after workspace reuse be bit-for-bit identical or toleranced? → A: Toleranced equivalence — backward error < 5e-11 and reconstruction error < 1e-12 are sufficient. Bit-for-bit identity is not required.
- Q: Should the simplicial fast path (US3) be included in this feature's scope or deferred? → A: Deferred to post-release. US1 workspace reuse addresses the overhead; a separate code path adds maintenance burden for diminishing returns.
- Q: How should per-worker workspaces be managed for parallel factorization? → A: Thread-local, following the existing Cell-based move semantics pattern used for g2l buffers in factor_tree_levelset.

## Assumptions

- The maximum front size is known after symbolic analysis and does not change during factorization (true by design: `AptpSymbolic` computes this)
- Pre-allocating workspace to `max_front_size` dimensions is acceptable in terms of memory; for most matrices, the largest front uses comparable memory to the total factor storage, so the workspace overhead is modest
- Zeroing a pre-allocated workspace before each supernode is cheaper than allocating + freeing, because the OS zero-page cache is not exhausted and memory remains cache-warm
- The parallel factorization architecture (Phase 8.2) already uses thread-local `g2l` buffers via `Cell`-based move semantics; workspace reuse should follow the same pattern
- faer's `Mat` type supports efficient in-place submatrix operations without requiring owned allocations for subviews

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: c-71 and c-big factorization time within 2x of SPRAL (currently 5.48x and 6.29x respectively)
- **SC-002**: Simplicial matrices (bloweybq, dixmaanl, mario001, etc.) within 1.5x of SPRAL (currently 2.0-2.5x)
- **SC-003**: Full SuiteSparse suite median factor time ratio at or below 1.0x vs SPRAL (currently ~1.13x)
- **SC-004**: Zero correctness regressions — all 65 SuiteSparse matrices maintain backward error < 5e-11
- **SC-005**: Zero performance regressions — no matrix in the suite becomes slower after optimization
- **SC-006**: Peak memory usage does not increase on any matrix in the test suite
