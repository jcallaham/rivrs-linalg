# Feature Specification: Assembly & Extraction Optimization (Phase 9.1c)

**Feature Branch**: `022-assembly-extraction-opt`
**Created**: 2026-02-22
**Status**: Complete
**Input**: Phase 9.1c from ssids-plan.md — optimize assembly/extraction overhead for matrices with high front sizes and low elimination counts.

## Context

Post-9.1b profiling of c-71 and c-big reveals that the remaining 4x gap vs SPRAL is dominated by assembly (37%) and extraction (35%), not the dense APTP kernel (23%). These matrices produce thousands of supernodes with front sizes of 2000-5300 but only 5-54 eliminated columns per supernode (elim/front ratio under 1%). Every supernode pays O(F^2) assembly and O((F-E)^2) extraction costs that are not amortized by useful factorization work.

Profiling of ncvxqp7 (which beats SPRAL at 0.74x) confirms the kernel-dominant pattern expected for well-structured matrices: assembly 11.7%, kernel 66.3%, extraction 15.6%.

SPRAL handles assembly differently through precomputed assembly maps (`amap`), two-phase assembly (pre-factor and post-factor), allocator-provided zeroing, and immediate contribution block cleanup. These optimizations reduce the per-supernode O(F^2) overhead that dominates our performance on c-71/c-big.

### Profiling Data

| Matrix   | Assembly | Kernel | Extraction | Ratio vs SPRAL |
|----------|----------|--------|------------|----------------|
| c-71     | 37.7%    | 23.0%  | 34.8%      | 4.06x          |
| c-big    | 36.6%    | 22.6%  | 37.0%      | 4.19x          |
| ncvxqp7  | 11.7%    | 66.3%  | 15.6%      | 0.74x          |
| dixmaanl | 5.4%     | 66.9%  | 4.3%       | 2.40x          |

### Algorithm References

- SPRAL assembly kernels: `assemble.hxx` (BSD-3) — `assemble_pre()`, `assemble_post()`, `add_a_block()`, `assemble_expected()`, `assemble_expected_contrib()`
- SPRAL workspace: `Workspace.hxx` (BSD-3) — thread-local reusable buffer
- SPRAL node data: `NumericNode.hxx` (BSD-3) — `alloc_contrib()`, pool allocator
- Hogg, Ovtchinnikov, Scott (2016): "A Sparse Symmetric Indefinite Direct Solver for GPU Architectures" — Sections 3.1, 3.3 on assembly and contribution management
- Reference Markdown files: `/workspace/rivrs-linalg/references/ssids/hogg2016.md`

## Clarifications

### Session 2026-02-22

- Q: Should scatter map memory growth be bounded by a budget, or is unbounded growth acceptable? → A: Unbounded — scatter map memory is proportional to nnz and accepted as part of symbolic analysis storage.

## User Scenarios & Testing

### User Story 1 — Reduce Assembly Overhead via Precomputed Scatter Maps (Priority: P1)

The solver precomputes scatter-map indices during symbolic analysis so that assembly of original matrix entries and child contributions into the frontal matrix avoids per-entry index arithmetic at factorization time. This eliminates the repeated permutation, global-to-local, and upper-triangle deduplication lookups that currently happen for every nonzero entry during assembly.

**Why this priority**: Assembly accounts for 37% of factor time on c-71/c-big. The per-entry cost of the current scatter function involves 5-7 index lookups per nonzero, and the extend-add function involves 3 indirections per contribution entry. Precomputing these mappings once during symbolic analysis (which is O(nnz) and amortized across factorizations with the same sparsity pattern) directly attacks the largest single overhead source.

**Independent Test**: Run `cargo test` (all unit tests pass) + `cargo test -- --ignored --test-threads=1` (65 SuiteSparse matrices, backward error < 5e-11). Run `profile_matrix` on c-71 and verify assembly percentage decreases.

**Acceptance Scenarios**:

1. **Given** a matrix with precomputed symbolic analysis, **When** factorization assembles original entries into a supernode's frontal matrix, **Then** assembly uses a precomputed index map rather than computing index translations per entry.
2. **Given** a child supernode's contribution block, **When** it is merged into the parent frontal matrix via extend-add, **Then** the scatter positions are determined by a precomputed map rather than per-entry global-to-local lookups.
3. **Given** any matrix in the 65-matrix SuiteSparse test suite, **When** factorized with precomputed scatter maps, **Then** backward error is identical to the current implementation (within floating-point tolerance).
4. **Given** the symbolic analysis result is reused for a second factorization with different numerical values, **When** the scatter maps are reused, **Then** the second factorization produces correct results.

---

### User Story 2 — Reduce Extraction Overhead via Bulk Copy Operations (Priority: P2)

The solver uses bulk memory operations (column-slice copies) instead of element-by-element loops when extracting the L factor, contribution block, and front factors from the factored frontal matrix. This reduces the per-supernode extraction cost from individually indexed writes to sequential memory operations that benefit from hardware prefetching.

**Why this priority**: Extraction accounts for 35-37% of factor time on c-71/c-big. The current contribution extraction and front factor extraction functions use nested loops with per-element indexing. Replacing these with column-wise bulk copies reduces overhead from indexing and bounds checking.

**Independent Test**: Run `cargo test` + SuiteSparse ignored tests. Verify backward errors are bit-exact identical to the pre-optimization baseline.

**Acceptance Scenarios**:

1. **Given** a factored frontal matrix, **When** the contribution block is extracted, **Then** the extraction uses contiguous column-slice copies rather than element-by-element indexing.
2. **Given** a factored frontal matrix, **When** L11 and L21 factors are extracted, **Then** column data is copied in bulk rather than element-by-element.
3. **Given** any matrix in the SuiteSparse suite, **When** factorized with bulk extraction, **Then** backward error is identical to the current implementation.

---

### User Story 3 — Optimize Frontal Matrix Zeroing (Priority: P3)

The solver reduces the cost of zeroing the frontal matrix workspace before each supernode by using bulk memory operations or by zeroing only the regions that will actually receive entries, rather than zeroing the full m x m lower triangle element-by-element.

**Why this priority**: Zeroing is part of the assembly phase timing and contributes to the O(F^2) per-supernode cost. For front=2501, the current implementation zeros ~3.1M entries via individual indexed assignments. Bulk zeroing or lazy zeroing could reduce this overhead. This is lower priority than scatter maps and bulk extraction because zeroing is a simpler operation with less per-element overhead.

**Independent Test**: Run `cargo test` + SuiteSparse ignored tests. Verify all backward errors unchanged.

**Acceptance Scenarios**:

1. **Given** a supernode with front size m, **When** the frontal workspace is prepared, **Then** the zeroing operation is faster than element-by-element indexed assignment.
2. **Given** any matrix in the SuiteSparse suite, **When** factorized with optimized zeroing, **Then** backward error is identical to the current implementation.
3. **Given** a supernode that uses only a portion of the workspace, **When** zeroing is applied, **Then** no more memory is zeroed than necessary for correctness.

---

### Edge Cases

- What happens when a supernode has zero contribution (root supernode, ne == m)? Extraction of contribution block is already skipped — scatter maps for this supernode have no extend-add component. No change needed.
- What happens when a supernode has zero children? Assembly only scatters original entries — extend-add scatter maps are empty but must not cause errors.
- What happens when delayed columns from a child change the expected scatter pattern? The scatter map for extend-add maps child row indices to parent front positions. Delayed columns (which are determined at factorization time, not symbolic time) shift which child rows appear in the contribution. The extend-add map must handle variable contribution sizes.
- What happens when amalgamated supernodes have non-contiguous owned column ranges? Scatter maps must encode the multi-range column ownership correctly, including cross-range upper-triangle deduplication logic.
- What happens when MC64 scaling is applied? Scaling values must be applied during scatter (not encoded in the map, since scaling values change per factorization while the map is reused).

## Requirements

### Functional Requirements

- **FR-001**: The solver MUST precompute assembly scatter maps during symbolic analysis that map original matrix entry positions to frontal matrix positions, eliminating per-entry permutation and deduplication lookups during factorization.
- **FR-002**: The solver MUST precompute extend-add scatter maps during symbolic analysis that map child contribution row indices to parent frontal matrix positions.
- **FR-003**: Scatter maps MUST correctly handle non-contiguous owned column ranges from supernode amalgamation.
- **FR-004**: Scatter maps MUST correctly handle upper-triangle deduplication to avoid double-counting symmetric entries.
- **FR-005**: The solver MUST use bulk memory copy operations for extracting L11, L21, and contribution block data from the factored frontal matrix, replacing element-by-element indexed loops.
- **FR-006**: The solver MUST use bulk memory operations for zeroing the frontal matrix workspace, replacing element-by-element indexed zeroing.
- **FR-007**: All optimizations MUST preserve numerical results — backward error must be identical (within floating-point tolerance) to the pre-optimization implementation on all 65 SuiteSparse test matrices.
- **FR-008**: Scatter maps computed during symbolic analysis MUST be reusable across multiple factorizations with the same sparsity pattern (consistent with the analyze-once, factor-many-times API).
- **FR-009**: The `diagnostic` feature timing instrumentation MUST continue to report assembly, kernel, and extraction times correctly after optimization.
- **FR-010**: Optimizations MUST work correctly with both sequential and parallel factorization paths.
- **FR-011**: Extend-add scatter maps MUST handle the variable number of delayed columns from child supernodes. Since delayed column counts are determined at factorization time (not symbolic time), the map for the non-delayed portion is precomputed while the delayed portion is resolved at factorization time.

### Key Entities

- **Assembly Map**: Per-supernode precomputed index pairs mapping original matrix entry positions to frontal matrix positions. Computed once during symbolic analysis. Stored as part of the symbolic analysis result.
- **Extend-Add Map**: Per-child-parent edge precomputed mapping from child contribution row indices to parent frontal matrix positions. Must handle delayed columns dynamically since their count varies at factorization time.
- **FactorizationWorkspace**: Existing pre-allocated workspace (from Phase 9.1b). Zeroing strategy is modified but the workspace structure is unchanged.

## Success Criteria

### Measurable Outcomes

- **SC-001**: c-71 factor time ratio vs SPRAL improves from 4.06x to below 3.0x.
- **SC-002**: c-big factor time ratio vs SPRAL improves from 4.19x to below 3.0x.
- **SC-003**: Assembly phase percentage of total factor time on c-71 decreases from 37.7% (measured via `profile_matrix` diagnostic).
- **SC-004**: Extraction phase percentage of total factor time on c-71 decreases from 34.8% (measured via `profile_matrix` diagnostic).
- **SC-005**: All 65 SuiteSparse matrices produce backward errors within 5e-11 with no regressions.
- **SC-006**: No matrix in the 65-matrix benchmark suite shows a factor time regression greater than 5% vs the post-9.1b baseline.
- **SC-007**: Median factor time ratio across the 65-matrix suite remains at or below 1.01x vs SPRAL.

## Assumptions

- The symbolic analysis phase can accommodate additional precomputation (scatter maps) without significantly increasing analysis time. Analysis is typically less than 10% of total solve time and runs once per sparsity pattern.
- Scatter map memory is unbounded and proportional to nnz (assembly maps) plus the sum of child contribution sizes (extend-add maps). This is accepted as part of symbolic analysis storage, matching SPRAL's approach. No memory budget or fallback path is needed.
- Column-major dense matrix storage provides contiguous column slices, enabling bulk copy operations on individual columns.
- The extend-add scatter map for contributions needs a hybrid approach: precompute the static portion (parent row indices) during symbolic analysis, then resolve the dynamic portion (delayed columns) at factorization time using the child's actual factorization result.
- Performance improvements primarily benefit matrices with high front-to-elimination ratios (c-71, c-big). Matrices that are already kernel-dominant (ncvxqp7, nd12k) may see negligible improvement.

## Out of Scope

- **Small leaf subtree fast path**: Simplicial matrix optimization via SmallLeafNumericSubtree is Phase 9.1d.
- **Pool allocator for contribution blocks**: SPRAL uses a pool allocator for contribution memory — a potential future optimization.
- **Amalgamation tuning** (nemin threshold changes): A separate investigation from assembly/extraction optimization.
- **Parallel assembly optimization**: Intra-node parallelism for assembly is out of scope — this feature targets sequential per-supernode overhead.
- **Two-phase assembly** (SPRAL's pre/post-factor assembly split): SPRAL assembles child contributions into the parent's contribution block in a post-factor step, bypassing the frontal matrix for the Schur complement portion. This architectural change would eliminate the contribution block copy from the parent's perspective but requires significant pipeline restructuring. Noted as a potential follow-up but out of scope for this feature.

## Dependencies

- Phase 9.1b (workspace reuse) — COMPLETE. This feature builds on the existing workspace infrastructure.
- Phase 9.1a (supernode amalgamation) — COMPLETE. Scatter maps must handle non-contiguous owned column ranges from amalgamation.
- Symbolic analysis result structure must be extended to store per-supernode scatter maps.
