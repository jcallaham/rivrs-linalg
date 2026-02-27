# Feature Specification: Direct GEMM into Contribution Buffer

**Feature Branch**: `024-direct-gemm-contrib`
**Created**: 2026-02-24 (revised)
**Status**: Draft
**Input**: Phase 9.1e from ssids-plan.md — eliminate the O(n²) contribution extraction copy by restructuring the BLAS-3 blocking loop to defer the Schur complement computation and write it directly into a pre-allocated contribution buffer, following SPRAL's zero-copy architecture.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Deferred Contribution GEMM (Priority: P1)

The multifrontal factorization's BLAS-3 blocking loop currently applies Schur complement updates to the **entire** trailing submatrix at every block iteration, including the non-fully-summed (NFS) region that will become the contribution block. Since no subsequent operation within the blocking loop reads the NFS×NFS region, these per-block updates to it are premature — the NFS×NFS Schur complement can be computed once, after the loop, in a single GEMM using all eliminated columns. This final GEMM should write its result directly into a pre-allocated contribution buffer, eliminating the post-factorization extraction copy entirely.

**Why this priority**: This addresses the architectural root cause identified in Phase 9.1c/9.1d. The profiling data showed that eliminating allocation overhead alone (9.1d pool experiment) netted zero improvement — the O(n²) copy itself is the bottleneck. Eliminating the copy requires changing where the GEMM writes, not just where the copy targets.

**Independent Test**: Can be verified by running the full SuiteSparse test suite (65 matrices) and confirming all backward errors remain below 5e-11. Diagnostic profiling should show ExtractContr near zero and the new contribution GEMM phase visible as a separate sub-phase.

**Acceptance Scenarios**:

1. **Given** a matrix with large fronts (e.g., c-71, max front ~2475, num_fully_summed ~26), **When** the solver factorizes it, **Then** per-block trailing updates within the blocking loop do not modify the NFS×NFS region (`A[p..m, p..m]` where `p = num_fully_summed`), and a single GEMM after the loop computes the NFS×NFS Schur complement directly into the contribution buffer.
2. **Given** a supernode with zero delayed columns (the common case), **When** the contribution is "extracted", **Then** no element-by-element or column-by-column data copy occurs — the contribution data is already in the buffer from the final GEMM, and only the row index array is built.
3. **Given** the full 65-matrix SuiteSparse test suite, **When** factorized with the restructured blocking loop, **Then** all backward errors remain within the same tolerance as Phase 9.1c (no regression).

---

### User Story 2 - Pre-Allocated Contribution Buffer (Priority: P1)

The solver should pre-allocate a contribution buffer as part of the factorization workspace, sized from symbolic analysis, and reuse it across supernodes. The final contribution GEMM writes directly into this buffer. This eliminates the thousands of per-supernode `Mat::zeros` allocations that cause 934K page faults and 3.1s of kernel (mmap/munmap) overhead.

**Why this priority**: The final GEMM in US1 needs a destination buffer. Pre-allocation ensures the buffer is warm in the TLB and avoids allocation overhead.

**Independent Test**: Can be verified by confirming zero per-supernode contribution allocations in the factorization loop (no `Mat::zeros` in `extract_contribution`), and by running `perf stat` to verify reduced page faults.

**Acceptance Scenarios**:

1. **Given** symbolic analysis results for any matrix, **When** the factorization workspace is constructed, **Then** a contribution buffer is pre-allocated with dimensions sufficient for the largest possible contribution block.
2. **Given** a factorization run, **When** monitored for allocations, **Then** no per-supernode contribution block allocations occur during the factorization loop (steady-state zero allocation via buffer swap between workspace and ContributionBlock).

---

### User Story 3 - Correct Handling of Delayed Columns (Priority: P1)

When the APTP kernel delays columns (pivots that fail the threshold check), the contribution block includes both delayed-column rows and NFS rows. The delayed-column portion and its cross-terms with NFS rows are computed during the blocking loop (they are in the fully-summed trailing region). Only the NFS×NFS portion is deferred. The contribution extraction must correctly assemble both portions into the contribution buffer: the small delayed portion copied from the workspace, and the NFS×NFS portion written by the final GEMM.

**Why this priority**: Delayed columns are the source of correctness complexity in this feature. The common case (zero delays) is straightforward, but the solver must handle arbitrary delay counts correctly.

**Independent Test**: Can be verified by factorizing matrices known to produce delayed columns (e.g., stokes128, bratu3d) and confirming backward error is unchanged.

**Acceptance Scenarios**:

1. **Given** a supernode where the APTP kernel delays `d` columns (ne < k), **When** the contribution is assembled, **Then** the contribution buffer contains: delayed×delayed interactions from the workspace, NFS×delayed cross-terms from the workspace, and NFS×NFS Schur complement from the final GEMM — all in the correct positions.
2. **Given** a supernode where all columns are delayed (ne = 0), **When** the contribution is assembled, **Then** the final GEMM is a no-op (rank-0 update), and the entire contribution is copied from the workspace (the assembled frontal matrix with no eliminations applied).

---

### User Story 4 - Parallel Path Compatibility (Priority: P2)

The parallel factorization path (rayon + thread-local workspaces) must continue to work correctly. Each thread's workspace should have its own pre-allocated contribution buffer. Cross-thread contribution transfer must continue to use owned data.

**Why this priority**: The parallel path is production-critical and must not regress. The core optimization (deferred GEMM) applies identically to both sequential and parallel paths since each thread processes its assigned supernodes independently.

**Independent Test**: Can be verified by running the full SuiteSparse suite with parallel factorization enabled.

**Acceptance Scenarios**:

1. **Given** the parallel factorization path, **When** a child supernode is factorized on thread A and the parent on thread B, **Then** the child's contribution (written by the final GEMM into thread A's workspace buffer, then moved into a ContributionBlock) is correctly transferred and assembled into the parent's frontal matrix on thread B.
2. **Given** the parallel path with thread-local workspaces, **When** multiple supernodes are factorized concurrently, **Then** each thread's contribution buffer is independent and no data races occur.

---

### Edge Cases

- What happens when a supernode has zero non-fully-summed rows (all rows are eliminated)? No contribution is produced; the final GEMM and contribution buffer are not used for that supernode.
- What happens when all columns are delayed (ne = 0)? The final GEMM has rank 0 (no eliminated columns), producing no output. The contribution is the assembled frontal matrix itself, copied from the workspace.
- What happens when the fully-summed count is very small relative to the front size (e.g., p=26, m=2475)? The per-block trailing updates become tiny (FS×FS is at most 26×26) while the deferred contribution GEMM is large (2449×2449). This is the primary performance win — the blocking loop becomes much cheaper.
- What happens when the fully-summed count equals the front size (leaf supernode with no pattern rows)? There is no NFS region, no contribution, and no deferred GEMM. The blocking loop behaves as before.
- What happens during the BLAS-3 blocking loop when columns are delayed (end_pos shrinks)? The `num_fully_summed` boundary `p` is fixed throughout the loop — it does not change when columns are delayed. `end_pos` shrinks but `p` stays the same. The NFS region `[p..m]` is stable.
- How does the deferred GEMM interact with the two-level factor path (`two_level_factor`) and the TPP small-front path (`tpp_factor_as_primary`)? Both paths produce the same workspace layout after completion (L21 in columns `[0..ne]`, D in MixedDiagonal). The deferred contribution GEMM is applied after whichever inner factorization path runs — it is orthogonal to the choice of pivoting strategy.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: Per-block trailing updates in the BLAS-3 blocking loop MUST skip the NFS×NFS region (`A[p..m, p..m]` where `p = num_fully_summed`), updating only the FS×FS region and cross-terms.
- **FR-002**: After the blocking loop completes, a single GEMM MUST compute the NFS×NFS Schur complement using the full L21 panel (rows `[p..m]`, all `ne` eliminated columns) and diagonal D, writing directly into the pre-allocated contribution buffer.
- **FR-003**: The solver MUST pre-allocate a contribution buffer on `FactorizationWorkspace`, sized from symbolic analysis to accommodate the largest contribution block across all supernodes.
- **FR-004**: The contribution extraction step MUST become index-only (constructing `row_indices` and `num_delayed`) for the NFS×NFS portion. When delayed columns exist, the small delayed portion and cross-terms MUST be copied from the workspace into the contribution buffer.
- **FR-005**: Delayed columns MUST be correctly positioned in the contribution buffer layout: rows/columns `0..num_delayed` for delayed pivots, `num_delayed..size` for the NFS Schur complement.
- **FR-006**: The parallel factorization path MUST use per-thread pre-allocated contribution buffers via the existing thread-local workspace mechanism.
- **FR-007**: All 65 SuiteSparse test matrices MUST continue to pass with backward error below 5e-11.
- **FR-008**: The solver MUST produce results within the same backward error tolerance (5e-11) as Phase 9.1c. Bit-identical results are not required — the GEMM decomposition change (many rank-`block_size` updates → one rank-`ne` update) changes floating-point operation ordering, which may produce last-ULP differences. Backward error is the meaningful correctness metric for a numerical solver.
- **FR-009**: The contribution buffer size MUST be predictable from symbolic analysis without requiring runtime reallocation during factorization (with dynamic resize as a fallback for delayed-column cascades, matching the frontal workspace pattern).
- **FR-010**: When a supernode produces no contribution (all rows eliminated), the solver MUST skip both the deferred GEMM and contribution extraction entirely.
- **FR-011**: The deferred contribution GEMM MUST work with all three factorization paths: `factor_inner` (standard blocking), `two_level_factor` (outer+inner blocking), and `tpp_factor_as_primary` (small fronts). In all cases the GEMM uses the workspace state after the path completes.
- **FR-012**: The original assembled NFS×NFS values (from original matrix scatter and child extend-add) MUST be correctly included in the contribution. The per-block trailing updates no longer modify the NFS×NFS region, so these values are preserved in the workspace. The deferred GEMM function copies them into the contribution buffer and applies the rank-`ne` update in-place, computing `contrib = assembled_NFS_NFS - L21_NFS * D * L21_NFS^T`.

### Key Entities

- **Contribution Buffer**: A pre-allocated dense matrix on `FactorizationWorkspace`, reused across supernodes via swap protocol. Sized to `max_front × max_front`. Receives the final Schur complement GEMM output for the NFS×NFS region.
- **Deferred Contribution GEMM**: A single symmetric rank-`ne` update computed after the blocking loop. Reads L21 and D from the workspace, writes to the contribution buffer. Replaces the many per-block rank-`block_size` updates to the NFS×NFS region.
- **NFS Boundary**: The fixed `num_fully_summed` value (`p`) that divides the frontal matrix into fully-summed columns (eligible for elimination) and non-fully-summed rows (contribution). This boundary is stable throughout the blocking loop.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The `ExtractContr` sub-phase (measured with `diagnostic` instrumentation) drops to near-zero as a fraction of factor time on c-71, confirming the copy is eliminated.
- **SC-002**: A new `ContribGEMM` sub-phase appears in diagnostic output, showing the direct contribution computation cost separately from the kernel.
- **SC-003**: System-level allocation overhead (page faults measured via `perf stat`) decreases significantly relative to Phase 9.1c baseline on c-71.
- **SC-004**: All 65 SuiteSparse test matrices pass with backward error below 5e-11, with no regressions from Phase 9.1c.
- **SC-005**: No new `unsafe` code is introduced (maintaining the all-safe-Rust property from Phase 8.2).

## Assumptions

- The NFS×NFS region of the trailing submatrix is not read by any operation within the BLAS-3 blocking loop (factor, TRSM, threshold check, cross-term update). This has been verified by analyzing the data dependencies in `factor_inner`.
- The final contribution GEMM produces the same Schur complement as the accumulated per-block updates (both compute `assembled_NFS_NFS - L21_NFS * D * L21_NFS^T`, just in different orderings).
- The contribution buffer sizing uses `max_front_size` as a conservative upper bound, accepting ~49 MB of additional memory for c-71-class problems.
- The two-level factor and TPP paths leave the workspace in the same state as `factor_inner` (L21 in columns `[0..ne]`, D in MixedDiagonal), so the deferred GEMM applies uniformly.

## Scope Boundary

**In scope**:
- Restructured per-block trailing update (skip NFS×NFS region)
- Final deferred contribution GEMM after blocking loop
- Pre-allocated contribution buffer on `FactorizationWorkspace`
- Index-only contribution extraction (plus small delayed-column copy)
- Sequential and parallel path compatibility
- Diagnostic instrumentation updates (new ContribGEMM sub-phase)

**Out of scope** (deferred to future phases):
- AppendAlloc / bump allocator for factor storage (saves 0.3% — not worth the lifetime propagation complexity now)
- BuddyAllocator for contributions (single swap buffer suffices)
- Split assembly into pre/post phases (SPRAL's `assemble_pre`/`assemble_post`)
- Task DAG scheduling (independent of memory architecture)
- Huge pages / `madvise(MADV_HUGEPAGE)` (orthogonal optimization)

## Dependencies

- Phase 9.1c (assembly/extraction optimization) — COMPLETE: provides `AssemblyMaps`, sub-phase timing, baseline measurements
- Phase 9.1d (architecture review) — COMPLETE: provides analysis showing pool-based reuse is insufficient; recommends direct GEMM; documents SPRAL's architecture
- Phase 8.2 (parallel factorization) — COMPLETE: provides thread-local workspace and parallel contribution transfer

## Algorithm References

The multifrontal method and contribution block management are described in:

- **Liu (1992)** — "The Multifrontal Method for Sparse Matrix Solution: Theory and Practice" (`/workspace/rivrs-linalg/references/ssids/liu1992.md`): Foundational theory of update matrices, frontal assembly via extend-add (Theorem 4.1), stack-based storage management (Section 5).
- **Hogg et al. (2016)** — "A Sparse Symmetric Indefinite Direct Solver for GPU Architectures" (`/workspace/rivrs-linalg/references/ssids/hogg2016.md`): SSIDS implementation. Generated element computation (lines 289-297): "the majority of floating-point operations typically derive from the formation of the generated element `C ← C - L2 D1 L2^T`". Sparse assembly of child contributions (lines 208-242).
- **Duff et al. (2005)** — "Scaling and Pivoting in the Multifrontal LDL^T Factorization (MA57)" (`/workspace/rivrs-linalg/references/ssids/duff2005.md`): Multifrontal framework (lines 96-108), delayed pivot handling (line 132), assembly as performance-dominant phase (line 604).
- **Ng (1993)** — "Supernodal Multifrontal Cholesky Factorization" (`/workspace/rivrs-linalg/references/ssids/ng1993.md`): Supernodal assembly loop (Figure 7, lines 148-160), dense relative concept for direct updates without indirect indexing (lines 203-218).
- **Duff and Reid (1984)** — "The Multifrontal Solution of Indefinite Sparse Symmetric Linear Equations" (`/workspace/rivrs-linalg/references/ssids/duff1984.md`): Original multifrontal indefinite method with stack-based contribution management.

SPRAL source code (BSD-3, freely consultable):
- `factor.hxx:92-103` — Direct GEMM into `node.contrib` (the architecture this feature adopts)
- `NumericNode.hxx` — Per-node memory layout (`lcol`, `contrib`, `perm`)
- `assemble.hxx:27-38` — Column-oriented scatter (`asm_col`, 4x unrolled)
- `AppendAlloc.hxx` — Bump allocator for factor storage (out of scope)
- `BuddyAllocator.hxx` — Buddy system for contributions (out of scope)

Internal references:
- `dev/phase9/phase-9.1c-profiling-report.md` — Sub-phase profiling identifying contribution extraction as 40.1% of factor time
- `dev/phase9/phase-9.1d-contrib-architecture-review.md` — Architecture comparison; pool-based reuse insufficient; SPRAL zero-copy analysis
