# Research: Two-Level APTP Factorization

**Feature**: 017-two-level-aptp
**Date**: 2026-02-16

## Research Topics & Decisions

### R1: Innermost Pivoting Strategy

**Decision**: Complete pivoting (Algorithm 4.1, Duff et al. 2020)

**Rationale**: For ib×ib blocks (ib=32), searching the entire remaining submatrix for the largest entry has negligible overhead since the block fits in L1 cache. Complete pivoting is provably at least as stable as threshold partial pivoting with u=0.25 (Duff et al. 2020, Section 4.1). It guarantees all L entries are bounded by 4 across three execution paths (1×1, successful 2×2, and failed-2×2 fallback to 1×1). Critically, complete pivoting never delays columns within the ib×ib block — it always finds a valid pivot unless the block is numerically singular. This reduces delay propagation to the outer level.

**Alternatives considered**:
- Current 1×1 threshold + 2×2 BK fallback: Simpler (reuses existing code) but can delay columns within inner blocks, increasing outer-level complexity. The threshold u=0.01 bounds L entries by 100, less tight than complete pivoting's bound of 4.
- Rook pivoting: Intermediate between partial and complete. More search than threshold but less than complete. Not used by SPRAL or the paper. No clear advantage for 32×32 blocks.

**References**:
- Duff, Hogg & Lopez (2020), Section 4, Algorithm 4.1 (`/workspace/rivrs-linalg/references/ssids/duff2020.md`)
- SPRAL: `block_ldlt()` in `ldlt_app.cxx` uses complete pivoting at innermost level

### R2: Backup Strategy

**Decision**: Per-block CopyBackup (save column of blocks before each outer block's Factor)

**Rationale**: Before factoring outer block j, save the diagonal block and all sub-diagonal blocks in that block column. If the a posteriori check reduces nelim, restore the failed columns from backup. This is SPRAL's `CopyBackup<T>` approach. Memory cost is O(nb × m) per outer block, where m is the frontal matrix dimension — small compared to the O(m²) frontal matrix.

**Alternatives considered**:
- Per-column backup (current): Backs up individual columns inside try_1x1_pivot/try_2x2_pivot. Incompatible with two-level APTP where the Factor phase processes the entire diagonal block speculatively before checking stability.
- Full-matrix backup (APP_AGGRESSIVE): Backs up the entire frontal matrix once, restores on any failure. Higher memory cost (2× frontal matrix) but simpler failure recovery. SPRAL supports this but CopyBackup is the default. Deferred unless benchmarking shows excessive per-block backup overhead.
- Pool backup (PoolBackup): Dynamic pool of block-sized allocations for multi-threaded use. Unnecessary for Phase 8.1 (sequential). Relevant for Phase 8.2.

**References**:
- SPRAL: `CopyBackup<T>` class (lines 423-597 in `ldlt_app.cxx`)
- SPRAL: `PoolBackup<T>` class (lines 606-760, for multi-threaded — Phase 8.2)

### R3: Block Size Defaults

**Decision**: Outer nb=256, inner ib=32 (matching SPRAL and the paper)

**Rationale**: These values were empirically chosen by the SPRAL authors for Intel Broadwell-era CPUs. Phase 7 profiling confirms 91% of supernodal matrices benefit from nb=256 blocking (max_front > 256). ib=32 fits in L1 cache (32×32 doubles = 8KB, typical L1 = 32-64KB). Block sizes are configurable via AptpOptions for future tuning.

**Alternatives considered**:
- Larger nb (512): Would reduce parallelism opportunities for Phase 8.2 but might improve GEMM efficiency. Can be explored via benchmarking.
- Smaller ib (16): Reduces inner-loop overhead but increases outer-block iterations. 32 is the sweet spot per SPRAL experiments.
- Dynamic nb based on front size: Adds complexity with unclear benefit. Fixed nb with automatic dispatch (two-level iff front > nb) is simpler.

**References**:
- Duff, Hogg & Lopez (2020), line 130: "We chose the values of 256 for nb and 32 for ib"
- SPRAL: `INNER_BLOCK_SIZE = 32` (line 41 in `ldlt_app.cxx`)

### R4: BLAS-3 Operation Mapping

**Decision**: Use faer's `matmul` and `solve_unit_lower_triangular_in_place` for Update and Apply phases

**Rationale**: faer provides high-performance dense operations that dispatch to optimized kernels. Using faer's APIs maintains transparent composition and avoids external BLAS dependencies.

**Specific API mapping**:
- **Apply (TRSM)**: `faer::linalg::triangular_solve::solve_unit_lower_triangular_in_place` for L21 = A21 × L11^{-T}. Then scale by D11^{-1} using MixedDiagonal's existing solve machinery.
- **Update (GEMM)**: `faer::linalg::matmul::matmul` for A22 -= L21 × D11 × L21^T. Compute W = L21 × D11 in workspace (O(nb × m) temporary), then GEMM.

**Alternatives considered**:
- Raw BLAS via FFI: Higher performance ceiling but adds external dependency (BLAS library). Violates constitution preference for pure Rust. faer's matmul is already well-optimized.
- Keep rank-1/rank-2 updates: No benefit from blocking. Defeats the purpose of Phase 8.1.

**References**:
- faer API: `faer::linalg::matmul::matmul` (already used in solve.rs)
- faer API: `faer::linalg::triangular_solve` module (already used in solve.rs)

### R5: Handling of Previously-Delayed Columns (UpdateNT/UpdateTN)

**Decision**: After each outer block j is factored (with nelim_j successful pivots), apply the rank-nelim_j update to all previously-delayed columns from blocks 0..j-1. This corresponds to UpdateNT and UpdateTN from Algorithm 3.1.

**Rationale**: Delayed columns from earlier blocks have not received the Schur complement contributions from block j's pivots. Without these updates, the delayed columns would be stale when eventually passed to the parent node or re-factored. SPRAL applies these updates as part of the block loop — they are small-granularity tasks but have low impact due to few delays in practice.

**Implementation**: After the Update phase (GEMM on trailing submatrix), iterate over each previously-delayed column region and apply the rank-nelim_j update. This uses the same L21 data from the Apply phase.

**References**:
- Duff, Hogg & Lopez (2020), Algorithm 3.1: UpdateNT, UpdateTN kernel descriptions

### R6: Interaction with AptpOptions and Fallback Strategy

**Decision**: The `AptpFallback` enum (BunchKaufman vs Delay) applies to the MIDDLE level of the two-level hierarchy — the inner APTP loop processing ib-sized sub-blocks within an nb-sized outer block. At the innermost level (ib×ib diagonal), complete pivoting is always used regardless of fallback setting. At the outer level, failure is handled by backup/restore.

**Rationale**: The three levels have different failure semantics:
- **Innermost (ib×ib)**: Complete pivoting — no failures (always finds a pivot unless singular)
- **Middle (ib sub-blocks within nb)**: APTP with configurable fallback — may delay sub-columns
- **Outer (nb blocks)**: A posteriori check on Apply — may reduce nelim and restore from backup

The existing `AptpFallback` setting naturally maps to the middle level, which is the closest analog to the current single-level behavior.

### R7: Testing Against Single-Level Baseline

**Decision**: Two-level and single-level will NOT produce bitwise-identical results, even on the same input. Tests compare tolerance-based metrics (reconstruction error < 1e-12, backward error < 5e-11) rather than exact numerical equality.

**Rationale**: BLAS-3 operations (GEMM, TRSM) perform floating-point arithmetic in a different order than column-by-column rank-1 updates. The results are equivalent to within machine epsilon accumulation but not bitwise identical. Additionally, complete pivoting at the innermost level may choose different pivot orderings than the current 1x1+BK strategy, leading to numerically equivalent but permutation-different factorizations.

**Testing strategy**:
1. Complete pivoting unit tests: verify Algorithm 4.1 on known matrices
2. Two-level kernel tests: verify reconstruction < 1e-12 on random matrices of size 64, 128, 256, 512, 1024
3. Integration tests: verify backward error no regression on CI SuiteSparse matrices
4. Benchmark: compare single-level vs two-level wall-clock time on isolated kernel calls
