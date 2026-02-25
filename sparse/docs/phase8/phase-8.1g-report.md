# Phase 8.1g: Sequential Profiling & Optimization Report

**Branch**: `018-sequential-profiling-optimization`
**Date**: 2026-02-20

## 1. Summary

Phase 8.1g instrumented the multifrontal factorization loop with per-supernode
timing, optimized allocation hotspots in the dense APTP kernel, collected
structured performance baselines for all 65 SuiteSparse matrices, and produced
a workload distribution analysis to inform Phase 8.2 parallelism strategy.

**Key finding**: 41/65 matrices (63%) are classified as **IntraNode** — the top
10% of supernodes by front size account for >80% of total factorization time.
Intra-node BLAS-3 parallelism (parallel dense operations within large fronts)
should be the **primary** parallelism strategy for Phase 8.2.

## 2. Per-Supernode Timing Instrumentation

Added `#[cfg(feature = "diagnostic")]` instrumentation to
`AptpNumeric::factor()` in `src/aptp/numeric.rs`:

- Three `Duration` fields on `PerSupernodeStats`: `assembly_time`,
  `kernel_time`, `extraction_time`
- Three aggregate `Duration` fields on `FactorizationStats`:
  `total_assembly_time`, `total_kernel_time`, `total_extraction_time`
- `ProfileSession` with top-level `factor_loop` section guard
- Per-supernode timing via `Instant::now()` (no per-supernode section guards —
  see Section 6 for rationale)
- `FinishedSession` stored on `AptpNumeric` for Chrome Trace export

**Zero overhead** when `diagnostic` feature is not enabled — all timing code
is behind `#[cfg(feature = "diagnostic")]` conditional compilation.

## 3. Allocation Optimization

Identified and fixed three allocation hotspots in `factor_inner`
(`src/aptp/factor.rs`):

| Hotspot | Before | After |
|---------|--------|-------|
| Panel row permutation buffer | `Vec<f64>` per row | Pre-allocated before loop, reused |
| Row permutation temp buffer | `vec![0.0; block_size]` per block | Pre-allocated before loop, reused |
| Column order slice copy | `.to_vec()` per block | Pre-allocated `Vec<usize>`, `copy_from_slice` |

**BlockBackup reuse** was evaluated and deferred to Phase 9.1 (arena memory) —
faer's `Mat` doesn't support in-place resize, and dimensions vary per block
iteration.

All 65/65 SuiteSparse matrices pass correctness checks after optimization.

## 4. Sequential Performance Baselines (CI Subset)

Collected via `examples/baseline_collection.rs` with `--ci-only` flag.
Ordering: MatchOrderMetis (MC64+METIS). Platform: linux aarch64.

| Matrix | n | nnz | ord_ms | fac_ms | slv_ms | backward_err |
|--------|---|-----|--------|--------|--------|-------------|
| bloweybq | 10,001 | 39,996 | 18.6 | 17.5 | 1.6 | 2.67e-18 |
| sparsine | 50,000 | 799,494 | 1,886 | 112,002 | 2,680 | 3.70e-15 |
| t2dal | 4,257 | 20,861 | 18.1 | 23.4 | 0.7 | 2.71e-18 |
| bratu3d | 27,792 | 88,627 | 121 | 530 | 10.7 | 2.47e-18 |
| cvxqp3 | 17,500 | 69,981 | 279 | 454 | 9.2 | 4.78e-14 |
| ncvxqp1 | 12,111 | 40,537 | 160 | 171 | 3.5 | 1.29e-16 |
| ncvxqp3 | 75,000 | 274,982 | 3,726 | 4,398 | 139 | 4.24e-14 |
| stokes128 | 49,666 | 295,938 | 311 | 203 | 10.1 | 3.74e-18 |
| cfd2 | 123,440 | 1,605,669 | 1,769 | 3,346 | 260 | 8.37e-19 |

**Observations**:
- sparsine dominates (112s factor) due to max_front=11,125 — cubic dense kernel
- Ordering (MC64+METIS) is a significant fraction for large matrices
  (ncvxqp3: 3.7s ordering vs 4.4s factor)
- All backward errors well below SPRAL's 5e-11 threshold
- Solve time is 1-5% of factor time across all matrices

## 5. Workload Distribution Analysis

Ran `examples/workload_analysis.rs` on all 65 SuiteSparse matrices.
Classification based on top-10% front size time fraction:
- **IntraNode** (>80%): 41 matrices (63%)
- **Mixed** (30-80%): 19 matrices (29%)
- **TreeLevel** (<30%): 5 matrices (8%)

### Classification Summary

**TreeLevel matrices** (5): bloweybq, dixmaanl, linverse, spmsrtls, blockqp1
- Characterized by many small supernodes (max_front 10-12), very even workload
- These are structurally trivial (diagonal or near-diagonal after ordering)
- Tree-level parallelism would help but these are already fast

**IntraNode matrices** (41): Most FEM/CFD/optimization problems
- Top 10% of fronts consume 80-99.8% of factorization time
- max_front ranges from 577 (stokes128) to 11,125 (sparsine)
- The single largest front accounts for 2-26% of total time
- These benefit most from parallel dense BLAS within large fronts

**Mixed matrices** (19): Moderate concentration (30-80% in top 10%)
- Includes large-n matrices with moderate front sizes (G3_circuit, ldoor, msdoor)
- Both tree-level and intra-node parallelism contribute
- Typically have many medium-size fronts (100-500) rather than few huge ones

### Front Size Distribution

Across all 65 matrices:
- **1-10**: Majority of supernodes by count (typically 40-80% of supernodes)
  but contribute <10% of factorization time
- **11-50**: Second largest count bucket, typically <10% of time
- **51-100**: Small count, moderate time contribution
- **101-500**: Moderate count, significant time (10-30% typical)
- **501-1k**: Small count, large time contribution
- **1k+**: Very few supernodes but dominant time (50-97% for IntraNode matrices)

**Key insight**: The distribution is extremely skewed. A handful of large fronts
dominate total factorization time. This makes intra-node BLAS-3 parallelism the
highest-impact optimization.

## 6. Profiling Infrastructure Findings

### O(N^2) Section Tree Bug (Fixed)

The `ProfileSession::finish()` method's `build_section_tree` had O(N^2)
complexity in `attach_children` — for each parent section, it scanned ALL events
to find children. With ~20K events (4 per supernode × 5K supernodes), this
caused the factorization to take minutes instead of milliseconds.

**Fix**: Removed per-supernode section guards (`supernode_N`, `assembly`,
`dense_kernel`, `extraction`) from the hot path. Per-supernode timing is
captured via `Instant::now()` and stored directly in `PerSupernodeStats`.
The `ProfileSession` retains only the top-level `factor_loop` guard, making
`finish()` O(1).

**Future improvement**: If Chrome Trace per-supernode output is needed, the
`build_section_tree` algorithm should be refactored to O(N) using a stack-based
approach rather than nested scanning. This is a Phase 9.1 item.

## 7. Phase 8.2 Parallelism Recommendation

### Primary Strategy: Intra-Node BLAS-3 Parallelism

**Evidence**: 41/65 matrices have >80% of factorization time concentrated in
the top 10% of supernodes by front size. For these matrices, parallelizing the
dense TRSM and GEMM operations within large fronts would directly accelerate
the bottleneck.

**Implementation**: Pass `faer::Par::rayon(nthreads)` to dense matmul and
triangular solve calls in `apply_and_check` and `update_trailing` within
`factor_inner`. faer's BLAS-3 kernels already support parallel execution via
the `Par` parameter — this is pure plumbing.

### Secondary Strategy: Tree-Level Parallelism

**Evidence**: 19 Mixed matrices have moderate workload distribution across
fronts. Tree-level parallelism (processing independent subtrees in parallel)
would help for these. The 5 TreeLevel matrices would also benefit but are
already fast (trivial structure).

**Implementation**: Level-set scheduling via Rayon for independent supernodes
at the same tree level. More complex than intra-node (requires dependency
tracking) but necessary for Mixed matrices.

### Recommended Phase 8.2 Order

1. **Intra-node first**: Add `Par` parameter to factor/solve, pass to faer
   dense operations. Low complexity, high impact on 63% of matrices.
2. **Tree-level second**: Add level-set scheduling for independent subtrees.
   Higher complexity, benefits remaining 29% of matrices.
3. **Benchmark both**: Compare speedup on representative matrices from each
   class (IntraNode: sparsine/nd12k, Mixed: G3_circuit/ldoor, TreeLevel:
   bloweybq/dixmaanl).

## 8. Tools Produced

- `examples/baseline_collection.rs`: Structured JSON baseline collection
  (per-phase timing, per-supernode stats, backward error, peak RSS).
  Supports `--ci-only` and `--compare <prev.json>` flags.
- `examples/workload_analysis.rs`: Workload distribution analysis and
  parallelism classification (TreeLevel/IntraNode/Mixed).
- JSON baselines: `target/benchmarks/baselines/baseline-<timestamp>.json`
