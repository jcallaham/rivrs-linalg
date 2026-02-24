# Research: Direct GEMM into Contribution Buffer

**Feature**: 024-direct-gemm-contrib
**Date**: 2026-02-24 (revised)

## R1: Deferred Contribution GEMM — Eliminating the Copy at Source

### Decision

Restructure the BLAS-3 blocking loop so that per-block trailing updates skip the NFS×NFS (non-fully-summed × non-fully-summed) region of the frontal matrix. After the loop, compute the entire NFS×NFS Schur complement in a single GEMM that writes directly to a pre-allocated contribution buffer. This eliminates both the O(n²) copy and the per-supernode allocation.

### Background: Why Pre-Allocation Alone Is Insufficient

Phase 9.1d measured what happens when allocation overhead is eliminated via a buffer pool:

| Metric | 9.1c (baseline) | 9.1d (pool) |
|--------|-----------------|-------------|
| sys time | 3.1s (32%) | 0.49s (7.3%) |
| page faults | 934K | 317K |
| **Factor time** | **5,920 ms** | **~6,050 ms** |

Allocation elimination was a wash. The pool's oversized buffers caused cache/TLB pollution that offset the allocation savings. But even with a properly-sized single buffer (avoiding the pool's pathology), the fundamental issue remains: copying 48 MB of data per large supernode (2449² elements for c-71's largest fronts) is a bandwidth-bound operation that cannot be eliminated by changing the destination.

The ssids-plan.md Section 9.1e explicitly recommends modifying `factor_inner` / `update_trailing` so the GEMM writes to the contribution buffer. The question is how.

### Architecture: Per-Block GEMM Decomposition

Each block's `update_trailing` currently computes a rank-`nelim` symmetric update on the entire trailing submatrix `A[block_end..m, block_end..m]`. This region decomposes into three sub-regions relative to the fixed `num_fully_summed` boundary (`p`):

```
                    block_end (ts)        p (num_fully_summed)        m
                        │                       │                     │
block_end (ts)     ─────┼───────────────────────┼─────────────────────┤
                        │                       │                     │
                        │     FS × FS           │     FS × NFS        │
                        │    (region 1)         │    (region 2)       │
                        │                       │                     │
p (num_fully_summed)────┼───────────────────────┼─────────────────────┤
                        │                       │                     │
                        │     NFS × FS          │     NFS × NFS       │
                        │    (region 2')        │    (region 3)       │
                        │                       │                     │
m                  ─────┼───────────────────────┼─────────────────────┘
```

- **Region 1** (FS×FS): Future fully-summed blocks. Subsequent iterations' factor and TRSM steps read these values. **Must be updated per-block.**
- **Region 2/2'** (cross-terms): NFS rows × FS columns and FS rows × NFS columns. Subsequent TRSM operations read the NFS rows of FS columns. **Must be updated per-block.**
- **Region 3** (NFS×NFS): The Schur complement proper. **No subsequent operation in the blocking loop reads from this region.** It is only consumed after factorization by `extract_contribution`.

SPRAL exploits this: the blocking loop updates regions 1+2 only, and a single GEMM after the loop computes region 3 directly into `node.contrib` (`factor.hxx:92-103`).

### Why Region 3 Can Be Deferred

The blocking loop in `factor_inner` iterates `while k < end_pos` where `end_pos` starts at `p` (num_fully_summed). Each iteration:

1. **Factor**: Operates on diagonal block `A[k..block_end, k..block_end]` — within FS region
2. **TRSM**: Computes L21 panel `A[block_end..m, k..block_end]` — reads FS and NFS rows of FS columns
3. **Trailing update**: `A[ts..m, ts..m] -= W * L21^T` — updates regions 1+2+3

The TRSM in step 2 reads from columns `[k..block_end]` (fully-summed), rows `[block_end..m]` (all trailing). These values include accumulated updates from previous iterations' trailing updates to regions 1+2. But the TRSM does **not** read from `A[p..m, p..m]` (region 3) — it reads columns in `[0..p]`, not `[p..m]`.

Therefore: skipping region 3's update during the loop does not affect any subsequent factor, TRSM, or threshold check within the loop. Region 3 can be computed once, after the loop, using the final accumulated L21 and D.

### Final Contribution GEMM

After the blocking loop, the factored workspace contains:
- `A[0..ne, 0..ne]`: L11, D11 (eliminated factors)
- `A[ne..m, 0..ne]`: L21 (from per-block TRSM, accumulated correctly because cross-terms were updated)
- `A[ne..p, ne..p]`: Delayed-column interactions (updated during the loop — region 1)
- `A[p..m, ne..p]`: NFS×delayed cross-terms (updated during the loop — region 2)
- `A[p..m, p..m]`: **Stale** — still contains original assembled values, no Schur complement updates applied

The final GEMM computes:
```
contrib[NFS×NFS] = assembled[NFS×NFS] - L21_NFS * D * L21_NFS^T
```
where `L21_NFS = A[p..m, 0..ne]` and D is the accumulated diagonal from the blocking loop. This is a single symmetric rank-`ne` update on the `(m-p) × (m-p)` contribution, writing directly to the pre-allocated buffer. For c-71's largest front (m=2475, p=26, ne=26): one GEMM of 2449×2449 with rank 26.

### Performance Characteristics

The single final GEMM replaces many per-block rank-`block_size` updates to the NFS×NFS region. Both approaches perform the same total FLOPs (the Schur complement has a fixed arithmetic cost regardless of blocking), but the single GEMM has:

- **Better BLAS-3 efficiency**: One large GEMM with rank `ne` (all eliminated columns at once) vs many small rank-`block_size` GEMMs. Larger rank means better cache reuse in the BLAS kernel.
- **Smaller per-block GEMMs**: Each block's trailing update shrinks from `(m-ts) × (m-ts)` to `(p-ts) × (p-ts)` plus a rectangular `(m-p) × (p-ts)` cross-term. For c-71 with p=26, the FS×FS region is at most 26×26 per block — negligible.
- **Zero copy**: The contribution data is written once, to its final location.

### Alternatives Considered

| Option | Approach | Outcome |
|--------|----------|---------|
| Pre-allocate + swap (original plan) | Eliminate allocation overhead, copy to warm memory | 9.1d data shows allocation elimination alone doesn't help. Retains O(n²) copy. |
| Per-block GEMM split to two destinations | Each block writes FS portion to workspace, NFS portion to contrib buffer | Same total FLOPs but splits each GEMM into two smaller ones. Less efficient than deferring. |
| **Deferred contribution GEMM (chosen)** | Skip NFS×NFS per-block, single final GEMM to contrib buffer | Zero copy, potentially faster GEMM (larger rank), cleaner separation of concerns |

## R2: Contribution Buffer Sizing from Symbolic Analysis

### Decision

Size the contribution buffer to `max_front_size × max_front_size` (same as the frontal workspace), reusing the existing `estimate_max_front_size()` computation.

### Rationale

The contribution dimension for a supernode is `m - ne`, which depends on runtime pivot outcomes. The absolute worst case is `max_front - 1`. Using `max_front_size` as the buffer dimension guarantees no runtime reallocation, matching the frontal workspace's sizing strategy. The workspace already has `ensure_capacity()` for the rare case where delayed column cascades exceed estimates.

Memory cost: one additional `max_front² × 8` bytes. For c-71 (max_front ~2475), this is ~49 MB — comparable to the existing frontal workspace allocation.

## R3: Delayed Column Handling

### Decision

The deferred GEMM handles the NFS×NFS region (`A[p..m, p..m]`). Delayed columns (`A[ne..p]`) and their cross-terms with NFS rows (`A[p..m, ne..p]`) are already correctly computed in the workspace during the blocking loop.

### Rationale

After the blocking loop:
- Delayed-column interactions `A[ne..p, ne..p]` received all rank-k updates during the loop (they were in the FS trailing region at each iteration)
- Cross-terms `A[p..m, ne..p]` were also updated per-block (region 2 in R1's decomposition)
- Only `A[p..m, p..m]` was skipped

For the contribution block, the data layout is:
- `[0..num_delayed, 0..num_delayed]`: delayed×delayed — **copy from workspace** (small: typically 0-10 rows)
- `[num_delayed..size, 0..num_delayed]`: NFS×delayed cross-term — **copy from workspace** (small)
- `[num_delayed..size, num_delayed..size]`: NFS×NFS Schur complement — **written by final GEMM** (dominant size)

In the common case (`ne = p`, no delayed columns), `num_delayed = 0` and the contribution is entirely NFS×NFS — true zero copy. When delayed columns exist, only the small delayed portion needs copying from the workspace. For c-71, the median number of delayed columns per supernode is ~0 with occasional outliers of a few columns.

## R4: ContributionBlock Ownership for Parallel Path

### Decision

Keep `ContributionBlock` as an owned type with `data: Mat<f64>`. The final GEMM writes into `workspace.contrib_buffer`, which is then moved into `ContributionBlock` (same swap protocol as the original plan). The parallel path uses per-thread workspace buffers via `Cell<FactorizationWorkspace>`.

### Rationale

The parallel path requires owned data for cross-thread transfer. The contribution buffer lives on `FactorizationWorkspace`, which is thread-local. After the final GEMM writes to the buffer, the buffer is moved into `ContributionBlock`. When the parent's extend-add consumes the contribution, the buffer is recycled back to the workspace.

For the sequential path, this creates a zero-allocation steady state: one buffer moves between workspace and ContributionBlock without ever being allocated or freed. For the parallel path, when contributions cross threads, the sending thread's buffer goes with the contribution and the thread lazily reallocates when needed.

## R5: Interaction with `extract_front_factors`

### Decision

The final contribution GEMM must execute while L21 and D are still in the frontal workspace — before the workspace is reused for the next supernode. The natural ordering is: blocking loop → final contribution GEMM → `extract_front_factors` → index-only `extract_contribution`. Since `extract_front_factors` reads from but does not modify the workspace, the contribution GEMM and factor extraction can happen in either order.

### Rationale

The final GEMM reads `L21_NFS = A[p..m, 0..ne]` and D from the blocking loop's `MixedDiagonal`. Both are in the workspace. `extract_front_factors` also reads from the workspace (copying L11, D, L21 into owned `FrontFactors`). Neither operation modifies the workspace, so they can be ordered flexibly. The only constraint is that both must complete before the workspace is zeroed for the next supernode.

## R6: Custom Allocators (AppendAlloc / BuddyAllocator)

### Decision

Neither custom allocator is needed for this feature.

### Rationale

- **AppendAlloc** (bump allocator for factor storage): Only saves 0.3% (`extract_front_factors`). Adopting it would require adding lifetime parameters to `FrontFactors`, `AptpNumeric`, and all downstream consumers — significant refactoring for minimal gain. Deferred.

- **BuddyAllocator** (for contribution blocks): SPRAL uses this to manage concurrent contribution buffers in the parallel case. Our simpler approach — one pre-allocated buffer per thread-local workspace, moved via swap — avoids the complexity of a buddy allocator while achieving the same goal (zero steady-state allocation). The buddy allocator's power-of-2 sizing and buddy merging are irrelevant when each thread uses a single reusable buffer.

Both allocators are independent optimizations that can be pursued later without affecting the contribution GEMM architecture.
