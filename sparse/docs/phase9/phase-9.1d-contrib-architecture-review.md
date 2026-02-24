# Phase 9.1d: Contribution Block Architecture Review

**Branch**: `023-contrib-workspace-reuse`
**Date**: 2026-02-23
**Status**: Feature abandoned after investigation. This document records the
findings and recommends an alternative approach based on SPRAL's architecture.

## 1. Problem Statement

Phase 9.1c profiling identified contribution block allocation as the dominant
bottleneck in multifrontal factorization:

| Sub-Phase       | Time (ms) | % of Factor | Notes |
|-----------------|-----------|-------------|-------|
| ExtractContr    | 2,374     | 40.1%       | Copy trailing submatrix to new `Mat<f64>` |
| Extend-add      | 1,974     | 33.3%       | Scatter child contribution into parent |
| Kernel          | 1,366     | 23.1%       | Dense APTP factorization |

Hardware counters showed:
- 934K page faults (first-touch on fresh `mmap` allocations)
- 3.1s sys time (32% of wall time — kernel mmap/munmap handling)
- 644B dTLB-load-misses (virtual address thrashing)

## 2. What We Tried: Pool-Based Contribution Reuse

Phase 9.1d implemented a free-list buffer pool on `FactorizationWorkspace`:

- **Pool mechanism**: `Vec<Mat<f64>>` free list with linear scan for first-fit.
  `take_from_pool` finds a buffer with `nrows >= size`, zeros the logical region,
  and returns it. `return_to_pool` pushes the buffer back after the parent's
  extend-add consumes it.

- **DFS traversal**: New iterative DFS postorder path (`factor_tree_dfs`) for
  the sequential case, with the pool enabling buffer recycling across siblings
  and cousins in the tree.

- **Level-set path**: Existing level-set traversal retained for parallel mode,
  with per-thread pools via `Cell<FactorizationWorkspace>`.

### Results

The pool eliminated allocation syscall overhead but did not improve factor time:

| Metric                | Phase 9.1c (before) | Phase 9.1d (pool) |
|-----------------------|---------------------|--------------------|
| sys time              | 3.1s (32%)          | 0.49s (7.3%)       |
| page faults           | 934K                | 317K                |
| Factor time (c-71)    | 5,920 ms            | ~6,050 ms           |
| ExtractContr          | 40.1%               | 19.3%               |
| Extend-add            | 33.3%               | 40.1%               |
| Kernel                | 23.1%               | 33.1%               |

ExtractContr dropped in absolute time (~2,374ms to ~1,172ms), but extend-add and
kernel times increased by comparable amounts, netting to zero improvement.

### Pool diagnostic data (c-71, DFS path)

| Metric | DFS | Level-Set |
|--------|-----|-----------|
| Hit rate | 77.3% | 18.4% |
| Waste (oversized buffers) | 11,902 MB | 3,943 MB |
| Max oversize ratio | **2374x** | 21x |

The DFS traversal causes catastrophic pool buffer oversizing: a buffer allocated
for a 2374-row contribution near the root is reused for a 1-row leaf contribution.
The physical pages remain resident, polluting caches and TLB. The level-set path
has poor hit rate (18.4%) because entire waves of contributions are consumed before
any buffers are returned to the pool.

### Conclusion

The pool addresses allocation overhead (syscall cost) but not the fundamental
problem: we are performing an O(n^2) copy of the Schur complement out of the
frontal workspace for every supernode. This copy exists because of an
architectural choice — using a single reusable frontal workspace — not because
of inefficient allocation.

## 3. How SPRAL Handles This

SPRAL avoids the contribution copy entirely through a different memory
architecture. The key insight is that **SPRAL does not have a reusable frontal
workspace**. Each node owns its own persistent memory.

### 3.1 Two-allocator design

SPRAL uses two specialized allocators:

1. **AppendAlloc** (factor storage): A bump/append-only allocator. Each node's
   `lcol` (factor matrix L and diagonal D) is allocated from a linked list of
   8MB+ pages via sequential pointer bump. Memory is calloc'd (zero-initialized)
   and never freed during factorization. Thread-safe via `#pragma omp critical`.

2. **BuddyAllocator** (contribution pool): A 16-level buddy-system allocator
   for transient contribution blocks. Allocations are power-of-2 sized with
   buddy merging on deallocation. Pre-sized to `max_front^2` to avoid growth.
   Thread-safe via OpenMP lock.

### 3.2 Per-node memory layout

Each `NumericNode` has three separately allocated regions:

```
node.lcol   → [L factor (m × n, leading dim = align_lda(m)) | D factor (2 × n)]
node.perm   → [permutation array (n integers)]
node.contrib → [Schur complement ((m-n) × (m-n), column-major, no padding)]
```

Where `m = nrow + ndelay_in`, `n = ncol + ndelay_in`.

`lcol` and `perm` are allocated from the AppendAlloc (permanent — never freed).
`contrib` is allocated from the BuddyAllocator (transient — freed after parent
consumes it).

### 3.3 Zero-copy factor storage

In SPRAL, the assembled frontal matrix IS `node.lcol`. Original matrix entries
and child contributions are scattered directly into `lcol` during `assemble_pre`.
After factorization, L11/L21/D remain in-place — no extraction needed. This
eliminates our `extract_front_factors` step (0.3% of factor time — small but
non-zero).

### 3.4 Direct GEMM into contribution buffer

This is the critical difference. In SPRAL's `factor_node_indef`
(`factor.hxx:92-103`), the final Schur complement update writes directly
into the contribution buffer:

```cpp
// Compute LD = L21 * D (into thread-local workspace)
calcLD<OP_N>(m-n, nelim, &lcol[nelim*ldl+n], ldl, &d[2*nelim], ld, ldld);

// GEMM: contrib = -L21 * LD^T  (output goes directly to node.contrib)
T rbeta = (nelim==0) ? 0.0 : 1.0;
host_gemm<T>(OP_N, OP_T, m-n, m-n, nelim,
    -1.0, &lcol[nelim*ldl+n], ldl, ld, ldld,
    rbeta, node.contrib, m-n);
```

The GEMM output destination is `node.contrib`, not the trailing submatrix of the
frontal workspace. There is no subsequent copy step. The Schur complement data is
written once, to its final location.

Compare with our implementation (`factor.rs:update_trailing`):

```rust
// GEMM: a22 -= W * L21^T  (output goes to trailing submatrix of frontal workspace)
tri_matmul::matmul_with_conj(a22, ..., w_ref, l21_ref.transpose(), -1.0, par);

// ... later, in extract_contribution:
// Copy trailing submatrix from workspace to contribution buffer
for j in 0..size {
    let src = &frontal_data.col_as_slice(ne + j)[ne + j..ne + j + col_len];
    data.col_as_slice_mut(j)[j..j + col_len].copy_from_slice(src);
}
```

Our GEMM writes to the workspace, and then we copy it out. SPRAL's GEMM writes
to the destination directly. For a front of 2400 rows with 26 eliminated, the
contribution is 2374^2 = 5.6M elements = 45 MB. This is 45 MB of reads + 45 MB
of writes that SPRAL avoids entirely.

### 3.5 Assembly pipeline

SPRAL splits assembly into two phases per node:

**assemble_pre** (before factorization):
1. Count incoming delays from children
2. Allocate `lcol` from AppendAlloc, `contrib` from BuddyAllocator
3. Scatter original matrix entries directly into `lcol`
4. For each child: copy delayed columns into `lcol`, scatter child
   `contrib` entries that map to fully-summed columns of this node into `lcol`

**assemble_post** (after factorization):
1. For each child: scatter remaining child `contrib` entries (those mapping
   to non-fully-summed rows) into this node's `contrib`
2. Free each child's `contrib` (returns to BuddyAllocator)

The split is significant: child contributions that map to fully-summed columns
go into `lcol` before factorization (assemble_pre), while contributions that
map to the Schur complement go into `contrib` after factorization (assemble_post).
This is possible because `lcol` and `contrib` are separate memory regions.

In our implementation, assemble_pre and assemble_post are merged: all child
contributions go into the single frontal workspace before factorization.

### 3.6 Extend-add implementation

SPRAL's extend-add uses column-oriented scatter with 4x unrolled inner loop
(`assemble.hxx:27-38`):

```cpp
void asm_col(int n, int const* idx, T const* src, T* dest) {
    int n2 = 4*(n/4);
    for(int j=0; j<n2; j+=4) {
        dest[idx[j+0]] += src[j+0];
        dest[idx[j+1]] += src[j+1];
        dest[idx[j+2]] += src[j+2];
        dest[idx[j+3]] += src[j+3];
    }
    for(int j=n2; j<n; j++)
        dest[idx[j]] += src[j];
}
```

Contributions are stored column-major without padding (leading dimension = number
of rows). The `idx` array maps child contribution row indices to parent local
indices, precomputed per child in a `Workspace`-allocated cache buffer.

Our extend-add (`extend_add_mapped` in `numeric.rs`) follows a similar pattern
but operates on the lower triangle only. The algorithmic complexity is the same.

### 3.7 Task scheduling

SPRAL uses OpenMP task DAG scheduling (`NumericSubtree.hxx:94-238`):

```cpp
#pragma omp task depend(inout: this_lcol[0:1]) depend(in: parent_lcol[0:1])
{
    assemble_pre(node);
    factor_node(node);
    assemble_post(node);
}
```

The `depend` clauses enforce bottom-up ordering without explicit barriers.
Children can execute in any order; the parent waits for all children via the
inout dependency. This is a natural fit for the tree structure and provides
automatic load balancing.

Our level-set approach achieves similar parallelism but with explicit wave
barriers. The DFS path added in this feature is purely sequential and does
not benefit from tree-level parallelism.

## 4. Architectural Comparison

| Aspect | SPRAL | Our implementation |
|--------|-------|--------------------|
| Frontal matrix | Per-node `lcol` (permanent) | Single reusable workspace |
| Factor extraction | Zero-copy (lcol IS the factor) | Copy L11/D11/L21 out |
| Contribution GEMM | Writes directly to `contrib` buffer | Writes to workspace trailing submatrix |
| Contribution extraction | No copy needed | O(n^2) copy from workspace |
| Factor allocator | AppendAlloc (bump, never freed) | Per-node `Vec` allocations in `FrontFactors` |
| Contribution allocator | BuddyAllocator (buddy-system pool) | Pool (`Vec<Mat<f64>>`) or fresh `Mat::zeros` |
| Memory usage | Higher (all nodes' lcol persist) | Lower (workspace reused, 24% RSS reduction) |
| Extend-add | Column-oriented scatter, 4x unrolled | Column-oriented scatter, precomputed maps |
| Parallelism | OpenMP task DAG | Level-set waves with rayon |

### Performance tradeoff

SPRAL's architecture uses more memory (all factor storage persists simultaneously)
but eliminates all copy overhead. The contribution copy accounts for 14-19% of
our factor time on c-71. The factor extraction copy adds another 0.3%.

For large problems, memory is typically abundant and speed is the priority. The
24% RSS reduction from workspace reuse (Phase 9.1b) is not worth the 14-19%
factor time penalty.

Additionally, SPRAL's per-node `lcol` allocations don't cause significant
cache/TLB overhead because:
- Large fronts exceed all cache levels regardless of allocation strategy
- AppendAlloc returns sequential addresses from contiguous pages (TLB-friendly)
- BuddyAllocator pre-sizes the pool to avoid growth and fragmentation

## 5. What Would Need to Change

To adopt SPRAL's architecture, the following changes would be required:

### 5.1 AppendAlloc for factor storage (moderate effort)

Replace per-node `FrontFactors` heap allocations with a bump allocator:

```rust
struct FactorAllocator {
    pages: Vec<Vec<f64>>,
    current_offset: usize,
}
```

Each node's L11, D11, L21 would be slices into the allocator's pages. The
allocator would be pre-sized from symbolic analysis (`predicted_nnz`).

**Rust consideration**: Returning slices from the allocator requires either
unsafe code (raw pointer + lifetime management) or an arena crate. The
`bumpalo` crate provides safe bump allocation but doesn't support `f64`
slices with faer's `MatRef` directly. A practical approach might be to
keep owned `Mat<f64>` per node for L/D (as now) but compute them in-place
rather than copying from a workspace.

### 5.2 Direct GEMM into contribution buffer (high impact, moderate effort)

The critical optimization. Instead of:

```
GEMM → workspace trailing submatrix → copy to contribution
```

Do:

```
GEMM → contribution buffer directly
```

This requires splitting `update_trailing` in `factor.rs` to handle the final
block specially. For intermediate blocks (within the fully-summed region), the
GEMM must still update the workspace in-place because subsequent blocks read the
accumulated updates. Only the final trailing update (rows/columns beyond
fully-summed) can target a separate contribution buffer.

Delayed columns add complexity: the contribution includes both delayed columns
(which are in the fully-summed region) and non-fully-summed rows. The final
GEMM must account for this split.

**Rust consideration**: The borrow checker makes it harder to have the GEMM write
to a separate buffer while reading from the frontal workspace. In C++, SPRAL
simply passes raw pointers. In Rust, we'd need to use submatrix views carefully
or pass raw slices. This is awkward but not impossible — `MatMut::col_as_slice_mut`
provides the necessary access patterns.

### 5.3 Split assembly into pre/post (moderate effort)

SPRAL's two-phase assembly (pre-factor and post-factor) is an elegant
consequence of having separate `lcol` and `contrib` regions. With a single
frontal workspace, we can't easily split assembly this way.

However, the split is not strictly necessary for the performance win. The
main benefit comes from the direct GEMM (5.2). The assembly split is an
optimization that allows child contributions to be freed earlier (after
assemble_pre rather than after the full factorization), reducing peak memory.

### 5.4 BuddyAllocator for contributions (low priority)

The buddy allocator is an elegant solution for contribution blocks but is not
the primary performance lever. A simpler approach — pre-allocating a single
`Mat<f64>` sized to `max_front^2` as a contribution workspace, with the GEMM
writing directly into it — would capture most of the benefit.

The buddy allocator becomes valuable when contributions from multiple children
must coexist (which they do in the parallel case). For the sequential case,
a single pre-allocated buffer suffices.

## 6. Recommended Path Forward

### Immediate: Roll back Phase 9.1d

The pool-based approach adds complexity without improving performance. The
diagnostic infrastructure (sub-phase timing, pool statistics, `--level-set`
flag) is useful and should be preserved, but the pool itself, the DFS
traversal, and the `frontal_data_alt` dead allocation should be removed.

### Next: Direct GEMM into contribution buffer

This is the highest-impact single change. It eliminates the 14-19%
`extract_contribution` cost without requiring the full AppendAlloc
refactoring. The approach:

1. Pre-allocate a contribution buffer on `FactorizationWorkspace`, sized
   to `max_contrib^2` (predictable from symbolic analysis).
2. Modify `factor_inner` / `update_trailing` so the final trailing GEMM
   (rows beyond fully-summed) writes to the contribution buffer instead
   of the frontal workspace.
3. Intermediate block updates (within fully-summed columns) continue to
   update the frontal workspace in-place.
4. `extract_contribution` becomes a no-op (or just builds `row_indices`).

This keeps the reusable frontal workspace (preserving the RSS benefit) while
eliminating the dominant copy.

### Future: AppendAlloc for factor storage

Lower priority. Replace per-node `FrontFactors` with slices into a bump
allocator. This eliminates the `extract_front_factors` copy (0.3%) and
reduces per-node heap allocation overhead. Worth doing but not urgent.

### Future: OpenMP-style task DAG scheduling

Replace level-set wave scheduling with dependency-driven task scheduling
(using rayon scoped tasks or similar). This would provide better load
balancing and avoid wave barriers. Independent of the memory architecture
changes.

## 7. Files Modified in Phase 9.1d (to be rolled back)

- `src/aptp/numeric.rs` — ContributionPool, DFS traversal, pool integration,
  ContribExtractResult, pool diagnostic fields
- `src/aptp/factor.rs` — `force_levelset` field on AptpOptions
- `src/aptp/solver.rs` — `force_levelset` field on FactorOptions
- `examples/profile_matrix.rs` — Pool waste display, `--level-set` flag

## 8. Diagnostic infrastructure to preserve

The following additions from Phase 9.1d are useful independent of the pool:

- `PerSupernodeStats.pool_buf_nrows` / `contrib_size` — buffer sizing diagnostics
- `FactorizationStats.pool_waste_elements` / `pool_max_oversize_ratio` — aggregate waste
- `--level-set` flag on `profile_matrix` — A/B testing of traversal strategies
- `force_levelset` on `AptpOptions` / `FactorOptions` — programmatic traversal control

These should be adapted for whatever allocation strategy replaces the pool.
