# Research: Contribution Workspace Reuse

**Feature**: 023-contrib-workspace-reuse
**Date**: 2026-02-23

## R1: How does SPRAL manage contribution block memory?

**Decision**: Use a contribution buffer pool (free-list of reusable `Mat<f64>`) — simpler than SPRAL's buddy allocator but addresses the same root cause.

**Findings**: SPRAL uses a two-allocator architecture:
- **FactorAllocator (AppendAlloc)**: Sequential append-only pool for persistent factor storage (L, D, permutations). Uses `calloc()` for zero-initialized pages. Never freed during factorization.
- **PoolAllocator (BuddyAllocator)**: Buddy system (binary tree splitting, free-list coalescing) for transient contribution blocks. Allocates per-node, freed immediately after assembly.

SPRAL does NOT reuse a single pre-allocated workspace for contributions. Instead, it allocates per-supernode from the BuddyAllocator, which efficiently recycles freed memory without OS mmap/munmap overhead. The buddy allocator handles variable-size blocks with O(log N) allocation and automatic coalescing.

Key SPRAL pattern (`NumericNode.hxx:37-51`):
- `alloc_contrib()`: Allocate from BuddyAllocator before assembly
- `free_contrib()`: Return to BuddyAllocator immediately after parent consumes it

**Rationale**: A contribution buffer pool (free-list) achieves the same effect as SPRAL's BuddyAllocator — memory is allocated from the OS only once and recycled. A simple free-list is sufficient because rivrs-sparse can size buffers to max_front (the workspace is already max_front × max_front). The buddy allocator's variable-size splitting is unnecessary when all buffers are the same maximum size.

**Alternatives considered**:
- Full buddy allocator: More complex, handles variable sizes elegantly but unnecessary given our max_front workspace sizing strategy
- Arena allocator: Would require knowing total contribution memory upfront; not practical with dynamic delayed columns
- Single pre-allocated contribution buffer: Lifetime conflict — contribution must survive until parent consumes it, but workspace is reused for next sibling

## R2: Can direct extend-add from frontal workspace work with the current traversal?

**Decision**: Restructure the sequential path from wave-based level-set to stack-based DFS postorder to enable direct extend-add for the last child of each parent.

**Findings**: The current `factor_tree_levelset` uses wave-based BFS: all leaves first, then their parents, etc. With this traversal, by the time a parent is processed, ALL children's frontal workspaces have been overwritten by subsequent supernodes. Direct extend-add is impossible because the child's data is gone.

Direct extend-add requires the child's frontal workspace to be valid when the extend-add happens. This is only possible if the extend-add happens immediately after the child is factored, before the workspace is reused.

In a DFS postorder traversal, the last child processed before a parent is always the parent's LAST child in the children list (because DFS processes the entire subtree rooted at each child before moving to the next sibling). This child's frontal data is still in the workspace when the parent starts processing.

**The dual-buffer mechanism**: Two frontal buffers (A and B) are needed so the parent's frontal matrix can be assembled in one buffer while the child is factored in the other. The child factors in buffer A, the parent's assembly accumulates in buffer B. After factoring the last child, extend-add reads directly from buffer A into buffer B. The parent then factors in buffer B.

**Multi-child handling**: For a parent with N children:
- Children 1..N-1: Extract contribution into pooled buffer (R1), immediately extend-add into parent's buffer B, return pool buffer
- Child N (last): Direct extend-add from buffer A into buffer B (zero-copy)

**Rationale**: DFS postorder is the natural traversal for the dual-buffer approach. It's also how SPRAL processes subtrees (`factor_subtree` in `factor_cpu.cxx`). The wave-based approach was chosen for parallelism (batching independent nodes), but the sequential path doesn't benefit from batching.

**Alternatives considered**:
- Keep wave-based traversal, use contribution workspace only: Simpler but misses the direct extend-add optimization (33% of factor time)
- Wave-based traversal with parent pre-allocation: Requires multiple parent buffers simultaneously (one per parent in the wave), which defeats the dual-buffer approach
- Recursive DFS: Natural but risks stack overflow on deep trees; iterative DFS with explicit stack is equivalent

## R3: How should the dual-buffer alternation work across tree depths?

**Decision**: Ping-pong between buffers based on whether we're factoring a child or its parent. Use a simple "active buffer" index that flips when moving from child processing to parent processing.

**Findings**: In a DFS postorder traversal, the processing sequence for a subtree rooted at P with children C1, C2 is:
```
[subtree(C1)], [subtree(C2)], P
```

Within `[subtree(C1)]`, the same dual-buffer ping-pong happens recursively. The key invariant: when we return from processing a subtree, the subtree root's factored data is in buffer X, and the parent should use buffer Y (the other buffer).

For the sequential path, this is managed by an explicit stack of "which buffer am I using":
- Push: entering a subtree, allocate parent's frontal in the other buffer
- Pop: leaving a subtree, parent's frontal is assembled, ready to factor

For multi-child parents where earlier children use pooled ContributionBlocks:
- Factor C1 in buffer A → extract to pool → extend-add pool into parent's buffer B → return pool buffer
- Factor C2 in buffer A → direct extend-add from buffer A into buffer B
- Factor P in buffer B → extract to pool (for P's parent) or direct extend-add into P's parent's buffer A

**Rationale**: The ping-pong approach is simple and correct. It requires only two buffers regardless of tree depth. Each buffer alternates between "factoring target" and "assembly accumulator" roles.

## R4: Delayed columns and parent pre-allocation timing

**Decision**: Pre-allocate the parent's frontal matrix at maximum estimated size (no delays). If delayed columns from children increase the actual size, grow the parent buffer after collecting all children's delay information.

**Findings**: The parent's front size `m = sn_ncols + delayed_from_children + pattern.len()` depends on how many columns children delayed. With MatchOrderMetis ordering, delays are rare (c-71: ~1 delay out of 6,350 supernodes with MatchOrderMetis). So the "estimated size" (without delays) is almost always exact.

For the rare delay case:
1. Pre-allocate parent buffer at estimated size (sn_ncols + pattern.len())
2. Process children, collecting delayed columns
3. If delays > 0: either grow buffer (if insufficient) or fall back to existing assembly pattern
4. The pre-allocation is not wasted — it's in a reusable workspace buffer

**Rationale**: Optimizing for the common case (no delays) and falling back for the rare case is the right tradeoff. The alternative — not pre-allocating until all children are done — would prevent direct extend-add entirely.

## R5: Parallel path strategy

**Decision**: Add contribution buffer pool to per-thread `FactorizationWorkspace`. Direct extend-add is NOT applicable for the parallel path (contributions cross thread boundaries). The parallel path uses the wave-based traversal unchanged.

**Findings**: In the parallel path:
- Tree-level parallelism: siblings processed on different rayon threads
- Each thread has a `Cell<FactorizationWorkspace>` via thread-local storage
- Contributions are owned `ContributionBlock` values transferred between threads via the `contributions[s]` vector

Direct extend-add requires the child's frontal workspace to be valid when the parent processes, which can't be guaranteed when child and parent are on different threads. The parallel path must continue using owned ContributionBlocks.

However, the contribution buffer pool still helps: each thread's workspace has its own pool, eliminating per-supernode allocation for all contributions produced within that thread's work items.

**Rationale**: The parallel path already works correctly with owned ContributionBlocks. Adding a per-thread pool addresses the allocation overhead without the complexity of cross-thread buffer management.

## R6: Liu (1992) stack-based workspace management

**Decision**: The DFS postorder traversal with LIFO contribution management follows Liu's Algorithm 6.1 directly.

**Findings**: Liu (1992), Section 6 describes the supernodal multifrontal method where update matrices (contribution blocks) flow through a LIFO stack:
1. For each supernode in postorder: pop children's update matrices from stack, extend-add into frontal matrix, factor, push update matrix onto stack
2. Stack discipline matches DFS postorder naturally — last child's update is on top

Section 7.4 discusses stack storage optimization: optimal child ordering minimizes peak stack depth. The formula `MinWstore(j) = max_child { max{MinWstore(c_k), |F_j|} + sum_{i<k} |U_{c_i}| }` quantifies the memory-throughput tradeoff.

**Rationale**: Our DFS postorder + dual-buffer approach is a modern realization of Liu's stack-based architecture, with the optimization that the "top of stack" (last child's update) doesn't need to be copied at all — it's read directly from the frontal workspace.

## R7: Duff (2020) SSIDS memory management

**Decision**: Adopt SSIDS's two-tier approach conceptually: persistent factor storage (existing `FrontFactors` with owned Mats) + transient contribution blocks (pool-allocated, immediately recycled).

**Findings**: Duff (2020), Section 5 describes SSIDS's memory management:
- Factors allocated via `calloc`-backed stack allocator (zero-on-alloc, never freed during factorization)
- Contribution blocks allocated via buddy system allocator (variable-size, freed after assembly)
- OS page cache exhaustion identified as performance bottleneck for large problems
- Contribution block zeroing folded into GEMM (`beta=0.0`) to avoid explicit zero loop

Key insight: the contribution block computation (trailing GEMM update) can zero-initialize as part of the GEMM call, avoiding a separate zeroing pass. This is relevant for our `extract_contribution` which currently zeroes via `Mat::zeros` then copies.

**Rationale**: Our buffer pool serves the same role as SSIDS's buddy allocator — recycling transient memory without OS interaction. The GEMM-based zeroing optimization is a potential follow-up but not in scope for this feature.

## R8: Root cause analysis — why rivrs-sparse allocation is 40% of factor time

**Decision**: The performance gap is dominated by **first-touch page faults** during faer's `Mat::zeros` zero-write pass on freshly mapped memory, NOT by mmap/munmap syscall overhead. The buffer pool approach addresses this by keeping pages physically resident across supernodes.

**Findings**: Diagnostic testing (2026-02-23) with custom allocator instrumentation and `MALLOC_MMAP_THRESHOLD_` override revealed the actual cost breakdown.

### Diagnostic results: mmap/munmap syscalls are NOT the bottleneck

Direct measurement with a custom global allocator wrapper timing every allocation ≥ 128KB:

| Matrix | Large Allocs | Total Bytes | alloc+dealloc time | % of Factor Time |
|--------|-------------|------------|-------------------|-----------------|
| SiNa (n=5743, max_front=2089) | 340 | 439 MB | 2.1 ms | 0.6% |
| bratu3d (n=27792, max_front=1496) | 425 | 298 MB | 0.84 ms | 0.3% |
| vibrobox (n=12328, max_front=882) | 182 | 98 MB | ~1.5 ms | ~1.8% |

The mmap/munmap syscalls themselves cost < 2ms total — negligible compared to the 50-128ms spent in ExtractContr.

`MALLOC_MMAP_THRESHOLD_=100000000` (forcing heap allocation) showed **no improvement** on any matrix:
- SiNa: median 305ms (default) → 335ms (threshold override) — *worse* due to heap fragmentation
- bratu3d: median 252ms → 251ms — no change
- vibrobox: median 83ms → 79ms — marginal improvement

### The actual bottleneck: first-touch page faults during zero-write + data copy

When glibc mmaps fresh pages, the `alloc()` call just sets up virtual memory mappings (cheap). But faer's `Mat::zeros(n, n)` then writes zeros to every element via `from_fn(|_, _| T::zero_impl())`, triggering a minor page fault per 4KB page on first touch. For a 10MB contribution block, that's ~2,500 page faults. This cost is incurred *during the write loop*, not during the alloc call — which is why the allocator wrapper measured only 2ms while ExtractContr measured 128ms.

Isolated buffer reuse benchmark confirms:

| Contribution Size | Fresh alloc+zero+copy | Buffer reuse+zero+copy | Speedup |
|-------------------|----------------------|------------------------|---------|
| 500×500 (2 MB) | 0.27 ms | 0.22 ms | 1.3x |
| 800×800 (5 MB) | 1.25 ms | 0.55 ms | 2.3x |
| 1000×1000 (8 MB) | 1.65 ms | 0.97 ms | 1.7x |
| 1200×1200 (12 MB) | 2.20 ms | 1.17 ms | 1.9x |
| 1500×1500 (18 MB) | 2.41 ms | 1.36 ms | 1.8x |

Reused buffers are **1.7-2.3x faster** because pages are already resident in the process page tables — no page faults on subsequent use.

### The dominant supernodes

A handful of supernodes with large fronts but few eliminated columns dominate ExtractContr time:

| Supernode | Front | Elim | Contrib Size | Contrib MB | ext_ms |
|-----------|-------|------|-------------|-----------|--------|
| 104 | 1199 | 62 | 1137 | 10.3 | 26.7 |
| 56 | 1146 | 82 | 1064 | 9.1 | 26.9 |
| 68 | 1315 | 101 | 1214 | 11.8 | 15.7 |

These ~5 supernodes produce contribution blocks > 10MB, each incurring ~2,500+ page faults. They account for the majority of ExtractContr time on SiNa.

### Revised cost model

| Factor | Estimated Cost | Previously Claimed |
|--------|---------------|-------------------|
| mmap/munmap syscalls | < 0.6% of factor time | "dominant bottleneck" |
| First-touch page faults (Mat::zeros) | **~15-25% of factor time** | "compounding factor" |
| TLB shootdowns | Negligible (consistent with low syscall cost) | "significant contributor" |
| Data copy (column-by-column scatter) | **~10-15% of factor time** | Not separately identified |

The 40.1% ExtractContr cost breaks down as: ~60% first-touch page faults during zeroing, ~40% data copy (column-by-column scatter from frontal matrix into contribution block). Both are addressed by the pool (page-fault-free reuse) and direct extend-add (eliminates copy entirely).

### Why SPRAL doesn't have this problem

SPRAL's BuddyAllocator keeps ALL contribution memory resident in a single contiguous virtual address range:
- One `calloc` at initialization, one `free` at cleanup
- Per-supernode: user-space pointer arithmetic only
- Pages stay physically allocated and TLB-warm between supernodes
- No first-touch page faults — pages are already resident from initial `calloc`

### Alternative quick fixes ruled out

1. **`MALLOC_MMAP_THRESHOLD_=100000000`**: Tested — no improvement. Heap fragmentation offsets any benefit from avoiding mmap. Not a viable workaround.

2. **jemalloc / mimalloc**: Would use `madvise(MADV_DONTNEED)` instead of `munmap` (retaining virtual mapping). But physical pages are still released and would re-fault on next use. Would NOT address the first-touch page fault bottleneck.

3. **faer upstream change to use `alloc_zeroed`/`calloc`**: Would let the kernel defer physical page allocation until non-zero writes, but doesn't help when we immediately write data into the buffer. Out of our control regardless.

**Rationale**: The contribution buffer pool is the correct optimization because it keeps pages **physically resident** across supernodes, eliminating the dominant first-touch page fault cost. Direct extend-add (Layer 2) is even more impactful — it eliminates both the zero-write AND the data copy for the last child, addressing the full 40% ExtractContr cost for those supernodes.
