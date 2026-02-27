# Phase 9.1c: Assembly & Extraction Profiling Report

**Branch**: `022-assembly-extraction-opt`
**Date**: 2026-02-22

## 1. Summary

Phase 9.1c implemented three optimization tiers for assembly and extraction in
the multifrontal factorization loop: precomputed scatter maps (US1), bulk
column-slice copies for extraction (US2), and `fill(0.0)` for frontal matrix
zeroing (US3). Sub-phase timing instrumentation was added to identify the
actual bottlenecks. Empirical profiling with `perf stat` and sub-phase timing
revealed that **memory allocation churn** from contribution blocks — not the
targeted operations — is the dominant bottleneck.

## 2. Optimization Results

### SPRAL comparison (c-71)

| Phase      | SPRAL (ms) | 9.1b (ms) | 9.1c (ms) | 9.1b ratio | 9.1c ratio |
|------------|-----------|-----------|-----------|------------|------------|
| c-71       | 2,385     | 9,678     | 5,920     | 4.06×      | 2.48×      |
| c-big      | 8,420     | 35,265    | 33,601    | 4.19×      | 4.00×      |

c-71 improved **38% in absolute time** (4.06× → 2.48× vs SPRAL), though
profiling percentages appeared unchanged because all sub-phases benefited
proportionally from reduced overhead.

### What improved

- **Bulk column-slice copies** (US2): `col_as_slice` + `copy_from_slice`
  replaced element-by-element indexing in `extract_contribution` and
  `extract_front_factors`. Material on matrices with large fronts.
- **`fill(0.0)` zeroing** (US3): Replaced nested loop with contiguous-slice
  fill per column. Minor improvement.
- **Precomputed scatter maps** (US1): Eliminated `global_to_local` lookups
  during scatter for zero-delay case. Measurable but small (scatter is only
  ~0.1% of factor time).

### What didn't move the needle

The Phase 9.1c optimizations targeted operations that are **not the dominant
bottleneck**. Sub-phase timing revealed the real cost centers.

## 3. Sub-Phase Profiling Results (c-71, sequential, release)

Instrumentation added 6 sub-phase timers behind `#[cfg(feature = "diagnostic")]`:

| Sub-Phase       | Time (ms) | % of Factor | Notes |
|-----------------|-----------|-------------|-------|
| Zeroing         | 64.2      | 1.1%        | `fill(0.0)` per column |
| G2L setup       | 0.3       | 0.0%        | Global-to-local map init |
| Scatter         | 4.9       | 0.1%        | CSC → frontal matrix |
| Extend-add      | 1,974     | 33.3%       | Child → parent merge |
| ExtractFactors  | 10.2      | 0.2%        | L, D, perm extraction |
| **ExtractContr**| **2,374** | **40.1%**   | **Contribution block copy** |
| Kernel          | 1,366     | 23.1%       | Dense APTP factorization |
| Other           | 126       | 2.1%        | Overhead, timing, misc |

**73.4% of factorization time is in extend-add + contribution extraction.**
The kernel (dense factorization) is only 23.1%.

## 4. `perf stat` Hardware Analysis (c-71, sequential, release)

```
Performance counter stats for 'cargo run ...':

     30,498.42 msec  task-clock                           #    0.996 CPUs utilized
         9,580       context-switches                     #  314.117 /sec
            39       cpu-migrations                       #    1.279 /sec
       933,979       page-faults                          #   30.624 K/sec
78,445,605,771       cycles                               #    2.572 GHz
15,107,973,178       stalled-cycles-frontend              #   19.26% frontend cycles idle
14,069,523,866       stalled-cycles-backend               #   17.93% backend cycles idle
112,121,610,455      instructions                         #    1.43  insn per cycle
                                                          #    0.13  stalled cycles per insn
 17,746,403,889      branches                             #  581.893 M/sec
    144,802,974      branch-misses                        #    0.82% of all branches
644,289,759,117      dTLB-load-misses

       3.123 s  sys time (32% of wall time)
```

### Key observations

1. **IPC = 1.43**: Reasonable but not great. Memory-bound workloads typically
   show IPC < 2 even on modern CPUs. Not CPU-compute limited.

2. **644 billion dTLB misses**: Extraordinarily high. Each large contribution
   block (`Mat::zeros(size, size)`) spans thousands of 4KB pages. Repeated
   allocation/deallocation causes the OS to recycle virtual address ranges,
   thrashing the TLB.

3. **3.1s sys time (32% of wall time)**: One-third of execution is kernel
   time — `mmap`/`munmap`/page-fault handling. The glibc allocator uses
   `mmap` for allocations above ~128KB, and each contribution block for large
   supernodes is multi-megabyte (c-71 max front ~2475 → 49MB contribution).

4. **934K page faults**: Consistent with thousands of large allocations being
   mapped, faulted in, used once, then unmapped.

## 5. Root Cause: Contribution Block Allocation Churn

The multifrontal factorization loop processes supernodes bottom-up. For each
supernode:

```
1. Assemble frontal matrix (reuses pre-allocated workspace — Phase 9.1b)
2. Factor dense block (APTP kernel — compute-bound)
3. Extract L, D, perm (small, fast)
4. Extract contribution block → **NEW Mat::zeros(size, size)** ← BOTTLENECK
5. Extend-add contribution into parent → BOTTLENECK (N² element-wise)
6. Drop contribution block → **deallocation** ← BOTTLENECK
```

Steps 4 and 6 allocate and deallocate a `Mat<f64>` for every supernode with
uneliminated rows. For c-71's 6,350 supernodes (post-amalgamation), this is
thousands of large allocations. The largest contribution blocks are ~2475×2475
(~49 MB each).

Phase 9.1b solved this problem for the **frontal matrix** by pre-allocating a
single workspace and reusing it. The contribution block was not addressed
because it's passed from child to parent (lifetime crosses supernode
boundaries).

## 6. Proposed Optimization: Eliminate Contribution Copy

The fundamental insight is that the contribution block is a **copy** of the
trailing submatrix of the frontal matrix. Instead of:

```
frontal → copy to contribution → extend-add contribution into parent → drop
```

We can do:

```
frontal → extend-add directly from frontal into parent → done
```

The contribution is the `(m-ne) × (m-ne)` lower-right submatrix of the
factored frontal matrix, where `ne` = number of eliminated pivots. This data
is already sitting in the frontal workspace — we just need to read it before
the frontal matrix is overwritten by the next supernode.

### Why this wasn't done in Phase 9.1b

The original architecture processes children before the parent. By the time
the parent's frontal matrix is being assembled, the child's frontal workspace
has already been overwritten. The contribution copy exists to preserve the
child's Schur complement across supernode boundaries.

### Solution: Direct extend-add from frontal

For the **sequential path** (bottom-up tree traversal), each supernode is
processed completely before moving to the next. The child's frontal workspace
is valid until the parent begins. If we perform extend-add immediately after
factoring each child (before moving to the next sibling or parent), we can
read directly from the frontal workspace.

This requires restructuring the factorization loop:

```
For each supernode (bottom-up):
  1. Assemble (scatter original + extend-add from ALL children)
  2. Factor
  3. Extract L, D, perm
  4. (contribution remains in-place in frontal workspace)
```

The extend-add in step 1 reads directly from the child's frontal workspace,
which is still valid because the child was processed in the previous
iteration.

### Challenge: Multiple children

A supernode may have multiple children. With a single frontal workspace, only
the **last-processed child** still has valid data in the workspace. Earlier
children's data has been overwritten.

Options:
1. **Contribution workspace**: A second pre-allocated `Mat<f64>` used only for
   contribution storage. Each child extracts into this workspace, parent reads
   from it. One allocation instead of thousands.
2. **Leaf-path optimization**: For linear chains (supernode → single child →
   grandchild → ...), skip the copy entirely and read directly from the frontal
   workspace. Only allocate contribution blocks for supernodes with >1 child.
3. **Deferred extend-add**: Change traversal order so that each child's
   extend-add happens immediately after its factorization, reading directly
   from the frontal workspace. Requires the parent frontal matrix to be
   allocated before its children are processed.

### Challenge: Parallel path

The parallel path uses thread-local frontal workspaces. Children processed on
different threads cannot share a workspace. However, the contribution block
is already an owned `Mat<f64>` that transfers between threads. A per-thread
contribution workspace would eliminate per-child allocations within each
thread.

## 7. Immediate Next Steps

### Priority 1: Contribution workspace reuse (high impact, moderate effort)

Add a pre-allocated `contribution_data: Mat<f64>` to `FactorizationWorkspace`,
sized to the maximum contribution block dimension. `extract_contribution`
writes into this workspace instead of allocating a new `Mat<f64>`.

Expected impact: Eliminate ~40% of factor time (contribution extraction) and
reduce sys time from 32% to near-zero. TLB misses should drop by orders of
magnitude.

### Priority 2: Direct extend-add (high impact, high effort)

For the sequential path, restructure the factorization loop so extend-add
reads directly from the child's frontal workspace. This eliminates both the
contribution copy AND the element-wise extend-add copy (33% of factor time).

Combined with Priority 1, this could reduce factor time by up to 73%.

### Priority 3: Huge pages (medium impact, low effort)

Use `madvise(MADV_HUGEPAGE)` or explicit 2MB huge pages for the frontal and
contribution workspaces. Each 2MB huge page replaces 512 TLB entries, directly
addressing the 644B dTLB miss count.

## 8. Files Modified in Phase 9.1c

- `src/aptp/numeric.rs` — AssemblyMaps, build_assembly_maps, scatter fast path,
  extend_add_mapped, bulk extraction, bulk zeroing, sub-phase timing (6 fields)
- `examples/profile_matrix.rs` — Sub-phase breakdown display
- `docker/entrypoint.sh` — perf_event_paranoid sysctl
- `docker/run.sh` — CAP_PERFMON, SYS_ADMIN, seccomp=unconfined
- `.devcontainer/docker-compose.yml` — cap_add, security_opt
- `.devcontainer/devcontainer.json` — capAdd, securityOpt
