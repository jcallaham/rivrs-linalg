# SPRAL SSIDS Code Review

An annotated architecture review of SPRAL's Sparse Symmetric Indefinite Direct
Solver (SSIDS) implementation. This document is designed so that a developer
unfamiliar with SPRAL can understand the three-phase architecture
(analyse, factor, solve) and where APTP logic lives within 30 minutes.

All file paths are relative to `references/spral/`.

---

## 1. Executive Summary

SSIDS solves sparse symmetric systems Ax = b via a multifrontal method. The
matrix A may be positive definite (Cholesky: A = PL(PL)^T) or indefinite
(LDL^T: A = PLD(PL)^T, where D is block diagonal with 1x1 and 2x2 blocks).
The solver is organized as a strict four-step pipeline: analyse, factor, solve,
free.

**Analyse (symbolic factorization).** The user supplies the matrix structure in
CSC format. SSIDS optionally cleans the data (removing duplicates and
out-of-range entries), computes a fill-reducing ordering (METIS by default, or
matching-based ordering), expands the lower-triangular structure to full,
constructs the supernodal elimination tree, and partitions the tree into
subtrees assigned to NUMA regions or GPUs. The output is an `ssids_akeep`
object that encodes the symbolic factorization and can be reused across
multiple numeric factorizations sharing the same sparsity pattern.

**Factor (numeric factorization).** Given the symbolic factorization in `akeep`
and the numerical values of A, SSIDS performs the numeric factorization. Leaf
subtrees are factored in parallel on their assigned CPU NUMA regions (or GPUs);
then all CPU resources cooperate on the root subtree. On the CPU, each
supernode is factored using either an a posteriori threshold pivoting (APTP)
kernel for indefinite systems or a Cholesky kernel for positive-definite
systems. Failed pivots are either retried with traditional threshold partial
pivoting (TPP) or delayed to parent nodes as contribution blocks. The output is
an `ssids_fkeep` object that stores the numeric factors.

**Solve (triangular solve).** Using the stored factors, SSIDS performs forward
substitution (PL), diagonal solve (D), and backward substitution ((PL)^T) to
compute x = A^{-1}b. The user may request partial solves (forward only,
diagonal only, backward only). Scaling is applied and reversed around the
solve. The solve phase is currently serial.

**Sources:**
`docs/Fortran/ssids.rst`,
`examples/Fortran/ssids.f90`,
`src/ssids/ssids.f90`

---

## 2. Module-by-Module Overview

The table below covers the top-level Fortran modules and the analyse layer.
CPU-specific files are covered separately in Section 3.

| File | Lines | Responsibility | Key Types / Exports |
|------|------:|----------------|---------------------|
| `src/ssids/ssids.f90` | 1479 | Top-level public API module. Defines generic interfaces for `ssids_analyse`, `ssids_factor`, `ssids_solve`, `ssids_free`, `ssids_enquire_*`, `ssids_alter`. Dispatches to internal routines. | `spral_ssids` module; re-exports `ssids_akeep`, `ssids_fkeep`, `ssids_options`, `ssids_inform` |
| `src/ssids/datatypes.f90` | 374 | Shared constants, error/warning codes, solver options, and internal data types used across all phases. | `ssids_options`, `node_type`, `smalloc_type`, `stack_type`, `thread_stats`, pivot method constants (`PIVOT_METHOD_APP_BLOCK`, etc.) |
| `src/ssids/akeep.f90` | 120 | Defines the symbolic factorization data structure returned by the analyse phase. Contains the elimination tree, supernode pointers, row lists, inverse permutation, subtree partitions, and cleaned matrix data. | `ssids_akeep` (with fields: `sptr`, `sparent`, `rptr`, `rlist`, `invp`, `nptr`, `nlist`, `part`, `subtree`, `contrib_ptr`, `contrib_idx`) |
| `src/ssids/fkeep.F90` | 463 | Defines the numeric factorization data structure. Contains the factored subtrees, scaling, and methods for `inner_factor`, `inner_solve`, `enquire_*`, and `alter`. CPU-only build; dispatches to subtree factor/solve routines. | `ssids_fkeep` (with fields: `scaling`, `pos_def`, `subtree` array, `inform`) |
| `src/ssids/inform.f90` | 216 | Information/statistics type returned to user. Tracks flags, pivot counts, flop counts, delays, rank, and more. Provides `reduce()` for combining stats from parallel threads. | `ssids_inform` (with fields: `flag`, `num_delay`, `num_factor`, `num_flops`, `num_neg`, `num_two`, `matrix_rank`, `maxfront`, `maxsupernode`) |
| `src/ssids/contrib.f90` | 77 | Contribution block type for passing uneliminated entries between subtrees. Contains the dense lower-triangular contribution matrix, row list, and delayed pivot data. | `contrib_type` (with fields: `n`, `val`, `ldval`, `rlist`, `ndelay`, `delay_perm`, `delay_val`, `owner`) |
| `src/ssids/subtree.f90` | 128 | Abstract base classes for symbolic and numeric subtrees. Defines the interface contract that CPU and GPU subtree implementations must satisfy: `factor()`, `get_contrib()`, `solve_fwd/diag/diag_bwd/bwd()`, `cleanup()`. | `symbolic_subtree_base` (abstract), `numeric_subtree_base` (abstract) |
| `src/ssids/anal.F90` | 1240 | Core analyse phase implementation. Calls `spral_core_analyse::basic_analyse` for the elimination tree and supernode construction, partitions the tree into subtrees for parallel factorization, assigns subtrees to NUMA regions and GPUs, and builds the mapping from original matrix entries to factor storage. | `analyse_phase`, `check_order`, `expand_pattern`, `expand_matrix` |

**Sources:**
`src/ssids/ssids.f90`,
`src/ssids/datatypes.f90`,
`src/ssids/akeep.f90`,
`src/ssids/fkeep.F90`,
`src/ssids/inform.f90`,
`src/ssids/contrib.f90`,
`src/ssids/subtree.f90`,
`src/ssids/anal.F90`

---

## 3. CPU Factorization Stack

The CPU factorization is implemented as a mixed Fortran/C++ stack. The Fortran
`cpu_subtree` module provides the Fortran-facing wrapper, while the core
factorization logic lives in C++ for performance and template flexibility.

### 3.1 CPU Subtree Layer

| File | Lines | Responsibility | Key Types |
|------|------:|----------------|-----------|
| `src/ssids/cpu/subtree.f90` | 425 | Fortran wrapper around C++ subtree. Defines `cpu_symbolic_subtree` and `cpu_numeric_subtree` extending the abstract base classes. Provides C interop bindings to create/destroy subtrees and dispatch solve operations. | `cpu_symbolic_subtree`, `cpu_numeric_subtree`, `construct_cpu_symbolic_subtree` |
| `src/ssids/cpu/cpu_iface.f90` | 175 | Fortran-C interop layer. Defines `cpu_factor_options` and `cpu_factor_stats` as `bind(C)` types. Provides wrapper subroutines for BLAS/LAPACK calls (`spral_c_dgemm`, `spral_c_dpotrf`, `spral_c_dsytrf`, `spral_c_dtrsm`, `spral_c_dsyrk`, `spral_c_dtrsv`, `spral_c_dgemv`). | `cpu_factor_options`, `cpu_factor_stats` |
| `src/ssids/cpu/cpu_iface.hxx` | 51 | C++ side of the interop. Defines `PivotMethod` enum (`app_aggressive`, `app_block`, `tpp`), `FailedPivotMethod` enum (`tpp`, `pass`), and `cpu_factor_options` struct. Also provides the `align_lda<T>()` template for AVX alignment. | `PivotMethod`, `FailedPivotMethod`, `cpu_factor_options` |

### 3.2 C++ Core

| File | Lines | Responsibility | Key Types |
|------|------:|----------------|-----------|
| `src/ssids/cpu/SymbolicNode.hxx` | 28 | POD struct for a single symbolic node: row/column counts, row list pointer, child linked list, parent index, assembly map. | `SymbolicNode` |
| `src/ssids/cpu/SymbolicSubtree.hxx` | 111 | Symbolic factorization of a CPU subtree. Builds the node linked-list tree, identifies small leaf subtrees (below `small_subtree_threshold` flops), and precomputes factor memory estimates. | `SymbolicSubtree` |
| `src/ssids/cpu/NumericNode.hxx` | 78 | Per-node numeric data: `nelim`, `ndelay_in`, `ndelay_out`, `lcol` (factor values), `perm`, `contrib` (contribution block), child/sibling linked-list pointers. | `NumericNode<T, PoolAllocator>` |
| `src/ssids/cpu/NumericSubtree.hxx` | 550 | Main numeric factorization driver for a CPU subtree. Allocates factor storage, iterates over nodes bottom-up, assembles contributions from children, calls `factor_node()` for each node, and builds the contribution block for the parent subtree. Uses OpenMP tasks for parallelism. | `NumericSubtree<posdef, T, PAGE_SIZE, FactorAlloc>` |
| `src/ssids/cpu/NumericSubtree.cxx` | 262 | Explicit template instantiations and C-callable entry points for `NumericSubtree`. Provides the `spral_ssids_cpu_create_num_subtree_dbl` function called from Fortran. | (instantiations) |
| `src/ssids/cpu/SmallLeafSymbolicSubtree.hxx` | 161 | Symbolic information for small leaf subtrees that are factored as a single serial task to reduce OpenMP scheduling overhead. | `SmallLeafSymbolicSubtree` |
| `src/ssids/cpu/SmallLeafNumericSubtree.hxx` | 507 | Serial factorization of small leaf subtrees. Mirrors `NumericSubtree` logic but without OpenMP task generation. | `SmallLeafNumericSubtree<posdef, T, FactorAlloc, PoolAlloc>` |
| `src/ssids/cpu/factor.hxx` | 187 | Node-level factorization dispatch. `factor_node_indef()` calls `ldlt_app_factor()` (APTP) and falls back to `ldlt_tpp_factor()` for failed pivots. `factor_node_posdef()` calls `cholesky_factor()`. | `factor_node_indef`, `factor_node_posdef`, `factor_node` |
| `src/ssids/cpu/ThreadStats.hxx` | 64 | Per-thread statistics accumulator: flag, num_delay, num_factor, num_flops, num_neg, num_two, num_zero, maxfront, not_first_pass, not_second_pass. | `ThreadStats` |

### 3.3 CPU Kernels

| File | Lines | Responsibility | Key Types / Functions |
|------|------:|----------------|----------------------|
| `src/ssids/cpu/kernels/ldlt_app.cxx` | 2591 | **Core APTP kernel.** Implements the a posteriori threshold pivoting LDL^T factorization. Contains `Column<T>` class, `ColumnData<T>`, `CopyBackup<T>`, `PoolBackup<T>`, `check_threshold()`, `apply_pivot()`, and the main `LDLT<...>::factor()` static method which dispatches to `run_elim_unpivoted` (aggressive) or `run_elim_pivoted` (block). See Section 5 for details. | `Column<T>`, `ColumnData<T>`, `ldlt_app_factor()`, `check_threshold()`, `apply_pivot()` |
| `src/ssids/cpu/kernels/ldlt_app.hxx` | 26 | Header for APTP kernel. Declares `ldlt_app_factor()`, `ldlt_app_solve_fwd()`, `ldlt_app_solve_diag()`, `ldlt_app_solve_bwd()`. | (declarations) |
| `src/ssids/cpu/kernels/ldlt_tpp.cxx` | 317 | Traditional Threshold Partial Pivoting (TPP) LDL^T kernel. Used as fallback when APTP fails, or when `pivot_method=3`. Serial, but numerically robust. | `ldlt_tpp_factor()` |
| `src/ssids/cpu/kernels/ldlt_nopiv.cxx` | 114 | Unpivoted LDL^T kernel for positive semi-definite blocks where pivoting is unnecessary. | `ldlt_nopiv_factor()` |
| `src/ssids/cpu/kernels/cholesky.cxx` | 215 | Blocked Cholesky factorization kernel using DPOTRF/DTRSM/DSYRK. | `cholesky_factor()` |
| `src/ssids/cpu/kernels/block_ldlt.hxx` | 415 | Block-level LDL^T factorization of small diagonal blocks (INNER_BLOCK_SIZE = 32). Includes Bunch-Kaufman style 1x1/2x2 pivot selection with threshold test. | `block_ldlt()`, `swap_cols()`, `find_maxloc()` |
| `src/ssids/cpu/kernels/assemble.hxx` | 445 | Assembly (scatter/gather) routines for building frontal matrices from child contributions and original matrix entries. | `assemble_node()`, `assemble_contrib()` |
| `src/ssids/cpu/kernels/calc_ld.hxx` | 120 | Computes L*D product needed for Schur complement updates. | `calcLD<op>()` |
| `src/ssids/cpu/kernels/common.hxx` | 143 | Shared enums (`operation`, `side`, `fill_mode`, `diagonal`) and small utility functions for the kernel layer. | `OP_N`, `OP_T`, `SIDE_LEFT`, `FILL_MODE_LWR`, etc. |
| `src/ssids/cpu/kernels/wrappers.cxx` | 92 | Thin C++ wrappers around BLAS/LAPACK calls (`host_gemm`, `host_trsm`, `host_syrk`, `host_trsv`, `gemv`). | `host_gemm<T>()`, `host_trsm<T>()`, `host_syrk<T>()` |
| `src/ssids/cpu/kernels/verify.hxx` | 183 | Debug verification routines for checking factorization correctness. Not used in production builds. | `Verify<T>` |
| `src/ssids/cpu/kernels/SimdVec.hxx` | 200 | SIMD vector abstraction for SSE2/AVX/AVX-512 used by `block_ldlt` for vectorized max-location searches. | `SimdVec<T>` |

**Sources:**
`src/ssids/cpu/subtree.f90`,
`src/ssids/cpu/cpu_iface.f90`,
`src/ssids/cpu/cpu_iface.hxx`,
`src/ssids/cpu/NumericSubtree.hxx`,
`src/ssids/cpu/NumericNode.hxx`,
`src/ssids/cpu/SymbolicSubtree.hxx`,
`src/ssids/cpu/factor.hxx`,
`src/ssids/cpu/kernels/ldlt_app.cxx`,
`src/ssids/cpu/kernels/ldlt_tpp.cxx`,
`src/ssids/cpu/kernels/cholesky.cxx`,
`src/ssids/cpu/kernels/block_ldlt.hxx`,
`src/ssids/cpu/kernels/assemble.hxx`

---

## 4. Data Flow Diagram

The following diagram shows data flow through the three phases, key data
structures created at each stage, and how information is passed between phases.

```
USER INPUT                              ANALYSE PHASE
+------------------+                    +------------------------------------+
| A (CSC format)   |    ssids_analyse   |                                    |
| - n, ptr, row    | =================>| 1. Clean/validate matrix data      |
| - val (optional) |                    | 2. Compute ordering (METIS/user)   |
| - order (opt.)   |                    | 3. Expand lower tri to full        |
| - options        |                    | 4. basic_analyse():                |
+------------------+                    |    - Elimination tree              |
                                        |    - Supernode detection           |
                                        |    - Supernode amalgamation        |
                                        | 5. Partition tree into subtrees    |
                                        | 6. Assign subtrees to NUMA/GPU    |
                                        | 7. Build A-to-L mapping (nlist)   |
                                        +------------------------------------+
                                                       |
                                                       v
                                        +------------------------------------+
                                        |         ssids_akeep                |
                                        | ---------------------------------- |
                                        | sptr, sparent  (elimination tree)  |
                                        | rptr, rlist    (row indices)       |
                                        | invp           (inverse perm)      |
                                        | nptr, nlist    (A-to-L map)        |
                                        | part, subtree  (tree partition)    |
                                        | contrib_ptr/idx (child contribs)   |
                                        | topology       (machine layout)    |
                                        +------------------------------------+
                                                       |
                                                       v
USER INPUT                              FACTOR PHASE
+------------------+                    +------------------------------------+
| A values (val)   |    ssids_factor    |                                    |
| - posdef flag    | =================>| 1. Optional scaling (MC64/Auction) |
| - options        |                    | 2. For each subtree partition:     |
| - scale (opt.)   |                    |    a. Allocate factor storage      |
+------------------+                    |    b. Bottom-up over nodes:        |
                                        |       - Assemble from A + children |
                                        |       - factor_node():            |
                                        |         posdef => cholesky_factor  |
                                        |         indef  => ldlt_app_factor  |
                                        |                   (APTP kernel)    |
                                        |         fallback => ldlt_tpp_factor|
                                        |       - Form contribution block   |
                                        |    c. Pass contrib to parent       |
                                        | 3. Leaf subtrees: parallel (OMP)  |
                                        | 4. Root subtree: all CPUs coop.   |
                                        +------------------------------------+
                                                       |
                                                       v
                                        +------------------------------------+
                                        |         ssids_fkeep                |
                                        | ---------------------------------- |
                                        | subtree[]  (numeric_subtree_ptr)   |
                                        |   -> nodes[].lcol  (L factors)     |
                                        |   -> nodes[].perm  (pivot order)   |
                                        |   -> nodes[].nelim (elim count)    |
                                        |   -> d[] stored after lcol in each |
                                        |      node (D^{-1} entries)         |
                                        | scaling[]  (diagonal scaling)      |
                                        | pos_def    (factorization type)    |
                                        +------------------------------------+
                                                       |
                                                       v
USER INPUT                              SOLVE PHASE
+------------------+                    +------------------------------------+
| b (rhs vector)   |    ssids_solve     |                                    |
| - job (opt.)     | =================>| 1. Permute & scale x:              |
+------------------+                    |    x2(i) = x(invp(i)) * scaling(i) |
                                        | 2. Forward solve (job=0,1):        |
                                        |    for part = 1..nparts:           |
                                        |      subtree.solve_fwd(x2)        |
                                        |    (L x2 = b, walking tree up)     |
                                        | 3. Diagonal solve (job=0,2,4):     |
                                        |    for part = 1..nparts:           |
                                        |      subtree.solve_diag(x2)       |
                                        |    (D x2 = b)                      |
                                        | 4. Backward solve (job=0,3,4):     |
                                        |    for part = nparts..1:           |
                                        |      subtree.solve_bwd(x2)        |
                                        |    (L^T x2 = b, walking tree down) |
                                        | 5. Unpermute & unscale:            |
                                        |    x(invp(i)) = x2(i) * scaling(i)|
                                        +------------------------------------+
                                                       |
                                                       v
                                        +------------------+
                                        | x (solution)     |
                                        +------------------+
```

### Contribution Block Flow Between Subtrees

```
Leaf Subtree A          Leaf Subtree B           Root Subtree
+-------------+         +-------------+         +-------------------+
| factor nodes|         | factor nodes|         |                   |
|  bottom-up  |         |  bottom-up  |         | assemble contribs |
|             |         |             |         | from children     |
+------+------+         +------+------+         |                   |
       |                       |                | factor remaining  |
       v                       v                | nodes bottom-up   |
  contrib_type            contrib_type          |                   |
  { val, rlist,           { val, rlist,    +--->|                   |
    ndelay,                 ndelay,        |    +-------------------+
    delay_perm,             delay_perm,    |
    delay_val }             delay_val }----+
       |                                   |
       +-----------------------------------+
               child_contrib[] array
```

**Sources:**
`src/ssids/ssids.f90` (analyse/factor/solve dispatch),
`src/ssids/fkeep.F90` (inner_factor_cpu, inner_solve_cpu),
`src/ssids/cpu/NumericSubtree.hxx` (node-level factorization loop),
`src/ssids/cpu/factor.hxx` (factor_node dispatch)

---

## 5. APTP Implementation Details

The A Posteriori Threshold Pivoting (APTP) algorithm is the distinguishing
feature of SSIDS's CPU factorization path. Unlike traditional threshold partial
pivoting (TPP), which selects pivots before performing elimination, APTP
proceeds optimistically and verifies stability after the fact. This enables a
Cholesky-like communication pattern with better data locality and parallelism.

### 5.1 Core File

The entire APTP implementation lives in a single file:

**`src/ssids/cpu/kernels/ldlt_app.cxx`** (2591 lines)

This file is by far the largest in the SSIDS codebase and contains the
`Column<T>` class, backup/restore infrastructure, threshold checking, pivot
application, and the main factorization loop in both task-parallel and serial
variants.

### 5.2 The `Column<T>` Class

```cpp
template<typename T>
class Column {
public:
   bool first_elim;  // True if first column with eliminations
   int nelim;        // Number of eliminated entries in this column
   T *d;             // Pointer to local D storage

   void init_passed(int passed);      // Initialize pass count
   void update_passed(int passed);    // Thread-safe min-reduction on pass count
   bool test_fail(int passed);        // Check if block failed threshold test
   void adjust(int& next_elim);       // Finalize nelim, handle split 2x2 pivots
   void move_back(...);               // Compact eliminated/failed entries
private:
   mutable spral::omp::Lock lock_;    // Lock for thread-safe npass updates
   int npass_ = 0;                    // Reduction variable for eliminated count
};
```

`Column<T>` tracks the elimination progress for one block column of the
factorization. The `npass_` counter records how many blocks in the column have
passed the a posteriori threshold test. All updates to `npass_` are
protected by a lock for thread safety. The `adjust()` method handles the edge
case of a split 2x2 pivot (where the last passed column is the first half of a
Bunch-Kaufman 2x2 pivot) by decrementing `npass_` to avoid an incomplete pivot.

The `ColumnData<T, IntAlloc>` wrapper manages an array of `Column<T>` objects
plus a local permutation vector. Its `calc_nelim()` method counts total
successful eliminations across all block columns by checking that every block in
each column passed.

### 5.3 `check_threshold()`

```cpp
template <enum operation op, typename T>
int check_threshold(int rfrom, int rto, int cfrom, int cto,
                    T u, T* aval, int lda);
```

This function checks whether the entries of a block satisfy the pivot threshold
condition. It scans entries `aval[rfrom:rto, cfrom:cto]` and returns the
index of the first row or column where `|aval[j*lda+i]| > 1/u`. The test is:
an entry in L must satisfy `|L_{ij}| <= 1/u` for the factorization to be
numerically stable with threshold parameter u.

The `op` template parameter determines whether the check is column-oriented
(`OP_N`, returns earliest failed column) or row-oriented (`OP_T`, returns
earliest failed row).

If all entries pass, the function returns `cto` (or `rto`), indicating full
success.

### 5.4 `apply_pivot()`

```cpp
template <enum operation op, typename T>
void apply_pivot(int m, int n, int from, const T *diag, const T *d,
                 const T small, T* aval, int lda);
```

This function performs the operation `L_{21} = A_{21} * L_{11}^{-T} * D_1^{-1}`
(for `op=OP_N`) or the transpose variant (for `op=OP_T`). It first applies the
triangular solve using `host_trsm`, then applies the diagonal scaling by D^{-1}.

The D storage convention is crucial (see Section 6.4). For 1x1 pivots, the
routine multiplies by `d[2*i]`. For 2x2 pivots (detected when `d[2*i+2]` is
not finite, i.e., Inf), it applies the 2x2 block:

```
d11 = d[2*i],  d21 = d[2*i+1],  d22 = d[2*i+3]
```

Zero pivots are handled specially: entries smaller than `small` are zeroed out,
while larger entries are set to infinity to propagate the singularity.

### 5.5 Three Pivot Methods

The `factor()` static method of the inner `LDLT` class dispatches to one of
two internal routines based on the `PivotMethod` enum:

#### APP_AGGRESSIVE (`PivotMethod::app_aggressive`, option value 1)

1. **Optimistic elimination**: `run_elim_unpivoted()` performs the entire
   factorization assuming all pivots are acceptable, using a Cholesky-like
   communication pattern (no pivot search, no row/column interchanges).
2. **Global a posteriori check**: After elimination, `check_threshold()` is
   applied to every block. If all pass, factorization is complete.
3. **On failure**: The factorization state is rolled back to the last
   fully-successful block column using `CopyBackup::restore_part()`, and
   the remaining columns are re-factored using `run_elim_pivoted()` (the
   block APTP method).

This method has the lowest overhead when pivoting is rarely needed (e.g.,
nearly positive-definite matrices) but requires a full copy of the matrix for
backup and a potentially expensive rollback on failure.

#### APP_BLOCK (`PivotMethod::app_block`, option value 2, **default**)

1. **Block-column-wise elimination**: `run_elim_pivoted()` processes one block
   column at a time. Within each block column:
   - Back up the current block column state.
   - Perform a small dense LDL^T factorization of the diagonal block using
     `block_ldlt()` (Bunch-Kaufman style, `INNER_BLOCK_SIZE = 32`).
   - Apply the pivot to below-diagonal blocks via `apply_pivot()`.
   - Check the threshold for each below-diagonal block via
     `check_threshold()`.
   - If a block fails, restore from backup and reduce `nelim` for that column.
2. **Column-level rollback**: A failed pivot only requires recalculation of
   entries within its own block column, not the entire matrix.
3. **Parallelism**: Updates to trailing blocks are issued as OpenMP tasks,
   allowing multiple block columns to be processed concurrently where
   dependencies permit.

This is the default method because it balances parallelism, numerical
robustness, and moderate memory overhead.

#### TPP (`PivotMethod::tpp`, option value 3)

Traditional Threshold Partial Pivoting, implemented in `ldlt_tpp.cxx`. This is
a serial algorithm that selects pivots before elimination. It is used:
- When explicitly requested via `options%pivot_method = 3`.
- As a fallback after APTP methods fail to eliminate all columns, controlled
  by `options%failed_pivot_method`:
  - `FAILED_PIVOT_METHOD_TPP (1)`: Retry failed columns with TPP.
  - `FAILED_PIVOT_METHOD_PASS (2)`: Pass failed columns directly to parent
    node as delays.

### 5.6 Failed Pivot Handling

When APTP fails to eliminate all n columns of a node, the remaining
`n - nelim` columns become "delayed pivots". The handling depends on context:

1. **Root nodes (m == n)**: Always fall back to TPP for remaining columns,
   since delays cannot be passed further up.
2. **Non-root nodes with `failed_pivot_method = tpp`**: Apply `ldlt_tpp_factor`
   to the remaining columns. If TPP also fails, the remaining columns are
   delayed to the parent.
3. **Non-root nodes with `failed_pivot_method = pass`**: Immediately delay the
   remaining columns to the parent node.

Delayed pivots increase the size of the parent node by `ndelay_out` rows and
columns, which are prepended to the parent's frontal matrix during assembly.
The `ndelay_in` / `ndelay_out` fields on `NumericNode` track this.

Statistics tracking:
- `not_first_pass`: Columns not eliminated by the primary APTP method.
- `not_second_pass`: Columns not eliminated by either APTP or TPP fallback.
- `num_delay`: Total delayed columns across all nodes.

**Sources:**
`src/ssids/cpu/kernels/ldlt_app.cxx` (lines 57-162: Column class; lines
303-322: check_threshold; lines 332-413: apply_pivot; lines 2273-2372:
factor dispatch; lines 2506-2534: ldlt_app_factor entry point),
`src/ssids/cpu/kernels/ldlt_tpp.cxx`,
`src/ssids/cpu/kernels/block_ldlt.hxx`,
`src/ssids/cpu/factor.hxx` (lines 36-136: factor_node_indef),
`src/ssids/cpu/cpu_iface.hxx` (PivotMethod enum)

---

## 6. Key Data Structures

### 6.1 `ssids_akeep` (Symbolic Factorization)

Defined in `src/ssids/akeep.f90`. This is the primary output of the analyse
phase and is passed unchanged to all subsequent calls.

| Field | Type | Description |
|-------|------|-------------|
| `n` | `integer` | Order of the matrix |
| `nnodes` | `integer` | Number of nodes (supernodes) in the elimination tree |
| `sptr(nnodes+1)` | `integer, allocatable` | Supernode pointers. Supernode i consists of columns `sptr(i)` through `sptr(i+1)-1` |
| `sparent(nnodes)` | `integer, allocatable` | Parent of node i in the elimination tree. `sparent(i) = nnodes+1` for a root. Parents are always numbered higher than children |
| `rptr(nnodes+1)` | `integer(long), allocatable` | Pointers into `rlist` for each node. Node i has row indices `rlist(rptr(i):rptr(i+1)-1)` |
| `rlist(*)` | `integer, allocatable` | Row indices for all nodes, stored contiguously. Within each node, indices are in elimination order |
| `invp(n)` | `integer(C_INT), allocatable` | Inverse permutation. `invp(i)` gives the position in the original matrix of the i-th variable in elimination order |
| `nptr(nnodes+1)` | `integer(long), allocatable` | Pointers into `nlist` for each node |
| `nlist(2,*)` | `integer(long), allocatable` | Map from original matrix values to factor storage. `nlist(1,j)` is the source index in `val`, `nlist(2,j)` is the destination index in `lcol` |
| `nparts` | `integer` | Number of subtree partitions |
| `part(nparts+1)` | `integer, allocatable` | Node ranges for each partition |
| `subtree(nparts)` | `symbolic_subtree_ptr, allocatable` | Symbolic subtree objects (polymorphic, CPU or GPU) |
| `contrib_ptr(nparts+1)` | `integer, allocatable` | Pointers into `contrib_idx` for child contribution lookup |
| `contrib_idx(nparts)` | `integer, allocatable` | Maps each subtree to its contribution target |
| `check` | `logical` | Whether matrix data was checked during analyse |
| `ptr(:)` | `integer(long), allocatable` | Cleaned column pointers (only if `check=.true.`) |
| `row(:)` | `integer, allocatable` | Cleaned row indices (only if `check=.true.`) |
| `map(:)` | `integer(long), allocatable` | Map from original to cleaned matrix (only if `check=.true.`) |
| `scaling(n)` | `real(wp), allocatable` | Scaling from matching-based ordering (only if `ordering=2`) |
| `topology(:)` | `numa_region, allocatable` | Machine topology for NUMA-aware factorization |

### 6.2 `ssids_fkeep` (Numeric Factorization)

Defined in `src/ssids/fkeep.F90`.

| Field | Type | Description |
|-------|------|-------------|
| `subtree(:)` | `numeric_subtree_ptr, allocatable` | Array of numeric subtree objects, one per partition. Each contains the factored nodes |
| `scaling(:)` | `real(wp), allocatable` | Scaling vector in elimination order. `scaling(i)` is the scale factor for the i-th eliminated variable |
| `pos_def` | `logical` | True if positive-definite factorization was performed |
| `inform` | `ssids_inform` | Copy of inform at end of factorization |

Each numeric subtree internally holds an array of `NumericNode` objects:

### 6.3 `NumericNode` (C++ Per-Node Numeric Data)

Defined in `src/ssids/cpu/NumericNode.hxx`.

| Field | Type | Description |
|-------|------|-------------|
| `symb` | `SymbolicNode const&` | Reference to corresponding symbolic node |
| `first_child` | `NumericNode*` | First child in linked list |
| `next_child` | `NumericNode*` | Next sibling in parent's child list |
| `ndelay_in` | `int` | Number of delayed pivots received from children |
| `ndelay_out` | `int` | Number of delayed pivots to pass to parent |
| `nelim` | `int` | Number of successfully eliminated columns |
| `lcol` | `T*` | Factor values. Layout: first `n` columns are L (column-major, `m x n`), followed by D storage (2n entries). Leading dimension is `align_lda<T>(m + ndelay_in)` |
| `perm` | `int*` | Permutation. `perm[i]` = original column index eliminated at local position i |
| `contrib` | `T*` | Contribution block. Dense lower-triangular `(m-n) x (m-n)` matrix |

Where `m = symb.nrow + ndelay_in` (total rows including delays) and
`n = symb.ncol + ndelay_in` (total columns including delays).

### 6.4 `node_type` (Fortran Per-Node Data)

Defined in `src/ssids/datatypes.f90`. Used by the Fortran layer.

| Field | Type | Description |
|-------|------|-------------|
| `nelim` | `integer` | Number of eliminated columns |
| `ndelay` | `integer` | Number of delays |
| `rdptr` | `integer(long)` | Entry into rebuilt rlist_direct |
| `ncpdb` | `integer` | Number of contributions to parent's diagonal block |
| `lcol(:)` | `real(wp), pointer` | Factor values (same layout as C++ `NumericNode.lcol`) |
| `perm(:)` | `integer, pointer` | Permutation of columns |
| `gpu_lcol` | `C_PTR` | GPU-side factor pointer |

### 6.5 `contrib_type` (Contribution Block Between Subtrees)

Defined in `src/ssids/contrib.f90`.

| Field | Type | Description |
|-------|------|-------------|
| `ready` | `logical` | Flag indicating contribution is ready to be consumed (used for synchronization) |
| `n` | `integer` | Size of the contribution block |
| `val(n,n)` | `real(C_DOUBLE), pointer` | Dense lower-triangular contribution matrix |
| `ldval` | `integer(C_INT)` | Leading dimension of `val` |
| `rlist(:)` | `integer(C_INT), pointer` | Row indices for the contribution (maps rows to original matrix variables) |
| `ndelay` | `integer` | Number of delayed pivots included |
| `delay_perm(:)` | `integer(C_INT), pointer` | Permutation for delayed entries |
| `delay_val(:)` | `real(C_DOUBLE), pointer` | Values for delayed entries |
| `lddelay` | `integer` | Leading dimension of delay_val |
| `owner` | `integer` | Cleanup owner: 0=CPU, 1=GPU |

### 6.6 Diagonal D Storage Convention

The block diagonal D (or rather D^{-1}) is stored as a flat array of `2*n`
entries adjacent to the L factor in each node's `lcol`. The storage convention:

```
D stored as d[0..2n-1]:

1x1 pivot at position i:
  d[2*i]   = D_{ii}^{-1}     (the inverse pivot value)
  d[2*i+1] = 0.0              (signals: not part of a 2x2 pivot)

2x2 pivot at positions i and i+1:
  d[2*i]   = element (1,1) of the 2x2 inverse block
  d[2*i+1] = element (2,1) of the 2x2 inverse block
  d[2*i+2] = Inf              (signals: this is the second half of a 2x2)
  d[2*i+3] = element (2,2) of the 2x2 inverse block

Zero pivot (singular):
  d[2*i]   = 0.0
  d[2*i+1] = 0.0
```

The key sentinel value is `d[2*i+2] = Inf`: when `std::isfinite(d[2*i+2])`
returns false, the code knows that positions i and i+1 form a 2x2 pivot.
This convention is used throughout `apply_pivot()`, `ldlt_app_solve_diag()`,
and `ldlt_app_solve_bwd()`.

**Sources:**
`src/ssids/akeep.f90`,
`src/ssids/fkeep.F90`,
`src/ssids/datatypes.f90` (node_type, lines 116-134),
`src/ssids/contrib.f90`,
`src/ssids/cpu/NumericNode.hxx`,
`src/ssids/cpu/kernels/ldlt_app.cxx` (D convention comment at lines 325-330)

---

## 7. External Dependencies and Configuration

### 7.1 BLAS/LAPACK Routines Used

SSIDS calls the following BLAS and LAPACK routines, wrapped through C-callable
interfaces in `src/ssids/cpu/cpu_iface.f90`:

| Routine | Purpose in SSIDS |
|---------|-----------------|
| `DGEMM` | Schur complement updates: `C -= L_{21} * D * L_{21}^T`. The dominant cost in factorization |
| `DSYRK` | Symmetric rank-k update for contribution blocks: `C -= L * L^T` |
| `DPOTRF` | Cholesky factorization of diagonal blocks (positive-definite path) |
| `DSYTRF` | Symmetric indefinite factorization of small diagonal blocks (used in block_ldlt fallback) |
| `DTRSM` | Triangular solve for off-diagonal blocks: `L_{21} = A_{21} * L_{11}^{-T}` |
| `DTRSV` | Triangular solve for single vectors (used in solve phase) |
| `DGEMV` | Matrix-vector product (used in solve phase for rectangular updates) |

These are called through thin C++ wrappers (`host_gemm<T>`, `host_trsm<T>`,
etc.) defined in `src/ssids/cpu/kernels/wrappers.cxx`.

### 7.2 Optional Dependencies

| Dependency | Purpose | Required? |
|-----------|---------|-----------|
| **METIS** | Fill-reducing ordering (`options%ordering = 1`, the default). Called through `spral_metis_wrapper`. Without METIS, user must supply their own ordering. | Strongly recommended |
| **hwloc** | Hardware topology detection for NUMA-aware subtree assignment. If not available, falls back to `OMP_NUM_THREADS`. | Optional |
| **CUDA / cuBLAS** | GPU factorization path. Enabled via `options%use_gpu = .true.` (default). Only activated if GPU hardware is detected and subtree has sufficient work (`min_gpu_work` threshold). | Optional |
| **OpenMP** | Thread parallelism for both subtree-level and node-level factorization. SSIDS uses OpenMP tasks within each subtree. Requires `OMP_CANCELLATION=true`. | Required for parallelism |

### 7.3 CPU vs GPU Path Separation

The assembly tree is partitioned into subtrees, each assigned to either a CPU
NUMA region or a GPU. The assignment is determined during the analyse phase
based on:

- `options%min_gpu_work`: Minimum flops threshold for GPU assignment (default:
  5 x 10^9).
- `options%gpu_perf_coeff`: Relative GPU-to-CPU performance ratio.
- `options%max_load_inbalance`: Maximum permissible load imbalance between
  NUMA regions.

Factorization proceeds in two stages:
1. **Leaf subtrees** are factored in parallel, each on its assigned resource
   (CPU NUMA region or GPU).
2. **Root subtree** (remaining nodes above all leaf subtrees) is factored
   cooperatively by all CPU threads.

CPU and GPU paths share the abstract interface defined in `subtree.f90` but
have completely separate implementations:
- CPU: `src/ssids/cpu/` (C++ templates)
- GPU: `src/ssids/gpu/` (CUDA kernels, not covered in this review)

Contribution blocks (`contrib_type`) provide the data exchange mechanism
between subtrees, regardless of whether they ran on CPU or GPU.

### 7.4 Key `ssids_options` Fields

Defined in `src/ssids/datatypes.f90`:

| Field | Default | Description |
|-------|---------|-------------|
| `pivot_method` | 2 (`APP_BLOCK`) | CPU pivot strategy: 1=APP_AGGRESSIVE, 2=APP_BLOCK, 3=TPP |
| `failed_pivot_method` | 1 (`TPP`) | What to do with APTP failures: 1=retry with TPP, 2=pass to parent |
| `u` | 0.01 | Relative pivot threshold. Range: [0, 0.5]. Larger values give more stability but more delays |
| `small` | 1e-20 | Absolute threshold below which an entry is treated as zero |
| `cpu_block_size` | 256 | Block size for task generation on larger nodes. Each block column of this width generates one factorization task |
| `nemin` | 32 | Supernode amalgamation threshold. Two neighbors in the elimination tree are merged if both involve fewer than `nemin` eliminations |
| `small_subtree_threshold` | 4 x 10^6 | Maximum flops in a subtree treated as a single serial task (reduces OpenMP scheduling overhead for small subtrees) |
| `ordering` | 1 | Ordering method: 0=user-supplied, 1=METIS, 2=matching-based |
| `scaling` | 0 | Scaling method: <=0=none/user, 1=Hungarian (MC64), 2=Auction, 3=from analyse, >=4=norm-equilibration (MC77) |
| `action` | .true. | Continue on singular matrix (with warning) if true; abort if false |
| `multiplier` | 1.1 | Factor memory over-allocation ratio to accommodate delayed pivots |
| `ignore_numa` | .true. | Treat entire machine as single NUMA region |
| `use_gpu` | .true. | Use GPU if available |
| `min_gpu_work` | 5 x 10^9 | Minimum flops for GPU subtree assignment |
| `max_load_inbalance` | 1.2 | Maximum load imbalance between NUMA regions |

**Sources:**
`src/ssids/cpu/cpu_iface.f90` (BLAS wrappers, lines 100-173),
`src/ssids/cpu/cpu_iface.hxx` (PivotMethod enum),
`src/ssids/datatypes.f90` (ssids_options, lines 188-284),
`src/ssids/cpu/kernels/wrappers.cxx`,
`src/ssids/anal.F90` (subtree partition and assignment),
`docs/Fortran/ssids.rst` (options documentation)

---

## Source Citations

Every section of this document was produced by reading the following SPRAL
source files (all paths relative to `references/spral/`):

**Section 1 (Executive Summary):**
- `docs/Fortran/ssids.rst` -- User documentation
- `examples/Fortran/ssids.f90` -- Usage example
- `src/ssids/ssids.f90` -- Top-level API module

**Section 2 (Module-by-Module Overview):**
- `src/ssids/ssids.f90` (1479 lines)
- `src/ssids/datatypes.f90` (374 lines)
- `src/ssids/akeep.f90` (120 lines)
- `src/ssids/fkeep.F90` (463 lines)
- `src/ssids/inform.f90` (216 lines)
- `src/ssids/contrib.f90` (77 lines)
- `src/ssids/subtree.f90` (128 lines)
- `src/ssids/anal.F90` (1240 lines)

**Section 3 (CPU Factorization Stack):**
- `src/ssids/cpu/subtree.f90` (425 lines)
- `src/ssids/cpu/cpu_iface.f90` (175 lines)
- `src/ssids/cpu/cpu_iface.hxx` (51 lines)
- `src/ssids/cpu/SymbolicNode.hxx` (28 lines)
- `src/ssids/cpu/SymbolicSubtree.hxx` (111 lines)
- `src/ssids/cpu/NumericNode.hxx` (78 lines)
- `src/ssids/cpu/NumericSubtree.hxx` (550 lines)
- `src/ssids/cpu/NumericSubtree.cxx` (262 lines)
- `src/ssids/cpu/SmallLeafSymbolicSubtree.hxx` (161 lines)
- `src/ssids/cpu/SmallLeafNumericSubtree.hxx` (507 lines)
- `src/ssids/cpu/factor.hxx` (187 lines)
- `src/ssids/cpu/ThreadStats.hxx` (64 lines)
- `src/ssids/cpu/kernels/ldlt_app.cxx` (2591 lines)
- `src/ssids/cpu/kernels/ldlt_app.hxx` (26 lines)
- `src/ssids/cpu/kernels/ldlt_tpp.cxx` (317 lines)
- `src/ssids/cpu/kernels/ldlt_nopiv.cxx` (114 lines)
- `src/ssids/cpu/kernels/cholesky.cxx` (215 lines)
- `src/ssids/cpu/kernels/block_ldlt.hxx` (415 lines)
- `src/ssids/cpu/kernels/assemble.hxx` (445 lines)
- `src/ssids/cpu/kernels/calc_ld.hxx` (120 lines)
- `src/ssids/cpu/kernels/common.hxx` (143 lines)
- `src/ssids/cpu/kernels/wrappers.cxx` (92 lines)
- `src/ssids/cpu/kernels/verify.hxx` (183 lines)
- `src/ssids/cpu/kernels/SimdVec.hxx` (200 lines)

**Section 4 (Data Flow Diagram):**
- `src/ssids/ssids.f90` (analyse/factor/solve dispatch)
- `src/ssids/fkeep.F90` (inner_factor_cpu, inner_solve_cpu)
- `src/ssids/cpu/NumericSubtree.hxx` (node factorization loop)
- `src/ssids/cpu/factor.hxx` (factor_node dispatch)

**Section 5 (APTP Implementation Details):**
- `src/ssids/cpu/kernels/ldlt_app.cxx` (core APTP kernel)
- `src/ssids/cpu/kernels/ldlt_tpp.cxx` (TPP fallback)
- `src/ssids/cpu/kernels/block_ldlt.hxx` (inner block factorization)
- `src/ssids/cpu/factor.hxx` (factor_node_indef)
- `src/ssids/cpu/cpu_iface.hxx` (PivotMethod enum)

**Section 6 (Key Data Structures):**
- `src/ssids/akeep.f90`
- `src/ssids/fkeep.F90`
- `src/ssids/datatypes.f90`
- `src/ssids/contrib.f90`
- `src/ssids/cpu/NumericNode.hxx`
- `src/ssids/cpu/kernels/ldlt_app.cxx` (D storage convention)

**Section 7 (External Dependencies and Configuration):**
- `src/ssids/cpu/cpu_iface.f90` (BLAS/LAPACK wrappers)
- `src/ssids/cpu/cpu_iface.hxx` (pivot method enums)
- `src/ssids/datatypes.f90` (ssids_options)
- `src/ssids/cpu/kernels/wrappers.cxx` (C++ BLAS wrappers)
- `src/ssids/anal.F90` (subtree assignment)
- `docs/Fortran/ssids.rst` (options documentation)
