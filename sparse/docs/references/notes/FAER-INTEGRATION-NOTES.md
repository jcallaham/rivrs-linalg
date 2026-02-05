# faer Integration Notes for SSIDS

This document maps faer's sparse infrastructure to the rivrs-sparse SSIDS
(Sparse Symmetric Indefinite Direct Solver) project. It identifies reusable
components, classifies them by integration strategy, and explains how faer's
existing Cholesky/Bunch-Kaufman implementation relates to the APTP (A Posteriori
Threshold Pivoting) algorithm we plan to implement.

---

## 1. Version and Overview

| Field        | Value                                                          |
|--------------|----------------------------------------------------------------|
| Library      | faer                                                           |
| Version      | 0.24.0                                                         |
| Commit       | `8dfccee` ("bump version to 0.24.0")                           |
| License      | MIT                                                            |
| Repository   | <https://codeberg.org/sarah-quinones/faer>                     |
| JOSS Paper   | El Kazdadi, S. (2023). "faer: A linear algebra library for the Rust programming language" |
| Rust Edition | 2021, MSRV 1.84.0                                             |

### Sparse Module Structure

faer's sparse infrastructure is organized under `faer/src/sparse/` with the
following layout:

```
faer/src/
  sparse/
    mod.rs            -- Top-level: CSC/CSR re-exports, FaerError, Pair/Triplet types
    csc/mod.rs        -- CSC storage (60 KB, ~2570 lines)
    csr/mod.rs        -- CSR storage (52 KB, ~2290 lines)
    ops.rs            -- Sparse matrix add/multiply operations
    utils.rs          -- Permutation of self-adjoint matrices, ghost indexing utilities
    solvers.rs        -- High-level solver wrappers
    linalg/
      mod.rs          -- SupernodalThreshold, SymbolicSupernodalParams, error types
      amd.rs          -- Approximate Minimum Degree ordering (24 KB, ~1050 lines)
      colamd.rs       -- Column Approximate Minimum Degree ordering (20 KB, ~760 lines)
      cholesky.rs     -- LLT, LDLT, LBLT factorization (160 KB, ~5585 lines)
      triangular_solve.rs -- Sparse triangular solve (28 KB, ~1003 lines)
      lu.rs           -- Sparse LU factorization (72 KB, ~2710 lines)
      qr.rs           -- Sparse QR factorization (92 KB)
      matmul.rs       -- Sparse matrix multiplication (12 KB)
  perm/
    mod.rs            -- Permutation types and row/column permutation (12 KB)
    permown.rs        -- Owned permutation (4 KB)
    permref.rs        -- Permutation reference (8 KB)
```

The library uses `dyn-stack` (v0.13.2) for stack-based workspace allocation
throughout all sparse algorithms, avoiding heap allocation in inner loops.

**Source**: `references/faer-rs/faer/src/sparse/mod.rs`,
`references/faer-rs/faer/src/sparse/linalg/mod.rs`,
`references/faer-rs/faer/Cargo.toml`

---

## 2. Component Reuse Classification Table

The following table classifies each faer component by its reuse strategy for the
SSIDS project. Classifications:

- **Direct use**: Use as a dependency via the faer crate; no modification needed.
- **Adapt**: The faer implementation provides the right structure but requires
  modification for APTP-specific behavior.
- **Reference**: Study the code for patterns and algorithms; reimplement
  independently.

| # | Component | faer Path | Size | Classification | APTP Notes |
|---|-----------|-----------|------|----------------|------------|
| 1 | CSC storage | `sparse/csc/mod.rs` | 60 KB | Direct use | Primary sparse format for SSIDS input and factor storage. |
| 2 | CSR storage | `sparse/csr/mod.rs` | 52 KB | Direct use | Needed for transpose operations and row-oriented access during solve. |
| 3 | AMD ordering | `sparse/linalg/amd.rs` | 24 KB | Direct use | Fill-reducing ordering for symmetric matrices; directly applicable. |
| 4 | COLAMD ordering | `sparse/linalg/colamd.rs` | 20 KB | Direct use | Column ordering; useful for unsymmetric or QR-based fallback paths. |
| 5 | Elimination tree | `sparse/linalg/cholesky.rs` (simplicial module) | Part of 160 KB | Direct use | `prefactorize_symbolic_cholesky` computes etree and column counts; reusable for SSIDS symbolic phase. |
| 6 | Triangular solve | `sparse/linalg/triangular_solve.rs` | 28 KB | Direct use | Lower/upper triangular solves with unit/generic diagonal; used directly for SSIDS solve phase. |
| 7 | Permutation utilities | `perm/mod.rs`, `perm/permown.rs`, `perm/permref.rs` | 24 KB total | Direct use | `PermRef`, `Perm`, `permute_rows`, `permute_cols`, `swap_rows_idx`; handles all reordering in SSIDS. |
| 8 | Workspace management | `dyn-stack` crate (v0.13.2) | External | Direct use | `MemStack`, `StackReq`, `MemBuffer` pattern; avoids per-call heap allocation. |
| 9 | Symbolic Cholesky structure | `sparse/linalg/cholesky.rs` (simplicial + supernodal) | Part of 160 KB | Adapt | `SymbolicSimplicialCholesky` and `SymbolicSupernodalCholesky` store factor sparsity; APTP needs additional fields for delayed columns and 2x2 pivot tracking. |
| 10 | Numeric Cholesky factorization | `sparse/linalg/cholesky.rs` (numeric functions) | Part of 160 KB | Adapt | Column-by-column update loop and `ereach` traversal are the right skeleton; must replace pivot logic with a posteriori stability check. |
| 11 | Dense Bunch-Kaufman (LBLT) | `linalg/cholesky/lblt/` (dense) + supernodal LBLT | Used in supernodal | Reference | Studies how faer handles 2x2 pivots (subdiagonal, permutation arrays); APTP differs fundamentally in pivot decision timing. |
| 12 | Sparse LU | `sparse/linalg/lu.rs` | 72 KB | Reference | Threshold pivoting and column elimination in an unsymmetric context; useful reference for pivot search strategies. |

**Source**: All paths relative to `references/faer-rs/faer/src/`.

---

## 3. Per-Component Deep Dive

This section provides details on the eight "direct use" components, including key
types, public API entry points, and integration notes for SSIDS.

### 3.1 CSC Storage (`sparse/csc/mod.rs`)

**Key Types**:
- `SymbolicSparseColMat<I>` -- Owned symbolic structure (col_ptr, row_idx, dimensions)
- `SymbolicSparseColMatRef<'a, I>` -- Borrowed symbolic view
- `SparseColMat<I, T>` -- Owned numeric CSC matrix
- `SparseColMatRef<'a, I, T>` -- Borrowed numeric CSC view
- `SparseColMatMut<'a, I, T>` -- Mutably borrowed numeric CSC view

**Public API**:
- `SparseColMat::try_new_from_triplets(nrows, ncols, &triplets)` -- Construction
  from triplet format with duplicate summation.
- `SymbolicSparseColMat::try_new_from_indices(nrows, ncols, &indices)` --
  Symbolic construction with separate `Argsort` for deferred value filling.
- `SparseColMatRef::row_idx_of_col(j)` / `val_of_col(j)` -- Column iteration.
- `col_ptr()`, `row_idx()`, `col_range(j)` -- Raw access to CSC arrays.
- `compute_nnz()` -- Total nonzero count.
- `to_row_major()` -- Conversion to CSR.

**Integration Notes**: CSC is the native format for SSIDS. The input matrix will
be provided as `SparseColMat<usize, f64>` (or `<u32, f64>` for smaller problems).
faer's CSC handles sorted and unsorted indices and supports the "uncompressed"
mode (`col_nnz` for per-column counts distinct from col_ptr gaps), which is useful
during assembly of the factor.

**Source**: `references/faer-rs/faer/src/sparse/csc/mod.rs`

### 3.2 AMD Ordering (`sparse/linalg/amd.rs`)

**Key Types**:
- `Control` -- Tuning parameters: `dense` threshold (default 10.0), `aggressive`
  absorption (default true).
- `FlopCount` -- Estimated flop counts: `n_div`, `n_mult_subs_ldl`, `n_mult_subs_lu`.

**Public API**:
- `order(perm, perm_inv, A, control, stack)` -- Compute AMD ordering of a sorted
  symmetric matrix. Returns `Result<FlopCount, FaerError>`.
- `order_maybe_unsorted(perm, perm_inv, A, control, stack)` -- Same but handles
  unsorted CSC input.
- `order_scratch::<I>(n, nnz_upper)` -- Workspace requirement for sorted input.
- `order_maybe_unsorted_scratch::<I>(n, nnz_upper)` -- Workspace requirement for
  unsorted input.

**Integration Notes**: SSIDS symbolic analysis begins with a fill-reducing
ordering. AMD is the default choice (SPRAL also uses AMD). faer's implementation
follows the classic Amestoy-Davis-Duff algorithm with aggressive absorption.
For larger problems, an external Metis binding could complement AMD, but AMD is
sufficient for the initial implementation. The `FlopCount` output is also useful
for the simplicial-vs-supernodal decision heuristic.

**Source**: `references/faer-rs/faer/src/sparse/linalg/amd.rs`

### 3.3 COLAMD Ordering (`sparse/linalg/colamd.rs`)

**Key Types**:
- `Control` -- Parameters: `dense_row`, `dense_col` thresholds (default 0.5),
  `aggressive` absorption.

**Public API**:
- `order(perm, perm_inv, A, control, stack)` -- Compute column ordering for QR.
- `order_scratch::<I>(nrows, ncols, A_nnz)` -- Workspace requirement.

**Integration Notes**: COLAMD is primarily for unsymmetric QR orderings. It is
not directly needed for symmetric SSIDS but could be useful if we later support
unsymmetric LU or for ordering auxiliary structures. Listed as "direct use"
because it is already available at no additional cost.

**Source**: `references/faer-rs/faer/src/sparse/linalg/colamd.rs`

### 3.4 Elimination Tree (`sparse/linalg/cholesky.rs`, simplicial module)

**Key Types**:
- `EliminationTreeRef<'a, I>` -- Wrapper around `&[I::Signed]` where each entry
  is either a non-negative parent index or `-1` (no parent / root).

**Public API**:
- `prefactorize_symbolic_cholesky(etree, col_counts, A, stack)` -- Computes the
  elimination tree and per-column nonzero counts from the upper triangular part
  of A. Returns `EliminationTreeRef`.
- `prefactorize_symbolic_cholesky_scratch::<I>(n, nnz)` -- Workspace: `O(n)`.

**Algorithm**: Processes columns left-to-right. For each column j, walks up the
etree from each row index i < j, linking unvisited nodes to j as parent and
incrementing column counts. This is the standard "row subtree" algorithm from
Davis (2006).

**Integration Notes**: The elimination tree is the foundation of the SSIDS
symbolic phase. faer's implementation is directly reusable -- call
`prefactorize_symbolic_cholesky` to get the etree, then use it for symbolic
factorization and `ereach` during numeric factorization.

The internal `ereach` function (elimination tree reach) computes the set of
columns in L that contribute to column k, traversed in topological order. This is
critical for the column-by-column numeric factorization loop and is reusable for
APTP.

**Source**: `references/faer-rs/faer/src/sparse/linalg/cholesky.rs`, lines
511--655

### 3.5 Triangular Solve (`sparse/linalg/triangular_solve.rs`)

**Public API** (all operate in-place on `rhs: MatMut<'_, T>`):
- `solve_lower_triangular_in_place(tril, conj, rhs, par)` -- Solve L*x = b.
- `solve_unit_lower_triangular_in_place(tril, conj, rhs, par)` -- Solve L*x = b
  (unit diagonal).
- `solve_lower_triangular_transpose_in_place(tril, conj, rhs, par)` -- Solve
  L^T*x = b.
- `solve_upper_triangular_in_place(triu, conj, rhs, par)` -- Solve U*x = b.
- Corresponding `_unit_` and `_transpose_` variants for upper triangular.

**Implementation Details**: Uses a 4-column unrolling strategy for multiple
right-hand sides, processing columns of the RHS in groups of 4, 3, 2, or 1.
Diagonal elements are assumed to be the first stored element per column (lower)
or last stored element (upper).

There is also an internal `ldlt_scale_solve_unit_lower_triangular_transpose_in_place_impl`
that combines the D^{-1} scaling with the L^T solve for LDLT factorizations.

**Integration Notes**: The SSIDS solve phase computes x = P^T L^{-T} D^{-1}
L^{-1} P b (for simplicial LDLT) or a variant with 2x2 block diagonal for APTP.
faer's triangular solves handle the L^{-1} and L^{-T} steps directly. For APTP,
we will need to add a block-diagonal solve step that handles mixed 1x1 and 2x2
pivots, but the triangular solve itself is used as-is.

**Source**: `references/faer-rs/faer/src/sparse/linalg/triangular_solve.rs`

### 3.6 Permutation Utilities (`perm/`)

**Key Types**:
- `Perm<I>` -- Owned permutation (forward + inverse arrays).
- `PermRef<'a, I>` -- Borrowed permutation view.

**Public API**:
- `PermRef::new_checked(fwd, inv, n)` -- Create from validated arrays.
- `permute_rows(dst, src, perm)` -- dst[i, :] = src[perm[i], :].
- `permute_cols(dst, src, perm)` -- dst[:, j] = src[:, perm[j]].
- `permute_rows_in_place(matrix, perm, stack)` -- In-place row permutation.
- `permute_cols_in_place(matrix, perm, stack)` -- In-place column permutation.
- `swap_rows_idx(mat, a, b)` / `swap_cols_idx(mat, a, b)` -- Single row/column swap.

**Integration Notes**: SSIDS uses permutations at three stages:
1. Fill-reducing ordering (AMD) before factorization.
2. Numeric pivoting permutation (APTP may permute columns within supernodes).
3. Solution phase: apply P and P^T to right-hand sides.

faer's permutation utilities handle all these cases. The
`permute_self_adjoint_to_unsorted` function in `sparse/utils.rs` is particularly
useful for applying the fill-reducing ordering to a symmetric matrix.

**Source**: `references/faer-rs/faer/src/perm/mod.rs`,
`references/faer-rs/faer/src/sparse/utils.rs`

### 3.7 Workspace Management (`dyn-stack` crate)

**Key Types** (from `dyn-stack` v0.13.2):
- `StackReq` -- A static description of memory requirements (size + alignment).
  Composable via `StackReq::all_of` (sum) and `StackReq::or` (max).
- `MemBuffer` -- A heap-allocated buffer satisfying a `StackReq`.
- `MemStack` -- A stack allocator backed by a `MemBuffer`. Sub-allocations are
  made via `stack.make_raw::<T>(n)` which returns `(&mut [T], &mut MemStack)`.
- `stack.collect(iter)` -- Allocate and fill from iterator.

**Pattern**:
```rust
// 1. Compute workspace requirement
let req = StackReq::all_of(&[
    StackReq::new::<usize>(n),
    StackReq::new::<f64>(nnz),
]);

// 2. Allocate buffer
let mut mem = MemBuffer::try_new(req)?;

// 3. Create stack and sub-allocate
let stack = MemStack::new(&mut mem);
let (workspace_a, stack) = unsafe { stack.make_raw::<usize>(n) };
let (workspace_b, _)     = unsafe { stack.make_raw::<f64>(nnz) };
```

**Integration Notes**: SSIDS should adopt this pattern throughout. Benefits:
- Single allocation per phase (symbolic analysis, numeric factorization, solve).
- No allocator calls in inner loops.
- Workspace requirements are statically composable, enabling the caller to
  pre-allocate a single buffer for the entire pipeline.
- The `StackReq::or` combinator is useful when phases share workspace.

The rivrs-sparse API should expose `*_scratch` functions following faer's
convention, so users can compute workspace requirements upfront.

**Source**: `references/faer-rs/faer/Cargo.toml` (dependency),
`references/faer-rs/faer/src/sparse/linalg/cholesky.rs` (usage throughout)

---

## 4. cholesky.rs Analysis

faer's `cholesky.rs` is the single largest file in the sparse module (160 KB,
~5585 lines). It implements three factorization variants in both simplicial
and supernodal modes, plus a unified high-level API.

### 4.1 What cholesky.rs Implements

**Factorization variants** (all for symmetric/Hermitian matrices):
- **LLT** ($A = LL^H$): Standard Cholesky for positive-definite matrices.
- **LDLT** ($A = LDL^H$): Diagonal pivoting for indefinite matrices; D is a
  real diagonal (1x1 pivots only); fails on zero pivots.
- **LBLT** ($A = LBL^H$): Bunch-Kaufman factorization; B is block-diagonal with
  1x1 and 2x2 blocks; handles indefinite matrices without failure.

**Structural variants**:
- **Simplicial**: Processes factor entries one at a time. Efficient when the
  factor L is very sparse (few nonzeros per column). Uses `ereach` for column
  dependencies and a column-by-column update loop.
- **Supernodal**: Groups columns with identical sparsity patterns into
  "supernodes" and uses dense BLAS operations on frontal matrices. Efficient
  when L is moderately dense (many nonzeros per column).

**Automatic selection**: `factorize_symbolic_cholesky` computes the
flop-to-nonzero ratio and selects simplicial or supernodal based on
`SupernodalThreshold` (default: ratio factor 40.0).

**Regularization support**: Both LLT and LDLT accept regularization parameters
(`LltRegularization`, `LdltRegularization`) with dynamic epsilon/delta for
handling near-singular matrices.

### 4.2 Bunch-Kaufman vs APTP: The Pivoting Difference

This is the central design distinction for SSIDS.

**Bunch-Kaufman (faer's LBLT)**:
- Pivot selection is **a priori** (decided before elimination of each column).
- At each step, the algorithm examines the current column and selects either a
  1x1 or 2x2 pivot based on the Bunch-Kaufman criterion (comparing the largest
  off-diagonal element to the diagonal).
- The pivot search requires accessing column entries that may not yet be fully
  updated, necessitating column lookups and potentially disrupting data locality.
- In the supernodal variant, faer calls
  `linalg::cholesky::lblt::factor::cholesky_in_place` on each dense frontal
  matrix, which performs dense Bunch-Kaufman within the supernode. The resulting
  `subdiag` array and `perm_forward`/`perm_inverse` arrays track the 2x2 blocks
  and local permutations.

**APTP (A Posteriori Threshold Pivoting)**:
- Pivot assessment is **a posteriori** (decided after a tentative elimination).
- The algorithm *optimistically* performs elimination assuming a 1x1 pivot.
- After elimination, it checks a stability criterion (typically
  |l_{ij}| <= 1/threshold for all entries in the column).
- If the criterion fails, the column is either:
  (a) paired with another column to form a 2x2 Bunch-Kaufman-style pivot, or
  (b) delayed to be processed later (pushed to a subsequent column/supernode).
- This approach has better data locality because it processes columns in order
  without needing to search for suitable pivots in advance.

**Why the difference matters**:
- APTP enables better vectorization and cache utilization because the "try then
  check" approach processes data sequentially.
- APTP naturally supports delayed pivots (columns that cannot be stably
  eliminated are deferred), which is essential for sparse indefinite systems
  where the fill-reducing ordering may place difficult pivots early.
- SPRAL's SSIDS uses APTP for exactly these reasons; the APTP approach was
  designed for the sparse supernodal context.

### 4.3 Reusable Patterns from cholesky.rs

Despite the pivot-logic difference, several structural patterns from cholesky.rs
transfer directly to an APTP implementation:

**1. Symbolic analysis pipeline** (lines 4627--4772):
The `factorize_symbolic_cholesky` function orchestrates:
ordering -> permute self-adjoint -> compute etree -> compute column counts ->
decide simplicial vs supernodal -> build symbolic structure.
This entire pipeline is reusable for SSIDS; only the symbolic structure types
need extension for delayed columns.

**2. Column-by-column update loop** (simplicial numeric):
The simplicial factorization iterates over columns, uses `ereach` to find
contributing ancestors, accumulates updates, then performs the diagonal operation.
For APTP, the same loop structure applies; only the "diagonal operation" changes
from a simple divide to a stability check with potential delay/pairing.

**3. `ereach` algorithm** (lines 624--655):
Computes the set of columns in L whose sparsity patterns overlap with column k,
in topological order. This is unchanged between Cholesky and APTP -- the
sparsity-driven traversal is identical.

**4. Supernodal assembly pattern** (lines 3470--3719):
The supernodal LBLT factorization shows how to:
- Map global row indices to local supernode indices (`global_to_local`).
- Assemble the original matrix entries into frontal matrices.
- Accumulate update contributions from descendant supernodes using dense
  matrix multiplication (`triangular::matmul` and `matmul::matmul`).
- Apply a dense pivot factorization to the assembled frontal matrix.
- Scale the off-diagonal block by D^{-1}.
This pattern maps directly to APTP supernodal factorization; only the dense
factorization kernel changes (from dense Bunch-Kaufman to dense APTP with delay
support).

**5. Postorder traversal** (lines 3721--3785):
`ghost_postorder` computes a postorder of the etree/supernode tree. This ensures
that all descendants of a node are processed before the node itself, which is
required for both Cholesky and APTP factorization.

### 4.4 What Must Be Built New for APTP

The following components have no direct equivalent in faer and must be
implemented from scratch (or from SPRAL reference):

1. **A posteriori stability check**: After tentatively eliminating a column,
   verify that all entries in the computed L column satisfy
   |l_{ij}| <= u^{-1} (where u is the pivot threshold, typically 0.01). This
   is the core APTP criterion.

2. **Column delay mechanism**: When a column fails the stability check and
   cannot form a 2x2 pivot, it must be "delayed" -- removed from the current
   supernode and appended to the parent supernode's frontal matrix. This
   requires:
   - Tracking delayed columns in the symbolic structure.
   - Expanding the frontal matrix to accommodate delayed columns.
   - Maintaining a mapping from original column indices to their current
     position.

3. **Hybrid 1x1/2x2 pivot tracking**: APTP produces a block diagonal D with
   mixed 1x1 and 2x2 blocks. The solve phase must handle this structure.
   faer's LBLT already stores `subdiag` and `perm_forward`/`perm_inverse` for
   exactly this purpose -- these data structures can be reused, but the logic
   for populating them differs.

4. **Modified symbolic structure**: `SymbolicSimplicialCholesky` and
   `SymbolicSupernodalCholesky` assume that column counts are fixed at symbolic
   analysis time. APTP delays can cause additional fill-in beyond the
   statically predicted pattern. The symbolic structure may need:
   - Conservative over-allocation for delayed columns.
   - Or a two-pass approach: optimistic allocation + reallocation on delay.

5. **Threshold parameter management**: APTP has a pivot threshold u (typically
   0.01--0.1) that controls the stability/fill-in trade-off. This needs to be
   exposed in the solver configuration API.

**Source**: `references/faer-rs/faer/src/sparse/linalg/cholesky.rs`

---

## 5. Integration Strategy and Workspace Patterns

### 5.1 Recommended Approach: faer as Direct Dependency

The recommended strategy is to use faer as a Cargo dependency for all
"direct use" components, rather than vendoring or reimplementing them.

**Rationale**:
- faer is MIT-licensed, compatible with Apache-2.0.
- The sparse infrastructure (CSC, ordering, etree, triangular solve,
  permutations) is mature and well-tested.
- Using faer directly avoids maintenance burden for generic index types,
  SIMD-accelerated dense kernels, and ordering algorithms.
- The `dyn-stack` workspace pattern integrates cleanly with our own code.

**Cargo dependency**:
```toml
[dependencies]
faer = { version = "0.24", default-features = false, features = ["sparse-linalg"] }
```

The `sparse-linalg` feature enables all sparse linear algebra modules. The
`default-features = false` avoids pulling in Rayon (parallel execution) unless
explicitly desired; parallelism can be added later.

### 5.2 MemStack/StackReq Workspace Pattern

All rivrs-sparse algorithms should follow faer's workspace convention:

1. **Every function that needs temporary memory** takes a `&mut MemStack`
   parameter and has a corresponding `*_scratch` function returning `StackReq`.

2. **Scratch functions are composable**:
   - `StackReq::all_of(&[a, b, c])` -- Total memory for sequential phases.
   - `StackReq::or(a, b)` -- Maximum of two mutually exclusive phases.

3. **The caller computes requirements and allocates once**:
   ```rust
   let req = StackReq::all_of(&[
       symbolic_analysis_scratch(n, nnz),
       numeric_factorize_scratch::<f64>(n, nnz),
       solve_scratch::<f64>(n, nrhs),
   ]);
   let mut mem = MemBuffer::try_new(req)?;
   let stack = MemStack::new(&mut mem);
   ```

4. **Inside functions**, sub-allocate with destructuring:
   ```rust
   let (workspace, stack) = unsafe { stack.make_raw::<f64>(n) };
   let (more_workspace, _) = unsafe { stack.make_raw::<usize>(n) };
   ```

This pattern ensures zero heap allocations during the hot path (numeric
factorization), which is critical for performance on large sparse systems.

### 5.3 Suggested Build Order for SSIDS Components

Based on the component analysis and dependency structure, the recommended
implementation order is:

**Phase 1 -- Symbolic Analysis** (mostly faer direct use):
1. Matrix input: accept `SparseColMat<I, f64>` from faer.
2. Fill-reducing ordering: call `amd::order_maybe_unsorted`.
3. Permute matrix: call `sparse::utils::permute_self_adjoint_to_unsorted`.
4. Elimination tree: call `simplicial::prefactorize_symbolic_cholesky`.
5. Symbolic factorization: adapt `factorize_simplicial_symbolic_cholesky` with
   APTP-aware column count estimation (conservative upper bounds for delays).

**Phase 2 -- Simplicial Numeric APTP** (adapt faer patterns):
1. Column-by-column factorization loop following faer's simplicial structure.
2. Use faer's `ereach` algorithm (or reimplement from the same academic source)
   for column dependency traversal.
3. Implement APTP stability check after each tentative column elimination.
4. Implement column delay mechanism (mark column as delayed, advance to next).
5. Implement 2x2 pivot pairing as fallback before delay.

**Phase 3 -- Solve Phase** (mostly faer direct use):
1. Forward permutation: `perm::permute_rows_in_place`.
2. Forward substitution: `triangular_solve::solve_unit_lower_triangular_in_place`.
3. Block-diagonal solve: custom implementation for mixed 1x1/2x2 D.
4. Back substitution: `triangular_solve::solve_unit_lower_triangular_transpose_in_place`.
5. Inverse permutation: `perm::permute_rows_in_place` with inverse perm.

**Phase 4 -- Supernodal Optimization** (adapt faer patterns):
1. Supernode detection from elimination tree (follow faer's amalgamation logic).
2. Dense frontal matrix assembly (follow faer's `global_to_local` pattern).
3. Dense APTP factorization kernel for frontal matrices.
4. Update matrix computation using dense BLAS via faer.
5. Delayed column propagation between supernodes.

---

## 6. Source Citations

All file paths are relative to `references/faer-rs/`.

| Section | Source File(s) |
|---------|---------------|
| Version and overview | `faer/Cargo.toml`, `paper.md` |
| CSC storage | `faer/src/sparse/csc/mod.rs` |
| CSR storage | `faer/src/sparse/csr/mod.rs` |
| AMD ordering | `faer/src/sparse/linalg/amd.rs` |
| COLAMD ordering | `faer/src/sparse/linalg/colamd.rs` |
| Elimination tree | `faer/src/sparse/linalg/cholesky.rs` (lines 511--655) |
| Triangular solve | `faer/src/sparse/linalg/triangular_solve.rs` |
| Permutation utilities | `faer/src/perm/mod.rs`, `faer/src/perm/permown.rs`, `faer/src/perm/permref.rs` |
| Workspace management | `faer/Cargo.toml` (dyn-stack dep), usage in `faer/src/sparse/linalg/cholesky.rs` |
| Symbolic Cholesky | `faer/src/sparse/linalg/cholesky.rs` (lines 673--770, 2392--2952, 4627--4772) |
| Numeric factorization | `faer/src/sparse/linalg/cholesky.rs` (lines 1033--1168, 3135--3719) |
| Dense LBLT (Bunch-Kaufman) | `faer/src/sparse/linalg/cholesky.rs` (supernodal LBLT, lines 3470--3719) |
| Sparse LU | `faer/src/sparse/linalg/lu.rs` |
| Sparse module structure | `faer/src/sparse/mod.rs`, `faer/src/sparse/linalg/mod.rs` |
| Self-adjoint permutation | `faer/src/sparse/utils.rs` |
