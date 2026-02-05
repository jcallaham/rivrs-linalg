# Research: Complete Phase 0.1 Reference Library

**Branch**: `001-reference-library` | **Date**: 2026-02-05

## Research Questions

This feature is documentation-only, so "research" means auditing existing materials and understanding their content well enough to synthesize deliverables. Three research threads were conducted in parallel.

## R1: SPRAL SSIDS Architecture Review

**Decision**: Document SPRAL's three-phase solver architecture with emphasis on APTP kernel location and data flow.

**Rationale**: SPRAL is the primary BSD-3 reference implementation. Understanding its module structure enables informed design of the Rust port.

**Key Findings**:

### Module Overview

| Module | File | Responsibility |
|--------|------|---------------|
| Main API | `src/ssids/ssids.f90` (~53K lines) | Public entry points: analyse, factor, solve, free, enquire, alter |
| Analysis keep | `src/ssids/akeep.f90` (~121 lines) | Symbolic factorization output: sptr, sparent, rptr, rlist |
| Factor keep | `src/ssids/fkeep.F90` (~15K lines) | Numeric factorization output: subtree array, scaling |
| Data types | `src/ssids/datatypes.f90` (~17K lines) | ssids_options, node_type, thread_stats |
| Inform | `src/ssids/inform.f90` (~217 lines) | Status reporting: flags, statistics, error codes |
| Contribution | `src/ssids/contrib.f90` (~78 lines) | Multifrontal contribution blocks between subtrees |
| Subtree base | `src/ssids/subtree.f90` (~129 lines) | Abstract interface for CPU/GPU subtrees |
| Core analysis | `src/core_analyse.f90` (~850+ lines) | Elimination tree, postorder, column counts, supernodes |
| Symbolic analysis | `src/ssids/anal.F90` (~44K lines) | Pattern expansion, analyse_phase dispatch |

### CPU Factorization Stack

| Layer | File | Role |
|-------|------|------|
| CPU subtree wrapper | `src/ssids/cpu/subtree.f90` (~6K lines) | Fortran-C++ interop |
| Numeric subtree | `src/ssids/cpu/NumericSubtree.hxx` (~21K lines) | C++ template; OMP task-based factorization |
| Numeric node | `src/ssids/cpu/NumericNode.hxx` (~2.5K lines) | Per-supernode factorization state |
| Symbolic subtree | `src/ssids/cpu/SymbolicSubtree.hxx` (~4.5K lines) | Symbolic info for CPU path |
| **APTP kernel** | `src/ssids/cpu/kernels/ldlt_app.cxx` (~2.6K lines) | **Core APTP implementation** |
| Block LDLT | `src/ssids/cpu/kernels/block_ldlt.hxx` (~14K lines) | Block-level dense operations |
| TPP fallback | `src/ssids/cpu/kernels/ldlt_tpp.cxx` (~11K lines) | Traditional threshold partial pivoting |
| No-pivot LDLT | `src/ssids/cpu/kernels/ldlt_nopiv.cxx` (~3.7K lines) | Positive definite fast path |
| Cholesky | `src/ssids/cpu/kernels/cholesky.cxx` (~8.3K lines) | LL^T for positive definite |
| Assembly | `src/ssids/cpu/kernels/assemble.hxx` (~16.6K lines) | Contribution assembly into parent front |
| LD computation | `src/ssids/cpu/kernels/calc_ld.hxx` (~4.3K lines) | Pivot value computation and testing |

### Three-Phase Data Flow

```
ANALYSE: Matrix A (CSC lower) + ordering options
  → expand to full symmetric → compute elimination tree
  → postorder → column counts → detect supernodes → partition into subtrees
  → OUTPUT: akeep {sptr, sparent, rptr, rlist, subtree[], invp}

FACTOR: akeep + matrix values + posdef flag
  → parallel over subtrees (OMP taskgroup)
  → per subtree: bottom-up tree traversal
    → per node: assemble children contributions → APTP or Cholesky factorization
    → generate contribution block for parent
  → OUTPUT: fkeep {subtree[] with factor values, perm, nelim per node}

SOLVE: fkeep + akeep + RHS x
  → forward substitution (leaves → root)
  → diagonal solve (D^{-1})
  → backward substitution (root → leaves)
  → OUTPUT: x (modified in-place)
```

### APTP Pivot Methods

- `PIVOT_METHOD_APP_AGGRESSIVE (1)`: Single-pass; on failure, restart entire node
- `PIVOT_METHOD_APP_BLOCK (2)` [default]: Block-column-level rollback
- `PIVOT_METHOD_TPP (3)`: Traditional serial pivoting (fallback only)

### Diagonal D Storage Convention

```
1x1 pivot: d[2k] = value,  d[2k+1] = 0.0
2x2 pivot: d[2k] = a,      d[2k+1] = b
            d[2k+2] = INF,  d[2k+3] = c
Zero pivot: d[2k] = 0.0,   d[2k+1] = 0.0
```

### External Dependencies

- BLAS: DGEMM, DSYRK, DTRSM, DTRSV
- LAPACK: DPOTRF, DSYTRF (fallback)
- Optional: METIS (ordering), hwloc (topology), CUDA (GPU path)

**Alternatives considered**: None — SPRAL is the only BSD-3 SSIDS implementation.

---

## R2: faer Sparse Infrastructure Assessment

**Decision**: Classify each faer component by reuse potential (direct use / adapt / reference only).

**Rationale**: faer provides ~70% of needed sparse infrastructure. Identifying exact reuse points prevents duplicating work.

**Key Findings**:

### faer Version

- Version: 0.24.0
- Commit: `8dfccee` ("bump version to 0.24.0")
- License: MIT

### Component Classification

| Component | File Path | Size | Classification | APTP Adaptation |
|-----------|-----------|------|----------------|-----------------|
| CSC storage | `sparse/csc/` | Multiple files | **Direct use** | None needed |
| CSR storage | `sparse/csr/` | Multiple files | **Direct use** | For input conversion |
| AMD ordering | `sparse/linalg/amd.rs` | 24 KB | **Direct use** | Call directly for initial ordering |
| COLAMD ordering | `sparse/linalg/colamd.rs` | 19 KB | **Direct use** | Available if needed |
| Elimination tree | In `cholesky.rs` | Part of 164 KB | **Direct use** | Extract etree computation |
| Triangular solve | `sparse/linalg/triangular_solve.rs` | — | **Direct use** | Forward/backward substitution |
| Permutation utils | `perm/` | — | **Direct use** | Forward/inverse permutations |
| Workspace mgmt | `MemStack`/`StackReq` | — | **Direct use** | Follow same allocation pattern |
| Sparse symbolic | `sparse/linalg/cholesky.rs` | Part of 164 KB | **Adapt** | Reuse ereach; adapt for indefinite |
| Sparse numeric | `sparse/linalg/cholesky.rs` | Part of 164 KB | **Adapt** | Replace Bunch-Kaufman with APTP |
| Dense LBLT | Dense module | — | **Reference** | Understand 2x2 block handling |
| Sparse LU | `sparse/linalg/lu.rs` | 72 KB | **Reference** | Etree patterns, supernode handling |

### Critical: cholesky.rs Analysis

faer's `cholesky.rs` (~164 KB, ~2700 lines) implements:
- **Both LLT and LDLT** factorizations (controlled by `is_llt` flag)
- **Simplicial** (column-by-column) and **supernodal** (blocked) modes
- **Bunch-Kaufman style LBLT** with 1x1 and 2x2 blocks
- **Regularization** via `eps`/`delta`/`signs` parameters

**Key difference from APTP**: faer's Bunch-Kaufman decides pivots **before** elimination. APTP decides pivots **after** (a posteriori) — this is the fundamental architectural difference. The column-by-column update pattern is reusable; the pivot decision logic must be replaced entirely.

### Workspace Pattern

faer uses `MemStack` (stack-based allocation from pre-allocated buffer) with `StackReq` for computing required sizes. Functions compute their scratch requirements via `*_scratch()` companion functions. This pattern should be adopted for SSIDS.

**Alternatives considered**: Building sparse infrastructure from scratch — rejected because faer already provides well-tested, SIMD-optimized implementations.

---

## R3: APTP Algorithm Extraction

**Decision**: Extract pseudocode for all key algorithms from the academic papers, organized by function.

**Rationale**: A unified pseudocode reference prevents developers from cross-referencing multiple papers during implementation.

**Key Findings**:

### Algorithms Extracted

1. **Algorithm 3.1** (Duff et al. 2020, §3): Main APTP factorization loop — 2D block partitioning with Factor/Apply/Adjust/Update kernels
2. **Algorithm 4.1** (Duff et al. 2020, §4): Complete pivoting for dense blocks — 1x1/2x2 selection with stability bounds
3. **Pivot acceptance test** (Duff et al. 2020, §3): `|l_ij| > u^{-1}` checked after applying pivot
4. **Column delay mechanism** (Duff et al. 2020, §3, §7): Fail-in-place with backup/restore
5. **GPU APTP variant** (Hogg et al. 2016, §3.2): Slightly different 2x2 test with Θ=100 bias
6. **Symbolic column counts** (Gilbert et al. 1994, §2): O(mα(m,n)) prediction of nonzero structure
7. **Multifrontal structure** (Liu 1992, Hogg et al. 2016, §2.3): Frontal matrix assembly and contribution blocks

### Pivot Decision Tree

```
For each column j in block:
  1. Attempt 1x1 pivot: use diagonal a_jj
  2. Apply to all blocks in column j
  3. Check: |l_ij| ≤ u^{-1} for all i?
     YES → Accept 1x1 pivot
     NO  → 4. Try 2x2 pivot with best partner
            Test: |Δ| ≥ (1/2)|a_mt|^2?
            YES → Accept 2x2 pivot
            NO  → 5. Delay column to parent node
                     Restore from backup
```

### Threshold Parameter Effects (Table 5, Duff et al. 2020)

| u | Stability | Fill-in | Delays | Speed |
|---|-----------|---------|--------|-------|
| 0.1 | Highest | Most | Most | Slowest |
| 0.01 (default) | Good | Moderate | Few | Fast |
| 0.001 | Minimum | Least | Fewest | Variable |

**Alternatives considered**: Using PARDISO's SBK approach (perturbed pivoting) — rejected because APTP provides proven stability bounds equivalent to TPP.
