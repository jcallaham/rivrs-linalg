# Sparse Module Restructuring Design

**Date**: 2026-02-25
**Status**: Design document (not yet implemented)
**Context**: The `sparse` crate currently contains a single solver (SSIDS-style multifrontal LDL^T with APTP pivoting). This document analyzes what's general-purpose vs algorithm-specific, proposes a restructuring to make room for additional solvers, and identifies high-value gaps in the Rust sparse LA ecosystem.

---

## 1. Current Architecture: General vs SSIDS-Specific

### Fully general infrastructure (reusable as-is)

| Module | Description |
|--------|-------------|
| `error.rs` | Error types (`SparseError`) |
| `validate.rs` | `backward_error()`, `sparse_backward_error()`, `validate_permutation()` |
| `profiling/` | Chrome Trace export, RSS tracking, section timing |
| `debug/` | Sparsity visualization, elimination tree display |
| `benchmarking/` | Criterion harness, baseline collection |
| `testing/` | `SolverTest` trait, `NumericalValidator`, test case loading |
| `lib.rs` | `SolverPhase` enum, module declarations |
| `io/registry.rs` | Matrix metadata registry |
| `aptp/ordering.rs` | METIS nested dissection |
| `aptp/matching.rs` | MC64 weighted bipartite matching + scaling |
| `aptp/perm.rs` | Permutation utilities |

### Generic multifrontal machinery (reusable with interface changes)

| Module | What's reusable |
|--------|-----------------|
| `aptp/numeric.rs` | `FactorizationWorkspace`, `AssemblyMaps`, extend-add scatter, column-oriented assembly, workspace reuse, small-leaf fast path, tree-level parallelism |
| `aptp/solve.rs` | Per-supernode gather/scatter/solve traversal pattern |
| `aptp/symbolic.rs` | `permute_symbolic_upper_triangle()`, faer `SymbolicCholesky` wrapper |
| `aptp/solver.rs` | Analyze/factor/solve API pattern, ordering dispatch |
| Amalgamation logic | Supernode merge predicates |

### Inherently APTP / symmetric-indefinite-specific

| Module | Why it's specific |
|--------|-------------------|
| `aptp/pivot.rs` | `PivotType` (1x1/2x2/Delayed), `Block2x2` — symmetric Bunch-Kaufman pivots |
| `aptp/diagonal.rs` | `MixedDiagonal` — D factor with mixed 1x1/2x2 block sizes |
| `aptp/factor.rs` | Dense APTP kernel, TPP fallback, threshold checking |
| `aptp/inertia.rs` | Eigenvalue sign counts (symmetric concept) |
| `io/mtx.rs` | Hardcoded symmetric format |
| `io/reference.rs` | `DBlock` reference factorization format |
| `validate.rs` (partial) | `reconstruction_error()` assumes LDL^T structure |

---

## 2. Module Naming: `aptp` vs Alternatives

`aptp` is an implementation detail — it tells you *how* the pivoting works, not *what* the solver does. From a user's perspective:

| Name | Pros | Cons |
|------|------|------|
| `aptp` | Precise, distinguishes from other LDL^T variants | Implementation detail; users care about "sparse symmetric indefinite solve" |
| `ldlt` | Describes the factorization | Too narrow if we add LL^T; faer already has LDL^T |
| `ssids` | Recognizable to SPRAL users | Too specific to one reference implementation |
| `symmetric` | Describes the matrix class | Correct scope — houses LDL^T, LL^T, Bunch-Kaufman |
| `direct` | Describes solver class | Too broad — encompasses unsymmetric too |
| `multifrontal` | Describes algorithm framework | Implementation detail, excludes simplicial |

**Recommendation**: `symmetric` as the public-facing module, with `aptp` becoming an internal submodule. The user-facing type `SparseLDLT` is already well-named regardless of module path.

---

## 3. Proposed Module Structure

```
sparse/src/
├── lib.rs
├── error.rs                    # General
├── validate.rs                 # General
│
├── io/                         # General infrastructure
│   ├── mtx.rs                  # Generalize: add format param (symmetric/general)
│   └── registry.rs             # General
│
├── ordering/                   # Extracted from aptp/
│   ├── metis.rs                # metis_ordering()
│   ├── matching.rs             # mc64_matching(), match_order_metis()
│   └── perm.rs                 # perm_from_forward()
│
├── multifrontal/               # Generic multifrontal machinery
│   ├── workspace.rs            # FactorizationWorkspace (parameterized by factor type)
│   ├── assembly.rs             # AssemblyMaps, extend_add, scatter_original
│   ├── symbolic.rs             # Symbolic analysis wrapping faer
│   ├── amalgamation.rs         # Supernode merge
│   └── tree.rs                 # Tree traversal, parallelism, small-leaf fast path
│
├── symmetric/                  # Symmetric solvers (current aptp/)
│   ├── solver.rs               # SparseLDLT (public API)
│   ├── factor.rs               # Multifrontal numeric (AptpNumeric)
│   ├── solve.rs                # Per-supernode solve
│   ├── aptp/                   # APTP-specific internals
│   │   ├── kernel.rs           # aptp_factor_in_place, TPP
│   │   ├── pivot.rs            # PivotType, Block2x2
│   │   └── diagonal.rs         # MixedDiagonal
│   └── inertia.rs              # Inertia
│
├── profiling/                  # General
├── debug/                      # General
├── testing/                    # General
└── benchmarking/               # General
```

---

## 4. Phased Migration Strategy

### Short term (before Phase 9.2 release)

Don't restructure yet. The public API type is `SparseLDLT`, which is already the right user-facing name. The `aptp` module path is an internal detail that doesn't leak into the user experience much.

### Medium term (when adding a second solver)

Extract `ordering/` as a sibling module to `aptp/` (or `symmetric/`). MC64 and METIS ordering are the most reusable, highest-value assets — they deserve to be first-class modules, not buried inside a solver-specific directory.

### Longer term (if building unsymmetric multifrontal)

Extract multifrontal assembly machinery (workspace, assembly maps, tree traversal, extend-add) into a generic `multifrontal/` module parameterized by factor type. This is where the real architectural payoff would be, but it's premature to design that interface without a concrete second consumer.

---

## 5. Gaps: High-Value Functionality Beyond faer

### What faer 0.22 provides in sparse

- **Matrix types**: CSC and CSR (owned, ref, mut, symbolic)
- **Symmetric factorization**: Simplicial + supernodal LL^T and LDL^T (with regularization, not APTP), supernodal intranode LBL^T (Bunch-Kaufman)
- **Unsymmetric factorization**: LU (partial pivoting), QR (no pivoting)
- **Ordering**: AMD, COLAMD
- **Triangular solve**: Sparse forward/backward substitution
- **Sparse matmul**: Symbolic + numeric sparse-sparse multiplication
- **Parallelism**: `Par` enum for tree-level and intra-node BLAS

### What faer does NOT have

**High-value gaps (strong candidates for rivrs-sparse):**

| Gap | Value | Notes |
|-----|-------|-------|
| MC64 matching + scaling | Very high | We already have this; benefits any solver, not just symmetric |
| METIS ordering | Very high | We have this via metis-sys; faer only has AMD/COLAMD |
| Robust indefinite symmetric solve | High | faer's LDL^T uses regularization, not APTP; ours handles genuinely indefinite systems |
| Iterative solvers (CG, GMRES, BiCGStab) | Very high | Workhorse for large systems; nothing in faer or the Rust ecosystem |
| Incomplete factorizations (ILU, ICC) | High | Preconditioners for iterative solvers; complements direct solvers |
| Sparse eigenvalue solver (Lanczos, LOBPCG) | High | Extremely common need; nothing in faer |
| Block-structured solvers (saddle point, KKT) | Medium | Common in optimization, PDE, coupled physics |

**Medium-value gaps:**

| Gap | Value | Notes |
|-----|-------|-------|
| Graph partitioning interface | Medium | Our METIS wrapper could be generalized |
| Sparse backward error | Medium | Our `sparse_backward_error()` is a useful diagnostic |
| Solver profiling infrastructure | Medium | Chrome Trace export for sparse solvers |

---

## 6. Unsymmetric Multifrontal (MUMPS-style) Reuse Assessment

If we were to add an unsymmetric multifrontal LU solver:

### Directly reusable (~40% of codebase)

- **Ordering**: `metis_ordering()`, `mc64_matching()` — MC64 is even *more* important for unsymmetric matrices (diagonal dominance improvement)
- **IO**: `registry.rs`, `mtx.rs` (with format generalization to accept `general` matrices)
- **Infrastructure**: profiling, debug, testing, benchmarking — all algorithm-agnostic
- **Validation**: `sparse_backward_error()`, `validate_permutation()`

### Reusable with refactoring (~25%)

- **`AssemblyMaps` and extend-add scatter**: Core multifrontal assembly pattern is identical, but frontal matrices are rectangular (m x n, m >= n) instead of square
- **`FactorizationWorkspace`**: Same workspace-reuse pattern, different buffer shapes
- **Tree traversal and parallelism**: Same postorder traversal, same rayon pattern
- **Small-leaf fast path**: Same subtree classification algorithm
- **Supernode amalgamation**: Similar merge predicates

### Must be rewritten (~35%)

- **Dense kernel**: LU with partial pivoting instead of APTP
- **Factor storage**: L + U factors instead of L + D (no `MixedDiagonal`)
- **Pivot tracking**: Row permutations instead of 1x1/2x2 block classification
- **Solve**: Forward L, then backward U (no diagonal middle step)
- **Symbolic**: Unsymmetric elimination tree (column etree), different fill estimates

---

## 7. Strategic Recommendations

### Highest-impact next solvers

1. **Iterative solvers + preconditioners** (CG + ICC for symmetric, GMRES + ILU for general) — These complement direct solvers perfectly, fill a large gap in the Rust ecosystem, and our direct solver could serve as a preconditioner for very large systems (domain decomposition, Schur complement).

2. **Sparse eigenvalue solver** (Lanczos for symmetric, Arnoldi for general) — Extremely common need, requires sparse matrix-vector product and a direct solver for shift-invert mode (which we already have).

3. **Unsymmetric multifrontal LU** — Significant reuse of ordering, assembly, and infrastructure code. Would validate the `multifrontal/` extraction and prove the architecture generalizes.

### Principle: extract on the second consumer

Don't prematurely generalize. Extract shared modules (`ordering/`, `multifrontal/`) when there's a concrete second consumer that exercises the interface. The current flat `aptp/` structure is fine for a single solver — premature extraction risks designing the wrong abstractions.

---

## 8. Shared-Memory vs Distributed-Memory Parallelism

### Summary

Distributed memory mostly **layers on top** of node-level solvers rather than requiring fundamental redesign — but the layering point depends heavily on the solver class.

### Direct solvers: distributed is a different beast

For multifrontal direct solvers, there are two levels where distribution could happen:

**Tree-level distribution** — assign subtrees of the elimination tree to different nodes. Independent subtrees factor in parallel, with communication only at merge points (extend-add across node boundaries). This is how MUMPS and PaStiX work at the coarse level. In principle our tree-level rayon parallelism could be replaced with MPI at this boundary, but...

**Root-front distribution** — the real problem. In 3D PDEs, most of the FLOPs concentrate in the large supernodes near the root of the elimination tree. For a cube mesh with n^3 unknowns, the root separator is O(n^2) and its dense frontal matrix is O(n^4) to factor. You can't parallelize this by giving subtrees to different nodes — you need 2D block-cyclic distribution of the dense frontal matrix *itself* across nodes (ScaLAPACK-style). This is what makes MUMPS and SuperLU_DIST so complex. It's fundamentally different from shared-memory BLAS.

A distributed multifrontal solver is a multi-year, multi-team effort (MUMPS has had 25+ years of development). It's not something to build into rivrs-sparse — it's a different project entirely.

### Iterative solvers: distributed layers on naturally

An iterative solver (CG, GMRES, BiCGStab) has a simple inner loop:

```
repeat:
    w = A * v          # sparse matrix-vector product
    z = M^{-1} * w     # preconditioner application
    update solution + convergence check
```

**SpMV distributes trivially** — partition rows across nodes, each node owns its chunk of the matrix, halo exchange for off-node column entries. The communication pattern is fixed (determined by sparsity structure).

**Preconditioner is the interesting part:**
- Block Jacobi / additive Schwarz: each node applies a *local* direct solve on its subdomain. Our `SparseLDLT` is exactly the right tool for this — runs on each node independently.
- Overlap Schwarz: same idea, slightly more communication.
- ILU/ICC: inherently sequential in their global form, but block variants work per-subdomain.

Key insight: **domain decomposition methods are designed so that the node-level work is a self-contained sparse direct solve**. Our shared-memory solver becomes a *component* inside a distributed framework, not something that needs to be distributed itself.

### What actually needs distributed memory?

**Clearly needs distributed:**
- Large-scale 3D PDEs (CFD, structural mechanics, geophysics) — O(n^2) fill makes direct solvers hit memory walls fast, problem sizes routinely exceed single-node RAM
- Coupled multiphysics (fluid-structure, thermo-mechanical) — very large saddle-point systems
- Time-dependent 3D simulations — must solve repeatedly, problem must fit in aggregate memory

**Usually fine on a single node:**
- 2D PDEs — O(n log n) fill, direct solvers handle millions of unknowns on a workstation
- Structural eigenproblems — shift-invert Lanczos with a direct solve
- Optimization (interior point) — KKT systems are structured; Schur complement approaches keep subproblems manageable
- Circuit simulation — large but very sparse, low fill

**Distributed for throughput, not memory:**
- Parametric studies / uncertainty quantification — many independent solves, embarrassingly parallel at the problem level
- Graph analytics — distributed SpMV but rarely direct solvers

The pattern: **distributed memory is compelling primarily when a single problem doesn't fit in a single node's RAM, which is mainly 3D discretizations.** Most other use cases either fit on a node or parallelize at a higher level.

### Recommendation

1. **Keep core solvers shared-memory** — the node-level solver is the building block.

2. **Design for composability** — the key API property is that our solver can be used as a subdomain solver inside a distributed framework. Requirements:
   - Accept matrices as input (don't own the mesh/discretization)
   - Support repeated factorization with the same sparsity (analyze once, factor many times)
   - Keep memory footprint predictable (so a distributed scheduler can plan)
   - We already have all of this.

3. **If we add iterative solvers**, the distributed extension point is clear: define a `LinearOperator` trait (matvec interface), implement distributed SpMV as one implementation, and let preconditioners use our local direct solver.

4. **If someone needs distributed direct solves**, the pragmatic answer is domain decomposition: partition the mesh, use our solver per subdomain, iterate via Schwarz/FETI/BDDC. This gets 90% of the scalability of MUMPS with 10% of the implementation complexity.

5. **Rust + MPI is still immature** — `rsmpi` exists but is a thin wrapper. The Rust ecosystem doesn't have a PETSc/Trilinos equivalent. Building that infrastructure is a separate project from building good solvers.

---

## 9. Release Readiness Assessment

### Performance vs SPRAL

- **41 out of 65 SuiteSparse matrices beat SPRAL** (63% win rate)
- **Median ratio: 0.98x** (rivrs is 2% faster on median)
- Flagship results: c-71 at 1.49x SPRAL, c-big at 1.37x SPRAL (improved from 2.16x and 2.30x respectively via Phase 9.1g column-oriented extend-add)
- Best performers: thread (0.41x), cvxqp3 (0.54x), ncvxqp5 (0.61x)
- Worst performers: a few matrices at 1.6-1.9x SPRAL (algorithm-level differences in pivoting strategy, not implementation quality)
- **All 65 matrices pass correctness** (backward error < 5e-11, typically ~1e-17)
- Remaining optimization targets (zeroing 8.8%, per-node storage ~9%) show diminishing returns

### API Design

**Strengths:**
- Clean three-phase API: analyze (reusable) -> factor (repeatable) -> solve
- Both allocating `solve()` and in-place `solve_in_place()` variants
- Transparent workspace via faer's `MemStack`
- `SparseLDLT` is a plain struct with no global state — multiple solvers can coexist
- Permutation, inertia, and factorization stats all queryable
- Parallelism controllable via `Par::Seq` / `Par::rayon(n)`
- faer types at the boundary (not custom wrappers)

**Gaps:**
- Single-column RHS only — no batched multi-RHS solve via dense matrix
- No condition number estimation
- No iterative refinement
- No serde support for solver options
- `OrderingStrategy::UserSupplied` exists but could be more discoverable

### Memory Efficiency

- Aggressive workspace reuse: frontal matrix, contribution buffer, g2l mapping all pre-allocated and reused across supernodes
- Swap-based contribution buffer protocol eliminates copy overhead
- Thread-local workspaces via `Cell` for parallel path
- Small-leaf fast path uses dedicated <=512KB workspace
- Peak RSS dominated by max_front^2 (inherent to multifrontal method, unavoidable)

### Code Quality

- Zero TODOs/FIXMEs/HACKs in source
- No `.unwrap()` in non-test code
- Rich error types with `thiserror`, proper `Result` propagation
- 401 unit tests, 65 SuiteSparse integration, 15 hand-constructed with analytical factorizations
- Comprehensive doc comments with academic citations on all public types

### Key Gaps for Release

**Not blockers but worth considering:**

| Feature | Impact | Effort | Notes |
|---------|--------|--------|-------|
| Multi-RHS solve | Medium | Low | `MatMut<'_, f64>` variant of solve; loop internally |
| Iterative refinement | Medium | Medium | One Newton step on residual; common in production solvers |
| Condition estimation | Low | Medium | Hager's 1-norm estimator; useful diagnostic |
| Sparse residual utility | Low | Low | `A*x - b` computation; exists as test utility but not public |
| Serde for options | Low | Low | Enables config-file-driven workflows |

### Verdict

The solver is production-ready. Performance is at parity with SPRAL, the API is clean and composable, memory is well-managed, and code quality is high. The missing features (multi-RHS, iterative refinement, condition estimation) are enhancements for post-1.0, not prerequisites. The most impactful pre-release work would be documentation polish and the multi-RHS solve variant.
