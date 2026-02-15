# Feature Specification: Multifrontal Numeric Factorization

**Feature Branch**: `015-multifrontal-factorization`
**Created**: 2026-02-15
**Status**: Draft
**Input**: Implement Phase 6 of ssids-plan.md — Multifrontal Numeric Factorization

## Clarifications

### Session 2026-02-15

- Q: How should Phase 6 handle matrices where `AptpSymbolic` produces simplicial (non-supernodal) decomposition? → A: Handle simplicial as trivial supernodes (each column = 1-column front), so all matrices pass through the same multifrontal code path.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Factor a Sparse Symmetric Indefinite Matrix (Priority: P1)

A solver developer provides a sparse symmetric matrix and its symbolic analysis (from Phase 3's `AptpSymbolic`) and receives a complete numeric factorization that decomposes the matrix into `P^T A P = L D L^T` via the multifrontal method. The factorization traverses the assembly tree in postorder, assembling and factoring frontal matrices at each supernode using Phase 5's dense APTP kernel, producing per-supernode factors (L11, D11, L21) and propagating contribution blocks to parent supernodes.

**Why this priority**: This is the core deliverable of Phase 6 — without end-to-end multifrontal factorization, none of the downstream phases (triangular solve, user-facing API) can proceed. This represents the transition from dense-only APTP to sparse-capable APTP.

**Independent Test**: Can be tested by factoring hand-constructed matrices with known sparsity patterns and verifying reconstruction error `||P^T A P - L D L^T|| / ||A|| < 10^-12`.

**Acceptance Scenarios**:

1. **Given** a sparse symmetric positive definite matrix and its `AptpSymbolic` analysis, **When** multifrontal factorization is invoked, **Then** the resulting factors satisfy reconstruction error < 10^-12 and all pivots are 1x1.
2. **Given** a sparse symmetric indefinite matrix (e.g., from the hand-constructed test set), **When** multifrontal factorization is invoked with default options, **Then** the factorization completes with correct 1x1 and 2x2 pivot decisions, and reconstruction error < 10^-12.
3. **Given** a sparse matrix and `AptpSymbolic` with supernodal structure, **When** factorization completes, **Then** the result contains per-supernode factors indexed by supernode ID, each with L11, D11 (as `MixedDiagonal`), and L21 blocks.
4. **Given** a sparse matrix, **When** factorization completes, **Then** factorization statistics report total 1x1 pivots, total 2x2 pivots, total delayed columns, and maximum front size.

---

### User Story 2 - Correctly Handle Delayed Pivots (Priority: P1)

When the dense APTP kernel cannot eliminate all fully-summed columns of a frontal matrix (due to numerical instability), the delayed columns are propagated to the parent supernode as additional fully-summed columns. This mechanism ensures that the factorization completes even for numerically challenging indefinite matrices — delayed columns get a second chance at elimination in the context of a larger frontal matrix.

**Why this priority**: Delayed pivot propagation is essential for correctness on indefinite systems. Without it, the factorization would fail on any matrix where a pivot is too small relative to its column at the current front but acceptable at an ancestor front. This is the mechanism that distinguishes multifrontal APTP from a simple block diagonal factorization.

**Independent Test**: Can be tested by constructing matrices where specific pivots fail at a child supernode but succeed at the parent, and verifying that the factorization completes with the correct number of delayed pivots reported.

**Acceptance Scenarios**:

1. **Given** a matrix where a specific column fails the APTP threshold at a child supernode, **When** factorization runs, **Then** the column is included in the parent's frontal matrix as an additional fully-summed column and is eventually eliminated.
2. **Given** a matrix with cascading delays (delayed from child to parent to grandparent), **When** factorization completes, **Then** all columns are eventually eliminated (zero delayed columns at the root) and reconstruction error < 10^-12.
3. **Given** a matrix where delays occur, **When** the factorization result is inspected, **Then** the `FactorizationStats` accurately report `total_delayed` reflecting the total number of delay events.

---

### User Story 3 - Assembly of Frontal Matrices (Priority: P1)

For each supernode in assembly-tree postorder, the multifrontal method assembles a dense frontal matrix by: (a) scattering original matrix entries from the sparse input into the appropriate positions, and (b) performing extend-add operations to merge contribution blocks from child supernodes. The assembled frontal matrix must have all fully-summed columns correctly populated before the APTP kernel is called.

**Why this priority**: Correct assembly is a prerequisite for correct factorization — if entries are placed in wrong positions or child contributions are lost, the factorization will silently produce wrong results. Assembly correctness is independently testable.

**Independent Test**: Can be tested by assembling frontal matrices for small matrices with known elimination trees and verifying that the assembled dense block matches the expected dense submatrix of the permuted sparse matrix.

**Acceptance Scenarios**:

1. **Given** a sparse matrix with known sparsity pattern, **When** a frontal matrix is assembled for a leaf supernode (no children), **Then** the assembled F11 and F21 blocks contain exactly the original matrix entries at the correct positions, with all other entries zero.
2. **Given** a supernode with child supernodes whose contribution blocks are known, **When** the frontal matrix is assembled, **Then** the extend-add operation correctly merges child contributions into the frontal matrix using row-index mapping.
3. **Given** a supernode with both original entries and child contributions, **When** assembly completes, **Then** the fully-summed portion of the frontal matrix equals the sum of scattered original entries and extend-added child contributions.

---

### User Story 4 - Schur Complement and Contribution Block Computation (Priority: P2)

After factoring the fully-summed columns (F11) of a frontal matrix, the solver computes the Schur complement to produce the contribution block (updated F22). This involves solving for L21 using triangular solve and updating F22 via a symmetric rank-k update. The contribution block is then available for extend-add into the parent supernode's frontal matrix.

**Why this priority**: This is a mathematically necessary step but is a direct application of dense linear algebra (TRSM + SYRK/GEMM) using faer. It carries less implementation risk than assembly or delayed pivots but must be correct for the factorization to work end-to-end.

**Independent Test**: Can be tested by factoring a single frontal matrix, extracting its contribution block, and verifying that re-assembling the parent front and factoring produces the same result as dense factorization of the full permuted matrix.

**Acceptance Scenarios**:

1. **Given** a factored frontal matrix with known L11, D11, **When** L21 is computed as `F21 * L11^{-T} * D11^{-1}`, **Then** the resulting L21 matches the expected subdiagonal factor block.
2. **Given** L21 and D11, **When** the Schur complement `F22 - L21 * D11 * L21^T` is computed, **Then** the contribution block matches the expected dense Schur complement to machine precision.
3. **Given** a multi-level assembly tree, **When** contribution blocks flow from leaves to root, **Then** the final reconstruction `P^T A P = L D L^T` is correct across all supernodes.

---

### Edge Cases

- What happens when a supernode has zero fully-summed columns after accounting for child structure? (Degenerate supernode — should be handled gracefully, possibly skipped.)
- What happens when all columns of a frontal matrix are delayed? (All columns propagate to parent; the front produces no factors for this supernode.)
- What happens when a matrix has only a single supernode (the entire matrix is one front)? (Reduces to dense APTP factorization — results must match Phase 5's output.)
- What happens when the root supernode receives delayed columns? (They must still be eliminated or reported as structurally/numerically singular.)
- How does the factorization handle 1x1 supernodes (simplicial columns)? (These should be treated as trivial fronts.)
- What happens when a 2x2 pivot straddles a supernode boundary? (The APTP kernel operates within a single frontal matrix, so this cannot happen — 2x2 pivots are always within one front's fully-summed columns.)

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST traverse the assembly tree in postorder (leaves before parents), processing each supernode's frontal matrix exactly once.
- **FR-002**: System MUST assemble each frontal matrix by scattering original sparse matrix entries into the correct dense positions using the global-to-local row index mapping.
- **FR-003**: System MUST perform the extend-add operation to merge contribution blocks from all child supernodes into the parent's frontal matrix, using row-index intersection to place entries correctly.
- **FR-004**: System MUST call Phase 5's `aptp_factor_in_place()` on the fully-summed portion (F11) of each frontal matrix, using the provided `AptpOptions`.
- **FR-005**: System MUST compute L21 by solving `L21 = F21 * L11^{-T} * D11^{-1}` using dense triangular solve and diagonal solve operations.
- **FR-006**: System MUST compute the Schur complement (contribution block) as `F22_updated = F22 - L21 * D11 * L21^T` using dense matrix operations.
- **FR-007**: System MUST propagate delayed columns from a child supernode to its parent's frontal matrix as additional fully-summed columns, expanding the parent's F11 block accordingly.
- **FR-008**: System MUST store per-supernode factors (L11, D11 as `MixedDiagonal`, L21) indexed by supernode ID in the output factorization result.
- **FR-009**: System MUST accept `AptpSymbolic` (Phase 3) and `SparseColMat<usize, f64>` as inputs, using the symbolic analysis to determine supernode structure, assembly tree, and fill-reducing permutation.
- **FR-010**: System MUST produce factorization statistics including total 1x1 pivots, total 2x2 pivots, total delayed pivot events, and maximum front size encountered.
- **FR-011**: System MUST return an error if the matrix is structurally singular (detected during assembly) or if the root supernode has unresolvable delayed columns.
- **FR-012**: System MUST use the fill-reducing permutation from `AptpSymbolic` to map between original matrix indices and permuted indices during assembly.
- **FR-013**: System MUST correctly handle the case where the `AptpSymbolic` analysis produced a supernodal decomposition (the expected case for Phase 6).
- **FR-014**: System MUST handle simplicial `AptpSymbolic` decompositions by treating each column as a trivial 1-column front (trivial supernode), ensuring all matrices pass through the same multifrontal code path regardless of whether faer chose simplicial or supernodal analysis.

### Key Entities

- **FrontalMatrix**: A dense matrix representing a supernode's contribution to the factorization. Partitioned into F11 (fully-summed, to be factored), F21 (subdiagonal block), and F22 (contribution block). Carries global row indices for local-to-global mapping.
- **FrontFactors**: The stored result of factoring one supernode — L11 (unit lower triangular), D11 (`MixedDiagonal` with 1x1/2x2 pivots), L21 (subdiagonal block), and the list of delayed column indices propagated to the parent.
- **AptpNumeric**: The complete numeric factorization result containing all per-supernode `FrontFactors`, aggregate statistics, and references back to the symbolic analysis.
- **ContributionBlock**: The Schur complement (updated F22) from a factored supernode, along with its global row indices. Consumed by the parent supernode during extend-add assembly.
- **FactorizationStats**: Aggregate statistics from the factorization — pivot counts, delay counts, front sizes, total floating-point operations performed.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Reconstruction error `||P^T A P - L D L^T|| / ||A||` is below 10^-12 for all 15 hand-constructed test matrices.
- **SC-002**: Reconstruction error is below 10^-12 for all 10 CI-subset SuiteSparse matrices.
- **SC-003**: On matrices small enough to also factor densely (n < 500), the multifrontal factorization produces the same reconstruction error (within 10x tolerance) as Phase 5's dense APTP factorization applied to the full dense matrix.
- **SC-004**: Delayed pivot events are correctly propagated — no matrix in the test suite produces a reconstruction error above tolerance due to mishandled delays.
- **SC-005**: Factorization of the full SuiteSparse "easy indefinite" test matrices (when run via `--ignored` tests) achieves reconstruction error below 10^-12.
- **SC-006**: Factorization statistics accurately reflect the total pivot counts and maximum front size (verified against hand-computed values for small matrices).
- **SC-007**: Inertia computed from the combined `MixedDiagonal` factors across all supernodes matches the expected inertia for all hand-constructed matrices with known eigenvalue sign distributions.

## Algorithm References

The multifrontal factorization algorithm is based on the following academic sources:

1. **Duff & Reid (1983)** — "The multifrontal solution of indefinite sparse symmetric linear equations." Foundational paper defining the multifrontal method, extend-add operator, contribution block propagation, and delayed pivot mechanism for indefinite systems.
   - Reference: `/workspace/rivrs-linalg/references/ssids/duff1984.md`

2. **Liu (1992)** — "The Multifrontal Method for Sparse Matrix Solution: Theory and Practice" (SIAM Review). Comprehensive survey formalizing frontal matrices, update matrices, elimination/assembly trees, supernodes, and the extend-add operator with mathematical definitions.
   - Reference: `/workspace/rivrs-linalg/references/ssids/liu1992.md`

3. **Hogg, Ovtchinnikov, Scott (2016)** — "A Sparse Symmetric Indefinite Direct Solver for GPU Architectures" (ACM TOMS). Describes SSIDS implementation including assembly (`assemble_pre`/`factor`/`assemble_post` pattern), threshold partial pivoting within multifrontal framework, and memory management.
   - Reference: `/workspace/rivrs-linalg/references/ssids/hogg2016.md`

4. **Hogg, Duff, Lopez (2020)** — "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting." Defines the APTP algorithm (Algorithm 3.1) and its integration with multifrontal factorization (Section 5), including the `assemble_pre`/`factor`/`assemble_post` three-phase pattern and failed-pivot handling.
   - Reference: `/workspace/rivrs-linalg/references/ssids/duff2020.md`

5. **Davis (2006/2016)** — Survey of direct methods for sparse linear systems. Additional context on multifrontal implementation techniques.
   - Reference: `/workspace/rivrs-linalg/references/ssids/davis2016.md`

## Assumptions

- **Unified code path for both decomposition types**: Phase 6 handles both supernodal and simplicial `AptpSymbolic` decompositions through a single multifrontal code path. When faer produces a simplicial decomposition (common for small matrices), each column is treated as a trivial 1-column front. The supernodal path remains the primary performance target.
- **Sequential execution**: Phase 6 implements single-threaded factorization. Parallelism (tree-level task parallelism, intra-node BLAS parallelism) is deferred to Phase 9.
- **Stack-based contribution storage**: Contribution blocks from child supernodes are consumed by their parent and can be deallocated after extend-add. A postorder traversal naturally supports stack-based lifetime management.
- **Dense operations via faer**: All dense linear algebra (matmul, TRSM, SYRK) uses faer's existing dense kernel APIs. No custom BLAS implementations.
- **Permutation handling**: The fill-reducing permutation from `AptpSymbolic::perm()` is applied to the sparse matrix before assembly. Within each front, `aptp_factor_in_place()` may produce an additional local permutation for pivoting, which must be tracked per-supernode.
- **AptpOptions shared across all fronts**: The same threshold and fallback strategy apply to every supernode's APTP factorization. Per-supernode tuning is not supported in this phase.

## Scope Boundaries

**In scope:**
- Multifrontal factorization loop (assembly tree traversal)
- Frontal matrix assembly (scatter + extend-add)
- Per-supernode APTP factorization via Phase 5's kernel
- Schur complement computation (L21 solve + contribution block update)
- Delayed pivot propagation to parent supernodes
- Per-supernode factor storage (L11, D11, L21)
- Factorization statistics collection
- Reconstruction-based correctness validation

**Out of scope:**
- Triangular solve and user-facing API (`SparseLDLT`) (Phase 7)
- Parallel factorization (tree-level or intra-node) (Phase 8)
- Two-level APTP blocking within large fronts (Phase 8.1)
- Workspace/arena memory optimization (Phase 9)
- Performance optimization for the simplicial-as-trivial-fronts path (functional correctness required, not performance parity with faer's native simplicial LDL^T)

## Dependencies

- **Phase 2** (APTP Data Structures): `PivotType`, `Block2x2`, `MixedDiagonal`, `Inertia`, `perm_from_forward`
- **Phase 3** (AptpSymbolic): `AptpSymbolic` providing supernodal structure, assembly tree, permutation, and symbolic pattern
- **Phase 5** (Dense APTP Kernel): `aptp_factor_in_place()`, `AptpFactorResult`, `AptpOptions`, `AptpFallback`
- **faer 0.22**: Dense matrix operations (`Mat`, `MatMut`, `MatRef`), triangular solve, matrix multiply, permutation types
