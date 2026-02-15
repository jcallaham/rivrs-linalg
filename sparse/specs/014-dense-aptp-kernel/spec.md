# Feature Specification: Dense APTP Factorization Kernel

**Feature Branch**: `014-dense-aptp-kernel`
**Created**: 2026-02-14
**Status**: Draft
**Input**: Phase 5 of SSIDS development plan — Dense A Posteriori Threshold Pivoting kernel for symmetric indefinite factorization

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Factor Dense Symmetric Indefinite Matrices (Priority: P1)

A sparse solver performing multifrontal factorization encounters a dense frontal matrix at each node of the elimination tree. The solver needs to partially factorize this dense symmetric indefinite matrix using an APTP strategy that preserves numerical stability while identifying columns that cannot be safely eliminated (delayed pivots).

The caller provides a dense symmetric matrix (the frontal matrix) and receives back the L and D factors for the successfully eliminated columns, along with a list of delayed columns to be propagated to the parent node in the elimination tree.

**Why this priority**: This is the core numerical kernel. Without stable dense factorization, the entire sparse solver cannot function. Every other feature (multifrontal assembly, triangular solve, end-to-end API) depends on this kernel producing correct results.

**Independent Test**: Can be fully tested by providing dense symmetric matrices (positive definite, indefinite, near-singular) and verifying that the factorization reconstructs the input matrix within tolerance. Delivers standalone dense LDL^T factorization capability.

**Acceptance Scenarios**:

1. **Given** a dense symmetric positive definite matrix, **When** APTP factorization is applied, **Then** all columns are eliminated with 1x1 pivots, no columns are delayed, and reconstruction error ||A - LDL^T|| / ||A|| < 10^-12
2. **Given** a dense symmetric indefinite matrix, **When** APTP factorization is applied, **Then** some columns use 2x2 Bunch-Kaufman pivots, the stability bound |l_ij| < 1/threshold holds for all entries, and reconstruction error < 10^-12
3. **Given** a dense frontal matrix where some columns fail both 1x1 and 2x2 stability checks, **When** APTP factorization is applied with fallback-to-delay, **Then** those columns are reported as delayed and the factorization of the remaining columns is numerically correct
4. **Given** a dense frontal matrix from a sparse problem preprocessed with MC64 scaling, **When** APTP factorization is applied, **Then** the number of delayed columns is reduced compared to unscaled input (MC64 improves diagonal dominance)

---

### User Story 2 — Configurable Stability Threshold and Fallback Strategy (Priority: P2)

A solver developer needs to tune the APTP kernel's behavior for different problem classes. For well-conditioned problems, a stricter threshold yields more accurate results. For highly indefinite problems from interior-point optimization, a more permissive threshold with robust fallback strategies reduces delayed pivots and improves solver convergence.

**Why this priority**: Configurability directly affects the solver's applicability to real-world problem classes. The default configuration must work well out-of-the-box, but expert users need control over threshold and fallback behavior.

**Independent Test**: Can be tested by factoring the same matrix with different threshold values and fallback strategies, verifying that all configurations produce valid factorizations with varying trade-offs between accuracy and delayed column counts.

**Acceptance Scenarios**:

1. **Given** a threshold value of 0.01 (permitting growth factor up to 100), **When** factorization is applied, **Then** the stability bound |l_ij| < 100 is enforced for all L entries
2. **Given** a stricter threshold value of 0.5 (permitting growth factor up to 2), **When** factorization is applied to an indefinite matrix, **Then** more columns are delayed compared to threshold 0.01, but the resulting factors are more accurate
3. **Given** the Bunch-Kaufman fallback strategy, **When** a 1x1 pivot fails the stability check, **Then** the kernel attempts a 2x2 Bunch-Kaufman pivot using the partner column before resorting to delay
4. **Given** the immediate-delay fallback strategy, **When** a 1x1 pivot fails, **Then** the column is immediately marked as delayed without attempting a 2x2 pivot

---

### User Story 3 — Factorization Diagnostics for Solver Monitoring (Priority: P3)

A solver developer needs insight into the factorization process to diagnose convergence problems, tune parameters, and validate the solver pipeline. The APTP kernel should report summary statistics (pivot counts, delayed columns, worst stability metric) and optionally a per-column pivot log.

**Why this priority**: Diagnostics are essential for solver development and debugging but do not affect numerical correctness. They enable informed tuning decisions and help identify when preprocessing (ordering, scaling) needs improvement.

**Independent Test**: Can be tested by factoring matrices of known structure and verifying that reported statistics match expected pivot patterns (e.g., a positive definite matrix reports zero 2x2 pivots and zero delays).

**Acceptance Scenarios**:

1. **Given** a completed factorization, **When** statistics are queried, **Then** the counts of 1x1 pivots, 2x2 pivot pairs, and delayed columns are reported and sum to the matrix dimension
2. **Given** a completed factorization, **When** the worst stability metric is queried, **Then** the maximum absolute L entry across all columns is reported and satisfies the threshold bound
3. **Given** a factorization with 2x2 pivots, **When** inertia is computed from the D factor, **Then** the reported eigenvalue signs (positive, negative, zero) are correct for the input matrix

---

### Edge Cases

- What happens when the input matrix is 1x1? The kernel must handle the trivial case (single diagonal element as a 1x1 pivot, or delay if zero/near-zero).
- What happens when the input matrix is entirely zero? All columns should be reported as delayed (zero pivots).
- What happens when a 2x2 pivot block has a near-zero determinant? The 2x2 pivot should be rejected and both columns delayed.
- What happens when ALL columns must be delayed? The kernel must return an empty factorization (zero eliminated columns) with all columns in the delayed list — this is not an error condition, as the parent node will handle them.
- What happens when the threshold is set to an extreme value (e.g., 0.0 or 1.0)? Threshold = 0.0 should accept all pivots (no stability checking); threshold values approaching 1.0 should delay most pivots. Both are valid configurations.
- What happens when the matrix has repeated eigenvalues or clustered spectrum? The kernel must still produce a valid factorization; pivot decisions may differ from other algorithms but correctness must hold.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST factorize a dense symmetric matrix into LDL^T form where L is unit lower triangular and D is block diagonal with 1x1 and 2x2 blocks
- **FR-002**: System MUST implement the A Posteriori Threshold Pivoting strategy: factor optimistically assuming 1x1 pivots, then check stability a posteriori
- **FR-003**: System MUST enforce the stability bound |l_ij| < 1/threshold for all entries of L, where threshold is a configurable parameter
- **FR-004**: System MUST support at least two fallback strategies when a 1x1 pivot fails: (a) 2x2 Bunch-Kaufman pivot attempt followed by delay, and (b) immediate delay
- **FR-005**: System MUST correctly identify and report columns that cannot be stably eliminated (delayed columns), returning them to the caller for handling
- **FR-006**: System MUST store the D factor using the existing MixedDiagonal type (from Phase 2), preserving mixed 1x1/2x2 block structure
- **FR-007**: System MUST store pivot classification using the existing PivotType enum (from Phase 2: OneByOne, TwoByTwo, Delayed)
- **FR-008**: System MUST compute the Schur complement update after eliminating each pivot to update the trailing submatrix
- **FR-009**: System MUST handle partial factorization — only the fully-summed portion of a frontal matrix is factorized; the contribution block is updated but not factorized
- **FR-010**: System MUST accept a configurable stability threshold with a sensible default value (0.01, following SPRAL convention)
- **FR-011**: System MUST provide factorization summary statistics: count of 1x1 pivots, 2x2 pivot pairs, delayed columns, and maximum absolute L entry
- **FR-012**: System MUST produce factorizations from which inertia can be correctly computed using the existing Inertia type and MixedDiagonal::compute_inertia()

### Key Entities

- **AptpFactorization**: The result of applying APTP to a dense matrix. Contains the L factor, D factor (MixedDiagonal), list of delayed columns, and factorization statistics. Represents a partial LDL^T decomposition where some columns may be uneliminated.
- **AptpOptions**: Configuration controlling the kernel's behavior. Contains the stability threshold and fallback strategy selection. Provides sensible defaults matching SPRAL conventions.
- **AptpStatistics**: Summary of what happened during factorization. Includes pivot type counts and worst-case stability metric. Used for solver monitoring and parameter tuning.
- **Pivot Log Record**: Per-column diagnostic recording what pivot type was used, the stability metric achieved, and whether fallback was invoked. Optional diagnostic output for development and debugging.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: All 15 hand-constructed test matrices factor with reconstruction error < 10^-12
- **SC-002**: 100+ randomly generated symmetric indefinite matrices factor with reconstruction error < 10^-12 and stability bound satisfied
- **SC-003**: All SuiteSparse CI-subset matrices (10 matrices) whose dense frontal submatrices are extracted: factorizations satisfy reconstruction error < 10^-12
- **SC-004**: Positive definite matrices produce zero delayed columns and zero 2x2 pivots in 100% of test cases
- **SC-005**: All fallback strategy variants produce valid factorizations on the same test suite — results may differ in delayed column counts but all satisfy reconstruction tolerance
- **SC-006**: Factorization statistics (pivot counts + delayed count) sum to the matrix dimension for every test case
- **SC-007**: Inertia computed from the D factor matches analytically known inertia for all hand-constructed test matrices
- **SC-008**: The worst stability metric (max |l_ij|) reported in statistics matches the actual maximum across all L entries, verified by independent scan on 100+ test matrices

## Algorithm References

The APTP algorithm is defined in the following academic literature, available as reference Markdown files:

| Reference                            | Location                                                     | Key Content                                                                                                                    |
|--------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------|
| Duff, Hogg & Lopez (2020)            | `/workspace/rivrs-linalg/references/ssids/duff2020.md`       | Core APTP algorithm (Algorithm 3.1), fail-in-place approach, 2D block partitioning, stability analysis (Section 4), complete pivoting inner kernel (Algorithm 4.1) |
| Hogg, Ovtchinnikov & Scott (2016)    | `/workspace/rivrs-linalg/references/ssids/hogg2016.md`       | SSIDS v1 design, threshold partial pivoting, delayed pivot mechanism, GPU factorization                                        |
| Duff & Reid (1983)                   | `/workspace/rivrs-linalg/references/ssids/duff1984.md`       | Foundational multifrontal method for indefinite systems, frontal matrix structure, contribution blocks                          |
| Liu (1992)                           | `/workspace/rivrs-linalg/references/ssids/liu1992.md`        | Multifrontal theory, elimination trees, supernode concepts                                                                     |
| Duff & Pralet (2005)                 | `/workspace/rivrs-linalg/references/ssids/duff2005.md`       | Strategies for scaling and pivoting in indefinite symmetric systems                                                            |
| Schenk, Gartner & Karypis (2006)     | `/workspace/rivrs-linalg/references/ssids/schenk2006.md`     | Matching-based preprocessing for LDL^T (context for MC64 integration)                                                          |

### APTP Algorithm Summary (from Duff, Hogg & Lopez 2020)

The APTP strategy applies to the dense partial factorization of each frontal matrix in a multifrontal solver. For a frontal matrix F partitioned as:

```
F = [ F_11  F_21^T ]
    [ F_21  F_22   ]
```

where F_11 (p x p) contains the fully-summed variables and F_21 ((m-p) x p) connects to non-fully-summed variables:

1. **Optimistic factorization**: Factor columns of F_11 assuming 1x1 pivots, producing L and D entries
2. **A posteriori stability check**: After factoring each column (or block of columns), verify that all entries satisfy |l_ij| < 1/threshold
3. **On failure**: Either attempt a 2x2 Bunch-Kaufman pivot with the adjacent column, or mark the column as "delayed" (fail-in-place)
4. **Schur complement update**: Update the trailing submatrix (F_22 and remaining F_11 columns) using the successfully eliminated pivots
5. **Delayed column handling**: Failed columns remain in place and are updated during subsequent elimination steps. After all passes, they are permuted to the end and returned to the caller.

The inner-block factorization (Algorithm 4.1 in the paper) uses complete pivoting on small blocks (order ib, typically 32), which provides stability equivalent to threshold partial pivoting with u = 0.25.

## Assumptions

- The input matrix is provided as a dense matrix (not sparse CSC). In the multifrontal context, frontal matrices are assembled as dense from sparse contributions before calling this kernel.
- The input matrix is symmetric. Only the lower triangle is read; the upper triangle may contain arbitrary values.
- The caller (Phase 6 multifrontal factorization) is responsible for applying MC64 scaling before calling this kernel. The kernel operates on the already-scaled matrix.
- The caller is responsible for propagating delayed columns to parent nodes in the elimination tree. The kernel's responsibility ends at identifying which columns are delayed.
- Phase 2 types (MixedDiagonal, PivotType, Block2x2, Inertia) are stable and will not change during Phase 5 development.
- faer 0.22 dense operations (matrix multiplication, triangular solve) are available and correct.
- Default threshold value of 0.01 matches SPRAL's default and provides good balance between stability and pivot acceptance for typical problems.
- The two-level blocking optimization (outer block size nb, inner block size ib as described in the paper) is deferred to Phase 9.1 (Two-Level APTP). Phase 5 implements single-level APTP where the entire frontal matrix is processed as one block.

## Dependencies

- **Phase 2 (APTP Data Structures)**: MixedDiagonal, PivotType, Block2x2, Inertia types — all complete
- **Phase 3 (AptpSymbolic)**: Symbolic analysis providing elimination tree structure — complete (used by Phase 6, not directly by this kernel)
- **Phase 4 (Ordering & Preprocessing)**: METIS ordering + MC64 scaling — complete (applied before this kernel is called)
- **Phase 0.4/1.1 (Test Infrastructure)**: NumericalValidator, test matrix registry, random generators — complete
- **faer 0.22**: Dense matrix types (Mat, MatRef, MatMut), dense operations (matrix multiplication, triangular solve), permutation types
