# Feature Specification: Supernode Amalgamation

**Feature Branch**: `020-supernode-amalgamation`
**Created**: 2026-02-21
**Status**: Draft
**Input**: User description: "Phase 9.1a — SPRAL-style nemin-based merge pass to reduce supernode count for narrow-supernode matrices"

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Solve Schenk Optimization Matrices at Competitive Speed (Priority: P1)

A user solving large Schenk-family optimization problems (c-71, c-big) observes that factorization takes 11-25x longer than SPRAL. The root cause is that faer's symbolic analysis produces ~35K narrow supernodes (avg 2.17 columns) where SPRAL produces ~8K amalgamated supernodes. The user needs the solver to automatically merge small supernodes so that assembly/extraction overhead does not dominate factorization time.

**Why this priority**: This is the primary motivating problem. c-71 and c-big are the worst-performing matrices in the entire SuiteSparse benchmark suite, and the root cause is directly traceable to excessive supernode count. Without amalgamation, 78% of factorization time is assembly + extraction overhead rather than the dense APTP kernel.

**Independent Test**: Can be tested by factorizing c-71 and c-big with the amalgamation pass enabled and measuring supernode count reduction (target: 3-5x fewer supernodes) and factor time improvement (target: 5-15x speedup on these matrices).

**Acceptance Scenarios**:

1. **Given** a Schenk optimization matrix (c-71: n=76638) with ~35K fundamental supernodes, **When** symbolic analysis runs with default settings, **Then** the amalgamation pass reduces supernode count to fewer than 12K supernodes.
2. **Given** the amalgamated supernode structure for c-71, **When** numeric factorization runs, **Then** factor time is within 2x of SPRAL's reference time (currently 24.5x).
3. **Given** any amalgamated factorization, **When** the solve phase completes, **Then** backward error remains below 5e-11 (no accuracy regression).

---

### User Story 2 — Transparent Amalgamation with No Regressions (Priority: P1)

A user solving a variety of sparse symmetric indefinite problems (the full 65-matrix SuiteSparse collection) needs amalgamation to improve narrow-supernode matrices without degrading performance or accuracy on any other matrix.

**Why this priority**: Equal priority with US1 because a regression on any matrix would block release. The amalgamation must be safe by default.

**Independent Test**: Can be tested by running the full 65-matrix SuiteSparse benchmark suite and verifying that no matrix has worse backward error or factor time than the pre-amalgamation baseline.

**Acceptance Scenarios**:

1. **Given** the full 65-matrix SuiteSparse suite, **When** factorized with amalgamation enabled (default), **Then** all 65 matrices produce backward error below 5e-11.
2. **Given** the full SuiteSparse suite, **When** factorized with amalgamation enabled, **Then** no matrix has factor time more than 10% slower than the pre-amalgamation baseline.
3. **Given** a matrix where faer's fundamental supernodes are already large (e.g., ldoor, inline_1), **When** amalgamation runs, **Then** the pass completes with minimal or no merges and negligible overhead.

---

### User Story 3 — Configurable Amalgamation Threshold (Priority: P2)

A library user with domain knowledge about their problem structure wants to control the amalgamation aggressiveness — for example, using a larger threshold for 3D FEM problems with bushy trees, or disabling amalgamation entirely for debugging or comparison purposes.

**Why this priority**: Important for power users and benchmarking, but the default (nemin=32, matching SPRAL) should work well for the vast majority of problems.

**Independent Test**: Can be tested by factorizing a matrix with different nemin values (e.g., 8, 32, 64) and observing that supernode count decreases and factor time changes predictably.

**Acceptance Scenarios**:

1. **Given** a user sets the amalgamation threshold to 64, **When** symbolic analysis runs, **Then** more aggressive merging occurs (fewer, larger supernodes) compared to the default of 32.
2. **Given** a user disables amalgamation (threshold = 1 or explicit disable), **When** symbolic analysis runs, **Then** the result matches the current (pre-amalgamation) behavior exactly.
3. **Given** a user sets the amalgamation threshold to 1 (disabled), **When** solving any matrix, **Then** results are bitwise identical to the pre-amalgamation solver.

---

### Edge Cases

- What happens when all supernodes are already larger than nemin? The amalgamation pass should be a no-op with negligible overhead.
- What happens with simplicial decomposition (every column is its own supernode, e.g., bloweybq)? Amalgamation should merge 1-column supernodes into larger groups, but this is not the primary target of Phase 9.1a (that is Phase 9.1b's workspace reuse). Amalgamation may still help reduce overhead for these matrices.
- What happens when a parent has many children (bushy assembly tree)? The merge pass must consider all children, not just the consecutively-numbered one (this is exactly faer's limitation that motivates this work).
- What happens when merging introduces significant fill-in? The nemin-based strategy merges unconditionally when both nodes are small, accepting some fill-in. SPRAL tracks explicit zeros (`ezero`) for statistics but does not use a fill-in limit for nemin merges.
- What happens with matrices that have a single root supernode with many children (star-shaped trees)? The merge pass should handle arbitrary tree shapes correctly.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The solver MUST perform a supernode amalgamation pass after faer's symbolic analysis and before numeric factorization, merging parent-child supernode pairs that are both smaller than a configurable threshold (`nemin`).
- **FR-002**: The amalgamation pass MUST consider ALL children of each parent supernode for merging (not restricted to consecutively-numbered children as faer's relaxed merging is).
- **FR-003**: Two supernodes MUST be merged when EITHER: (a) the parent has 1 eliminated column and the parent's column count equals the child's column count minus 1 (structural match — these are supernodal by the fundamental definition), OR (b) both parent and child have fewer than `nemin` eliminated columns (the nemin condition).
- **FR-004**: The default value of `nemin` MUST be 32, matching SPRAL's default (`nemin_default = 32` in `datatypes.f90:21`).
- **FR-005**: The `nemin` parameter MUST be exposed on `AnalyzeOptions` (and propagated through `SolverOptions`) so users can tune or disable amalgamation.
- **FR-006**: The amalgamation pass MUST produce valid supernode metadata: begin/end column ranges, row patterns (with fill entries for merged structure), parent pointers, and column counts that are consistent with the merged structure.
- **FR-007**: The amalgamation pass MUST maintain the postorder property of the assembly tree (children processed before parents in factorization).
- **FR-008**: The row pattern for a merged supernode MUST be the union of the constituent supernodes' row patterns, with additional fill entries as needed for the rectangular block structure.
- **FR-009**: The numeric factorization MUST consume the amalgamated supernode structure without changes to the dense APTP kernel or solve path.
- **FR-010**: Backward error on all 65 SuiteSparse matrices MUST remain below 5e-11 after amalgamation.
- **FR-011**: Factor time on c-71 and c-big MUST improve by at least 5x compared to the pre-amalgamation baseline.
- **FR-012**: Factor time on all other SuiteSparse matrices MUST not regress by more than 10% compared to the pre-amalgamation baseline.

### Key Entities

- **AmalgamatedSupernode**: A supernode that may contain columns from multiple faer fundamental supernodes that have been merged. Characterized by: column range (begin, end), row pattern (union of constituent patterns with fill), parent pointer, number of eliminated columns.
- **Amalgamation configuration (nemin)**: The minimum eliminated column count threshold. Parent-child pairs where both have fewer than `nemin` eliminated columns are merged unconditionally. Default: 32.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: c-71 supernode count reduced from ~35K to fewer than 12K (matching the 3-5x reduction ratio observed in SPRAL's 7,697 supernodes).
- **SC-002**: c-71 and c-big factor times within 2x of SPRAL reference times (from current 11-25x).
- **SC-003**: All 65 SuiteSparse matrices maintain backward error below 5e-11.
- **SC-004**: No SuiteSparse matrix factor time regresses by more than 10% vs pre-amalgamation baseline.
- **SC-005**: Amalgamation pass overhead is less than 5% of total symbolic analysis time for any matrix.

## Assumptions

- faer's `SymbolicCholesky` result (fundamental supernodes, elimination tree, column counts, row patterns) provides sufficient information to implement a post-hoc amalgamation pass without modifying faer.
- The existing `build_supernode_info()` function in `numeric.rs` is the natural integration point — amalgamation transforms the `Vec<SupernodeInfo>` before factorization consumes it.
- The dense APTP kernel (`aptp_factor_in_place`) and solve path (`aptp_solve`) operate on per-supernode data and do not need to know whether supernodes were amalgamated.
- SPRAL's `nemin=32` default is appropriate for the target workloads (Schenk optimization matrices from 3D nested dissection).
- The fill-in increase from nemin-based merging is acceptable for performance-critical matrices (SPRAL uses the same strategy and the approach is well-validated in the literature).

## Algorithm References

The following academic references and source code inform this feature:

### Primary Algorithm Reference
- **SPRAL `core_analyse.f90:528-853`** — `find_supernodes()` and `do_merge()` functions. The `do_merge()` predicate (lines 806-822) implements the two-condition merge: structural match OR both-nodes-small. The `merge_nodes()` subroutine (lines 827-853) handles the actual merge bookkeeping including explicit zero tracking. Available at `/workspace/rivrs-linalg/references/spral/src/core_analyse.f90`.

### faer Limitation (Motivating This Work)
- **faer `cholesky.rs:2461`** — The consecutivity check `child_index + 1 == parent_index` blocks most merges for bushy 3D assembly trees. Documented in ssids-plan.md Phase 9.1a analysis.

### Academic Papers (Markdown summaries in `/workspace/rivrs-linalg/references/ssids/`)
- **Liu (1992)** — "The Multifrontal Method for Sparse Matrix Solution: Theory and Practice" (SIAM Review). Defines fundamental supernodes, assembly trees, and the extend-add operation. File: `liu1992.md`
- **Hogg, Scott & Sherwood-Jones (2016)** — "A New Sparse LDLT Solver using A Posteriori Threshold Pivoting". Describes SSIDS's supernodal structure and `nemin` parameter (default 32 for GPU, 64 tested for CPU). File: `hogg2016.md`
- **Duff, Hogg & Lopez (2020)** — "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting" (SIAM J. Sci. Comput.). Two-level APTP architecture and its interaction with supernodal structure. File: `duff2020.md`
- **Schenk & Gartner (2006)** — "On fast factorization pivoting methods for sparse symmetric indefinite systems". Discusses relaxed supernodes and the cost-benefit of amalgamation for BLAS-3 utilization. File: `schenk2006.md`
- **Gilbert, Ng & Peyton (1992)** — "An efficient algorithm to compute row and column counts for sparse Cholesky factorization". Column count prediction used in fill analysis. File: `gilbert1992.md`

### SPRAL Configuration Reference
- `datatypes.f90:21` — `nemin_default = 32`
- `datatypes.f90:213` — `options%nemin` field definition
- `anal.F90:981-988` — nemin validation and propagation to `find_supernodes()`
