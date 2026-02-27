# Research: Robustness — Testing & Hardening

**Feature**: 026-robustness-hardening
**Date**: 2026-02-25

## Decision 1: SPRAL Test Scenario Inventory

### Decision
SPRAL's kernel test suites contain **77+ deterministic tests** and **3000 torture test instances** across two files:

**ldlt_app.cxx** (APP/complete pivoting kernel):
- 48 deterministic tests in 8 categories (A–H): small block, aggressive, two-level blocking, delays, square, singular diagonal blocks
- Matrix sizes: (2,1) through (128,128), both rectangular (m>n) and square (m=n)
- 2 torture test suites: 500 instances each at (128,128) and (128,48)
- Perturbation helpers: `cause_delays()`, `make_singular()`, `make_dblk_singular()`
- Backward error threshold: 1e-13
- Key assertion: q1+q2 = m (all rows fully processed)

**ldlt_tpp.cxx** (TPP/scalar kernel):
- 23 deterministic tests in 4 categories (A–D): simple, delays, singular, singular+delays
- Matrix sizes: (1,1) through (500,500)
- 2 torture test suites: 1000 instances each at (100,100) and (100,50)
- Perturbation helpers: `cause_delays()`, `make_singular()` (no `make_dblk_singular`)
- Backward error threshold: 2e-13 (slightly relaxed vs APP)
- Extra assertion: L growth bound ≤ 1/u

### Rationale
Complete inventory needed before mapping to rivrs-sparse coverage. The two files test different kernels that rivrs-sparse maps to `aptp_factor_in_place` (complete pivoting path for blocks ≥ 32) and `tpp_factor_as_primary` (TPP path for blocks < 32).

### Alternatives Considered
- Auditing only torture tests (insufficient — deterministic tests cover specific edge cases)
- Auditing full SPRAL test suite beyond kernels (out of scope — multifrontal/solve tests are integration-level, covered by SuiteSparse suite)

## Decision 2: Perturbation Helper Design

### Decision
Three perturbation helpers operating on dense `Mat<f64>`, matching SPRAL's conceptual approach:

1. **`cause_delays(matrix, rng)`**: Select n/8 random rows, multiply by 1000. Also multiply n/8 random individual entries by 1000. Ensures first oversized row is past first block (BLOCK_SIZE rows) when n > BLOCK_SIZE. Forces pivot threshold failures → delayed columns.

2. **`make_singular(matrix, col1, col2)`**: Make col2 a scaled copy of col1 (read col1, scale, write to col2, symmetrize). Creates exact rank deficiency → tests delayed column handling and inertia tracking.

3. **`make_dblk_singular(matrix, block_row, block_size)`**: Select first and last columns of a diagonal block, call `make_singular`. Creates singular diagonal block → tests 2×2 pivot detection within blocks.

### Rationale
SPRAL's perturbation helpers (ldlt_app.cxx:78-134) are the proven approach for creating adversarial test matrices that exercise the full pivot decision tree. Our helpers follow the same *concept* but are implemented from scratch in Rust (clean room).

### Alternatives Considered
- Random perturbation (too unpredictable, doesn't target specific code paths)
- Hand-constructed adversarial matrices only (insufficient coverage breadth)

## Decision 3: Proptest Integration Approach

### Decision
Add `proptest = "1.4"` as a dev-dependency. Create two new modules:

- `src/testing/strategies.rs`: Custom proptest strategies for generating symmetric matrices (PD and indefinite) at configurable sizes/densities. Wraps existing `generate_random_symmetric` infrastructure.
- `src/testing/perturbations.rs`: The three perturbation helpers (shared between proptest properties and torture tests).

Property tests go in `src/aptp/factor.rs` `#[cfg(test)]` module. End-to-end property tests via `SparseLDLT` go in a new integration test file.

### Rationale
Proptest is the Rust ecosystem standard for property-based testing. Existing random matrix generators provide a solid foundation. Custom strategies wrap them for proptest's shrinking/replay infrastructure.

### Alternatives Considered
- `quickcheck` (less flexible shrinking, smaller ecosystem)
- Manual random testing without framework (no automatic shrinking on failure, no replay from seeds)

## Decision 4: Torture Test vs Property Test Separation

### Decision
**Torture tests** (SPRAL-style, FR-002–FR-004):
- Dense kernel level only (`aptp_factor_in_place`, `tpp_factor_as_primary`)
- Fixed seed for reproducibility, 500+ instances per config
- `#[ignore]` (long-running, run explicitly)
- Test probabilistic perturbation mix: 70% delays, 20% singular, 10% dblk_singular
- Verify: no panics, backward error < 5e-11 for non-singular, correct inertia for singular

**Property tests** (proptest, FR-005):
- End-to-end through `SparseLDLT` (sparse matrices, full pipeline)
- Proptest-managed seeds with automatic shrinking
- Run with normal `cargo test` (not `#[ignore]`)
- Properties: backward error, inertia sum = n, permutation validity, no panics
- Size range 5–500, configurable case count (default 256 for CI speed, 1000+ for thorough)

### Rationale
Torture tests match SPRAL's approach (kernel-level, high instance count, specific perturbations). Property tests complement with broader parameter space and automatic shrinking. Different `#[ignore]` policies reflect different runtime characteristics.

### Alternatives Considered
- Combining both into one framework (loses the targeted kernel-level coverage of torture tests)
- Property tests at kernel level only (misses integration-level issues)

## Decision 5: Adversarial Input Handling

### Decision
Test the following adversarial inputs through `SparseLDLT::solve_full` and `SparseLDLT::analyze_with_matrix`:

| Input | Expected Behavior |
|-------|------------------|
| 0×0 empty matrix | Clean error (not panic) |
| 1×1 with nonzero diagonal | Correct solve |
| 1×1 with zero diagonal | Clean error or zero-inertia factorization |
| Pure diagonal matrix | Trivial correct factorization |
| Arrowhead matrix | Correct factorization |
| Near-overflow values (1e308) | Clean error or correct result |
| Near-underflow values (1e-308) | Clean error or correct result |
| Matrix with NaN entries | Clean error (not panic) |
| Matrix with Inf entries | Clean error (not panic) |
| Non-symmetric sparsity pattern | Clean error |

Some of these may require new guard code in the solver (e.g., NaN/Inf detection). The spec requires zero panics (FR-007).

### Rationale
Production robustness requires graceful handling of all inputs. Current solver may panic on some of these (e.g., 0×0 matrix likely triggers index-out-of-bounds). These tests will drive defensive code additions.

### Alternatives Considered
- Testing only well-formed inputs (insufficient for production robustness)
- Fuzzing with cargo-fuzz (complementary but doesn't replace structured adversarial tests — could be added post-release)

## Decision 6: Test Suite Pruning Criteria

### Decision
During the audit, evaluate each existing test against these criteria for removal:

**Remove if**:
- Test is redundant: another test covers the exact same code path with the same or broader assertions
- Test is TDD scaffolding: tests an intermediate internal detail (e.g., a helper function's return shape) that is now fully covered by higher-level tests
- Test has no independent assertion value: it only checks "doesn't crash" on a case where a neighboring test already checks correctness

**Keep if**:
- Test covers a meaningful scenario even if SPRAL doesn't test it (e.g., MC64 scaling validation, parallel determinism, workspace reuse)
- Test is the only coverage for a specific code path or edge case
- Test provides regression protection for a previously-found bug
- Test validates a structural property (inertia consistency, permutation validity)

### Rationale
SPRAL parity is a floor, not a ceiling (per spec clarification). The goal is a lean, high-signal suite — not blindly matching or trimming to SPRAL's scope.

### Alternatives Considered
- No pruning (keeps dead weight from TDD phases)
- Aggressive pruning to match SPRAL exactly (loses valuable rivrs-specific coverage)

## SPRAL-to-rivrs Architecture Mapping

Understanding how SPRAL's test structure maps to rivrs-sparse's architecture:

| SPRAL Concept | rivrs-sparse Equivalent | Notes |
|---------------|------------------------|-------|
| `ldlt_app` (APP block kernel) | `aptp_factor_in_place` with complete pivoting (num_fully_summed ≥ inner_block_size) | rivrs uses `factor_inner` with BLAS-3 pipeline |
| `ldlt_tpp` (TPP scalar kernel) | `tpp_factor_as_primary` (num_fully_summed < inner_block_size) | Called from `aptp_factor_in_place` |
| `BLOCK_SIZE` template param | `inner_block_size` in `AptpOptions` (default 32) | rivrs uses runtime config, not template |
| `app_block` vs `app_aggressive` | Not directly mapped | rivrs only implements `app_block` equivalent; aggressive (Cholesky-like) not implemented — N/A |
| `outer_block_size` (two-level) | `outer_block_size` in `AptpOptions` (default 256) | Same concept |
| `do_update` (Schur complement) | `compute_contribution_gemm` in `numeric.rs` | Different architecture (multifrontal vs standalone kernel test) |
| `permute_rows` (row permutation) | Physical row/column swaps in `factor_inner` | rivrs swaps in-place rather than permute-on-read |
| Backward error check | `dense_reconstruction_error()` + `sparse_backward_error()` | rivrs uses reconstruction as primary oracle |
