# Implementation Plan: Parallel Factorization & Solve (Phase 8.2)

**Branch**: `019-parallel-factorization-solve` | **Date**: 2026-02-21 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/019-parallel-factorization-solve/spec.md`

## Summary

Add shared-memory parallelism to the SSIDS multifrontal LDL^T solver, targeting two complementary strategies: (1) intra-node BLAS-3 parallelism for large dense fronts via faer's `Par` enum, and (2) tree-level parallelism for independent subtrees via Rayon's `scope`-based fork-join. The `Par` enum is added to `FactorOptions` and `SolverOptions` with `Par::Seq` as default for full backward compatibility. Phase 8.1g profiling shows 63% of SuiteSparse matrices are IntraNode-dominated (top 10% of fronts account for >80% of factor time), making intra-node BLAS parallelism the highest-impact first step. Tree-level scheduling benefits the remaining 29% (Mixed class) and 8% (TreeLevel class).

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (dense LA + sparse infrastructure), rayon 1.x (NEW — tree-level parallelism), metis-sys 0.3 (ordering)
**Storage**: N/A (in-memory computation)
**Testing**: cargo test (358 unit tests + 10 CI SuiteSparse integration tests)
**Target Platform**: Linux (Docker dev environment, >=4 cores)
**Project Type**: Single Rust library crate
**Performance Goals**: >=3x speedup on 4 cores for IntraNode-class matrices (SC-001/SC-002), >=1.5x for Mixed-class (SC-003), <5% overhead for small matrices (SC-007)
**Constraints**: Bitwise determinism at fixed thread count (SC-004), backward error < 5e-11 for all 65 SuiteSparse matrices (SC-006)
**Scale/Scope**: 65 SuiteSparse test matrices (n: 420–923,136, nnz: 1,572–7,309,820). Factorization time ranges from 17ms to 112s sequential.

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Parallel results validated against same backward error threshold; all 65 matrices must pass. Determinism at fixed thread count ensures reproducibility. |
| II. Clean Room | PASS | Parallelism patterns from academic papers (Duff 2020 Section 5-6) and SPRAL source (BSD-3). Rayon is MIT-licensed. No restricted sources consulted. |
| III. TDD | PASS | Tests will verify: backward error < 5e-11, determinism (run twice compare), regression against sequential, correctness at each step. |
| IV. Documentation | PASS | Algorithm references documented in spec. All parallel patterns will cite academic sources and SPRAL equivalents. |
| V. Numerical Stability | PASS | Parallelism does not change the APTP algorithm — same pivot decisions, same factorization (at fixed thread count). Relaxed cross-thread-count identity acknowledges BLAS summation order differences. |
| VI. Structured Development | PASS | Phase 8.2 follows Phase 8.1g completion. Implementation order matches profiling recommendations. |
| VII. Code Quality | PASS | Uses faer's `Par` enum (transparent composition). Rayon is idiomatic Rust. No unsafe code needed. |

**Post-design re-check**: PASS — all principles satisfied. No violations.

## Project Structure

### Documentation (this feature)

```text
specs/019-parallel-factorization-solve/
├── plan.md              # This file
├── spec.md              # Feature specification
├── research.md          # Phase 0 research findings
├── data-model.md        # Data model changes
├── quickstart.md        # Implementation quickstart guide
├── contracts/           # API contract changes
│   └── api-changes.md   # Public and internal API modifications
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code (modified files)

```text
src/
├── aptp/
│   ├── solver.rs        # Add par to FactorOptions, SolverOptions; propagate through SparseLDLT
│   ├── factor.rs        # Add par to AptpOptions; replace Par::Seq in TRSM (L1281) and GEMM (L1426)
│   ├── numeric.rs       # Refactor postorder loop → recursive rayon::scope; per-task global_to_local
│   └── solve.rs         # Add par parameter; parallel diagonal/forward/backward solve
├── Cargo.toml           # Add rayon = "1"
examples/
└── parallel_scaling.rs  # NEW: benchmark tool for parallel scaling report
```

**Structure Decision**: All changes are within the existing `src/aptp/` module. No new modules needed. One new example binary for benchmarking.

## Implementation Phases

### Phase A: API Surface + Rayon Dependency (Low risk, foundation)

**Goal**: Add `par: Par` to all option structs and propagate through the call chain without changing behavior. All existing tests pass with `Par::Seq` default.

**Changes**:
1. Add `rayon = "1"` to `Cargo.toml`
2. Add `par: Par` field to `FactorOptions`, `SolverOptions`, `AptpOptions` (default `Par::Seq`)
3. Add `par: Par` parameter to `SparseLDLT::solve_in_place` and `aptp_solve`
4. Thread `par` from `FactorOptions` → `AptpOptions` in `SparseLDLT::factor`
5. Thread `par` from `aptp_solve` → per-supernode solve functions
6. Update all call sites (tests, examples, benchmarks) with new signatures

**Tests**: All existing 358 tests pass. No behavioral change.
**Risk**: Low — pure refactoring, no parallelism activated.

### Phase B: Intra-Node BLAS-3 Parallelism (Medium risk, highest impact)

**Goal**: Replace hardcoded `Par::Seq` in dense APTP kernel with `options.par`, gated by front size threshold.

**Changes**:
1. In `AptpNumeric::factor()`: compute `effective_par` per supernode — `Par::Seq` if front dimension < 256, else `options.par`
2. In `apply_and_check()` (factor.rs): replace `Par::Seq` with `par` parameter at TRSM call (line 1278)
3. In `update_trailing()` (factor.rs): replace `Par::Seq` with `par` parameter at GEMM call (line 1426)
4. Thread `par` through `factor_inner()` → `two_level_factor()` → `apply_and_check()` + `update_trailing()`
5. In per-supernode solve functions: apply same threshold logic for dense triangular solve and matmul

**Tests**:
- Regression: all 358 tests pass with `Par::Seq`
- New: SuiteSparse CI subset with `Par::rayon(4)` — backward error < 5e-11
- New: determinism test — factor twice with same Par, compare bitwise
- New: small-matrix overhead test — time Par::Seq vs Par::Rayon on small problems

**Risk**: Medium — faer's parallel BLAS may produce different ULP-level results. Must verify backward error, not bitwise identity across thread counts.

### Phase C: Tree-Level Parallel Factorization (Higher risk, broader coverage)

**Goal**: Refactor sequential postorder loop into recursive function with `rayon::scope` at tree forks.

**Changes**:
1. Refactor `AptpNumeric::factor()` loop body into `fn factor_node(s, ...)` helper
2. Add recursive tree traversal: `factor_subtree(root, symbolic, ...)` that spawns children in `rayon::scope`
3. Replace shared `global_to_local` with per-task allocation (each parallel supernode builds its own)
4. Replace sequential `front_factors_vec.push()` with pre-allocated `Vec<Option<FrontFactors>>` indexed by supernode
5. Replace sequential `per_sn_stats.push()` with pre-allocated `Vec<Option<PerSupernodeStats>>` indexed by supernode
6. Compute aggregate `FactorizationStats` from collected per-supernode stats after tree traversal
7. When `par` is `Par::Seq`: use sequential recursive traversal (no rayon overhead)
8. When `par` is `Par::Rayon`: use `rayon::scope` for children with >1 child

**Key design decisions**:
- **Return-value pattern for contributions**: `factor_supernode()` returns `(FrontFactors, Option<ContributionBlock>, PerSupernodeStats)` as owned values. The recursive `factor_subtree()` collects children's returned contributions after the scope barrier, then passes them to the parent's assembly. No shared mutable `contributions` Vec is needed — all data flows through return values and scope barriers. This avoids `unsafe` code entirely and leverages Rust's ownership system for correctness.
- `global_to_local` cannot be shared (each supernode builds/resets its own range) → per-task `Vec<usize>` allocation
- Root supernodes (no parent): no contribution block, just front factors

**Tests**:
- Regression: sequential recursive path matches original postorder loop
- New: all 65 SuiteSparse matrices with `Par::rayon(4)` — backward error < 5e-11
- New: determinism — run twice, compare bitwise at fixed thread count
- New: Mixed/TreeLevel class matrices show measurable speedup

**FR-008 (no data races) coverage**: All parallel code uses safe Rust — no `unsafe`, no `UnsafeCell`, no raw pointers. Rust's type system statically prevents data races at compile time. The return-value pattern for contributions avoids shared mutable state entirely. A dedicated ThreadSanitizer task is unnecessary; compilation is sufficient proof.

**Risk**: Higher — Rust's ownership system will enforce correctness of shared data, but the refactoring of the factorization loop is significant. The return-value pattern eliminates shared mutable state for contributions, reducing the risk surface.

### Phase D: Parallel Solve (Medium risk, moderate impact)

**Goal**: Parallelize the three-phase solve (forward, diagonal, backward).

**Changes**:
1. **Diagonal solve**: Use `rayon::par_chunks` or `par_iter` over supernodes — each reads/writes only its own `col_indices` (no overlaps). Per-task workspace.
2. **Forward solve**: Recursive tree traversal with `rayon::scope`, matching factorization's tree structure. Children scatter to `rhs[row_indices]` before parent reads — scope barrier ensures ordering.
3. **Backward solve**: Reverse tree traversal with `rayon::scope`. Parent updates `rhs[col_indices]` before children read — reverse scope barrier.
4. Per-task workspace allocation (work/work2 buffers) for parallel supernodes.

**Key design decisions**:
- Forward solve scatter: siblings may have overlapping `row_indices` (shared parent's column space) → scope processes children first (independently), then parent (after barrier). No sibling conflicts because siblings don't share row_indices.
- Workspace: allocate per-task rather than trying to pool (simplicity over allocation optimization)

**Tests**:
- Regression: solve with `Par::Seq` matches existing output
- New: solve with `Par::rayon(4)` — backward error < 5e-11
- New: timing comparison for diagonal solve (should show improvement)
- New: timing comparison for forward/backward on wide-tree matrices

**Risk**: Medium — the RHS vector is shared mutable state. Must verify no data races via scope barrier correctness.

### Phase E: Benchmark Tool & Scaling Report (Low risk, validation)

**Goal**: Build `examples/parallel_scaling.rs` that produces a structured parallel scaling report.

**Changes**:
1. New example: `parallel_scaling.rs` (requires `diagnostic` feature)
2. Run factorization and solve at 1, 2, 4, 8 threads per matrix
3. Compute speedup = T_seq / T_par, efficiency = speedup / n_threads
4. Classify by workload type (IntraNode/Mixed/TreeLevel)
5. Output structured JSON report + human-readable summary
6. Compare against Phase 8.1g sequential baselines

**Tests**: Verify report structure contains required fields (per SC-008).
**Risk**: Low — benchmarking infrastructure, no numerical changes.

## Complexity Tracking

No constitution violations — no complexity tracking needed.
