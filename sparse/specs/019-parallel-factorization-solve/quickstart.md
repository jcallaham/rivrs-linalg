# Quickstart: Parallel Factorization & Solve (Phase 8.2)

## Implementation Order

### Step 1: Add `par` to options structs and Rayon dependency
- Add `rayon = "1"` to `Cargo.toml`
- Add `par: Par` field to `FactorOptions`, `SolverOptions`, `AptpOptions`
- Default to `Par::Seq` everywhere
- Update `SparseLDLT::factor` and `solve_in_place` to accept and propagate `par`
- **Tests**: Existing tests should pass unchanged (all using default `Par::Seq`)

### Step 2: Thread `par` through dense APTP kernel
- Replace `Par::Seq` at TRSM (factor.rs:1281) and GEMM (factor.rs:1426) with `options.par`
- Apply front-size threshold: fronts < 256 get `Par::Seq` regardless
- **Tests**: Run existing SuiteSparse suite with `Par::Seq` (regression), then with `Par::rayon(4)` (correctness)

### Step 3: Tree-level parallel factorization
- Refactor sequential postorder loop into recursive function
- Add `rayon::scope` at each tree fork (when `par` is `Par::Rayon`)
- Handle per-task `global_to_local` allocation
- Pre-allocate output vectors as `Vec<Option<_>>` for index-based writes
- Compute aggregate stats after tree traversal completes
- **Tests**: Correctness (backward error < 5e-11), determinism (run twice, compare bitwise)

### Step 4: Thread `par` through solve
- Add `par` parameter to `aptp_solve` and per-supernode solve functions
- Replace `Par::Seq` at all 4 sites in solve.rs with `par` (threshold-gated)
- Implement diagonal solve with `rayon::par_iter` (or par_chunks)
- Implement forward/backward solve with recursive `rayon::scope` matching factorization tree
- Per-task workspace allocation for parallel solve
- **Tests**: Compare parallel solve output against sequential (backward error equivalence)

### Step 5: Benchmark tool and scaling report
- Extend `baseline_collection.rs` or create new `parallel_scaling.rs` example
- Run factorization and solve at 1, 2, 4, 8 threads
- Compute speedup, parallel efficiency per matrix
- Classify by workload type (IntraNode/Mixed/TreeLevel)
- **Validation**: SC-001 through SC-008

## Key Files to Modify

| File | Changes |
|------|---------|
| `Cargo.toml` | Add `rayon = "1"` |
| `src/aptp/solver.rs` | Add `par` to FactorOptions, SolverOptions; propagate through SparseLDLT |
| `src/aptp/factor.rs` | Add `par` to AptpOptions; replace Par::Seq in TRSM/GEMM |
| `src/aptp/numeric.rs` | Refactor postorder loop to recursive + rayon::scope |
| `src/aptp/solve.rs` | Add `par` parameter; parallel diagonal/forward/backward |
| `examples/parallel_scaling.rs` | New: benchmark tool for scaling report |

## Smoke Test

After Step 2 (intra-node parallelism):
```bash
# Quick correctness check
cargo test -- --ignored test_sparsine --test-threads=1

# Timing comparison
cargo run --example solve_timing --release -- --matrix sparsine
```

After Step 3 (tree-level parallelism):
```bash
# Full suite correctness
cargo test -- --ignored --test-threads=1

# Parallel benchmark
cargo run --example parallel_scaling --features diagnostic --release -- --ci-only
```
