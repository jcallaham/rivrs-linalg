# Changelog

## 0.1.0

Initial release of rivrs-sparse: a sparse symmetric indefinite direct solver
(SSIDS) for Rust, implementing multifrontal LDL^T with A Posteriori Threshold
Pivoting (APTP).

### Features

- **Multifrontal LDL^T factorization** with two-level APTP pivoting (TPP for
  small fronts, complete pivoting with BLAS-3 blocking for large fronts)
- **MC64 matching & scaling** + METIS nested dissection ordering for
  hard-indefinite matrices (KKT, saddle-point, optimization)
- **Three-phase API**: analyze (symbolic) -> factor (numeric) -> solve, with
  symbolic reuse across refactorizations
- **Parallel factorization & solve** via rayon (tree-level) and faer `Par`
  (intra-node BLAS)
- **MatrixMarket I/O** with automatic symmetric mirroring
- **Backward error validation** (`sparse_backward_error`)
- **Inertia computation** (eigenvalue sign counts from the D factor)
- **Supernode amalgamation** for reduced dispatch overhead
- **Diagnostic instrumentation** (optional `diagnostic` feature) with Chrome
  Trace export

### Performance

Benchmarked on 65 SuiteSparse matrices:

- Sequential: median 0.94x SPRAL factor time (36% faster, 38% comparable, 24% slower)
- Parallel (8 threads): median 0.89x SPRAL factor time
- vs MUMPS sequential: ~2x faster (median 0.46x)
- All 65 matrices: backward error < 5e-11

### Correctness

- 546 unit tests + 23 ignored (full SuiteSparse suite)
- 65/65 SuiteSparse matrices pass with backward error < 5e-11
- SPRAL-style torture tests (4500 factorizations, zero panics)
- Property-based tests via proptest (7 properties x 256 cases)
- Adversarial edge-case tests (14 tests: 0x0, NaN, Inf, near-overflow)
