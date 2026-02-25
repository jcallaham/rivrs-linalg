# Test API Contracts: Robustness — Testing & Hardening

**Feature**: 026-robustness-hardening
**Date**: 2026-02-25

## Perturbation Helper API

### `cause_delays`

```
Input:  mutable dense symmetric matrix (n×n), RNG
Output: matrix modified in-place
Effect: n/8 random rows multiplied by 1000, n/8 random entries multiplied by 1000
        First oversized row index > block_size when n > block_size
Post:   Matrix remains symmetric
```

### `make_singular`

```
Input:  mutable dense symmetric matrix (n×n), col1 index, col2 index
Output: matrix modified in-place
Effect: col2 becomes scaled copy of col1; corresponding rows updated for symmetry
Post:   Matrix is rank-deficient (rank decreased by 1), remains symmetric
```

### `make_dblk_singular`

```
Input:  mutable dense symmetric matrix (n×n), block starting row, block size
Output: matrix modified in-place
Effect: First and last columns of specified diagonal block made linearly dependent
Post:   Specified diagonal block is singular, matrix remains symmetric
```

## Proptest Strategy API

### `arb_symmetric_pd`

```
Input:  size range (min..=max)
Output: proptest Strategy yielding Mat<f64>
Post:   Generated matrix is symmetric positive definite (diagonal dominance)
Shrinks: toward smaller sizes
```

### `arb_symmetric_indefinite`

```
Input:  size range (min..=max)
Output: proptest Strategy yielding Mat<f64>
Post:   Generated matrix is symmetric indefinite (mixed diagonal signs)
Shrinks: toward smaller sizes
```

### `arb_sparse_symmetric`

```
Input:  size range (min..=max), density range (min..=max)
Output: proptest Strategy yielding SparseColMat<usize, f64>
Post:   Generated matrix is sparse, symmetric (full CSC storage), nonzero
Shrinks: toward smaller sizes and lower density
```

## Torture Test Entry Points

### `ldlt_app_torture_test`

```
Input:  AptpOptions, matrix size (m, n), num_instances, seed
Output: test pass/fail
Calls:  aptp_factor_in_place with complete pivoting path
Assert: zero panics, backward error < 5e-11 for non-singular instances,
        correct inertia for singular instances (n_zero > 0)
```

### `ldlt_tpp_torture_test`

```
Input:  AptpOptions, matrix size (m, n), num_instances, seed
Output: test pass/fail
Calls:  aptp_factor_in_place with TPP path (num_fully_summed < inner_block_size)
Assert: zero panics, backward error < 5e-11 for non-singular instances,
        L growth bound ≤ 1/threshold
```

## Property Test Contracts

### Property: Reconstruction Valid

```
For all symmetric matrices M (size 5–500):
  let result = aptp_factor(M)
  if result.is_ok():
    assert reconstruction_error(M, result) < 1e-12
```

### Property: Inertia Consistent

```
For all symmetric matrices M (size 5–500):
  let result = aptp_factor(M)
  if result.is_ok():
    let inertia = result.inertia()
    assert inertia.n_positive + inertia.n_negative + inertia.n_zero == M.nrows()
```

### Property: Permutation Valid

```
For all symmetric matrices M (size 5–500):
  let result = aptp_factor(M)
  if result.is_ok():
    let perm = result.permutation()
    assert perm is a valid permutation of 0..n
    assert perm.fwd() and perm.inv() are inverses
```

### Property: No Panics on Adversarial

```
For all adversarial matrices M (from perturbation helpers):
  aptp_factor_in_place(M) does not panic
  (may return Ok or Err)
```

## Adversarial Input Contracts

### Edge Case Tests

```
For each edge case input:
  SparseLDLT::analyze_with_matrix(input) either:
    - returns Ok (valid analysis) followed by factor/solve success, OR
    - returns Err(SparseError::*) with appropriate variant
  NEVER panics
```

| Input | Expected Error Variant (if error) |
|-------|----------------------------------|
| 0×0 matrix | InvalidInput or DimensionMismatch |
| 1×1 zero diagonal | NumericalSingularity or successful zero-inertia |
| NaN entries | InvalidInput |
| Inf entries | InvalidInput |
| Non-symmetric pattern | InvalidInput |
