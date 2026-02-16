# API Contract: Internal Solve Functions

**Feature**: 016-triangular-solve-api
**Module**: `rivrs_sparse::aptp::solve`

## Core Solve Function

```rust
/// Per-supernode triangular solve through the multifrontal factor structure.
///
/// Solves L D L^T x = b in-place, where b is provided in `rhs` (permuted
/// coordinate system — the caller handles P and S). Traverses the assembly
/// tree performing per-supernode gather/solve/scatter operations.
///
/// # Algorithm
///
/// 1. Forward solve (postorder): L y = b
/// 2. Diagonal solve: D z = y
/// 3. Backward solve (reverse postorder): L^T x = z
///
/// # Arguments
///
/// - `symbolic`: Symbolic analysis (supernode structure, traversal order)
/// - `numeric`: Numeric factors (per-supernode L11, D11, L21)
/// - `rhs`: Right-hand side, modified in-place to contain solution
/// - `stack`: Pre-allocated workspace (from `aptp_solve_scratch()`)
///
/// # Errors
///
/// Returns `SparseError::DimensionMismatch` if `rhs.len() != symbolic.nrows()`.
///
/// # References
///
/// - Duff, Hogg & Lopez (2020), Section 3: APTP solve algorithm
/// - Liu (1992), Sections 4-5: multifrontal solve traversal
/// - faer `SupernodalIntranodeLbltRef::solve_in_place_no_numeric_permute_with_conj`
///   (reference for per-supernode gather/scatter pattern)
pub(crate) fn aptp_solve(
    symbolic: &AptpSymbolic,
    numeric: &AptpNumeric,
    rhs: &mut [f64],
    stack: &mut MemStack,
) -> Result<(), SparseError>;

/// Workspace requirement for `aptp_solve`.
///
/// Returns `StackReq` for the temporary per-supernode buffer.
/// Size is determined by `max_front_size` from factorization stats.
pub(crate) fn aptp_solve_scratch(
    numeric: &AptpNumeric,
    rhs_ncols: usize,
) -> StackReq;
```

## Per-Supernode Helper Functions

```rust
/// Forward solve for a single supernode.
///
/// 1. Gather rhs[col_indices[i]] → local buffer
/// 2. Solve L11 * y = local (unit lower triangular, in-place)
/// 3. Write local → rhs[col_indices[i]]
/// 4. Scatter: rhs[row_indices[j]] -= (L21 * local)[j]
///
/// # Arguments
///
/// - `ff`: Per-supernode factors (L11, D11, L21, col_indices, row_indices)
/// - `rhs`: Global RHS vector (permuted coordinates), modified in-place
/// - `work`: Workspace buffer of length >= ne + r (sliced from MemStack)
pub(crate) fn forward_solve_supernode(
    ff: &FrontFactors,
    rhs: &mut [f64],
    work: &mut [f64],
);

/// Diagonal solve for a single supernode.
///
/// 1. Gather rhs[col_indices[i]] → local buffer
/// 2. D11.solve_in_place(local) (handles 1x1, 2x2, zero pivots)
/// 3. Write local → rhs[col_indices[i]]
pub(crate) fn diagonal_solve_supernode(
    ff: &FrontFactors,
    rhs: &mut [f64],
    work: &mut [f64],
);

/// Backward solve for a single supernode.
///
/// 1. Gather rhs[col_indices[i]] → local buffer
/// 2. Gather rhs[row_indices[j]] → tmp buffer (if r > 0)
/// 3. local -= L21^T * tmp (if r > 0)
/// 4. Solve L11^T * z = local (unit upper triangular via transpose)
/// 5. Write local → rhs[col_indices[i]]
pub(crate) fn backward_solve_supernode(
    ff: &FrontFactors,
    rhs: &mut [f64],
    work: &mut [f64],
);
```

## Validation Addition

```rust
/// Module: rivrs_sparse::validate

/// Compute backward error using sparse matrix-vector multiply.
///
/// Returns ||Ax - b|| / (||A||_F * ||x|| + ||b||).
/// Uses `sparse_dense_matmul` to avoid O(n^2) dense conversion.
///
/// **Symmetric handling**: Since A is stored as lower triangle in CSC,
/// the full symmetric product A*x is computed via two passes:
///   A*x = A_lower*x + A_lower^T*x - diag(A)*x
/// This is O(nnz) work and O(n) extra memory (diagonal extraction).
///
/// For ||A||_F, iterates over sparse values directly with
/// 2x multiplier for off-diagonal entries (stored once in lower triangle).
pub fn sparse_backward_error(
    a: &SparseColMat<usize, f64>,
    x: &Col<f64>,
    b: &Col<f64>,
) -> f64;
```

## Modified Existing Functions

```rust
/// Module: rivrs_sparse::aptp::numeric

/// MODIFIED: Add optional scaling parameter.
///
/// When `scaling` is Some, matrix entries are scaled during scatter:
///   scaled_val = scaling[perm_inv[orig_row]] * val * scaling[perm_inv[orig_col]]
///
/// The APTP kernel, extend-add, and all other logic remain unchanged.
pub fn factor(
    symbolic: &AptpSymbolic,
    matrix: &SparseColMat<usize, f64>,
    options: &AptpOptions,
    scaling: Option<&[f64]>,  // NEW parameter
) -> Result<AptpNumeric, SparseError>;
```

```rust
/// Module: rivrs_sparse::aptp::diagonal

/// MODIFIED: Handle zero pivots gracefully.
///
/// For 1x1 pivots with d == 0.0: sets x[col] = 0.0.
/// For 2x2 blocks with det == 0.0: sets both components to 0.0.
/// Previously panicked on these conditions via debug_assert.
pub fn solve_in_place(&self, x: &mut [f64]);
```
