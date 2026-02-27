# Module Contracts: 004-repo-setup

**Date**: 2026-02-06

This feature is a Rust library (not a web service), so contracts are expressed as module-level public API signatures.

## Module: `src/io/mtx.rs` (Matrix Market Parser)

```rust
/// Parse a Matrix Market file into a faer sparse symmetric matrix.
///
/// Supports `coordinate real symmetric` format only. Entries are
/// 1-indexed in the file and converted to 0-indexed internally.
/// The returned matrix stores the full symmetric matrix (both triangles).
///
/// # Errors
/// - File not found or unreadable
/// - Unsupported Matrix Market format (not coordinate real symmetric)
/// - Malformed header, size line, or data lines
/// - Row/col indices out of bounds
pub fn load_mtx(path: &Path) -> Result<SparseColMat<usize, f64>, SparseError>;
```

## Module: `src/io/reference.rs` (Reference Factorization Loader)

```rust
/// Load a reference factorization from a companion JSON file.
///
/// # Errors
/// - File not found or unreadable
/// - Invalid JSON structure
/// - Inconsistent data (e.g., l_entry indices out of bounds for given permutation size)
pub fn load_reference(path: &Path) -> Result<ReferenceFactorization, SparseError>;
```

## Module: `src/io/registry.rs` (Test Matrix Registry)

```rust
/// Load the test matrix registry from metadata.json.
///
/// # Errors
/// - metadata.json not found at expected location
/// - Invalid JSON structure
pub fn load_registry() -> Result<Vec<MatrixMetadata>, SparseError>;

/// Load a specific test matrix by name, including its sparse matrix
/// and optional reference factorization.
///
/// Returns None if the matrix's .mtx file does not exist on disk
/// (e.g., gitignored SuiteSparse matrix not extracted).
///
/// # Errors
/// - Matrix name not found in registry
/// - .mtx file exists but fails to parse
/// - .json file exists but fails to parse
pub fn load_test_matrix(name: &str) -> Result<Option<TestMatrix>, SparseError>;
```

## Module: `src/validate.rs` (Numerical Validation)

```rust
/// Compute the reconstruction error: ||P^T A P - L D L^T||_F / ||A||_F
///
/// Converts all inputs to dense for computation. Appropriate for
/// small matrices (hand-constructed, up to ~20×20). For large matrices,
/// a sparse-aware validation will be needed in later phases.
///
/// # Arguments
/// - `a`: The original sparse symmetric matrix
/// - `reference`: The reference factorization (L, D, P)
///
/// # Returns
/// The relative Frobenius-norm reconstruction error (dimensionless scalar).
pub fn reconstruction_error(
    a: SparseColMatRef<usize, f64>,
    reference: &ReferenceFactorization,
) -> f64;

/// Compute the scaled backward error: ||Ax - b|| / (||A|| ||x|| + ||b||)
///
/// Uses sparse matrix-vector multiply for A*x.
///
/// # Arguments
/// - `a`: The coefficient matrix (sparse)
/// - `x`: The computed solution vector
/// - `b`: The right-hand side vector
///
/// # Returns
/// The scaled backward error (dimensionless scalar).
pub fn backward_error(
    a: SparseColMatRef<usize, f64>,
    x: &[f64],
    b: &[f64],
) -> f64;

/// Check that computed inertia exactly matches expected inertia.
pub fn check_inertia(computed: &Inertia, expected: &Inertia) -> bool;
```

## Module: `src/io.rs` (Re-exports)

```rust
pub mod mtx;
pub mod reference;
pub mod registry;
```

## Test Infrastructure

### `tests/hand_constructed.rs`

Integration test that:
1. Loads all 15 hand-constructed matrices via registry
2. Verifies dimensions and nnz match metadata
3. For each matrix with a reference factorization, computes reconstruction error
4. Asserts reconstruction error < 10^-12

### `tests/suitesparse_ci.rs`

Integration test that:
1. Loads all 10 CI-subset matrices via registry
2. Verifies dimensions and nnz match metadata
3. Verifies matrix is loadable and has correct structure

### `benches/matrix_loading.rs`

Criterion benchmark that:
1. Benchmarks loading `arrow-10-indef.mtx` from disk
2. Benchmarks loading + property verification
