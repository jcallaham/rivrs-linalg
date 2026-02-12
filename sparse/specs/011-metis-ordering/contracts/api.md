# API Contract: METIS Ordering

**Feature**: 011-metis-ordering
**Date**: 2026-02-12

## Public API Surface

### Function: `metis_ordering`

**Module**: `crate::aptp::ordering` (new module, re-exported from `crate::aptp`)

**Signature** (notional):
```rust
/// Compute a METIS nested dissection fill-reducing ordering for a symmetric sparse matrix.
///
/// Extracts the adjacency graph from the matrix's sparsity pattern and calls METIS
/// to compute a fill-reducing permutation via multilevel nested dissection.
///
/// The returned permutation is suitable for use with
/// `AptpSymbolic::analyze(matrix, SymmetricOrdering::Custom(perm.as_ref()))`.
///
/// # Arguments
///
/// * `matrix` — Symbolic sparsity pattern of a symmetric sparse matrix.
///   Only the structural pattern is used; numerical values are ignored.
///
/// # Returns
///
/// A fill-reducing permutation as `Perm<usize>`.
///
/// # Errors
///
/// Returns `SparseError` if:
/// - Matrix dimension exceeds `i32::MAX` (METIS limitation)
/// - METIS reports an internal error
///
/// # Algorithm References
///
/// - Karypis & Kumar (1998), "A Fast and High Quality Multilevel Scheme for
///   Partitioning Irregular Graphs", SIAM J. Sci. Comput.
/// - George (1973), "Nested Dissection of a Regular Finite Element Mesh"
pub fn metis_ordering(
    matrix: SymbolicSparseColMatRef<'_, usize>,
) -> Result<Perm<usize>, SparseError>
```

**Input contract**:
- `matrix` is the symbolic structure of a symmetric sparse matrix
- Matrix may store upper triangle, lower triangle, or full structure
- Dimension must fit in `i32` (< 2,147,483,647)

**Output contract**:
- Returns a valid `Perm<usize>` of dimension `n` (= number of rows/columns)
- Forward and inverse arrays are both valid permutations of `{0, ..., n-1}`

**Edge case behavior**:
- Dimension 0: returns empty permutation (identity)
- Dimension 1: returns trivial 1-element permutation
- Diagonal matrix (no off-diagonal entries): returns identity permutation (METIS has no edges to partition)
- Disconnected graph: METIS handles internally; valid permutation returned

### Re-export from `crate::aptp`

```rust
// In src/aptp/mod.rs, add:
pub mod ordering;
pub use ordering::metis_ordering;
```

## Internal Functions (not public API)

### `extract_adjacency`

Extracts the undirected graph adjacency structure from a symmetric CSC matrix.

```rust
/// Extract CSR adjacency arrays (xadj, adjncy) from a symmetric sparse matrix.
///
/// Excludes diagonal entries. Ensures symmetry (both (i,j) and (j,i) present).
fn extract_adjacency(
    matrix: SymbolicSparseColMatRef<'_, usize>,
) -> (Vec<i32>, Vec<i32>)
```

### Integration with existing API

No changes to `AptpSymbolic::analyze()`. Usage pattern:

```rust
use rivrs_sparse::aptp::{AptpSymbolic, metis_ordering};
use faer::sparse::linalg::cholesky::SymmetricOrdering;

let perm = metis_ordering(matrix.symbolic())?;
let symbolic = AptpSymbolic::analyze(
    matrix.symbolic(),
    SymmetricOrdering::Custom(perm.as_ref()),
)?;
```
