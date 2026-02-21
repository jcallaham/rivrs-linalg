# API Contract: SparseLDLT Public API

**Feature**: 016-triangular-solve-api
**Module**: `rivrs_sparse::aptp::solver`

## SparseLDLT

```rust
/// High-level sparse symmetric indefinite solver.
///
/// Wraps the APTP multifrontal factorization (Phases 2-6) with a
/// three-phase API: analyze → factor → solve. The symbolic analysis
/// is reusable across factorizations with the same sparsity pattern.
///
/// # SPRAL Equivalent
/// Corresponds to the `ssids_akeep` (symbolic) + `ssids_fkeep` (numeric)
/// data structures and the `ssids_analyse`, `ssids_factor`, `ssids_solve`
/// subroutines.
pub struct SparseLDLT { /* fields private */ }

impl SparseLDLT {
    /// Symbolic analysis phase.
    ///
    /// Computes the fill-reducing ordering, elimination tree, and
    /// supernode structure. Optionally computes MC64 scaling factors
    /// when `MatchOrderMetis` ordering is requested.
    ///
    /// # Errors
    /// - `SparseError::AnalysisFailure` if ordering or symbolic analysis fails
    pub fn analyze(
        matrix: SymbolicSparseColMatRef<'_, usize>,
        options: &AnalyzeOptions,
    ) -> Result<Self, SparseError>;

    /// Numeric factorization phase.
    ///
    /// Computes P^T A P = L D L^T using multifrontal APTP.
    /// If scaling factors are present (from MC64 ordering), entries
    /// are scaled at assembly time: scaled_value = s[i] * a[i][j] * s[j].
    ///
    /// # Errors
    /// - `SparseError::DimensionMismatch` if matrix dimensions differ from analysis
    pub fn factor(
        &mut self,
        matrix: &SparseColMat<usize, f64>,
        options: &FactorOptions,
    ) -> Result<(), SparseError>;

    /// Refactor with same sparsity pattern.
    ///
    /// Equivalent to calling `factor()` again — provided for API clarity
    /// and to mirror SPRAL's `ssids_factor` with existing `akeep`.
    pub fn refactor(
        &mut self,
        matrix: &SparseColMat<usize, f64>,
        options: &FactorOptions,
    ) -> Result<(), SparseError>;

    /// Solve Ax = b (allocating).
    ///
    /// Returns the solution vector. Applies permutation and optional
    /// scaling around the core triangular solve.
    ///
    /// # Errors
    /// - `SparseError::SolveBeforeFactor` if `factor()` has not been called
    /// - `SparseError::DimensionMismatch` if rhs length != matrix dimension
    pub fn solve(
        &self,
        rhs: &Col<f64>,
        stack: &mut MemStack,
    ) -> Result<Col<f64>, SparseError>;

    /// Solve Ax = b (in-place).
    ///
    /// Modifies `rhs` to contain the solution. Applies permutation
    /// and optional scaling around the core triangular solve.
    ///
    /// # Errors
    /// - `SparseError::SolveBeforeFactor` if `factor()` has not been called
    /// - `SparseError::DimensionMismatch` if rhs length != matrix dimension
    pub fn solve_in_place(
        &self,
        rhs: &mut Col<f64>,
        stack: &mut MemStack,
    ) -> Result<(), SparseError>;

    /// Workspace requirement for solve.
    ///
    /// Returns the `StackReq` needed by `solve_in_place()` for the
    /// given number of RHS columns (Phase 7: always 1).
    pub fn solve_scratch(&self, rhs_ncols: usize) -> StackReq;

    /// One-shot: analyze + factor + solve.
    ///
    /// Convenience method that allocates workspace internally.
    pub fn solve_full(
        matrix: &SparseColMat<usize, f64>,
        rhs: &Col<f64>,
        options: &SolverOptions,
    ) -> Result<Col<f64>, SparseError>;

    /// Get inertia from the most recent factorization.
    ///
    /// Returns `None` if `factor()` has not been called.
    pub fn inertia(&self) -> Option<Inertia>;

    /// Get factorization statistics.
    ///
    /// Returns `None` if `factor()` has not been called.
    pub fn stats(&self) -> Option<&FactorizationStats>;

    /// Matrix dimension (from symbolic analysis).
    pub fn n(&self) -> usize;
}
```

## Option Types

```rust
/// Configuration for the symbolic analysis phase.
#[derive(Debug, Clone)]
pub struct AnalyzeOptions {
    /// Fill-reducing ordering strategy.
    pub ordering: OrderingStrategy,
}

/// Fill-reducing ordering strategy.
#[derive(Debug, Clone)]
pub enum OrderingStrategy {
    /// AMD ordering (faer built-in).
    Amd,
    /// METIS ordering (via metis-sys). Default.
    Metis,
    /// MC64 matching + METIS ordering. Produces scaling factors.
    MatchOrderMetis,
    /// User-supplied ordering permutation.
    UserSupplied(Perm<usize>),
}

/// Configuration for the numeric factorization phase.
#[derive(Debug, Clone)]
pub struct FactorOptions {
    /// APTP pivot threshold (u parameter). Default: 0.01.
    pub threshold: f64,
    /// Fallback strategy for failed 1x1 pivots. Default: BunchKaufman.
    pub fallback: AptpFallback,
}

/// Configuration for the one-shot solve.
#[derive(Debug, Clone)]
pub struct SolverOptions {
    /// Fill-reducing ordering strategy.
    pub ordering: OrderingStrategy,
    /// APTP pivot threshold. Default: 0.01.
    pub threshold: f64,
    /// Fallback strategy. Default: BunchKaufman.
    pub fallback: AptpFallback,
}
```

## Error Additions

```rust
/// New variant added to SparseError:
SolveBeforeFactor {
    context: String,
},
```

## Re-exports from `aptp/mod.rs`

```rust
pub use solver::{SparseLDLT, AnalyzeOptions, FactorOptions, SolverOptions, OrderingStrategy};
```
