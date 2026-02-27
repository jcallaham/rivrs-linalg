//! User-facing sparse symmetric indefinite solver API.
//!
//! Provides [`SparseLDLT`], the high-level solver struct wrapping the
//! three-phase analyze → factor → solve pipeline. The symbolic analysis
//! is reusable across factorizations with the same sparsity pattern.
//!
//! # References
//!
//! - Duff, Hogg & Lopez (2020), "A New Sparse LDL^T Solver Using A Posteriori
//!   Threshold Pivoting", SIAM J. Sci. Comput. 42(4)
//! - Liu (1992), "The Multifrontal Method for Sparse Matrix Solution"

use faer::Col;
use faer::Par;
use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::perm::Perm;
use faer::sparse::linalg::cholesky::SymmetricOrdering;
use faer::sparse::{SparseColMat, SymbolicSparseColMatRef};

use super::factor::{AptpFallback, AptpOptions};
use super::inertia::Inertia;
use super::numeric::{AptpNumeric, FactorizationStats};
use super::solve::{aptp_solve, aptp_solve_scratch};
use super::symbolic::AptpSymbolic;
use crate::error::SparseError;
use crate::ordering::{match_order_metis, metis_ordering};

/// Fill-reducing ordering strategy.
///
/// The ordering is a user-configurable option rather than an automatic
/// heuristic. Guidance from Duff, Hogg & Lopez (2020):
///
/// - **Easy indefinite** (structural FEM, thermal, acoustic, model reduction,
///   quantum chemistry): use [`Metis`](Self::Metis) (the default).
/// - **Hard indefinite** (KKT/saddle-point, optimal control, power networks,
///   mixed FEM / Stokes): use [`MatchOrderMetis`](Self::MatchOrderMetis).
#[derive(Debug, Clone)]
pub enum OrderingStrategy {
    /// AMD ordering (faer built-in).
    Amd,
    /// METIS ordering (via metis-sys). Default.
    /// Best for easy-indefinite problems with naturally dominant diagonal.
    Metis,
    /// MC64 matching + METIS ordering. Produces scaling factors.
    /// Best for hard-indefinite problems (KKT, saddle-point, zero diagonal blocks).
    MatchOrderMetis,
    /// User-supplied ordering permutation.
    UserSupplied(Perm<usize>),
}

/// Configuration for the symbolic analysis phase.
#[derive(Debug, Clone)]
pub struct AnalyzeOptions {
    /// Fill-reducing ordering strategy.
    pub ordering: OrderingStrategy,
}

impl Default for AnalyzeOptions {
    fn default() -> Self {
        Self {
            ordering: OrderingStrategy::MatchOrderMetis,
        }
    }
}

/// Configuration for the numeric factorization phase.
#[derive(Debug, Clone)]
pub struct FactorOptions {
    /// APTP pivot threshold (u parameter). Default: 0.01.
    pub threshold: f64,
    /// Fallback strategy for failed 1x1 pivots. Default: BunchKaufman.
    pub fallback: AptpFallback,
    /// Outer block size (nb) for two-level APTP. Default: 256.
    pub outer_block_size: usize,
    /// Inner block size (ib) for two-level APTP. Default: 32.
    pub inner_block_size: usize,
    /// Parallelism control for dense BLAS operations within supernodes
    /// and tree-level scheduling. Default: `Par::Seq`.
    pub par: Par,
    /// Minimum supernode size for amalgamation. Default: 32.
    /// Setting `nemin = 1` disables amalgamation.
    pub nemin: usize,
    /// Front-size threshold for small-leaf subtree fast path. Supernodes with
    /// front size strictly less than this value may be grouped into small-leaf
    /// subtrees and processed via a streamlined code path that reuses a single
    /// small workspace, avoiding per-supernode dispatch overhead.
    ///
    /// Default: 256. Set to 0 to disable the fast path entirely.
    ///
    /// # Relationship to `INTRA_NODE_THRESHOLD`
    ///
    /// This threshold determines which subtrees use the fast path (workspace
    /// reuse, sequential processing). `INTRA_NODE_THRESHOLD` (also 256) controls
    /// whether intra-node BLAS uses parallel or sequential execution. They are
    /// independent but share the same default because fronts below 256 benefit
    /// from neither parallel BLAS nor the general-path overhead.
    pub small_leaf_threshold: usize,
}

impl Default for FactorOptions {
    fn default() -> Self {
        Self {
            threshold: 0.01,
            fallback: AptpFallback::BunchKaufman,
            outer_block_size: 256,
            inner_block_size: 32,
            par: Par::Seq,
            nemin: 32,
            small_leaf_threshold: 256,
        }
    }
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
    /// Parallelism control for factorization and solve. Default: `Par::Seq`.
    pub par: Par,
    /// Minimum supernode size for amalgamation. Default: 32.
    /// Setting `nemin = 1` disables amalgamation.
    pub nemin: usize,
    /// Front-size threshold for small-leaf subtree fast path.
    /// Default: 256. Set to 0 to disable. See [`FactorOptions::small_leaf_threshold`].
    pub small_leaf_threshold: usize,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            ordering: OrderingStrategy::MatchOrderMetis,
            threshold: 0.01,
            fallback: AptpFallback::BunchKaufman,
            par: Par::Seq,
            nemin: 32,
            small_leaf_threshold: 256,
        }
    }
}

/// High-level sparse symmetric indefinite solver.
///
/// Wraps the APTP multifrontal factorization with a
/// three-stage API: analyze → factor → solve. The symbolic analysis
/// is reusable across factorizations with the same sparsity pattern.
///
/// # Examples
///
/// ```
/// use faer::sparse::{SparseColMat, Triplet};
/// use faer::Col;
/// use rivrs_sparse::symmetric::{SparseLDLT, SolverOptions};
///
/// let triplets = vec![
///     Triplet::new(0, 0, 4.0),
///     Triplet::new(0, 1, 1.0),
///     Triplet::new(1, 0, 1.0),
///     Triplet::new(1, 1, 3.0),
/// ];
/// let a = SparseColMat::try_new_from_triplets(2, 2, &triplets).unwrap();
/// let b = Col::from_fn(2, |i| [5.0, 4.0][i]);
/// let x = SparseLDLT::solve_full(&a, &b, &SolverOptions::default()).unwrap();
/// // A = [[4, 1], [1, 3]], b = [5, 4] => x = [1, 1]
/// assert!((x[0] - 1.0).abs() < 1e-12);
/// assert!((x[1] - 1.0).abs() < 1e-12);
/// ```
pub struct SparseLDLT {
    symbolic: AptpSymbolic,
    numeric: Option<AptpNumeric>,
    scaling: Option<Vec<f64>>,
    /// Cached forward permutation (perm_fwd[new] = old). Computed once during analyze.
    perm_fwd: Vec<usize>,
}

impl SparseLDLT {
    /// Build a `SparseLDLT` with cached permutation vectors from symbolic analysis.
    fn new_with_cached_perm(symbolic: AptpSymbolic, scaling: Option<Vec<f64>>) -> Self {
        let (perm_fwd, _) = symbolic.perm_vecs();
        SparseLDLT {
            symbolic,
            numeric: None,
            scaling,
            perm_fwd,
        }
    }

    /// Symbolic analysis phase.
    ///
    /// Computes the fill-reducing ordering, elimination tree, and
    /// supernode structure. Optionally computes MC64 scaling factors
    /// when `MatchOrderMetis` ordering is requested.
    ///
    /// Note: `MatchOrderMetis` ordering requires numeric matrix values;
    /// use [`analyze_with_matrix`](Self::analyze_with_matrix) instead.
    ///
    /// # Errors
    ///
    /// - `SparseError::AnalysisFailure` if ordering or symbolic analysis fails
    pub fn analyze(
        matrix: SymbolicSparseColMatRef<'_, usize>,
        options: &AnalyzeOptions,
    ) -> Result<Self, SparseError> {
        let n = matrix.nrows();

        match &options.ordering {
            OrderingStrategy::Amd => {
                let symbolic = AptpSymbolic::analyze(matrix, SymmetricOrdering::Amd)?;
                Ok(Self::new_with_cached_perm(symbolic, None))
            }
            OrderingStrategy::Metis => {
                // Build a dummy matrix for METIS (needs SymbolicSparseColMatRef)
                let col_ptrs = matrix.col_ptr();
                let row_indices = matrix.row_idx();
                let mut triplets = Vec::new();
                for j in 0..n {
                    for &i in &row_indices[col_ptrs[j]..col_ptrs[j + 1]] {
                        triplets.push(faer::sparse::Triplet::new(i, j, 1.0f64));
                    }
                }
                let dummy_matrix =
                    SparseColMat::try_new_from_triplets(n, n, &triplets).map_err(|e| {
                        SparseError::AnalysisFailure {
                            reason: format!("Failed to construct matrix for METIS ordering: {}", e),
                        }
                    })?;

                let perm = metis_ordering(dummy_matrix.symbolic())?;
                let symbolic =
                    AptpSymbolic::analyze(matrix, SymmetricOrdering::Custom(perm.as_ref()))?;
                Ok(Self::new_with_cached_perm(symbolic, None))
            }
            OrderingStrategy::MatchOrderMetis => {
                Err(SparseError::AnalysisFailure {
                    reason: "MatchOrderMetis requires numeric matrix values; use analyze_with_matrix() instead".to_string(),
                })
            }
            OrderingStrategy::UserSupplied(perm) => {
                let symbolic =
                    AptpSymbolic::analyze(matrix, SymmetricOrdering::Custom(perm.as_ref()))?;
                Ok(Self::new_with_cached_perm(symbolic, None))
            }
        }
    }

    /// Symbolic analysis with full matrix (required for MatchOrderMetis).
    ///
    /// When `MatchOrderMetis` is requested, this method performs MC64
    /// matching on the numeric matrix to compute scaling factors and
    /// a fill-reducing ordering. Other orderings delegate to [`analyze`](Self::analyze).
    pub fn analyze_with_matrix(
        matrix: &SparseColMat<usize, f64>,
        options: &AnalyzeOptions,
    ) -> Result<Self, SparseError> {
        let n = matrix.nrows();

        match &options.ordering {
            OrderingStrategy::MatchOrderMetis => {
                let result = match_order_metis(matrix)?;
                let ordering_perm = result.ordering;
                let symbolic = AptpSymbolic::analyze(
                    matrix.symbolic(),
                    SymmetricOrdering::Custom(ordering_perm.as_ref()),
                )?;

                // Transform scaling to elimination order
                let (perm_fwd, _) = symbolic.perm_vecs();
                let elim_scaling: Vec<f64> = (0..n).map(|i| result.scaling[perm_fwd[i]]).collect();

                Ok(Self::new_with_cached_perm(symbolic, Some(elim_scaling)))
            }
            // All other orderings don't need numeric values
            _ => Self::analyze(matrix.symbolic(), options),
        }
    }

    /// Numeric factorization phase.
    ///
    /// Computes P^T A P = L D L^T using multifrontal APTP.
    /// If scaling factors are present (from MC64 ordering), entries
    /// are scaled at assembly time.
    ///
    /// # Errors
    ///
    /// - `SparseError::DimensionMismatch` if matrix dimensions differ from analysis
    pub fn factor(
        &mut self,
        matrix: &SparseColMat<usize, f64>,
        options: &FactorOptions,
    ) -> Result<(), SparseError> {
        let aptp_options = AptpOptions {
            threshold: options.threshold,
            fallback: options.fallback,
            outer_block_size: options.outer_block_size,
            inner_block_size: options.inner_block_size,
            par: options.par,
            nemin: options.nemin,
            small_leaf_threshold: options.small_leaf_threshold,
            ..AptpOptions::default()
        };

        // Drop old numeric before building new one to avoid holding both
        // in memory simultaneously. For H2O (max_front=9258), old+new
        // front_factors can exceed 7 GB.
        self.numeric = None;
        let numeric = AptpNumeric::factor(
            &self.symbolic,
            matrix,
            &aptp_options,
            self.scaling.as_deref(),
        )?;
        self.numeric = Some(numeric);
        Ok(())
    }

    /// Refactor with same sparsity pattern but different numeric values.
    ///
    /// Reuses the symbolic analysis from [`analyze_with_matrix`](Self::analyze_with_matrix),
    /// avoiding the cost of re-ordering and re-computing the elimination tree.
    /// This is the inner-loop operation for applications where the sparsity
    /// pattern is fixed but values change — interior-point iterations,
    /// time-stepping schemes, and parameter sweeps.
    ///
    /// Equivalent to calling [`factor`](Self::factor) again.
    pub fn refactor(
        &mut self,
        matrix: &SparseColMat<usize, f64>,
        options: &FactorOptions,
    ) -> Result<(), SparseError> {
        self.factor(matrix, options)
    }

    /// Solve Ax = b (allocating).
    ///
    /// Returns the solution vector.
    ///
    /// # Errors
    ///
    /// - `SparseError::SolveBeforeFactor` if `factor()` has not been called
    /// - `SparseError::DimensionMismatch` if rhs length != matrix dimension
    pub fn solve(
        &self,
        rhs: &Col<f64>,
        stack: &mut MemStack,
        par: Par,
    ) -> Result<Col<f64>, SparseError> {
        let mut result = rhs.to_owned();
        self.solve_in_place(&mut result, stack, par)?;
        Ok(result)
    }

    /// Solve Ax = b (in-place).
    ///
    /// Modifies `rhs` to contain the solution.
    ///
    /// # Errors
    ///
    /// - `SparseError::SolveBeforeFactor` if `factor()` has not been called
    /// - `SparseError::DimensionMismatch` if rhs length != matrix dimension
    pub fn solve_in_place(
        &self,
        rhs: &mut Col<f64>,
        stack: &mut MemStack,
        par: Par,
    ) -> Result<(), SparseError> {
        let numeric = self
            .numeric
            .as_ref()
            .ok_or_else(|| SparseError::SolveBeforeFactor {
                context: "factor() must be called before solve()".to_string(),
            })?;

        let n = self.symbolic.nrows();
        if rhs.nrows() != n {
            return Err(SparseError::DimensionMismatch {
                expected: (n, 1),
                got: (rhs.nrows(), 1),
                context: "RHS length must match matrix dimension".to_string(),
            });
        }

        if n == 0 {
            return Ok(());
        }

        let perm_fwd = &self.perm_fwd;

        // 1. Permute: rhs_perm[new] = rhs[perm_fwd[new]]
        //    perm_fwd[new] = old, so this gathers original RHS into permuted order.
        let mut rhs_perm = vec![0.0f64; n];
        for new in 0..n {
            rhs_perm[new] = rhs[perm_fwd[new]];
        }

        // 2. Scale (if scaling present): rhs_perm[i] *= scaling[i]
        if let Some(ref scaling) = self.scaling {
            for i in 0..n {
                rhs_perm[i] *= scaling[i];
            }
        }

        // 3. Core solve: aptp_solve(symbolic, numeric, &mut rhs_perm, stack, par)
        aptp_solve(&self.symbolic, numeric, &mut rhs_perm, stack, par)?;

        // 4. Unscale: rhs_perm[i] *= scaling[i] (symmetric: same S applied)
        if let Some(ref scaling) = self.scaling {
            for i in 0..n {
                rhs_perm[i] *= scaling[i];
            }
        }

        // 5. Unpermute: rhs[perm_fwd[new]] = rhs_perm[new]
        for new in 0..n {
            rhs[perm_fwd[new]] = rhs_perm[new];
        }

        Ok(())
    }

    /// Workspace requirement for solve.
    ///
    /// Returns the [`StackReq`] needed by [`solve`](Self::solve) or
    /// [`solve_in_place`](Self::solve_in_place). Allocate once with
    /// [`MemBuffer::new`](faer::dyn_stack::MemBuffer::new) and reuse across
    /// multiple solves to avoid repeated allocation.
    ///
    /// `rhs_ncols` is the number of right-hand side columns (1 for a single
    /// vector solve).
    pub fn solve_scratch(&self, rhs_ncols: usize) -> StackReq {
        if let Some(ref numeric) = self.numeric {
            aptp_solve_scratch(numeric, rhs_ncols)
        } else {
            StackReq::empty()
        }
    }

    /// One-shot: analyze + factor + solve.
    ///
    /// Convenience method that runs the full pipeline in a single call.
    /// For repeated solves (multiple RHS or refactorization), use the
    /// three-phase API instead to amortize the symbolic analysis cost.
    ///
    /// # Examples
    ///
    /// ```
    /// use faer::sparse::{SparseColMat, Triplet};
    /// use faer::Col;
    /// use rivrs_sparse::symmetric::{SparseLDLT, SolverOptions};
    ///
    /// let triplets = vec![
    ///     Triplet::new(0, 0, 4.0),
    ///     Triplet::new(0, 1, 1.0), Triplet::new(1, 0, 1.0),
    ///     Triplet::new(1, 1, 3.0),
    /// ];
    /// let a = SparseColMat::try_new_from_triplets(2, 2, &triplets).unwrap();
    /// let b = Col::from_fn(2, |i| [5.0, 4.0][i]);
    /// let x = SparseLDLT::solve_full(&a, &b, &SolverOptions::default()).unwrap();
    /// assert!((x[0] - 1.0).abs() < 1e-12);
    /// ```
    pub fn solve_full(
        matrix: &SparseColMat<usize, f64>,
        rhs: &Col<f64>,
        options: &SolverOptions,
    ) -> Result<Col<f64>, SparseError> {
        let analyze_opts = AnalyzeOptions {
            ordering: options.ordering.clone(),
        };
        let factor_opts = FactorOptions {
            threshold: options.threshold,
            fallback: options.fallback,
            par: options.par,
            nemin: options.nemin,
            small_leaf_threshold: options.small_leaf_threshold,
            ..FactorOptions::default()
        };

        let mut solver = Self::analyze_with_matrix(matrix, &analyze_opts)?;
        solver.factor(matrix, &factor_opts)?;

        let scratch = solver.solve_scratch(1);
        let mut mem = MemBuffer::new(scratch);
        let stack = MemStack::new(&mut mem);
        solver.solve(rhs, stack, options.par)
    }

    /// Get inertia (eigenvalue sign counts) from the most recent factorization.
    ///
    /// Returns the number of positive, negative, and zero eigenvalues of A.
    /// Useful for checking definiteness, verifying KKT optimality conditions
    /// (expected inertia `(n, m, 0)` for `n` primal and `m` dual variables),
    /// and detecting numerical singularity (zero eigenvalues).
    ///
    /// Returns `None` if [`factor`](Self::factor) has not been called.
    pub fn inertia(&self) -> Option<Inertia> {
        self.numeric.as_ref().map(|numeric| {
            let mut inertia = Inertia {
                positive: 0,
                negative: 0,
                zero: 0,
            };
            for ff in numeric.front_factors() {
                let local_inertia = ff.d11().compute_inertia();
                inertia.positive += local_inertia.positive;
                inertia.negative += local_inertia.negative;
                inertia.zero += local_inertia.zero;
            }
            inertia
        })
    }

    /// Get factorization statistics from the most recent factorization.
    ///
    /// Returns aggregate statistics including pivot counts (`total_1x1_pivots`,
    /// `total_2x2_pivots`), delayed pivot count (`total_delayed`), maximum
    /// front size, and supernode count. With the `diagnostic` feature enabled,
    /// also includes per-phase timing breakdowns.
    ///
    /// Returns `None` if [`factor`](Self::factor) has not been called.
    pub fn stats(&self) -> Option<&FactorizationStats> {
        self.numeric.as_ref().map(|n| n.stats())
    }

    /// Get per-supernode diagnostic statistics from the most recent factorization.
    ///
    /// Returns one [`PerSupernodeStats`](super::numeric::PerSupernodeStats) per
    /// supernode with front size, pivot counts, and delay counts. With the
    /// `diagnostic` feature enabled, each entry also includes sub-phase timing
    /// (assembly, kernel, extraction, extend-add).
    ///
    /// Returns `None` if [`factor`](Self::factor) has not been called.
    pub fn per_supernode_stats(&self) -> Option<&[super::numeric::PerSupernodeStats]> {
        self.numeric.as_ref().map(|n| n.per_supernode_stats())
    }

    /// Matrix dimension (from symbolic analysis).
    pub fn n(&self) -> usize {
        self.symbolic.nrows()
    }
}
