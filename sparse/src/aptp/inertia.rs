//! Eigenvalue sign classification (inertia) for symmetric matrices.
//!
//! Provides [`Inertia`], which records the count of positive, negative, and zero
//! eigenvalues of a symmetric matrix. Used both for reference factorization
//! validation and for live inertia computation from [`super::MixedDiagonal`].
//!
//! Reference: Hogg, Duff & Lopez (2020), "A New Sparse LDL^T Solver Using
//! A Posteriori Threshold Pivoting", SIAM J. Sci. Comput. 42(4), Section 2.

use serde::Deserialize;

/// Eigenvalue sign classification of a symmetric matrix.
///
/// Records the number of positive, negative, and zero eigenvalues. This is a
/// fundamental property of symmetric matrices that is preserved under congruence
/// transformations (Sylvester's Law of Inertia) and can be determined from the
/// D factor in an LDL^T factorization without computing actual eigenvalues.
///
/// # References
///
/// - Hogg, Duff & Lopez (2020), Section 2: inertia as a correctness check
///   for indefinite factorizations
/// - Bunch & Kaufman (1977): inertia computation from pivot classifications
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, Deserialize)]
pub struct Inertia {
    /// Count of positive eigenvalues.
    pub positive: usize,
    /// Count of negative eigenvalues.
    pub negative: usize,
    /// Count of zero eigenvalues.
    pub zero: usize,
}

impl Inertia {
    /// Total matrix dimension (positive + negative + zero).
    pub fn dimension(&self) -> usize {
        self.positive + self.negative + self.zero
    }
}
