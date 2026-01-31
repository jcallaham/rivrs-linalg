//! Condition estimation for Sylvester equations.
//!
//! Provides separation estimation for detecting ill-conditioned or singular
//! Sylvester equations. The separation `sep(A, B)` measures how close
//! the eigenvalues of A are to the negatives of the eigenvalues of B
//! (continuous case) or how close the product of eigenvalues is to -1
//! (discrete case).
//!
//! # Separation Definition
//!
//! For continuous-time (AX + XB = C):
//!   `sep(A, B) = min_{||X||_F = 1} ||AX + XB||_F`
//!
//! For discrete-time (AXB + X = C):
//!   `sep_d(A, B) = min_{||X||_F = 1} ||AXB + X||_F`
//!
//! When sep is zero, A and -B share an eigenvalue and the equation is singular.
//! When sep is small, the equation is ill-conditioned.
//!
//! # Algorithm
//!
//! Uses a cheap lower bound estimator based on the eigenvalues of the
//! Schur forms. For the continuous case, this is:
//!   `sep_est ≥ min_{i,j} |λ_i(A) + λ_j(B)|`
//!
//! For the discrete case:
//!   `sep_est ≥ min_{i,j} |1 + λ_i(A) * λ_j(B)|`
//!
//! This is a lower bound on the true separation, which is sufficient
//! for detecting singular or near-singular cases.
//!
//! # References
//!
//! - Golub & Van Loan (2013), "Matrix Computations" (4th Ed), Section 7.6
//! - Varah (1979), "On the separation of two matrices", SIAM J. Numer. Anal. 16:216-222
//! - Byers (1984), "A LINPACK-style condition estimator for the equation AX-XB=C"

use faer::prelude::*;

use super::triangular::detect_blocks;
use super::types::EquationType;

/// Threshold for declaring a Sylvester equation as singular.
/// If the estimated separation is below this value, the equation
/// is considered too ill-conditioned to solve reliably.
pub const SEPARATION_THRESHOLD: f64 = 1e-14;

/// Threshold for marking an equation as near-singular.
/// If separation is below this but above SEPARATION_THRESHOLD,
/// the solution may have poor accuracy.
pub const NEAR_SINGULAR_THRESHOLD: f64 = 1e-10;

/// Result of separation estimation.
#[derive(Debug, Clone)]
pub struct SeparationEstimate {
    /// Estimated separation value (lower bound).
    pub separation: f64,
    /// Whether the equation appears to be singular (separation ≈ 0).
    pub is_singular: bool,
    /// Whether the equation appears to be near-singular.
    pub is_near_singular: bool,
}

/// Estimates the separation for a Sylvester equation given Schur forms.
///
/// For continuous-time (AX + XB = C), estimates:
///   `sep(T, S) ≥ min_{i,j} |λ_i(T) + λ_j(S)|`
///
/// For discrete-time (AXB + X = C), estimates:
///   `sep_d(T, S) ≥ min_{i,j} |1 + λ_i(T) * λ_j(S)|`
///
/// The Schur forms T and S are quasi-triangular, so eigenvalues are
/// read from 1x1 diagonal blocks (real eigenvalues) and 2x2 diagonal
/// blocks (complex conjugate pairs).
///
/// # Arguments
///
/// - `schur_a`: Matrix A in real Schur form (quasi-triangular)
/// - `schur_b`: Matrix B in real Schur form (quasi-triangular)
/// - `equation_type`: Whether to estimate for continuous or discrete case
///
/// # References
///
/// - Golub & Van Loan (2013), Section 7.6
/// - Varah (1979), "On the separation of two matrices"
pub fn estimate_separation(
    schur_a: MatRef<'_, f64>,
    schur_b: MatRef<'_, f64>,
    equation_type: EquationType,
) -> SeparationEstimate {
    let eig_a = extract_eigenvalues(schur_a);
    let eig_b = extract_eigenvalues(schur_b);

    if eig_a.is_empty() || eig_b.is_empty() {
        return SeparationEstimate {
            separation: f64::INFINITY,
            is_singular: false,
            is_near_singular: false,
        };
    }

    let mut min_sep = f64::INFINITY;

    match equation_type {
        EquationType::Continuous => {
            // sep ≥ min |λ_i(A) + λ_j(B)|
            for &(a_re, a_im) in &eig_a {
                for &(b_re, b_im) in &eig_b {
                    // |λ_i(A) + λ_j(B)| for complex eigenvalues
                    // Need to consider all combinations of conjugate pairs
                    let sep1 = ((a_re + b_re).powi(2) + (a_im + b_im).powi(2)).sqrt();
                    let sep2 = ((a_re + b_re).powi(2) + (a_im - b_im).powi(2)).sqrt();
                    let sep = sep1.min(sep2);
                    min_sep = min_sep.min(sep);
                }
            }
        }
        EquationType::Discrete => {
            // sep ≥ min |1 + λ_i(A) * λ_j(B)|
            for &(a_re, a_im) in &eig_a {
                for &(b_re, b_im) in &eig_b {
                    // (a_re + i*a_im) * (b_re + i*b_im) = (a_re*b_re - a_im*b_im) + i*(a_re*b_im + a_im*b_re)
                    let prod_re = a_re * b_re - a_im * b_im;
                    let prod_im = a_re * b_im + a_im * b_re;
                    let sep1 = ((1.0 + prod_re).powi(2) + prod_im.powi(2)).sqrt();

                    // Also with conjugate: (a_re + i*a_im) * (b_re - i*b_im)
                    let prod_re2 = a_re * b_re + a_im * b_im;
                    let prod_im2 = -a_re * b_im + a_im * b_re;
                    let sep2 = ((1.0 + prod_re2).powi(2) + prod_im2.powi(2)).sqrt();

                    let sep = sep1.min(sep2);
                    min_sep = min_sep.min(sep);
                }
            }
        }
    }

    SeparationEstimate {
        separation: min_sep,
        is_singular: min_sep < SEPARATION_THRESHOLD,
        is_near_singular: min_sep < NEAR_SINGULAR_THRESHOLD,
    }
}

/// Extracts eigenvalues from a quasi-triangular (real Schur) matrix.
///
/// Returns eigenvalues as (real, imaginary) pairs.
/// - 1x1 diagonal blocks give real eigenvalues: (d, 0)
/// - 2x2 diagonal blocks give complex conjugate pairs: (re, im) and (re, -im),
///   but we only return the positive imaginary part since both are considered.
fn extract_eigenvalues(schur: MatRef<'_, f64>) -> Vec<(f64, f64)> {
    let n = schur.nrows();
    if n == 0 {
        return Vec::new();
    }

    let blocks = detect_blocks(schur, n);
    let mut eigenvalues = Vec::with_capacity(n);

    for (k1, k2) in blocks {
        if k1 == k2 {
            // 1x1 block: real eigenvalue
            eigenvalues.push((schur[(k1, k1)], 0.0));
        } else {
            // 2x2 block: complex conjugate pair
            // [[a, b], [c, d]] has eigenvalues ((a+d)/2) ± sqrt(((a-d)/2)^2 + b*c)
            let a = schur[(k1, k1)];
            let b = schur[(k1, k2)];
            let c = schur[(k2, k1)];
            let d = schur[(k2, k2)];

            let re = (a + d) / 2.0;
            let disc = ((a - d) / 2.0).powi(2) + b * c;

            if disc < 0.0 {
                // Complex conjugate pair
                let im = (-disc).sqrt();
                eigenvalues.push((re, im));
            } else {
                // Two real eigenvalues (unusual for 2x2 block in Schur form)
                let sq = disc.sqrt();
                eigenvalues.push((re + sq, 0.0));
                eigenvalues.push((re - sq, 0.0));
            }
        }
    }

    eigenvalues
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_extract_eigenvalues_diagonal() {
        let t = mat![[1.0, 0.0], [0.0, 3.0f64]];
        let eigs = extract_eigenvalues(t.as_ref());
        assert_eq!(eigs.len(), 2);
        assert!((eigs[0].0 - 1.0).abs() < 1e-14);
        assert!((eigs[0].1).abs() < 1e-14);
        assert!((eigs[1].0 - 3.0).abs() < 1e-14);
    }

    #[test]
    fn test_extract_eigenvalues_2x2_block() {
        // 2x2 block with complex eigenvalues: 1 ± 2i
        // [[1, 4], [-1, 1]] has eigenvalues 1 ± sqrt(4*(-1)) = 1 ± 2i
        let t = mat![[1.0, 4.0], [-1.0, 1.0f64]];
        let eigs = extract_eigenvalues(t.as_ref());
        assert_eq!(eigs.len(), 1); // One pair (re, im)
        assert!((eigs[0].0 - 1.0).abs() < 1e-14);
        assert!((eigs[0].1 - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_separation_continuous_well_separated() {
        // A has eigenvalue 1, B has eigenvalue 5
        // sep = |1 + 5| = 6
        let t = mat![[1.0f64]];
        let s = mat![[5.0f64]];
        let est = estimate_separation(t.as_ref(), s.as_ref(), EquationType::Continuous);
        assert!((est.separation - 6.0).abs() < 1e-14);
        assert!(!est.is_singular);
        assert!(!est.is_near_singular);
    }

    #[test]
    fn test_separation_continuous_common_eigenvalues() {
        // A has eigenvalue 2, B has eigenvalue -2
        // sep = |2 + (-2)| = 0
        let t = mat![[2.0f64]];
        let s = mat![[-2.0f64]];
        let est = estimate_separation(t.as_ref(), s.as_ref(), EquationType::Continuous);
        assert!(est.separation < 1e-14);
        assert!(est.is_singular);
    }

    #[test]
    fn test_separation_continuous_near_singular() {
        // A has eigenvalue 1, B has eigenvalue -1 + 1e-12
        // sep = |1 + (-1 + 1e-12)| = 1e-12
        let t = mat![[1.0f64]];
        let s = mat![[-1.0 + 1e-12f64]];
        let est = estimate_separation(t.as_ref(), s.as_ref(), EquationType::Continuous);
        assert!(est.separation < NEAR_SINGULAR_THRESHOLD);
        assert!(est.is_near_singular);
    }

    #[test]
    fn test_separation_continuous_complex_eigenvalues() {
        // A has eigenvalues 1 ± 2i (from 2x2 block)
        // B has eigenvalue 3 (1x1 block)
        // sep = min(|1+2i + 3|, |1-2i + 3|) = min(|(4+2i)|, |(4-2i)|) = sqrt(16+4) = sqrt(20)
        let t = mat![[1.0, 4.0], [-1.0, 1.0f64]]; // eigenvalues 1±2i
        let s = mat![[3.0f64]]; // eigenvalue 3
        let est = estimate_separation(t.as_ref(), s.as_ref(), EquationType::Continuous);
        let expected = (20.0f64).sqrt();
        assert!((est.separation - expected).abs() < 1e-12);
        assert!(!est.is_singular);
    }

    #[test]
    fn test_separation_discrete_well_separated() {
        // A has eigenvalue 0.5, B has eigenvalue 0.3
        // sep = |1 + 0.5*0.3| = |1.15| = 1.15
        let t = mat![[0.5f64]];
        let s = mat![[0.3f64]];
        let est = estimate_separation(t.as_ref(), s.as_ref(), EquationType::Discrete);
        assert!((est.separation - 1.15).abs() < 1e-14);
        assert!(!est.is_singular);
    }

    #[test]
    fn test_separation_discrete_singular() {
        // A has eigenvalue 2, B has eigenvalue -0.5
        // sep = |1 + 2*(-0.5)| = |1 - 1| = 0
        let t = mat![[2.0f64]];
        let s = mat![[-0.5f64]];
        let est = estimate_separation(t.as_ref(), s.as_ref(), EquationType::Discrete);
        assert!(est.separation < 1e-14);
        assert!(est.is_singular);
    }

    #[test]
    fn test_separation_empty_matrices() {
        let t = Mat::<f64>::zeros(0, 0);
        let s = Mat::zeros(0, 0);
        let est = estimate_separation(t.as_ref(), s.as_ref(), EquationType::Continuous);
        assert!(est.separation.is_infinite());
        assert!(!est.is_singular);
    }

    #[test]
    fn test_separation_multi_eigenvalue() {
        // A has eigenvalues 1, 3; B has eigenvalues -1, 5
        // Continuous sep = min(|1+(-1)|, |1+5|, |3+(-1)|, |3+5|) = min(0, 6, 2, 8) = 0
        let t = mat![[1.0, 0.0], [0.0, 3.0f64]];
        let s = mat![[-1.0, 0.0], [0.0, 5.0f64]];
        let est = estimate_separation(t.as_ref(), s.as_ref(), EquationType::Continuous);
        assert!(est.separation < 1e-14);
        assert!(est.is_singular);
    }
}
