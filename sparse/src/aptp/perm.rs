//! Permutation construction utilities for APTP factorization.
//!
//! Provides [`perm_from_forward`] to bridge ordering algorithm output (forward
//! permutation arrays) to faer's [`Perm<usize>`](faer::perm::Perm) type by
//! computing the inverse mapping automatically.
//!
//! Reference: Hogg, Duff & Lopez (2020), "A New Sparse LDL^T Solver Using
//! A Posteriori Threshold Pivoting", SIAM J. Sci. Comput. 42(4) — permutation
//! handling in the context of fill-reducing orderings (Section 2.1).

use faer::perm::Perm;

use crate::error::SparseError;
use crate::validate;

/// Construct a [`Perm<usize>`] from a forward permutation array, computing
/// the inverse automatically.
///
/// Takes ownership of the forward array to avoid an extra allocation when
/// converting to `Box<[usize]>` for faer's `Perm::new_checked`.
///
/// # Errors
///
/// Returns [`SparseError::InvalidInput`] if the array contains duplicate
/// indices or out-of-bounds values.
///
/// # Examples
///
/// ```
/// use rivrs_sparse::aptp::perm_from_forward;
///
/// let fwd = vec![2, 0, 1];
/// let perm = perm_from_forward(fwd).unwrap();
/// let (fwd_arr, inv_arr) = perm.as_ref().arrays();
/// assert_eq!(fwd_arr, &[2, 0, 1]);
/// assert_eq!(inv_arr, &[1, 2, 0]);
/// ```
pub fn perm_from_forward(fwd: Vec<usize>) -> Result<Perm<usize>, SparseError> {
    let n = fwd.len();
    validate::validate_permutation(&fwd, n)?;

    // Compute inverse: inv[fwd[i]] = i
    let mut inv = vec![0usize; n];
    for (i, &f) in fwd.iter().enumerate() {
        inv[f] = i;
    }

    let fwd_box = fwd.into_boxed_slice();
    let inv_box = inv.into_boxed_slice();
    Ok(Perm::new_checked(fwd_box, inv_box, n))
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- perm_from_forward unit tests ----

    #[test]
    fn identity_permutation() {
        let fwd = vec![0, 1, 2, 3, 4];
        let perm = perm_from_forward(fwd).unwrap();
        let (fwd_arr, inv_arr) = perm.as_ref().arrays();
        assert_eq!(fwd_arr, &[0, 1, 2, 3, 4]);
        assert_eq!(inv_arr, &[0, 1, 2, 3, 4]);
    }

    #[test]
    fn nontrivial_permutation_round_trip() {
        let fwd = vec![3, 1, 4, 0, 2];
        let perm = perm_from_forward(fwd).unwrap();
        let (fwd_arr, inv_arr) = perm.as_ref().arrays();

        // Verify round-trip: inv[fwd[i]] == i for all i
        for i in 0..5 {
            assert_eq!(inv_arr[fwd_arr[i]], i);
        }
        // And fwd[inv[i]] == i for all i
        for i in 0..5 {
            assert_eq!(fwd_arr[inv_arr[i]], i);
        }
    }

    #[test]
    fn composition_of_two_permutations() {
        let p1 = perm_from_forward(vec![1, 0, 2]).unwrap();
        let p2 = perm_from_forward(vec![2, 1, 0]).unwrap();
        let (p1_fwd, _) = p1.as_ref().arrays();
        let (p2_fwd, _) = p2.as_ref().arrays();

        // Compose: apply p1 then p2 → composed_fwd[i] = p2_fwd[p1_fwd[i]]
        let composed: Vec<usize> = (0..3).map(|i| p2_fwd[p1_fwd[i]]).collect();
        let composed_perm = perm_from_forward(composed).unwrap();
        let (comp_fwd, comp_inv) = composed_perm.as_ref().arrays();

        // Verify: composition of swap(0,1) then reverse → should be [2, 0, 1]
        // p1: [1, 0, 2] — swap 0 and 1
        // p2: [2, 1, 0] — reverse
        // composed[0] = p2[p1[0]] = p2[1] = 1
        // composed[1] = p2[p1[1]] = p2[0] = 2
        // composed[2] = p2[p1[2]] = p2[2] = 0
        assert_eq!(comp_fwd, &[1, 2, 0]);

        // Verify forward/inverse relationship
        for i in 0..3 {
            assert_eq!(comp_inv[comp_fwd[i]], i);
        }
    }

    #[test]
    fn empty_array_dimension_0() {
        let perm = perm_from_forward(vec![]).unwrap();
        let (fwd_arr, inv_arr) = perm.as_ref().arrays();
        assert!(fwd_arr.is_empty());
        assert!(inv_arr.is_empty());
    }

    #[test]
    fn single_element() {
        let perm = perm_from_forward(vec![0]).unwrap();
        let (fwd_arr, inv_arr) = perm.as_ref().arrays();
        assert_eq!(fwd_arr, &[0]);
        assert_eq!(inv_arr, &[0]);
    }

    #[test]
    fn error_on_duplicate_indices() {
        let result = perm_from_forward(vec![0, 1, 1]);
        assert!(result.is_err());
    }

    #[test]
    fn error_on_out_of_bounds() {
        let result = perm_from_forward(vec![0, 5, 2]);
        assert!(result.is_err());
    }

    #[test]
    fn forward_inverse_relationship() {
        let fwd = vec![4, 2, 0, 3, 1];
        let perm = perm_from_forward(fwd).unwrap();
        let (fwd_arr, inv_arr) = perm.as_ref().arrays();

        for i in 0..5 {
            assert_eq!(inv_arr[fwd_arr[i]], i, "inv[fwd[{}]] != {}", i, i);
            assert_eq!(fwd_arr[inv_arr[i]], i, "fwd[inv[{}]] != {}", i, i);
        }
    }

    #[test]
    fn scale_test_n_10000() {
        // perm_from_forward at n=10,000, verify round-trip identity
        let n = 10_000;
        // Create a permutation: rotate by 1
        let fwd: Vec<usize> = (0..n).map(|i| (i + 1) % n).collect();
        let perm = perm_from_forward(fwd).unwrap();
        let (fwd_arr, inv_arr) = perm.as_ref().arrays();

        // Verify round-trip
        for i in 0..n {
            assert_eq!(inv_arr[fwd_arr[i]], i);
            assert_eq!(fwd_arr[inv_arr[i]], i);
        }
    }
}
