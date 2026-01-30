use faer::prelude::*;

use crate::error::SylvesterError;

/// Validates that matrices A, B, C have compatible dimensions for a Sylvester equation.
///
/// Requirements:
/// - A must be square (n x n)
/// - B must be square (m x m)
/// - C must be n x m
pub fn validate_dimensions(
    a: MatRef<'_, f64>,
    b: MatRef<'_, f64>,
    c: MatRef<'_, f64>,
) -> Result<(), SylvesterError> {
    if a.nrows() != a.ncols() {
        return Err(SylvesterError::NotSquare {
            matrix: 'A',
            dims: (a.nrows(), a.ncols()),
        });
    }

    if b.nrows() != b.ncols() {
        return Err(SylvesterError::NotSquare {
            matrix: 'B',
            dims: (b.nrows(), b.ncols()),
        });
    }

    if c.nrows() != a.nrows() || c.ncols() != b.ncols() {
        return Err(SylvesterError::DimensionMismatch {
            expected: (a.nrows(), b.ncols()),
            got: (c.nrows(), c.ncols()),
            context: "Matrix C must be n x m where A is n x n and B is m x m".into(),
        });
    }

    Ok(())
}

/// Validates that all elements of a matrix are finite (no NaN or Inf).
pub fn validate_finite(matrix: MatRef<'_, f64>, name: char) -> Result<(), SylvesterError> {
    for j in 0..matrix.ncols() {
        for i in 0..matrix.nrows() {
            let val = matrix[(i, j)];
            if !val.is_finite() {
                return Err(SylvesterError::InvalidInput {
                    reason: format!(
                        "Matrix {} contains non-finite value at ({}, {})",
                        name, i, j
                    ),
                });
            }
        }
    }
    Ok(())
}

/// Validates that a matrix is in quasi-triangular (real Schur) form.
///
/// A quasi-triangular matrix is upper triangular except for possible
/// 2x2 blocks on the diagonal that represent complex conjugate eigenvalue pairs.
///
/// The 2x2 blocks have the property that the sub-diagonal element is nonzero,
/// and consecutive 2x2 blocks cannot overlap.
pub fn validate_quasi_triangular(
    matrix: MatRef<'_, f64>,
    name: char,
) -> Result<(), SylvesterError> {
    let n = matrix.nrows();
    if n <= 1 {
        return Ok(());
    }

    let tol = f64::EPSILON * 100.0;

    let mut i = 0;
    while i < n {
        // Check if this starts a 2x2 block
        let is_2x2_block = if i + 1 < n {
            matrix[(i + 1, i)].abs() > tol * (matrix[(i, i)].abs() + matrix[(i + 1, i + 1)].abs())
        } else {
            false
        };

        if is_2x2_block {
            // For a 2x2 block, check that everything below row i+1 in columns i and i+1 is zero
            for k in (i + 2)..n {
                if matrix[(k, i)].abs() > tol {
                    return Err(SylvesterError::NotQuasiTriangular {
                        matrix: name,
                        location: (k, i),
                    });
                }
                if matrix[(k, i + 1)].abs() > tol {
                    return Err(SylvesterError::NotQuasiTriangular {
                        matrix: name,
                        location: (k, i + 1),
                    });
                }
            }
            i += 2; // Skip the 2x2 block
        } else {
            // 1x1 block - check below diagonal in column i
            for k in (i + 1)..n {
                if matrix[(k, i)].abs() > tol {
                    return Err(SylvesterError::NotQuasiTriangular {
                        matrix: name,
                        location: (k, i),
                    });
                }
            }
            i += 1;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_validate_dimensions_valid() {
        let a = mat![[1.0, 2.0], [3.0, 4.0]];
        let b = mat![[5.0, 6.0, 7.0], [8.0, 9.0, 10.0], [11.0, 12.0, 13.0]];
        let c = Mat::zeros(2, 3);
        assert!(validate_dimensions(a.as_ref(), b.as_ref(), c.as_ref()).is_ok());
    }

    #[test]
    fn test_validate_dimensions_a_not_square() {
        let a = Mat::<f64>::zeros(2, 3);
        let b = Mat::zeros(3, 3);
        let c = Mat::zeros(2, 3);
        let err = validate_dimensions(a.as_ref(), b.as_ref(), c.as_ref()).unwrap_err();
        match err {
            SylvesterError::NotSquare { matrix, .. } => assert_eq!(matrix, 'A'),
            _ => panic!("Expected NotSquare error for A"),
        }
    }

    #[test]
    fn test_validate_dimensions_b_not_square() {
        let a = Mat::<f64>::zeros(2, 2);
        let b = Mat::zeros(2, 3);
        let c = Mat::zeros(2, 3);
        let err = validate_dimensions(a.as_ref(), b.as_ref(), c.as_ref()).unwrap_err();
        match err {
            SylvesterError::NotSquare { matrix, .. } => assert_eq!(matrix, 'B'),
            _ => panic!("Expected NotSquare error for B"),
        }
    }

    #[test]
    fn test_validate_dimensions_c_mismatch() {
        let a = Mat::<f64>::zeros(2, 2);
        let b = Mat::zeros(3, 3);
        let c = Mat::zeros(3, 3); // Wrong: should be 2x3
        let err = validate_dimensions(a.as_ref(), b.as_ref(), c.as_ref()).unwrap_err();
        match err {
            SylvesterError::DimensionMismatch { expected, got, .. } => {
                assert_eq!(expected, (2, 3));
                assert_eq!(got, (3, 3));
            }
            _ => panic!("Expected DimensionMismatch error"),
        }
    }

    #[test]
    fn test_validate_finite_ok() {
        let m = mat![[1.0, 2.0], [3.0, 4.0]];
        assert!(validate_finite(m.as_ref(), 'A').is_ok());
    }

    #[test]
    fn test_validate_finite_nan() {
        let m = mat![[1.0, f64::NAN], [3.0, 4.0]];
        assert!(validate_finite(m.as_ref(), 'A').is_err());
    }

    #[test]
    fn test_validate_finite_inf() {
        let m = mat![[1.0, 2.0], [f64::INFINITY, 4.0]];
        assert!(validate_finite(m.as_ref(), 'B').is_err());
    }

    #[test]
    fn test_validate_quasi_triangular_diagonal() {
        let m = mat![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        assert!(validate_quasi_triangular(m.as_ref(), 'T').is_ok());
    }

    #[test]
    fn test_validate_quasi_triangular_upper_triangular() {
        let m = mat![[1.0, 2.0, 3.0], [0.0, 4.0, 5.0], [0.0, 0.0, 6.0]];
        assert!(validate_quasi_triangular(m.as_ref(), 'T').is_ok());
    }

    #[test]
    fn test_validate_quasi_triangular_with_2x2_block() {
        // 2x2 block at top-left representing complex eigenvalues
        let m = mat![
            [1.0, 2.0, 3.0],
            [-0.5, 1.0, 4.0],
            [0.0, 0.0, 5.0]
        ];
        assert!(validate_quasi_triangular(m.as_ref(), 'T').is_ok());
    }

    #[test]
    fn test_validate_quasi_triangular_not_valid() {
        // Lower triangular part has nonzero entry not part of a 2x2 block
        let m = mat![[1.0, 2.0, 3.0], [0.0, 4.0, 5.0], [1.0, 0.0, 6.0]];
        assert!(validate_quasi_triangular(m.as_ref(), 'T').is_err());
    }

    #[test]
    fn test_validate_empty_matrices() {
        let a = Mat::<f64>::zeros(0, 0);
        let b = Mat::zeros(0, 0);
        let c = Mat::zeros(0, 0);
        assert!(validate_dimensions(a.as_ref(), b.as_ref(), c.as_ref()).is_ok());
    }
}
