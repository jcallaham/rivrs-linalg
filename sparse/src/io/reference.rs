//! Reference factorization data loader (JSON format).
//!
//! Loads companion `.json` files for hand-constructed test matrices that contain
//! analytically known LDL^T factorizations (L factor, D diagonal, permutation,
//! and inertia).

use std::path::Path;

use serde::Deserialize;

use crate::error::SparseError;

/// Eigenvalue sign classification of a symmetric matrix.
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

/// Single entry in the strict lower triangle of L.
#[derive(Debug, Clone, Deserialize)]
pub struct LEntry {
    /// Row index (0-indexed).
    pub row: usize,
    /// Column index (0-indexed), must satisfy col < row.
    pub col: usize,
    /// Entry value.
    pub value: f64,
}

/// One block of the block diagonal D.
///
/// Either a 1×1 scalar pivot or a 2×2 symmetric pivot block.
#[derive(Debug, Clone)]
pub enum DBlock {
    /// 1×1 scalar pivot.
    OneByOne { value: f64 },
    /// 2×2 symmetric pivot block (row-major).
    TwoByTwo { values: [[f64; 2]; 2] },
}

/// The known-correct LDL^T factorization of a hand-constructed matrix.
#[derive(Debug, Clone, Deserialize)]
pub struct ReferenceFactorization {
    /// Matrix name (must match the MatrixMetadata name).
    pub matrix_name: String,
    /// Pivot permutation (0-indexed).
    pub permutation: Vec<usize>,
    /// Strict lower-triangular entries of L.
    pub l_entries: Vec<LEntry>,
    /// Block diagonal D (1×1 or 2×2 blocks).
    pub d_blocks: Vec<DBlock>,
    /// Eigenvalue sign counts.
    pub inertia: Inertia,
    /// Human-readable description.
    #[serde(default)]
    pub notes: String,
}

/// Load a reference factorization from a companion JSON file.
///
/// # Errors
///
/// - File not found or unreadable
/// - Invalid JSON structure
/// - Inconsistent data (l_entry indices out of bounds, invalid permutation)
pub fn load_reference(path: &Path) -> Result<ReferenceFactorization, SparseError> {
    let content = std::fs::read_to_string(path).map_err(|e| SparseError::IoError {
        source: e.to_string(),
        path: path.display().to_string(),
    })?;

    let refdata: ReferenceFactorization =
        serde_json::from_str(&content).map_err(|e| SparseError::ParseError {
            reason: e.to_string(),
            path: path.display().to_string(),
            line: None,
        })?;

    // Validate l_entries: col < row
    for (i, entry) in refdata.l_entries.iter().enumerate() {
        if entry.col >= entry.row {
            return Err(SparseError::ParseError {
                reason: format!(
                    "l_entry[{}] has col ({}) >= row ({}); must be strict lower triangle",
                    i, entry.col, entry.row
                ),
                path: path.display().to_string(),
                line: None,
            });
        }
    }

    // Validate permutation: must be a valid permutation of 0..n
    let n = refdata.permutation.len();
    let mut seen = vec![false; n];
    for (i, &p) in refdata.permutation.iter().enumerate() {
        if p >= n {
            return Err(SparseError::ParseError {
                reason: format!("permutation[{}] = {} is out of bounds (n = {})", i, p, n),
                path: path.display().to_string(),
                line: None,
            });
        }
        if seen[p] {
            return Err(SparseError::ParseError {
                reason: format!("permutation has duplicate index {}", p),
                path: path.display().to_string(),
                line: None,
            });
        }
        seen[p] = true;
    }

    Ok(refdata)
}

// Custom serde for DBlock: the JSON uses {"size": 1, "values": [...]} or {"size": 2, "values": [[...], [...]]}
impl<'de> Deserialize<'de> for DBlock {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::Error;

        let raw: serde_json::Value = Deserialize::deserialize(deserializer)?;
        let obj = raw
            .as_object()
            .ok_or_else(|| D::Error::custom("d_block must be an object"))?;

        let size = obj
            .get("size")
            .and_then(|v| v.as_u64())
            .ok_or_else(|| D::Error::custom("d_block must have integer 'size' field"))?;

        let values = obj
            .get("values")
            .ok_or_else(|| D::Error::custom("d_block must have 'values' field"))?;

        match size {
            1 => {
                // values is an array with one element: [scalar]
                let arr = values
                    .as_array()
                    .ok_or_else(|| D::Error::custom("1x1 d_block values must be an array"))?;
                if arr.len() != 1 {
                    return Err(D::Error::custom(format!(
                        "1x1 d_block values must have exactly 1 element, got {}",
                        arr.len()
                    )));
                }
                let value = arr[0]
                    .as_f64()
                    .ok_or_else(|| D::Error::custom("1x1 d_block value must be a number"))?;
                Ok(DBlock::OneByOne { value })
            }
            2 => {
                // values is a 2x2 array: [[a, b], [c, d]]
                let arr = values
                    .as_array()
                    .ok_or_else(|| D::Error::custom("2x2 d_block values must be an array"))?;
                if arr.len() != 2 {
                    return Err(D::Error::custom(format!(
                        "2x2 d_block values must have exactly 2 rows, got {}",
                        arr.len()
                    )));
                }
                let mut vals = [[0.0f64; 2]; 2];
                for (i, row) in arr.iter().enumerate() {
                    let row_arr = row.as_array().ok_or_else(|| {
                        D::Error::custom(format!("2x2 d_block row {} must be an array", i))
                    })?;
                    if row_arr.len() != 2 {
                        return Err(D::Error::custom(format!(
                            "2x2 d_block row {} must have exactly 2 elements, got {}",
                            i,
                            row_arr.len()
                        )));
                    }
                    for (j, val) in row_arr.iter().enumerate() {
                        vals[i][j] = val.as_f64().ok_or_else(|| {
                            D::Error::custom(format!(
                                "2x2 d_block value at ({}, {}) must be a number",
                                i, j
                            ))
                        })?;
                    }
                }
                Ok(DBlock::TwoByTwo { values: vals })
            }
            _ => Err(D::Error::custom(format!(
                "d_block size must be 1 or 2, got {}",
                size
            ))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test-data")
    }

    #[test]
    fn load_arrow_5_pd_reference() {
        let path = test_data_dir().join("hand-constructed/arrow-5-pd.json");
        let refdata = load_reference(&path).expect("failed to load arrow-5-pd.json");
        assert_eq!(
            refdata.inertia,
            Inertia {
                positive: 5,
                negative: 0,
                zero: 0
            }
        );
        assert_eq!(refdata.l_entries.len(), 10);
        assert_eq!(refdata.permutation.len(), 5);
        assert_eq!(refdata.d_blocks.len(), 5);
        // All d_blocks should be 1x1 for this PD matrix
        for block in &refdata.d_blocks {
            assert!(matches!(block, DBlock::OneByOne { .. }));
        }
    }

    #[test]
    fn load_stress_delayed_pivots_2x2_blocks() {
        let path = test_data_dir().join("hand-constructed/stress-delayed-pivots.json");
        let refdata = load_reference(&path).expect("failed to load stress-delayed-pivots.json");
        assert_eq!(refdata.d_blocks.len(), 5);
        // All d_blocks should be 2x2 for this stress test
        for block in &refdata.d_blocks {
            assert!(matches!(block, DBlock::TwoByTwo { .. }));
        }
        assert_eq!(
            refdata.inertia,
            Inertia {
                positive: 5,
                negative: 5,
                zero: 0
            }
        );
    }

    #[test]
    fn invalid_json_returns_error() {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test-tmp");
        std::fs::create_dir_all(&dir).ok();
        let path = dir.join("invalid.json");
        std::fs::write(&path, "{ not valid json }").unwrap();
        let result = load_reference(&path);
        assert!(result.is_err());
    }
}
