//! Test matrix registry backed by metadata.json.
//!
//! Provides functions to load the test matrix catalog and individual matrices
//! by name. Uses `CARGO_MANIFEST_DIR` to locate the `test-data/` directory.

use std::path::{Path, PathBuf};

use faer::sparse::SparseColMat;
use serde::Deserialize;

use crate::error::SparseError;
use crate::io::mtx;
use crate::io::reference::{self, ReferenceFactorization};

/// Structural and numerical properties of a test matrix.
#[derive(Debug, Clone, Deserialize)]
pub struct MatrixProperties {
    /// Always true for this project.
    pub symmetric: bool,
    /// Whether matrix is positive definite.
    #[serde(default)]
    pub positive_definite: bool,
    /// Whether matrix is indefinite.
    #[serde(default)]
    pub indefinite: bool,
    /// Difficulty classification: "trivial", "easy", "hard".
    #[serde(default)]
    pub difficulty: String,
    /// Structure type: "arrow", "tridiagonal", "block-diagonal", etc.
    #[serde(default)]
    pub structure: Option<String>,
    /// Kind of problem (SuiteSparse matrices).
    #[serde(default)]
    pub kind: Option<String>,
    /// Expected delayed pivot behavior (SuiteSparse matrices).
    #[serde(default)]
    pub expected_delayed_pivots: Option<String>,
}

/// Metadata for a single test matrix from metadata.json.
#[derive(Debug, Clone, Deserialize)]
pub struct MatrixMetadata {
    /// Matrix identifier (e.g., "arrow-10-indef").
    pub name: String,
    /// Origin: "hand-constructed" or "suitesparse".
    pub source: String,
    /// Classification: "hand-constructed", "easy-indefinite", "hard-indefinite", "positive-definite".
    pub category: String,
    /// Relative path from test-data/ to .mtx file.
    pub path: String,
    /// Matrix dimension (n for n×n).
    pub size: usize,
    /// Number of stored nonzeros (lower triangle for symmetric).
    pub nnz: usize,
    /// Whether the .mtx file is committed to git.
    #[serde(default)]
    pub in_repo: bool,
    /// Whether this matrix is in the CI subset.
    #[serde(default)]
    pub ci_subset: bool,
    /// Structural/numerical properties.
    pub properties: MatrixProperties,
    /// Academic paper references.
    #[serde(default)]
    pub paper_references: Vec<String>,
    /// Reference solver results (currently empty, reserved for future use).
    #[serde(default)]
    pub reference_results: serde_json::Value,
    /// Relative path to companion .json factorization file (hand-constructed only).
    #[serde(default)]
    pub factorization_path: Option<String>,
}

/// A loaded test matrix with metadata and optional reference factorization.
#[derive(Debug)]
pub struct TestMatrix {
    /// Metadata from metadata.json.
    pub metadata: MatrixMetadata,
    /// The sparse matrix loaded from .mtx file.
    pub matrix: SparseColMat<usize, f64>,
    /// Reference factorization from .json file (hand-constructed only).
    pub reference: Option<ReferenceFactorization>,
}

/// Top-level structure of metadata.json.
///
/// The `schema_version`, `generated`, and `total_count` fields are parsed
/// for forward-compatibility but not currently used by the library.
#[derive(Debug, Deserialize)]
struct MetadataFile {
    #[allow(dead_code)]
    schema_version: String,
    #[allow(dead_code)]
    generated: String,
    #[allow(dead_code)]
    total_count: usize,
    matrices: Vec<MatrixMetadata>,
}

/// Return the absolute path to the test-data/ directory.
fn test_data_dir() -> PathBuf {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    Path::new(manifest_dir).join("test-data")
}

/// Load the test matrix registry from metadata.json.
///
/// # Errors
///
/// - metadata.json not found at expected location
/// - Invalid JSON structure
pub fn load_registry() -> Result<Vec<MatrixMetadata>, SparseError> {
    let path = test_data_dir().join("metadata.json");
    let path_str = path.display().to_string();
    let content = std::fs::read_to_string(&path).map_err(|e| SparseError::IoError {
        source: e.to_string(),
        path: path_str.clone(),
    })?;
    let metadata: MetadataFile =
        serde_json::from_str(&content).map_err(|e| SparseError::ParseError {
            reason: e.to_string(),
            path: path_str,
            line: None,
        })?;
    Ok(metadata.matrices)
}

/// Resolve the .mtx file path for a matrix entry.
///
/// For CI-subset matrices, prefers the `suitesparse-ci/<category>/<name>.mtx`
/// copy (committed to git) over the gitignored `suitesparse/` path.
fn resolve_mtx_path(entry: &MatrixMetadata) -> PathBuf {
    if entry.ci_subset {
        if let Some(ci_path) = ci_subset_path(entry) {
            if ci_path.exists() {
                return ci_path;
            }
        }
    }
    test_data_dir().join(&entry.path)
}

/// Build the suitesparse-ci/ path for a CI-subset matrix.
///
/// Metadata paths are `suitesparse/<category>/<name>/<name>.mtx`.
/// CI copies are at `suitesparse-ci/<category>/<name>.mtx`.
fn ci_subset_path(entry: &MatrixMetadata) -> Option<PathBuf> {
    let rest = entry.path.strip_prefix("suitesparse/")?;
    let category = rest.split('/').next()?;
    let file_name = Path::new(&entry.path).file_name()?;
    Some(
        test_data_dir()
            .join("suitesparse-ci")
            .join(category)
            .join(file_name),
    )
}

/// Load a test matrix from a pre-resolved registry entry.
///
/// This avoids re-parsing metadata.json when loading multiple matrices
/// from an already-loaded registry. Returns `Ok(None)` if the .mtx file
/// does not exist on disk (e.g., gitignored SuiteSparse matrix not extracted).
///
/// # Errors
///
/// - `.mtx` file exists but fails to parse
/// - `.json` file exists but fails to parse
/// - Reference factorization permutation length doesn't match matrix dimension
pub fn load_test_matrix_from_entry(
    entry: &MatrixMetadata,
) -> Result<Option<TestMatrix>, SparseError> {
    let mtx_path = resolve_mtx_path(entry);

    if !mtx_path.exists() {
        return Ok(None);
    }

    let matrix = mtx::load_mtx(&mtx_path)?;

    // Load reference factorization if available
    let reference = if let Some(ref fact_path) = entry.factorization_path {
        let json_path = test_data_dir().join(fact_path);
        if json_path.exists() {
            let refdata = reference::load_reference(&json_path)?;
            // Cross-check: permutation length must match matrix dimension
            if refdata.permutation.len() != matrix.nrows() {
                return Err(SparseError::ParseError {
                    reason: format!(
                        "reference factorization permutation length ({}) != matrix dimension ({})",
                        refdata.permutation.len(),
                        matrix.nrows()
                    ),
                    path: json_path.display().to_string(),
                    line: None,
                });
            }
            Some(refdata)
        } else {
            None
        }
    } else {
        None
    };

    Ok(Some(TestMatrix {
        metadata: entry.clone(),
        matrix,
        reference,
    }))
}

/// Load a specific test matrix by name, including its sparse matrix
/// and optional reference factorization.
///
/// Returns `Ok(None)` if the matrix's `.mtx` file does not exist on disk
/// (e.g., gitignored SuiteSparse matrix not extracted).
///
/// # Errors
///
/// - Matrix name not found in registry
/// - `.mtx` file exists but fails to parse
/// - `.json` file exists but fails to parse
pub fn load_test_matrix(name: &str) -> Result<Option<TestMatrix>, SparseError> {
    let registry = load_registry()?;
    let entry =
        registry
            .iter()
            .find(|m| m.name == name)
            .ok_or_else(|| SparseError::MatrixNotFound {
                name: name.to_string(),
            })?;

    load_test_matrix_from_entry(entry)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_arrow_5_pd_returns_some() {
        let test = load_test_matrix("arrow-5-pd")
            .expect("registry error")
            .expect("matrix should exist on disk");
        assert_eq!(test.matrix.nrows(), 5);
        assert_eq!(test.matrix.ncols(), 5);
        assert!(test.reference.is_some());
    }

    #[test]
    fn nonexistent_matrix_returns_error() {
        let result = load_test_matrix("nonexistent-matrix-name");
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SparseError::MatrixNotFound { .. }
        ));
    }

    #[test]
    fn missing_mtx_file_returns_none() {
        // Construct a fake entry pointing to a nonexistent file
        let fake_entry = MatrixMetadata {
            name: "fake-missing-matrix".to_string(),
            source: "test".to_string(),
            category: "test".to_string(),
            path: "nonexistent/path/fake.mtx".to_string(),
            size: 5,
            nnz: 10,
            in_repo: false,
            ci_subset: false,
            properties: MatrixProperties {
                symmetric: true,
                positive_definite: false,
                indefinite: false,
                difficulty: "trivial".to_string(),
                structure: None,
                kind: None,
                expected_delayed_pivots: None,
            },
            paper_references: vec![],
            reference_results: serde_json::Value::Null,
            factorization_path: None,
        };
        let result =
            load_test_matrix_from_entry(&fake_entry).expect("should not error for missing file");
        assert!(result.is_none(), "missing .mtx file should return None");
    }

    #[test]
    fn load_via_entry_matches_load_by_name() {
        let registry = load_registry().expect("failed to load registry");
        let entry = registry.iter().find(|m| m.name == "arrow-5-pd").unwrap();

        let by_entry = load_test_matrix_from_entry(entry)
            .expect("entry load error")
            .expect("should exist");
        let by_name = load_test_matrix("arrow-5-pd")
            .expect("name load error")
            .expect("should exist");

        assert_eq!(by_entry.matrix.nrows(), by_name.matrix.nrows());
        assert_eq!(by_entry.matrix.ncols(), by_name.matrix.ncols());
        assert_eq!(by_entry.matrix.compute_nnz(), by_name.matrix.compute_nnz());
    }
}
