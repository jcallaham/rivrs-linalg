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
pub struct TestMatrix {
    /// Metadata from metadata.json.
    pub metadata: MatrixMetadata,
    /// The sparse matrix loaded from .mtx file.
    pub matrix: SparseColMat<usize, f64>,
    /// Reference factorization from .json file (hand-constructed only).
    pub reference: Option<ReferenceFactorization>,
}

/// Top-level structure of metadata.json.
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
    let content = std::fs::read_to_string(&path).map_err(|e| SparseError::IoError {
        source: e.to_string(),
        path: path.display().to_string(),
    })?;
    let metadata: MetadataFile =
        serde_json::from_str(&content).map_err(|e| SparseError::ParseError {
            reason: e.to_string(),
            path: path.display().to_string(),
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
    let entry = registry
        .into_iter()
        .find(|m| m.name == name)
        .ok_or_else(|| SparseError::ParseError {
            reason: format!("matrix '{}' not found in registry", name),
            path: test_data_dir().join("metadata.json").display().to_string(),
            line: None,
        })?;

    // Resolve .mtx path. CI-subset matrices have copies in suitesparse-ci/
    // which is committed to git, while their metadata path points to the
    // gitignored suitesparse/ directory.
    let mtx_path = resolve_mtx_path(&entry);

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
        metadata: entry,
        matrix,
        reference,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_registry_returns_82_entries() {
        let registry = load_registry().expect("failed to load registry");
        assert_eq!(registry.len(), 82);
    }

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
    }

    #[test]
    fn gitignored_matrix_returns_none() {
        // Find a matrix in the registry that has in_repo: true but lives in the
        // gitignored suitesparse/ directory. We'll use a known SuiteSparse name.
        let registry = load_registry().expect("failed to load registry");
        let suitesparse_entry = registry
            .iter()
            .find(|m| m.path.starts_with("suitesparse/") && !m.path.starts_with("suitesparse-ci/"));

        if let Some(entry) = suitesparse_entry {
            let mtx_path = test_data_dir().join(&entry.path);
            if !mtx_path.exists() {
                // This matrix is gitignored and not present — should return None
                let result = load_test_matrix(&entry.name).expect("registry error");
                assert!(result.is_none(), "gitignored matrix should return None");
            }
            // If the file happens to exist (extracted), we can't test the None case
        }
    }
}
