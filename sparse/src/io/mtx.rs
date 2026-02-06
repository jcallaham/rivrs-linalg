//! Matrix Market (.mtx) file parser.
//!
//! Parses `coordinate real symmetric` format files into faer sparse matrices.
//! This is a lightweight custom parser scoped to the format used by all test
//! matrices in this project.

use std::path::Path;

use faer::sparse::{SparseColMat, Triplet};

use crate::error::SparseError;

/// Parse a Matrix Market file into a faer sparse symmetric matrix.
///
/// Supports `coordinate real symmetric` format only. Entries are
/// 1-indexed in the file and converted to 0-indexed internally.
/// The returned matrix stores the full symmetric matrix (both triangles).
///
/// # Errors
///
/// - File not found or unreadable
/// - Unsupported Matrix Market format (not coordinate real symmetric)
/// - Malformed header, size line, or data lines
/// - Row/col indices out of bounds
pub fn load_mtx(path: &Path) -> Result<SparseColMat<usize, f64>, SparseError> {
    let path_str = path.display().to_string();
    let content = std::fs::read_to_string(path).map_err(|e| SparseError::IoError {
        source: e.to_string(),
        path: path_str.clone(),
    })?;

    let mut lines = content.lines().enumerate();

    // Parse header line
    let (_, header) = lines.next().ok_or_else(|| SparseError::ParseError {
        reason: "empty file".to_string(),
        path: path_str.clone(),
        line: Some(1),
    })?;

    let header_lower = header.to_lowercase();
    if !header_lower.starts_with("%%matrixmarket matrix coordinate real symmetric") {
        return Err(SparseError::ParseError {
            reason: format!(
                "unsupported format: expected '%%MatrixMarket matrix coordinate real symmetric', got '{}'",
                header
            ),
            path: path_str.clone(),
            line: Some(1),
        });
    }

    // Skip comment lines
    let mut size_line = None;
    for (line_idx, line) in &mut lines {
        if !line.starts_with('%') {
            size_line = Some((line_idx, line));
            break;
        }
    }

    // Parse size line: nrows ncols nnz
    let (size_line_idx, size_line_str) = size_line.ok_or_else(|| SparseError::ParseError {
        reason: "missing size line".to_string(),
        path: path_str.clone(),
        line: None,
    })?;

    let size_parts: Vec<&str> = size_line_str.split_whitespace().collect();
    if size_parts.len() != 3 {
        return Err(SparseError::ParseError {
            reason: format!(
                "size line must have 3 values (nrows ncols nnz), got {}",
                size_parts.len()
            ),
            path: path_str.clone(),
            line: Some(size_line_idx + 1),
        });
    }

    let nrows: usize = size_parts[0].parse().map_err(|_| SparseError::ParseError {
        reason: format!("invalid nrows: '{}'", size_parts[0]),
        path: path_str.clone(),
        line: Some(size_line_idx + 1),
    })?;
    let ncols: usize = size_parts[1].parse().map_err(|_| SparseError::ParseError {
        reason: format!("invalid ncols: '{}'", size_parts[1]),
        path: path_str.clone(),
        line: Some(size_line_idx + 1),
    })?;
    let _nnz: usize = size_parts[2].parse().map_err(|_| SparseError::ParseError {
        reason: format!("invalid nnz: '{}'", size_parts[2]),
        path: path_str.clone(),
        line: Some(size_line_idx + 1),
    })?;

    // Parse triplets and symmetrize
    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
    for (line_idx, line) in lines {
        let line = line.trim();
        if line.is_empty() || line.starts_with('%') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() != 3 {
            return Err(SparseError::ParseError {
                reason: format!(
                    "data line must have 3 values (row col val), got {}",
                    parts.len()
                ),
                path: path_str.clone(),
                line: Some(line_idx + 1),
            });
        }

        let row: usize = parts[0].parse().map_err(|_| SparseError::ParseError {
            reason: format!("invalid row index: '{}'", parts[0]),
            path: path_str.clone(),
            line: Some(line_idx + 1),
        })?;
        let col: usize = parts[1].parse().map_err(|_| SparseError::ParseError {
            reason: format!("invalid col index: '{}'", parts[1]),
            path: path_str.clone(),
            line: Some(line_idx + 1),
        })?;
        let val: f64 = parts[2].parse().map_err(|_| SparseError::ParseError {
            reason: format!("invalid value: '{}'", parts[2]),
            path: path_str.clone(),
            line: Some(line_idx + 1),
        })?;

        // Convert from 1-indexed to 0-indexed
        if row == 0 || col == 0 {
            return Err(SparseError::ParseError {
                reason: "Matrix Market indices are 1-based; got 0".to_string(),
                path: path_str.clone(),
                line: Some(line_idx + 1),
            });
        }
        let row = row - 1;
        let col = col - 1;

        if row >= nrows || col >= ncols {
            return Err(SparseError::ParseError {
                reason: format!(
                    "index ({}, {}) out of bounds for {}x{} matrix",
                    row, col, nrows, ncols
                ),
                path: path_str.clone(),
                line: Some(line_idx + 1),
            });
        }

        triplets.push(Triplet { row, col, val });
        // Symmetrize: add (col, row, val) for off-diagonal entries
        if row != col {
            triplets.push(Triplet {
                row: col,
                col: row,
                val,
            });
        }
    }

    SparseColMat::try_new_from_triplets(nrows, ncols, &triplets).map_err(|e| {
        SparseError::ParseError {
            reason: format!("failed to construct sparse matrix: {:?}", e),
            path: path_str,
            line: None,
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test-data")
    }

    #[test]
    fn load_arrow_5_pd() {
        let path = test_data_dir().join("hand-constructed/arrow-5-pd.mtx");
        let mat = load_mtx(&path).expect("failed to load arrow-5-pd.mtx");
        assert_eq!(mat.nrows(), 5);
        assert_eq!(mat.ncols(), 5);
        // File has 9 entries (lower triangle including diagonal).
        // After symmetrization: 5 diagonal + 4 off-diagonal pairs = 13 stored entries.
        assert_eq!(mat.compute_nnz(), 13);
    }

    #[test]
    fn malformed_input_returns_parse_error() {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test-tmp");
        std::fs::create_dir_all(&dir).ok();
        let path = dir.join("malformed.mtx");
        std::fs::write(&path, "this is not a matrix market file\n").unwrap();
        let result = load_mtx(&path);
        assert!(result.is_err());
        match result.unwrap_err() {
            SparseError::ParseError { .. } => {}
            e => panic!("expected ParseError, got: {:?}", e),
        }
    }

    #[test]
    fn unsupported_format_returns_error() {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test-tmp");
        std::fs::create_dir_all(&dir).ok();
        let path = dir.join("unsupported.mtx");
        std::fs::write(
            &path,
            "%%MatrixMarket matrix array real general\n2 2\n1.0\n2.0\n3.0\n4.0\n",
        )
        .unwrap();
        let result = load_mtx(&path);
        assert!(result.is_err());
        match result.unwrap_err() {
            SparseError::ParseError { .. } => {}
            e => panic!("expected ParseError, got: {:?}", e),
        }
    }

    #[test]
    fn nonexistent_file_returns_io_error() {
        let result = load_mtx(Path::new("/nonexistent/path/matrix.mtx"));
        assert!(result.is_err());
        match result.unwrap_err() {
            SparseError::IoError { .. } => {}
            e => panic!("expected IoError, got: {:?}", e),
        }
    }
}
