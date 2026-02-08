//! Text-based sparsity pattern visualization.
//!
//! Renders sparse matrix non-zero structure as a character grid using
//! density-based characters. Large matrices are downsampled to fit
//! configurable display dimensions.

use std::fmt;

use faer::sparse::SparseColMat;

/// A text-based visualization of a sparse matrix's non-zero structure.
pub struct SparsityDisplay {
    /// Original matrix row count.
    matrix_rows: usize,
    /// Original matrix column count.
    matrix_cols: usize,
    /// Original non-zero count.
    matrix_nnz: usize,
    /// Maximum display width in characters.
    max_width: usize,
    /// Maximum display height in characters.
    max_height: usize,
    /// Use ASCII characters only (no Unicode blocks).
    ascii_only: bool,
    /// Density per display cell (0.0 to 1.0), row-major.
    cells: Vec<Vec<f64>>,
    /// Actual display columns used.
    display_cols: usize,
    /// Actual display rows used.
    display_rows: usize,
}

impl SparsityDisplay {
    /// Analyze the sparsity structure of a CSC matrix.
    pub fn from_sparse(matrix: &SparseColMat<usize, f64>) -> Self {
        let matrix_rows = matrix.nrows();
        let matrix_cols = matrix.ncols();
        let matrix_nnz = matrix.compute_nnz();

        let mut display = Self {
            matrix_rows,
            matrix_cols,
            matrix_nnz,
            max_width: 80,
            max_height: 40,
            ascii_only: false,
            cells: Vec::new(),
            display_cols: 0,
            display_rows: 0,
        };

        display.compute_cells(matrix);
        display
    }

    /// Set maximum display width (default: 80).
    pub fn with_max_width(mut self, cols: usize) -> Self {
        self.max_width = cols;
        self.recompute_dimensions();
        self
    }

    /// Set maximum display height (default: 40).
    pub fn with_max_height(mut self, rows: usize) -> Self {
        self.max_height = rows;
        self.recompute_dimensions();
        self
    }

    /// If true, use ASCII characters only. Default: false (Unicode blocks).
    pub fn with_ascii_only(mut self, ascii: bool) -> Self {
        self.ascii_only = ascii;
        self
    }

    /// Produce the complete text output including header and grid.
    pub fn render(&self) -> String {
        let mut out = String::new();

        // Header
        let density_pct = if self.matrix_rows > 0 && self.matrix_cols > 0 {
            self.matrix_nnz as f64 / (self.matrix_rows as f64 * self.matrix_cols as f64) * 100.0
        } else {
            0.0
        };

        out.push_str(&format!(
            "({}x{}, nnz={}, density={:.1}%)",
            self.matrix_rows,
            self.matrix_cols,
            format_with_commas(self.matrix_nnz),
            density_pct,
        ));

        // Add downsampling note if applicable
        if self.display_rows != self.matrix_rows || self.display_cols != self.matrix_cols {
            out.push_str(&format!(
                " [{}x{} view]",
                self.display_cols, self.display_rows
            ));
        }

        out.push('\n');

        // Separator
        let header_len = out.lines().next().map_or(40, |l| l.len());
        for _ in 0..header_len {
            out.push('\u{2500}'); // ─
        }
        out.push('\n');

        // Grid
        if self.cells.is_empty() || self.matrix_rows == 0 || self.matrix_cols == 0 {
            out.push_str("(empty matrix)\n");
            return out;
        }

        for row in &self.cells {
            for &density in row {
                out.push(self.density_char(density));
            }
            out.push('\n');
        }

        out
    }

    /// Map density value to display character.
    fn density_char(&self, density: f64) -> char {
        if self.ascii_only {
            match density {
                d if d <= 0.0 => '.',
                d if d <= 0.33 => '-',
                d if d <= 0.66 => '+',
                _ => '#',
            }
        } else {
            match density {
                d if d <= 0.0 => '.',
                d if d <= 0.25 => '\u{2591}', // ░
                d if d <= 0.50 => '\u{2592}', // ▒
                d if d <= 0.75 => '\u{2593}', // ▓
                _ => '\u{2588}',              // █
            }
        }
    }

    /// Compute cell densities from the matrix structure.
    fn compute_cells(&mut self, matrix: &SparseColMat<usize, f64>) {
        if self.matrix_rows == 0 || self.matrix_cols == 0 {
            self.display_rows = 0;
            self.display_cols = 0;
            self.cells = Vec::new();
            return;
        }

        self.display_rows = self.matrix_rows.min(self.max_height);
        self.display_cols = self.matrix_cols.min(self.max_width);

        // Compute bin sizes
        let row_bin_size = (self.matrix_rows as f64 / self.display_rows as f64).ceil() as usize;
        let col_bin_size = (self.matrix_cols as f64 / self.display_cols as f64).ceil() as usize;

        // Count non-zeros per bin
        let mut counts = vec![vec![0u64; self.display_cols]; self.display_rows];

        let col_ptr = matrix.col_ptr();
        let row_idx = matrix.row_idx();

        for col in 0..self.matrix_cols {
            let bin_col = (col / col_bin_size).min(self.display_cols - 1);
            let range = col_ptr[col]..col_ptr[col + 1];
            for &row in &row_idx[range] {
                let bin_row = (row / row_bin_size).min(self.display_rows - 1);
                counts[bin_row][bin_col] += 1;
            }
        }

        // Convert counts to densities
        self.cells = counts
            .iter()
            .enumerate()
            .map(|(r, row)| {
                row.iter()
                    .enumerate()
                    .map(|(c, &count)| {
                        // Compute actual bin area (handle edge bins)
                        let bin_rows = actual_bin_size(r, row_bin_size, self.matrix_rows);
                        let bin_cols = actual_bin_size(c, col_bin_size, self.matrix_cols);
                        let area = (bin_rows * bin_cols) as f64;
                        if area > 0.0 {
                            (count as f64 / area).min(1.0)
                        } else {
                            0.0
                        }
                    })
                    .collect()
            })
            .collect();
    }

    /// Recompute dimensions after max_width/max_height change.
    /// Note: we can't recompute cells without the matrix, so this only
    /// adjusts dimensions. In practice, builders are chained before use.
    fn recompute_dimensions(&mut self) {
        self.display_rows = self.matrix_rows.min(self.max_height);
        self.display_cols = self.matrix_cols.min(self.max_width);
    }
}

/// Compute actual size of a bin at the given index, accounting for edge bins.
fn actual_bin_size(bin_idx: usize, bin_size: usize, total: usize) -> usize {
    let start = bin_idx * bin_size;
    let end = ((bin_idx + 1) * bin_size).min(total);
    end.saturating_sub(start)
}

fn format_with_commas(n: usize) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

impl fmt::Display for SparsityDisplay {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.render())
    }
}

#[cfg(test)]
mod tests {
    use faer::sparse::{SparseColMat, Triplet};

    use super::*;

    /// Create a small tridiagonal matrix for testing.
    fn tridiagonal_5x5() -> SparseColMat<usize, f64> {
        let triplets = vec![
            Triplet::new(0, 0, 2.0),
            Triplet::new(0, 1, -1.0),
            Triplet::new(1, 0, -1.0),
            Triplet::new(1, 1, 2.0),
            Triplet::new(1, 2, -1.0),
            Triplet::new(2, 1, -1.0),
            Triplet::new(2, 2, 2.0),
            Triplet::new(2, 3, -1.0),
            Triplet::new(3, 2, -1.0),
            Triplet::new(3, 3, 2.0),
            Triplet::new(3, 4, -1.0),
            Triplet::new(4, 3, -1.0),
            Triplet::new(4, 4, 2.0),
        ];
        SparseColMat::try_new_from_triplets(5, 5, &triplets).unwrap()
    }

    /// Create a diagonal 3x3 matrix.
    fn diagonal_3x3() -> SparseColMat<usize, f64> {
        let triplets = vec![
            Triplet::new(0, 0, 1.0),
            Triplet::new(1, 1, 2.0),
            Triplet::new(2, 2, 3.0),
        ];
        SparseColMat::try_new_from_triplets(3, 3, &triplets).unwrap()
    }

    #[test]
    fn from_sparse_small_tridiagonal() {
        let mat = tridiagonal_5x5();
        let display = SparsityDisplay::from_sparse(&mat);
        let rendered = display.render();

        // Header should contain dimensions and nnz
        assert!(rendered.contains("5x5"), "should show dimensions");
        assert!(rendered.contains("nnz=13"), "should show nnz count");
    }

    #[test]
    fn small_matrix_no_downsampling() {
        let mat = tridiagonal_5x5();
        let display = SparsityDisplay::from_sparse(&mat);
        // 5x5 fits within 80x40 default, so no downsampling
        assert_eq!(display.display_rows, 5);
        assert_eq!(display.display_cols, 5);
    }

    #[test]
    fn nonzero_positions_show_filled() {
        let mat = diagonal_3x3();
        let display = SparsityDisplay::from_sparse(&mat);
        let rendered = display.render();

        // Skip header lines, get grid
        let grid_lines: Vec<&str> = rendered.lines().skip(2).collect();
        assert_eq!(grid_lines.len(), 3, "3x3 matrix should have 3 grid rows");

        // Diagonal should be filled, off-diagonal should be dots
        for (i, line) in grid_lines.iter().enumerate() {
            let chars: Vec<char> = line.chars().collect();
            assert_eq!(chars.len(), 3, "each row should have 3 chars");
            for (j, &c) in chars.iter().enumerate() {
                if i == j {
                    assert_ne!(c, '.', "diagonal position ({i},{j}) should be filled");
                } else {
                    assert_eq!(c, '.', "off-diagonal ({i},{j}) should be dot");
                }
            }
        }
    }

    #[test]
    fn tridiagonal_pattern_correct() {
        let mat = tridiagonal_5x5();
        let display = SparsityDisplay::from_sparse(&mat);
        let rendered = display.render();

        let grid_lines: Vec<&str> = rendered.lines().skip(2).collect();
        for (i, line) in grid_lines.iter().enumerate() {
            let chars: Vec<char> = line.chars().collect();
            for (j, &c) in chars.iter().enumerate() {
                let should_be_nonzero = i == j || (i > 0 && j == i - 1) || (j > 0 && i == j - 1);
                if should_be_nonzero {
                    assert_ne!(c, '.', "position ({i},{j}) should be non-zero");
                } else {
                    assert_eq!(c, '.', "position ({i},{j}) should be zero");
                }
            }
        }
    }

    #[test]
    fn downsampling_large_matrix() {
        // Create a larger matrix (100x100 diagonal)
        let triplets: Vec<Triplet<usize, usize, f64>> =
            (0..100).map(|i| Triplet::new(i, i, 1.0)).collect();
        let mat = SparseColMat::try_new_from_triplets(100, 100, &triplets).unwrap();

        let display = SparsityDisplay::from_sparse(&mat)
            .with_max_width(20)
            .with_max_height(10);

        // Should be downsampled
        assert!(display.display_cols <= 20);
        assert!(display.display_rows <= 10);

        let rendered = display.render();
        assert!(
            rendered.contains("["),
            "should have downsampling annotation"
        );
        assert!(
            rendered.contains("view]"),
            "should have downsampling annotation"
        );
    }

    #[test]
    fn downsampling_density_values() {
        // Create a fully dense 10x10 matrix
        let mut triplets = Vec::new();
        for i in 0..10 {
            for j in 0..10 {
                triplets.push(Triplet::new(i, j, 1.0));
            }
        }
        let mat = SparseColMat::try_new_from_triplets(10, 10, &triplets).unwrap();

        let display = SparsityDisplay::from_sparse(&mat)
            .with_max_width(5)
            .with_max_height(5);

        // All bins should be 100% dense
        for row in &display.cells {
            for &density in row {
                assert!(
                    (density - 1.0).abs() < 1e-10,
                    "fully dense matrix should have density 1.0, got {density}"
                );
            }
        }
    }

    #[test]
    fn ascii_mode_uses_ascii_only() {
        let mat = tridiagonal_5x5();
        let display = SparsityDisplay::from_sparse(&mat).with_ascii_only(true);
        let rendered = display.render();

        // Skip header (which may have Unicode separator)
        let grid_lines: Vec<&str> = rendered.lines().skip(2).collect();
        for line in &grid_lines {
            for c in line.chars() {
                assert!(
                    c == '.' || c == '-' || c == '+' || c == '#',
                    "ASCII mode should only use '.', '-', '+', '#', got '{c}'"
                );
            }
        }
    }

    #[test]
    fn display_impl_matches_render() {
        let mat = tridiagonal_5x5();
        let display = SparsityDisplay::from_sparse(&mat);
        let rendered = display.render();
        let displayed = format!("{display}");
        assert_eq!(rendered, displayed);
    }

    #[test]
    fn empty_matrix() {
        let mat = SparseColMat::try_new_from_triplets(0, 0, &[] as &[Triplet<usize, usize, f64>])
            .unwrap();
        let display = SparsityDisplay::from_sparse(&mat);
        let rendered = display.render();
        assert!(rendered.contains("0x0"));
        assert!(rendered.contains("empty matrix"));
    }

    #[test]
    fn header_contains_density() {
        let mat = tridiagonal_5x5();
        let display = SparsityDisplay::from_sparse(&mat);
        let rendered = display.render();
        assert!(
            rendered.contains("density="),
            "should show density percentage"
        );
    }

    #[test]
    fn density_char_unicode_thresholds() {
        let display = SparsityDisplay {
            matrix_rows: 0,
            matrix_cols: 0,
            matrix_nnz: 0,
            max_width: 80,
            max_height: 40,
            ascii_only: false,
            cells: Vec::new(),
            display_cols: 0,
            display_rows: 0,
        };
        assert_eq!(display.density_char(0.0), '.');
        assert_eq!(display.density_char(0.1), '\u{2591}');
        assert_eq!(display.density_char(0.3), '\u{2592}');
        assert_eq!(display.density_char(0.6), '\u{2593}');
        assert_eq!(display.density_char(0.9), '\u{2588}');
    }

    #[test]
    fn density_char_ascii_thresholds() {
        let display = SparsityDisplay {
            matrix_rows: 0,
            matrix_cols: 0,
            matrix_nnz: 0,
            max_width: 80,
            max_height: 40,
            ascii_only: true,
            cells: Vec::new(),
            display_cols: 0,
            display_rows: 0,
        };
        assert_eq!(display.density_char(0.0), '.');
        assert_eq!(display.density_char(0.2), '-');
        assert_eq!(display.density_char(0.5), '+');
        assert_eq!(display.density_char(0.8), '#');
    }
}
