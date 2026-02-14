//! Mixed 1x1/2x2 block diagonal storage for the D factor in LDL^T.
//!
//! Provides [`MixedDiagonal`], the data structure representing the D factor in
//! P^T A P = L D L^T where D contains a mix of 1x1 scalar pivots and 2x2
//! symmetric Bunch-Kaufman pivot blocks. Supports incremental construction,
//! query, solve, and inertia computation.
//!
//! Reference: Hogg, Duff & Lopez (2020), "A New Sparse LDL^T Solver Using
//! A Posteriori Threshold Pivoting", SIAM J. Sci. Comput. 42(4), Section 3.

use super::inertia::Inertia;
use super::pivot::{Block2x2, PivotType};

/// The D factor in P^T A P = L D L^T with mixed 1x1 and 2x2 blocks.
///
/// Stores a block diagonal matrix where each block is either a 1x1 scalar
/// pivot or a 2x2 symmetric Bunch-Kaufman pivot block. The structure is built
/// incrementally during factorization (each column starts as [`PivotType::Delayed`]
/// and is set to 1x1 or 2x2 as pivots are decided).
///
/// # Storage layout
///
/// Uses parallel arrays indexed by column for O(1) access during solve:
/// - `pivot_map[col]`: pivot classification per column
/// - `diag[col]`: diagonal value — for 1x1 pivots the scalar value, for 2x2
///   blocks the `a` value at the owner column and `c` at the partner column
/// - `off_diag[col]`: off-diagonal `b` for 2x2 blocks (at the owner column);
///   0.0 at 1x1 and delayed columns
///
/// [`Block2x2`] serves as the API type for input ([`set_2x2`](Self::set_2x2))
/// and output ([`get_2x2`](Self::get_2x2)) but is not the storage format.
///
/// # References
///
/// - Hogg, Duff & Lopez (2020), Section 3: mixed diagonal D storage in APTP
/// - Bunch & Kaufman (1977): 2x2 pivot block structure
#[derive(Debug)]
pub struct MixedDiagonal {
    pivot_map: Vec<PivotType>,
    diag: Vec<f64>,
    off_diag: Vec<f64>,
    n: usize,
}

impl MixedDiagonal {
    /// Create a new `MixedDiagonal` of dimension `n`.
    ///
    /// All columns start as [`PivotType::Delayed`] (unset).
    pub fn new(n: usize) -> Self {
        Self {
            pivot_map: vec![PivotType::Delayed; n],
            diag: vec![0.0; n],
            off_diag: vec![0.0; n],
            n,
        }
    }

    /// Set column `col` as a 1x1 pivot with the given diagonal value.
    ///
    /// # Panics (debug only)
    ///
    /// - `col >= n` (bounds check)
    /// - Column is not currently [`PivotType::Delayed`] (cannot overwrite a set pivot)
    pub fn set_1x1(&mut self, col: usize, value: f64) {
        debug_assert!(
            col < self.n,
            "set_1x1: col {} out of bounds (n = {})",
            col,
            self.n
        );
        debug_assert!(
            self.pivot_map[col] == PivotType::Delayed,
            "set_1x1: col {} is already set ({:?})",
            col,
            self.pivot_map[col]
        );
        self.pivot_map[col] = PivotType::OneByOne;
        self.diag[col] = value;
    }

    /// Set a 2x2 pivot block starting at `block.first_col`.
    ///
    /// Marks both `first_col` and `first_col + 1` as [`PivotType::TwoByTwo`].
    ///
    /// # Panics (debug only)
    ///
    /// - `first_col + 1 >= n` (block must fit)
    /// - Either column is not currently [`PivotType::Delayed`]
    pub fn set_2x2(&mut self, block: Block2x2) {
        let col = block.first_col;
        debug_assert!(
            col + 1 < self.n,
            "set_2x2: first_col {} + 1 out of bounds (n = {})",
            col,
            self.n
        );
        debug_assert!(
            self.pivot_map[col] == PivotType::Delayed,
            "set_2x2: col {} is already set ({:?})",
            col,
            self.pivot_map[col]
        );
        debug_assert!(
            self.pivot_map[col + 1] == PivotType::Delayed,
            "set_2x2: col {} is already set ({:?})",
            col + 1,
            self.pivot_map[col + 1]
        );
        self.pivot_map[col] = PivotType::TwoByTwo { partner: col + 1 };
        self.pivot_map[col + 1] = PivotType::TwoByTwo { partner: col };
        self.diag[col] = block.a;
        self.diag[col + 1] = block.c;
        self.off_diag[col] = block.b;
    }

    /// Matrix dimension.
    pub fn dimension(&self) -> usize {
        self.n
    }

    /// Pivot type for the given column.
    ///
    /// # Panics (debug only)
    ///
    /// `col >= n`
    pub fn get_pivot_type(&self, col: usize) -> PivotType {
        debug_assert!(
            col < self.n,
            "get_pivot_type: col {} out of bounds (n = {})",
            col,
            self.n
        );
        self.pivot_map[col]
    }

    /// Diagonal value for a 1x1 pivot column.
    ///
    /// # Panics (debug only)
    ///
    /// Column is not [`PivotType::OneByOne`].
    pub fn get_1x1(&self, col: usize) -> f64 {
        debug_assert!(
            self.pivot_map[col] == PivotType::OneByOne,
            "get_1x1: col {} is not OneByOne ({:?})",
            col,
            self.pivot_map[col]
        );
        self.diag[col]
    }

    /// Block data for a 2x2 pivot (by the lower-indexed column).
    ///
    /// Constructs a [`Block2x2`] from the inline parallel arrays. This is O(1).
    ///
    /// # Panics (debug only)
    ///
    /// Column is not the owner (lower-indexed) of a [`PivotType::TwoByTwo`] block.
    pub fn get_2x2(&self, first_col: usize) -> Block2x2 {
        debug_assert!(
            matches!(self.pivot_map[first_col], PivotType::TwoByTwo { partner } if partner > first_col),
            "get_2x2: col {} is not a 2x2 block owner ({:?})",
            first_col,
            self.pivot_map[first_col]
        );
        Block2x2 {
            first_col,
            a: self.diag[first_col],
            b: self.off_diag[first_col],
            c: self.diag[first_col + 1],
        }
    }

    /// Number of columns still marked as [`PivotType::Delayed`].
    pub fn num_delayed(&self) -> usize {
        self.pivot_map
            .iter()
            .filter(|p| **p == PivotType::Delayed)
            .count()
    }

    /// Number of 1x1 pivots.
    pub fn num_1x1(&self) -> usize {
        self.pivot_map
            .iter()
            .filter(|p| **p == PivotType::OneByOne)
            .count()
    }

    /// Number of 2x2 pivot pairs.
    pub fn num_2x2_pairs(&self) -> usize {
        self.pivot_map
            .iter()
            .enumerate()
            .filter(|(i, p)| matches!(p, PivotType::TwoByTwo { partner } if *partner > *i))
            .count()
    }

    /// Solve D x = b in place, where `x` initially contains the right-hand side b.
    ///
    /// For 1x1 pivots: `x[i] /= d[i]`.
    /// For 2x2 blocks: solves `[[a, b], [b, c]] * [x1, x2]^T = [r1, r2]^T`
    /// via Cramer's rule (analytical 2x2 inverse).
    ///
    /// # Panics (debug only)
    ///
    /// - Any column is still [`PivotType::Delayed`]
    /// - Any 1x1 pivot value is zero
    /// - Any 2x2 block has zero determinant
    /// - `x.len() != n`
    ///
    /// # References
    ///
    /// - Cramer's rule for 2x2 symmetric systems
    /// - Hogg, Duff & Lopez (2020), Section 3: D-solve in APTP context
    pub fn solve_in_place(&self, x: &mut [f64]) {
        debug_assert_eq!(
            x.len(),
            self.n,
            "solve_in_place: x.len() = {} != n = {}",
            x.len(),
            self.n
        );
        debug_assert!(
            self.num_delayed() == 0,
            "solve_in_place: {} delayed columns remain",
            self.num_delayed()
        );

        let mut col = 0;
        while col < self.n {
            match self.pivot_map[col] {
                PivotType::OneByOne => {
                    let d = self.diag[col];
                    debug_assert!(d != 0.0, "solve_in_place: zero 1x1 pivot at col {}", col);
                    x[col] /= d;
                    col += 1;
                }
                PivotType::TwoByTwo { partner } => {
                    if partner > col {
                        // This is the owner (lower-indexed column) of the 2x2 block.
                        // Read a, b, c directly from parallel arrays — O(1).
                        let a = self.diag[col];
                        let b = self.off_diag[col];
                        let c = self.diag[partner];
                        let det = a * c - b * b;
                        debug_assert!(
                            det != 0.0,
                            "solve_in_place: zero determinant 2x2 block at col {}",
                            col
                        );
                        let r1 = x[col];
                        let r2 = x[partner];
                        // Cramer's rule: [[a,b],[b,c]]^-1 = (1/det) * [[c,-b],[-b,a]]
                        x[col] = (c * r1 - b * r2) / det;
                        x[partner] = (a * r2 - b * r1) / det;
                    }
                    // Skip partner column if we've already processed this pair
                    col += 1;
                }
                PivotType::Delayed => {
                    unreachable!("solve_in_place: delayed column at {}", col);
                }
            }
        }
    }

    /// Compute eigenvalue sign counts from stored pivots.
    ///
    /// For 1x1 pivots, the sign of the diagonal value determines the eigenvalue
    /// sign. For 2x2 blocks, the trace and determinant classify the eigenvalue
    /// signs without computing actual eigenvalues:
    ///
    /// | Condition | Eigenvalue signs |
    /// |-----------|-----------------|
    /// | det > 0, trace > 0 | both positive |
    /// | det > 0, trace < 0 | both negative |
    /// | det < 0 | one positive, one negative |
    /// | det = 0, trace > 0 | one positive, one zero |
    /// | det = 0, trace < 0 | one negative, one zero |
    /// | det = 0, trace = 0 | both zero |
    ///
    /// # Panics (debug only)
    ///
    /// Any column is still [`PivotType::Delayed`].
    ///
    /// # References
    ///
    /// - Standard eigenvalue sign classification from trace/determinant
    /// - Hogg, Duff & Lopez (2020), Section 2: inertia in APTP context
    /// - Bunch & Kaufman (1977): inertia from pivot classifications
    pub fn compute_inertia(&self) -> Inertia {
        debug_assert!(
            self.num_delayed() == 0,
            "compute_inertia: {} delayed columns remain",
            self.num_delayed()
        );

        let mut positive = 0usize;
        let mut negative = 0usize;
        let mut zero = 0usize;

        let mut col = 0;
        while col < self.n {
            match self.pivot_map[col] {
                PivotType::OneByOne => {
                    let d = self.diag[col];
                    if d > 0.0 {
                        positive += 1;
                    } else if d < 0.0 {
                        negative += 1;
                    } else {
                        zero += 1;
                    }
                    col += 1;
                }
                PivotType::TwoByTwo { partner } => {
                    if partner > col {
                        // Read a, b, c directly from parallel arrays — O(1).
                        let a = self.diag[col];
                        let b = self.off_diag[col];
                        let c = self.diag[partner];
                        let det = a * c - b * b;
                        let trace = a + c;

                        if det > 0.0 {
                            if trace > 0.0 {
                                positive += 2;
                            } else {
                                // trace < 0 (trace == 0 impossible with det > 0 for real symmetric)
                                negative += 2;
                            }
                        } else if det < 0.0 {
                            positive += 1;
                            negative += 1;
                        } else {
                            // det == 0
                            if trace > 0.0 {
                                positive += 1;
                                zero += 1;
                            } else if trace < 0.0 {
                                negative += 1;
                                zero += 1;
                            } else {
                                zero += 2;
                            }
                        }
                    }
                    col += 1;
                }
                PivotType::Delayed => {
                    unreachable!("compute_inertia: delayed column at {}", col);
                }
            }
        }

        Inertia {
            positive,
            negative,
            zero,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aptp::pivot::{Block2x2, PivotType};

    // ---- T009: MixedDiagonal construction and query ----

    #[test]
    fn new_creates_all_delayed() {
        let diag = MixedDiagonal::new(5);
        assert_eq!(diag.dimension(), 5);
        for col in 0..5 {
            assert_eq!(diag.get_pivot_type(col), PivotType::Delayed);
        }
        assert_eq!(diag.num_delayed(), 5);
        assert_eq!(diag.num_1x1(), 0);
        assert_eq!(diag.num_2x2_pairs(), 0);
    }

    #[test]
    fn set_1x1_marks_correct_pivot_type() {
        let mut diag = MixedDiagonal::new(4);
        diag.set_1x1(0, 3.5);
        diag.set_1x1(2, -1.0);

        assert_eq!(diag.get_pivot_type(0), PivotType::OneByOne);
        assert_eq!(diag.get_pivot_type(1), PivotType::Delayed);
        assert_eq!(diag.get_pivot_type(2), PivotType::OneByOne);
        assert_eq!(diag.get_pivot_type(3), PivotType::Delayed);

        assert_eq!(diag.get_1x1(0), 3.5);
        assert_eq!(diag.get_1x1(2), -1.0);

        assert_eq!(diag.num_1x1(), 2);
        assert_eq!(diag.num_delayed(), 2);
    }

    #[test]
    fn set_2x2_marks_both_columns() {
        let mut diag = MixedDiagonal::new(6);
        let block = Block2x2 {
            first_col: 2,
            a: 2.0,
            b: 0.5,
            c: -3.0,
        };
        diag.set_2x2(block);

        assert_eq!(diag.get_pivot_type(2), PivotType::TwoByTwo { partner: 3 });
        assert_eq!(diag.get_pivot_type(3), PivotType::TwoByTwo { partner: 2 });
        assert_eq!(diag.get_2x2(2), block);
        assert_eq!(diag.num_2x2_pairs(), 1);
        assert_eq!(diag.num_delayed(), 4);
    }

    #[test]
    fn mixed_pivots_correct_counts() {
        let mut diag = MixedDiagonal::new(6);
        diag.set_2x2(Block2x2 {
            first_col: 0,
            a: 2.0,
            b: 0.5,
            c: -3.0,
        });
        diag.set_1x1(2, 4.0);
        diag.set_1x1(3, -1.0);
        diag.set_1x1(4, 7.0);
        diag.set_1x1(5, 2.0);

        assert_eq!(diag.num_2x2_pairs(), 1);
        assert_eq!(diag.num_1x1(), 4);
        assert_eq!(diag.num_delayed(), 0);
        assert_eq!(diag.dimension(), 6);
    }

    // ---- T010: solve_in_place ----

    #[test]
    fn solve_all_1x1() {
        // D = diag(2, 4, -1, 5)
        // b = [6, 12, -3, 20]
        // x = [3, 3, 3, 4]
        let mut diag = MixedDiagonal::new(4);
        diag.set_1x1(0, 2.0);
        diag.set_1x1(1, 4.0);
        diag.set_1x1(2, -1.0);
        diag.set_1x1(3, 5.0);

        let mut x = vec![6.0, 12.0, -3.0, 20.0];
        let b = x.clone();
        diag.solve_in_place(&mut x);

        assert_eq!(x, vec![3.0, 3.0, 3.0, 4.0]);

        // Verify D*x = b (relative error)
        let dx: Vec<f64> = vec![2.0 * x[0], 4.0 * x[1], -x[2], 5.0 * x[3]];
        let norm_b: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        let norm_diff: f64 = dx
            .iter()
            .zip(b.iter())
            .map(|(d, bi)| (d - bi).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(norm_diff / norm_b < 1e-14);
    }

    #[test]
    fn solve_all_2x2() {
        // D = [[2, 0.5], [0.5, -3]]  (one 2x2 block)
        // b = [4.5, -0.5] → det = 2*(-3) - 0.5^2 = -6.25
        // x1 = ((-3)*4.5 - 0.5*(-0.5)) / (-6.25) = (-13.5 + 0.25) / (-6.25) = -13.25/-6.25 = 2.12
        // x2 = (2*(-0.5) - 0.5*4.5) / (-6.25) = (-1 - 2.25) / (-6.25) = -3.25/-6.25 = 0.52
        let mut diag = MixedDiagonal::new(2);
        diag.set_2x2(Block2x2 {
            first_col: 0,
            a: 2.0,
            b: 0.5,
            c: -3.0,
        });

        let b = vec![4.5, -0.5];
        let mut x = b.clone();
        diag.solve_in_place(&mut x);

        // Verify: D*x should equal b
        let dx0 = 2.0 * x[0] + 0.5 * x[1];
        let dx1 = 0.5 * x[0] + (-3.0) * x[1];
        let norm_b: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        let norm_diff = ((dx0 - b[0]).powi(2) + (dx1 - b[1]).powi(2)).sqrt();
        assert!(
            norm_diff / norm_b < 1e-14,
            "relative error: {:.2e}",
            norm_diff / norm_b
        );
    }

    #[test]
    fn solve_mixed_1x1_and_2x2() {
        // Dimension 6: 2x2 block at [0,1], then 1x1 at [2,3,4,5]
        let mut diag = MixedDiagonal::new(6);
        diag.set_2x2(Block2x2 {
            first_col: 0,
            a: 2.0,
            b: 0.5,
            c: -3.0,
        });
        diag.set_1x1(2, 4.0);
        diag.set_1x1(3, -1.0);
        diag.set_1x1(4, 7.0);
        diag.set_1x1(5, 2.0);

        let b = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut x = b.clone();
        diag.solve_in_place(&mut x);

        // Verify: reconstruct D*x and compare to b
        // 2x2 block: D*x at [0,1]
        let dx0 = 2.0 * x[0] + 0.5 * x[1];
        let dx1 = 0.5 * x[0] + (-3.0) * x[1];
        // 1x1 blocks: D*x at [2,3,4,5]
        let dx2 = 4.0 * x[2];
        let dx3 = -x[3];
        let dx4 = 7.0 * x[4];
        let dx5 = 2.0 * x[5];

        let dx = [dx0, dx1, dx2, dx3, dx4, dx5];
        let norm_b: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        let norm_diff: f64 = dx
            .iter()
            .zip(b.iter())
            .map(|(d, bi)| (d - bi).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(
            norm_diff / norm_b < 1e-14,
            "relative error: {:.2e}",
            norm_diff / norm_b
        );
    }

    #[test]
    fn solve_dimension_0_is_noop() {
        let diag = MixedDiagonal::new(0);
        let mut x: Vec<f64> = vec![];
        diag.solve_in_place(&mut x);
        assert!(x.is_empty());
    }

    // ---- T011: Edge case tests ----

    #[test]
    fn dimension_0() {
        let diag = MixedDiagonal::new(0);
        assert_eq!(diag.dimension(), 0);
        assert_eq!(diag.num_delayed(), 0);
        assert_eq!(diag.num_1x1(), 0);
        assert_eq!(diag.num_2x2_pairs(), 0);
    }

    #[test]
    fn dimension_1_single_1x1() {
        let mut diag = MixedDiagonal::new(1);
        diag.set_1x1(0, 5.0);
        assert_eq!(diag.get_pivot_type(0), PivotType::OneByOne);
        assert_eq!(diag.get_1x1(0), 5.0);
        assert_eq!(diag.num_1x1(), 1);
        assert_eq!(diag.num_delayed(), 0);
    }

    #[test]
    fn dimension_2_single_2x2() {
        let mut diag = MixedDiagonal::new(2);
        let block = Block2x2 {
            first_col: 0,
            a: 1.0,
            b: 0.0,
            c: 1.0,
        };
        diag.set_2x2(block);
        assert_eq!(diag.num_2x2_pairs(), 1);
        assert_eq!(diag.num_delayed(), 0);
    }

    #[test]
    fn all_2x2_even_n() {
        let mut diag = MixedDiagonal::new(4);
        diag.set_2x2(Block2x2 {
            first_col: 0,
            a: 1.0,
            b: 0.0,
            c: 1.0,
        });
        diag.set_2x2(Block2x2 {
            first_col: 2,
            a: 2.0,
            b: 0.5,
            c: 3.0,
        });
        assert_eq!(diag.num_2x2_pairs(), 2);
        assert_eq!(diag.num_delayed(), 0);
    }

    #[test]
    #[should_panic]
    fn solve_panics_on_delayed_columns() {
        let mut diag = MixedDiagonal::new(3);
        diag.set_1x1(0, 1.0);
        // columns 1, 2 still delayed
        let mut x = vec![1.0, 2.0, 3.0];
        diag.solve_in_place(&mut x); // should panic in debug mode
    }

    #[test]
    #[should_panic]
    fn set_2x2_at_last_column_odd_n_panics() {
        let mut diag = MixedDiagonal::new(3);
        // first_col = 2, but 2+1 = 3 which is NOT < 3, so debug-assert should fire
        diag.set_2x2(Block2x2 {
            first_col: 2,
            a: 1.0,
            b: 0.0,
            c: 1.0,
        });
    }

    // ---- T017: compute_inertia unit tests ----

    #[test]
    fn inertia_all_positive_1x1() {
        let mut diag = MixedDiagonal::new(4);
        for i in 0..4 {
            diag.set_1x1(i, (i + 1) as f64);
        }
        let inertia = diag.compute_inertia();
        assert_eq!(
            inertia,
            Inertia {
                positive: 4,
                negative: 0,
                zero: 0
            }
        );
    }

    #[test]
    fn inertia_mixed_sign_1x1() {
        let mut diag = MixedDiagonal::new(5);
        diag.set_1x1(0, 3.0); // +
        diag.set_1x1(1, -2.0); // -
        diag.set_1x1(2, 1.0); // +
        diag.set_1x1(3, -0.5); // -
        diag.set_1x1(4, 0.0); // zero
        let inertia = diag.compute_inertia();
        assert_eq!(
            inertia,
            Inertia {
                positive: 2,
                negative: 2,
                zero: 1
            }
        );
    }

    #[test]
    fn inertia_2x2_det_negative_one_plus_one_minus() {
        // [[2, 0.5], [0.5, -3]] → det = -6.25 < 0 → one +, one -
        let mut diag = MixedDiagonal::new(2);
        diag.set_2x2(Block2x2 {
            first_col: 0,
            a: 2.0,
            b: 0.5,
            c: -3.0,
        });
        let inertia = diag.compute_inertia();
        assert_eq!(
            inertia,
            Inertia {
                positive: 1,
                negative: 1,
                zero: 0
            }
        );
    }

    #[test]
    fn inertia_2x2_det_positive_trace_positive() {
        // [[5, 1], [1, 3]] → det = 15-1 = 14 > 0, trace = 8 > 0 → both positive
        let mut diag = MixedDiagonal::new(2);
        diag.set_2x2(Block2x2 {
            first_col: 0,
            a: 5.0,
            b: 1.0,
            c: 3.0,
        });
        let inertia = diag.compute_inertia();
        assert_eq!(
            inertia,
            Inertia {
                positive: 2,
                negative: 0,
                zero: 0
            }
        );
    }

    #[test]
    fn inertia_2x2_det_positive_trace_negative() {
        // [[-5, 1], [1, -3]] → det = 15-1 = 14 > 0, trace = -8 < 0 → both negative
        let mut diag = MixedDiagonal::new(2);
        diag.set_2x2(Block2x2 {
            first_col: 0,
            a: -5.0,
            b: 1.0,
            c: -3.0,
        });
        let inertia = diag.compute_inertia();
        assert_eq!(
            inertia,
            Inertia {
                positive: 0,
                negative: 2,
                zero: 0
            }
        );
    }

    #[test]
    fn inertia_mixed_1x1_and_2x2() {
        // quickstart.md example:
        // 2x2 [[2, 0.5], [0.5, -3]] → det < 0 → 1+, 1-
        // 1x1: 4(+), -1(-), 7(+), 2(+)
        // Total: 4+, 2-
        let mut diag = MixedDiagonal::new(6);
        diag.set_2x2(Block2x2 {
            first_col: 0,
            a: 2.0,
            b: 0.5,
            c: -3.0,
        });
        diag.set_1x1(2, 4.0);
        diag.set_1x1(3, -1.0);
        diag.set_1x1(4, 7.0);
        diag.set_1x1(5, 2.0);
        let inertia = diag.compute_inertia();
        assert_eq!(
            inertia,
            Inertia {
                positive: 4,
                negative: 2,
                zero: 0
            }
        );
    }

    #[test]
    fn scale_test_n_10000() {
        // SC-001: MixedDiagonal at n=10,000 with random mixed pivot pattern
        let n = 10_000;
        let mut diag = MixedDiagonal::new(n);

        // Alternate: 2x2 block, then 1x1, then 2x2, then 1x1...
        // Pattern: pairs at (0,1), 1x1 at 2, pairs at (3,4), 1x1 at 5, ...
        let mut col = 0;
        while col < n {
            if col + 1 < n && col % 3 != 2 {
                diag.set_2x2(Block2x2 {
                    first_col: col,
                    a: 2.0 + (col as f64) * 0.001,
                    b: 0.1,
                    c: 3.0 + (col as f64) * 0.001,
                });
                col += 2;
            } else {
                diag.set_1x1(col, 1.0 + (col as f64) * 0.001);
                col += 1;
            }
        }

        assert_eq!(diag.num_delayed(), 0);
        assert_eq!(diag.dimension(), n);

        // Solve round-trip: set b = [1, 2, 3, ...], solve D*x = b, verify D*x ≈ b
        let b: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        let mut x = b.clone();
        diag.solve_in_place(&mut x);

        // Reconstruct D*x
        let mut dx = vec![0.0; n];
        for i in 0..n {
            match diag.get_pivot_type(i) {
                PivotType::OneByOne => {
                    dx[i] = diag.get_1x1(i) * x[i];
                }
                PivotType::TwoByTwo { partner } => {
                    if i < partner {
                        let block = diag.get_2x2(i);
                        dx[i] = block.a * x[i] + block.b * x[partner];
                        dx[partner] = block.b * x[i] + block.c * x[partner];
                    }
                    // Skip partner column (already handled)
                }
                PivotType::Delayed => unreachable!(),
            }
        }

        let norm_b: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        let norm_diff: f64 = dx
            .iter()
            .zip(b.iter())
            .map(|(d, bi)| (d - bi).powi(2))
            .sum::<f64>()
            .sqrt();
        let rel_err = norm_diff / norm_b;
        assert!(
            rel_err < 1e-14,
            "SC-001 scale test: relative error {:.2e} exceeds 1e-14",
            rel_err
        );
    }
}
