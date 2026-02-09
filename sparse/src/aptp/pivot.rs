//! Pivot classification types for APTP factorization.
//!
//! Provides [`PivotType`] for classifying column pivot decisions and [`Block2x2`]
//! for storing 2x2 symmetric diagonal blocks. These types represent the pivot
//! decisions made during the A Posteriori Threshold Pivoting algorithm described
//! in Hogg, Duff & Lopez (2020), "A New Sparse LDL^T Solver Using A Posteriori
//! Threshold Pivoting", SIAM J. Sci. Comput. 42(4).

/// Classification of a single column's pivot decision during APTP factorization.
///
/// Each column in the elimination is classified as one of:
/// - [`OneByOne`](PivotType::OneByOne): standard 1x1 scalar pivot
/// - [`TwoByTwo`](PivotType::TwoByTwo): 2x2 Bunch-Kaufman pivot with an adjacent partner column
/// - [`Delayed`](PivotType::Delayed): column failed the APTP stability check and is deferred
///   to an ancestor node in the elimination tree
///
/// # References
///
/// - Hogg, Duff & Lopez (2020), Section 3: pivot classification in APTP
/// - Bunch & Kaufman (1977): 2x2 pivot strategy for symmetric indefinite systems
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PivotType {
    /// Standard 1x1 scalar pivot.
    OneByOne,
    /// 2x2 Bunch-Kaufman pivot. `partner` identifies the paired column.
    TwoByTwo { partner: usize },
    /// Column failed APTP stability check; deferred to ancestor node.
    Delayed,
}

/// Storage for a single 2x2 symmetric diagonal block `[[a, b], [b, c]]`.
///
/// Represents a 2x2 block on the diagonal of the D factor in P^T A P = L D L^T.
/// The block is symmetric, so only three independent values are stored (the
/// off-diagonal element `b` appears in both positions `D[i, i+1]` and `D[i+1, i]`).
///
/// # References
///
/// - Bunch & Kaufman (1977): 2x2 symmetric pivot block structure in
///   "Some Stable Methods for Calculating Inertia and Solving Symmetric
///   Linear Systems", Math. Comp.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Block2x2 {
    /// Index of the first (lower-indexed) column of the 2x2 block.
    pub first_col: usize,
    /// Top-left element: `D[i, i]`.
    pub a: f64,
    /// Off-diagonal element: `D[i, i+1] = D[i+1, i]`.
    pub b: f64,
    /// Bottom-right element: `D[i+1, i+1]`.
    pub c: f64,
}

impl Block2x2 {
    /// Determinant of the 2x2 block: `ac - b²`.
    pub fn determinant(&self) -> f64 {
        self.a * self.c - self.b * self.b
    }

    /// Trace of the 2x2 block: `a + c`.
    pub fn trace(&self) -> f64 {
        self.a + self.c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- T007: PivotType tests ----

    #[test]
    fn pivot_type_one_by_one() {
        let p = PivotType::OneByOne;
        assert_eq!(p, PivotType::OneByOne);
    }

    #[test]
    fn pivot_type_two_by_two_with_partner() {
        let p = PivotType::TwoByTwo { partner: 5 };
        assert_eq!(p, PivotType::TwoByTwo { partner: 5 });
        assert_ne!(p, PivotType::TwoByTwo { partner: 3 });
    }

    #[test]
    fn pivot_type_delayed() {
        let p = PivotType::Delayed;
        assert_eq!(p, PivotType::Delayed);
        assert_ne!(p, PivotType::OneByOne);
    }

    #[test]
    fn pivot_type_copy_clone() {
        let p = PivotType::TwoByTwo { partner: 2 };
        let p2 = p; // Copy
        #[allow(clippy::clone_on_copy)]
        let p3 = p.clone(); // Clone — intentionally tests Clone impl
        assert_eq!(p, p2);
        assert_eq!(p, p3);
    }

    #[test]
    fn pivot_type_debug_formatting() {
        let one = PivotType::OneByOne;
        let two = PivotType::TwoByTwo { partner: 7 };
        let delayed = PivotType::Delayed;
        // Just verify Debug doesn't panic and produces non-empty output
        assert!(!format!("{:?}", one).is_empty());
        assert!(format!("{:?}", two).contains("7"));
        assert!(!format!("{:?}", delayed).is_empty());
    }

    #[test]
    fn pivot_type_eq_is_structural() {
        // Verify Eq: all variants compare correctly
        assert_eq!(PivotType::OneByOne, PivotType::OneByOne);
        assert_eq!(PivotType::Delayed, PivotType::Delayed);
        assert_ne!(PivotType::OneByOne, PivotType::Delayed);
        assert_ne!(
            PivotType::TwoByTwo { partner: 1 },
            PivotType::TwoByTwo { partner: 2 }
        );
    }

    // ---- T008: Block2x2 tests ----

    #[test]
    fn block2x2_field_storage() {
        let block = Block2x2 {
            first_col: 3,
            a: 2.0,
            b: 0.5,
            c: -3.0,
        };
        assert_eq!(block.first_col, 3);
        assert_eq!(block.a, 2.0);
        assert_eq!(block.b, 0.5);
        assert_eq!(block.c, -3.0);
    }

    #[test]
    fn block2x2_determinant() {
        // [[2, 0.5], [0.5, -3]] → det = 2*(-3) - 0.5^2 = -6.25
        let block = Block2x2 {
            first_col: 0,
            a: 2.0,
            b: 0.5,
            c: -3.0,
        };
        assert_eq!(block.determinant(), -6.25);
    }

    #[test]
    fn block2x2_determinant_identity() {
        // [[1, 0], [0, 1]] → det = 1
        let block = Block2x2 {
            first_col: 0,
            a: 1.0,
            b: 0.0,
            c: 1.0,
        };
        assert_eq!(block.determinant(), 1.0);
    }

    #[test]
    fn block2x2_trace() {
        let block = Block2x2 {
            first_col: 0,
            a: 2.0,
            b: 0.5,
            c: -3.0,
        };
        assert_eq!(block.trace(), -1.0);
    }

    #[test]
    fn block2x2_trace_positive() {
        let block = Block2x2 {
            first_col: 0,
            a: 5.0,
            b: 1.0,
            c: 3.0,
        };
        assert_eq!(block.trace(), 8.0);
    }

    #[test]
    fn block2x2_copy_clone() {
        let block = Block2x2 {
            first_col: 0,
            a: 1.0,
            b: 2.0,
            c: 3.0,
        };
        let b2 = block; // Copy
        #[allow(clippy::clone_on_copy)]
        let b3 = block.clone(); // Clone — intentionally tests Clone impl
        assert_eq!(block, b2);
        assert_eq!(block, b3);
    }

    #[test]
    fn block2x2_partial_eq() {
        let b1 = Block2x2 {
            first_col: 0,
            a: 1.0,
            b: 2.0,
            c: 3.0,
        };
        let b2 = Block2x2 {
            first_col: 0,
            a: 1.0,
            b: 2.0,
            c: 3.0,
        };
        let b3 = Block2x2 {
            first_col: 1,
            a: 1.0,
            b: 2.0,
            c: 3.0,
        };
        assert_eq!(b1, b2);
        assert_ne!(b1, b3);
    }
}
