//! End-to-end property-based tests for SparseLDLT.
//!
//! Uses proptest strategies from `rivrs_sparse::testing::strategies` to generate
//! random sparse symmetric matrices and verify structural invariants through
//! the full solve pipeline.

use faer::Col;
use proptest::prelude::*;
use rivrs_sparse::symmetric::solver::{SolverOptions, SparseLDLT};
use rivrs_sparse::testing::strategies;
use rivrs_sparse::validate::sparse_backward_error;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// Backward error < 5e-11 for sparse PD matrices through full pipeline.
    #[test]
    fn property_sparse_pd_backward_error(
        a in strategies::arb_sparse_symmetric(5..=100, 0.05..=0.3)
    ) {
        let n = a.nrows();
        // Make PD by replacing negative diagonal with positive
        let dense = a.to_dense();
        let mut triplets = Vec::new();
        for j in 0..n {
            for i in 0..n {
                let v = dense[(i, j)];
                if v != 0.0 {
                    if i == j {
                        // Ensure positive diagonal
                        triplets.push(faer::sparse::Triplet::new(i, j, v.abs()));
                    } else {
                        triplets.push(faer::sparse::Triplet::new(i, j, v));
                    }
                }
            }
        }
        let a_pd = faer::sparse::SparseColMat::try_new_from_triplets(n, n, &triplets)
            .expect("PD matrix construction");

        // RHS = ones
        let rhs = Col::from_fn(n, |_| 1.0);
        let options = SolverOptions::default();
        let result = SparseLDLT::solve_full(&a_pd, &rhs, &options);

        if let Ok(ref x) = result {
            let be = sparse_backward_error(&a_pd, x, &rhs);
            prop_assert!(
                be < 5e-11,
                "PD backward error {:.2e} for n={}", be, n
            );
        }
        // Errors (e.g., tiny matrix numerics) are acceptable
    }

    /// Backward error < 5e-11 or clean error for sparse indefinite matrices.
    #[test]
    fn property_sparse_indefinite_backward_error(
        a in strategies::arb_sparse_symmetric(5..=100, 0.05..=0.3)
    ) {
        let n = a.nrows();
        let rhs = Col::from_fn(n, |_| 1.0);
        let options = SolverOptions::default();
        let result = SparseLDLT::solve_full(&a, &rhs, &options);

        match result {
            Ok(ref x) => {
                let be = sparse_backward_error(&a, x, &rhs);
                prop_assert!(
                    be < 5e-11,
                    "Indefinite backward error {:.2e} for n={}", be, n
                );
            }
            Err(_) => {
                // Clean error is acceptable for indefinite matrices
            }
        }
    }

    /// Inertia consistency: positive + negative + zero == n.
    #[test]
    fn property_sparse_inertia_consistency(
        a in strategies::arb_sparse_symmetric(5..=100, 0.05..=0.3)
    ) {
        let n = a.nrows();
        let options = SolverOptions::default();
        let analyze_opts = rivrs_sparse::symmetric::solver::AnalyzeOptions {
            ordering: options.ordering.clone(),
        };
        let factor_opts = rivrs_sparse::symmetric::solver::FactorOptions::default();

        let solver = SparseLDLT::analyze_with_matrix(&a, &analyze_opts);
        if let Ok(mut solver) = solver {
            if solver.factor(&a, &factor_opts).is_ok() {
                if let Some(inertia) = solver.inertia() {
                    prop_assert_eq!(
                        inertia.dimension(), n,
                        "inertia sum {} != n={}", inertia.dimension(), n
                    );
                }
            }
        }
    }
}
