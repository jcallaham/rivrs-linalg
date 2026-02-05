//! Example: Demonstrating error handling for singular Sylvester equations.
//!
//! When A and -B share eigenvalues, the Sylvester equation AX + XB = C
//! is singular and cannot be solved. This example demonstrates how the
//! solver detects and reports such cases.

use faer::prelude::*;
use rivrs_control::error::SylvesterError;
use rivrs_control::sylvester::solve_continuous;

fn main() {
    println!("=== Singular Case Detection ===\n");

    // Case 1: Common eigenvalues (exactly singular)
    // A has eigenvalue 2, B has eigenvalue -2
    // λ(A) + λ(B) = 2 + (-2) = 0 => singular
    println!("Case 1: A=[2], B=[-2] (common eigenvalues)");
    let a = mat![[2.0f64]];
    let b = mat![[-2.0f64]];
    let c = mat![[1.0f64]];

    match solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()) {
        Ok(_) => println!("  Unexpected success"),
        Err(SylvesterError::CommonEigenvalues {
            separation,
            threshold,
        }) => {
            println!("  Detected: common eigenvalues");
            println!(
                "  Separation: {:.2e} (threshold: {:.2e})",
                separation, threshold
            );
        }
        Err(e) => println!("  Other error: {}", e),
    }

    // Case 2: Dimension mismatch
    println!("\nCase 2: Dimension mismatch (A is 2x2, B is 3x3, C is 2x2)");
    let a = mat![[1.0, 0.0], [0.0, 2.0f64]];
    let b = mat![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0f64]];
    let c = mat![[1.0, 0.0], [0.0, 1.0f64]];

    match solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()) {
        Ok(_) => println!("  Unexpected success"),
        Err(e) => println!("  Error: {}", e),
    }

    // Case 3: NaN input
    println!("\nCase 3: NaN in input matrix");
    let a = mat![[f64::NAN, 0.0], [0.0, 1.0]];
    let b = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let c = mat![[1.0, 0.0], [0.0, 1.0f64]];

    match solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()) {
        Ok(_) => println!("  Unexpected success"),
        Err(e) => println!("  Error: {}", e),
    }

    // Case 4: Well-conditioned problem (for comparison)
    println!("\nCase 4: Well-conditioned problem (A=[1], B=[5])");
    let a = mat![[1.0f64]];
    let b = mat![[5.0f64]];
    let c = mat![[12.0f64]];

    match solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()) {
        Ok(result) => {
            let x = &result.solution * (1.0 / result.scale);
            println!("  Solution: x = {:.6}", x[(0, 0)]);
            println!("  Residual: {:.2e}", result.residual_norm);
            println!("  Near-singular: {}", result.near_singular);
        }
        Err(e) => println!("  Unexpected error: {}", e),
    }
}
