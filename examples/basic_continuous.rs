//! Example: Solving a continuous-time Sylvester equation AX + XB = C.
//!
//! Demonstrates basic usage of the `solve_continuous` function.

use rivrs_linalg::sylvester::solve_continuous;
use faer::prelude::*;

fn main() {
    // Define the matrices for AX + XB = C
    let a = mat![[1.0, 2.0, 0.0], [0.0, 3.0, 1.0], [0.0, 0.0, 5.0f64]];
    let b = mat![[2.0, 1.0], [0.0, 4.0f64]];
    let c = mat![[10.0, 20.0], [30.0, 40.0], [50.0, 60.0f64]];

    println!("Solving continuous-time Sylvester equation: AX + XB = C");
    println!("A = {:?}", &a);
    println!("B = {:?}", &b);
    println!("C = {:?}", &c);

    match solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()) {
        Ok(result) => {
            let x = &result.solution * (1.0 / result.scale);
            println!("\nSolution X:");
            for i in 0..x.nrows() {
                for j in 0..x.ncols() {
                    print!("  {:.6}", x[(i, j)]);
                }
                println!();
            }
            println!("\nScale factor: {}", result.scale);
            println!(
                "Residual norm ||AX + XB - C||: {:.2e}",
                result.residual_norm
            );
            println!("Near-singular: {}", result.near_singular);
        }
        Err(e) => {
            eprintln!("Error: {}", e);
        }
    }
}
