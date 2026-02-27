//! Example: Solving a discrete-time Sylvester equation AXB + X = C.
//!
//! Demonstrates basic usage of the `solve_discrete` function.

use faer::prelude::*;
use rivrs_control::sylvester::solve_discrete;

fn main() {
    // Define the matrices for AXB + X = C
    let a = mat![[0.5, 0.1, 0.0], [0.0, 0.8, 0.2], [0.0, 0.0, 0.3f64]];
    let b = mat![[0.6, 0.1], [0.0, 0.9f64]];
    let c = mat![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0f64]];

    println!("Solving discrete-time Sylvester equation: AXB + X = C");
    println!("A = {:?}", &a);
    println!("B = {:?}", &b);
    println!("C = {:?}", &c);

    match solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()) {
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
                "Residual norm ||AXB + X - C||: {:.2e}",
                result.residual_norm
            );
            println!("Near-singular: {}", result.near_singular);
        }
        Err(e) => {
            eprintln!("Error: {}", e);
        }
    }
}
