//! Compare our solver's backward error against SPRAL SSIDS with identical inputs.
//!
//! Exports our matrix + ordering + scaling + RHS to SPRAL, runs both solvers,
//! and compares backward error and factorization statistics.
//!
//! # Prerequisites
//!
//! Build SPRAL library and driver:
//! ```sh
//! tools/build_spral.sh
//! METIS_LIB=$(find target -name "libmetis.a" | head -1)
//! gfortran -O2 -I /tmp/spral_ssids -o /tmp/spral_full_solve \
//!   tools/spral_full_solve.f90 \
//!   -Wl,--whole-archive /tmp/spral_ssids/libspral.a -Wl,--no-whole-archive \
//!   $METIS_LIB -lopenblas -lstdc++ -lm
//! ```
//!
//! # Usage
//!
//! ```sh
//! # Compare on specific matrix files
//! cargo run --release --example spral_solve_comparison -- \
//!   test-data/suitesparse/easy-indefinite/sparsine/sparsine.mtx
//!
//! # Compare on multiple matrices
//! cargo run --release --example spral_solve_comparison -- \
//!   test-data/suitesparse/easy-indefinite/copter2/copter2.mtx \
//!   test-data/suitesparse/easy-indefinite/dawson5/dawson5.mtx
//! ```

use std::io::Write;
use std::path::Path;
use std::process::Command;

use faer::Col;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::SparseColMat;

use rivrs_sparse::aptp::{
    AnalyzeOptions, FactorOptions, OrderingStrategy, SparseLDLT, match_order_metis,
};
use rivrs_sparse::io::mtx::load_mtx;

const SPRAL_SOLVE_BIN: &str = "/tmp/spral_full_solve";

/// Extract lower triangle from a full symmetric CSC matrix.
/// SPRAL expects lower-triangle-only CSC input.
fn extract_lower_triangle(matrix: &SparseColMat<usize, f64>) -> (Vec<i32>, Vec<i32>, Vec<f64>) {
    let n = matrix.nrows();
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();

    let mut lower_ptr = vec![0i32; n + 1]; // 1-indexed
    let mut lower_row = Vec::new();
    let mut lower_val = Vec::new();

    for j in 0..n {
        lower_ptr[j] = lower_row.len() as i32 + 1; // 1-indexed
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for idx in start..end {
            let i = row_indices[idx];
            if i >= j {
                // Lower triangle (including diagonal)
                lower_row.push(i as i32 + 1); // 1-indexed
                lower_val.push(values[idx]);
            }
        }
    }
    lower_ptr[n] = lower_row.len() as i32 + 1;

    (lower_ptr, lower_row, lower_val)
}

/// Format input for SPRAL full solve driver.
fn format_spral_input(
    matrix: &SparseColMat<usize, f64>,
    ordering_inv: &[usize], // inv[old_idx] = new_pos (SPRAL convention: order(i) = position)
    scaling: &[f64],
    rhs: &[f64],
) -> String {
    let n = matrix.nrows();
    let (lower_ptr, lower_row, lower_val) = extract_lower_triangle(matrix);
    let nnz = lower_row.len();

    let mut out = String::new();
    out.push_str(&format!("{} {}\n", n, nnz));

    // Column pointers (1-indexed)
    for i in 0..=n {
        out.push_str(&format!("{}\n", lower_ptr[i]));
    }

    // Row indices and values
    for k in 0..nnz {
        out.push_str(&format!("{} {:.17e}\n", lower_row[k], lower_val[k]));
    }

    // Ordering (1-indexed: order(i) = position of variable i)
    out.push_str("ORDERING\n");
    for &pos in ordering_inv {
        out.push_str(&format!("{}\n", pos + 1)); // Convert to 1-indexed
    }

    // Scaling factors
    out.push_str("SCALING\n");
    for &s in scaling {
        out.push_str(&format!("{:.17e}\n", s));
    }

    // RHS
    out.push_str("RHS\n");
    for &b in rhs {
        out.push_str(&format!("{:.17e}\n", b));
    }

    out
}

struct SpralResult {
    backward_error: f64,
    num_delay: i64,
    num_two: i64,
    not_first_pass: i64,
    not_second_pass: i64,
    maxfront: i64,
}

fn parse_spral_output(stdout: &str) -> Option<SpralResult> {
    let mut result = SpralResult {
        backward_error: f64::NAN,
        num_delay: 0,
        num_two: 0,
        not_first_pass: 0,
        not_second_pass: 0,
        maxfront: 0,
    };

    for line in stdout.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("BACKWARD_ERROR ") {
            result.backward_error = rest.trim().parse().ok()?;
        } else if let Some(rest) = line.strip_prefix("num_delay ") {
            result.num_delay = rest.trim().parse().ok()?;
        } else if let Some(rest) = line.strip_prefix("num_two ") {
            result.num_two = rest.trim().parse().ok()?;
        } else if let Some(rest) = line.strip_prefix("not_first_pass ") {
            result.not_first_pass = rest.trim().parse().ok()?;
        } else if let Some(rest) = line.strip_prefix("not_second_pass ") {
            result.not_second_pass = rest.trim().parse().ok()?;
        } else if let Some(rest) = line.strip_prefix("maxfront ") {
            result.maxfront = rest.trim().parse().ok()?;
        }
    }

    Some(result)
}

fn compare_matrix(path: &Path) -> Option<()> {
    let name = if let Some(parent) = path.parent().and_then(|p| p.file_name()) {
        format!(
            "{}/{}",
            parent.to_str().unwrap_or(""),
            path.file_stem()?.to_str()?
        )
    } else {
        path.file_stem()?.to_str()?.to_string()
    };

    // Load matrix
    let matrix = match load_mtx(path) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("  {:<40} LOAD FAILED: {}", name, e);
            return None;
        }
    };
    let n = matrix.nrows();

    // Run our match_order_metis
    let mo_result = match match_order_metis(&matrix) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("  {:<40} match_order_metis FAILED: {}", name, e);
            return None;
        }
    };

    let (_, ordering_inv) = mo_result.ordering.as_ref().arrays();

    // Generate RHS: A * ones(n) so we know exact solution
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();

    let mut rhs = vec![0.0f64; n];
    for j in 0..n {
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for idx in start..end {
            let i = row_indices[idx];
            rhs[i] += values[idx]; // A * ones = sum of each row
        }
    }

    // ---- Run our solver ----
    let analyze_opts = AnalyzeOptions {
        ordering: OrderingStrategy::MatchOrderMetis,
    };
    let factor_opts = FactorOptions::default();

    let our_be = match (|| -> Result<f64, Box<dyn std::error::Error>> {
        let mut solver = SparseLDLT::analyze_with_matrix(&matrix, &analyze_opts)?;
        solver.factor(&matrix, &factor_opts)?;

        let rhs_col = Col::from_fn(n, |i| rhs[i]);
        let scratch = solver.solve_scratch(1);
        let mut mem = MemBuffer::new(scratch);
        let stack = MemStack::new(&mut mem);
        let x_col = solver.solve(&rhs_col, stack)?;

        // Compute backward error
        let mut ax = vec![0.0f64; n];
        for j in 0..n {
            let start = col_ptrs[j];
            let end = col_ptrs[j + 1];
            for idx in start..end {
                let i = row_indices[idx];
                ax[i] += values[idx] * x_col[j];
            }
        }
        let norm_r: f64 = ax
            .iter()
            .zip(&rhs)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        let norm_a: f64 = values.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let norm_x: f64 = (0..n).map(|i| x_col[i].abs()).fold(0.0, f64::max);
        let norm_b: f64 = rhs.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let denom = norm_a * norm_x + norm_b;
        Ok(if denom > 0.0 { norm_r / denom } else { 0.0 })
    })() {
        Ok(be) => be,
        Err(e) => {
            eprintln!("  {:<40} OUR SOLVER FAILED: {}", name, e);
            return None;
        }
    };

    // ---- Run SPRAL solver ----
    let input = format_spral_input(&matrix, ordering_inv, &mo_result.scaling, &rhs);

    let mut child = match Command::new(SPRAL_SOLVE_BIN)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
    {
        Ok(c) => c,
        Err(e) => {
            eprintln!("  {:<40} SPRAL SPAWN FAILED: {}", name, e);
            return None;
        }
    };

    if let Some(stdin) = child.stdin.as_mut() {
        if stdin.write_all(input.as_bytes()).is_err() {
            eprintln!("  {:<40} SPRAL STDIN WRITE FAILED", name);
            return None;
        }
    }

    let output = match child.wait_with_output() {
        Ok(o) => o,
        Err(e) => {
            eprintln!("  {:<40} SPRAL WAIT FAILED: {}", name, e);
            return None;
        }
    };

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        eprintln!(
            "  {:<40} SPRAL FAILED (exit {}): {}",
            name, output.status, stderr
        );
        return None;
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr_str = String::from_utf8_lossy(&output.stderr);
    let spral = match parse_spral_output(&stdout) {
        Some(r) => r,
        None => {
            eprintln!("  {:<40} SPRAL PARSE FAILED", name);
            eprintln!("    stderr: {}", stderr_str);
            return None;
        }
    };

    // ---- Report ----
    let ratio = if spral.backward_error > 0.0 {
        our_be / spral.backward_error
    } else {
        f64::NAN
    };

    eprintln!(
        "  {:<35} {:>7} {:>12.2e} {:>12.2e} {:>8.1} {:>8} {:>8} {:>8} {:>8} {:>8}",
        name,
        n,
        our_be,
        spral.backward_error,
        ratio,
        spral.num_delay,
        spral.num_two,
        spral.not_first_pass,
        spral.not_second_pass,
        spral.maxfront,
    );

    Some(())
}

fn main() {
    if !Path::new(SPRAL_SOLVE_BIN).exists() {
        eprintln!(
            "Error: {} not found.\nRun tools/build_spral.sh then compile tools/spral_full_solve.f90.",
            SPRAL_SOLVE_BIN
        );
        std::process::exit(1);
    }

    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        eprintln!("Usage: spral_solve_comparison <matrix.mtx> [matrix2.mtx] ...");
        std::process::exit(1);
    }

    eprintln!(
        "\n  {:<35} {:>7} {:>12} {:>12} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "Matrix", "n", "our_BE", "spral_BE", "ratio", "delays", "2x2", "!1st", "!2nd", "maxfrt"
    );
    eprintln!("  {:-<131}", "");

    let mut pass = 0;
    let mut fail = 0;

    for path_str in &args {
        let path = Path::new(path_str);
        if !path.exists() {
            eprintln!("  {:<35} FILE NOT FOUND", path_str);
            fail += 1;
            continue;
        }

        if compare_matrix(path).is_some() {
            pass += 1;
        } else {
            fail += 1;
        }
    }

    eprintln!("\n  {}/{} compared successfully", pass, pass + fail);
}
