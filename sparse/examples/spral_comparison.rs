//! Compare our match_order_metis against SPRAL's on SuiteSparse matrices.
//!
//! Compares:
//! - Scaling vectors (should match closely since MC64 is deterministic)
//! - Fill quality via AptpSymbolic (predicted nnz within reasonable bounds)
//!
//! # Prerequisites
//!
//! Build the SPRAL Fortran driver (see tools/spral_match_order.f90 for instructions):
//!
//! ```sh
//! mkdir -p /tmp/spral_mo
//! gfortran -O2 -c -o /tmp/spral_mo/matrix_util.o /opt/references/spral/src/matrix_util.f90 -J /tmp/spral_mo
//! gfortran -O2 -c -o /tmp/spral_mo/scaling.o /opt/references/spral/src/scaling.f90 -J /tmp/spral_mo -I /tmp/spral_mo
//! gfortran -O2 -c -o /tmp/spral_mo/metis5_wrapper.o -DSPRAL_HAVE_METIS_H=0 /opt/references/spral/src/metis5_wrapper.F90 -J /tmp/spral_mo -I /tmp/spral_mo
//! gfortran -O2 -c -o /tmp/spral_mo/match_order.o /opt/references/spral/src/match_order.f90 -J /tmp/spral_mo -I /tmp/spral_mo
//! METIS_LIB=$(find target -name "libmetis.a" | head -1)
//! gfortran -O2 -I /tmp/spral_mo -o /tmp/spral_match_order \
//!   /tmp/spral_mo/matrix_util.o /tmp/spral_mo/scaling.o /tmp/spral_mo/metis5_wrapper.o /tmp/spral_mo/match_order.o \
//!   tools/spral_match_order.f90 $METIS_LIB -lm
//! ```
//!
//! # Usage
//!
//! ```sh
//! # Compare on CI-subset matrices
//! cargo run --release --example spral_comparison
//!
//! # Compare on specific matrix files
//! cargo run --release --example spral_comparison -- test-data/suitesparse/GHS_indef/sparsine.mtx
//!
//! # Compare on all matrices in a directory
//! cargo run --release --example spral_comparison -- test-data/suitesparse-ci/**/*.mtx
//! ```

use std::io::Write;
use std::path::Path;
use std::process::Command;

use faer::sparse::SparseColMat;
use faer::sparse::linalg::cholesky::SymmetricOrdering;

use rivrs_sparse::aptp::{AptpSymbolic, match_order_metis};
use rivrs_sparse::io::mtx::load_mtx;
use rivrs_sparse::testing::{TestCaseFilter, load_test_cases};

const SPRAL_MATCH_ORDER_BIN: &str = "/tmp/spral_match_order";

/// Export a full symmetric CSC matrix in the format expected by the SPRAL
/// match_order_metis Fortran driver (1-indexed, full both triangles).
fn export_full_csc_for_spral(matrix: &SparseColMat<usize, f64>) -> String {
    let n = matrix.nrows();
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();
    let nnz = col_ptrs[n];

    let mut out = String::new();
    out.push_str(&format!("{} {}\n", n, nnz));

    // Column pointers (1-indexed)
    for &p in &col_ptrs[..=n] {
        out.push_str(&format!("{}\n", p + 1));
    }

    // Row indices and values (1-indexed rows)
    for k in 0..nnz {
        out.push_str(&format!("{} {:.17e}\n", row_indices[k] + 1, values[k]));
    }

    out
}

/// Parse SPRAL match_order_metis output.
/// Returns (scaling, ordering) where ordering[i] = position of variable i (0-indexed).
fn parse_spral_match_order_output(stdout: &str) -> Option<(Vec<f64>, Vec<usize>)> {
    let lines = stdout.lines();
    let mut scaling = Vec::new();
    let mut ordering = Vec::new();
    let mut n = 0usize;
    let mut section = "";

    for line in lines {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("N ") {
            n = rest.trim().parse().ok()?;
        } else if let Some(rest) = line.strip_prefix("FLAG ") {
            let flag: i32 = rest.trim().parse().ok()?;
            if flag < 0 {
                return None; // SPRAL error
            }
        } else if line == "SCALING" {
            section = "scaling";
        } else if line == "ORDERING" {
            section = "ordering";
        } else if section == "scaling" && scaling.len() < n {
            scaling.push(line.parse::<f64>().ok()?);
        } else if section == "ordering" && ordering.len() < n {
            // SPRAL ordering is 1-indexed
            let pos: usize = line.trim().parse().ok()?;
            ordering.push(pos - 1); // Convert to 0-indexed
        }
    }

    if scaling.len() == n && ordering.len() == n {
        Some((scaling, ordering))
    } else {
        None
    }
}

/// Run SPRAL's match_order_metis on a matrix and return (scaling, ordering).
fn run_spral_match_order(matrix: &SparseColMat<usize, f64>) -> Option<(Vec<f64>, Vec<usize>)> {
    let input = export_full_csc_for_spral(matrix);

    let mut child = Command::new(SPRAL_MATCH_ORDER_BIN)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .ok()?;

    child.stdin.as_mut()?.write_all(input.as_bytes()).ok()?;

    let output = child.wait_with_output().ok()?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        eprintln!("SPRAL match_order_metis failed: {}", stderr);
        return None;
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_spral_match_order_output(&stdout)
}

/// Compare a single matrix against SPRAL.
fn compare_matrix(name: &str, matrix: &SparseColMat<usize, f64>) -> bool {
    let n = matrix.nrows();

    // Run our pipeline
    let our_result = match match_order_metis(matrix) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("  {:<40} OUR FAILED: {}", name, e);
            return false;
        }
    };

    // Run SPRAL pipeline
    let (spral_scaling, spral_ordering) = match run_spral_match_order(matrix) {
        Some(s) => s,
        None => {
            eprintln!("  {:<40} SPRAL FAILED", name);
            return false;
        }
    };

    // Compare scaling vectors
    if our_result.scaling.len() != spral_scaling.len() {
        eprintln!(
            "  {:<40} SCALING LENGTH MISMATCH: {} vs {}",
            name,
            our_result.scaling.len(),
            spral_scaling.len()
        );
        return false;
    }

    let mut max_diff = 0.0_f64;
    let mut sum_sq_diff = 0.0_f64;
    for (ours, &spral_s) in our_result.scaling.iter().zip(&spral_scaling) {
        let diff = (ours - spral_s).abs();
        let rel_diff = if spral_s.abs() > 1e-300 {
            diff / spral_s.abs()
        } else {
            diff
        };
        max_diff = max_diff.max(rel_diff);
        sum_sq_diff += rel_diff * rel_diff;
    }
    let rms_diff = (sum_sq_diff / n as f64).sqrt();

    // Compare fill quality via AptpSymbolic — our ordering
    let our_nnz = match AptpSymbolic::analyze(
        matrix.symbolic(),
        SymmetricOrdering::Custom(our_result.ordering.as_ref()),
    ) {
        Ok(s) => s.predicted_nnz(),
        Err(e) => {
            eprintln!("  {:<40} OUR SYMBOLIC FAILED: {}", name, e);
            return false;
        }
    };

    // SPRAL ordering: order(i) = position of variable i = inv[old_idx]
    let spral_inv = spral_ordering;
    let mut spral_fwd = vec![0usize; n];
    for (old_idx, &new_pos) in spral_inv.iter().enumerate() {
        if new_pos >= n {
            eprintln!(
                "  {:<40} SPRAL ORDERING OUT OF RANGE: order[{}] = {}",
                name, old_idx, new_pos
            );
            return false;
        }
        spral_fwd[new_pos] = old_idx;
    }

    let spral_perm = faer::perm::Perm::new_checked(
        spral_fwd.into_boxed_slice(),
        spral_inv.into_boxed_slice(),
        n,
    );
    let spral_nnz = match AptpSymbolic::analyze(
        matrix.symbolic(),
        SymmetricOrdering::Custom(spral_perm.as_ref()),
    ) {
        Ok(s) => s.predicted_nnz(),
        Err(e) => {
            eprintln!("  {:<40} SPRAL SYMBOLIC FAILED: {}", name, e);
            return false;
        }
    };

    let fill_ratio = our_nnz as f64 / spral_nnz.max(1) as f64;

    let scaling_ok = max_diff < 1e-10;
    let fill_ok = fill_ratio < 2.0 && fill_ratio > 0.5;
    let status = if scaling_ok && fill_ok {
        "OK"
    } else {
        "MISMATCH"
    };

    eprintln!(
        "  {:<40} {:>8} {:>12.2e} {:>12.2e} {:>12} {:>12} {:>10.3}  {}",
        name, n, max_diff, rms_diff, our_nnz, spral_nnz, fill_ratio, status,
    );

    if !scaling_ok {
        eprintln!(
            "    WARNING: scaling max relative diff = {:.2e} (threshold 1e-10)",
            max_diff
        );
    }
    if !fill_ok {
        eprintln!(
            "    WARNING: fill ratio {:.3} outside [0.5, 2.0]",
            fill_ratio
        );
    }

    scaling_ok && fill_ok
}

fn main() {
    if !Path::new(SPRAL_MATCH_ORDER_BIN).exists() {
        eprintln!(
            "Error: {} not found.\nSee tools/spral_match_order.f90 for build instructions.",
            SPRAL_MATCH_ORDER_BIN
        );
        std::process::exit(1);
    }

    let args: Vec<String> = std::env::args().skip(1).collect();

    eprintln!(
        "\n  {:<40} {:>8} {:>12} {:>12} {:>12} {:>12} {:>10}",
        "Matrix", "n", "max_s_diff", "rms_s_diff", "our_nnzL", "spral_nnzL", "fill_ratio"
    );
    eprintln!("  {:-<112}", "");

    let mut pass = 0;
    let mut fail = 0;

    if args.is_empty() {
        // Default: run on CI-subset matrices
        let cases = load_test_cases(&TestCaseFilter::ci_subset())
            .expect("failed to load CI-subset matrices");

        let suitesparse: Vec<_> = cases
            .iter()
            .filter(|c| c.properties.source == "suitesparse")
            .collect();

        for case in &suitesparse {
            if compare_matrix(&case.name, &case.matrix) {
                pass += 1;
            } else {
                fail += 1;
            }
        }
    } else {
        // Run on specified matrix files
        for path_str in &args {
            let path = Path::new(path_str);
            if !path.exists() {
                eprintln!("  {:<40} FILE NOT FOUND", path_str);
                fail += 1;
                continue;
            }

            let matrix = match load_mtx(path) {
                Ok(m) => m,
                Err(e) => {
                    eprintln!("  {:<40} LOAD FAILED: {}", path_str, e);
                    fail += 1;
                    continue;
                }
            };

            // Use filename as label
            let name = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or(path_str);

            // Try to include parent directory for context
            let label = if let Some(parent) = path.parent().and_then(|p| p.file_name()) {
                format!("{}/{}", parent.to_str().unwrap_or(""), name)
            } else {
                name.to_string()
            };

            if compare_matrix(&label, &matrix) {
                pass += 1;
            } else {
                fail += 1;
            }
        }
    }

    eprintln!("\n  {}/{} passed", pass, pass + fail);

    if fail > 0 {
        std::process::exit(1);
    }
}
