//! Compare dense APTP kernel accuracy across different inner block sizes.
//!
//! Creates a dense symmetric indefinite matrix and factors it with different
//! inner_block_size values. Compares reconstruction error to verify that
//! blocking does not introduce accuracy regressions.
//!
//! Usage:
//!   cargo run --release --example dense_aptp_blocking [matrix_size]
//!
//! Default matrix size is 256.

use faer::Mat;
use rivrs_sparse::aptp::factor::{AptpOptions, aptp_factor};
use rivrs_sparse::aptp::pivot::PivotType;

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(256);

    let block_sizes = [16, 32, 64, 128, 256, 100000];

    eprintln!("=== Dense APTP blocking comparison, n={} ===\n", n);
    eprintln!(
        "{:<12} {:>6}/{:<6} {:>12} {:>8} {:>6}",
        "config", "ne", "n", "recon_err", "delays", "2x2"
    );
    eprintln!("{}", "-".repeat(60));

    // Deterministic symmetric indefinite matrix with small diagonal (forces 2x2 pivots)
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    let a = Mat::from_fn(n, n, |i, j| {
        let (lo, hi) = if i >= j { (j, i) } else { (i, j) };
        let mut h = DefaultHasher::new();
        (lo, hi, n).hash(&mut h);
        let v = (h.finish() as f64 / u64::MAX as f64) * 2.0 - 1.0;
        if i == j { v * 0.01 } else { v }
    });

    for &ib in &block_sizes {
        if ib > n * 2 && ib != 100000 {
            continue;
        }
        let options = AptpOptions {
            inner_block_size: ib.min(100000),
            outer_block_size: ib.clamp(256, 100000),
            ..AptpOptions::default()
        };

        let label = if ib >= 100000 {
            "unblocked".to_string()
        } else {
            format!("ib={}", ib)
        };

        match aptp_factor(a.as_ref(), &options) {
            Ok(result) => {
                let l = &result.l;
                let ne = l.ncols();
                if ne < 2 {
                    eprintln!("{:<12} {:>6}/{:<6} — too few eliminated", label, ne, n);
                    continue;
                }
                let d = &result.d;
                let perm = &result.perm;

                // Build dense L and D
                let mut l_dense = Mat::zeros(n, ne);
                for j in 0..ne {
                    l_dense[(j, j)] = 1.0;
                    for i in (j + 1)..n {
                        l_dense[(i, j)] = l[(i, j)];
                    }
                }

                let mut d_dense = Mat::zeros(ne, ne);
                let mut col = 0;
                while col < ne {
                    match d.get_pivot_type(col) {
                        PivotType::OneByOne => {
                            d_dense[(col, col)] = d.get_1x1(col);
                            col += 1;
                        }
                        PivotType::TwoByTwo { partner } if partner > col && col + 1 < ne => {
                            let blk = d.get_2x2(col);
                            d_dense[(col, col)] = blk.a;
                            d_dense[(col + 1, col)] = blk.b;
                            d_dense[(col, col + 1)] = blk.b;
                            d_dense[(col + 1, col + 1)] = blk.c;
                            col += 2;
                        }
                        _ => {
                            col += 1;
                        }
                    }
                }

                // Compute P^T A P vs L D L^T
                let (fwd, _) = perm.as_ref().arrays();
                let mut pap = Mat::zeros(n, n);
                for i in 0..n {
                    for j in 0..n {
                        pap[(i, j)] = a[(fwd[i], fwd[j])];
                    }
                }

                let ld = &l_dense * &d_dense;
                let ldlt = &ld * l_dense.transpose();

                let mut max_err = 0.0f64;
                let mut norm_a = 0.0f64;
                for i in 0..ne {
                    for j in 0..ne {
                        max_err = max_err.max((pap[(i, j)] - ldlt[(i, j)]).abs());
                        norm_a = norm_a.max(pap[(i, j)].abs());
                    }
                }
                let rel_err = if norm_a > 0.0 {
                    max_err / norm_a
                } else {
                    max_err
                };

                eprintln!(
                    "{:<12} {:>6}/{:<6} {:>12.2e} {:>8} {:>6}",
                    label, ne, n, rel_err, result.stats.num_delayed, result.stats.num_2x2
                );
            }
            Err(e) => {
                eprintln!("{:<12} FAILED: {}", label, e);
            }
        }
    }
}
