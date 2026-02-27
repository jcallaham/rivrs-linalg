//! Profiling integration tests for MC64 matching and scaling.
//!
//! - `profile_mc64_ci_subset`: quick CI sanity check with ProfileSession
//! - `profile_mc64_suitesparse_full`: full collection with wall-clock timing,
//!   throughput, RSS tracking, and red flag detection

use std::time::Instant;

use rivrs_sparse::benchmarking::read_current_rss_kb;
use rivrs_sparse::io::registry;
use rivrs_sparse::ordering::{Mc64Job, mc64_matching};
use rivrs_sparse::profiling::ProfileSession;
use rivrs_sparse::testing::{TestCaseFilter, load_test_cases};

#[test]
fn profile_mc64_ci_subset() {
    let cases =
        load_test_cases(&TestCaseFilter::ci_subset()).expect("failed to load CI-subset matrices");

    let suitesparse: Vec<_> = cases
        .iter()
        .filter(|c| c.properties.source == "suitesparse")
        .collect();

    assert!(!suitesparse.is_empty(), "expected SuiteSparse CI matrices");

    let session = ProfileSession::new();

    for case in &suitesparse {
        let _guard = session.enter_section(&case.name);
        mc64_matching(&case.matrix, Mc64Job::MaximumProduct)
            .unwrap_or_else(|e| panic!("MC64 failed for '{}': {}", case.name, e));
    }

    let finished = session.finish();
    let report = finished.summary_report();

    // Verify profiling produced meaningful output
    assert!(
        report.contains("Profile Summary"),
        "report should have header"
    );
    for case in &suitesparse {
        assert!(
            report.contains(&case.name),
            "report should contain matrix '{}'",
            case.name
        );
    }

    eprintln!("\n{}", report);
}

/// Maximum nnz for MC64 profiling (same cap as mc64_matching tests).
const MAX_NNZ_FOR_MC64: usize = 10_000_000;

/// Minimum number of SuiteSparse matrices to consider the full collection present.
const MIN_FULL_COLLECTION_SIZE: usize = 20;

struct MatrixTiming {
    name: String,
    n: usize,
    nnz: usize,
    time_ms: f64,
    throughput_nnz_per_ms: f64,
    rss_kb: Option<u64>,
    rss_delta_kb: Option<i64>,
}

#[test]
#[ignore = "requires full SuiteSparse collection"]
fn profile_mc64_suitesparse_full() {
    let all_meta = registry::load_registry().expect("failed to load metadata.json");
    let mut entries: Vec<_> = all_meta
        .into_iter()
        .filter(|m| m.source == "suitesparse")
        .collect();

    // Sort by nnz ascending for progressive profiling
    entries.sort_by_key(|m| m.nnz);

    eprintln!(
        "Profiling MC64 on {} SuiteSparse entries (nnz cap: {})",
        entries.len(),
        MAX_NNZ_FOR_MC64
    );

    let session = ProfileSession::new();
    let mut timings: Vec<MatrixTiming> = Vec::new();
    let mut skipped_large = 0usize;
    let mut skipped_parse = 0usize;
    let mut loaded = 0usize;

    for meta in &entries {
        if meta.nnz > MAX_NNZ_FOR_MC64 {
            skipped_large += 1;
            continue;
        }

        let test_matrix = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(tm)) => tm,
            Ok(None) => continue,
            Err(_) => {
                skipped_parse += 1;
                continue;
            }
        };
        loaded += 1;

        let rss_before = read_current_rss_kb();

        let start = Instant::now();
        {
            let _guard = session.enter_section(&meta.name);
            mc64_matching(&test_matrix.matrix, Mc64Job::MaximumProduct)
                .unwrap_or_else(|e| panic!("MC64 failed for '{}': {}", meta.name, e));
        }
        let elapsed = start.elapsed();
        let time_ms = elapsed.as_secs_f64() * 1000.0;

        let rss_after = read_current_rss_kb();
        let rss_delta_kb = match (rss_before, rss_after) {
            (Some(before), Some(after)) => Some(after as i64 - before as i64),
            _ => None,
        };

        let throughput = if time_ms > 0.0 {
            meta.nnz as f64 / time_ms
        } else {
            f64::INFINITY
        };

        timings.push(MatrixTiming {
            name: meta.name.clone(),
            n: meta.size,
            nnz: meta.nnz,
            time_ms,
            throughput_nnz_per_ms: throughput,
            rss_kb: rss_after,
            rss_delta_kb,
        });

        // Drop matrix data before loading next
        drop(test_matrix);
    }

    if loaded < MIN_FULL_COLLECTION_SIZE {
        eprintln!(
            "Only {} SuiteSparse matrices loaded (need >= {}). \
             Extract the full collection to test-data/suitesparse/.",
            loaded, MIN_FULL_COLLECTION_SIZE,
        );
        return;
    }

    // Print formatted table
    eprintln!();
    eprintln!(
        "  {:<30} {:>8} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "Matrix", "n", "nnz", "Time(ms)", "nnz/ms", "RSS(KB)", "Delta(KB)"
    );
    eprintln!("  {:-<98}", "");

    for t in &timings {
        let rss = t
            .rss_kb
            .map(|v| format!("{}", v))
            .unwrap_or_else(|| "N/A".to_string());
        let delta = t
            .rss_delta_kb
            .map(|v| format!("{:+}", v))
            .unwrap_or_else(|| "N/A".to_string());
        eprintln!(
            "  {:<30} {:>8} {:>10} {:>10.2} {:>10.0} {:>10} {:>10}",
            t.name, t.n, t.nnz, t.time_ms, t.throughput_nnz_per_ms, rss, delta
        );
    }

    // Summary statistics
    let total_time_ms: f64 = timings.iter().map(|t| t.time_ms).sum();
    let total_nnz: usize = timings.iter().map(|t| t.nnz).sum();
    let overall_throughput = if total_time_ms > 0.0 {
        total_nnz as f64 / total_time_ms
    } else {
        0.0
    };

    let slowest = timings
        .iter()
        .max_by(|a, b| a.time_ms.total_cmp(&b.time_ms));
    let lowest_throughput = timings
        .iter()
        .min_by(|a, b| a.throughput_nnz_per_ms.total_cmp(&b.throughput_nnz_per_ms));
    let max_rss_delta = timings
        .iter()
        .filter_map(|t| t.rss_delta_kb)
        .max_by_key(|d| d.unsigned_abs());

    eprintln!();
    eprintln!("  Summary ({} matrices profiled):", timings.len());
    eprintln!("    Total time:        {:.1} ms", total_time_ms);
    eprintln!("    Total nnz:         {}", total_nnz);
    eprintln!("    Overall throughput: {:.0} nnz/ms", overall_throughput);

    if let Some(s) = slowest {
        eprintln!(
            "    Slowest:           {} ({:.2} ms, {} nnz)",
            s.name, s.time_ms, s.nnz
        );
    }
    if let Some(lt) = lowest_throughput {
        eprintln!(
            "    Lowest throughput: {} ({:.0} nnz/ms, {} nnz)",
            lt.name, lt.throughput_nnz_per_ms, lt.nnz
        );
    }
    if let Some(delta) = max_rss_delta {
        eprintln!("    Max RSS delta:     {:+} KB", delta);
    }

    if skipped_large > 0 {
        eprintln!(
            "    Skipped (nnz > {}): {}",
            MAX_NNZ_FOR_MC64, skipped_large
        );
    }
    if skipped_parse > 0 {
        eprintln!("    Skipped (parse error): {}", skipped_parse);
    }

    // Red flag detection
    let mut red_flags = Vec::new();

    for t in &timings {
        if t.throughput_nnz_per_ms < 100.0 {
            red_flags.push(format!(
                "LOW THROUGHPUT: {} — {:.0} nnz/ms (threshold: 100)",
                t.name, t.throughput_nnz_per_ms
            ));
        }
        if let Some(delta) = t.rss_delta_kb {
            if delta > 100 * 1024 {
                red_flags.push(format!(
                    "HIGH RSS DELTA: {} — {:+} KB ({:.1} MB)",
                    t.name,
                    delta,
                    delta as f64 / 1024.0
                ));
            }
        }
    }

    if !red_flags.is_empty() {
        eprintln!();
        eprintln!("  RED FLAGS ({}):", red_flags.len());
        for flag in &red_flags {
            eprintln!("    {}", flag);
        }
    }

    // Finish profiling session and print summary
    let finished = session.finish();
    eprintln!("\n{}", finished.summary_report());
}
