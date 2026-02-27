# Quickstart: Profiling and Debug Tools

**Feature**: 008-profiling-debug-tools

## Prerequisites

- Rust 1.87+ (edition 2024)
- `rivrs-sparse` with `test-util` feature enabled

All profiling and debug tools are gated behind the `test-util` feature flag:

```toml
[dependencies]
rivrs-sparse = { path = "../sparse", features = ["test-util"] }
```

In tests and benchmarks, the self-referencing dev-dependency already enables this.

## 1. Profile Solver Components

Wrap code sections with the profiler to get a hierarchical timing breakdown:

```rust
use rivrs_sparse::profiling::ProfileSession;

let mut session = ProfileSession::new();

{
    let _guard = session.enter_section("analyze");
    // ... symbolic analysis work ...
    {
        let _inner = session.enter_section("build_etree");
        // ... elimination tree construction ...
    }
}

{
    let _guard = session.enter_section("factor");
    for node in nodes {
        let _kernel = session.enter_section("dense_kernel");
        // ... factorize node ...
    }
}

let finished = session.finish();

// Text summary to terminal
println!("{}", finished.summary_report());

// Export for Chrome Trace viewer (chrome://tracing or Perfetto)
std::fs::write("trace.json", finished.export_chrome_trace()).unwrap();
```

### Summary Report Output

```
Profile Summary (total: 45.2 ms)
──────────────────────────────────────────────────────────────────
  Section           Total      Calls   Mean       Min        Max      Parent%
  analyze           12.3 ms    1       12.3 ms    12.3 ms    12.3 ms  27.2%
    build_etree      8.1 ms    1        8.1 ms     8.1 ms     8.1 ms  65.9%
  factor            32.9 ms    1       32.9 ms    32.9 ms    32.9 ms  72.8%
    dense_kernel    28.4 ms    150      0.19 ms    0.12 ms    0.45 ms  86.3%
```

### View in Chrome Trace Viewer

1. Open `chrome://tracing` in Chrome (or use Perfetto UI)
2. Click "Load" and select the `trace.json` file
3. Navigate the hierarchical timeline view

## 2. Track Memory Usage

Record RSS snapshots at key points to understand memory allocation patterns:

```rust
use rivrs_sparse::profiling::memory::MemoryTracker;

let mut tracker = MemoryTracker::new();
tracker.snapshot("start");

// ... load matrix ...
tracker.snapshot("matrix_loaded");

// ... factorize ...
tracker.snapshot("factorized");

// ... solve ...
tracker.snapshot("solved");

let report = tracker.report();
println!("{}", report.display_report());
```

### Memory Report Output

```
Memory Report
────────────────────────────────────────────────────────
  Label                 RSS (KB)   Peak (KB)   Delta (KB)
  start                 12,340     12,340      —
  matrix_loaded         14,560     14,560      +2,220
  factorized            28,900     28,900      +14,340
  solved                29,100     29,100      +200
────────────────────────────────────────────────────────
  Peak RSS: 29,100 KB (28.4 MB)
```

**Note**: Memory tracking uses `/proc/self/status` and is only available on Linux. On other platforms, RSS values will be reported as "unavailable".

## 3. Visualize Sparsity Patterns

Inspect matrix structure in the terminal:

```rust
use rivrs_sparse::debug::SparsityDisplay;

let matrix = /* load a SparseColMat<usize, f64> */;

// Default: 80-column Unicode display
println!("{}", SparsityDisplay::from_sparse(&matrix));

// Customized: narrower, ASCII-only
let display = SparsityDisplay::from_sparse(&matrix)
    .with_max_width(40)
    .with_ascii_only(true);
println!("{}", display.render());
```

## 4. Inspect Elimination Trees

Visualize tree structure and compute statistics:

```rust
use rivrs_sparse::debug::ETreeDisplay;

let parent_array: &[usize] = /* from symbolic analysis */;
let etree = ETreeDisplay::from_parent_array(parent_array);

// For small trees: text tree diagram
println!("{}", etree.render_tree());

// For any size: summary statistics
println!("{}", etree.render_stats());

// Programmatic access
let stats = etree.stats();
println!("Tree depth: {}, Leaves: {}", stats.depth, stats.num_leaves);
```

## Feature Gating

All profiling and debug tools are excluded from production builds. When `test-util` is not enabled:
- Profiling types exist as zero-sized types
- All methods compile to no-ops
- No runtime overhead
- No additional dependencies pulled in
