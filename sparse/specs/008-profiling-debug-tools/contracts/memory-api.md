# API Contract: Memory Tracking Module

**Module path**: `rivrs_sparse::profiling::memory`
**Feature gate**: `test-util`

## Core Types

### MemoryTracker

```
MemoryTracker
  ::new() -> MemoryTracker
  .snapshot(label: &str)
  .report() -> MemoryReport
```

- `new()`: Creates a tracker, optionally taking an initial "baseline" snapshot.
- `snapshot(label)`: Records the current process RSS and peak RSS with the given label. On non-Linux platforms, records `None` for both values.
- `report()`: Computes deltas between consecutive snapshots and returns a `MemoryReport`.

### MemoryReport

```
MemoryReport
  .peak_rss_kb() -> Option<u64>
  .snapshots() -> &[MemorySnapshot]
  .deltas() -> &[MemoryDelta]
  .display_report() -> String
```

- `peak_rss_kb()`: Maximum peak RSS observed across all snapshots.
- `snapshots()`: Ordered list of snapshots.
- `deltas()`: Computed differences between consecutive snapshots.
- `display_report()`: Formatted text table suitable for terminal output.

## RSS Reading Functions (Shared Utility)

```
read_current_rss_kb() -> Option<u64>
read_peak_rss_kb() -> Option<u64>
```

- `read_current_rss_kb()`: Reads `VmRSS` from `/proc/self/status` (Linux only).
- `read_peak_rss_kb()`: Reads `VmHWM` from `/proc/self/status` (Linux only). Factored out from existing `benchmarking::rss` module.
- Both return `None` on non-Linux platforms.

## Display Report Format

```
Memory Report
────────────────────────────────────────────────────────
  Label                 RSS (KB)   Peak (KB)   Delta (KB)
  baseline              12,340     12,340      —
  after_analyze         14,560     14,560      +2,220
  after_factor          28,900     28,900      +14,340
  after_solve           29,100     29,100      +200
────────────────────────────────────────────────────────
  Peak RSS: 29,100 KB (28.4 MB)
```

## Usage Pattern

```
let mut tracker = MemoryTracker::new();
tracker.snapshot("baseline");

// ... symbolic analysis ...
tracker.snapshot("after_analyze");

// ... numeric factorization ...
tracker.snapshot("after_factor");

// ... solve ...
tracker.snapshot("after_solve");

let report = tracker.report();
println!("{}", report.display_report());
```
