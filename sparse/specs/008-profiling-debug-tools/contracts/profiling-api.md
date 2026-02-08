# API Contract: Profiling Module

**Module path**: `rivrs_sparse::profiling`
**Feature gate**: `test-util`

## Core Types

### ProfileSession

```
ProfileSession
  ::new() -> ProfileSession
  .enter_section(name: &str) -> SectionGuard
  .finish() -> FinishedSession
```

- `new()`: Creates a session and records the start timestamp.
- `enter_section(name)`: Begins timing a named section. Returns a guard that stops timing on drop. Supports nesting — calling `enter_section` while another guard is alive creates a child section.
- `finish()`: Consumes the session, merges thread-local events, computes aggregations. Returns a read-only `FinishedSession`.

**Thread safety**: `ProfileSession` is `Send + Sync`. Each thread records events into thread-local storage. `finish()` collects all thread-local buffers.

### SectionGuard

```
SectionGuard (RAII)
  drop(): Records the section's duration and depth.
```

- Drop-based timing ensures sections are closed even on early return or panic.
- Must not be sent across threads (not `Send`).

### FinishedSession

```
FinishedSession
  .summary_report() -> String
  .export_chrome_trace() -> String
  .sections() -> &[ProfileSection]
  .total_duration() -> Duration
```

- `summary_report()`: Hierarchical text table showing name, total time, call count, min/max/mean, parent %.
- `export_chrome_trace()`: JSON string in Chrome Trace Event format (Complete Events, type "X").
- `sections()`: Programmatic access to aggregated section tree.
- `total_duration()`: Wall-clock time from session start to finish.

## Zero-Cost No-Op Pattern

When `test-util` feature is disabled:

```
ProfileSession  → zero-sized type
SectionGuard    → zero-sized type
enter_section() → returns ZST guard, no timing
finish()        → returns FinishedSession with empty data
summary_report() → returns empty string
export_chrome_trace() → returns empty string
```

All methods exist with the same signatures but compile to no-ops. Instrumented code compiles without `cfg` guards at call sites.

## Usage Pattern

```
let mut session = ProfileSession::new();

{
    let _outer = session.enter_section("factorize");
    {
        let _inner = session.enter_section("dense_kernel");
        // ... work ...
    } // _inner dropped, inner section timed
} // _outer dropped, outer section timed

let finished = session.finish();
println!("{}", finished.summary_report());
std::fs::write("trace.json", finished.export_chrome_trace());
```

## Chrome Trace Export Format

```json
{
  "traceEvents": [
    {
      "name": "factorize",
      "cat": "profiling",
      "ph": "X",
      "ts": 0.0,
      "dur": 1500.3,
      "pid": 1,
      "tid": 1,
      "args": { "call_count": 1 }
    },
    {
      "name": "dense_kernel",
      "cat": "profiling",
      "ph": "X",
      "ts": 10.2,
      "dur": 1200.1,
      "pid": 1,
      "tid": 1,
      "args": { "call_count": 50 }
    }
  ]
}
```

- Timestamps in microseconds relative to session start.
- Individual invocations preserved (not aggregated) in trace export.
- Aggregated summary available via `summary_report()`.
