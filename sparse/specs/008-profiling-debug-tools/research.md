# Research: Profiling and Debug Tools

**Feature**: 008-profiling-debug-tools
**Date**: 2026-02-08

## R1: Chrome Trace Event Format for Profiling Export

**Decision**: Use Chrome Trace Event Format with Complete Events ("X" type).

**Rationale**: Chrome Trace Event format is widely supported (Chrome `chrome://tracing`, Perfetto, Trace Compass) and requires no external dependencies — it's just JSON. Complete Events (type "X") are simpler to emit than Begin/End pairs and represent the full duration in a single record.

**Key format details**:
- JSON object with `traceEvents` array
- Required fields per event: `name` (string), `ph` ("X"), `ts` (microseconds since trace start), `dur` (microseconds), `pid` (process ID), `tid` (thread ID)
- Optional: `cat` (category string for filtering), `args` (key-value metadata)
- Hierarchy is implicit: viewer infers nesting from overlapping time ranges on the same tid
- Timestamps in microseconds (float allowed for sub-µs)

**Alternatives considered**:
- **`tracing` crate integration**: Powerful but heavyweight dependency; adds async runtime concerns and ecosystem coupling not needed for this use case
- **Custom binary format**: Better performance but requires custom viewer tooling
- **Flamegraph SVG**: Good for static analysis but loses temporal information

## R2: Thread-Safe Profiling Architecture

**Decision**: Thread-local storage with synchronous merge at session end.

**Rationale**: Thread-local storage (`thread_local!` + `RefCell`) provides zero-contention recording — critical for sub-microsecond overhead. The solver's profiling sessions are well-bounded (one matrix factorization = one session), making synchronous merge at session end simple and sufficient. No background thread or channel infrastructure needed.

**Architecture**:
- Each thread records `ProfileEvent` structs into a thread-local `Vec`
- Events use `std::time::Instant` for monotonic, cross-thread-comparable timestamps
- At session end, the session owner collects all thread-local buffers and merges by timestamp
- Hierarchy reconstructed from overlapping time ranges (same logic as Chrome Trace viewer)
- Aggregation (total time, call count, min/max/mean) computed during merge

**Alternatives considered**:
- **Shared `Mutex<Vec>`**: High contention at sub-µs recording frequency; serializes all threads
- **Lock-free channels (`crossbeam`)**: Adds external dependency; overkill when sessions are bounded and merge is synchronous
- **Single-threaded only**: Rejected per clarification — must be thread-safe from day one

**Zero-cost pattern**: Use `#[cfg(feature = "test-util")]` to gate the profiling module entirely. When disabled, profiling structs and methods don't exist. Instrumentation call sites use a macro that expands to nothing when the feature is off.

**Timing**: `std::time::Instant` is monotonic, `Send + Sync`, and comparable across threads. Use `checked_duration_since()` for robustness against rare platform anomalies.

## R3: RSS Memory Reading — Existing Infrastructure

**Decision**: Extend `benchmarking/rss.rs` with current-RSS reading; import from `crate::benchmarking::rss` in the memory tracker. No code move needed.

**Rationale**: The existing `read_peak_rss_kb()` in `src/benchmarking/rss.rs` already implements the Linux `/proc/self/status` reading pattern with proper `#[cfg(target_os = "linux")]` gating. The memory tracker needs both peak RSS (`VmHWM`) and current RSS (`VmRSS`). Since both `benchmarking` and `profiling` are behind the same `test-util` feature flag, the memory tracker can import directly from `crate::benchmarking::rss` without factoring into a separate shared module. A shared utility can be extracted later if the dependency direction becomes awkward.

**Existing implementation** (`benchmarking/rss.rs`):
- `pub fn read_peak_rss_kb() -> Option<u64>` — reads `VmHWM` from `/proc/self/status`
- Linux-only with `#[cfg(target_os = "linux")]`; returns `None` on other platforms
- No external dependencies (pure `std::fs`)
- Currently used in `benches/solver_benchmarks.rs` for before/after snapshots

**Extension needed**:
- Add `read_current_rss_kb() -> Option<u64>` reading `VmRSS` line in `benchmarking/rss.rs`
- Memory tracker imports from `crate::benchmarking::rss` directly

**Alternatives considered**:
- **Global allocator wrapper** (`#[global_allocator]`): Tracks exact allocations but invasive, affects all code, complex to attribute to labeled sections
- **`jemalloc` statistics**: Requires jemalloc dependency; not portable
- **`libc::getrusage`**: Less granular than `/proc/self/status`; doesn't provide current RSS on Linux

## R4: Text-Based Sparsity Pattern Visualization

**Decision**: Implement from scratch using character grid with density-based Unicode block characters and configurable downsampling.

**Rationale**: No existing Rust crate provides a lightweight text-based spy plot. The implementation is ~150-200 LOC and avoids adding dependencies. Unicode block characters (`█▓▒░`) provide density-level visualization in downsampled cells, are widely supported in modern terminals, and degrade gracefully to ASCII if needed.

**Approach**:
- **Small matrices** (n ≤ display width): 1:1 mapping, `#` for non-zero, `.` for zero
- **Large matrices** (n > display width): Bin rows and columns into display cells. Count non-zeros per bin. Map density to Unicode blocks:
  - 0%: `.` (space/dot)
  - 1-25%: `░` (light)
  - 26-50%: `▒` (medium)
  - 51-75%: `▓` (dark)
  - 76-100%: `█` (full)
- **Header**: Matrix name, dimensions (n×m), nnz count, density percentage, downsampling ratio if applicable
- **Configurable**: Max display width (default 80), max display height (default 40), ASCII-only mode fallback

**Alternatives considered**:
- **Braille characters** (U+2800-28FF): Higher resolution (2×4 pixels per char) but less readable for density visualization
- **External crate** (`plotters`, `textplots`): Adds dependency for a simple feature; plotters is image-focused, not terminal-focused
- **SVG/image output**: Deferred per spec assumption — text output is sufficient for this phase
