# Research: Sequential Profiling & Optimization (Phase 8.1g)

**Date**: 2026-02-20

## R1: Profiling Integration Strategy

**Decision**: Use existing `ProfileSession`/`SectionGuard` infrastructure with `diagnostic` feature gating. Instrument at per-supernode granularity within `AptpNumeric::factor()`.

**Rationale**: The profiling module (`src/profiling/`) is already feature-complete — thread-safe hierarchical timing, Chrome Trace export, RAII guards. No new infrastructure needed; just need to call the existing API from the right places.

**Integration Points** (in `AptpNumeric::factor()` postorder loop, lines 380-488 of numeric.rs):
- `"delayed_collection"` — lines 383-389 (gather delayed columns from children)
- `"assembly"` — lines 412-435 (scatter original entries + extend-add child contributions)
- `"dense_kernel"` — line 438 (call `aptp_factor_in_place`)
- `"extraction"` — lines 468-470 (call `extract_front_factors`)

**Alternatives considered**:
- External profiler only (cargo-flamegraph, perf): Rejected — doesn't provide per-supernode attribution or structured output
- Custom timing framework: Rejected — duplicates existing ProfileSession capabilities
- Always-on timing (no feature gate): Rejected — adds overhead to production builds

---

## R2: Allocation Hotspot Analysis

**Decision**: Target the top 3 allocation hotspots by frequency, all in `factor_inner` (factor.rs):

1. **Panel row permutation (line 1664)** — ~15,000+ allocations per large front. A `Vec<f64>` of size `block_size` is allocated for *every panel row* in *every block iteration*. Fix: hoist a single reusable buffer outside the loop.

2. **Row permutation temp buffer (line 1692)** — one `vec![0.0; block_size]` per block. Fix: hoist outside the block loop.

3. **Column order slice copy (line 1704)** — one `to_vec()` per block. Fix: use a pre-allocated buffer or stack array.

**Secondary targets** (lower frequency but larger size):
- `BlockBackup::create` (line 1146) — one `Mat::zeros(rows, block_cols)` per block. Could pre-allocate and reuse.
- `compute_ld` workspace (line 1440) — 0-2 per block. Could pre-allocate.

**Not targeted** (per-supernode, variable size, hard to pool):
- Frontal matrix `Mat::zeros(m, m)` (numeric.rs line 414) — each supernode has different `m`
- `frontal_rows` vector (numeric.rs line 396) — different pattern per supernode

**Rationale**: Line 1664 accounts for ~95% of allocations per large front. Fixing just this one hotspot eliminates ~15,000 allocations per `factor_inner` call on a 1000×1000 front.

**Alternatives considered**:
- Arena allocator for all frontal data: Rejected for Phase 8.1g — too invasive, deferred to Phase 9.1
- SmallVec for block_size buffers: Rejected — a simple hoisted `Vec` is cleaner and sufficient
- Global workspace pool across supernodes: Rejected — variable supernode sizes make pooling complex

---

## R3: Baseline Collection Strategy

**Decision**: Create a new `examples/baseline_collection.rs` tool that extends the `solve_timing.rs` pattern with profiling, per-supernode stats, and structured JSON export.

**Rationale**: The existing infrastructure provides all building blocks:
- `BenchmarkResult` / `BenchmarkSuiteResult` for phase-level timing (in `src/benchmarking/`)
- `PerSupernodeStats` already collected during factorization (accessible via `SparseLDLT::per_supernode_stats()`)
- `MemoryTracker` for RSS snapshots
- `ProfileSession` for Chrome Trace export
- `Baseline` struct with save/load and `detect_regressions()`

Phase 8.1g adds a new `PerformanceBaseline` struct that bundles all per-matrix data (timing per phase, per-supernode stats, peak RSS) and exports to JSON for Phase 8.2 comparison.

**Storage**: `target/benchmarks/baselines/` (already used by Criterion baseline framework)

**Alternatives considered**:
- Extend Criterion benchmarks directly: Rejected — Criterion is designed for statistical micro-benchmarks, not systematic single-run collection with profiling
- External profiling tools only: Rejected — doesn't produce structured data for programmatic analysis

---

## R4: Workload Distribution Analysis Approach

**Decision**: Compute workload distribution from per-supernode timing data collected during baseline runs. Classify matrices into three categories based on front size distribution.

**Matrix classes** (to inform Phase 8.2 parallelism):
1. **Small/sparse** — many small fronts, tree parallelism matters most
2. **Medium** — mixed front sizes, both parallelism types may help
3. **Large/dense fronts** — few large fronts dominate, intra-node BLAS-3 parallelism matters most

**Key metric**: fraction of total factorization time in top 10% of fronts (by size). If >80%, intra-node parallelism is critical. If <30%, tree-level parallelism dominates.

**Existing infrastructure**: `examples/front_sizes.rs` already does symbolic-only front size analysis. `examples/supernode_stats.rs` already outputs per-supernode pivot details. Phase 8.1g combines these with actual timing data.

**Alternatives considered**:
- Symbolic-only analysis (front sizes without timing): Rejected — front size is a proxy for time but doesn't account for varying pivot complexity and delay rates
- Per-matrix individual reports: Rejected for primary output — too granular. Summary report with matrix classes is more actionable for Phase 8.2.

---

## R5: Feature Flag for Instrumentation

**Decision**: Use existing `diagnostic` feature flag.

**Rationale**: Already in Cargo.toml, semantically correct for profiling/diagnostic instrumentation. Avoids adding overhead to `test-util` builds (which enable random generators for testing). The `diagnostic` flag is currently unused in the codebase, making it a clean fit.

**Alternatives considered**:
- `test-util`: Rejected — would add profiling overhead to all test builds
- New `profiling` feature: Rejected — unnecessary proliferation of feature flags

---

## R6: PerSupernodeStats Timing Extension

**Decision**: Add three `Duration`-based timing fields to `PerSupernodeStats`, conditionally compiled behind `diagnostic` feature.

**Fields to add**:
- `assembly_time`: Time for scatter + extend-add
- `kernel_time`: Time for dense APTP kernel
- `extraction_time`: Time for front factor extraction

**Rationale**: These map directly to the three main phases within each supernode's factorization. Using `Duration` (not `u64` nanoseconds) is idiomatic Rust and avoids unit confusion.

**Conditional compilation**: Fields present only when `diagnostic` is enabled. Default `PerSupernodeStats` (production) remains unchanged. This avoids bloating the per-supernode storage for non-profiling builds.

**Alternatives considered**:
- HashMap<String, Duration>: Rejected — runtime overhead, not statically typed
- Separate TimingStats struct: Rejected — unnecessary indirection when three fixed fields suffice
- Always-present timing fields: Rejected — adds memory overhead in production builds where timing is not collected
