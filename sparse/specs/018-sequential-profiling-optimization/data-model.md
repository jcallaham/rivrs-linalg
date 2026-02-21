# Data Model: Sequential Profiling & Optimization (Phase 8.1g)

**Date**: 2026-02-20

## Entities

### Extended Entity: PerSupernodeStats

Existing per-supernode statistics collected during multifrontal factorization. Phase 8.1g extends this with optional timing fields.

**Existing fields** (unchanged):
- `snode_id: usize` — Supernode index in postorder
- `front_size: usize` — Total frontal matrix dimension (m)
- `num_fully_summed: usize` — Fully-summed rows (k)
- `num_eliminated: usize` — Columns successfully factored (ne)
- `num_delayed: usize` — Delayed columns (k - ne)
- `num_1x1: usize` — 1x1 pivots accepted
- `num_2x2: usize` — 2x2 pivot pairs accepted
- `max_l_entry: f64` — Stability metric (max |L| entry)

**New fields** (conditional on `diagnostic` feature):
- `assembly_time: Duration` — Wall-clock time for scatter + extend-add
- `kernel_time: Duration` — Wall-clock time for dense APTP kernel
- `extraction_time: Duration` — Wall-clock time for front factor extraction

**Invariants**:
- `num_eliminated + num_delayed == num_fully_summed`
- `num_1x1 + 2 * num_2x2 == num_eliminated`
- When timing is collected: `assembly_time + kernel_time + extraction_time ≈ total supernode time` (modulo overhead)

---

### Extended Entity: FactorizationStats

Existing aggregate statistics. Phase 8.1g extends with total timing breakdown.

**Existing fields** (unchanged):
- `total_1x1_pivots: usize`
- `total_2x2_pivots: usize`
- `total_delayed: usize`
- `zero_pivots: usize`
- `max_front_size: usize`

**New fields** (conditional on `diagnostic` feature):
- `total_assembly_time: Duration` — Sum of assembly times across all supernodes
- `total_kernel_time: Duration` — Sum of kernel times across all supernodes
- `total_extraction_time: Duration` — Sum of extraction times across all supernodes

**Invariants**:
- `total_1x1_pivots == sum(per_snode.num_1x1)`
- `total_assembly_time == sum(per_snode.assembly_time)`

---

### New Entity: PerformanceBaseline

A structured record capturing per-matrix performance data for one complete solve pipeline.

**Fields**:
- `matrix_name: String` — SuiteSparse matrix identifier
- `matrix_dim: usize` — Matrix dimension (n)
- `matrix_nnz: usize` — Number of nonzeros
- `ordering_time: Duration` — Time for ordering computation (METIS/MC64+METIS)
- `symbolic_time: Duration` — Time for symbolic analysis
- `factor_time: Duration` — Time for numeric factorization
- `solve_time: Duration` — Time for triangular solve
- `total_time: Duration` — End-to-end time (ordering + symbolic + factor + solve)
- `peak_rss_kb: Option<u64>` — Peak RSS during factorization (platform-dependent)
- `backward_error: f64` — Solution backward error
- `num_supernodes: usize` — Number of supernodes in assembly tree
- `max_front_size: usize` — Largest frontal matrix dimension
- `factorization_stats: FactorizationStats` — Aggregate pivot/delay statistics
- `per_supernode_stats: Vec<PerSupernodeStats>` — Per-front breakdown (with timing if diagnostic)

**Invariants**:
- `total_time >= ordering_time + symbolic_time + factor_time + solve_time` (overhead from setup/teardown)
- `backward_error > 0.0` (non-negative)
- `per_supernode_stats.len() == num_supernodes`

---

### New Entity: BaselineSuite

Collection of baselines for the full SuiteSparse test suite, with metadata.

**Fields**:
- `timestamp: String` — ISO-8601 collection timestamp
- `platform: String` — Machine/OS description
- `solver_version: String` — Git commit hash or tag
- `ordering_strategy: String` — Ordering used (e.g., "MatchOrderMetis")
- `baselines: Vec<PerformanceBaseline>` — Per-matrix baselines

**Invariants**:
- `baselines.len() == 65` (full SuiteSparse suite) or subset for CI

---

### New Entity: WorkloadProfile

Analysis of workload distribution for a single matrix, derived from per-supernode timing data.

**Fields**:
- `matrix_name: String`
- `total_factor_time: Duration`
- `top_10pct_time_fraction: f64` — Fraction of total time in top 10% of fronts (by size)
- `top_1_front_time_fraction: f64` — Fraction of total time in single largest front
- `front_size_histogram: Vec<(usize, usize)>` — (front_size_bucket, count) pairs
- `time_by_front_size: Vec<(usize, Duration)>` — (front_size_bucket, total_time) pairs
- `parallelism_recommendation: ParallelismClass` — Tree-level / Intra-node / Mixed

**ParallelismClass variants**:
- `TreeLevel` — Many small independent fronts; tree-level parallelism most beneficial
- `IntraNode` — Few large fronts dominate; intra-node BLAS-3 parallelism most beneficial
- `Mixed` — Significant time in both small and large fronts; both strategies needed

---

## Relationships

```
BaselineSuite
  └── 1:N ── PerformanceBaseline (one per matrix)
                ├── 1:1 ── FactorizationStats
                ├── 1:N ── PerSupernodeStats (one per supernode)
                └── derives ── WorkloadProfile

WorkloadProfile
  └── derived from PerformanceBaseline.per_supernode_stats timing data
```
