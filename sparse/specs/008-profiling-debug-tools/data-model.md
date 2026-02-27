# Data Model: Profiling and Debug Tools

**Feature**: 008-profiling-debug-tools
**Date**: 2026-02-08

## Entities

### ProfileSession

A bounded profiling run collecting hierarchical timing data, potentially from multiple threads.

| Field | Type | Description |
|-------|------|-------------|
| start_time | Instant | Monotonic start of the session |
| thread_events | Map<ThreadId, Vec\<ProfileEvent\>> | Raw events per thread (pre-merge) |
| merged_sections | Vec\<ProfileSection\> | Aggregated hierarchy (post-merge) |

**Lifecycle**: Created → Recording (events appended per-thread) → Finished (merge + aggregate) → Exported/Reported

**Identity**: One session per profiling scope (e.g., one factorization run). Not persisted.

### ProfileEvent

A raw timing record from a single thread. Intermediate representation — aggregated into ProfileSection during merge.

| Field | Type | Description |
|-------|------|-------------|
| name | &str / String | Section label |
| start | Instant | When section entered |
| duration | Duration | Wall-clock time spent |
| thread_id | ThreadId | Originating thread |
| depth | usize | Nesting depth (0 = top level) |

**Note**: Events are append-only during recording. Consumed during merge.

### ProfileSection

An aggregated view of one or more invocations of a named section at a specific nesting position.

| Field | Type | Description |
|-------|------|-------------|
| name | String | Section label |
| total_duration | Duration | Sum of all invocations |
| call_count | u64 | Number of times entered |
| min_duration | Duration | Shortest invocation |
| max_duration | Duration | Longest invocation |
| mean_duration | Duration | Average (total / count) |
| parent_pct | f64 | Percentage of parent section's total time |
| children | Vec\<ProfileSection\> | Nested subsections (recursive) |

**Aggregation rule**: Sections with the same name at the same depth under the same parent are aggregated. Cross-thread sections with the same label are aggregated together in the merged view.

### MemorySnapshot

A point-in-time RSS measurement annotated with a label.

| Field | Type | Description |
|-------|------|-------------|
| label | String | Developer-provided section name |
| timestamp | Instant | When snapshot was taken |
| rss_kb | Option\<u64\> | Current RSS in KB (None if unavailable) |
| peak_rss_kb | Option\<u64\> | Peak RSS (VmHWM) in KB (None if unavailable) |

**Platform behavior**: On Linux, both fields populated from `/proc/self/status`. On other platforms, both are `None`.

### MemoryReport

Aggregated view of a sequence of memory snapshots.

| Field | Type | Description |
|-------|------|-------------|
| snapshots | Vec\<MemorySnapshot\> | Ordered sequence of snapshots |
| peak_rss_kb | Option\<u64\> | Maximum peak_rss_kb across all snapshots |
| deltas | Vec\<MemoryDelta\> | Computed differences between consecutive snapshots |

### MemoryDelta

Computed difference between two consecutive snapshots.

| Field | Type | Description |
|-------|------|-------------|
| from_label | String | Label of the earlier snapshot |
| to_label | String | Label of the later snapshot |
| rss_delta_kb | Option\<i64\> | Change in current RSS (can be negative) |

### SparsityDisplay

A rendered text representation of a matrix's non-zero structure.

| Field | Type | Description |
|-------|------|-------------|
| matrix_rows | usize | Original matrix row count |
| matrix_cols | usize | Original matrix column count |
| matrix_nnz | usize | Original non-zero count |
| display_rows | usize | Rendered grid rows |
| display_cols | usize | Rendered grid columns |
| cells | Vec\<Vec\<f64\>\> | Density per cell (0.0 to 1.0) |

**Rendering**: Density mapped to characters — `.` (0%), `░` (1-25%), `▒` (26-50%), `▓` (51-75%), `█` (76-100%). ASCII fallback: `.` (0%), `-` (1-33%), `+` (34-66%), `#` (67-100%).

### EliminationTreeStats

Summary statistics for an elimination tree.

| Field | Type | Description |
|-------|------|-------------|
| num_nodes | usize | Total nodes |
| num_leaves | usize | Leaf count (nodes with no children) |
| depth | usize | Maximum root-to-leaf path length |
| branching_min | usize | Minimum children per non-leaf |
| branching_max | usize | Maximum children per non-leaf |
| branching_mean | f64 | Average children per non-leaf |
| subtree_sizes | Vec\<usize\> | Size of subtree rooted at each node |

## Relationships

```
ProfileSession 1──* ProfileEvent     (raw, pre-merge)
ProfileSession 1──* ProfileSection   (aggregated, post-merge)
ProfileSection 1──* ProfileSection   (parent-children hierarchy)

MemoryReport   1──* MemorySnapshot   (ordered sequence)
MemoryReport   1──* MemoryDelta      (computed from consecutive snapshots)
```

## State Transitions

### ProfileSession Lifecycle

```
Created ──[enter_section()]──► Recording
Recording ──[enter/exit_section()]──► Recording  (repeated)
Recording ──[finish()]──► Finished
Finished ──[summary_report()]──► Finished  (read-only)
Finished ──[export_chrome_trace()]──► Finished  (read-only)
```

### MemoryReport Lifecycle

```
Created ──[snapshot()]──► Tracking
Tracking ──[snapshot()]──► Tracking  (repeated)
Tracking ──[report()]──► Complete
Complete ──[display()]──► Complete  (read-only)
```
