# Data Model: Direct GEMM into Contribution Buffer

**Feature**: 024-direct-gemm-contrib
**Date**: 2026-02-24 (revised)

## Modified Entities

### FactorizationWorkspace (modified)

Existing pre-allocated workspace for multifrontal factorization, extended with a contribution buffer.

| Field | Type | Description |
|-------|------|-------------|
| `frontal_data` | Dense matrix (m×m) | Existing. Reusable frontal matrix, max_front × max_front |
| `frontal_row_indices` | Integer array | Existing. Global row indices for current frontal matrix |
| `delayed_cols_buf` | Integer array | Existing. Collects delayed columns from children |
| `global_to_local` | Integer array | Existing. Global-to-local index mapping, length n |
| **`contrib_buffer`** | **Dense matrix** | **NEW. Pre-allocated contribution buffer, max_front × max_front. Receives the deferred contribution GEMM output. Reused across supernodes via swap protocol.** |

**Invariants**:
- `contrib_buffer` is sized ≥ max contribution dimension for any supernode
- `contrib_buffer` may be temporarily empty (moved into ContributionBlock) and refilled via swap
- When empty (moved out), the workspace cannot receive a deferred GEMM until a buffer is swapped back in
- The deferred GEMM writes the NFS×NFS Schur complement directly into this buffer

### ContributionBlock (unchanged interface, different data source)

The Schur complement from a factored supernode. Interface unchanged; the NFS×NFS data is now written directly by the deferred GEMM rather than copied from the frontal workspace.

| Field | Type | Description |
|-------|------|-------------|
| `data` | Dense matrix (size×size) | Schur complement data (lower triangle). NFS×NFS portion from deferred GEMM; delayed portion copied from workspace. |
| `row_indices` | Integer array | Global permuted column indices for rows/columns |
| `num_delayed` | Integer | Number of delayed columns at start of block |

**Layout within `data`**:
```
[0..num_delayed, 0..num_delayed]     → delayed × delayed interactions (from workspace)
[num_delayed..size, 0..num_delayed]  → NFS × delayed cross-terms (from workspace)
[num_delayed..size, num_delayed..size] → NFS × NFS Schur complement (from deferred GEMM)
```

### FrontalMatrix (unchanged)

Temporary borrowed view into the workspace. No changes.

### AptpKernelWorkspace (unchanged)

Dense APTP kernel workspace buffers (ld_buf, copy_buf). No changes to buffer structure, though their usage changes: the ld_buf and copy_buf may be reused by the deferred contribution GEMM.

## Per-Block Trailing Update Decomposition

The key structural change: the trailing submatrix `A[ts..m, ts..m]` is decomposed relative to the fixed `num_fully_summed` boundary (`p`):

```
                    ts (block_end)          p (num_fully_summed)           m
                        │                         │                        │
ts (block_end)     ─────┼─────────────────────────┼────────────────────────┤
                        │                         │                        │
                        │    Region 1: FS × FS    │   Region 2: FS × NFS  │
                        │    (lower-triangular)    │   (rectangular)        │
                        │    UPDATED PER-BLOCK     │   UPDATED PER-BLOCK    │
                        │                         │                        │
p                  ─────┼─────────────────────────┼────────────────────────┤
                        │                         │                        │
                        │   Region 2': NFS × FS   │   Region 3: NFS × NFS │
                        │   (rectangular)          │   (lower-triangular)   │
                        │   UPDATED PER-BLOCK      │   ** DEFERRED **       │
                        │                         │                        │
m                  ─────┼─────────────────────────┼────────────────────────┘
```

- **Region 1** (FS×FS): Updated per-block. Read by subsequent factor + TRSM.
- **Region 2/2'** (cross-terms): Updated per-block. NFS rows of FS columns read by subsequent TRSM.
- **Region 3** (NFS×NFS): **Deferred** to single post-loop GEMM → writes to `contrib_buffer`.

## State Transitions

### Contribution Buffer Lifecycle (Sequential Path)

```
State 1: WORKSPACE_HOLDS_BUFFER
  workspace.contrib_buffer = allocated (max_front × max_front)

  → Blocking loop runs (update_trailing skips NFS×NFS)
  → Deferred GEMM writes NFS×NFS Schur complement into contrib_buffer
  → Small delayed portion (if any) copied from workspace into contrib_buffer
  → contrib_buffer MOVED into ContributionBlock

State 2: BUFFER_IN_CONTRIBUTION
  workspace.contrib_buffer = empty (moved out)
  contributions[s] = Some(ContributionBlock { data: <the buffer>, ... })

  → Parent's extend_add() reads from ContributionBlock
  → After extend_add, ContributionBlock consumed
  → Its data (Mat<f64>) RETURNED to workspace.contrib_buffer

State 3: WORKSPACE_HOLDS_BUFFER (recycled)
  workspace.contrib_buffer = recycled buffer (may be oversized — no shrinking)

  → If recycled buffer too small for next contribution: reallocate (rare)
  → Otherwise: reuse as-is
  → Back to State 1 for next supernode
```

### Contribution Buffer Lifecycle (Parallel Path)

```
State 1: THREAD_LOCAL_WORKSPACE_HOLDS_BUFFER
  thread_workspace.contrib_buffer = allocated

  → factor_single_supernode: blocking loop + deferred GEMM + extraction
  → contrib_buffer MOVED into ContributionBlock result

State 2: BUFFER_IN_CROSS_THREAD_CONTRIBUTION
  ContributionBlock moved to parent's thread via Vec<Option<ContributionBlock>>
  Thread's workspace.contrib_buffer = empty

  → Thread picks up next supernode from wave
  → If workspace.contrib_buffer is empty: allocate a new one

State 3: PARENT_THREAD_CONSUMES
  Parent's extend_add reads ContributionBlock
  Returned buffer swapped into parent's workspace.contrib_buffer
```

## Relationships

```
AptpSymbolic ──predicts──> max_front_size ──sizes──> FactorizationWorkspace
                                                        │
                                                        ├── frontal_data (existing)
                                                        └── contrib_buffer (NEW)
                                                              ▲
                                                              │ WRITE (deferred GEMM)
                                                              │
                                              factor_inner ───┘
                                              (blocking loop skips NFS×NFS,
                                               single GEMM after loop writes
                                               NFS×NFS to contrib_buffer)
                                                              │
                                                              │ MOVE (swap into ContributionBlock)
                                                              ▼
                                                        ContributionBlock
                                                              │
                                                              │ READ (extend_add)
                                                              ▼
                                                        FrontalMatrix (parent)
                                                              │
                                                              │ RECYCLE (return buffer)
                                                              ▼
                                                        FactorizationWorkspace.contrib_buffer
```
