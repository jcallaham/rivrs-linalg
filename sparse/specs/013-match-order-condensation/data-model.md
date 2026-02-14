# Data Model: Match-Order Condensation Pipeline

**Feature**: 013-match-order-condensation
**Date**: 2026-02-13

## Entities

### MatchOrderResult (public output)

The final result of the combined matching-ordering pipeline.

| Field | Type | Description | Constraints |
|-------|------|-------------|-------------|
| ordering | Perm\<usize\> | Fill-reducing permutation with matched pair adjacency | Valid permutation of {0..n-1}; 2-cycle pairs consecutive |
| scaling | Vec\<f64\> | MC64 symmetric scaling factors (linear domain) | All positive; \|s_i * a_ij * s_j\| <= 1 for matched entries |
| matched | usize | Number of matched entries from MC64 | 0 <= matched <= n; equals n for nonsingular |
| condensed_dim | usize | Dimension of the condensed graph | condensed_dim <= n; strictly < n when 2-cycles exist |
| singletons | usize | Count of singleton nodes in cycle decomposition | singletons + 2 * two_cycles + unmatched = n |
| two_cycles | usize | Count of 2-cycle pairs in cycle decomposition | two_cycles >= 1 implies condensed_dim < n |

**Relationships**: Contains the MC64 scaling/matched data (from Mc64Result) plus the fill-reducing ordering and condensation diagnostics. Consumed by AptpSymbolic::analyze() (ordering) and future Phase 5 numeric factorization (scaling).

**Lifecycle**: Created once per pipeline invocation. Immutable after construction.

### CycleDecomposition (internal intermediate)

Result of splitting an MC64 matching permutation into singletons, 2-cycles, and unmatched indices.

| Field | Type | Description | Constraints |
|-------|------|-------------|-------------|
| partner | Vec\<isize\> | Per-node classification | -2=unmatched, -1=singleton, >0=partner index |
| old_to_new | Vec\<usize\> | Original index → condensed super-node index | Many-to-one for 2-cycle pairs |
| new_to_old | Vec\<usize\> | Condensed super-node index → first original index | One-to-one; length = condensed_dim |
| condensed_dim | usize | Number of condensed super-nodes (matched only) | condensed_dim <= n |
| singletons | usize | Count of singletons | |
| two_cycles | usize | Count of 2-cycle pairs | |

**Relationships**: Derived from Mc64Result.matching and Mc64Result.is_matched. The `is_matched` field distinguishes true singletons (matched, fwd[i]==i) from unmatched indices. Consumed by condensed graph construction and ordering expansion.

**State transitions**: None (computed once, used twice: for graph condensation and ordering expansion).

### CondensedGraph (internal intermediate)

Compressed sparse adjacency graph with one node per singleton/2-cycle super-node.

| Field | Type | Description | Constraints |
|-------|------|-------------|-------------|
| xadj | Vec\<i32\> | CSR row pointers (METIS convention) | Length = condensed_dim + 1; 0-based |
| adjncy | Vec\<i32\> | CSR adjacency indices (METIS convention) | No self-loops; no duplicates; symmetric |
| dim | usize | Graph dimension (= condensed_dim) | Fits in i32 |

**Relationships**: Built from the original matrix sparsity pattern + CycleDecomposition mappings. Consumed by METIS_NodeND.

**Validation rules**:
- No self-loops (diagonal entries excluded)
- Symmetric: if (i,j) is an edge, then (j,i) is an edge
- No duplicate edges (enforced by marker-array deduplication)
- All indices in [0, condensed_dim)

## Data Flow

```
Input Matrix (SparseColMat<usize, f64>)
    │
    ▼
MC64 Matching ──► Mc64Result { matching, scaling, matched, is_matched }
    │
    ▼
Cycle Splitting ──► CycleDecomposition { partner, old_to_new, new_to_old, ... }
  (uses matching.fwd + is_matched to distinguish singletons from unmatched)
    │
    ├──► Condensed Graph ──► CondensedGraph { col_ptr, adj, dim }
    │                              │
    │                              ▼
    │                        METIS_NodeND ──► condensed_order: Vec<usize>
    │                              │
    └──────────────────────────────┘
                    │
                    ▼
            Ordering Expansion ──► final_ordering: Perm<usize>
                    │
                    ▼
            MatchOrderResult { ordering, scaling, matched, condensed_dim, ... }
                    │
                    ▼
            AptpSymbolic::analyze(matrix, SymmetricOrdering::Custom(ordering))
```
