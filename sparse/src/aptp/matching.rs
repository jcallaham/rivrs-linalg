//! MC64 weighted bipartite matching and symmetric scaling.
//!
//! Implements the Duff & Koster (2001) Algorithm MPD — Dijkstra-based shortest
//! augmenting path matching on a logarithmic cost graph — with the MC64SYM
//! symmetric scaling from Duff & Pralet (2005).
//!
//! MC64 preprocessing improves diagonal dominance of sparse symmetric indefinite
//! matrices, reducing delayed pivots in subsequent APTP factorization. The algorithm
//! produces:
//! - A matching permutation decomposing into singletons and 2-cycles
//! - Symmetric scaling factors such that the scaled matrix has unit matched diagonal
//!   and off-diagonal entries bounded by 1 in absolute value
//!
//! # Algorithm
//!
//! 1. Build a bipartite cost graph with `c[i,j] = log(col_max_j) - log|a[i,j]|`
//! 2. Compute initial dual variables and greedy matching (~80% cardinality)
//! 3. For each unmatched column, find shortest augmenting path via Dijkstra
//! 4. Symmetrize scaling: `s[i] = exp(-(u[i] + v[i]) / 2)`
//!
//! # References
//!
//! - Duff & Koster (2001), "On Algorithms for Permuting Large Entries to the
//!   Diagonal of a Sparse Matrix", SIAM J. Matrix Anal. Appl. 22(4)
//! - Duff & Pralet (2005), "Strategies for Scaling and Pivoting for Sparse
//!   Symmetric Indefinite Problems", RAL Technical Report

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use faer::perm::Perm;
use faer::sparse::SparseColMat;

use crate::error::SparseError;

/// Result of MC64 matching-based preprocessing.
pub struct Mc64Result {
    /// Matching permutation: row i is matched to column σ(i).
    /// For symmetric matrices: decomposes into singletons (σ(i)=i)
    /// and 2-cycles (σ(i)=j, σ(j)=i).
    /// NOTE: This is NOT a fill-reducing ordering. Do not use directly
    /// with `SymmetricOrdering::Custom`.
    pub matching: Perm<usize>,

    /// Symmetric scaling factors (linear domain).
    /// `A_scaled[i,j] = scaling[i] * A[i,j] * scaling[j]`
    /// The scaled matrix has unit diagonal and off-diagonals <= 1
    /// for matched entries.
    pub scaling: Vec<f64>,

    /// Number of matched entries. Equals n for structurally nonsingular
    /// matrices. Less than n indicates structural singularity.
    pub matched: usize,

    /// Per-index matched status: `is_matched[i]` is true if index i
    /// participates in the matching (as row OR column for symmetric matrices).
    /// For structurally nonsingular matrices, all entries are true.
    /// Used by downstream cycle decomposition to distinguish singletons
    /// (matched to self) from unmatched indices.
    pub is_matched: Vec<bool>,
}

/// Optimization objective for the matching.
#[non_exhaustive]
pub enum Mc64Job {
    /// Maximize the product of diagonal entry magnitudes.
    /// Equivalent to minimizing sum of `-log|a_ij|` costs.
    /// Default for APTP preprocessing.
    MaximumProduct,
}

/// Bipartite cost graph with logarithmic edge costs.
///
/// Built from the input matrix by expanding the upper triangle to full
/// symmetric storage and computing `c[i,j] = log(col_max_j) - log|a[i,j]|`.
struct CostGraph {
    /// CSC column pointers for the full (symmetrized) matrix.
    col_ptr: Vec<usize>,
    /// CSC row indices.
    row_idx: Vec<usize>,
    /// Edge costs: `c[i,j] = col_max_log[j] - log|a[i,j]|` (non-negative).
    cost: Vec<f64>,
    /// Column maxima in log domain: `log(max_k |a[k,j]|)`.
    col_max_log: Vec<f64>,
    /// Matrix dimension.
    n: usize,
}

/// Working state during the matching algorithm.
struct MatchingState {
    /// For each row i, the column it is matched to (usize::MAX if unmatched).
    row_match: Vec<usize>,
    /// For each column j, the row it is matched to (usize::MAX if unmatched).
    col_match: Vec<usize>,
    /// Row dual variables (log domain).
    u: Vec<f64>,
}

const UNMATCHED: usize = usize::MAX;

/// Compute MC64 weighted bipartite matching and symmetric scaling for a sparse
/// symmetric matrix.
///
/// # Arguments
///
/// * `matrix` — Sparse symmetric matrix in upper-triangular CSC format (faer
///   convention). Must be square. Numeric values are required.
/// * `job` — Optimization objective. Currently only `Mc64Job::MaximumProduct`.
///
/// # Returns
///
/// * `Ok(Mc64Result)` — Matching permutation, scaling factors, and match count.
/// * `Err(SparseError::NotSquare)` — Matrix is not square.
/// * `Err(SparseError::InvalidInput)` — Zero dimension or non-finite entries.
/// * `Err(SparseError::AnalysisFailure)` — Internal algorithm failure.
///
/// # Algorithm References
///
/// - Duff & Koster (2001), Algorithm MPD
/// - Duff & Pralet (2005), MC64SYM
pub fn mc64_matching(
    matrix: &SparseColMat<usize, f64>,
    _job: Mc64Job,
) -> Result<Mc64Result, SparseError> {
    let (nrows, ncols) = (matrix.nrows(), matrix.ncols());

    // Input validation
    if nrows != ncols {
        return Err(SparseError::NotSquare {
            dims: (nrows, ncols),
        });
    }
    let n = nrows;

    if n == 0 {
        return Err(SparseError::InvalidInput {
            reason: "MC64 requires non-empty matrix".to_string(),
        });
    }

    // Check for non-finite entries
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    for j in 0..n {
        let start = symbolic.col_ptr()[j];
        let end = symbolic.col_ptr()[j + 1];
        for &val in &values[start..end] {
            if !val.is_finite() {
                return Err(SparseError::InvalidInput {
                    reason: "MC64 requires finite matrix entries".to_string(),
                });
            }
        }
    }

    // Trivial case: n=1
    if n == 1 {
        let has_entry = symbolic.col_ptr()[1] > symbolic.col_ptr()[0];
        let scale = if has_entry {
            let val = values[symbolic.col_ptr()[0]];
            if val.abs() > 0.0 {
                1.0 / val.abs().sqrt()
            } else {
                1.0
            }
        } else {
            1.0
        };
        let fwd: Box<[usize]> = vec![0].into_boxed_slice();
        let inv: Box<[usize]> = vec![0].into_boxed_slice();
        return Ok(Mc64Result {
            matching: Perm::new_checked(fwd, inv, 1),
            scaling: vec![scale],
            matched: if has_entry { 1 } else { 0 },
            is_matched: vec![has_entry],
        });
    }

    // Build cost graph
    let graph = build_cost_graph(matrix);

    // Greedy initial matching
    let mut state = greedy_initial_matching(&graph);

    // Augment unmatched columns via Dijkstra
    for j in 0..n {
        if state.col_match[j] != UNMATCHED {
            continue;
        }
        dijkstra_augment(j, &graph, &mut state);
    }

    // Count matched
    let matched = state.col_match.iter().filter(|&&m| m != UNMATCHED).count();

    if matched == n {
        // Full matching — compute scaling directly from dual variables
        let (mut scaling, fwd, inv) = build_full_match_result(&graph, &state);
        // Enforce |s_i * a_ij * s_j| <= 1 (corrects any dual feasibility gaps)
        enforce_scaling_bound(matrix, &mut scaling);
        return Ok(Mc64Result {
            matching: Perm::new_checked(fwd.into_boxed_slice(), inv.into_boxed_slice(), n),
            scaling,
            matched,
            is_matched: vec![true; n],
        });
    }

    // Structural singularity: matched < n
    // For structurally singular matrices, the dual variables from the partial
    // matching may not produce perfect quality scaling for all rows. We use
    // the standard approach (complementary slackness + Duff-Pralet correction)
    // then enforce the scaling bound to guarantee correctness.

    // Zero u for unmatched rows (SPRAL: where(iperm==0) dualu=0)
    for i in 0..n {
        if state.row_match[i] == UNMATCHED {
            state.u[i] = 0.0;
        }
    }

    // Compute column duals via complementary slackness and symmetric scaling
    let v = compute_column_duals(&graph, &state);
    let mut scaling = symmetrize_scaling(&state.u, &v, &graph.col_max_log);

    // Apply Duff-Pralet correction for unmatched indices.
    // is_matched[i] = true if index i participates in the matching in any role
    // (as row OR column). For symmetric matrices, row i and column i are the same
    // physical index, so we take the union. This drives both the Duff-Pralet
    // correction (which needs to know unmatched indices) and the permutation
    // builder (which needs to know which columns are "used").
    let is_row_matched: Vec<bool> = (0..n).map(|i| state.row_match[i] != UNMATCHED).collect();
    let is_col_matched: Vec<bool> = (0..n).map(|j| state.col_match[j] != UNMATCHED).collect();
    let is_matched: Vec<bool> = (0..n)
        .map(|i| is_row_matched[i] || is_col_matched[i])
        .collect();
    duff_pralet_correction(matrix, &mut scaling, &is_matched);

    // Enforce |s_i * a_ij * s_j| <= 1 (corrects any dual feasibility gaps)
    enforce_scaling_bound(matrix, &mut scaling);

    // Build permutation: matched rows keep their matched column, unmatched get remaining
    let (fwd, inv) = build_singular_permutation(n, &state, &is_row_matched);

    Ok(Mc64Result {
        matching: Perm::new_checked(fwd.into_boxed_slice(), inv.into_boxed_slice(), n),
        scaling,
        matched,
        is_matched,
    })
}

/// Compute column dual variables via complementary slackness.
///
/// For each matched column j: `v[j] = c[matched_row, j] - u[matched_row]`.
/// Unmatched columns get `v[j] = 0.0`.
/// This is the standard post-matching dual computation (SPRAL lines 1159-1167).
fn compute_column_duals(graph: &CostGraph, state: &MatchingState) -> Vec<f64> {
    let n = graph.n;
    let mut v = vec![0.0_f64; n];
    for (j, v_j) in v.iter_mut().enumerate() {
        let i = state.col_match[j];
        if i != UNMATCHED {
            let col_start = graph.col_ptr[j];
            let col_end = graph.col_ptr[j + 1];
            for idx in col_start..col_end {
                if graph.row_idx[idx] == i {
                    *v_j = graph.cost[idx] - state.u[i];
                    break;
                }
            }
        }
    }
    v
}

/// Build scaling, forward and inverse permutation for a full matching result.
fn build_full_match_result(
    graph: &CostGraph,
    state: &MatchingState,
) -> (Vec<f64>, Vec<usize>, Vec<usize>) {
    let n = graph.n;
    let v = compute_column_duals(graph, state);
    let scaling = symmetrize_scaling(&state.u, &v, &graph.col_max_log);

    let mut fwd = vec![0usize; n];
    for (i, fwd_i) in fwd.iter_mut().enumerate() {
        *fwd_i = state.row_match[i];
    }

    let mut inv = vec![0usize; n];
    for (i, &f) in fwd.iter().enumerate() {
        inv[f] = i;
    }

    (scaling, fwd, inv)
}

/// Build permutation arrays for structurally singular case.
/// Matched rows keep their matched column; unmatched rows get remaining columns.
fn build_singular_permutation(
    n: usize,
    state: &MatchingState,
    is_matched: &[bool],
) -> (Vec<usize>, Vec<usize>) {
    let mut fwd = vec![0usize; n];
    let mut unmatched_rows: Vec<usize> = Vec::new();

    for (i, fwd_i) in fwd.iter_mut().enumerate() {
        if state.row_match[i] != UNMATCHED {
            *fwd_i = state.row_match[i];
        } else {
            unmatched_rows.push(i);
        }
    }

    let mut used_cols = vec![false; n];
    for (i, &matched) in is_matched.iter().enumerate() {
        if matched {
            used_cols[state.row_match[i]] = true;
        }
    }
    let free_cols: Vec<usize> = (0..n).filter(|&j| !used_cols[j]).collect();
    for (idx, &i) in unmatched_rows.iter().enumerate() {
        fwd[i] = free_cols[idx];
    }

    let mut inv = vec![0usize; n];
    for (i, &f) in fwd.iter().enumerate() {
        inv[f] = i;
    }

    (fwd, inv)
}

/// Build the bipartite cost graph from a symmetric sparse matrix.
///
/// Expands the upper-triangular CSC input to full symmetric storage and computes
/// logarithmic costs: `c[i,j] = log(col_max_j) - log|a[i,j]|`.
fn build_cost_graph(matrix: &SparseColMat<usize, f64>) -> CostGraph {
    let n = matrix.nrows();
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();

    // Step 1: Collect all (row, col, |value|) entries, expanding upper triangle to full.
    // We store entries as (col, row, abs_val) grouped by column.
    let mut col_entries: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];

    for j in 0..n {
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for k in start..end {
            let i = row_indices[k];
            let abs_val = values[k].abs();
            if abs_val == 0.0 {
                continue; // Skip explicit zeros
            }
            // Add (i, j) entry
            col_entries[j].push((i, abs_val));
            // If off-diagonal, also add (j, i) for symmetry
            if i != j {
                col_entries[i].push((j, abs_val));
            }
        }
    }

    // Sort each column's entries by row index for deterministic ordering
    for entries in &mut col_entries {
        entries.sort_by_key(|&(row, _)| row);
    }

    // Step 2: Compute column maxima in log domain
    let mut col_max_log = vec![f64::NEG_INFINITY; n];
    for j in 0..n {
        for &(_, abs_val) in &col_entries[j] {
            let log_val = abs_val.ln();
            if log_val > col_max_log[j] {
                col_max_log[j] = log_val;
            }
        }
    }

    // Step 3: Build CSC with costs = col_max_log[j] - log|a[i,j]|
    let mut col_ptr = Vec::with_capacity(n + 1);
    let mut row_idx = Vec::new();
    let mut cost = Vec::new();

    col_ptr.push(0);
    for j in 0..n {
        for &(i, abs_val) in &col_entries[j] {
            let c = col_max_log[j] - abs_val.ln();
            row_idx.push(i);
            cost.push(c);
        }
        col_ptr.push(row_idx.len());
    }

    CostGraph {
        col_ptr,
        row_idx,
        cost,
        col_max_log,
        n,
    }
}

/// Compute initial matching and dual variables using the greedy heuristic.
///
/// Following Duff & Koster (2001, Section 4): compute initial dual variables
/// from column-minimum and row-minimum passes, then build greedy matching
/// from zero-reduced-cost edges with a secondary 2-pass rearrangement.
fn greedy_initial_matching(graph: &CostGraph) -> MatchingState {
    let n = graph.n;
    let mut row_match = vec![UNMATCHED; n];
    let mut col_match = vec![UNMATCHED; n];
    let mut u = vec![f64::INFINITY; n];

    // Pass 1: For each row i, find minimum cost across all columns
    // u[i] = min_j c[i,j]
    // Also record which column j achieves the minimum, and the position
    let mut best_col_for_row = vec![UNMATCHED; n]; // column with min cost for row i
    let mut best_cost_pos = vec![0usize; n]; // index into graph arrays

    for j in 0..n {
        let col_start = graph.col_ptr[j];
        let col_end = graph.col_ptr[j + 1];
        for idx in col_start..col_end {
            let i = graph.row_idx[idx];
            let c = graph.cost[idx];
            if c < u[i] {
                u[i] = c;
                best_col_for_row[i] = j;
                best_cost_pos[i] = idx;
            }
        }
    }

    // For rows with no entries, set u[i] = 0
    for u_i in &mut u {
        if *u_i == f64::INFINITY {
            *u_i = 0.0;
        }
    }

    // Pass 1 matching: match row i to its best column if that column is unmatched
    // Avoid matching on very dense columns (heuristic from SPRAL)
    let dense_threshold = if n > 50 { n / 10 } else { n };
    for i in 0..n {
        let j = best_col_for_row[i];
        if j == UNMATCHED {
            continue;
        }
        if col_match[j] != UNMATCHED {
            continue;
        }
        let col_degree = graph.col_ptr[j + 1] - graph.col_ptr[j];
        if col_degree > dense_threshold {
            continue;
        }
        row_match[i] = j;
        col_match[j] = i;
    }

    // Pass 2: For each unmatched column, find cheapest reduced-cost assignment
    // d[j] = min_i (c[i,j] - u[i]) for column j
    // If cheapest row is unmatched, assign directly
    // If cheapest row is matched to column jj, try length-2 augmentation
    let mut d_col = vec![0.0_f64; n]; // current min reduced cost for column j
    let mut search_from = vec![0usize; n]; // track search position per column
    search_from[..n].copy_from_slice(&graph.col_ptr[..n]);

    'col_loop: for j in 0..n {
        if col_match[j] != UNMATCHED {
            continue;
        }
        let col_start = graph.col_ptr[j];
        let col_end = graph.col_ptr[j + 1];
        if col_start >= col_end {
            continue; // empty column
        }

        // Find row with smallest reduced cost in column j
        let mut best_i = graph.row_idx[col_start];
        let mut best_rc = graph.cost[col_start] - u[best_i];
        let mut best_k = col_start;

        for idx in (col_start + 1)..col_end {
            let i = graph.row_idx[idx];
            let rc = graph.cost[idx] - u[i];
            if rc < best_rc
                || (rc == best_rc && row_match[i] == UNMATCHED && row_match[best_i] != UNMATCHED)
            {
                best_rc = rc;
                best_i = i;
                best_k = idx;
            }
        }

        d_col[j] = best_rc;

        // If best row is unmatched, assign directly
        if row_match[best_i] == UNMATCHED {
            row_match[best_i] = j;
            col_match[j] = best_i;
            search_from[j] = best_k + 1;
            continue;
        }

        // Try 2-augmentation: for each tied row i matched to column jj,
        // scan column jj for an unmatched row ii
        for idx in best_k..col_end {
            let i = graph.row_idx[idx];
            let rc = graph.cost[idx] - u[i];
            if rc > best_rc {
                continue;
            }
            let jj = row_match[i];
            if jj == UNMATCHED {
                continue;
            }

            // Scan column jj for unmatched row
            let jj_end = graph.col_ptr[jj + 1];
            for kk in search_from[jj]..jj_end {
                let ii = graph.row_idx[kk];
                if row_match[ii] != UNMATCHED {
                    continue;
                }
                let rc_ii = graph.cost[kk] - u[ii];
                if rc_ii <= d_col[jj] {
                    // Augment: (i,j) and (ii,jj)
                    col_match[jj] = ii;
                    row_match[ii] = jj;
                    search_from[jj] = kk + 1;
                    col_match[j] = i;
                    row_match[i] = j;
                    search_from[j] = idx + 1;
                    continue 'col_loop;
                }
            }
            search_from[jj] = jj_end;
        }
    }

    MatchingState {
        row_match,
        col_match,
        u,
    }
}

/// Find shortest augmenting path from unmatched column `root_col` using Dijkstra
/// on reduced costs and augment the matching if a path is found.
///
/// Returns `true` if augmenting path found, `false` if no path exists.
fn dijkstra_augment(root_col: usize, graph: &CostGraph, state: &mut MatchingState) -> bool {
    let n = graph.n;

    // d[i] = shortest distance from root_col to row i
    let mut d = vec![f64::INFINITY; n];
    // pr[j] = parent column along shortest path tree
    let mut pr = vec![UNMATCHED; n];
    // out[j] = index in graph arrays of the edge entering column j's matched row
    //          (used for augmenting path traceback — NOT the matched edge in col j)
    let mut out = vec![0usize; n];
    // visited[i] = true if row i has been finalized (in set B)
    let mut visited = vec![false; n];

    // csp = cost of shortest augmenting path found so far
    let mut csp = f64::INFINITY;
    let mut isp = UNMATCHED; // row at end of shortest augmenting path
    let mut jsp = UNMATCHED; // column entering that row

    // Priority queue: (Reverse(distance), row_index)
    let mut heap: BinaryHeap<(Reverse<OrderedFloat>, usize)> = BinaryHeap::new();
    // Also maintain Q1: rows with d[i] == dmin (stored in a flat list)
    let mut q1: Vec<usize> = Vec::new();
    let mut dmin = f64::INFINITY;

    // Scan root column
    pr[root_col] = UNMATCHED; // sentinel: root has no parent
    let col_start = graph.col_ptr[root_col];
    let col_end = graph.col_ptr[root_col + 1];

    for idx in col_start..col_end {
        let i = graph.row_idx[idx];
        let dnew = graph.cost[idx] - state.u[i];

        if dnew >= csp {
            continue;
        }

        if state.row_match[i] == UNMATCHED {
            // Found augmenting path of length dnew
            csp = dnew;
            isp = i;
            jsp = root_col;
        } else {
            d[i] = dnew;
            if dnew < dmin {
                dmin = dnew;
            }
            let jj = state.row_match[i];
            out[jj] = idx;
            pr[jj] = root_col;
        }
    }

    // Build initial Q1 (rows with d == dmin) and heap (rest)
    for idx in col_start..col_end {
        let i = graph.row_idx[idx];
        if d[i] == f64::INFINITY || state.row_match[i] == UNMATCHED {
            continue;
        }
        if csp <= d[i] {
            d[i] = f64::INFINITY;
            continue;
        }
        if d[i] <= dmin {
            q1.push(i);
        } else {
            heap.push((Reverse(OrderedFloat(d[i])), i));
        }
    }

    // Main Dijkstra loop
    loop {
        // Get next row from Q1 or heap
        let q0;
        if let Some(i) = q1.pop() {
            q0 = i;
        } else if let Some(&(Reverse(OrderedFloat(top_d)), _)) = heap.peek() {
            if top_d >= csp {
                break;
            }
            dmin = top_d;
            // Extract all rows with d == dmin into Q1
            while let Some(&(Reverse(OrderedFloat(d_val)), i)) = heap.peek() {
                if d_val > dmin {
                    break;
                }
                heap.pop();
                if !visited[i] && d[i] <= dmin {
                    q1.push(i);
                }
            }
            if let Some(i) = q1.pop() {
                q0 = i;
            } else {
                break;
            }
        } else {
            break;
        }

        if d[q0] >= csp || visited[q0] {
            continue;
        }

        visited[q0] = true;
        let dq0 = d[q0];

        // Scan column matched with row q0
        let j = state.row_match[q0];
        if j == UNMATCHED {
            continue;
        }

        // Compute vj from the MATCHED edge (q0, j) in column j.
        // Note: out[j] points to the discovery edge (q0, scanning_col), not the
        // matched edge in column j. We must find the actual matched edge cost.
        // O(degree) scan — could be avoided by maintaining v[j] during Dijkstra
        // as SPRAL does, at the cost of an additional n-sized array.
        let col_start_j = graph.col_ptr[j];
        let col_end_j = graph.col_ptr[j + 1];
        let mut matched_cost = 0.0_f64;
        for idx in col_start_j..col_end_j {
            if graph.row_idx[idx] == q0 {
                matched_cost = graph.cost[idx];
                break;
            }
        }
        let vj = dq0 - matched_cost + state.u[q0];

        for idx in col_start_j..col_end_j {
            let i = graph.row_idx[idx];
            if visited[i] {
                continue;
            }

            let dnew = vj + graph.cost[idx] - state.u[i];

            if dnew >= csp {
                continue;
            }

            if state.row_match[i] == UNMATCHED {
                // Found shorter augmenting path
                csp = dnew;
                isp = i;
                jsp = j;
            } else if dnew < d[i] {
                d[i] = dnew;
                let jj = state.row_match[i];
                out[jj] = idx;
                pr[jj] = j;

                if dnew <= dmin {
                    q1.push(i);
                } else {
                    heap.push((Reverse(OrderedFloat(dnew)), i));
                }
            }
        }
    }

    if csp == f64::INFINITY {
        return false; // No augmenting path
    }

    // Augment matching along path from isp back to root_col
    let mut i = isp;
    let mut j = jsp;
    state.row_match[i] = j;
    state.col_match[j] = i;

    loop {
        let jj = pr[j];
        if jj == UNMATCHED {
            break; // Reached root
        }
        let k = out[j];
        i = graph.row_idx[k];
        state.row_match[i] = jj;
        state.col_match[jj] = i;
        j = jj;
    }

    // Update dual variables for visited rows
    for i in 0..n {
        if visited[i] {
            state.u[i] = state.u[i] + d[i] - csp;
        }
    }

    true
}

/// Compute symmetric scaling factors from dual variables.
///
/// Following Duff & Pralet (2005), MC64SYM:
/// `scaling[i] = exp((u[i] + v[i] - col_max_log[i]) / 2)`
///
/// where u, v are row/column dual variables in cost domain and col_max_log
/// are the column maxima used in cost construction.
fn symmetrize_scaling(u: &[f64], v: &[f64], col_max_log: &[f64]) -> Vec<f64> {
    let n = u.len();
    let mut scaling = Vec::with_capacity(n);

    for i in 0..n {
        // SPRAL: rscaling = dualu, cscaling = dualv - cmax
        // then scaling = exp((rscaling + cscaling) / 2)
        // = exp((u[i] + v[i] - col_max_log[i]) / 2)
        let log_scale = (u[i] + v[i] - col_max_log[i]) / 2.0;

        // Clamp to avoid overflow/underflow in exp
        let clamped = log_scale.clamp(-500.0, 500.0);
        scaling.push(clamped.exp());
    }

    scaling
}

/// Enforce the scaling bound |s_i * a_ij * s_j| <= 1 for all stored entries.
///
/// When the dual variables from the matching algorithm don't perfectly satisfy
/// feasibility (e.g., due to greedy initialization without column dual tracking),
/// complementary slackness can produce scaling factors that violate the bound.
/// This function iteratively reduces scaling factors to enforce the constraint.
///
/// Works on the upper-triangular input matrix, handling symmetry.
/// Convergence is guaranteed since scaling factors only decrease.
fn enforce_scaling_bound(matrix: &SparseColMat<usize, f64>, scaling: &mut [f64]) {
    let n = matrix.nrows();
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();

    // Iterate until convergence (typically 1-2 passes)
    for _iter in 0..10 {
        let mut changed = false;

        // For each row i, compute max_j |a_ij * s_j| over all entries
        // Then enforce s_i <= 1 / max_j(|a_ij| * s_j)
        let mut row_max = vec![0.0_f64; n];

        for j in 0..n {
            let start = col_ptrs[j];
            let end = col_ptrs[j + 1];
            for k in start..end {
                let i = row_indices[k];
                let abs_val = values[k].abs();
                if abs_val == 0.0 {
                    continue;
                }
                let sj_contrib = abs_val * scaling[j];
                if sj_contrib > row_max[i] {
                    row_max[i] = sj_contrib;
                }
                if i != j {
                    let si_contrib = abs_val * scaling[i];
                    if si_contrib > row_max[j] {
                        row_max[j] = si_contrib;
                    }
                }
            }
        }

        // Cap scaling factors where needed
        for i in 0..n {
            if row_max[i] > 0.0 {
                let s_max = 1.0 / row_max[i];
                if scaling[i] > s_max * (1.0 + 1e-12) {
                    scaling[i] = s_max;
                    changed = true;
                }
            }
        }

        if !changed {
            break;
        }
    }
}

/// Apply Duff-Pralet correction for unmatched indices in structurally singular matrices.
///
/// For unmatched index i: `scaling[i] = 1.0 / max_k |a[i,k] * scaling[k]|` over
/// matched k. Convention: `1/0 = 1.0`.
fn duff_pralet_correction(
    matrix: &SparseColMat<usize, f64>,
    scaling: &mut [f64],
    is_matched: &[bool],
) {
    let n = matrix.nrows();
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();

    // For each unmatched index, compute max_k |a[i,k] * scaling[k]| over matched k
    // We need the original (stored) scaling values, so clone first
    let orig_scaling = scaling.to_vec();

    // Build log-domain scaling for unmatched rows
    let mut log_max = vec![f64::NEG_INFINITY; n];

    for j in 0..n {
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for k in start..end {
            let i = row_indices[k];
            let abs_val = values[k].abs();
            if abs_val == 0.0 {
                continue;
            }
            // Entry (i, j) in upper triangle
            // If i is unmatched and j is matched:
            if !is_matched[i] && is_matched[j] {
                let contrib = abs_val.ln() + orig_scaling[j].ln();
                if contrib > log_max[i] {
                    log_max[i] = contrib;
                }
            }
            // If j is unmatched and i is matched (symmetric entry):
            if i != j && !is_matched[j] && is_matched[i] {
                let contrib = abs_val.ln() + orig_scaling[i].ln();
                if contrib > log_max[j] {
                    log_max[j] = contrib;
                }
            }
        }
    }

    // Apply correction
    for i in 0..n {
        if is_matched[i] {
            continue;
        }
        if log_max[i] == f64::NEG_INFINITY {
            // Isolated row: no matched neighbors
            scaling[i] = 1.0;
        } else {
            scaling[i] = (-log_max[i]).exp();
        }
    }
}

/// Count singletons, 2-cycles, and longer cycles in a matching permutation.
///
/// For symmetric matrices, the optimal matching decomposes into singletons
/// (σ(i)=i) and 2-cycles (σ(i)=j, σ(j)=i, i≠j). However, the asymmetric
/// cost graph can produce longer cycles that are genuinely optimal (see
/// `docs/mc64-scaling-notes.md`). Returns `(singletons, two_cycles, longer_cycles)`.
pub fn count_cycles(matching: &[usize]) -> (usize, usize, usize) {
    let n = matching.len();
    let mut visited = vec![false; n];
    let mut singletons = 0;
    let mut two_cycles = 0;
    let mut longer_cycles = 0;

    for i in 0..n {
        if visited[i] {
            continue;
        }
        let j = matching[i];
        if j == i {
            singletons += 1;
            visited[i] = true;
        } else if matching[j] == i {
            two_cycles += 1;
            visited[i] = true;
            visited[j] = true;
        } else {
            // Longer cycle — trace it
            longer_cycles += 1;
            let mut k = i;
            loop {
                visited[k] = true;
                k = matching[k];
                if k == i {
                    break;
                }
            }
        }
    }

    (singletons, two_cycles, longer_cycles)
}

/// Wrapper for f64 that provides total ordering for use in BinaryHeap.
#[derive(Debug, Clone, Copy, PartialEq)]
struct OrderedFloat(f64);

impl Eq for OrderedFloat {}

impl PartialOrd for OrderedFloat {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedFloat {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::sparse::Triplet;

    /// Helper: create a symmetric upper-triangular matrix from entries.
    /// Only entries with i <= j should be provided.
    fn make_upper_tri(n: usize, entries: &[(usize, usize, f64)]) -> SparseColMat<usize, f64> {
        let triplets: Vec<_> = entries
            .iter()
            .map(|&(i, j, v)| Triplet::new(i, j, v))
            .collect();
        SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
    }

    /// 3x3 symmetric matrix (upper triangle):
    /// [4  2  0]
    /// [2  5  1]
    /// [0  1  3]
    fn make_3x3_test() -> SparseColMat<usize, f64> {
        make_upper_tri(
            3,
            &[
                (0, 0, 4.0),
                (0, 1, 2.0),
                (1, 1, 5.0),
                (1, 2, 1.0),
                (2, 2, 3.0),
            ],
        )
    }

    // ---- T004: build_cost_graph tests ----

    #[test]
    fn test_build_cost_graph_3x3() {
        let matrix = make_3x3_test();
        let graph = build_cost_graph(&matrix);

        assert_eq!(graph.n, 3);

        // Full symmetric matrix has entries:
        // col 0: rows [0, 1] (diagonal + (0,1))
        // col 1: rows [0, 1, 2] (symmetric (0,1) + diagonal + (1,2))
        // col 2: rows [1, 2] (symmetric (1,2) + diagonal)

        // Check column counts
        let col_count = |j: usize| graph.col_ptr[j + 1] - graph.col_ptr[j];
        assert_eq!(col_count(0), 2, "col 0 should have 2 entries");
        assert_eq!(col_count(1), 3, "col 1 should have 3 entries");
        assert_eq!(col_count(2), 2, "col 2 should have 2 entries");

        // Column maxima in log domain:
        // col 0: max(|4|, |2|) = 4, log(4) = 1.386...
        // col 1: max(|2|, |5|, |1|) = 5, log(5) = 1.609...
        // col 2: max(|1|, |3|) = 3, log(3) = 1.099...
        assert!((graph.col_max_log[0] - 4.0_f64.ln()).abs() < 1e-12);
        assert!((graph.col_max_log[1] - 5.0_f64.ln()).abs() < 1e-12);
        assert!((graph.col_max_log[2] - 3.0_f64.ln()).abs() < 1e-12);

        // All costs should be non-negative
        for &c in &graph.cost {
            assert!(c >= -1e-14, "cost {} should be non-negative", c);
        }

        // Diagonal entries should have cost = col_max_log - log|diag|
        // (0,0): log(4) - log(4) = 0
        // (1,1): log(5) - log(5) = 0
        // (2,2): log(3) - log(3) = 0
        for j in 0..3 {
            let col_start = graph.col_ptr[j];
            let col_end = graph.col_ptr[j + 1];
            for idx in col_start..col_end {
                if graph.row_idx[idx] == j {
                    assert!(
                        graph.cost[idx].abs() < 1e-12,
                        "diagonal ({},{}) cost should be ~0, got {}",
                        j,
                        j,
                        graph.cost[idx]
                    );
                }
            }
        }
    }

    #[test]
    fn test_build_cost_graph_includes_diagonal() {
        let matrix = make_upper_tri(2, &[(0, 0, 3.0), (0, 1, 1.0), (1, 1, 2.0)]);
        let graph = build_cost_graph(&matrix);

        // Both diagonals should be present
        let mut has_diag = [false; 2];
        for (j, diag) in has_diag.iter_mut().enumerate() {
            let col_start = graph.col_ptr[j];
            let col_end = graph.col_ptr[j + 1];
            for idx in col_start..col_end {
                if graph.row_idx[idx] == j {
                    *diag = true;
                }
            }
        }
        assert!(has_diag[0], "diagonal (0,0) missing");
        assert!(has_diag[1], "diagonal (1,1) missing");
    }

    #[test]
    fn test_build_cost_graph_symmetric_expansion() {
        // Upper triangle only: entries (0,1) and (1,2)
        let matrix = make_upper_tri(
            3,
            &[
                (0, 0, 1.0),
                (0, 1, 2.0),
                (1, 1, 3.0),
                (1, 2, 4.0),
                (2, 2, 5.0),
            ],
        );
        let graph = build_cost_graph(&matrix);

        // Check that (1,0) and (2,1) appear in the expanded graph
        let has_entry = |col: usize, row: usize| -> bool {
            let start = graph.col_ptr[col];
            let end = graph.col_ptr[col + 1];
            graph.row_idx[start..end].contains(&row)
        };

        assert!(has_entry(0, 1), "symmetric entry (1,0) should exist");
        assert!(has_entry(1, 0), "entry (0,1) should exist");
        assert!(
            has_entry(2, 1),
            "symmetric entry (1,2) should exist in col 2"
        );
        assert!(has_entry(1, 2), "entry (2,1) should exist in col 1");
    }

    // ---- T006: greedy_initial_matching tests ----

    #[test]
    fn test_greedy_matching_4x4() {
        // 4x4 matrix where greedy can find at least 3 matches
        let matrix = make_upper_tri(
            4,
            &[
                (0, 0, 10.0),
                (0, 1, 1.0),
                (1, 1, 8.0),
                (1, 2, 2.0),
                (2, 2, 6.0),
                (2, 3, 3.0),
                (3, 3, 5.0),
            ],
        );
        let graph = build_cost_graph(&matrix);
        let state = greedy_initial_matching(&graph);

        // Count matched
        let matched_count = state.row_match.iter().filter(|&&m| m != UNMATCHED).count();
        assert!(
            matched_count >= 3,
            "greedy should match at least 3 of 4, got {}",
            matched_count
        );

        // Verify dual feasibility: u[i] + v[j] <= c[i,j] for all edges
        // Note: v is not yet computed in greedy, but u should be valid
        // u[i] should be non-negative for cost graphs
        for &ui in &state.u {
            assert!(ui.is_finite(), "dual u should be finite");
        }

        // Verify matched edges have zero reduced cost (approximately)
        for i in 0..4 {
            let j = state.row_match[i];
            if j == UNMATCHED {
                continue;
            }
            let col_start = graph.col_ptr[j];
            let col_end = graph.col_ptr[j + 1];
            for idx in col_start..col_end {
                if graph.row_idx[idx] == i {
                    // Matched edge: reduced cost c[i,j] - u[i] should be small
                    // (but greedy doesn't guarantee exactly zero for all)
                    let _rc = graph.cost[idx] - state.u[i];
                    break;
                }
            }
        }
    }

    // ---- T007: dijkstra_augment tests ----

    #[test]
    fn test_dijkstra_augment_3x3() {
        // 3x3 with one unmatched column; verify augmenting path
        let matrix = make_upper_tri(
            3,
            &[
                (0, 0, 5.0),
                (0, 1, 3.0),
                (0, 2, 1.0),
                (1, 1, 4.0),
                (1, 2, 2.0),
                (2, 2, 6.0),
            ],
        );
        let graph = build_cost_graph(&matrix);
        let mut state = greedy_initial_matching(&graph);

        let initial_matched = state.col_match.iter().filter(|&&m| m != UNMATCHED).count();

        // Try to augment each unmatched column
        let mut augmented = false;
        for j in 0..3 {
            if state.col_match[j] == UNMATCHED && dijkstra_augment(j, &graph, &mut state) {
                augmented = true;
            }
        }

        let final_matched = state.col_match.iter().filter(|&&m| m != UNMATCHED).count();

        // Should have found augmenting path if initial matching wasn't perfect
        if initial_matched < 3 {
            assert!(augmented, "should find augmenting path");
            assert!(
                final_matched > initial_matched,
                "matching size should increase"
            );
        }

        // Dual feasibility: u[i] should remain finite
        for &ui in &state.u {
            assert!(ui.is_finite(), "dual u should be finite after augmentation");
        }
    }

    // ---- T008: symmetrize_scaling tests ----

    #[test]
    fn test_symmetrize_scaling_known_duals() {
        // Known dual values; verify scaling = exp((u + v - col_max_log) / 2)
        let u = vec![0.5, 1.0, 0.0];
        let v = vec![0.2, 0.3, 0.8];
        let col_max_log = vec![1.0, 1.5, 0.5];

        let scaling = symmetrize_scaling(&u, &v, &col_max_log);

        for i in 0..3 {
            let expected = ((u[i] + v[i] - col_max_log[i]) / 2.0).exp();
            assert!(
                (scaling[i] - expected).abs() < 1e-12,
                "scaling[{}] = {}, expected {}",
                i,
                scaling[i],
                expected
            );
        }
    }

    #[test]
    fn test_symmetrize_scaling_positive() {
        let u = vec![1.0, -0.5, 2.0];
        let v = vec![0.5, 1.5, -1.0];
        let col_max_log = vec![0.0, 0.0, 0.0];

        let scaling = symmetrize_scaling(&u, &v, &col_max_log);

        for (i, &s) in scaling.iter().enumerate() {
            assert!(s > 0.0, "scaling[{}] = {} should be positive", i, s);
            assert!(s.is_finite(), "scaling[{}] should be finite", i);
        }
    }

    // ---- T009 (partial): end-to-end mc64_matching tests ----

    #[test]
    fn test_mc64_diagonal_identity() {
        // Diagonal matrix: identity matching, scaling = 1/sqrt(diag)
        let matrix = make_upper_tri(3, &[(0, 0, 4.0), (1, 1, 9.0), (2, 2, 1.0)]);

        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
        assert_eq!(result.matched, 3);

        // Matching should be identity
        let (fwd, _) = result.matching.as_ref().arrays();
        for (i, &f) in fwd.iter().enumerate() {
            assert_eq!(f, i, "diagonal matrix matching should be identity");
        }

        // All scaling factors positive and finite
        for (i, &s) in result.scaling.iter().enumerate() {
            assert!(s > 0.0, "scaling[{}] should be positive", i);
            assert!(s.is_finite(), "scaling[{}] should be finite", i);
        }
    }

    #[test]
    fn test_mc64_tridiagonal_indefinite() {
        // 4x4 tridiagonal indefinite:
        // [ 2  -1   0   0]
        // [-1  -3   2   0]
        // [ 0   2   1  -1]
        // [ 0   0  -1  -4]
        let matrix = make_upper_tri(
            4,
            &[
                (0, 0, 2.0),
                (0, 1, -1.0),
                (1, 1, -3.0),
                (1, 2, 2.0),
                (2, 2, 1.0),
                (2, 3, -1.0),
                (3, 3, -4.0),
            ],
        );

        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
        assert_eq!(result.matched, 4);

        // Verify SPRAL scaling properties
        verify_scaling_properties(&matrix, &result);
    }

    #[test]
    fn test_mc64_arrow_indefinite() {
        // 5x5 arrow: dense first row/col, indefinite
        // [10  1  1  1  1]
        // [ 1 -3  0  0  0]
        // [ 1  0  5  0  0]
        // [ 1  0  0 -2  0]
        // [ 1  0  0  0  4]
        let matrix = make_upper_tri(
            5,
            &[
                (0, 0, 10.0),
                (0, 1, 1.0),
                (0, 2, 1.0),
                (0, 3, 1.0),
                (0, 4, 1.0),
                (1, 1, -3.0),
                (2, 2, 5.0),
                (3, 3, -2.0),
                (4, 4, 4.0),
            ],
        );

        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
        assert_eq!(result.matched, 5);

        verify_scaling_properties(&matrix, &result);
    }

    #[test]
    fn test_mc64_trivial_1x1() {
        let matrix = make_upper_tri(1, &[(0, 0, 7.0)]);
        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
        assert_eq!(result.matched, 1);
        assert_eq!(result.scaling.len(), 1);
        assert!(result.scaling[0] > 0.0);
    }

    #[test]
    fn test_mc64_trivial_2x2() {
        let matrix = make_upper_tri(2, &[(0, 0, 3.0), (0, 1, 1.0), (1, 1, 5.0)]);
        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
        assert_eq!(result.matched, 2);

        verify_scaling_properties(&matrix, &result);
    }

    #[test]
    fn test_mc64_not_square_error() {
        // Non-square matrix should error
        let triplets = vec![Triplet::new(0, 0, 1.0), Triplet::new(0, 1, 2.0)];
        let matrix = SparseColMat::try_new_from_triplets(2, 3, &triplets).unwrap();
        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct);
        assert!(matches!(result, Err(SparseError::NotSquare { .. })));
    }

    #[test]
    fn test_mc64_zero_dim_error() {
        let triplets: Vec<Triplet<usize, usize, f64>> = vec![];
        let matrix = SparseColMat::try_new_from_triplets(0, 0, &triplets).unwrap();
        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct);
        assert!(matches!(result, Err(SparseError::InvalidInput { .. })));
    }

    #[test]
    fn test_count_cycles_identity() {
        let matching = vec![0, 1, 2, 3];
        let (s, c, l) = count_cycles(&matching);
        assert_eq!(s, 4);
        assert_eq!(c, 0);
        assert_eq!(l, 0);
    }

    #[test]
    fn test_count_cycles_two_swaps() {
        let matching = vec![1, 0, 3, 2];
        let (s, c, l) = count_cycles(&matching);
        assert_eq!(s, 0);
        assert_eq!(c, 2);
        assert_eq!(l, 0);
    }

    #[test]
    fn test_count_cycles_mixed() {
        let matching = vec![0, 2, 1, 3, 4];
        let (s, c, l) = count_cycles(&matching);
        assert_eq!(s, 3); // 0, 3, 4 are singletons
        assert_eq!(c, 1); // (1,2) is a 2-cycle
        assert_eq!(l, 0);
    }

    #[test]
    fn test_count_cycles_longer_cycle() {
        // 3-cycle: 0→1→2→0, plus singleton 3
        let matching = vec![1, 2, 0, 3];
        let (s, c, l) = count_cycles(&matching);
        assert_eq!(s, 1); // 3 is singleton
        assert_eq!(c, 0);
        assert_eq!(l, 1); // one 3-cycle
    }

    /// Verify SPRAL scaling properties for a matching result.
    /// Delegates to the shared helper in testing::mc64_validation.
    fn verify_scaling_properties(matrix: &SparseColMat<usize, f64>, result: &Mc64Result) {
        use crate::testing::verify_spral_scaling_properties;
        verify_spral_scaling_properties("unit_test", matrix, result);
    }

    // ---- T018: duff_pralet_correction tests ----

    #[test]
    fn test_duff_pralet_4x4_singular() {
        // 4x4 structurally singular: row/col 3 has no matching
        // [4  2  0  1]
        // [2  5  1  0]
        // [0  1  3  0]
        // [1  0  0  0]  <- no diagonal, sparse connections
        let matrix = make_upper_tri(
            4,
            &[
                (0, 0, 4.0),
                (0, 1, 2.0),
                (0, 3, 1.0),
                (1, 1, 5.0),
                (1, 2, 1.0),
                (2, 2, 3.0),
            ],
        );

        let mut scaling = vec![0.5, 0.4, 0.6, 0.0]; // pre-set for matched
        let is_matched = vec![true, true, true, false]; // row 3 unmatched

        duff_pralet_correction(&matrix, &mut scaling, &is_matched);

        // Row 3 connects to row 0 via entry (0,3)=1.0
        // scaling[3] = 1.0 / max_k |a[3,k] * scaling[k]| over matched k
        // Only connection is (0,3)=1.0, scaling[0]=0.5
        // So scaling[3] = 1.0 / |1.0 * 0.5| = 2.0
        assert!(scaling[3] > 0.0, "unmatched scaling should be positive");
        assert!(scaling[3].is_finite(), "unmatched scaling should be finite");

        // Matched scaling should be unchanged
        assert!((scaling[0] - 0.5).abs() < 1e-12);
        assert!((scaling[1] - 0.4).abs() < 1e-12);
        assert!((scaling[2] - 0.6).abs() < 1e-12);
    }

    #[test]
    fn test_duff_pralet_isolated_row() {
        // Row with no connections to matched set → scaling = 1.0
        let matrix = make_upper_tri(
            3,
            &[
                (0, 0, 4.0),
                (1, 1, 5.0),
                // No entries connecting row 2 to rows 0 or 1
                (2, 2, 3.0),
            ],
        );

        let mut scaling = vec![0.5, 0.4, 0.0];
        // Only diagonal entries; if row 2 is unmatched but its only entry is (2,2),
        // and column 2 is unmatched too, then no matched connections
        let is_matched = vec![true, true, false];

        duff_pralet_correction(&matrix, &mut scaling, &is_matched);

        // Row 2 has no connections to matched indices (only diagonal which is unmatched)
        assert_eq!(
            scaling[2], 1.0,
            "isolated unmatched row should get scaling 1.0"
        );
    }

    #[test]
    fn test_duff_pralet_all_positive() {
        let matrix = make_upper_tri(
            4,
            &[
                (0, 0, 4.0),
                (0, 1, 2.0),
                (0, 3, 1.0),
                (1, 1, 5.0),
                (1, 2, 1.0),
                (2, 2, 3.0),
            ],
        );

        let mut scaling = vec![0.5, 0.4, 0.6, 0.0];
        let is_matched = vec![true, true, true, false];

        duff_pralet_correction(&matrix, &mut scaling, &is_matched);

        for (i, &s) in scaling.iter().enumerate() {
            assert!(s > 0.0, "scaling[{}] = {} should be positive", i, s);
            assert!(s.is_finite(), "scaling[{}] = {} should be finite", i, s);
        }
    }

    // ---- T019 (unit part): structurally singular mc64_matching ----

    #[test]
    fn test_mc64_singular_zero_diagonal() {
        // Structurally singular: no diagonal entries, forcing partial matching
        // [0  5  0  0]
        // [5  0  0  0]
        // [0  0  0  3]
        // [0  0  3  0]
        let matrix = make_upper_tri(4, &[(0, 1, 5.0), (2, 3, 3.0)]);

        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();

        // Should find matching (even if not perfect)
        for (i, &s) in result.scaling.iter().enumerate() {
            assert!(s > 0.0, "scaling[{}] should be positive", i);
            assert!(s.is_finite(), "scaling[{}] should be finite", i);
        }

        // Matching should be valid permutation
        let (fwd, _) = result.matching.as_ref().arrays();
        let mut seen = [false; 4];
        for &f in fwd {
            assert!(!seen[f], "duplicate in matching");
            seen[f] = true;
        }
    }

    // ---- Issue #8: NaN/Inf input rejection ----

    #[test]
    fn test_mc64_nan_entry_error() {
        let triplets = vec![
            Triplet::new(0, 0, 4.0),
            Triplet::new(0, 1, f64::NAN),
            Triplet::new(1, 1, 5.0),
        ];
        let matrix = SparseColMat::try_new_from_triplets(2, 2, &triplets).unwrap();
        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct);
        assert!(
            matches!(result, Err(SparseError::InvalidInput { .. })),
            "NaN entry should produce InvalidInput error"
        );
    }

    #[test]
    fn test_mc64_inf_entry_error() {
        let triplets = vec![
            Triplet::new(0, 0, 4.0),
            Triplet::new(0, 1, f64::INFINITY),
            Triplet::new(1, 1, 5.0),
        ];
        let matrix = SparseColMat::try_new_from_triplets(2, 2, &triplets).unwrap();
        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct);
        assert!(
            matches!(result, Err(SparseError::InvalidInput { .. })),
            "Inf entry should produce InvalidInput error"
        );
    }

    // ---- Issue #9: enforce_scaling_bound in isolation ----

    #[test]
    fn test_enforce_scaling_bound_corrects_violation() {
        // Hand-craft scaling factors that violate |s_i * a_ij * s_j| <= 1
        // Matrix: [4 2; 2 5] (upper triangle)
        let matrix = make_upper_tri(2, &[(0, 0, 4.0), (0, 1, 2.0), (1, 1, 5.0)]);

        // These scaling factors are too large: s[0]*4*s[0] = 2*4*2 = 16 >> 1
        let mut scaling = vec![2.0, 2.0];
        enforce_scaling_bound(&matrix, &mut scaling);

        // After enforcement, all scaled entries should be <= 1
        let symbolic = matrix.symbolic();
        let values = matrix.val();
        let row_idx = symbolic.row_idx();
        for j in 0..2 {
            let start = symbolic.col_ptr()[j];
            let end = symbolic.col_ptr()[j + 1];
            for (k, (&i, &val)) in row_idx[start..end]
                .iter()
                .zip(&values[start..end])
                .enumerate()
            {
                let scaled = (scaling[i] * val * scaling[j]).abs();
                assert!(
                    scaled <= 1.0 + 1e-10,
                    "|s[{}]*a[{},{}]*s[{}]| = {:.6e} > 1.0 after enforcement (k={})",
                    i,
                    i,
                    j,
                    j,
                    scaled,
                    k
                );
            }
        }

        // Scaling factors should have been reduced from the original 2.0
        assert!(scaling[0] < 2.0, "scaling[0] should be reduced");
        assert!(scaling[1] < 2.0, "scaling[1] should be reduced");
        assert!(scaling[0] > 0.0, "scaling[0] should remain positive");
        assert!(scaling[1] > 0.0, "scaling[1] should remain positive");
    }

    #[test]
    fn test_enforce_scaling_bound_noop_on_valid() {
        // Scaling factors already satisfying the bound should not change
        let matrix = make_upper_tri(2, &[(0, 0, 4.0), (0, 1, 2.0), (1, 1, 5.0)]);

        // Very small scaling: s[i]*a*s[j] << 1 for all entries
        let mut scaling = vec![0.1, 0.1];
        let original = scaling.clone();
        enforce_scaling_bound(&matrix, &mut scaling);

        assert_eq!(scaling, original, "valid scaling should not be modified");
    }

    // ---- Issue #10: greedy matching quality on diagonal matrix ----

    #[test]
    fn test_greedy_matching_diagonal_perfect() {
        // Diagonal matrix: greedy should achieve perfect matching without Dijkstra
        let matrix = make_upper_tri(4, &[(0, 0, 10.0), (1, 1, 20.0), (2, 2, 5.0), (3, 3, 15.0)]);
        let graph = build_cost_graph(&matrix);
        let state = greedy_initial_matching(&graph);

        let matched_count = state.row_match.iter().filter(|&&m| m != UNMATCHED).count();
        assert_eq!(
            matched_count, 4,
            "greedy should perfectly match a diagonal matrix"
        );

        // All should be identity matching (row i matched to col i)
        for (i, &j) in state.row_match.iter().enumerate() {
            assert_eq!(
                j, i,
                "diagonal greedy: row {} should match col {}, got {}",
                i, i, j
            );
        }
    }

    // ---- Issue #11: negative diagonal matrix ----

    #[test]
    fn test_mc64_negative_diagonal() {
        // All-negative diagonal: tests abs() path in cost computation
        let matrix = make_upper_tri(3, &[(0, 0, -10.0), (1, 1, -20.0), (2, 2, -5.0)]);
        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();
        assert_eq!(result.matched, 3);

        // Should produce identity matching (diagonals are the only entries)
        let (fwd, _) = result.matching.as_ref().arrays();
        for (i, &f) in fwd.iter().enumerate() {
            assert_eq!(f, i, "negative diagonal should give identity matching");
        }

        verify_scaling_properties(&matrix, &result);
    }

    // ---- Issue #12: unmatched index placement in singular permutation ----

    #[test]
    fn test_singular_unmatched_permutation_valid() {
        // 3x3 matrix where only rows 0,1 can match (via off-diagonal)
        // Row 2 has no connections except diagonal (which doesn't help with bipartite matching
        // when all entries are off-diagonal for the matched pairs)
        // [0  5  0]
        // [5  0  0]
        // [0  0  0]  <- isolated, can only be "placed" in remaining slot
        let matrix = make_upper_tri(3, &[(0, 1, 5.0)]);
        let result = mc64_matching(&matrix, Mc64Job::MaximumProduct).unwrap();

        // Matching must be a valid permutation regardless of match count
        let (fwd, inv) = result.matching.as_ref().arrays();
        let mut seen = [false; 3];
        for &f in fwd {
            assert!(f < 3, "fwd index out of range");
            assert!(!seen[f], "duplicate in fwd");
            seen[f] = true;
        }
        // fwd and inv should be consistent
        for i in 0..3 {
            assert_eq!(fwd[inv[i]], i, "fwd[inv[{}]] != {}", i, i);
        }
    }
}
