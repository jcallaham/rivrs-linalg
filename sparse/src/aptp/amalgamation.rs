//! Supernode amalgamation pass for the multifrontal factorization.
//!
//! Merges small parent-child supernode pairs after faer's symbolic analysis
//! to reduce assembly and extraction overhead. This implements SPRAL's
//! two-condition merge predicate:
//!
//! 1. **Structural match**: parent has 1 eliminated column and column count
//!    matches child's minus 1 (zero fill-in merge)
//! 2. **Both small**: both parent and child have fewer than `nemin` eliminated
//!    columns (bounded fill-in merge)
//!
//! # References
//!
//! - Liu (1992), "The Multifrontal Method for Sparse Matrix Solution: Theory
//!   and Practice" — supernode definitions, assembly trees
//! - SPRAL `core_analyse.f90:528-853` — `find_supernodes`, `do_merge`,
//!   `merge_nodes` (BSD-3-Clause)

use std::ops::Range;

use super::numeric::SupernodeInfo;

/// Check whether a child supernode should merge into its parent.
///
/// Implements SPRAL's `do_merge` predicate (`core_analyse.f90:806-822`):
///
/// - **Condition (a)**: structural match — parent has exactly 1 eliminated
///   column and its column count equals child's minus 1. This produces
///   zero fill-in.
/// - **Condition (b)**: both small — both parent and child have fewer than
///   `nemin` eliminated columns.
///
/// # Arguments
///
/// - `parent_nelim`: number of eliminated columns in parent
/// - `parent_cc`: total column count of parent (nelim + pattern.len())
/// - `child_nelim`: number of eliminated columns in child
/// - `child_cc`: total column count of child (nelim + pattern.len())
/// - `nemin`: minimum supernode size threshold
fn do_merge(
    parent_nelim: usize,
    parent_cc: usize,
    child_nelim: usize,
    child_cc: usize,
    nemin: usize,
) -> bool {
    // Condition (a): structural match (zero fill-in)
    if parent_nelim == 1 && parent_cc == child_cc.saturating_sub(1) {
        return true;
    }
    // Condition (b): both nodes are small
    if parent_nelim < nemin && child_nelim < nemin {
        return true;
    }
    false
}

/// Compute the sorted union of two sorted slices, excluding indices in a range.
///
/// Returns `sorted_union(a, b) \ {x : exclude_range.contains(x)}`.
/// Both input slices must be sorted in ascending order.
///
/// Used during supernode merging to compute the merged pattern:
/// `merged_pattern = sorted_union(parent.pattern, child.pattern)` minus
/// the fully-summed columns of the merged supernode.
fn sorted_union_excluding(a: &[usize], b: &[usize], exclude_range: Range<usize>) -> Vec<usize> {
    let mut result = Vec::with_capacity(a.len() + b.len());
    let mut ia = 0;
    let mut ib = 0;

    while ia < a.len() && ib < b.len() {
        let va = a[ia];
        let vb = b[ib];
        let next = if va < vb {
            ia += 1;
            va
        } else if vb < va {
            ib += 1;
            vb
        } else {
            // equal — take one, advance both
            ia += 1;
            ib += 1;
            va
        };
        if !exclude_range.contains(&next) {
            result.push(next);
        }
    }
    // Drain remaining
    for &v in &a[ia..] {
        if !exclude_range.contains(&v) {
            result.push(v);
        }
    }
    for &v in &b[ib..] {
        if !exclude_range.contains(&v) {
            result.push(v);
        }
    }
    result
}

/// Amalgamate supernodes by merging small parent-child pairs.
///
/// Processes the assembly tree in postorder, merging child supernodes into
/// their parents when SPRAL's merge predicate is satisfied. Returns a
/// compacted `Vec<SupernodeInfo>` with renumbered parent pointers.
///
/// # Arguments
///
/// - `supernodes`: fundamental supernodes from `build_supernode_info()`
/// - `nemin`: minimum supernode size threshold. `nemin = 1` disables amalgamation.
///
/// # SPRAL Equivalent
///
/// Corresponds to the amalgamation logic in `find_supernodes`
/// (`core_analyse.f90:618-641`) with `do_merge` predicate
/// (`core_analyse.f90:806-822`) and `merge_nodes` operation
/// (`core_analyse.f90:827-853`).
pub(crate) fn amalgamate(mut supernodes: Vec<SupernodeInfo>, nemin: usize) -> Vec<SupernodeInfo> {
    let n = supernodes.len();
    if n <= 1 {
        return supernodes;
    }

    // Track which supernodes are deleted (merged into parent)
    let mut deleted = vec![false; n];

    // Track accumulated nelim per supernode (starts as col_end - col_begin)
    let mut nelim: Vec<usize> = supernodes
        .iter()
        .map(|sn| sn.col_end - sn.col_begin)
        .collect();

    // Build children lists from parent pointers
    let mut children = vec![Vec::new(); n];
    for (s, sn) in supernodes.iter().enumerate() {
        if let Some(p) = sn.parent {
            children[p].push(s);
        }
    }

    // Process parents in ascending order (postorder: children have lower indices)
    for p in 0..n {
        if deleted[p] {
            continue;
        }
        // Iterate over children of p, checking merge predicate for each.
        // We need to handle the fact that merging changes nelim/cc of p,
        // affecting subsequent merge decisions for later children.
        // Take the children list to avoid borrow issues.
        let p_children = std::mem::take(&mut children[p]);
        for &c in &p_children {
            if deleted[c] {
                continue;
            }
            let parent_nelim = nelim[p];
            let parent_cc = nelim[p] + supernodes[p].pattern.len();
            let child_nelim = nelim[c];
            let child_cc = nelim[c] + supernodes[c].pattern.len();

            if do_merge(parent_nelim, parent_cc, child_nelim, child_cc, nemin) {
                // Merge child c into parent p
                let new_col_begin = supernodes[p].col_begin.min(supernodes[c].col_begin);
                let new_col_end = supernodes[p].col_end.max(supernodes[c].col_end);
                let exclude = new_col_begin..new_col_end;

                // Take child pattern out to avoid double borrow
                let child_pattern = std::mem::take(&mut supernodes[c].pattern);
                let parent_pattern = std::mem::take(&mut supernodes[p].pattern);
                let merged_pattern =
                    sorted_union_excluding(&parent_pattern, &child_pattern, exclude);

                // Merge owned_ranges: parent inherits child's owned column ranges
                let child_owned = std::mem::take(&mut supernodes[c].owned_ranges);
                supernodes[p].owned_ranges.extend(child_owned);

                supernodes[p].col_begin = new_col_begin;
                supernodes[p].col_end = new_col_end;
                supernodes[p].pattern = merged_pattern;
                nelim[p] += nelim[c];

                // Reparent c's children to p
                let c_children = std::mem::take(&mut children[c]);
                for &gc in &c_children {
                    if !deleted[gc] {
                        supernodes[gc].parent = Some(p);
                    }
                }
                children[p].extend(c_children);

                // Mark child as deleted
                deleted[c] = true;
            }
        }
        // Restore undeleted children + any newly reparented ones.
        // children[p] already has the reparented ones from merges above.
        // We also need to add back the original children that weren't merged.
        for &c in &p_children {
            if !deleted[c] {
                children[p].push(c);
            }
        }
    }

    // Compact: remove deleted supernodes, renumber parent pointers
    let mut old_to_new = vec![0usize; n];
    let mut new_idx = 0;
    for s in 0..n {
        if !deleted[s] {
            old_to_new[s] = new_idx;
            new_idx += 1;
        }
    }

    let mut result = Vec::with_capacity(new_idx);
    for s in 0..n {
        if !deleted[s] {
            let mut sn = std::mem::replace(
                &mut supernodes[s],
                SupernodeInfo {
                    col_begin: 0,
                    col_end: 0,
                    pattern: Vec::new(),
                    parent: None,
                    owned_ranges: Vec::new(),
                    in_small_leaf: false,
                },
            );
            sn.parent = sn.parent.map(|p| old_to_new[p]);
            result.push(sn);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to construct a SupernodeInfo with owned_ranges = [col_begin..col_end].
    #[allow(clippy::single_range_in_vec_init)]
    fn sn(
        col_begin: usize,
        col_end: usize,
        pattern: Vec<usize>,
        parent: Option<usize>,
    ) -> SupernodeInfo {
        SupernodeInfo {
            col_begin,
            col_end,
            pattern,
            parent,
            owned_ranges: vec![col_begin..col_end],
            in_small_leaf: false,
        }
    }

    // -----------------------------------------------------------------------
    // Phase 2: Foundational — do_merge and sorted_union_excluding tests
    // -----------------------------------------------------------------------

    /// T008: Structural match — parent with 1 eliminated col and cc(parent) == cc(child) - 1.
    #[test]
    fn test_do_merge_structural_match() {
        // Parent: 1 eliminated col, cc = 5.  Child: 2 eliminated cols, cc = 6.
        // Condition (a): nelim(parent)==1 AND cc(parent)==cc(child)-1 → true
        assert!(do_merge(1, 5, 2, 6, 32));
    }

    /// T009: Both nodes small (nemin condition).
    #[test]
    fn test_do_merge_nemin_both_small() {
        // Both have nelim < nemin=32
        assert!(do_merge(4, 20, 8, 30, 32));
    }

    /// T010: One node is large (>= nemin), should NOT merge.
    #[test]
    fn test_do_merge_one_large() {
        // Parent has nelim=32 >= nemin=32, child has nelim=4 < nemin
        assert!(!do_merge(32, 50, 4, 20, 32));
        // Child has nelim=32 >= nemin=32, parent has nelim=4 < nemin
        assert!(!do_merge(4, 20, 32, 50, 32));
    }

    /// T011: Both nodes large, should NOT merge.
    #[test]
    fn test_do_merge_both_large() {
        assert!(!do_merge(40, 60, 35, 50, 32));
    }

    /// T012: Sorted union of disjoint sets.
    #[test]
    fn test_sorted_union_disjoint() {
        let result = sorted_union_excluding(&[1, 3, 5], &[2, 4, 6], 0..0);
        assert_eq!(result, vec![1, 2, 3, 4, 5, 6]);
    }

    /// T013: Sorted union of overlapping sets.
    #[test]
    fn test_sorted_union_overlapping() {
        let result = sorted_union_excluding(&[1, 3, 5], &[3, 5, 7], 0..0);
        assert_eq!(result, vec![1, 3, 5, 7]);
    }

    /// T014: Sorted union with exclusion range.
    #[test]
    fn test_sorted_union_with_exclusion() {
        let result = sorted_union_excluding(&[5, 8, 10], &[3, 5, 7], 3..6);
        assert_eq!(result, vec![7, 8, 10]);
    }

    // -----------------------------------------------------------------------
    // Phase 3: US1 — Core amalgamation algorithm tests
    // -----------------------------------------------------------------------

    /// T019: All supernodes large (> nemin) — no merges, output identical to input.
    #[test]
    fn test_no_merges_large_supernodes() {
        let supernodes = vec![
            sn(0, 40, vec![200, 201], Some(4)),
            sn(40, 80, vec![200, 202], Some(4)),
            sn(80, 120, vec![200, 203], Some(4)),
            sn(120, 160, vec![200, 204], Some(4)),
            sn(160, 200, vec![], None),
        ];
        let result = amalgamate(supernodes, 32);
        assert_eq!(
            result.len(),
            5,
            "no merges expected — all supernodes are large"
        );
    }

    /// T020: Simple parent-child pair, both small — should merge.
    #[test]
    fn test_nemin_merge_simple_pair() {
        // child: cols [0,4), pattern [4,5,10]
        // parent: cols [4,8), pattern [10]
        // Both nelim=4 < nemin=32 → merge
        let supernodes = vec![sn(0, 4, vec![4, 5, 10], Some(1)), sn(4, 8, vec![10], None)];
        let result = amalgamate(supernodes, 32);
        assert_eq!(result.len(), 1, "pair should merge into one");
        assert_eq!(result[0].col_begin, 0);
        assert_eq!(result[0].col_end, 8);
        // Pattern: union of [4,5,10] and [10] minus [0..8) = [10]
        assert_eq!(result[0].pattern, vec![10]);
        assert!(result[0].parent.is_none());
    }

    /// T021: Structural match merge — parent with 1 col, cc matches child.
    #[test]
    fn test_structural_match_merge() {
        // child: cols [0,3), pattern [3,10,20] → cc=6
        // parent: cols [3,4), pattern [10,20,30,40] → cc=5
        // Condition (a): nelim(parent)==1, cc(parent)==5==cc(child)-1=5 → merge
        let supernodes = vec![
            sn(0, 3, vec![3, 10, 20], Some(1)),
            sn(3, 4, vec![10, 20, 30, 40], None),
        ];
        let result = amalgamate(supernodes, 32);
        assert_eq!(result.len(), 1, "structural match should merge");
        assert_eq!(result[0].col_begin, 0);
        assert_eq!(result[0].col_end, 4);
        // Pattern: union of [3,10,20] and [10,20,30,40] minus [0..4) = [10,20,30,40]
        assert_eq!(result[0].pattern, vec![10, 20, 30, 40]);
    }

    /// T022: Chain of 5 small supernodes — progressive merge in postorder.
    #[test]
    fn test_chain_merge() {
        // s0→s1→s2→s3→s4 (chain), all nelim=2
        let supernodes = vec![
            sn(0, 2, vec![2, 3, 4, 5, 6, 7, 8, 9], Some(1)),
            sn(2, 4, vec![4, 5, 6, 7, 8, 9], Some(2)),
            sn(4, 6, vec![6, 7, 8, 9], Some(3)),
            sn(6, 8, vec![8, 9], Some(4)),
            sn(8, 10, vec![], None),
        ];
        let result = amalgamate(supernodes, 32);
        // With nemin=32, all nodes have nelim<32, so every child merges into parent.
        // s0 merges into s1, s1 merges into s2, etc.
        // Eventually all collapse into one supernode.
        assert_eq!(
            result.len(),
            1,
            "chain of 5 small supernodes should merge to 1"
        );
        assert_eq!(result[0].col_begin, 0);
        assert_eq!(result[0].col_end, 10);
        assert!(result[0].pattern.is_empty());
    }

    /// T023: Parent with 4 small children — all merge into parent.
    #[test]
    fn test_bushy_tree_merge() {
        // 4 children (indices 0-3), all small, all point to parent (index 4)
        let supernodes = vec![
            sn(0, 2, vec![8, 9], Some(4)),
            sn(2, 4, vec![8, 9], Some(4)),
            sn(4, 6, vec![8, 9], Some(4)),
            sn(6, 8, vec![8, 9], Some(4)),
            sn(8, 10, vec![], None),
        ];
        let result = amalgamate(supernodes, 32);
        // All children are small (<32), parent is small (<32).
        // After first merge: parent nelim becomes 2+2=4.
        // After second: 4+2=6. Third: 6+2=8. Fourth: 8+2=10.
        // All still < 32, so all merge.
        assert_eq!(
            result.len(),
            1,
            "all 4 small children should merge into parent"
        );
        assert_eq!(result[0].col_begin, 0);
        assert_eq!(result[0].col_end, 10);
    }

    /// T024: Parent with 3 children — 2 small, 1 large — only small merge.
    #[test]
    fn test_partial_merge_mixed_sizes() {
        let supernodes = vec![
            // small child 1
            sn(0, 4, vec![100, 140], Some(3)),
            // small child 2
            sn(4, 8, vec![100, 140], Some(3)),
            // large child (40 cols >= nemin=32)
            sn(8, 48, vec![100, 140], Some(3)),
            // parent (small, nelim=4)
            sn(100, 104, vec![140], Some(4)),
            // root (nelim=50, large)
            sn(140, 190, vec![], None),
        ];
        let result = amalgamate(supernodes, 32);
        // Small children (0,1) merge into parent (3) → parent nelim becomes 4+4+4=12.
        // Large child (2, nelim=40) cannot merge with parent (nelim=12<32 but child=40>=32).
        // Parent (nelim=12) cannot merge with root (nelim=50>=32).
        // Result: 3 supernodes (large child, merged parent, root)
        assert_eq!(
            result.len(),
            3,
            "only 2 small children should merge, large child stays"
        );
        // The large child should be unchanged
        assert_eq!(result[0].col_end - result[0].col_begin, 40);
    }

    /// T025: When child C merges into parent P, C's grandchildren become P's children.
    #[test]
    fn test_parent_reparenting() {
        // gc0, gc1 → child → parent (root)
        // child is small, parent is small → merge
        // After merge: gc0, gc1 → parent
        let supernodes = vec![
            // gc0
            sn(0, 2, vec![4, 5, 8, 9], Some(2)),
            // gc1
            sn(2, 4, vec![4, 5, 8, 9], Some(2)),
            // child (merges into parent)
            sn(4, 6, vec![8, 9], Some(3)),
            // parent (root)
            sn(8, 10, vec![], None),
        ];
        let result = amalgamate(supernodes, 32);
        // child merges into parent. gc0, gc1 are reparented.
        // Then gc0, gc1 are also small and can merge into parent.
        // All end up as 1 supernode.
        assert_eq!(result.len(), 1, "all small nodes should eventually merge");
        assert_eq!(result[0].col_begin, 0);
        assert_eq!(result[0].col_end, 10);
    }

    /// T026: Pattern union on merge — verify merged pattern is correct.
    #[test]
    fn test_pattern_union_on_merge() {
        // child: cols [0,2), pattern [2,5,10,20]
        // parent: cols [2,4), pattern [5,15,20]
        // Merged: cols [0,4), pattern = union([2,5,10,20],[5,15,20]) \ [0..4) = [5,10,15,20]
        let supernodes = vec![
            sn(0, 2, vec![2, 5, 10, 20], Some(1)),
            sn(2, 4, vec![5, 15, 20], None),
        ];
        let result = amalgamate(supernodes, 32);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].pattern, vec![5, 10, 15, 20]);
    }

    /// T027: After amalgamation, all parent indices > child indices (postorder).
    #[test]
    fn test_postorder_preserved() {
        let supernodes = vec![
            sn(0, 2, vec![4, 8], Some(2)),
            sn(2, 4, vec![4, 8], Some(2)),
            sn(4, 6, vec![8], Some(3)),
            sn(8, 12, vec![], None),
        ];
        // Use nemin=4 so that supernodes with nelim=2 merge but
        // supernodes with nelim=4 don't merge with other large ones
        let result = amalgamate(supernodes, 4);
        // Verify postorder: for all non-root supernodes, parent index > self index
        for (i, sn) in result.iter().enumerate() {
            if let Some(p) = sn.parent {
                assert!(
                    p > i,
                    "postorder violation: supernode {} has parent {} (should be > {})",
                    i,
                    p,
                    i
                );
            }
        }
    }

    /// T028: Single supernode (root) — returned unchanged.
    #[test]
    fn test_single_supernode_passthrough() {
        let supernodes = vec![sn(0, 10, vec![], None)];
        let result = amalgamate(supernodes, 32);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].col_begin, 0);
        assert_eq!(result[0].col_end, 10);
    }

    /// T028b: 100 single-column supernodes (simulating simplicial like bloweybq).
    #[test]
    fn test_simplicial_many_single_column_supernodes() {
        let n = 100;
        // Build a chain: s0→s1→s2→...→s99 (each 1 column)
        // Pattern for si = [i+1, i+2, ..., min(i+5, n-1)] (up to 5 pattern entries)
        let supernodes: Vec<SupernodeInfo> = (0..n)
            .map(|i| {
                let pattern: Vec<usize> = ((i + 1)..n.min(i + 6)).collect();
                let parent = if i + 1 < n { Some(i + 1) } else { None };
                sn(i, i + 1, pattern, parent)
            })
            .collect();
        let result = amalgamate(supernodes, 32);
        // With nemin=32, single-column supernodes should merge aggressively.
        // The exact count depends on pattern compatibility, but should be
        // significantly less than 100.
        assert!(
            result.len() < 20,
            "expected significant reduction from 100, got {}",
            result.len()
        );
        // Verify postorder
        for (i, sn) in result.iter().enumerate() {
            if let Some(p) = sn.parent {
                assert!(p > i, "postorder violation at {}: parent={}", i, p);
            }
        }
        // Verify col_begin < col_end for all
        for (i, sn) in result.iter().enumerate() {
            assert!(
                sn.col_begin < sn.col_end,
                "empty supernode at {}: [{}, {})",
                i,
                sn.col_begin,
                sn.col_end
            );
        }
        // Verify no pattern entry is within the supernode's own range
        for sn in &result {
            for &r in &sn.pattern {
                assert!(
                    r < sn.col_begin || r >= sn.col_end,
                    "pattern entry {} is within [{}, {})",
                    r,
                    sn.col_begin,
                    sn.col_end
                );
            }
        }
    }

    /// T028c: Root with 20 small children — accumulated nelim may prevent later merges.
    #[test]
    fn test_star_tree_many_children() {
        // 20 children (each nelim=2), parent (nelim=4), nemin=32
        // After merging some children, parent's accumulated nelim grows.
        // Once accumulated nelim >= 32, further nemin merges stop.
        let mut supernodes: Vec<SupernodeInfo> = (0..20)
            .map(|i| sn(i * 2, i * 2 + 2, vec![40, 41, 42, 43], Some(20)))
            .collect();
        supernodes.push(sn(40, 44, vec![], None));

        let result = amalgamate(supernodes, 32);

        // The parent starts with nelim=4. Merging children adds 2 each.
        // After 14 merges: 4 + 14*2 = 32 >= nemin → stop.
        // So ~14 children merge (up to 15 since check is <nemin), remaining ~5-6 stay.
        // Exact count depends on merge order.
        assert!(
            result.len() < 21,
            "at least some children should merge, got {} supernodes",
            result.len()
        );
        // Parent's accumulated nelim should be bounded by nemin considerations
        // Verify postorder
        for (i, sn) in result.iter().enumerate() {
            if let Some(p) = sn.parent {
                assert!(p > i, "postorder violation at {}: parent={}", i, p);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Phase 5 / US3: Configurable nemin threshold tests
    // -----------------------------------------------------------------------

    /// T042: nemin=1 disables amalgamation — same input as test_nemin_merge_simple_pair
    /// should produce no merges, output identical to input.
    #[test]
    fn test_nemin_1_disables_amalgamation() {
        // Parent-child pair, both with nelim=4
        let supernodes = vec![sn(0, 4, vec![4, 5, 6], Some(1)), sn(4, 8, vec![5, 6], None)];
        let result = amalgamate(supernodes.clone(), 1);

        assert_eq!(
            result.len(),
            2,
            "nemin=1 should disable amalgamation: got {} supernodes",
            result.len()
        );
        // Verify structure preserved
        assert_eq!(result[0].col_begin, 0);
        assert_eq!(result[0].col_end, 4);
        assert_eq!(result[1].col_begin, 4);
        assert_eq!(result[1].col_end, 8);
    }

    /// T043: nemin=64 merges more aggressively — supernodes with nelim in 32-63
    /// should merge with nemin=64 but not with nemin=32.
    #[test]
    fn test_nemin_64_more_aggressive() {
        // Parent-child pair, both with nelim=40 (> 32, < 64)
        let supernodes = vec![
            sn(0, 40, vec![40, 50, 60], Some(1)),
            sn(40, 80, vec![50, 60], None),
        ];

        // nemin=32 → no merge (both >= 32)
        let result_32 = amalgamate(supernodes.clone(), 32);
        assert_eq!(
            result_32.len(),
            2,
            "nemin=32 should NOT merge supernodes with nelim=40: got {} supernodes",
            result_32.len()
        );

        // nemin=64 → merge (both < 64)
        let result_64 = amalgamate(supernodes, 64);
        assert_eq!(
            result_64.len(),
            1,
            "nemin=64 should merge supernodes with nelim=40: got {} supernodes",
            result_64.len()
        );
        assert_eq!(result_64[0].col_begin, 0);
        assert_eq!(result_64[0].col_end, 80);
    }
}
