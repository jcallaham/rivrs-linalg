//! Elimination tree visualization and statistics.
//!
//! Provides text-based tree diagrams for small trees and summary statistics
//! for any size, supporting debugging of symbolic analysis phases.

use std::fmt::Write;

/// Summary statistics for an elimination tree.
#[derive(Debug, Clone)]
pub struct EliminationTreeStats {
    /// Total number of nodes.
    pub num_nodes: usize,
    /// Number of leaf nodes (no children).
    pub num_leaves: usize,
    /// Maximum root-to-leaf path length.
    pub depth: usize,
    /// Minimum children per non-leaf node.
    pub branching_min: usize,
    /// Maximum children per non-leaf node.
    pub branching_max: usize,
    /// Average children per non-leaf node.
    pub branching_mean: f64,
    /// Subtree size rooted at each node.
    pub subtree_sizes: Vec<usize>,
}

/// Text-based elimination tree display.
pub struct ETreeDisplay {
    parent: Vec<usize>,
    children: Vec<Vec<usize>>,
    roots: Vec<usize>,
}

impl ETreeDisplay {
    /// Construct from a standard parent-pointer array.
    ///
    /// `parent[i]` is the parent of node `i`. A root node has `parent[i] == i`
    /// or `parent[i] >= n` (sentinel).
    pub fn from_parent_array(parent: &[usize]) -> Self {
        let n = parent.len();
        let mut children = vec![vec![]; n];
        let mut roots = Vec::new();

        for (i, &p) in parent.iter().enumerate() {
            if p == i || p >= n {
                roots.push(i);
            } else {
                children[p].push(i);
            }
        }

        Self {
            parent: parent.to_vec(),
            children,
            roots,
        }
    }

    /// For small trees (n < 20), render a text tree diagram.
    /// For larger trees, falls back to stats view.
    pub fn render_tree(&self) -> String {
        let n = self.parent.len();
        if n == 0 {
            return "Elimination Tree (n=0)\n(empty)\n".to_string();
        }
        if n >= 20 {
            return self.render_stats();
        }

        let mut out = String::new();
        writeln!(out, "Elimination Tree (n={n})").unwrap();

        for (idx, &root) in self.roots.iter().enumerate() {
            let is_last = idx == self.roots.len() - 1;
            self.render_node(&mut out, root, "", is_last);
        }

        out
    }

    /// Render summary statistics.
    pub fn render_stats(&self) -> String {
        let stats = self.stats();
        let mut out = String::new();

        writeln!(out, "Elimination Tree Statistics (n={})", stats.num_nodes).unwrap();
        writeln!(out, "  Depth:            {}", stats.depth).unwrap();
        writeln!(
            out,
            "  Leaves:           {} ({:.1}%)",
            stats.num_leaves,
            if stats.num_nodes > 0 {
                stats.num_leaves as f64 / stats.num_nodes as f64 * 100.0
            } else {
                0.0
            }
        )
        .unwrap();

        if stats.num_nodes > stats.num_leaves {
            writeln!(
                out,
                "  Branching factor:  min={}, max={}, mean={:.1}",
                stats.branching_min, stats.branching_max, stats.branching_mean
            )
            .unwrap();
        }

        if !stats.subtree_sizes.is_empty() {
            let mut sorted = stats.subtree_sizes.clone();
            sorted.sort_unstable();
            let median = sorted[sorted.len() / 2];
            writeln!(
                out,
                "  Subtree sizes:    min={}, max={}, median={}",
                sorted[0],
                sorted[sorted.len() - 1],
                median
            )
            .unwrap();
        }

        out
    }

    /// Compute and return statistics for programmatic access.
    pub fn stats(&self) -> EliminationTreeStats {
        let n = self.parent.len();
        if n == 0 {
            return EliminationTreeStats {
                num_nodes: 0,
                num_leaves: 0,
                depth: 0,
                branching_min: 0,
                branching_max: 0,
                branching_mean: 0.0,
                subtree_sizes: Vec::new(),
            };
        }

        // Count leaves
        let num_leaves = self.children.iter().filter(|c| c.is_empty()).count();

        // Compute depth via DFS
        let depth = self
            .roots
            .iter()
            .map(|&r| self.node_depth(r))
            .max()
            .unwrap_or(0);

        // Branching factor stats (non-leaf nodes only)
        let non_leaf_children: Vec<usize> = self
            .children
            .iter()
            .filter(|c| !c.is_empty())
            .map(|c| c.len())
            .collect();

        let (branching_min, branching_max, branching_mean) = if non_leaf_children.is_empty() {
            (0, 0, 0.0)
        } else {
            let min = *non_leaf_children.iter().min().unwrap();
            let max = *non_leaf_children.iter().max().unwrap();
            let sum: usize = non_leaf_children.iter().sum();
            let mean = sum as f64 / non_leaf_children.len() as f64;
            (min, max, mean)
        };

        // Compute subtree sizes
        let mut subtree_sizes = vec![0usize; n];
        for &root in &self.roots {
            self.compute_subtree_sizes(root, &mut subtree_sizes);
        }

        EliminationTreeStats {
            num_nodes: n,
            num_leaves,
            depth,
            branching_min,
            branching_max,
            branching_mean,
            subtree_sizes,
        }
    }

    /// Compute depth of the subtree rooted at `node`.
    fn node_depth(&self, node: usize) -> usize {
        if self.children[node].is_empty() {
            1
        } else {
            1 + self.children[node]
                .iter()
                .map(|&c| self.node_depth(c))
                .max()
                .unwrap_or(0)
        }
    }

    /// Compute subtree sizes bottom-up.
    fn compute_subtree_sizes(&self, node: usize, sizes: &mut [usize]) {
        let mut size = 1;
        for &child in &self.children[node] {
            self.compute_subtree_sizes(child, sizes);
            size += sizes[child];
        }
        sizes[node] = size;
    }

    /// Render a single node with tree-drawing characters.
    fn render_node(&self, out: &mut String, node: usize, prefix: &str, is_last: bool) {
        let connector = if is_last {
            "\u{2514}\u{2500}\u{2500} "
        } else {
            "\u{251C}\u{2500}\u{2500} "
        };
        writeln!(out, "{prefix}{connector}{node}").unwrap();

        let child_prefix = if is_last {
            format!("{prefix}    ")
        } else {
            format!("{prefix}\u{2502}   ")
        };

        for (idx, &child) in self.children[node].iter().enumerate() {
            let child_is_last = idx == self.children[node].len() - 1;
            self.render_node(out, child, &child_prefix, child_is_last);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Linear chain: 0 → 1 → 2 → 3 → 4 (4 is root)
    fn linear_chain_5() -> Vec<usize> {
        vec![1, 2, 3, 4, 4]
    }

    /// Balanced binary tree:
    ///       6
    ///      / \
    ///     4   5
    ///    / \
    ///   0   1
    ///  2   3
    /// Actually: let's make a simpler one:
    ///     6
    ///    / \
    ///   2   5
    ///  / \ / \
    /// 0  1 3  4
    fn balanced_binary_7() -> Vec<usize> {
        // parent[0]=2, parent[1]=2, parent[2]=6, parent[3]=5, parent[4]=5, parent[5]=6, parent[6]=6 (root)
        vec![2, 2, 6, 5, 5, 6, 6]
    }

    /// Star topology: all nodes point to root (node 4)
    fn star_5() -> Vec<usize> {
        vec![4, 4, 4, 4, 4]
    }

    // --- EliminationTreeStats tests ---

    #[test]
    fn stats_linear_chain() {
        let parent = linear_chain_5();
        let etree = ETreeDisplay::from_parent_array(&parent);
        let stats = etree.stats();

        assert_eq!(stats.num_nodes, 5);
        assert_eq!(stats.num_leaves, 1); // only node 0 is a leaf
        assert_eq!(stats.depth, 5); // chain of 5
        assert_eq!(stats.branching_min, 1);
        assert_eq!(stats.branching_max, 1);
        assert!((stats.branching_mean - 1.0).abs() < 0.01);
    }

    #[test]
    fn stats_balanced_binary() {
        let parent = balanced_binary_7();
        let etree = ETreeDisplay::from_parent_array(&parent);
        let stats = etree.stats();

        assert_eq!(stats.num_nodes, 7);
        assert_eq!(stats.num_leaves, 4); // nodes 0, 1, 3, 4
        assert_eq!(stats.depth, 3);
        assert_eq!(stats.branching_min, 2);
        assert_eq!(stats.branching_max, 2);
        assert!((stats.branching_mean - 2.0).abs() < 0.01);
    }

    #[test]
    fn stats_star_topology() {
        let parent = star_5();
        let etree = ETreeDisplay::from_parent_array(&parent);
        let stats = etree.stats();

        assert_eq!(stats.num_nodes, 5);
        assert_eq!(stats.num_leaves, 4); // nodes 0-3
        assert_eq!(stats.depth, 2); // root + one level
        assert_eq!(stats.branching_min, 4);
        assert_eq!(stats.branching_max, 4);
    }

    #[test]
    fn stats_subtree_sizes_linear() {
        let parent = linear_chain_5();
        let etree = ETreeDisplay::from_parent_array(&parent);
        let stats = etree.stats();

        // Node 0: size 1, node 1: size 2, ..., node 4: size 5
        assert_eq!(stats.subtree_sizes[0], 1);
        assert_eq!(stats.subtree_sizes[1], 2);
        assert_eq!(stats.subtree_sizes[4], 5);
    }

    #[test]
    fn stats_subtree_sizes_star() {
        let parent = star_5();
        let etree = ETreeDisplay::from_parent_array(&parent);
        let stats = etree.stats();

        for i in 0..4 {
            assert_eq!(stats.subtree_sizes[i], 1);
        }
        assert_eq!(stats.subtree_sizes[4], 5);
    }

    #[test]
    fn stats_empty_tree() {
        let etree = ETreeDisplay::from_parent_array(&[]);
        let stats = etree.stats();
        assert_eq!(stats.num_nodes, 0);
        assert_eq!(stats.num_leaves, 0);
        assert_eq!(stats.depth, 0);
    }

    // --- Render tests ---

    #[test]
    fn render_tree_small_has_box_drawing() {
        let parent = balanced_binary_7();
        let etree = ETreeDisplay::from_parent_array(&parent);
        let rendered = etree.render_tree();

        assert!(rendered.contains("Elimination Tree (n=7)"));
        assert!(
            rendered.contains("\u{251C}\u{2500}\u{2500}")
                || rendered.contains("\u{2514}\u{2500}\u{2500}"),
            "should contain box-drawing characters ├── or └──"
        );
    }

    #[test]
    fn render_tree_shows_parent_child_relationships() {
        // Simple tree: 0 → 2, 1 → 2, 2 is root
        let parent = vec![2, 2, 2];
        let etree = ETreeDisplay::from_parent_array(&parent);
        let rendered = etree.render_tree();

        assert!(rendered.contains("2"), "should show root node 2");
        assert!(rendered.contains("0"), "should show child node 0");
        assert!(rendered.contains("1"), "should show child node 1");
    }

    #[test]
    fn render_tree_large_falls_back_to_stats() {
        // 20 nodes — should fall back to stats
        let parent: Vec<usize> = (0..20).map(|i| if i == 19 { 19 } else { i + 1 }).collect();
        let etree = ETreeDisplay::from_parent_array(&parent);
        let rendered = etree.render_tree();

        assert!(
            rendered.contains("Elimination Tree Statistics"),
            "n>=20 should fall back to stats view"
        );
    }

    #[test]
    fn render_stats_format() {
        let parent = balanced_binary_7();
        let etree = ETreeDisplay::from_parent_array(&parent);
        let stats_str = etree.render_stats();

        assert!(stats_str.contains("Elimination Tree Statistics (n=7)"));
        assert!(stats_str.contains("Depth:"));
        assert!(stats_str.contains("Leaves:"));
        assert!(stats_str.contains("Branching factor:"));
        assert!(stats_str.contains("Subtree sizes:"));
    }

    #[test]
    fn render_tree_empty() {
        let etree = ETreeDisplay::from_parent_array(&[]);
        let rendered = etree.render_tree();
        assert!(rendered.contains("n=0"));
        assert!(rendered.contains("empty"));
    }

    #[test]
    fn render_tree_single_root() {
        let parent = vec![0]; // single node, root
        let etree = ETreeDisplay::from_parent_array(&parent);
        let rendered = etree.render_tree();
        assert!(rendered.contains("Elimination Tree (n=1)"));
        assert!(rendered.contains("0"), "should show the single node");
    }

    #[test]
    fn multiple_roots() {
        // Two disconnected components: 0→2, 1→2 (root), 3→4 (root)
        let parent = vec![2, 2, 2, 4, 4];
        let etree = ETreeDisplay::from_parent_array(&parent);
        let stats = etree.stats();

        // Two roots: 2 and 4
        assert_eq!(etree.roots.len(), 2);
        assert_eq!(stats.depth, 2);
    }
}
