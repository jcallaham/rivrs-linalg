//! Debug visualization tools for sparse solver development.
//!
//! Provides text-based visualizations for sparsity patterns and elimination
//! trees, useful for inspecting matrix structure during solver development.
//!
//! Gated behind the `test-util` Cargo feature flag. Not included in production builds.
//!
//! # Examples
//!
//! ```
//! use rivrs_sparse::debug::ETreeDisplay;
//!
//! // Visualize a small elimination tree
//! let parent = vec![2, 2, 2]; // nodes 0,1 → root 2
//! let etree = ETreeDisplay::from_parent_array(&parent);
//! let rendered = etree.render_tree();
//! assert!(rendered.contains("Elimination Tree"));
//! ```

mod etree;
mod sparsity;

pub use etree::{ETreeDisplay, EliminationTreeStats};
pub use sparsity::SparsityDisplay;
