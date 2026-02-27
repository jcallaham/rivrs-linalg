# API Contract: Debug Visualization Module

**Module path**: `rivrs_sparse::debug`
**Feature gate**: `test-util`

## Sparsity Visualization

### SparsityDisplay

```
SparsityDisplay
  ::from_sparse(matrix: &SparseColMat<usize, f64>) -> SparsityDisplay
  .with_max_width(cols: usize) -> Self
  .with_max_height(rows: usize) -> Self
  .with_ascii_only(ascii: bool) -> Self
  .render() -> String
```

- `from_sparse(matrix)`: Analyzes the sparsity structure of a CSC matrix.
- `with_max_width(cols)`: Sets maximum display width (default: 80).
- `with_max_height(rows)`: Sets maximum display height (default: 40).
- `with_ascii_only(ascii)`: If true, uses ASCII characters only (`#.+-`). Default: false (Unicode blocks `█▓▒░`).
- `render()`: Produces the complete text output including header and grid.

### Implements `fmt::Display`

`SparsityDisplay` implements `Display`, delegating to `render()`.

### Output Format

```
bcsstk05 (153×153, nnz=2,423, density=10.3%)
─────────────────────────────────────────────
██▒░..........░▒
█▓▒░..........░▒
▒▓█▓▒░........░.
░▒▓█▓▒.........
..░▒█▓▒░.......
...░▒▓█▓▒.....
.....░▒█▓▒░...
......░▒▓█▓▒..
.......░▒█▓▒░.
........░▒▓█▓▒
░........░▒▓█▓
▒░........░▒▓█
```

When downsampled:
```
bcsstk05 (153×153, nnz=2,423, density=10.3%) [20×10 view]
```

### Density Character Mapping

**Unicode mode** (default):
| Density | Character |
|---------|-----------|
| 0%      | `.`       |
| 1-25%   | `░`       |
| 26-50%  | `▒`       |
| 51-75%  | `▓`       |
| 76-100% | `█`       |

**ASCII mode**:
| Density | Character |
|---------|-----------|
| 0%      | `.`       |
| 1-33%   | `-`       |
| 34-66%  | `+`       |
| 67-100% | `#`       |

## Elimination Tree Visualization

### ETreeDisplay

```
ETreeDisplay
  ::from_parent_array(parent: &[usize]) -> ETreeDisplay
  .render_tree() -> String
  .render_stats() -> String
  .stats() -> EliminationTreeStats
```

- `from_parent_array(parent)`: Constructs from a standard parent-pointer array where `parent[i]` is the parent of node `i` (root has `parent[i] == i` or a sentinel value).
- `render_tree()`: For small trees (n < 20), produces a text tree showing parent-child relationships. For larger trees, returns the stats view instead.
- `render_stats()`: Produces a summary statistics block.
- `stats()`: Returns the `EliminationTreeStats` struct for programmatic access.

### Tree Output Format (small trees)

```
Elimination Tree (n=10)
├── 0
│   ├── 1
│   └── 2
├── 3
│   ├── 4
│   │   └── 5
│   └── 6
└── 7
    ├── 8
    └── 9
```

### Stats Output Format

```
Elimination Tree Statistics (n=840)
  Depth:            23
  Leaves:           412 (49.0%)
  Branching factor:  min=1, max=8, mean=2.3
  Subtree sizes:    min=1, max=840, median=3
```

## Usage Patterns

### Sparsity Pattern

```
use rivrs_sparse::debug::SparsityDisplay;

let matrix = load_matrix("bcsstk05");
let display = SparsityDisplay::from_sparse(&matrix)
    .with_max_width(60);
println!("{}", display.render());
```

### Side-by-Side Comparison

```
let before = SparsityDisplay::from_sparse(&original)
    .with_max_width(40)
    .with_max_height(20);
let after = SparsityDisplay::from_sparse(&reordered)
    .with_max_width(40)
    .with_max_height(20);
println!("Before:                                  After:");
// ... interleave lines ...
```

### Elimination Tree

```
use rivrs_sparse::debug::ETreeDisplay;

let etree = ETreeDisplay::from_parent_array(&parent);
println!("{}", etree.render_stats());
```
