# API Contract: Internal Types (pub(crate))

**Module**: `src/aptp/numeric.rs` (internal types)
**Visibility**: `pub(crate)` unless noted

## SupernodeInfo

Precomputed unified supernode descriptor. Built once from `AptpSymbolic` before the factorization loop.

```
pub(crate) struct SupernodeInfo {
    col_begin: usize,
    col_end: usize,
    pattern: Vec<usize>,
    parent: Option<usize>,
}
```

### Construction

```
pub(crate) fn build_supernode_info(
    symbolic: &AptpSymbolic,
) -> Vec<SupernodeInfo>
```

- For supernodal: delegates to `AptpSymbolic::supernode_begin/end/pattern/parent`
- For simplicial: derives from `etree()` and simplicial CSC structure
- Returns one entry per supernode, indexed by supernode ID
- Postcondition: `info[s].parent.map_or(true, |p| p > s)` (postorder)

## FrontalMatrix

Temporary dense matrix for one supernode's assembly and factorization.

```
pub(crate) struct FrontalMatrix {
    data: Mat<f64>,
    row_indices: Vec<usize>,
    num_fully_summed: usize,
}
```

### Construction

```
pub(crate) fn assemble_frontal_matrix(
    sn: &SupernodeInfo,
    delayed_from_children: &[(Vec<usize>, /* delayed global indices */
                              Mat<f64>    /* delayed column data */)],
    matrix: &SparseColMat<usize, f64>,
    perm_fwd: &[usize],
    perm_inv: &[usize],
    child_contributions: &[ContributionBlock],
    children: &[usize],
    global_to_local: &mut [usize],
) -> FrontalMatrix
```

### Steps

1. Compute frontal row index set: supernode columns + delayed columns + pattern
2. Set up global-to-local mapping
3. Allocate Mat::zeros(m, m)
4. Scatter original matrix entries using permutation
5. Extend-add each child's contribution block
6. Return assembled FrontalMatrix

## ContributionBlock

```
pub(crate) struct ContributionBlock {
    data: Mat<f64>,
    row_indices: Vec<usize>,
    num_delayed: usize,
}
```

### Extraction

```
pub(crate) fn extract_contribution(
    frontal: &FrontalMatrix,
    result: &AptpFactorResult,
) -> ContributionBlock
```

- Copies trailing `(m - ne) x (m - ne)` submatrix from factored frontal matrix
- Computes row indices: delayed from `result.delayed_cols`, non-fully-summed from `frontal.row_indices[k..m]`
- Sets `num_delayed = k - ne`

## Extend-Add Operation

```
pub(crate) fn extend_add(
    parent_frontal: &mut FrontalMatrix,
    child_contribution: &ContributionBlock,
    global_to_local: &[usize],
)
```

- For each entry in child_contribution at local position (i, j):
  - Map to global indices via `child_contribution.row_indices`
  - Map to parent local indices via `global_to_local`
  - Add value to parent frontal matrix at mapped position
- Only lower triangle is processed (symmetric)

## Factor Extraction

```
pub(crate) fn extract_front_factors(
    frontal: &FrontalMatrix,
    result: &AptpFactorResult,
    sn: &SupernodeInfo,
) -> FrontFactors
```

- Extracts L11 from `frontal.data[0..ne, 0..ne]`
- Takes D11 from `result.d` (first ne entries)
- Extracts L21 from `frontal.data[k..m, 0..ne]`
- Records `local_perm = result.perm[0..k]`
- Computes `col_indices` from local_perm mapping through `frontal.row_indices`
- Sets `row_indices = frontal.row_indices[k..m]`
