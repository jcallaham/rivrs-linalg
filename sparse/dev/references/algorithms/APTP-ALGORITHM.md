# APTP Algorithm Reference

A unified pseudocode reference for the **A Posteriori Threshold Pivoting** algorithm,
sufficient for implementing a simplicial LDL^T factorization without consulting multiple papers.

---

## Table of Contents

1. [APTP Overview](#1-aptp-overview)
2. [Main APTP Factorization](#2-main-aptp-factorization)
3. [Dense Block Kernels](#3-dense-block-kernels)
4. [Complete Pivoting Algorithm](#4-complete-pivoting-algorithm)
5. [Pivot Acceptance and Column Delay](#5-pivot-acceptance-and-column-delay)
6. [Symbolic Analysis](#6-symbolic-analysis)
7. [Multifrontal Structure](#7-multifrontal-structure)
8. [Source Citations](#8-source-citations)

---

## 1. APTP Overview

**Source**: Duff, Hogg, Lopez (2020) Sections 1--3; Hogg, Ovtchinnikov, Scott (2016) Section 2.3;
Schenk & Gartner (2006) Section 2.

### What is APTP?

A Posteriori Threshold Pivoting (APTP) is a pivoting strategy for the LDL^T factorization of
sparse symmetric indefinite matrices. It was introduced by Duff, Hogg, and Lopez (2020) as an
evolution of the GPU-oriented a posteriori pivoting from the original SSIDS solver (Hogg,
Ovtchinnikov, Scott 2016). The key insight is that instead of deciding whether a pivot is
acceptable *before* performing elimination (as in traditional approaches), APTP:

1. **Performs elimination optimistically**, assuming all pivots in a block will succeed.
2. **Checks stability after the fact**, by verifying that no entry of L exceeds the
   threshold bound u^{-1}.
3. **Falls back gracefully** when a pivot fails: the block is restored from a backup
   ("fail-in-place"), and the failed columns are handled later.

This strategy enables **2D block partitioning** of the dense frontal matrix at each node of
the assembly tree, yielding significantly more parallelism than the 1D block-column partitioning
required by traditional threshold partial pivoting (TPP).

### Comparison: TPP vs SBK vs APTP

| Property | TPP (MA86/MA97) | SBK (PARDISO) | APTP (SSIDS) |
|---|---|---|---|
| **Pivot decision** | A priori: check column max before elimination | A priori: Bunch-Kaufman within supernode diagonal block | A posteriori: eliminate, then verify |
| **Stability guarantee** | Bounded L entries: \|l_{ij}\| < u^{-1} | No guarantee; uses perturbation + iterative refinement | Bounded L entries: \|l_{ij}\| < u^{-1} |
| **Block partitioning** | 1D block-column only | 2D (static supernodal structure) | 2D block partitioning |
| **Parallelism** | Limited by column dependencies | High (static structure, no delayed pivots) | High (speculative execution, fail-in-place) |
| **Failed pivots** | Delayed to parent node | Perturbed: set to sign(a_{ii}) * epsilon * \|\|A\|\|_1 | Fail-in-place, then fallback (re-factor or delay) |
| **Accuracy** | Robust; backward error near machine precision | May require iterative refinement; can fail on hard problems | Same robustness as TPP |
| **Inertia preservation** | Yes | No (perturbation changes inertia) | Yes |

**Source**: Duff et al. (2020), Section 2 describes TPP and SBK as background. Section 3 introduces
APTP. Schenk & Gartner (2006), Algorithm 2.1 defines the SBK pivot selection.

### Why APTP over TPP?

In TPP, the pivot decision for column j requires the maximum absolute value in column j
*below the diagonal*. When using 2D block partitioning, all subdiagonal blocks in that column
must be up-to-date before the pivot can be selected. This serializes the factorization within
each block column. APTP removes this bottleneck: the diagonal block is factored independently,
and the pivot check is distributed across all subdiagonal blocks, each of which can operate in
parallel.

### Why APTP over SBK?

SBK restricts pivoting to within the diagonal block of a supernode and perturbs the system
when no acceptable pivot is found. This means the factorization is *not* backward stable:
the computed factors satisfy P(A + E)P^T = LDL^T where E may be non-negligible. Iterative
refinement is needed but may not converge for hard indefinite systems. APTP maintains the
same backward stability as TPP while achieving comparable parallel performance to SBK.

**Source**: Duff et al. (2020), Section 2 (paragraphs on PARDISO); Tables 3--4 and Section 7
for experimental comparison showing PARDISO failures on stokes128 and cvxqp3.

---

## 2. Main APTP Factorization

**Source**: Duff, Hogg, Lopez (2020), Algorithm 3.1, Section 3.

### Parameters

- **nb**: Outer block size for 2D partitioning (typical value: 256).
- **ib**: Inner block size used within the Factor kernel (typical value: 32).
- **u**: Threshold parameter in (0, 0.5]. Default: 0.01. The stability bound is |l_{ij}| <= u^{-1}.
- **nelim_j**: Number of successfully eliminated pivots in block column j. Initialized to the
  block size and reduced atomically by Apply and Adjust kernels.
- **nblk**: Number of block columns = ceil(n / nb).
- **mblk**: Number of block rows = ceil(m / nb).

The frontal matrix F at a node has m rows and n fully summed columns (n <= m). It is
partitioned into nblk block columns and mblk block rows of size nb (the last block may be
smaller).

### Algorithm 3.1: APTP LDL^T Factorization

```
Algorithm 3.1: A posteriori threshold pivoting LDL^T factorization.
  [Duff, Hogg, Lopez (2020), Section 3, Algorithm 3.1]

  Input: m x n frontal matrix A (2D block partitioned, nb x nb blocks)
  Output: L, D factors; indices of delayed (uneliminated) columns

  Let nblk = ceil(n / nb)    -- number of block columns
  Let mblk = ceil(m / nb)    -- number of block rows

  for j = 1 to nblk do

      -- Step 1: Factorize the diagonal block
      Factor(A[j,j], nelim_j)

      -- Step 2: Apply factorization to superdiagonal blocks (blocks in row j, left of diagonal)
      for i = 1 to j-1 do
          ApplyT(A[j,i], nelim_j; A[j,j])
      end for

      -- Step 3: Apply factorization to subdiagonal blocks (blocks below diagonal in column j)
      for i = j+1 to mblk do
          ApplyN(A[i,j], nelim_j; A[j,j])
      end for

      -- Step 4: Synchronize and finalize nelim_j
      Adjust(nelim_j; all blocks in column j and row j)

      -- Step 5: Update previously failed columns (from iterations 1..j-1)
      for k = 1 to j-1 do
          for i = k to j-1 do
              UpdateTN(A[i,k]; A[j,i], A[j,k], nelim_j)
          end for
          for i = j to mblk do
              UpdateNT(A[i,k]; A[i,j], A[j,k], nelim_j)
          end for
      end for

      -- Step 6: Update trailing submatrix (columns j..nblk)
      for k = j to nblk do
          for i = k to mblk do
              UpdateNN(A[i,k]; A[i,j], A[k,j], nelim_j)
          end for
      end for

  end for

  -- Post-processing: permute failed columns to back, optionally re-factorize
```

### Task Dependencies

Each subroutine call represents a task. Parameters before the semicolon are **inout** (updated);
parameters after the semicolon are **in** (read-only). The variable nelim_j is updated via
**atomic reductions** (atomicMin), not through explicit task dependencies. This allows Apply
tasks to run in parallel while independently reducing nelim_j.

**Source**: Duff et al. (2020), Algorithm 3.1, paragraphs immediately following it.

---

## 3. Dense Block Kernels

**Source**: Duff, Hogg, Lopez (2020), Section 3 (kernel descriptions following Algorithm 3.1),
Section 4 (complete pivoting within Factor).

### Factor(A[j,j], nelim_j)

Factorizes the nb x nb diagonal block. Uses a **two-level hierarchical** approach:

```
Factor(A_jj, nelim_j):
  [Duff et al. (2020), Section 3 -- Factor kernel; Section 4 for inner pivoting]

  1. Store backup of block A_jj.
  2. Perform pivoted LDL^T factorization:
       P_j * A_jj * P_j^T = L_jj * D_jj * L_jj^T
     using APTP with inner block size ib:
       - Partition A_jj into inner blocks of size ib.
       - For each inner diagonal block, use Complete Pivoting (Algorithm 4.1).
       - Apply inner block's factors to remaining inner blocks (serial APTP loop).
  3. Initialize nelim_j = block_size(A_jj).
     (Will be reduced by Apply kernels if any entry of L exceeds u^{-1}.)
```

The two-level structure is illustrated in the blocking strategy figure (Duff et al. 2020, Figure 2):
the outer level uses block size nb with parallel APTP, and the inner level uses block size ib
with serial APTP. The innermost ib x ib diagonal blocks are factored using complete pivoting
(Algorithm 4.1).

### ApplyN(A[i,j], nelim_j; A[j,j])

Applies the diagonal block's factorization to a **subdiagonal** block:

```
ApplyN(A_ij, nelim_j; A_jj):
  [Duff et al. (2020), Section 3 -- ApplyN kernel]

  1. Store backup of block A_ij.
  2. Compute: L_ij = A_ij * P_j * (L_jj * D_jj)^{-T}
     That is, apply the permutation P_j to columns, then solve
     the triangular system to obtain the L entries.
  3. Scan columns of L_ij left to right. Find the first column index
     nelim_ij such that some entry |l_{pq}| > u^{-1}.
  4. Atomic reduction: nelim_j = min(nelim_j, nelim_ij).
```

### ApplyT(A[j,i], nelim_j; A[j,j])

Applies the diagonal block's factorization to a **superdiagonal** block (i.e., a block in
row j to the left of the diagonal, stored as transpose):

```
ApplyT(A_ji, nelim_j; A_jj):
  [Duff et al. (2020), Section 3 -- ApplyT kernel]

  1. Store backup of block A_ji.
  2. Compute: L_ji = (L_jj * D_jj)^{-1} * P_j^T * A_ji
     Only on columns corresponding to previously failed pivots.
  3. Scan rows of L_ji top to bottom. Find the first row index
     nelim_ji such that some entry |l_{pq}| > u^{-1}.
  4. Atomic reduction: nelim_j = min(nelim_j, nelim_ji).

  Note: If block A_ji contains no failed entries, this is a no-op.
```

### Adjust(nelim_j; all blocks A[:,j] and A[j,:])

Synchronization barrier that finalizes the value of nelim_j:

```
Adjust(nelim_j; all blocks in column j and row j):
  [Duff et al. (2020), Section 3 -- Adjust kernel]

  1. Wait until all ApplyN and ApplyT tasks for block column j have completed
     their atomic reductions.
  2. If nelim_j would accept only the first column of a 2x2 pivot
     (i.e., the pivot at position nelim_j is the first column of a
     2x2 Bunch-Kaufman pivot), decrement nelim_j by one.
     This ensures we never split a 2x2 pivot across the accepted/failed boundary.
```

### UpdateNN(A[i,k]; A[i,j], A[k,j], nelim_j)

Updates a block in the **trailing submatrix** (to the right and below the current block column):

```
UpdateNN(A_ik; A_ij, A_kj, nelim_j):
  [Duff et al. (2020), Section 3 -- UpdateNN kernel]

  1. If k == j (block is in the same column as the pivotal column):
       Restore failed entries from backup.
  2. Perform the rank-nelim_j update:
       A_ik = A_ik - L_ij * D_jj * L_kj^T
     where only the first nelim_j columns of L_ij and L_kj participate.
     If k == j, update only the uneliminated (restored) entries.
     Otherwise, update the entire block.
```

### UpdateNT(A[i,k]; A[i,j], A[j,k], nelim_j)

Updates a block at the **bottom-left** of the eliminated block column (row i >= j, column k < j):

```
UpdateNT(A_ik; A_ij, A_jk, nelim_j):
  [Duff et al. (2020), Section 3 -- UpdateNT kernel]

  1. If i == j (block is in the same row as the eliminated column):
       Restore failed entries from backup.
  2. Perform the update on uneliminated entries:
       A_ik = A_ik - L_ij * D_jj * L_jk
     where L_ij uses the first nelim_j columns and L_jk is from the
     superdiagonal application.
```

### UpdateTN(A[i,k]; A[j,i], A[j,k], nelim_j)

Updates a block at the **top-left** of the eliminated block column (row i < j, column k < j):

```
UpdateTN(A_ik; A_ji, A_jk, nelim_j):
  [Duff et al. (2020), Section 3 -- UpdateTN kernel]

  Perform the update on uneliminated entries:
    A_ik = A_ik - L_ji^T * D_jj * L_jk
  where L_ji is from ApplyT and L_jk is from ApplyT on another block.
```

---

## 4. Complete Pivoting Algorithm

**Source**: Duff, Hogg, Lopez (2020), Section 4, Algorithm 4.1 and stability analysis (Section 4.1).

At the finest granularity, ib x ib blocks on the diagonal are factored using a **complete pivoting**
strategy. This is chosen because:
- The block fits in cache, so searching the entire block for the best pivot is cheap.
- It causes failed pivots to occur as late as possible within each block (maximizing nelim).

### Algorithm 4.1: Complete Pivoting on a Dense Block

```
Algorithm 4.1: Complete pivoting algorithm.
  [Duff et al. (2020), Section 4, Algorithm 4.1]

  Input: ib x ib dense symmetric block A (partially factored)
  Output: L, D factors; pivot sequence; count of successful pivots
  Parameters: delta (singularity threshold, e.g. 1e-20 * max|A|)

  while uneliminated columns remain do

      1. Find the uneliminated entry with maximum absolute value.
         Let its position be (t, m) with value a_{tm}.

      2. if |a_{tm}| < delta then
             -- All remaining pivots are effectively zero.
             Record all remaining columns as zero pivots.
             break

      3. else if m == t then
             -- Maximum entry is on the diagonal.
             Use a_{mm} as a 1x1 pivot.
             Eliminate column m.

      4. else
             -- Maximum entry is off-diagonal at (t, m).
             Compute: Delta = a_{mm} * a_{tt} - a_{mt}^2

             if |Delta| >= (1/2) * |a_{mt}|^2 then
                 -- 2x2 pivot is stable.
                 Use (m, t) as a 2x2 pivot.
                 Eliminate columns m and t.

             else
                 -- 2x2 pivot would be unstable.
                 -- Fall back to 1x1 on the larger diagonal.
                 Use max(|a_{mm}|, |a_{tt}|) as a 1x1 pivot.
                 Eliminate the corresponding column.

             end if

      end if

  end while
```

### Stability Analysis

**Source**: Duff et al. (2020), Section 4.1.

Algorithm 4.1 is at least as stable as traditional threshold partial pivoting with u = 0.25.
The three cases are:

**Case 1: Immediate 1x1 pivot** (m == t, diagonal entry is the global maximum).

    |l_{i1}| = |a_{i1}| / |a_{11}| <= 1 < 4 = (0.25)^{-1}

**Case 2: Successful 2x2 pivot** (|Delta| >= (1/2)|a_{mt}|^2).

Without loss of generality, assume |a_{21}| >= |a_{11}| >= |a_{22}| and the condition
|Delta| >= (1/2)|a_{21}|^2 holds. Then:

    |l_{i1}|, |l_{i2}| <= 2 / (|a_{21}|^2) * (|a_{21}|^2 + |a_{21}|^2) = 4 = (0.25)^{-1}

This uses the triangle inequality and the fact that |a_{21}| >= |a_{i1}|, |a_{i2}| for all i.

**Case 3: Failed 2x2 pivot** (|Delta| < (1/2)|a_{mt}|^2, fall back to 1x1).

The failure of the 2x2 condition implies:

    (1/2)|a_{21}|^2 < |a_{11}| * |a_{22}| < |a_{11}|^2

Therefore 1/|a_{11}| < sqrt(2)/|a_{21}|, giving:

    |l_{i1}| = |a_{i1}|/|a_{11}| <= |a_{21}|/|a_{11}| <= sqrt(2) < 4

All three paths yield |l_{ij}| < 4 = (0.25)^{-1}, equivalent to u = 0.25 threshold partial pivoting.

---

## 5. Pivot Acceptance and Column Delay

**Source**: Duff, Hogg, Lopez (2020), Section 3 (fail-in-place mechanism),
Section 7 (threshold parameter guide, Table 5).

### Threshold Test

After the Factor kernel processes diagonal block A[j,j] and each Apply kernel computes L
entries, the threshold test is:

```
For each entry l_{ij} in the computed L block:
    if |l_{ij}| > u^{-1} then
        FAIL: this column and all subsequent columns in the block are marked as failed.
        nelim_j = min(nelim_j, index_of_first_failed_column)
```

The parameter u controls the trade-off between stability and performance:
- Larger u (e.g., 0.1): stricter test, more failed columns, more delayed pivots, better
  stability, slower factorization.
- Smaller u (e.g., 0.001): more permissive, fewer failures, faster factorization,
  slightly reduced stability (but still bounded).

### Fail-in-Place Mechanism

**Source**: Duff et al. (2020), Section 3 (paragraphs on fail-in-place and speculative execution).

The fail-in-place approach is central to APTP's efficiency:

```
Fail-in-Place Protocol:

  BEFORE elimination (in Factor, ApplyN, ApplyT):
    1. Store a backup copy of the block's entries.

  AFTER elimination:
    2. Check threshold: |l_{ij}| <= u^{-1} for all entries.
    3. If check PASSES: discard backup. Entries remain as computed.
    4. If check FAILS at column c:
         a. Record nelim = c - 1 (columns 0..c-1 are accepted).
         b. Atomic reduction: nelim_j = min(nelim_j, nelim).
         c. The failed columns REMAIN IN PLACE (not moved).

  DURING UpdateNN for k == j:
    5. Restore failed entries from backup before applying updates.
       (The failed columns now contain their original values, updated
       by all successfully eliminated columns from this block.)
```

This avoids the expensive row/column permutations that would be needed if failed columns
were moved out of the block immediately. Failed columns are simply left in their original
positions and updated incrementally as later blocks are processed.

### Column Delay to Parent Node

After the APTP loop completes for the entire frontal matrix:

```
Post-APTP Processing:
  [Duff et al. (2020), Section 3, final paragraphs; Section 5]

  1. Permute all failed (uneliminated) columns to the back of the matrix.
  2. OPTION A (default in SSIDS): Re-factorize failed columns using
     sequential TPP. This is effective because the failed columns have
     been updated by all successfully eliminated columns and may now
     have acceptable pivots.
  3. OPTION B: Pass uneliminated columns directly to the parent node
     in the assembly tree as "delayed pivots." The parent's frontal
     matrix is enlarged to accommodate them.
  4. Any columns still uneliminated after the fallback are delayed
     to the parent node.
```

### Threshold Parameter Guide

**Source**: Duff et al. (2020), Section 7, Table 5.

| u value | Behavior | Use case |
|---|---|---|
| 0.1 | Strict: many failures and delays, highest stability | Numerically very challenging problems |
| **0.01** (default) | **Balanced: few failures, good stability** | **General-purpose default** |
| 0.001 | Permissive: very few failures, fastest factorization | Easy indefinite or nearly positive-definite problems |

Experimental data from Table 5 of Duff et al. (2020) shows the impact:

- **pct20stif**: u=0.1 gives 1974 failed columns, 239 delays; u=0.001 gives 5 failed, 2 delays.
- **sparsine**: u=0.1 gives 9547 failed columns (61.3s); u=0.001 gives 104 failed (5.38s).
- **Ga10As10H30**: u=0.1 gives 1157 failed; u=0.001 gives 0 failed. The factorization time
  difference (22.69s vs 18.69s) is driven entirely by the fallback TPP cost.

The backward error degrades by roughly one order of magnitude when going from u=0.1 to u=0.001,
which is acceptable for most applications.

---

## 6. Symbolic Analysis

**Source**: Gilbert, Ng, Peyton (1994), Sections 2--3; Liu (1992), Sections 3--4.

Before numeric factorization, the symbolic analysis phase determines the sparsity structure
of the factor L. This consists of three key computations:

### 6.1 Elimination Tree Computation

**Source**: Liu (1992), Section 3.1; Gilbert et al. (1994), Section 1.2.

The elimination tree T(A) has n nodes {1, ..., n} where node p is the parent of j if and only if:

```
p = min{ i > j | l_{ij} != 0 }
```

That is, the parent of column j is the row index of the first off-diagonal nonzero in column j
of L. The elimination tree can be computed in O(m * alpha(m, n)) time using the algorithm
of Liu (1990), which uses disjoint set union (UNION-FIND) operations.

```
Elimination Tree Algorithm:
  [Liu (1990/1992), Section 3.1]

  Input: Sparse symmetric matrix A (n x n, m off-diagonal nonzeros)
  Output: parent[j] for j = 1..n

  Initialize: ancestor[j] = j for all j   -- disjoint set forest

  for j = 1 to n do
      for each i in {k > j : a_{kj} != 0} do    -- higher-numbered neighbors
          r = FIND(i)                            -- find representative
          if r != j then
              parent[r] = j                      -- set parent in elimination tree
              ancestor[r] = j                    -- union by attaching to j
          end if
      end for
  end for
```

**Properties of the elimination tree** (Liu 1992, Theorems 3.1--3.2):
- If k is a descendant of j in T, then the structure of column k (below row j) is
  contained in the structure of column j.
- If l_{jk} != 0 with k < j, then k is a descendant of j in T.

### 6.2 Row and Column Counts

**Source**: Gilbert, Ng, Peyton (1994), Section 2 (Figures 2--3), Section 3.

The column count cc(v) -- the number of nonzeros in column v of L -- determines the front
sizes, the total nonzero count |L|, and the operation count (sum of squared column counts).
Gilbert, Ng, and Peyton (1994) give an O(m * alpha(m, n)) algorithm that computes both
row and column counts without forming L explicitly.

**Key idea**: Decompose each row subtree T_r[u] into paths using the Least Common Ancestor (LCA)
of consecutive lower-numbered neighbors. The row count for vertex u is:

```
rc(u) = 1 + sum_{i=1}^{k-1} ( level(p_i) - level(lca(p_i, p_{i+1})) )
```

where p_1 < p_2 < ... < p_{k-1} are the lower-numbered neighbors of u in G(A) union T(A),
p_k = u, and level(v) is the depth of v in T.

For column counts, define weights wt(v) such that cc(v) = sum of wt(x) over all descendants
x of v. The algorithm computes these weights incrementally:

```
Row and Column Count Algorithm:
  [Gilbert, Ng, Peyton (1994), Figure 3]

  Input: Sparse symmetric matrix A, elimination tree T with postordering
  Output: rc[u] (row counts), cc[v] (column counts) for all vertices

  Compute level(u) for all u (depth in T from root)
  Compute fst_desc(u) for all u (first descendant in postorder)
  Initialize: prev_p[u] = 0, prev_nbr[u] = 0 for all u
  Initialize: rc[u] = 1 for all u
  Initialize: wt[u] = 1 if u is a leaf, 0 otherwise

  for p = 1 to n do                       -- postorder traversal
      if p != root then
          wt[parent(p)] -= 1
      end if

      for each u in hadj[p] do            -- higher-numbered neighbors
          if fst_desc(p) > prev_nbr[u] then    -- skeleton graph test
              wt[p] += 1
              p' = prev_p[u]
              if p' == 0 then
                  rc[u] += level(p) - level(u)
              else
                  q = FIND(p')             -- q = lca(p', p)
                  rc[u] += level(p) - level(q)
                  wt[q] -= 1
              end if
              prev_p[u] = p
          end if
          prev_nbr[u] = p
      end for

      UNION(p, parent(p))
  end for

  -- Compute column counts by summing weights bottom-up
  cc[v] = wt[v] for all v
  for v = 1 to n-1 do
      cc[parent(v)] += cc[v]
  end for
```

**Complexity**: O(m * alpha(m, n)) where m = |E(G(A))| and alpha is the inverse Ackermann
function. In practice, alpha(m, n) <= 4 for all feasible problem sizes, so this is effectively
linear.

The "supernodal" version (with the skeleton graph test on lines marked with `fst_desc`)
reduces the number of FIND operations from m to m^{-} (the number of skeleton graph edges),
which is often much smaller.

### 6.3 Fill-Reducing Ordering

Before computing the elimination tree, a fill-reducing permutation P is computed. Common choices:

- **AMD** (Approximate Minimum Degree): O(m) heuristic, good for 2D problems.
- **Nested dissection** (e.g., METIS): O(m log m), better for 3D problems and parallelism.

The ordering determines P in PAP^T = LDL^T. The symbolic analysis then operates on the
reordered matrix. The faer library provides AMD and COLAMD implementations that can be
leveraged directly.

---

## 7. Multifrontal Structure

**Source**: Liu (1992), Sections 3--5 (formal definitions of frontal/update matrices, assembly
tree, extend-add operator); Duff, Hogg, Lopez (2020), Section 5 (SSIDS multifrontal
implementation); Hogg, Ovtchinnikov, Scott (2016), Section 2.3 (factorize phase).

### 7.1 Three-Phase Solver Structure

The complete solver follows the standard three-phase approach:

```
Phase 1: ANALYSE (symbolic analysis)
  [Liu (1992), Sections 3--4; Gilbert et al. (1994)]

  Input: Sparsity pattern of A
  Output: Elimination tree, assembly tree, column counts, memory layout

  1. Compute fill-reducing ordering P (AMD or nested dissection).
  2. Compute elimination tree T(PAP^T).
  3. Compute column counts using Gilbert-Ng-Peyton algorithm.
  4. Detect supernodes (groups of columns with identical structure).
  5. Amalgamate small supernodes to increase dense block sizes.
  6. Allocate storage for factors based on column counts.

  The result of the analyse phase is reusable for any matrix with the
  same sparsity pattern.


Phase 2: FACTORIZE (numeric factorization)
  [Duff et al. (2020), Section 5; Liu (1992), Algorithm 4.1]

  Input: Numeric values of A, symbolic analysis result
  Output: L, D factors, permutation P

  Traverse the assembly tree from leaves to root:
  for each node f in topological order do
      1. assemble_pre(f):
           - Allocate memory for fully summed columns and contribution block.
           - Copy delayed columns from children.
           - Assemble contributions from children into fully summed columns.

      2. factor(f):
           - Apply APTP algorithm (Algorithm 3.1) to the fully summed columns.
           - Compute the contribution block.

      3. assemble_post(f):
           - Assemble contributions from children into the contribution block.
  end for


Phase 3: SOLVE (triangular solve)
  [Liu (1992), implicit in the factorization structure; Hogg et al. (2016), Section 2.4]

  Input: L, D factors, right-hand side b
  Output: Solution x

  Forward solve (leaves to root):
      For each node, solve L_node * y_node = b_node (triangular solve)
      then apply subdiagonal block (matrix-vector multiply).

  Diagonal solve:
      z = D^{-1} * y   (trivial for 1x1 blocks; 2x2 blocks solved directly)

  Backward solve (root to leaves):
      For each node, solve L_node^T * x_node = z_node (transpose triangular solve).
```

### 7.2 Frontal Matrix Assembly

**Source**: Liu (1992), Section 4, Theorem 4.1.

At each node j of the assembly tree, the frontal matrix F_j is assembled from:
1. Original matrix entries A_{*j} (the column of A corresponding to node j).
2. Update (contribution) matrices U_{c1}, ..., U_{cs} from child nodes c1, ..., cs.

The assembly uses the **extend-add** operator: given two dense matrices with potentially
different index sets, extend each to the union of the index sets (filling with zeros) and add:

```
F_j = (original entries from A_{*j}) extend-add U_{c1} extend-add ... extend-add U_{cs}
```

The frontal matrix has the block structure:

```
         n_fs columns    n_cb columns
       +-------------+-------------+
       |             |             |
n_fs   |    F_11     |   F_21^T    |   Fully summed part
rows   |             |             |
       +-------------+-------------+
       |             |             |
n_cb   |    F_21     |    F_22     |   Contribution block (not fully summed)
rows   |             |             |
       +-------------+-------------+

  F_11: n_fs x n_fs symmetric, fully summed -- pivots chosen from here
  F_21: n_cb x n_fs, fully summed -- updated but not eliminated here
  F_22: n_cb x n_cb, partially summed -- contribution block passed to parent
```

Pivots can only be chosen from F_11 (the fully summed part). After eliminating q <= n_fs
pivots, the contribution block C = F_22 - F_21 * D_1^{-1} * F_21^T is passed to the parent.

### 7.3 How APTP Integrates into the Multifrontal Framework

**Source**: Duff et al. (2020), Section 5, Figure 3.

At each node, the `factor(f)` routine applies APTP (Algorithm 3.1) to the fully summed
columns of the frontal matrix. The integration works as follows:

```
APTP within Multifrontal:

  1. The fully summed part [F_11; F_21] has m rows and p = n_fs columns.

  2. APTP partitions this m x p region into nb x nb blocks and runs
     Algorithm 3.1. The diagonal blocks come from F_11; the subdiagonal
     blocks come from F_21.

  3. After APTP, q <= p pivots have been successfully eliminated.
     The remaining p - q columns are "delayed":
       a. First, attempt re-factorization of failed columns using TPP fallback.
       b. Any columns that still fail are passed to the parent as delayed pivots.

  4. The contribution block is formed:
       C = F_22 - L_2 * D_1 * L_2^T
     where L_2 contains the subdiagonal factor entries corresponding to the
     contribution block rows, and D_1 is the block diagonal of the q
     successfully eliminated pivots.

  5. Delayed columns are appended to the parent's fully summed part,
     enlarging the parent's frontal matrix.
```

The parallel implementation uses OpenMP tasks (Duff et al. 2020, Figure 3):

```
Parallel Multifrontal Factorization:
  [Duff et al. (2020), Section 5, Figure 3]

  #pragma omp taskgroup
  {
      for f = 0 to nfronts-1 do
          #pragma omp task depend(inout: front[f]) depend(in: parent[f])
          {
              assemble_pre(f, children[f])
              factor(f)          // uses APTP with nested OpenMP tasks
              assemble_post(f, children[f])
          }
      end for
  }
```

Tree parallelism arises because independent branches can execute simultaneously. Node
parallelism arises from the 2D block structure within APTP, where Factor, Apply, and Update
tasks at the same node can overlap.

### 7.4 Update Matrix Management

**Source**: Liu (1992), Section 5.

Update matrices are managed using a stack when the assembly tree is traversed in postorder.
Each update matrix U_j is pushed onto the stack after it is computed. When forming frontal
matrix F_j, the children's update matrices are popped from the stack. Postorder traversal
guarantees LIFO (last-in-first-out) access, minimizing working memory.

---

## 8. Source Citations

Every pseudocode block and algorithm in this document is traceable to specific sections
of the following academic papers:

### Primary References

| Paper | Citation | Sections Used |
|---|---|---|
| Duff, Hogg, Lopez (2020) | "A new sparse LDL^T solver using a posteriori threshold pivoting." SIAM J. Sci. Comput. 42(2), C23--C42. DOI: 10.1137/18M1225963 | **Sections 2--5**: Algorithm 3.1 (main APTP loop), all kernel descriptions, Algorithm 4.1 (complete pivoting), stability proofs, multifrontal integration, Table 5 (threshold guide) |
| Hogg, Ovtchinnikov, Scott (2016) | "A sparse symmetric indefinite direct solver for GPU architectures." ACM Trans. Math. Softw. 42(1), Article 1. DOI: 10.1145/2756548 | **Sections 2.3--3.2**: Original a posteriori pivoting on GPU, frontal matrix structure, tile-column factorization, threshold test, backup/restore mechanism |
| Gilbert, Ng, Peyton (1994) | "An efficient algorithm to compute row and column counts for sparse Cholesky factorization." SIAM J. Matrix Anal. Appl. 15(4), 1075--1091. | **Section 2**: Row count via path decomposition and LCA; **Section 3**: Column count via subtree weight summation; Figure 3 (combined algorithm) |
| Liu (1992) | "The multifrontal method for sparse matrix solution: theory and practice." SIAM Review 34(1), 82--109. | **Sections 3--5**: Elimination tree, frontal/update matrix definitions, extend-add operator, assembly tree, Algorithm 4.1 (Cholesky via frontal matrices), stack management |
| Schenk, Gartner (2006) | "On fast factorization pivoting methods for sparse symmetric indefinite systems." Electron. Trans. Numer. Anal. 23, 158--179. | **Section 2**: SBK algorithm (Algorithm 2.1), pivot perturbation strategy, comparison point for APTP |

### Cross-Reference: Pseudocode to Source

| Pseudocode | Source |
|---|---|
| Algorithm 3.1 (APTP main loop) | Duff et al. (2020), Section 3, Algorithm 3.1 |
| Factor kernel | Duff et al. (2020), Section 3 (Factor description); Section 4 (inner pivoting) |
| ApplyN kernel | Duff et al. (2020), Section 3, items 4--6 |
| ApplyT kernel | Duff et al. (2020), Section 3, items 7--9 |
| Adjust kernel | Duff et al. (2020), Section 3, Adjust description |
| UpdateNN kernel | Duff et al. (2020), Section 3, items 10--11 |
| UpdateNT kernel | Duff et al. (2020), Section 3, items 12--13 |
| UpdateTN kernel | Duff et al. (2020), Section 3, UpdateTN description |
| Algorithm 4.1 (complete pivoting) | Duff et al. (2020), Section 4, Algorithm 4.1 |
| Stability bounds for Algorithm 4.1 | Duff et al. (2020), Section 4.1, Equations (4.1)--(4.7) |
| Fail-in-place protocol | Duff et al. (2020), Section 3 (paragraphs 1--2); Hogg et al. (2016), Section 3.2 |
| Threshold parameter guide (Table 5) | Duff et al. (2020), Section 7, Table 5 |
| Elimination tree algorithm | Liu (1992), Section 3.1; Gilbert et al. (1994), Section 1.2 |
| Row/column count algorithm (Figure 3) | Gilbert, Ng, Peyton (1994), Section 2 (row counts), Section 2.3 (column counts), Figure 3 |
| Extend-add assembly | Liu (1992), Section 4, Theorem 4.1 |
| Multifrontal factorization loop | Liu (1992), Algorithm 4.1; Duff et al. (2020), Section 5, Figure 3 |
| Three-phase solver structure | Hogg et al. (2016), Section 2; Duff et al. (2020), Section 5 |
| SBK algorithm (for comparison) | Schenk & Gartner (2006), Algorithm 2.1 |
