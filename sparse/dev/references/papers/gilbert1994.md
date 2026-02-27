# AN EFFICIENT ALGORITHM TO COMPUTE ROW AND COLUMN COUNTS FOR SPARSE CHOLESKY FACTORIZATION* 

JOHN R. GILBERT ${ }^{\dagger}$, ESMOND G. NG ${ }^{\ddagger}$, AND BARRY W. PEYTON ${ }^{\ddagger}$


#### Abstract

Let an undirected graph $G$ be given, along with a specified depth-first spanning tree $T$. Almost-linear-time algorithms are given to solve the following two problems. First, for every vertex $v$, compute the number of descendants $w$ of $v$ for which some descendant of $w$ is adjacent (in $G$ ) to $v$. Second, for every vertex $v$, compute the number of ancestors of $v$ that are adjacent (in $G$ ) to at least one descendant of $v$.

These problems arise in Cholesky and $Q R$ factorizations of sparse matrices. The authors' algorithms can be used to determine the number of nonzero entries in each row and column of the triangular factor of a matrix from the zero/nonzero structure of the matrix. Such a prediction makes storage allocation for sparse matrix factorizations more efficient. The authors' algorithms run in time linear in the size of the input times a slowly growing inverse of Ackermann's function. The best previously known algorithms for these problems ran in time linear in the sum of the nonzero counts, which is usually much larger. Experimental results are given demonstrating the practical efficiency of the new algorithms.


Key words. sparse Cholesky factorization, sparse $Q R$ factorization, symbolic factorization, graph algorithms, chordal graph completion, disjoint set union, column counts, row counts

AMS subject classifications. 65F50, 68Q20

1. Introduction. Direct solution of a sparse symmetric positive definite linear system requires four steps [7], [15]: reordering, symbolic factorization, sparse Cholesky factorization, and sparse triangular solutions. Let $A$ be the $n \times n$ coefficient matrix of the linear system after it has been reordered to reduce fill, and let $L$ be the lower triangular Cholesky factor of $A$. This paper presents improved algorithms for computing the number of nonzero entries in each row and column of $L$ prior to the symbolic factorization step. We refer to these parameters as the row counts and column counts of $L$.

In least squares computations, $A$ is $m \times n$, with $m \geq n$. It is often necessary to compute the orthogonal factorization $A=Q R$. Our algorithms can be used also to predict upper bounds on the row counts and column counts of the upper triangular factor $R$, since the structure of $R$ is always contained in the structure of the Cholesky factor of $A^{T} A$ [12].

Throughout the paper we assume familiarity with graphs, trees, and such basic techniques as depth-first search [24]. We also assume a basic knowledge of the four steps in solving sparse systems by Cholesky factorization, and with the use of graphs in these algorithms [15]. More specifically, we assume familiarity with elimination trees [19], skeleton graphs [18], postorderings, supernodes [1], [2], [16], [20], [21], and the subscript compression scheme for $L$ [15], [25].

[^0]1.1. Applications. Here we survey some of the sparse matrix settings in which it is useful to precompute the row counts, the column counts, or the total number of nonzeros in the Cholesky factor of a sparse matrix.

Either the row or column counts can be used to compute $|L|$, the total number of nonzeros in the factor. (We write $|X|$ for the number of nonzeros in a matrix $X$ or the number of elements in a set $X$.) Knowing $|L|$ before the numeric factorization step makes it possible to allocate storage all at once instead of dynamically. In sparse Cholesky factorization, the time required to compute $|L|$ by existing methods is dominated by the time required for numerical factorization; but there are at least two settings in which it is valuable to be able to compute $|L|$ as fast as possible.

First, some methods for large-scale numerical optimization use Cholesky factorization on a Hessian matrix [5], [6]. If the Hessian is indefinite, Cholesky factorization will abort, but the partial factorization contains enough information to help determine a good descent direction containing negative curvature information. In this case, the symbolic factorization time may dominate the time spent on the numeric factorization before it aborts. Thus it may be more efficient to skip the symbolic phase and to build the data structure for $L$ during the numeric factorization. However, for this to be efficient, we still need to find $|L|$ (and perhaps the column counts) before starting the factorization.

Second, much research remains to be done on the issue of how best to reorder the initial matrix to reduce fill, i.e., to reduce $|L|$. It is sometimes useful to compute $|L|$ for many different orderings of the same matrix, both in experiments with reordering algorithms and when trying to optimize an ordering for a specific matrix. Our new algorithms make this much faster.

Besides fill, there are several other measures of the quality of a reordering. Some of them can be computed from the column counts; for example, the total number of arithmetic operations is the sum of the squares of the column counts, and the maximum front size is equal to the largest column count. The smallest maximum front size, over all reorderings of a graph, is one more than the graph's treewidth [3]. Thus the fast column count algorithm may also be useful in experimental studies of treewidth.

Two applications related to the supernodal structure of $L$ also require the column counts. Supernodes are clusters of columns with related nonzero patterns, which can be exploited to use fast dense matrix computation kernels in sparse factorization; §3 describes them in more detail. First, there is a simple, flexible $O(n)$ scheme for computing supernode partitions [2], [17] that takes the column counts and the elimination tree as input. This algorithm is more versatile and faster than the $O(|A|)$ algorithm of Liu, Ng, and Peyton [20], which takes the original matrix and its elimination tree as input. The latter algorithm computes the so-called fundamental supernode partition. Given a fast algorithm to compute column counts, the more flexible scheme could be used efficiently to compute coarser supernode partitions [2], which trade extra fill for a simpler sparsity structure that can be used to improve efficiency on vector supercomputers or to reduce synchronization overhead on shared-memory multiprocessors.

The second supernodal application of the column counts is to compute the storage required for indexing information for $L$ in the usual compressed format generated by the symbolic factorization step [25]. Current software packages [4], [9] do not precompute the space needed for this compressed symbolic factorization, because it is too expensive using the currently known algorithms. The storage required for the other three steps in the solution process is usually computed in advance; we believe
that the new algorithms introduced here are efficient enough to be used by a software package to precompute the storage requirement of the symbolic factorization step as well.

Finally, we know of only one application that specifically requires the row counts rather than the column counts. The row counts are the numbers of column modifications (sparse SAXPYs) required to complete each column in sparse Cholesky factorization algorithms. Some parallel implementations [13], [14] need the row counts to determine when all the modifications have arrived for each column.
1.2. Previous work. Like many combinatorial algorithms in sparse matrix factorization, all the efficient algorithms for row and column counts begin by computing the elimination tree of the matrix (defined in the next section). The fastest known elimination tree algorithm is due to Liu [19]. The time complexity for this algorithm is dominated by disjoint set union operations, which take time $O(m \alpha(m, n))$, where $A$ is $n \times n$ and has $2 m$ off-diagonal nonzeros. Here $\alpha(m, n)$ is a slowly growing inverse of Ackermann's function defined by Tarjan [27]; for all values of $m$ and $n$ less than the number of elementary particles in the observable universe, $\alpha(m, n) \leq 4$. Thus a function that is $O(m \alpha(m, n))$ is often called "almost linear."

The fastest previously known algorithm for computing row and column counts is also due to Liu [19]. It first computes the elimination tree of $A$ and then traverses each "row subtree" of the elimination tree (defined in the next section). The total size of the row subtrees is the number of nonzeros in the factor, so the running time of this step is $O(|L|)$. Unless the factor is extremely sparse, the subtree traversals dominate the time to find the elimination tree. To put this in perspective, suppose $A$ is the matrix of an $n$-node finite difference mesh ordered by nested dissection. Then $m$ is $O(n)$, and $|L|$ is $O(n \log n)$ in two dimensions or $O\left(n^{4 / 3}\right)$ in three dimensions.

The algorithm in this paper also takes $A$ and the elimination tree as input but runs in almost-linear time $O(m \alpha(m, n))$; the time complexity for the new algorithm is dominated by disjoint set union operations. Thus it computes the row and column counts in the same asymptotic time needed to find the elimination tree. As we will see in §4, this asymptotic efficiency is also reflected in practice.
1.3. Outline. Section 2 presents the row and column count algorithm from a graph-theoretic point of view. Here it is convenient to think of the input not as the graph $G(A)$ of a matrix, but as the graph $G(A) \cup T(A)$ that has edges both for the matrix nonzeros and for the elimination tree. (The elimination tree $T(A)$ usually has edges not contained in $G(A)$.) The elimination tree is a depth-first spanning tree of the graph $G(A) \cup T(A)$; thus for the purpose of the high-level view in §2, the input is just an undirected graph with a specified depth-first spanning tree. In this setting, we suspect that our results may be useful in efficient algorithms involving chordal graphs, chordal completion, and treewidth.

In §3 we return to the matrix-computation point of view, and discuss details of the implementation in the sparse matrix setting. Two points of practical importance arise here: we modify the algorithm slightly to make only one pass over its input, and we take advantage of supernodal structure to compute only with a subgraph called the skeleton graph. We show how to organize the entire computation, including the skeleton graph reduction, within the framework of the fundamental supernode algorithm of Liu, Ng, and Peyton [20].

Section 4 contains experimental results. We experiment with both the nodal and supernodal versions of the algorithm, as well as with several implementations of the disjoint set union operations (UNION and FIND) that dominate the asymptotic
running time. The best version is the supernodal algorithm with path-halving and no union by rank (definitions are in §3.3); it performs well enough that we argue it should be a standard part of high-performance sparse factorization codes. Finally, §5 contains concluding remarks.

## 2. The algorithm.

2.1. Definitions and problem statement. Let $G=(V, E)$ be a connected undirected graph with $n$ vertices and $m$ edges, and let $T$ be a specific depth-first spanning tree for $G$ (e.g., $G=G(A) \cup T(A)$ and $T=T(A)$ ). We call vertices $v$ and $w$ adjacent if they are joined by an edge in $G$; that is, if $(v, w) \in E$. We say that vertex $v$ is an ancestor of vertex $w$ if $v$ is on the path in $T$ from $w$ to the root of $T$. Vertex $v$ is a descendant of $w$ if $w$ is an ancestor of $v$. Note that a vertex is its own ancestor and its own descendant; a proper ancestor or descendant is one that is different from the vertex itself. We write $T[v]$ for the set of descendants of $v$ and also for the subtree of $T$ (rooted at $v$ ) that those vertices induce.

Since $T$ is a depth-first spanning tree, every edge of $G$ (whether or not it is an edge of $T$ ) joins an ancestor in $T$ to a descendant in $T$ [24].

To simplify notation, we assume that the vertices of $G$ are the integers 1 through $n$. We also assume that the vertex numbers are a postorder on $T$; that is, that for every vertex $v$, the vertices of $T[v]$ are numbered consecutively, with $v$ numbered last. Thus vertex $n$ is the root of $T$.

The level of vertex $v$, which we write level(v), is its distance in $T$ from the root. The least common ancestor of vertices $v$ and $w$, which we write $l c a(v, w)$, is the ancestor of $v$ and $w$ with the smallest postorder number (or the largest level). Both a postorder numbering and the vertex levels for an arbitrary tree can be computed in linear time by depth-first search [26]. Given a set of $k$ pairs $\{v, w\}$ of vertices, the $k$ least common ancestors $l c a(v, w)$ can be computed in $O(k \alpha(k, n))$ time, where $\alpha$ is the very slowly growing inverse of Ackermann's function mentioned above [28]. We describe these algorithms in more detail in §3.

We consider the following two problems.
Problem 1 (row counts). For every node $u \in V$, let row $[u]$ be the set of descendants $v$ of $u$ for which either $v=u$ or there exists an edge $(u, w)$ with $w \in T[v]$. The problem is to compute $r c(u)=|r o w[u]|$ for every $u$.

Problem 2 (column counts). For every node $v \in V$, let $\operatorname{col}[v]$ be the set of ancestors $u$ of $v$ for which either $u=v$ or there exists an edge $(u, w)$ with $w \in T[v]$. The problem is to compute $c c(v)=|\operatorname{col}[v]|$ for every $v$.

Note that $v \in \operatorname{row}[u]$ if and only if $u \in \operatorname{col}[v]$, and that $u$ is an element of both $\operatorname{row}[u]$ and $\operatorname{col}[u]$. For each $u$, the subgraph of $T$ induced by row $[u]$, denoted by $T_{r}[u]$ and referred to as the row subtree of $u$, is connected; it is a "pruned subtree" rooted at $u$. The subgraph of $T$ induced by $\operatorname{col}[v]$ may not be connected.

We conclude by briefly describing the relationship between these problems and sparse Cholesky factorization. It may seem a bit confusing that we include the elimination tree edges in the graph $G$ in the graph problem but not in the matrix problem; however, the answer is the same in either case.

Let an $n \times n$ symmetric, positive definite matrix $A$ be given, and let $G(A)$ be its undirected graph (whose vertices are the integers 1 through $n$ ). Let $G^{+}(A)$ be the filled graph of $G(A)$ [22] obtained by adding to $G(A)$ edge ( $v, w$ ) whenever there is a path in $G(A)$ from $v$ to $w$ whose intermediate vertices are all smaller than both $v$ and $w$. The graph $G^{+}(A)$ is chordal, and (ignoring numerical cancellation) is the graph of $L+L^{T}$, where $L$ is the Cholesky factor of $A$ [23].

![](https://cdn.mathpix.com/cropped/10fc0455-f575-497d-8a59-f8116a7c4f80-05.jpg?height=571&width=634&top_left_y=197&top_left_x=519)
Fig. 1. Example of path decomposition.

The elimination tree of $A$, denoted $T(A)$, has vertices 1 through $n$, and the parent of vertex $v$ is the smallest $w>v$ such that $(v, w)$ is an edge of $G^{+}(A)$. Liu [19] surveys the uses and properties of this structure. It is a forest with one tree for each connected component of $G(A)$; if $A$ is irreducible then $T(A)$ is a tree. The elimination tree may not be a subgraph of $G(A)$, but it is a subgraph of $G^{+}(A)$, and in fact it is a depth-first spanning tree of that graph. If $A^{\prime}$ is a matrix whose graph is $G\left(A^{\prime}\right)=G(A) \cup T(A)$, it is straightforward that $G^{+}\left(A^{\prime}\right)=G^{+}(A)$ and $T\left(A^{\prime}\right)=T(A)$.

Now consider problems (1) and (2) above for $G=G\left(A^{\prime}\right)$ and $T=T\left(A^{\prime}\right)$. It is easy to show [19] that the edges of $G^{+}(A)=G^{+}\left(A^{\prime}\right)$ are exactly those ( $u, v$ ) for which $v \neq u$ and $v \in \operatorname{row}[u]$ (or $u \in \operatorname{col}[v]$ ). Thus $r c(u)$ is the number of nonzeros in row $u$ of the Cholesky factor $L$ of $A$, and $c c(v)$ is the number of nonzeros in column $v$ of $L$.
2.2. Row counts. We count the vertices in $r o w[u]$ by counting the edges in the pruned subtree $T_{r}[u]$ of $T$ that row $[u]$ induces. The following lemma lets us partition those edges into paths.

Lemma 2.1. Let $p_{1}<p_{2}<\cdots<p_{k}$ be some of the vertices of a rooted tree $R$ (where $<$ is postorder), and suppose all the leaves and the root of $R$ are among the $p_{i}$ 's. Let $q_{i}$ be the least common ancestor of $p_{i}$ and $p_{i+1}$ for $1 \leq i<k$. Then each edge $(s, t)$ of the tree is on the tree path from $p_{j}$ to $q_{j}$ for exactly one $j$.

Proof. Suppose $t$ is the parent of $s$ in $R$. The descendants of $s$ include at least one leaf, so they include at least one $p_{i}$. Let $p_{j}$ be the largest $p_{i}$ among the descendants of $s$. Then $p_{j} \leq s<p_{j+1}$. (There must be a $p_{j+1}$ - that is, we cannot have $j=k$ because $p_{k}$ is the root, which is a proper ancestor of $s$.) Since $s$ is an ancestor of $p_{j}$ but not of $p_{j+1}$, the least common ancestor $q_{j}$ of $p_{j}$ and $p_{j+1}$ is a proper ancestor of $s$, and hence an ancestor of $t$. Therefore ( $s, t$ ) is on the path from $p_{j}$ to $q_{j}$.

Now consider an $i \neq j$. If $s$ is not an ancestor of $p_{i}$, then ( $s, t$ ) is not on the path from $p_{i}$ to its ancestor $q_{i}$. If $s$ is an ancestor of $p_{i}$, then $p_{i} \leq s$, and $i \neq j$ implies $p_{i} \leq p_{i+1} \leq s$. Since postorder assigns consecutive numbers to the vertices in a subtree, this means that $s$ is also an ancestor of $p_{i+1}$, and hence of the least common ancestor $q_{i}$. Thus ( $s, t$ ) is not on the path from $p_{i}$ to $q_{i}$.

Figure 1 shows an example of the path decomposition.

Recall that $T$ is a depth-first spanning tree of $G$ and hence every edge of $G$ joins an ancestor in $T$ to a descendant in $T$. Now consider a vertex $u$ of $G$. If the lowernumbered neighbors of $u$ in $G$ are $p_{1}<p_{2}<\cdots<p_{k-1}$, and if $p_{k}=u$, then the pruned subtree $R=T_{r}[u]$ induced by row[u] satisfies the hypotheses of Lemma 2.1. Thus the number of edges in $T_{r}[u]$ is the sum of the lengths of the paths in the lemma. The length of the path from $p_{i}$ to its ancestor $q_{i}$ is the difference of their levels. The number of vertices in row $[u]$ is one more than the number of edges, so

$$
r c(u)=1+\sum_{1 \leq i<k}\left(\operatorname{level}\left(p_{i}\right)-\operatorname{level}\left(\operatorname{lca}\left(p_{i}, p_{i+1}\right)\right)\right) .
$$

(Here $l c a$ and level are taken in $T$ rather than $T_{r}[u]$, but it is clear that for any two vertices in row $[u]$ the least common ancestor and the difference in levels are the same in either tree.)

Let ladj $[u]$ be the lower numbered neighbors of $u$ in $G$. The algorithm to compute $r c(u)$ for all $u$ first sorts each set $l a d j[u] \cup\{u\}$ by postorder, then computes all the necessary least common ancestors, and finally computes the sum above for each $u$. Computing level numbers (and the postorder itself if necessary) takes linear time, and sorting the sets $\operatorname{ladj}[u] \cup\{u\}$ into postorder takes linear time by a lexicographic bucket sort. There is one least-common-ancestor computation for each edge of $G$, so the dominant term in the algorithm's time complexity is $O(m \alpha(m, n))$.
2.3. Column counts. Because $u \in \operatorname{col}[v]$ if and only if $v \in r o w[u]$, the column count $c c(v)$ is equal to the number of row subtrees $T_{r}[u]$ that contain $v$. We could compute $c c(v)$ by traversing each row subtree in turn, and counting the number of times each vertex was traversed [19]. This, however, would take time proportional to $\sum_{v} c c(v)$.

To get a faster algorithm, we define weights $w t(v)$ on the vertices of $G$ in such a way that the column count for vertex $v$ turns out to be the sum of the weights of the descendants of $v$. The key observation is that we can compute these weights as a sum of contributions from each row subtree, and that the row subtree contributions can be computed efficiently using the same least common ancestors as in the row count algorithm.

Here are the details. For each vertex $u$, define $\chi_{u}$ to be the characteristic function of $r o w[u]$, so that $\chi_{u}(v)=1$ if $v \in r o w[u]$ and $\chi_{u}(v)=0$ otherwise. Define $w t_{u}$ by

$$
w t_{u}(v)=\chi_{u}(v)-\sum_{\text {children } y \text { of } v} \chi_{u}(y)
$$

These weights may be positive, negative, or zero. This definition implies that

$$
\chi_{u}(v)=\sum_{x \in T[v]} w t_{u}(x)
$$

In a sense, $w t_{u}$ is a "first difference" down the tree of the characteristic function of row $[u]$. Finally, define

$$
w t(v)=\sum_{u \in V} w t_{u}(v)
$$

Now we prove three lemmas relating the column counts to the weights, the weights to the sets row $[u]$, and finally the row $[u]$, once more, to the least common ancestors.

Lemma 2.2. For every vertex $v$,

$$
c c(v)=\sum_{x \in T[v]} w t(x)
$$

Proof. Because $v \in \operatorname{row}[u]$ if and only if $u \in \operatorname{col}[v]$, we have

$$
c c(v)=|\operatorname{col}[v]|=\sum_{u \in V} \chi_{u}(v)
$$

Equation (2) states that this is equal to

$$
\sum_{u \in V} \sum_{x \in T[v]} w t_{u}(x) .
$$

The result follows by reversing the order of summation and using (3).
Lemma 2.2 implies that we can compute the column counts easily and efficiently from the weights by traversing the tree in postorder and summing the weights of the subtrees. It remains to describe how to compute the weights.

Lemma 2.3. Let $u$ and $v$ be vertices. Suppose that $d$ of the children of $v$ are vertices of row $[u]$. Then

$$
w t_{u}(v)=\left\{\begin{array}{cl}
1-d & \text { if } v \in \text { row }[u] \\
-1 & \text { if } v \text { is the parent of } u \\
0 & \text { otherwise }
\end{array}\right.
$$

Proof. This is immediate from (1) and the definition of $\chi_{u}$.
Lemma 2.3 implies that the only vertices $v$ for which $w t_{u}(v)$ is nonzero are the leaves of the pruned row subtree $T_{r}[u]$, the internal vertices of $T_{r}[u]$ that have more than one child in $T_{r}[u]$, and the parent of $u$. The following lemma allows us to compute $w t_{u}(v)$ for each $v$ from the same $p_{i}$ 's and $q_{i}$ 's we used in the row count algorithm.

Lemma 2.4. Let $p_{1}<p_{2}<\cdots<p_{k}$ be some of the vertices of a rooted tree $R$ (where $<$ is postorder), and suppose all the leaves and the root of $R$ are among the $p_{i}$ 's. Let $q_{i}$ be the least common ancestor of $p_{i}$ and $p_{i+1}$, for $1 \leq i<k$. Then for each vertex $v$ of $R$, the number of children of $v$ in $R$ is

$$
\left|\left\{i: q_{i}=v\right\}\right|-\left|\left\{i: p_{i}=v\right\}\right|+1
$$

Proof. Let $Q=\left|\left\{i: q_{i}=v\right\}\right|$, let $P=\left|\left\{i: p_{i}=v\right\}\right|$, and let $d$ be the number of children of $v$ in $R$. Consider the set of directed paths from $p_{i}$ to $q_{i}$ in $R$, for $1 \leq i<k$. For any collection of directed paths, each path that includes vertex $v$ either begins at $v$ or enters $v$ along edges from other vertices. Similarly, each path that includes vertex $v$ either ends at $v$ or leaves $v$ along edges to other vertices. Consequently:

The number of paths that either begin at $v$ or enter $v$ along edges from other vertices must be equal to the number of paths that either end at $v$ or leave $v$ along edges to other vertices.
(This is essentially Kirchoff's law for a flow of unit size from $p_{i}$ to $q_{i}$ for each $i$.) Lemma 2.1 says that every edge of $R$ is on exactly one of these paths. Therefore one path enters $v$ from each of the $d$ children of $v$; exactly one path leaves $v$, to its parent, unless $v$ is the root; one path begins at $v$ for each $i$ such that $p_{i}=v$ (except for $i=k$

```
Sort the vertices and their lists of neighbors by a postorder of $T$;
Compute level ( $u$ ) as the distance from $u$ to $n$ (the root), for all $u$;
Compute $l c a\left(p, p^{\prime}\right)$ for every $p$ and its successor $p^{\prime}$ in $l a d j[u] \cup\{u\}$, for all $u$;
$r c(u) \leftarrow 1$, for all $u$;
$w t(u) \leftarrow 1$, for all $u$;
for $u \leftarrow 1$ to $n$ do
    if $u \neq n$ then
        $w t($ parent $(u)) \leftarrow w t($ parent $(u))-1$;
    end if
    for $p \in \operatorname{ladj}[u]$ (in order) do
        $w t(p) \leftarrow w t(p)+1 ;$
        $p^{\prime} \leftarrow$ the successor of $p$ in ladj $[u] \cup\{u\}$;
        $q \leftarrow l c a\left(p, p^{\prime}\right) ;$
        $r c(u) \leftarrow r c(u)+\operatorname{level}(p)-\operatorname{level}(q) ;$
        $w t(q) \leftarrow w t(q)-1 ;$
    end for
end for
$c c(v) \leftarrow w t(v)$, for all $v ;$
for $v \leftarrow 1$ to $n-1$ do
    $c c(\operatorname{parent}(v)) \leftarrow c c(\operatorname{parent}(v))+c c(v) ;$
end for
```

Fig. 2. Algorithm to compute row and column counts.
if $v$ is the root); and one path ends at $v$ for each $i$ such that $q_{i}=v$. A trivial path with $p_{i}=q_{i}=v$ both starts and ends at $v$, but does not enter or leave $v$. Thus the relation above is

$$
P+d=Q+1
$$

if $v$ is not the root of $R$, or

$$
(P-1)+d=Q+0
$$

if $v$ is the root. In either case, we have $d=Q-P+1$ as desired. $\square$

Now consider a vertex $u$ of $G$. If the vertices of ladj $[u]$ are $p_{1}<p_{2}<\cdots<p_{k-1}$, and if $p_{k}=u$, then the pruned subtree $R=T_{r}[u]$ induced by $r o w[u]$ satisfies the hypotheses of Lemma 2.4. Therefore, using Lemma 2.3, if $v$ is a vertex of row[u] then $w t_{u}(v)=\left|\left\{i: p_{i}=v\right\}\right|-\left|\left\{i: q_{i}=v\right\}\right|$. Thus we could compute $w t_{u}(v)$ for all $v$ by initializing each weight to zero, setting the weight of the parent of $u$ to -1 , and then adding one to the weight of each $p_{i}$ and subtracting one from the weight of each $q_{i}$.

In fact we do not need to compute $w t_{u}(v)$ separately for each $u$; we can compute $w t(v)=\sum_{u} w t_{u}(v)$ all at once. The algorithm begins, like the row count algorithm, by sorting each set ladj $[u] \cup\{u\}$ in postorder and computing all the necessary least common ancestors. It initializes $w t(u)$ to one for each $u$. Then, for each $u$, it subtracts one from the weight of the parent of $u$, adds one to $w t(p)$ for each $p \in \operatorname{ladj}[u]$, and subtracts one from $w t(q)$ for the least common ancestor $q$ of each pair $p$ and $p^{\prime}$ of consecutive members (in postorder) of $\operatorname{ladj}[u] \cup\{u\}$. Finally, the algorithm computes $c c(v)$ for all $v$ by summing the weights of each subtree in postorder. Figure 2 sketches the algorithm to compute both row and column counts. The only step
that takes more than linear time is the least common ancestor computation, and the dominant term in the algorithm's time complexity is $O(m \alpha(m, n))$.
3. Implementation. The discussion in the previous section was in a general graph-theoretic setting. However, to obtain the most efficient implementation of the new algorithm for our applications, we need to switch back to a sparse matrix setting.

Consider a symmetric matrix $A$ and its graph $G(A)$. Assume that the elimination tree $T(A)$, the postordering, and the values level(u) (with respect to $T(A)$ ) have been computed, as required in Fig. 2. Two other requirements must be met to obtain a practical and efficient implementation of the new algorithm.

First, we must reorganize the computation to avoid sorting the adjacency lists by postorder and precomputing all the least common ancestors. Indeed, direct implementation of the algorithm in Fig. 2 would require that $G(A)$ be processed three times, and we doubt that any multiple-pass implementation will come close to realizing the practical efficiency of the single-pass implementation presented in this section.

Second, we must discard some edges of $G(A)$ that do not affect the result. Recall from Liu [18] that the skeleton graph $G^{-}(A)$ is obtained from $G(A)$ by removing every edge ( $u, v$ ) for which $v<u$ and the vertex $v$ is not a leaf of $T_{r}[u]$. The skeleton graph is the smallest subgraph of $G(A)$ whose filled graph is identical with that of $G(A)$. Consequently, the new algorithm produces the same results when applied to $G^{-}(A)$ as when applied to $G(A)$. Indeed, if $G=G^{-}(A) \cup T(A)$ rather than $G=G(A) \cup T(A)$ in Lemmas 2.1 and 2.4, then every vertex $p_{1}, p_{2}, \ldots, p_{k-1}$ is a leaf in the tree $R$. This reduces the number of edges searched and least common ancestors computed by the new algorithm to the minimum possible. Since $G^{-}(A)$ often has far fewer edges than $G(A)$ in practice, an implementation that processes $G^{-}(A)$ rather than $G(A)$ promises to be substantially faster; we see in §4 that this is indeed the case.

The skeleton graph $G^{-}(A)$ is closely related to fundamental supernodes of $A$, and can be computed efficiently in linear time by a simple modification of the algorithm of Liu, Ng, and Peyton [20] to find fundamental supernodes. Indeed, that algorithm is a good framework for implementing our new algorithm, whether the skeleton graph is exploited or not. We can combine the two algorithms to obtain an efficient singlepass implementation. As this implementation processes the edges of $G(A)$, it discards edges not in the skeleton graph, and uses only the skeleton edges to compute the data for the row and column counts. If $m^{-}$is the number of edges in $G^{-}(A)$, then this scheme runs in $O\left(m+m^{-} \alpha\left(m^{-}, n\right)\right)$ time.

Section 3.1 below reviews the material we need from Liu, Ng, and Peyton [20]. Section 3.2 presents a detailed version of the new combined implementation. Section 3.3 briefly describes our implementation of the disjoint set union algorithm for computing the least common ancestors, upon which the time complexity of our algorithm depends.
3.1. A fast algorithm for finding supernodes. Liu, Ng, and Peyton [20] introduced an $O(|A|)$ algorithm to compute a fundamental supernode partition. Their algorithm assumes that the elimination tree $T(A)$ has been computed and that the vertices are numbered by a postordering of $T(A)$. Let the higher adjacency set of $v$, denoted by hadj[v], be the set of neighbors of $v$ in $G(A)$ that are numbered higher than $v$, and let $h a d j^{+}[v]$ be the higher adjacency set of $v$ in $G^{+}(A)$. Ashcraft and Grimes [2] defined a fundamental supernode as a maximal contiguous set of vertices $\{v, v+1, \ldots, v+s\}$ such that $v+i$ is the only child of $v+i+1$ in the elimination
tree (for $i=0,1, \ldots, s-1$ ) and

$$
\operatorname{hadj}^{+}[v]=\operatorname{hadj}^{+}[v+s] \cup\{v+1, v+2, \ldots, v+s\}
$$

The fundamental supernodes partition the vertices of $G(A)$.
In matrix terms, a supernode is any group of consecutive columns in $L$ with a full diagonal block and with identical column patterns below the diagonal block. A fundamental supernode is maximal subject to the following condition: every column of the supernode except the last is an only child in the elimination tree. Liu, Ng, and Peyton [20] give several reasons why fundamental supernodes are the most appropriate choice of supernodes for most applications, one of which is that they are independent of the choice of postordering for $T(A)$.

Finding the set of fundamental supernodes is equivalent to finding the first vertex of each supernode. These "first vertices" are characterized by the following result.

Theorem 3.1 (Liu, Ng, and Peyton [20]). Vertex $v$ is the first vertex in a fundamental supernode if and only if vertex $v$ has two or more children in the elimination tree, or $v$ is a leaf of some row subtree of $T(A)$.

The key observation is that the vertices required by the row/column count algorithm (the $p_{i}$ 's and $q_{i}$ 's) are in fact first vertices of fundamental supernodes. It follows from the discussion immediately after Lemma 2.3 in § 2.3 that the vertex pairs $p_{i}, p_{i+1}$ whose least common ancestors must be found can be restricted to vertices that are leaves of some row subtree of $T(A)$. This is equivalent to restricting the algorithm in Fig. 2 to the skeleton graph $G^{-}(A)$. Furthermore, when the $p_{i}$ 's are restricted in this manner, it is clear that every least common ancestor $q_{i}=l c a\left(p_{i}, p_{i+1}\right)$ has two or more children. Consequently, the Liu, Ng, and Peyton algorithm is an excellent vehicle for an efficient implementation of our new algorithm.
3.2. Detailed implementation of the new algorithm. The details of our single-pass, column-oriented implementation are given in Fig. 3. Note that it traverses the higher adjacency sets hadj[p] rather than the lower adjacency sets used by the algorithm in Fig. 2. Again, the vertices are numbered by a postorder of the tree $T(A)$, but here no assumption is made concerning the order of the vertices in hadj[p], nor are the least common ancestors computed in advance. Consequently, this implementation makes only a single pass through $G(A)$.

The vector of markers prev_ $p(u)$ stores the most recently visited vertex $p^{\prime}$ that is a leaf in $T_{r}[u]$. The pairs $p, p^{\prime}$ produced by the algorithm are precisely the multiset consisting of every consecutive pair of leaves in every row subtree $T_{r}[u]$. The reason for this is that one of the if tests in the algorithm screens out all edges in $G(A)$ except those in the skeleton graph $G^{-}(A)$. The lines marked with asterisks have been added to the algorithm solely for this purpose. Of these, the key line is the test for whether or not the first (i.e., lowest numbered) descendant of $p\left(\operatorname{fst\_ desc}(p)\right)$ is greater than the most recently visited vertex in ladj $[u]$, namely the vertex stored in the marker variable prev_nbr (u). It is not difficult to verify that when the condition holds true, no descendant of $p$ is adjacent to $u$ in $G(A)$; hence $p$ is indeed a leaf in $T_{r}[u]$. For full details of this test, see Liu, Ng, and Peyton [20].

The implementation is correct with or without the starred lines. We have implemented both versions: we call the one with the starred lines the supernodal version, and the one without these lines the nodal version. ${ }^{1}$ We experiment with both versions of the algorithm in our tests in §4.

[^1]Sort the vertices by a postorder of $T(A)$;
Compute level ( $u$ ) as the distance from $u$ to $n$ (the root), for all $u$;

* Compute fst_desc $(u)$ as the first (least) descendant of $u$ in $T(A)$, for all $u$; prev_p $(u) \leftarrow 0$, for all $u$;
* $\quad$ prev_nbr $(u) \leftarrow 0$, for all $u$;
$r c(u) \leftarrow 1$, for all $u$;
$w t(u) \leftarrow 0$, for all nonleaves $u$ in $T(A)$;
$w t(u) \leftarrow 1$, for all leaves $u$ in $T(A)$;
for $p \leftarrow 1$ to $n$ do
if $p \neq n$ then
$w t(\operatorname{parent}(p)) \leftarrow w t(\operatorname{parent}(p))-1 ;$
end if
for $u \in h a d j[p]$ do
if $f s t \_$desc $(p)>$ prev_nbr $(u)$ then
$w t(p) \leftarrow w t(p)+1 ;$
$p^{\prime} \leftarrow$ prev_p $(u)$;
if $p^{\prime}=0$ then
$r c(u) \leftarrow r c(u)+\operatorname{level}(p)-\operatorname{level}(u) ;$
else
$q \leftarrow \operatorname{FIND}\left(p^{\prime}\right) ;$
$r c(u) \leftarrow r c(u)+\operatorname{level}(p)-\operatorname{level}(q) ;$
$w t(q) \leftarrow w t(q)-1 ;$
end if
prev_ $p(u) \leftarrow p$;
end if
prev_nbr $(u) \leftarrow p$;
end for
UNION $(p$, parent $(p))$;
end for
$c c(v) \leftarrow w t(v)$, for all $v ;$
for $v \leftarrow 1$ to $n-1$ do
$c c(\operatorname{parent}(v)) \leftarrow c c(\operatorname{parent}(v))+c c(v) ;$
end for

Fig. 3. Implementation of algorithm to compute row and column counts.
3.3. Disjoint set union. To compute least common ancestors, the algorithm in Fig. 3 must manipulate disjoint sets of vertices, each of which induces a subtree of the elimination tree. The highest numbered vertex in each set (the root of the subtree) is used to "name" the set, and is called the representative vertex of the set. Initially each vertex $p$ from 1 to $n$ is a singleton set. As the algorithm proceeds, it executes a sequence of FIND and UNION operations which are defined as follows.
$\operatorname{FIND}(p)$ : return the representative vertex of the unique set that contains $p$.
UNION $(u, v)$ : combine the two distinct sets represented by $u$ and $v$ into a single set, which will be represented by the larger of $u$ and $v$.
It is not hard to show that the call to $\operatorname{FIND}\left(p^{\prime}\right)$ in our algorithm returns $l c a\left(p^{\prime}, p\right)$; see Tarjan [28] for details.

Each disjoint set is implemented as a tree stored using a parent vector (not to be confused with the parent vector in the elimination tree). The operation UNION( $u, v)$ joins the two distinct trees represented by $u$ and $v$ together by making one of the roots a child of the other root. Consequently, UNION is a constant-time operation. This is not the case for FIND. The operation FIND( $p$ ) traces the find path from $p$ to the root of $p$ 's tree. This root either is the representative vertex or contains a pointer to the representative vertex, depending on the implementation of UNION.

Tarjan [29] describes several techniques to shorten the find paths and thus reduce the amount of work spent on the FIND operations. Union by rank makes the shorter tree's root a child of the taller tree's root in UNION, which tends to keep the trees short and bushy. With no other enhancements, union by rank ensures that find paths are no longer than $O\left(\log _{2}(n)\right)$. This is usually combined with one of two techniques for shortening the find path during a FIND operation. The first of these is path compression, which, after finding the root, makes the parent for each vertex on the find path point to the root during a second pass along the path. Alternatively, path halving resets the parent pointer for every other vertex on the find path to point to its grandparent. Path compression shortens the find path more, but requires two passes over the find path; path halving needs only one pass.

Tarjan [27], [29] showed that when union by rank is combined with either path compression or path halving, any sequence of $n$ UNION's and $m$ FIND's takes only $O(m \alpha(m, n))$ time. Tarjan [28] pointed out how to use the disjoint set union algorithm to find the least common ancestors of an arbitrary set of pairs of vertices from the same tree; our implementation of the row and column count algorithm uses the same method. Consequently, we can implement the nodal version of our algorithm to run in $O(m \alpha(m, n))$ time, and similarly we can implement the supernodal version to run in $O\left(m+m^{-} \alpha\left(m^{-}, n\right)\right)$ time.

Gabow and Tarjan [10] showed that if the order of the UNION operations is known in advance (as is the case in our problem), then disjoint set union can be implemented so that a sequence of $n$ UNION's and $m(\geq n)$ FIND's takes only $O(m)$ time. Their sophisticated hybrid algorithm partitions the vertices into microsets and performs all the operations in a hierarchical fashion, using table look-up to answer queries within the microsets, and using the standard disjoint set union algorithm on the microsets themselves. We did not implement this algorithm; we believe its increased overhead would wipe out the difference between $O(m \alpha(m, n))$ and $O(m)$ in our application.

We implemented and tested the following six combinations.

1. No union by rank, no path compression or halving.
2. No union by rank, path compression.
3. No union by rank, path halving.
4. Union by rank, no path compression or halving.
5. Union by rank, path compression.
6. Union by rank, path halving.

We found surprisingly little difference in performance among the various options. Far more important is whether or not the row/column count processing is limited to the skeleton graph, as we see in the next section. We found that any gains due to union by rank were more than offset by the additional overhead required for its implementation. The third option-no union by rank, path halving-performed slightly better on most machines we tried. Path halving was clearly superior to path compression when the skeleton adjacency structure was not exploited. Consequently, we recommend path halving to those implementing the method, and in the next section all our timings
were obtained using path halving and no union by rank.
4. Experimental results. We ran the new algorithms on several problems from the Harwell-Boeing sparse matrix collection [8]. Table 1 lists our test problems, and Table 2 contains the problem statistics that have a bearing on the observed performance of our algorithms. Throughout this section supcnt refers to the "supernodal" version of the algorithm (Fig. 3 with the starred lines), which identifies the edges of the skeleton graph $G^{-}(A)$ and uses only those edges in its row and column count calculations, and nodcnt refers to the "nodal" version of the algorithm (Fig. 3 without the starred lines), which uses all the edges of $G(A)$.
4.1. Performance of the disjoint set union options. The primary purpose of Table 3 is to explain two things we observed in our tests: (i) why exploiting the skeleton graph is so beneficial and (ii) why the various disjoint set union (DSU) implementation options have so little influence on the performance of our code. The number of FIND operations required by nodent and supcnt is bounded above by $m$ and $m^{-}$, respectively, and bounded below by $m-n$ and $m^{-}-n$. Thus, the huge difference between the number of FIND's required by nodent and the number of FIND's required by supcnt (see Table 3) simply reflects the fact that the skeleton graph of $A$ is typically much sparser than the graph of $A$ (see Table 2).

Each FIND( $p$ ) operation traverses the find path in $p$ 's tree beginning at $p$ and ending at the root of the tree. The average number of vertices on these find paths is reported for each DSU implementation. We tested only two options for nodcnt: path compression and path halving, both without union by rank. Note that the average number of vertices on a find path ranges from 2 to 2.7 , with path compression faring slightly better than path halving. The performance of path compression suffers, however, because the find path must be traversed twice, compared with once for path halving. Our tests indicate that path halving does indeed substantially outperform path compression, and in nodcnt, where the number of FIND's is large, the gain in efficiency is substantial.

We tried all six options mentioned in §3.3 in our implementations of supcnt and, as noted earlier, we saw little difference in performance from one option to the next. The primary explanation for this phenomenon is the small proportion of supcnt's total work devoted to DSU operations. The number of FIND operations is small relative to $m$, and the average number of vertices on a find path is small (from 1.4 to 2.6) for five of the six options tested. For the sixth option (no DSU enhancements), the average number of vertices on a find path is still quite modest (from 3.6 to 5.8), with less work required for each vertex visited. Consequently, even this option is competitive in our tests.

When path compression or path halving is used, union by rank obtains only modest reductions in the average number of nodes visited. The overhead costs associated with union by rank more than offset any advantages conferred by the technique. Comparing path compression and path halving with no union by rank, the same observations made previously for nodcnt hold for supcnt also. The primary difference is that the total work associated with DSU operations in supcnt is so small that the performance edge of path halving over path compression is quite small. Nonetheless, path halving with no union by rank has proven most effective overall and has the added advantage of simplicity. Finally, note that for our chosen option the total number of vertices visited by FIND operations is much less than $m$ for most of the test problems.

Table 1
List of test problems.
| Problem | Brief description |
| :--- | :--- |
| NASA1824 | Structure from NASA Langley, 1824 degrees of freedom |
| NASA2910 | Structure from NASA Langley, 2910 degrees of freedom |
| NASA4704 | Structure from NASA Langley, 4704 degrees of freedom |
| BCSSTK13 | Stiffness matrix-fluid flow generalized eigenvalues |
| BCSSTK14 | Stiffness matrix-roof of Omni Coliseum, Atlanta |
| BCSSTK15 | Stiffness matrix-module of an offshore platform |
| BCSSTK16 | Stiffness matrix-Corps of Engineers dam |
| BCSSTK17 | Stiffness matrix-elevated pressure vessel |
| BCSSTK18 | Stiffness matrix-R. E. Ginna nuclear power station |
| BCSSTK23 | Stiffness matrix-portion of a 3D globally triangular building |
| BCSSTK24 | Stiffness matrix-winter sports arena |


Table 2
Problem statistics.
| Problem | Dimension $n$ | Edges in $G(A) m$ | Edges in $G^{-}(A) m^{-}$ | Edges in $G^{+}(A) m^{+}$ |
| :--- | :--- | :--- | :--- | :--- |
| NASA1824 | 1824 | 18692 | 3565 | 71875 |
| NASA2910 | 2910 | 85693 | 8113 | 201493 |
| NASA4704 | 4704 | 50026 | 9672 | 276768 |
| BCSSTK13 | 2003 | 40940 | 5598 | 269668 |
| BCSSTK14 | 1806 | 30824 | 4352 | 110461 |
| BCSSTK15 | 3948 | 56934 | 13186 | 647274 |
| BCSSTK16 | 4884 | 142747 | 11665 | 736294 |
| BCSSTK17 | 10974 | 208838 | 24569 | 994885 |
| BCSSTK18 | 11948 | 68571 | 23510 | 650777 |
| BCSSTK23 | 3134 | 21022 | 8500 | 417177 |
| BCSSTK24 | 3562 | 78174 | 6977 | 275360 |


Table 3
Average number of vertices on a find path for DSU implementation options: PC is path compression, PH is path halving, R is union by rank, and NR is no union by rank.
| Problem | nodcnt |  |  | supcnt |  |  |  |  |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | vertices path |  | FIND's | vertices |  |  |  |  |  | FIND's |
|  | PC | PH |  | none |  | PC |  | PH |  |  |
|  | NR |  |  | NR | R | NR | R | NR | R |  |
| NASA1824 | 2.1 | 2.3 | 17050 | 4.1 | 1.9 | 2.3 | 1.6 | 2.5 | 1.6 | 1923 |
| NASA2920 | 2.0 | 2.1 | 83071 | 3.6 | 1.6 | 2.1 | 1.4 | 2.2 | 1.4 | 5491 |
| NASA4704 | 2.2 | 2.3 | 45809 | 4.1 | 1.9 | 2.2 | 1.6 | 2.4 | 1.6 | 5455 |
| BCSSTK13 | 2.1 | 2.2 | 39125 | 5.5 | 2.2 | 2.2 | 1.8 | 2.4 | 1.8 | 3783 |
| BCSSTK14 | 2.1 | 2.2 | 29200 | 4.2 | 1.7 | 2.2 | 1.5 | 2.3 | 1.5 | 2728 |
| BCSSTK15 | 2.2 | 2.3 | 53468 | 4.7 | 2.0 | 2.2 | 1.6 | 2.3 | 1.7 | 9720 |
| BCSSTK16 | 2.1 | 2.1 | 138121 | 4.4 | 2.0 | 2.1 | 1.7 | 2.2 | 1.7 | 7039 |
| BCSSTK17 | 2.1 | 2.1 | 199092 | 4.1 | 1.9 | 2.2 | 1.6 | 2.2 | 1.6 | 14823 |
| BCSSTK18 | 2.3 | 2.5 | 59624 | 5.5 | 2.2 | 2.5 | 1.8 | 2.8 | 1.9 | 14563 |
| BCSSTK23 | 2.4 | 2.7 | 18419 | 5.8 | 2.4 | 2.4 | 1.9 | 2.6 | 1.9 | 5897 |
| BCSSTK24 | 2.0 | 2.1 | 74762 | 3.8 | 1.7 | 2.1 | 1.6 | 2.2 | 1.6 | 3565 |


Table 4
Run times in seconds on an IBM RS/6000 (model 320).
| Problem | E-tree | Postordering | Row/column counts |  |  | Supernodes |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  |  | Liu's lnzcnt | New |  |  |
|  |  |  |  | nodcnt | supcnt |  |
| NASA1824 | . 035 | . 006 | . 076 | . 047 | . 038 | . 031 |
| NASA2920 | . 156 | . 009 | . 256 | . 198 | . 144 | . 128 |
| NASA4704 | . 096 | . 016 | . 261 | . 128 | . 104 | . 085 |
| BCSSTK13 | . 078 | . 006 | . 238 | . 098 | . 074 | . 064 |
| BCSSTK14 | . 057 | . 005 | . 118 | . 074 | . 056 | . 048 |
| BCSSTK15 | . 108 | . 013 | . 513 | . 142 | . 113 | . 091 |
| BCSSTK16 | . 262 | . 016 | . 691 | . 331 | . 239 | . 216 |
| BCSSTK17 | . 391 | . 037 | . 965 | . 500 | . 408 | . 329 |
| BCSSTK18 | . 144 | . 040 | . 549 | . 197 | . 181 | . 141 |
| BCSSTK23 | . 044 | . 010 | . 310 | . 059 | . 054 | . 039 |
| BCSSTK24 | . 143 | . 012 | . 295 | . 184 | . 134 | . 120 |


4.2. Performance of the row and column count algorithm. We coded nodcnt and supcnt in Fortran 77 and ran our tests on an IBM RS/6000 (model 320). We used the standard Fortran compiler and compiler optimization flag (xlf -0). We used a high-resolution timer (readrtc) to obtain our timings on this machine, repeating each run ten times in succession and returning the average elapsed time. The results are shown in Table 4. We used path halving and no union by rank in the implementation of the disjoint set union algorithm for both nodcnt and supcnt. The time required to compute the elimination tree and postordering are of interest for two reasons. First, they must be computed before the row/column counts can be computed. Second, the algorithm for computing the elimination tree is, like nodcnt and supcnt, a single-pass $O(m \alpha(m, n))$ algorithm that relies on efficient implementation of the disjoint set union operations for efficiency. Thus it is interesting to compare its performance with that of the new algorithms.

Both nodent and supcnt are much more efficient than lnzcnt, the $O(|L|)$ algorithm from Liu [19]. Algorithm nodcnt is 1.29 to 5.25 times faster than lnzcnt, while supcnt is, in turn, 1.08 to 1.39 times faster than nodcnt. For every problem but one, supcnt is at least twice as fast as lnzcnt. (For NASA2920, supcnt is 1.77 times faster than lnzcnt.) For four of the problems, supcnt is more than three times faster than lnzcnt. For BCSSTK15 supcnt is 4.54 times faster, and for BCSSTK23 supcnt is 5.74 times faster.

Finally, it is interesting to compare the timings for the elimination tree algorithm [19] and the supernode algorithm [20] with those for supcnt. First, supcnt can be viewed as an extension of the supernode algorithm, and consequently the time for supcnt should be bounded below by the time for the supernode algorithm. Though there are some differences in the amount and kind of $O(n)$ work performed by the two algorithms before and after the main loop, the difference in the two timings can nevertheless be viewed as a crude measure of the cost of adding the instructions necessary to compute row and column counts to the supernode algorithm. Clearly, this cost is quite small, especially considering the simplicity and demonstrated practical efficiency of the supernode algorithm. Note also that the timings for supcnt and the elimination tree algorithm closely track each other. From these observations, we conclude that it is probably not possible to improve the performance of supcnt much beyond what we are currently observing.
5. Conclusion. We have considered in this paper the problem of predicting the row counts and column counts in the Cholesky factor $L$ of a sparse symmetric positive definite matrix $A$, given the zero/nonzero structure of $A$ and the elimination tree $T(A)$. We have presented new algorithms for determining the counts, the complexities of which are linear in $|A|$ times a slowly growing inverse of Ackermann's function; the previously known algorithms ran in $O(|L|)$ time. The key to the new algorithms is the computation of least common ancestors in a tree using the disjoint set union algorithm. We have investigated different ways of implementing the disjoint set union operations in our algorithms. Based on our experimental results, we conclude that path halving with no union by rank is the best technique for an efficient implementation of the disjoint set union algorithm.

We have further improved our new algorithms by exploiting the skeleton graph of $A$. We have demonstrated that the supernodal version is faster than the nodal version in all of the problems we tested. Moreover, both the nodal and supernodal versions are much more efficient than the previously known $O(|L|)$-time algorithms. We expect the algorithms described in this paper to be of practical use in a wide range of sparse matrix computations.

## REFERENCES

[1] C. C. Ashcraft, A Vector Implementation of the Multifrontal Method for Large Sparse, Symmetric Positive Definite Linear Systems, Tech. Report ETA-TR-51, Engineering Technology Applications Division, Boeing Computer Services, Seattle, WA, 1987.
[2] C. C. Ashcraft and R. G. Grimes, The influence of relaxed supernode partitions on the multifrontal method, ACM Trans. Math. Software, 15 (1989), pp. 291-309.
[3] H. Bodlaender, J. R. Gilbert, H. Hafsteinsson, and T. Kloks, Approximating treewidth, pathwidth, frontsize, and minimum elimination tree height, J. Algorithms, to appear.
[4] E. C. H. Chu, A. George, J. W-H. Liu, and E. G-Y. Ng, User's Guide for SPARSPAK-A: Waterloo Sparse Linear Equations Package, Tech. Report CS-84-36, University of Waterloo, Waterloo, Ontario, 1984.
[5] T. Coleman and Y. Li, On the Convergence of Reflective Newton Methods for Large-Scale Nonlinear Minimization Subject to Bounds, Tech. Report 92-1314, Cornell University Computer Science Department, Ithaca, NY, 1992.
[6] -, A Reflective Newton Method for Minimizing a Quadratic Function Subject to Bounds on the Variables, Tech. Report 92-1315, Cornell University Computer Science Department, Ithaca, NY, 1992.
[7] I. S. Duff, A. M. Erisman, and J. K. Reid, Direct Methods for Sparse Matrices, Oxford University Press, Oxford, England, 1987.
[8] I. S. Duff, R. G. Grimes, and J. G. Lewis, Sparse matrix test problems, ACM Trans. Math. Software, 15 (1989), pp. 1-14.
[9] S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman, The Yale sparse matrix package I: The symmetric codes, Internat. J. Numer. Meth. Engrg., 18 (1982), pp. 1145-1151.
[10] H. N. Gabow and R. E. Tarjan, A linear time algorithm for a special case of disjoint set union, J. Comput. Syst. Sci., 30 (1985), pp. 209-221.
[11] F. Gavril, The intersection graphs of subtrees in trees are exactly the chordal graphs, J. Combinatorial Theory B, 16 (1974), pp. 47-56.
[12] A. George and M. T. Heath, Solution of sparse linear least squares problems using Givens rotations, Linear Algebra Appl., 34 (1980), pp. 69-83.
[13] A. George, M. T. Heath, J. W-H. Liu, and E. G-Y. Ng, Solution of sparse positive definite systems on a shared memory multiprocessor, Internat. J. Parallel Programming, 15 (1986), pp. 309-325.
[14] -, Sparse Cholesky factorization on a local-memory multiprocessor, SIAM J. Sci. Statist. Comput., 9 (1988), pp. 327-340.
[15] A. George and J. W-H. Liu, Computer Solution of Large Sparse Positive Definite Systems, Prentice-Hall, Englewood Cliffs, NJ, 1981.
[16] M. T. Heath, E. NG, and B. W. Peyton, Parallel algorithms for sparse linear systems, SIAM Rev., 33 (1991), pp. 420-460.
[17] J. G. Lewis, B. W. Peyton, and A. Pothen, A fast algorithm for reordering sparse matrices for parallel factorization, SIAM J. Sci. Statist. Comput., 10 (1989), pp. 1156-1173.
[18] J. W-H. Liu, A compact row storage scheme for Cholesky factors using elimination trees, ACM Trans. Math. Software, 12 (1986), pp. 127-148.
[19] [-, The role of elimination trees in sparse factorization, SIAM J. Matrix Anal. Appl., 11 (1990), pp. 134-172.
[20] J. W-H. Liu, E. Ng, and B. W. Peyton, On finding supernodes for sparse matrix computations, SIAM J. Matrix Anal. Appl., 14 (1993), pp. 242-252.
[21] A. Pothen, Simplicial Cliques, Shortest Elimination Trees, and Supernodes in Sparse Cholesky Factorization, Tech. Report CS-88-13, Department of Computer Science, The Pennsylvania State University, University Park, PA, 1988.
[22] D. J. Rose, A graph-theoretic study of the numerical solution of sparse positive definite systems of linear equations, in Graph Theory and Computing, R. C. Read, ed., Academic Press, 1972, pp. 183-217.
[23] D. J. Rose, R. E. Tarjan, and G. S. Lueker, Algorithmic aspects of vertex elimination on graphs, SIAM J. Comput., 5 (1976), pp. 266-283.
[24] R. Sedgewick, Algorithms, Addison-Wesley, Reading, MA, 1983.
[25] A. H. Sherman, On the efficient solution of sparse systems of linear and nonlinear equations, Ph.D. thesis, Yale University, New Haven, CT, 1975.
[26] R. E. Tarjan, Depth-first search and linear graph algorithms, SIAM J. Comput., 1 (1972), pp. 146-160.
[27] -, Efficiency of a good but not linear set union algorithm, J. ACM, 22 (1975), pp. 215-225.
[28] -, Applications of path compression on balanced trees, J. ACM, 26 (1979), pp. 690-715.
[29] -, Data Structures and Network Algorithms, CBMS-NSF Regional Conference Series in Applied Math, Society for Industrial and Applied Mathematics, Philadelphia, PA, 1983.


[^0]:    * Received by the editors September 14, 1992; accepted for publication (in revised form) May 14, 1993. This work was supported in part by the Applied Mathematical Sciences Research Program, Office of Energy Research, United States Department of Energy contract DE-AC05-84OR21400 with Martin Marietta Energy Systems, Incorporated, and in part by the Institute for Mathematics and Its Applications with funds provided by the National Science Foundation.
    ${ }^{\dagger}$ Xerox Palo Alto Research Center, 3333 Coyote Hill Road, Palo Alto, California 94304-1314 (gilbert@parc.xerox.com). Copyright © 1992, 1993 by Xerox Corporation. All rights reserved.
    ${ }^{\ddagger}$ Mathematical Sciences Section, Oak Ridge National Laboratory, P. O. Box 2008, Oak Ridge, Tennessee 37831-6367 (esmond@msr.epm.ornl.gov, peyton@msr.epm.ornl.gov).

[^1]:    ${ }^{1}$ In the nodal version, $\operatorname{prev\_ } p(u)$ functions precisely as $p r e v \_n b r(u)$ does in the supernodal version.

