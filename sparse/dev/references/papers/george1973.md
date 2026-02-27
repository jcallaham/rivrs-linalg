# Nested Dissection of a Regular Finite Element Mesh 

Alan George

SIAM Journal on Numerical Analysis, Vol. 10, No. 2 (Apr., 1973), 345-363.

## Stable URL:

http://links.jstor.org/sici?sici=0036-1429\(197304\)10\%3A2\<345\%3ANDOARF\>2.0.CO\%3B2-D
SIAM Journal on Numerical Analysis is currently published by Society for Industrial and Applied Mathematics.

Your use of the JSTOR archive indicates your acceptance of JSTOR's Terms and Conditions of Use, available at http://www.jstor.org/about/terms.html. JSTOR's Terms and Conditions of Use provides, in part, that unless you have obtained prior permission, you may not download an entire issue of a journal or multiple copies of articles, and you may use content in the JSTOR archive only for your personal, non-commercial use.

Please contact the publisher regarding any further use of this work. Publisher contact information may be obtained at http://www.jstor.org/journals/siam.html.

Each copy of any part of a JSTOR transmission must contain the same copyright notice that appears on the screen or printed page of such transmission.

JSTOR is an independent not-for-profit organization dedicated to creating and preserving a digital archive of scholarly journals. For more information regarding JSTOR, please contact support@jstor.org.

# NESTED DISSECTION OF A REGULAR FINITE ELEMENT MESH* 

## ALAN GEORGE †


#### Abstract

Let $M$ be a mesh consisting of $n^{2}$ squares called elements, formed by subdividing the unit square $(0,1) \times(0,1)$ into $n^{2}$ small squares of side $1 / h$, and having a node at each of the $(n+1)^{2}$ grid points. With $M$ we associate the $N \times N$ symmetric positive definite system $A x=b$, where $N=(n+1)^{2}$, each $x_{i}$ is associated with a node of $M$, and $A_{i j} \neq 0$ if and only if $x_{i}$ and $x_{j}$ are associated with nodes of the same element. If we solve the equations via the standard symmetric factorization of $A$, then $O\left(n^{4}\right)$ arithmetic operations are required if the usual row by row (banded) numbering scheme is used, and the storage required is $O\left(n^{3}\right)$. In this paper we present an unusual numbering of the mesh (unknowns) and show that if we avoid operating on zeros, the $L D L^{T}$ factorization of $A$ can be computed using the same standard algorithm in $O\left(n^{3}\right)$ arithmetic operations. Furthermore, the storage required is only $O\left(n^{2} \log _{2} n\right)$. Finally, we prove that all orderings of the mesh must yield an operation count of at least $O\left(n^{3}\right)$, provided we use the standard factorization algorithm.


1. Introduction. It is well known that if we avoid operating on and storing zeros, the way we number or order the unknowns of a sparse system of equations can drastically affect the amount of computation and storage required for their direct solution. In this paper we consider this ordering problem for a finite element system of equations associated with a regular $n \times n$ mesh or grid. This sparse system has a very definite structure which we show can be exploited to considerable advantage. The problem we consider is special, but as we indicate in our concluding remarks, the techniques we use can be applied in more general situations.

We first define a finite element system of equations. Let $M$ be any mesh formed by subdividing a planar region $R$ with boundary $\partial R$ by a number of arcs, all of which terminate on an arc or on $\partial R$. The mesh so formed consists of the union of subregions which we call elements. We require that $M$ have a node at each vertex in the mesh, and it may have nodes on edges and in the interior of some or all of the elements. (A vertex of $M$ is any point in $R \cup \partial R$ having more than two arcs emanating from it.) We refer to these nodes as vertex, edge, and interior nodes respectively. We call such a mesh a finite element mesh.

Let $M$ have $N$ nodes, numbered in some way from 1 to $N$. Associating an unknown $x_{i}$ with the $i$ th node leads us to the following.

Definition. A finite element system of equations associated with the finite element mesh $M$ is any $N \times N$ symmetric positive definite system $A x=b$ having the property that $A_{i j} \neq 0$ if and only if $x_{i}$ and $x_{j}$ are associated with nodes of the same element of $M$.

The reader familiar with the use of finite element techniques may wonder why we allow $M$ to be more general than the meshes usually employed in finite element methods. Usually $M$ is the union of triangular and/or quadrilateral elements with adjacent elements having a common side or vertex. Our intention

[^0]is to introduce a sequence of meshes which correspond (in the sense above) to a sequence of matrices which arise during the computation. These derived meshes have a less restricted topology than the original $M$.

Since there is a $1-1$ correspondence between nodes and unknowns, we do not distinguish between them in the sequel.

In this paper we consider the ordering problem for the finite element system

$$
A x=b,
$$

associated with the mesh $M$ formed by subdividing the unit square $R=(0,1) \times(0,1)$ into $n^{2}$ small squares (elements) having side length $1 / h$. The mesh has a node at each of the $N=(n+1)^{2}$ vertices. The method we use to solve (1.1) is direct; we compute the symmetric factorization $L D L^{T}$ of $A$, where $L$ is unit lower triangular and $D$ is a positive diagonal matrix. We then obtain $x$ by solving $L y=b, D z=y$, and finally $L^{T} x=z$.

Our measure of computational difficulty for the solution of (1.1) is $\theta$, the number of multiplicative operations (multiplications and divisions) required to factor $A$. We regard this as a reasonable measure, since the required number of additions and subtractions is about the same, and the factorization is typically the major portion of the computation. In the following, "operations" will mean multiplicative operations. We measure storage requirements by $\eta$, the number of nonzero off-diagonal components in $L$.

Consider any matrix $A^{(5)}$ obtained from a matrix $A$ in our class by setting $A_{i j}$ to zero unless $x_{i}$ and $x_{j}$ are associated with nodes on the extremity of the same element edge. Such a matrix arises when we apply the usual 5 -point difference operator in connection with solving self-adjoint elliptic boundary problems with rectangular domains [3]. Hoffman, Martin and Rose [6] show that symmetric positive definite matrices having the structure of $A^{(5)}$ require at least $O\left(n^{3}\right)$ operations for their factorization, and the corresponding lower triangular factor must have at least $O\left(n^{2} \log _{2} n\right)$ nonzero components. Since $A$ is obtained from $A^{(5)}$ by adding nonzero components, it follows that these results also hold for $A$. It is important to appreciate that these results apply for the particular algorithm described in § 2, and, in this context, we conclude that the ordering we present here is optimal, in the order of magnitude sense.

Since we do not make explicit use of the actual numerical values of $A$, our upper bounds for $\theta$ and $\eta$ hold regardless of whether or not unknowns associated with the same element are indeed connected. However, to avoid tedious qualifications in the discourse that follows, we assume that if unknowns $x_{i}$ and $x_{j}$ are associated with the same element, then $A_{i j} \neq 0$. If this is not the case, our upper bounds will simply not be sharp.

An ordering similar but somewhat inferior to the one we present here appears in [4]. The bound on $\theta$ we obtain here halves the one obtained in that article.

In § 2 and § 3 we introduce quantities and a model which allow us to conveniently determine $\theta$ and $\eta$ for a given ordering (numbering) of $M$. In § 4 we present an ordering of the mesh which results in $\theta=O\left(n^{3}\right)$ and $\eta=O\left(n^{2} \log _{2} n\right)$, and in § 5 we show that if we apply the symmetric decomposition algorithm, any ordering of the mesh must result in $\theta=O\left(n^{3}\right)$.
2. Symmetric elimination. Using the outer product formulation employed by Rose [8], we describe the factorization of $A$ into $L D L^{T}$ by the following equations. Setting $A=A_{0}=B_{0}$, we have

$$
\begin{aligned}
A_{0} & =\left(\begin{array}{cc}
d_{1} & v_{1}^{T} \\
v_{1} & B_{1}^{\prime}
\end{array}\right)=\left(\begin{array}{cc}
1 & 0 \\
\frac{v_{1}}{d_{1}} & I_{N-1}
\end{array}\right)\left(\begin{array}{cc}
d_{1} & 0 \\
0 & B_{1}^{\prime}-\frac{v_{1} v_{1}^{T}}{d_{1}}
\end{array}\right)\left(\begin{array}{cc}
1 & v_{1}^{T} / d_{1} \\
0 & I_{N-1}
\end{array}\right) \\
& =L_{1}\left(\begin{array}{cc}
d_{1} & 0 \\
0 & B_{1}
\end{array}\right) L_{1}^{T}=L_{1} A_{1} L_{1}^{T}, \\
A_{1} & =\left(\begin{array}{ccc}
d_{1} & 0 & \\
0 & d_{2} & v_{2}^{T} \\
& v_{2} & B_{2}^{\prime}
\end{array}\right) \\
& =\left(\begin{array}{ccc}
1 & 0 \\
0 & \frac{v_{2}}{d_{2}} & I_{N-2}
\end{array}\right)\left(\begin{array}{ccc}
d_{1} & d_{2} & 0 \\
0 & B_{2}^{\prime}-\frac{v_{2} v_{2}^{T}}{d_{2}}
\end{array}\right)\left(\begin{array}{ccc}
1 & 0 \\
0 & 1 & v_{2}^{T} / d_{2} \\
0 & I_{N-2}
\end{array}\right) \\
& =L_{2} A_{2} L_{2}^{T}, \\
& \vdots \\
A_{N-1} & =D .
\end{aligned}
$$

Here $d_{k}$ is a positive scalar, $v_{k}$ is a vector of length $N-k$, and $B_{k}^{\prime}$ is an $(N-k) \times(N-k)$ symmetric positive definite matrix. In the sequel $B_{k}=B_{k}^{\prime}-v_{k} v_{k}^{T} / d_{k}$ is referred to as "the part of $A$ remaining to be factored" after the first $k$ steps of the factorization have been performed. Following Rose [8], we refer to performing the $k$ th step of the factorization as "eliminating variable $x_{k}$."

Since $A$ is sparse, the vectors $v_{k}$ will usually also have some zeros. Define $v_{k}$ to be the number of nonzero components in $v_{k}$. Then we have the following lemma.

Lemma 2.1 (Rose [8]). Provided we avoid operating on zeros, the number of multiplicative operations required to factor $A$ into $L D L^{T}$ is

$$
\theta=\sum_{k=1}^{N-1} \frac{v_{k}\left(v_{k}+3\right)}{2}
$$

and the number of nonzero off-diagonal components in $L$ is

$$
\eta=\sum_{i=1}^{N-1} v_{k}
$$

Proof. It is straightforward to verify from (2.1) that the cost $\theta_{k}$ of performing the $k$ th step of the factorization is $v_{k}\left(v_{k}+3\right) / 2$. Summing over $k$ and noting that $v_{N} \equiv 0$ yields (2.2).

Equation (2.3) can be derived by observing first that (2.1) implies that $A=L_{1} L_{2} \cdots L_{N-1} D L_{N-1}^{T} L_{N-2}^{T} \cdots L_{1}^{T}$, and it is easy to show that

$$
L=\sum_{k=1}^{N-1} L_{k}-(N-2) I
$$

Thus, the $k$ th column of $L$ is precisely the $k$ th column of $L_{k}$, which immediately implies (2.3). This concludes the proof.

Un-eliminated variables $x_{j}$ and $x_{k}$ are referred to as being connected if their corresponding off-diagonal components in $B_{i}(i<j, k)$ are nonzero. Obviously, unconnected variables can become connected as the factorization proceeds. Specifically, if we assume that we do not create zeros through cancellation (which is a reasonable assumption in the presence of roundoff error), then the following holds.

Lemma 2.2 (Parter [7]). The elimination of variable $x_{k}$ pairwise connects all variables $x_{i}, i>k$, to which $x_{k}$ was connected at the point of its elimination.

Proof. Referring to equations (2.1), we note that eliminating $x_{k}$ modifies $B_{k}^{\prime}$ by subtracting the rank-one matrix $v_{k} v_{k}^{T}$ from it, forming $B_{k}$. The matrix $v_{k} v_{k}^{T}$ has nonzeros in position $(i, j)$ for all $i$ and $j$ corresponding to nonzero components in $v_{k}$. Assuming no cancellation in the subtraction, $B_{k}$ must have nonzeros in the same positions, proving the lemma.

Rose [8] refers to sets of unknowns which are pairwise connected (every unknown in the set is connected to every other unknown in the set) as cliques. Thus, eliminating an unknown $x_{k}$ which is connected to a set of variables $S$ renders that set a clique.

The following lemma is a direct application of Lemmas 2.1 and 2.2.
Lemma 2.3. Let $Q=\left\{x_{i_{1}}, x_{i_{2}}, \cdots, x_{i_{q}}\right\}$ and $R=\left\{x_{j_{1}}, x_{j_{2}}, \cdots, x_{j_{r}}\right\}$ be two sets of unknowns, with unknowns $x_{k} \in Q$ connected to every unknown $x_{l} \in Q \cup R$ and no others. Then if $i_{k}<j_{l}$ for $1 \leqq k \leqq q$ and $1 \leqq l \leqq r$, the contribution to $\theta$ from elimination of the unknowns in $Q$ is

$$
m(q, r)=m_{1}(q)+m_{2}(q, r)
$$

where

$$
m_{1}(q)=\frac{q^{3}}{6}+\frac{q^{2}}{2}-\frac{2}{3} q
$$

and

$$
m_{2}(q, r)=\frac{q r}{2}(q+r+2)
$$

The storage required for columns $i_{1}, i_{2}, \cdots, i_{q}$ of $L$ is

$$
s(q, r)=\frac{1}{2} q(q+2 r-1)
$$

Proof. Since we assume variables in $Q$ are connected only among themselves and to those in $R$, we can without loss of generality assume $i_{k}=k$ and $j_{k}=q+k$. We then need only consider the cost of performing the first $q$ steps of the factorization of a $(q+r) \times(q+r)$ matrix whose first $q$ columns are dense. This implies $v_{k}=q+r-k, k=1,2, \cdots, q$, and using (2.2), we have

$$
\begin{aligned}
m(q, r) & =\sum_{k=1}^{q} \frac{1}{2}(q+r-k)(q+r-k+3) \\
& =\sum_{k=1}^{q} \frac{1}{2}(q-k)(q-k+3)+\sum_{k=1}^{q} \frac{1}{2} r(r+2 q-2 k+3) \\
& =m_{1}(q)+m_{2}(q, r)
\end{aligned}
$$

Similarly, using (2.3) along with the above formula for $v_{k}$, we have

$$
s(q, r)=\sum_{k=1}^{q}(q+r-k)=\frac{1}{2} q(q+2 r-1) .
$$

3. A mesh model for the analysis of the factorization. In the Introduction we defined the zero-nonzero structure of finite element systems of equations in terms of a planar mesh. We now establish a correspondence between elimination of certain sets of unknowns in (1.1) and corresponding changes in $M$. Denoting $M$ by $M_{0}$, we arrange that after the $k$ th step of the factorization, we have a mesh $M_{k}$ which corresponds to $B_{k}$ in the sense that unknowns still to be eliminated are connected only if they are associated with the same element of $M_{k}$. In other words, $B_{k}$ is a finite element matrix corresponding to $M_{k}$.

For example, consider the mesh $M_{0}$ depicted in Fig. 3.1, consisting of rectangular elements, and numbered as indicated. Unnumbered nodes are assumed to have numbers greater than 1 .

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-06.jpg?height=438&width=440&top_left_y=1008&top_left_x=606)
FIG. 3.1. The mesh $M_{0}$

By our definition of finite element matrices, $x_{1}$ is connected to all of the 8 other unknowns which are associated with the same elements as $x_{1}$. Keeping in mind that only unknowns associated with the same element may be connected, and recalling Lemma 2.2, the mesh $M_{1}$ which reflects the structure of $B_{1}$ must have the 8 unknowns mentioned above all associated with the same element. Thus, $M_{1}$ must be as shown in Fig. 3.2, where a new element has been formed by merging four elements of $M$.

Consider a second mesh example (Fig. 3.3), which has some edge nodes as well as vertex nodes. Here we assume unnumbered nodes have numbers greater than 2 . Which mesh $M_{1}$ correctly reflects the structure of $B_{1}$, the part of the matrix remaining after the first step of the factorization is complete? By the same arguments as before, the 8 unknowns $x_{k}, k>2$, associated with the two elements sharing nodes 1 and 2 must be associated with the same element in $M_{1}$. Thus, these two adjacent elements must coalesce. Node (unknown) 2 is already connected to all unknowns of the two elements, and eliminating $x_{1}$ does not connect it to any others. Thus, $B_{1}$ 's structure is correctly described by the mesh in Fig. 3.4, where node 2 is now interior to the newly formed element. Since elimination of

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-07.jpg?height=437&width=444&top_left_y=290&top_left_x=613)
Fig. 3.2. Derived mesh $M_{1}$

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-07.jpg?height=440&width=442&top_left_y=893&top_left_x=611)
Fig. 3.3. Mesh M ${ }_{0}$

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-07.jpg?height=436&width=440&top_left_y=1484&top_left_x=611)
Fig. 3.4. Mesh $M_{1}$ derived from the mesh $M_{0}$ of Fig. 3.3

variable $x_{2}$ only involves variables to which $x_{2}$ is already connected, the mesh $M_{2}$ derived from $M_{1}$ and corresponding to $B_{2}$ would simply be the mesh of Fig. 3.4 with node 2 removed.

Finally, to completely fix these ideas in mind, consider the mesh $M_{0}$ depicted in Fig. 3.5. Unnumbered nodes are assumed to have numbers greater than 9 . Using our examples developed above, the reader should now be able to verify

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-08.jpg?height=454&width=589&top_left_y=311&top_left_x=516)
Fig. 3.5. Mesh $M_{0}$

that the meshes $M_{1}, M_{2}, \cdots, M_{9}$ shown in Fig. 3.6 correspond to $B_{1}, B_{2}, \cdots, B_{9}$ respectively.

Summarizing, we have developed the following mesh transformation rules.
Rule 1. The elimination of a variable associated with an interior node corresponds to removal of that node from the mesh; for example, the transformation from $M_{7}$ to $M_{8}$ in Fig. 3.6.

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-08.jpg?height=986&width=896&top_left_y=1233&top_left_x=362)
Fig. 3.6. Meshes $M_{1}$ to $M_{9}$ corresponding to $B_{1}$ through $B_{9}$

Rule 2. The elimination of a variable associated with an edge node on an edge $e$ corresponds to removal of the node and edge from the mesh, but other nodes on $e$ remain, and become interior nodes of a new element; for example, the transformation of $M_{6}$ to $M_{7}$ in Fig. 3.6.

Rule 3. The elimination of a variable associated with a vertex node corresponds to removal of the vertex node and all its incident element edges. Other nodes lying on these edges remain, and become interior nodes of a newly formed element.

In order to preserve our element structure when eliminating nodes on the boundary, we must interpret our rules somewhat carefully. First observe that as far as the matrix structure is concerned, the two meshes in Fig. 3.7 are equivalent.

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-09.jpg?height=247&width=743&top_left_y=855&top_left_x=322)
Fig. 3.7. Two meshes having equivalent matrix structure

When eliminating unknowns associated with nodes on the boundary, we shall assume that edge nodes have been moved to the interior of an element and vertex nodes have been moved onto an edge as indicated in Fig. 3.7. We then apply our rules as before to this perturbed mesh, which has no nodes on $\partial R$. An example of the applications of the rules in this situation is given in Fig. 3.8.

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-09.jpg?height=702&width=917&top_left_y=1490&top_left_x=231)
FIG. 3.8. Modification of the mesh when variables associated with boundary nodes are eliminated

Thus, using this model, every numbering of the mesh determines a sequence of finite element meshes which reflects the changing structure of the part of the matrix remaining to be factored. The last member in the sequence $M_{0}, M_{1}, \cdots, M_{N}$ consists of a single boundary element whose boundary is $\partial R$, with all nodes removed.

The reader familiar with the graph theory model of elimination extensively employed by Rose [8] for analyzing sparse matrix computations will recognize the close relationship between our model and the graph model. By simply joining every pair of nodes associated with the same element of $M$, we obtain the (undirected) graph corresponding to $A$. Doing the same for our $M_{k}, k=1,2, \cdots, N$, yields a sequence of elimination graphs. Elements correspond to cliques. Thus, our model can be viewed as a graph theory model where many of the edges of the graph are not explicitly drawn, but are made implicit through the introduction of elements. For our particular application, our model has the advantage of being somewhat easier to draw and interpret. For less structured problems and in other contexts, the graph model is more appropriate.
4. A nested dissection ordering of $M$ and some crude bounds. In this section we describe an ordering of our $n \times n$ mesh $M$ which yields $\theta=O\left(n^{3}\right)$ and $\eta=O\left(n^{2} \log _{2} n\right)$. To simplify the presentation, we do not attempt to find the constants involved; the interested reader will find a careful analysis in the Appendix.

We begin by defining some sets of unknowns (nodes), where $x_{i j}$ denotes the unknown associated with node (ih, $j h$ ). We assume initially that $n=2^{l}$. For $i=1, \cdots, n-1$ define

$$
\pi(i)=p+1 \quad \text { if } i=2^{p}(2 q+1) .
$$

Furthermore, define $\pi(0)=1$ and $\pi(n)=1$. For example, when $n=16$, the values of $\pi(i), i=0, \cdots, 16$, are given in the first row of Fig. 4.1. For $k=1, \cdots, l$ define the set of nodes $P_{k}$ by

$$
P_{k}=\left\{x_{i j} \mid \max (\pi(i), \pi(j))=k\right\} .
$$

Denoting membership in $P_{k}$ by the number $k$, these node sets are depicted in Fig. 4.1 for $n=2^{4}=16$. To aid the description we have put lines around $P_{4}$ and a few subsets of $P_{1}, P_{2}$ and $P_{3}$.

Notice that $P_{4}$ subdivides the nodes (unknowns) into 4 subsets which are mutually independent in the sense that if $x_{i}$ and $x_{j}$ are in different subsets, then $A_{i j}=0$. In the same way, $P_{3}$ subdivides each of these subsets into 4 mutually independent subsets, and so on. Hence the name "nested dissection." ${ }^{1}$ Each of the $P_{i}$ themselves consists of independent sets of nodes which increase in size with $i$. For $i>1$ these subsets are "+" shaped. The subsets of each $P_{i}$ might be appropriately named "separating sets."

Now our overall strategy is to number the unknowns in $P_{1}$, followed by those in $P_{2}$ and so on, finally numbering the unknowns in $P_{l}$. The results we establish in this section are independent of the way each $P_{k}$ is numbered.

[^1]| 1 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 4 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 1 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 4 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 1 |
| 2 | 2 | 2 | 2 | 3 | 2 | 2 | 2 | 4 | 2 | 2 | 2 | 3 | 2 | 2 | 2 | 2 |
| 1 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 4 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 1 |
| 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 4 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 |
| 1 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 4 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 1 |
| 2 | 2 | 2 | 2 | 3 | 2 | 2 | 2 | 4 | 2 | 2 | 2 | 3 | 2 | 2 | 2 | 2 |
| 1 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 4 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 1 |
| 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
| 1 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 4 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 1 |
| 2 | 2 | 2 | 2 | 3 | 2 | 2 | 2 | 4 | 2 | 2 | 2 | 3 | 2 | 2 | 2 | 2 |
| 1 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 4 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 1 |
| 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 4 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 |
| 1 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 4 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 1 |
| 2 | 2 | 2 | 2 | 3 | 2 | 2 | 2 | 4 | 2 | 2 | 2 | 3 | 2 | 2 | 2 | 2 |
| 1 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 4 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 1 |
| 1 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 4 | 1 | 2 | 1 | 3 | 1 | 2 | 1 | 1 |

Fig. 4.1. The node sets $P_{k}, k=1,2,3,4$

Using our example above ( $n=16$ ), and the model we introduced in § 3, it is not difficult to see that the meshes $M_{100}, M_{196}$, and $M_{256}$ are as indicated in Figs. 4.2, 4.3 and 4.4 respectively. They correspond to the structure of the matrix remaining to be factored after unknowns in $P_{1}, P_{1} \cup P_{2}$, and $P_{1} \cup P_{2} \cup P_{3}$ respectively, have been eliminated.

Theorem 4.1. Let $n=2^{l}$ and define the node sets $P_{k}, k=1,2, \cdots, l$, by (4.1). Number the nodes in increasing order beginning with those in $P_{1}$,followed by those in $P_{2}$, etc., finally numbering those in $P_{l}$. Then there exist constants $C_{1}$ and $C_{2}$

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-11.jpg?height=699&width=707&top_left_y=1477&top_left_x=346)
FIG. 4.2. The mesh $M_{100}$ formed by elimination of unknowns in $P_{1}$

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-12.jpg?height=688&width=695&top_left_y=335&top_left_x=492)
Fig. 4.3. The mesh $M_{196}$ formed by eliminating unknowns in $P_{1}$ and $P_{2}$

such that

$$
\theta<C_{1} n^{3}
$$

and

$$
\eta<C_{2} n^{2} \log _{2} n .
$$

Proof. First observe that $P_{k}$ consists of $n^{2} / 2^{2 k}$ independent sets of unknowns, and note that they remain independent during the elimination. That is, unknowns in different subsets of $P_{k}$ never become connected during the elimination. Each independent set has no more than $2^{k+1}$ unknowns in it, and each unknown in the

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-12.jpg?height=684&width=688&top_left_y=1532&top_left_x=484)
FIG. 4.4. The mesh $M_{256}$ formed by eliminating unknowns in $P_{1}, P_{2}$ and $P_{3}$

set is connected to no more than $6 \cdot 2^{k}-3$ unknowns at the point of its elimination. For example, unknowns in $P_{2}$ are connected to at most 20 unknowns before they are eliminated. Thus, using (2.2), we have

$$
\begin{aligned}
\theta & <\sum_{k=1}^{l}\left(\frac{n^{2}}{2^{2 k}}\right) 2^{k+1}\left(6 \cdot 2^{k}\right)^{2} \\
& =C_{1}^{\prime} n^{2} \sum_{k=1}^{l} 2^{k} \leqq C_{1} n^{3}
\end{aligned}
$$

Similarly, using (2.3), we have

$$
\begin{aligned}
\eta & <\sum_{k=1}^{l} \frac{n^{2}}{2^{2 k}} \cdot 2^{k+1} \cdot 6 \cdot 2^{k}=C_{2} n^{2} \sum_{k=1}^{l} 1 \\
& =C_{2} n^{2} \log _{2} n
\end{aligned}
$$

Corollary 4.2. For any $n>2$, there exist constants $C_{3}$ and $C_{4}$ such that

$$
\theta<C_{3} n^{3}
$$

and

$$
\eta<C_{4} n^{2} \log _{2} n .
$$

Proof. Let $2^{l-1}<n<2^{l}=\bar{n}$, let $M$ and $\bar{M}$ be meshes corresponding to $n$ and $\bar{n}$, and let $A$ correspond to $M$. Now augment our system of equations (1.1) by adding trivial equations of the form $1 \cdot x_{i}=1, i>(n+1)^{2}$ corresponding to nodes of $\bar{M}$ that are not nodes of $M$, so that the dimension of our expanded coefficient matrix is $(\bar{n}+1)^{2}$. Now number $\bar{M}$ as in Theorem 4.1 and solve our expanded system. A subvector of the solution will be the solution of the original problem. Let $\alpha=\bar{n} / n<2$. Then by Theorem 4.1,

$$
\theta<C_{1} \bar{n}^{3}=C_{1} \alpha^{3} n^{3}<8 C_{1} n^{3}=C_{3} n^{3}
$$

and

$$
\begin{aligned}
\eta & <C_{2} \tilde{n}^{2} \log _{2} \bar{n}<C_{2} \alpha^{2} n^{2}\left(\log _{2} n+\log _{2} \alpha\right) \\
& <C_{2} \alpha^{2} n^{2}\left(\log _{2} n+1\right)<8 C_{2} n^{2} \log _{2} n \\
& =C_{4} n^{2} \log _{2} n
\end{aligned}
$$

As one might expect, using this crude analysis yields constants $C_{1}, C_{2}, C_{3}$ and $C_{4}$ which are, from a practical viewpoint, discouragingly large. However, in the Appendix we show that if $n=2^{l}$, then $\theta<10 n^{3}$ and $\eta<8 \ln ^{2}$.
5. A lower bound for $\theta$. In this section we show that for the algorithm described in § 2, any order of $M$ must lead to $\theta=O\left(n^{3}\right)$ or greater. Following a suggestion of D. J. Rose, our strategy is to show that some member $M_{i}$ of the mesh sequence $M_{k}, k=0,1,2, \cdots, N$, must contain an element $T$ having $n+1$ nodes associated with it. This implies that $B_{i}$ must have a dense $n+1$ submatrix in it.

Lemma 5.1. Let $M=M_{0}$ be the regular $n \times n$ finite element mesh described in the Introduction, and let $M_{k}, k=1,2, \cdots, N=(n+1)^{2}$, be the mesh sequence generated by an arbitrary ordering of $M$.

Then at least one element $T$ having $n+1$ or more unknowns associated with it appears in the mesh sequence.

Proof. Let $Q=\left(x_{i}, y_{i}\right)$ be the first node of $M$ to be removed which completely vacates a row or column of the mesh. Since all nodes are eventually removed, this situation must arise.

Suppose the removal of $Q$ vacates exactly one row \{column\} $L$. Then $L$ is contained in an element $T$, since by Rule $2, \S 3$, removal of the nodes on $L$ removes their incident edges. Moreover, the $n+1$ columns \{rows\} orthogonal to $L$ necessarily have one or more nodes remaining on them, so proceeding along them in one direction or the other from $L$ we must collide with a node lying on the boundary of $T$, proving the theorem.

Suppose the removal of $Q$ simultaneously vacates both row $R$ and column $C$ of $M$. Since all other nodes on $R$ and $C$ have been removed, $Q$ must necessarily be an interior node of an element $T$ containing both $R$ and $C$. By the same argument as above, $T$ is an element having one or more nodes in each of the $n+1$ rows and columns, again proving the theorem.

Theorem 5.2. The number of multiplicative operations required to complete the symmetric factorization of $A$ is greater than $n^{3} / 6$, regardless of the way the equations are numbered.

Proof. By Lemma 5.1, during the decomposition we necessarily create one or more elements $T_{i}$ having at least $n+1$ unknowns associated with it, and they are all connected by virtue of belonging to $T_{i}$. Thus, $B_{i}$ has a dense $(n+1) \times(n+1)$ submatrix whose symmetric factorization alone using the algorithm of § 2 requires $n^{3} / 6+n^{2} / 2-2 n / 3$ multiplicative operations.

COROLLARY 5.3. Every ordering of $M$ results in a bandwidthm $=\max _{A_{i j} \neq 0}|i-j|$ satisfying $m \geqq n$.

Proof. The $B_{i}$ of Theorem 5.2 has a dense $(n+1) \times(n+1)$ submatrix, which means the bandwidth of $B_{i}$ is at least $n$, which in turn implies the bandwidth of $A$ is at least $n$.

Corollary 5.3 indicates how unreliable bandwidth can be as a measure of computational complexity, since computation estimates based on bandwidth would suggest $\theta$ must be at least $O\left(n^{4}\right)$. It is interesting in this connection to observe that for our ordering, $m \simeq n^{2}$ !
6. Concluding remarks. We have presented an ordering of a system of $N=(n+1)^{2}$ equations $A x=b$ derived from a regular $n \times n$ mesh and have shown that using the given ordering, $A$ can always be factored in $O\left(n^{3}\right)$ multiplicative operations, and the number of nonzero components in $L$ is $O\left(n^{2} \log _{2} n\right)$. These results, combined with lower bounds of the same form due to Hoffman, Martin and Rose [6], leads us to conclude that our ordering is optimal in the order of magnitude sense (provided we use the algorithm described in § 2) with respect to both $\theta$ and $\eta$.

The class of matrices we have studied includes some well-known special cases, notably those which arise when the usual five-point or nine-point difference operator [3] is applied in connection with solving self-adjoint elliptic boundary value problems with rectangular domains. In certain special circumstances it is well known that systems with coefficient matrices having the structure of $A$ can
be solved in $O\left(n^{3}\right)$ or even $O\left(n^{2} \log _{2} n\right)$ operations [1], but these algorithms are special in the sense that they exploit the actual values of the components of $A$.

However, if the matrix has no special characteristics other than being symmetric, positive definite, and having the zero-nonzero structure corresponding to $M$, then it has been assumed that factorization of $A$ required $O\left(n^{4}\right)$ operations and $O\left(n^{3}\right)$ storage locations. (That is, the row by row ordering was used.) Since it can be shown that some iterative schemes, under suitable circumstances, will reduce the error in the solution of (1.1) below a specified tolerance in $O\left(n^{3}\right)$ operations, it is often argued that iterative schemes are more efficient than direct methods for solving these sparse systems [2], [9]. We feel that Theorem 4.1 indicates that a re-comparison of direct and iterative methods is due. Of course, iterative schemes in general require much less storage than direct methods, even when we use our ordering, but the development of large memories and fast peripheral storage devices makes storage a less important factor than it has been in the past.

The extension of our upper bound results to the regular unit cube mesh containing $n^{3}$ cuboid elements is straightforward, although a little tedious. One now numbers sets of independent unknowns enclosed in increasingly large cubes rather than squares, finally numbering three intersecting planes of nodes, each pair having a common line of nodes, and the three having one common vertex at the centroid of the mesh. For this ordering of the mesh, $\theta=O\left(n^{6}\right)$ and $\eta=O\left(n^{4}\right)$. We have been unable to extend our lower bound proof to the three-dimensional problem.

Another straightforward generalization of our results is the case where $M$ has nodes on edges and in the interior of each element, and more than one unknown is associated with each node. The results are essentially unchanged, except for some adjustments in constants. See [4] for some results in this direction, for a slightly different ordering than the one presented here.

A closely related matrix problem, arising in connection with the use of splines, has the property that grid points (unknowns) $p$ and $q$ are connected provided that the maximum difference in their $x$ - and $y$-cordinates is bounded by some number $d$ (which depends on the degree of the spline). We must now define our $P_{k}$ 's to consist of strips of horizontal and vertical grid lines, each strip consisting of $d$ parallel grid lines. In this way nodes on opposite sides of a strip cannot be connected. We then proceed as before, obtaining $\theta=O\left(n^{3}\right)$ and $\eta=O\left(n^{2} \log _{2} n\right)$.

Finally, it should be obvious that our results apply even if $A$ is unsymmetric, provided that its zero-nonzero structure is symmetric and we do not have to pivot during the decomposition to maintain numerical stability. Such matrix problems arise in connection with the use of finite element methods for solving non-self-adjoint elliptic boundary value problems. Depending on the way the problems are formulated, the matrix may or may not be symmetric, but the matrices always have symmetric structure. Furthermore, they are often diagonally dominant, which implies that pivoting for numerical stability is not required [10].

Appendix. We now decribe in detail a particular ordering using the general strategy described in §4. We assume $n=2^{l}$, and obtain sharp bounds for $\theta$ and $\eta$.

Our numbering strategy is most easily described in terms of the mesh sequence $M_{i}, i=1,2, \cdots, N$. Recall that numbering the first $k$ nodes corresponds to
carrying out the first $k$ steps of the elimination, which leaves us with a matrix $B_{k}$ remaining to be factored which has structure represented by $M_{k}$. Within each set $P_{i}$, our numbering (elimination) strategy is to number (eliminate) the node (unknown) in $M_{k}\left(B_{k}\right)$ which is connected to the fewest nodes (unknowns). When this strategy is used solely, it is known as the minimum degree algorithm. Here its application is restricted to the sequence of sets $P_{1}, P_{2}, \cdots, P_{l}$. Since the subsets of $P_{i}, 0<i \leqq l$, are independent, their relative numbering is immaterial. Therefore, we can apply the above strategy to each subset of $P_{i}$ in turn, or apply it globally to the entire set $P_{i}$. In the example in Fig. A.1, we take the former course.

It is convenient at this point to define a set $P_{0}=\left\{x_{i j} \mid i, j=0, n\right\}$, and redefine $P_{1}$ to be $P_{1}-P_{0}$. Then each subset of $P_{k}, 0<k \leqq l$, consists of a " + " shaped set of vertices, which we refer to as a hub with four spokes. (For $k=1$, most of the spokes are null.)

Now the application of the minimum degree strategy described above to each subset leads to a numbering best described by the diagrams in Figs. A.2, A. 3 and A.4, where the hashed lines indicate the boundary $\partial R$. The actual relative numbering on each edge is immaterial. The value of $p$ is undetermined and not important here.

Obviously, the number of nodes on the edges in Figs. A.2-A. 4 depend on $k$. To distinguish between the three types of subsets, we refer to those depicted in Figs. A.2-A. 4 as interior, boundary, and corner subsets respectively.

We now make repeated use of Lemma 2.3 to compute $\theta$ and $\eta$.
Lemma A.1. The contribution to $\theta$ from the elimination of the unknowns in $P_{0}$ is 36 , and the number of nonzeros in the first four columns of $L$ is 12 .

| 81 82 |  | 173 | 85 | 231 | 87 | 178 | 89 | 265 | 91 | 184 | 93 | 242 | 95 | 190 | 98 | 97 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | 83 <br> 84 | 174 | 86 | 232 | 88 | 179 | 90 | 266 | 92 | 185 | 94 | 243 | 96 | 191 | 100 | 99 |
| 171 | 172 | 176 | 175 | 233 | 181 | 182 | 183 | 267 | 187 | 188 | 189 | 244 | 196 | 195 | 193 | 192 |
| 71 | 72 | 177 | 73 | 234 | 74 | 180 | 75 | 268 | 76 | 186 | 77 | 245 | 78 | 194 | 80 | 79 |
| 227 | 228 | 229 | 230 | 238 | 237 | 236 | 235 | 269 | 256 | 255 | 254 | 253 | 249 | 248 | 247 | 246 |
| 61 | 62 | 152 | 63 | 239 | 64 | 157 | 65 | 270 | 66 | 162 | 67 | 252 | 68 | 168 | 70 | 69 |
| 149 | 150 | 153 | 151 | 240 | 155 | 158 | 156 | 271 | 160 | 163 | 161 | 251 | 167 | 169 | 166 | 165 |
| 51 | 52 | 154 | 53 | 241 | 54 | 159 | 55 | 272 | 56 | 164 | 57 | 250 | 58 | 170 | 60 | 59 |
| 273 | 274 | 275 | 276 | 277 | 278 | 279 | 280 | 281 | 282 | 283 | 284 | 285 | 286 | 287 | 288 | 289 |
| 41 | 42 | 130 | 43 | 250 | 44 | 135 | 45 | 264 | 46 | 140 | 47 | 226 | 48 | 146 | 50 | 49 |
| 127 | 128 | 131 | 129 | 206 | 133 | 136 | 134 | 263 | 138 | 141 | 139 | 225 | 145 | 147 | 144 | 143 |
| 31 | 32 | 132 | 33 | 207 | 34 | 137 | 35 | 262 | 36 | 142 | 37 | 224 | 38 | 148 | 40 | 39 |
| 201 | 202 | 203 | 204 | 208 | 209 | 210 | 211 | 261 | 220 | 221 | 222 | 233 | 219 | 218 | 217 | 216 |
| 21 | 22 | 105 | 23 | 200 | 24 | 110 | 25 | 260 | 26 | 116 | 27 | 215 | 28 | 126 | 30 | 29 |
| 103 | 104 | 106 | 107 | 199 | 111 | 112 | 113 | 259 | 117 | 118 | 119 | 214 | 124 | 125 | 123 | 122 |
| 3 | 4 | 102 | 6 | 198 | 8 | 109 | 10 | 258 | 12 | 115 | 14 | 213 | 16 | 121 | 20 | 19 |
| 1 | 2 | 101 | 5 | 197 | 7 | 108 | 9 | 257 | 11 | 114 | 13 | 212 | 15 | 120 | 18 | 17 |

FIG. A.1. Detailed numbering of $M$ for $n=16$. The set $P_{4}$ and some subsets of $P_{1}, P_{2}$ and $P_{3}$ have been outlined.

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-17.jpg?height=295&width=1062&top_left_y=341&top_left_x=168)
Fig. A.2. Mesh sequence corresponding to elimination of unknowns in a subset of $P_{k}$ which has no nodes on $\partial R$

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-17.jpg?height=283&width=1052&top_left_y=781&top_left_x=172)
FIG. A.3. Mesh sequence corresponding to elimination of unknowns of a subset of $P_{k}$ which has one node on $\partial R$

![](https://cdn.mathpix.com/cropped/a54c4909-4f6e-4732-962b-bfb3edaba055-17.jpg?height=269&width=1035&top_left_y=1237&top_left_x=173)
FIG. A.4. Mesh sequence corresponding to elimination of unknowns of a subset of $P_{k}$ having two nodes on $\partial R$

Proof. Obviously each corner node is connected to only 3 other unknowns, which means $q=1$ and $r=3$ in Lemma 2.3. Thus, the contribution to $\theta$ is $4 \cdot m(1,3)=36$, and the contribution to $\eta$ is $4 \cdot s(1,3)=12$.

Lemma A.2. The contribution to $\theta$ from the elimination of interior subsets in $P_{k}$, $0<k<l-1$, is

$$
\theta_{I}^{k}=\left(\frac{n}{2^{k}}-2\right)^{2}\left(\frac{371}{24} \cdot 2^{3 k}-17 \cdot 2^{2 k}-\frac{22}{3} \cdot 2^{k}+3\right)
$$

and the number of nonzeros in the corresponding columns of $L$ is

$$
\eta_{I}^{k}=\left(\frac{n}{2^{k}}-2\right)^{2}\left(\frac{31 \cdot 2^{2 k}}{4}-13 \cdot 2^{k}+3\right) .
$$

Proof. First observe that there are $\left(n / 2^{k}-2\right)^{2}$ such interior subsets in $P_{k}$, $0<k<l-1$. The sets $P_{0}, P_{l-1}$ and $P_{l}$ have none. The first two spokes eliminated
contribute $m\left(2^{k-1}, 3 \cdot 2^{k}\right)$, since each spoke has $2^{k-1}$ nodes on it, and the number of nodes not on the spoke but connected to unknowns on the spoke is $3 \cdot 2^{k}$. Elimination of the final two spokes and hub contributes $m\left(2^{k}-1,2^{k+2}\right)$.

Similarly, the storage required for the corresponding columns of $L$ is $2 \cdot s\left(2^{k-1}, 3 \cdot 2^{k}\right)+s\left(2^{k}-1,2^{k+2}\right)$. Expanding and simplifying these functions yields (A.1) and (A.2).

Lemma A.3. The contribution to $\theta$ from the elimination of unknowns in boundary subsets of $P_{k}, 0<k<l-1$, is

$$
\theta_{B}^{k}=4\left(\frac{n}{2^{k}}-2\right)\left(\frac{121}{12} \cdot 2^{3 k}-\frac{13}{14} \cdot 2^{2 k}-\frac{41}{6} \cdot 2^{k}+1\right)
$$

and the number of nonzeros in the corresponding columns of $L$ is

$$
\eta_{B}^{k}=4\left(\frac{n}{2^{k}}-2\right)\left(\frac{25}{4} \cdot 2^{2 k}-7 \cdot 2^{k}+1\right)
$$

Proof. Elimination of variables on the boundary spoke and the ones on the opposite spoke requires $m\left(2^{k-1}, 2^{k+1}+1\right)$ and $m\left(2^{k-1}-1,3 \cdot 2^{k}\right)$ operations respectively. The remaining variables require $m\left(2^{k}-1,3 \cdot 2^{k}+1\right)$ operations. Using this, along with the fact that there are $4\left(n / 2^{k}-2\right)$ boundary subsets in $P_{k}$, yields (A.3).

Similarly, the number of nonzeros in the corresponding columns of $L$ for each subset is

$$
s\left(2^{k-1}, 2^{k+1}+1\right)+s\left(2^{k-1}-1,3 \cdot 2^{k}\right)+s\left(2^{k}-1,3 \cdot 2^{k}+1\right)
$$

Expanding, simplifying and multiplying by $4\left(n / 2^{k}-2\right)$ yields (A.4). Note that $P_{0}, P_{l-1}$ and $P_{l}$ have no boundary subsets.

Lemma A.4. The contribution to $\theta$ from the elimination of variables in the four corner subsets of $P_{k}, 0<k<l$, is given by

$$
\theta_{C}^{k}=\frac{125}{6} \cdot 2^{3 k}+18 \cdot 2^{2 k}-\frac{34}{3} \cdot 2^{k}
$$

and the number of nonzero components in the corresponding columns of $L$ is

$$
\eta_{C}^{k}=18 \cdot 2^{2 k}-8 \cdot 2^{k}
$$

Proof. Elimination of variables of the first spoke contribute

$$
m\left(2^{k-1}, 3 \cdot 2^{k-1}+1\right)
$$

and $s\left(2^{k-1}, 3 \cdot 2^{k-1}+1\right)$ to $\theta$ and $\eta$ respectively. The second spoke contributes $m\left(2^{k-1}, 2^{k+1}+1\right)$ and $s\left(2^{k-1}, 2^{k+1}+1\right)$ to $\theta$ and $\eta$ respectively. Finally, the remaining $2^{k-1}$ unknowns contribute $m\left(2^{k}-1,2^{k+1}+1\right)$ operations and $s\left(2^{k}-1,2^{k+1}+1\right)$ nonzeros respectively. The set $P_{l}$ has no corner subsets.

Lemma A.5. The number of arithmetic operations required to eliminate the unknowns in $P_{l}$ is

$$
\theta^{l}=\frac{23}{24} n^{3}+\frac{7}{2} n^{2}+\frac{5}{8} n
$$

and the number of nonzeros in the corresponding columns of $L$ is

$$
\eta^{l}=\frac{7}{4} n^{2}+n .
$$

Proof. Elimination of the first two spokes requires $2 \cdot m(n / 2, n+1)$ operations, and the final $n+1$ unknowns require $m(n+1,0)$ operations. Similarly, the number of nonzeros is given by $2 \cdot s(n / 2, n+1)+s(n+1,0)$.

Theorem A.6. The number of multiplicative operations required to factor the finite element matrix $A$ associated with a regular $n \times n$ mesh, numbered as described above, and with $n=2^{l}$, is given by

$$
\theta=\frac{267}{28} n^{3}-17 n^{2} \log _{2} n+\frac{847}{28} n^{2}+O\left(n \log _{2} n\right),
$$

and the number of nonzero components in $L$ is

$$
\eta=\frac{93}{12} n^{2} \log _{2} n-\frac{73}{3} n^{2}+24 n \log _{2} n+O(n) .
$$

Proof. The proof consists of merely summing the quantities given in Lemmas A.1-A. 5 over the appropriate ranges. Normally, the algebra would be extremely tedious, but fortunately the author had access to the symbolic algebra system Altran [5], so the quantities were summed by machine. The expansions of the $m$ 's and $s$ 's in the above lemmas were also checked by machine.

Briefly, we compute

$$
\sum_{k=1}^{l-1} \theta_{C}^{k}+\sum_{k=1}^{l-2}\left(\theta_{I}^{k}+\theta_{B}^{k}\right)
$$

and then add $\theta^{l}$ for $P_{l}$ and 36 for $P_{0}$ yielding (A.9).
Similarly, to obtain $\eta$ we compute

$$
\sum_{k=1}^{l-1} \eta_{C}^{k}+\sum_{k=1}^{l-2}\left(\eta_{I}^{k}-\eta_{B}^{k}\right)
$$

and then add $\eta^{l}$ for $P_{l}$ and 12 for $P_{0}$.
Table A. 1 compares $\theta$ and $\eta$ for this ordering with the corresponding values ( $\bar{\theta}$ and $\bar{\eta}$ ) which result when the natural row by row numbering scheme is used.

Table A. 1
| $n$ | $\boldsymbol{N}$ | $\theta$ | $\eta$ | $\bar{\theta}$ | $\bar{\eta}$ |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 4 | 25 | 376 | 100 | 504 | 120 |
| 8 | 81 | 3,172 | 572 | 4,496 | 720 |
| 16 | 289 | 28,664 | 3,340 | 50,336 | 4,896 |
| 32 | 1,089 | 257,036 | 18,828 | 669,792 | 36,459 |


Acknowledgment. The author thanks Professors Donald J. Rose, Cleve B. Moler and R. B. Simpson for many helpful suggestions. The referees' comments are also gratefully acknowledged.

## REFERENCES

[1] F. W. Dorr, The direct solution of the discrete Poisson equation on a rectangle, SIAM Rev., 12 (1970), pp. 248-263.
[2] George Fix and Kate Larsen, Iterative methods for finite element approximations to elliptic boundary value problems, this Journal, 8 (1971), pp. 536-547.
[3] G. E. Forsythe and W. R. Wasow, Finite-Difference Methods for Partial Differential Equations, John Wiley, New York, 1960.
[4] J. A. George, Block elimination on finite element systems of equations, Sparse Matrices and their Applications, D. J. Rose and R. A. Willoughby, eds., Plenum Press, New York, 1972.
[5] A. D. Hall, The Altran system for rational function manipulation-a survey, Comm. ACM, 14 (1971), pp. 517-522.
[6] A. J. Hoffman, M. S. Martin and D. J. Rose, Complexity bounds for regular finite difference and finite element grids, this Journal, 10 (1973), pp. 364-369.
[7] S. V. Parter, The use of linear graphs in Gauss elimination, SIAM Rev., 3 (1961), pp. 119-130.
[8] D. J. Rose, A graph-theoretic study of the numerical solution of sparse positive definitive systems of linear equations, Graph Theory and Computing, R. Read, ed., Academic Press, New York, 1972.
[9] R. S. Varga, Matrix Iterative Analysis, Prentice-Hall, Englewood Cliffs, N.J., 1962.
[10] J. H. Wilkinson, The Algebraic Eigenvalue Problem, Clarendon Press, Oxford, 1965.


[^0]:    * Received by the editors March 17, 1972, and in revised form July 25, 1972.
    † This work was supported in part by a University of Waterloo research grant, and in part by Canadian National Research Council Grant A8111. Alan George received his Ph.D. in Computer Science in 1971 from Stanford University under the joint direction of Professor Forsythe and Dr. Fred Dorr. He is now Assistant Professor of Computer Science at the University of Waterloo, Waterloo, Ontario, Canada.

[^1]:    ${ }^{1}$ The author is indebted to Professor Garret Birkhoff of Harvard University for coining this very appropriate term.

