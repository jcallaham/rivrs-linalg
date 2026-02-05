## Chapter 7

## Unsymmetric Eigenvalue Problems

7.1 Properties and Decompositions
7.2 Perturbation Theory
7.3 Power Iterations
7.4 The Hessenberg and Real Schur Forms
7.5 The Practical QR Algorithm
7.6 Invariant Subspace Computations
7.7 The Generalized Eigenvalue Problem
7.8 Hamiltonian and Product Eigenvalue Problems
7.9 Pseudospectra

Having discussed linear equations and least squares, we now direct our attention to the third major problem area in matrix computations, the algebraic eigenvalue problem. The unsymmetric problem is considered in this chapter and the more agreeable symmetric case in the next.

Our first task is to present the decompositions of Schur and Jordan along with the basic properties of eigenvalues and invariant subspaces. The contrasting behavior of these two decompositions sets the stage for §7.2 in which we investigate how the eigenvalues and invariant subspaces of a matrix are affected by perturbation. Condition numbers are developed that permit estimation of the errors induced by roundoff.

The key algorithm of the chapter is the justly famous QR algorithm. This procedure is one of the most complex algorithms presented in the book and its development is spread over three sections. We derive the basic QR iteration in §7.3 as a natural generalization of the simple power method. The next two sections are devoted to making this basic iteration computationally feasible. This involves the introduction of the Hessenberg decomposition in §7.4 and the notion of origin shifts in §7.5.

The QR algorithm computes the real Schur form of a matrix, a canonical form that displays eigenvalues but not eigenvectors. Consequently, additional computations
usually must be performed if information regarding invariant subspaces is desired. In §7.6, which could be subtitled, "What to Do after the Real Schur Form is Calculated," we discuss various invariant subspace calculations that can be performed after the QR algorithm has done its job.

The next two sections are about Schur decomposition challenges. The generalized eigenvalue problem $A x=\lambda B x$ is the subject of §7.7. The challenge is to compute the Schur decomposition of $B^{-1} A$ without actually forming the indicated inverse or the product. The product eigenvalue problem is similar, only arbitrarily long sequences of products are considered. This is treated in §7.8 along with the Hamiltonian eigenproblem where the challenge is to compute a Schur form that has a special 2-by-2 block structure.

In the last section the important notion of pseudospectra is introduced. It is sometimes the case in unsymmetric matrix problems that traditional eigenvalue analysis fails to tell the "whole story" because the eigenvector basis is ill-conditioned. The pseudospectra framework effectively deals with this issue.

We mention that it is handy to work with complex matrices and vectors in the more theoretical passages that follow. Complex versions of the QR factorization, the singular value decomposition, and the CS decomposition surface in the discussion.

## Reading Notes

Knowledge of Chapters 1-3 and §§5.1-§5.2 are assumed. Within this chapter there are the following dependencies:
![](https://cdn.mathpix.com/cropped/59af40f2-d7a0-4751-8bf2-c2837a3551d5-02.jpg?height=132&width=1081&top_left_y=1138&top_left_x=216)

Excellent texts for the dense eigenproblem include Chatelin (EOM), Kressner (NMSE), Stewart (MAE), Stewart and Sun (MPA), Watkins (MEP), and Wilkinson (AEP).

### 7.1 Properties and Decompositions

In this section the background necessary to develop and analyze the eigenvalue algorithms that follow are surveyed. For further details, see Horn and Johnson (MA).

### 7.1.1 Eigenvalues and Invariant Subspaces

The eigenvalues of a matrix $A \in \mathbb{C}^{n \times n}$ are the $n$ roots of its characteristic polynomial $p(z)=\operatorname{det}(z I-A)$. The set of these roots is called the spectrum of $A$ and is denoted by

$$
\lambda(A)=\{z: \operatorname{det}(z I-A)=0\} .
$$

If $\lambda(A)=\left\{\lambda_{1}, \ldots, \lambda_{n}\right\}$, then

$$
\operatorname{det}(A)=\lambda_{1} \lambda_{2} \cdots \lambda_{n}
$$

and

$$
\operatorname{tr}(A)=\lambda_{1}+\cdots+\lambda_{n}
$$

where the trace function, introduced in §6.4.1, is the sum of the diagonal entries, i.e.,

$$
\operatorname{tr}(A)=\sum_{i=1}^{n} a_{i i}
$$

These characterizations of the determinant and the trace follow by looking at the constant term and the coefficient of $z^{n-1}$ in the characteristic polynomial.

Four other attributes associated with the spectrum of $A \in \mathbb{C}^{n \times n}$ include the

$$
\begin{aligned}
\text { Spectral Radius: } & \rho(A)=\max _{\lambda \in \lambda(A)}|\lambda|, \\
\text { Spectral Abscissa: } & \alpha(A)=\max _{\lambda \in \lambda(A)} \operatorname{Re}(\lambda), \\
\text { Numerical Radius: } & r(A)=\max _{\lambda \in \lambda(A)}\left\{\left|x^{H} A x\right|:\|x\|_{2}=1\right\}, \\
\text { Numerical Range }: & W(A)=\left\{x^{H} A x:\|x\|_{2}=1\right\} .
\end{aligned}
$$

The numerical range, which is sometimes referred to as the field of values, obviously includes $\lambda(A)$. It can be shown that $W(A)$ is convex.

If $\lambda \in \lambda(A)$, then the nonzero vectors $x \in \mathbb{C}^{n}$ that satisfy $A x=\lambda x$ are eigenvectors. More preciscly, $x$ is a right eigenvector for $\lambda$ if $A x=\lambda x$ and a left eigenvector if $x^{H} A=\lambda x^{H}$. Unless otherwise stated, "eigenvector" means "right eigenvector."

An eigenvector defines a 1 -dimensional subspace that is invariant with respect to premultiplication by $A$. A subspace $S \subseteq \mathbb{C}^{n}$ with the property that

$$
x \in S \Longrightarrow A x \in S
$$

is said to be invariant (for $A$ ). Note that if

$$
A X=X B, \quad B \in \mathbb{C}^{k \times k}, X \in \mathbb{C}^{n \times k}
$$

then $\operatorname{ran}(X)$ is invariant and $B y=\lambda y \Rightarrow A(X y)=\lambda(X y)$. Thus, if $X$ has full column rank, then $A X=X B$ implies that $\lambda(B) \subseteq \lambda(A)$. If $X$ is square and nonsingular, then $A$ and $B=X^{-1} A X$ are similar, $X$ is a similarity transformation, and $\lambda(A)=\lambda(B)$.

### 7.1.2 Decoupling

Many eigenvalue computations involve breaking the given problem down into a collection of smaller eigenproblems. The following result is the basis for these reductions.

Lemma 7.1.1. If $T \in \mathbb{C}^{n \times n}$ is partitioned as follows,

$$
\left.T=\underset{\substack{T_{11} \\ 0 \\ p}}{\underset{q}{T_{22}}}\right]_{q}
$$

then $\lambda(T)=\lambda\left(T_{11}\right) \cup \lambda\left(T_{22}\right)$.

Proof. Suppose

$$
T x=\left[\begin{array}{cc}
T_{11} & T_{12} \\
0 & T_{22}
\end{array}\right]\left[\begin{array}{l}
x_{1} \\
x_{2}
\end{array}\right]=\lambda\left[\begin{array}{l}
x_{1} \\
x_{2}
\end{array}\right]
$$

where $x_{1} \in \mathbb{C}^{p}$ and $x_{2} \in \mathbb{C}^{q}$. If $x_{2} \neq 0$, then $T_{22} x_{2}=\lambda x_{2}$ and so $\lambda \in \lambda\left(T_{22}\right)$. If $x_{2}=0$, then $T_{11} x_{1}=\lambda x_{1}$ and so $\lambda \in \lambda\left(T_{11}\right)$. It follows that $\lambda(T) \subset \lambda\left(T_{11}\right) \cup \lambda\left(T_{22}\right)$. But since both $\lambda(T)$ and $\lambda\left(T_{11}\right) \cup \lambda\left(T_{22}\right)$ have the same cardinality, the two sets are equal.

### 7.1.3 Basic Unitary Decompositions

By using similarity transformations, it is possible to reduce a given matrix to any one of several canonical forms. The canonical forms differ in how they display the eigenvalues and in the kind of invariant subspace information that they provide. Because of their numerical stability we begin by discussing the reductions that can be achieved with unitary similarity.

Lemma 7.1.2. If $A \in \mathbb{C}^{n \times n}, B \in \mathbb{C}^{p \times p}$, and $X \in \mathbb{C}^{n \times p}$ satisfy

$$
A X=X B, \quad \operatorname{rank}(X)=p
$$

then there exists a unitary $Q \in \mathbb{C}^{n \times n}$ such that

$$
Q^{H} A Q=T=\left[\begin{array}{cc}
T_{11} & T_{12} \\
0 & T_{22} \\
p & n-p
\end{array}\right]_{n-p}^{p}
$$

and $\lambda\left(T_{11}\right)=\lambda(A) \cap \lambda(B)$.
Proof. Let

$$
X=Q\left[\begin{array}{c}
R_{1} \\
0
\end{array}\right], \quad Q \in \mathbb{C}^{n \times n}, R_{1} \in \mathbb{C}^{p \times p}
$$

be a QR factorization of $X$. By substituting this into (7.1.5) and rearranging we have

$$
\left[\begin{array}{ll}
T_{11} & T_{12} \\
T_{21} & T_{22}
\end{array}\right]\left[\begin{array}{c}
R_{1} \\
0
\end{array}\right]=\left[\begin{array}{c}
R_{1} \\
0
\end{array}\right] B
$$

where

$$
Q^{H} A Q=\underset{p \quad n-p}{\left[\begin{array}{ll}
T_{11} & T_{12} \\
T_{21} & T_{22}
\end{array}\right]_{n-p}^{p}} \underset{n-p}{p} .
$$

By using the nonsingularity of $R_{1}$ and the equations $T_{21} R_{1}=0$ and $T_{11} R_{1}=R_{1} B$, we can conclude that $T_{21}=0$ and $\lambda\left(T_{11}\right)=\lambda(B)$. The lemma follows because from Lemma 7.1.1 we have $\lambda(A)=\lambda(T)=\lambda\left(T_{11}\right) \cup \lambda\left(T_{22}\right)$.

Lemma 7.1.2 says that a matrix can be reduced to block triangular form using unitary similarity transformations if we know one of its invariant subspaces. By induction we can readily establish the decomposition of Schur (1909).

Theorem 7.1.3 (Schur Decomposition). If $A \in \mathbb{C}^{n \times n}$, then there exists a unitary $Q \in \mathbb{C}^{n \times n}$ such that

$$
Q^{H} A Q=T=D+N
$$

where $D=\operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{n}\right)$ and $N \in \mathbb{C}^{n \times n}$ is strictly upper triangular. Furthermore, $Q$ can be chosen so that the eigenvalues $\lambda_{i}$ appear in any order along the diagonal.

Proof. The theorem obviously holds if $n=1$. Suppose it holds for all matrices of order $n-1$ or less. If $A x=\lambda x$ and $x \neq 0$, then by Lemma 7.1.2 (with $B=(\lambda)$ ) there exists a unitary $U$ such that

$$
U^{H} A U=\left[\begin{array}{cc}
\lambda & w^{H} \\
0 & C \\
1 & n-1
\end{array}\right]_{n-1}^{1}
$$

By induction there is a unitary $\tilde{U}$ such that $\tilde{U}^{H} C \tilde{U}$ is upper triangular. Thus, if $Q=U \cdot \operatorname{diag}(1, \tilde{U})$, then $Q^{H} A Q$ is upper triangular.

If $Q=\left[q_{1}|\cdots| q_{n}\right]$ is a column partitioning of the unitary matrix $Q$ in (7.1.7), then the $q_{i}$ are referred to as Schur vectors. By equating columns in the equations $A Q=Q T$, we see that the Schur vectors satisfy

$$
A q_{k}=\lambda_{k} q_{k}+\sum_{i=1}^{k-1} n_{i k} q_{i}, \quad k=1: n
$$

From this we conclude that the subspaces

$$
S_{k}=\operatorname{span}\left\{q_{1}, \ldots, q_{k}\right\}, \quad k=1: n
$$

are invariant. Moreover, it is not hard to show that if $Q_{k}=\left[q_{1}|\cdots| q_{k}\right]$, then $\lambda\left(Q_{k}^{H} A Q_{k}\right)=\left\{\lambda_{1}, \ldots, \lambda_{k}\right\}$. Since the eigenvalues in (7.1.7) can be arbitrarily ordered, it follows that there is at least one $k$-dimensional invariant subspace associated with each subset of $k$ eigenvalues. Another conclusion to be drawn from (7.1.8) is that the Schur vector $q_{k}$ is an eigenvector if and only if the $k$ th column of $N$ is zero. This turns out to be the case for $k=1: n$ whenever $A^{H} A=A A^{H}$. Matrices that satisfy this property are called normal.

Corollary 7.1.4. $A \in \mathbb{C}^{n \times n}$ is normal if and only if there exists a unitary $Q \in \mathbb{C}^{n \times n}$ such that $Q^{H} A Q=\operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{n}\right)$.

Proof. See P7.1.1.
Note that if $Q^{H} A Q=T=\operatorname{diag}\left(\lambda_{i}\right)+N$ is a Schiur decomposition of a general $n$-by- $n$ matrix $A$, then $\|N\|_{F}$ is independent of the choice of $Q$ :

$$
\|N\|_{F}^{2}=\|A\|_{F}^{2}-\sum_{i=1}^{n}\left|\lambda_{i}\right|^{2} \equiv \Delta^{2}(A) .
$$

This quantity is referred to as $A$ 's departure from normality. Thus, to make $T$ "more diagonal," it is necessary to rely on nonunitary similarity transformations.

### 7.1.4 Nonunitary Reductions

To see what is involved in nonunitary similarity reduction, we consider the block diagonalization of a 2 -by- 2 block triangular matrix.

Lemma 7.1.5. Let $T \in \mathbb{C}^{n \times n}$ be partitioned as follows:

$$
T=\underset{p}{\left[\begin{array}{cc}
T_{11} & T_{12} \\
0 & T_{22} \\
p
\end{array}\right]_{q}^{p} .}
$$

Define the linear transformation $\phi: \mathbb{C}^{p \times q} \rightarrow \mathbb{C}^{p \times q}$ by

$$
\phi(X)=T_{11} X-X T_{22}
$$

where $X \in \mathbb{C}^{p \times q}$. Then $\phi$ is nonsingular if and only if $\lambda\left(T_{11}\right) \cap \lambda\left(T_{22}\right)=\emptyset$. If $\phi$ is nonsingular and $Y$ is defined by

$$
Y=\left[\begin{array}{rr}
I_{p} & Z \\
0 & I_{q}
\end{array}\right]
$$

where $\phi(Z)=-T_{12}$, then $Y^{-1} T Y=\operatorname{diag}\left(T_{11}, T_{22}\right)$.
Proof. Suppose $\phi(X)=0$ for $X \neq 0$ and that

$$
U^{H} X V=\underset{r}{\left[\begin{array}{cc}
\Sigma_{r} & 0 \\
0 & 0
\end{array}\right]_{p-r}} \underset{p-r}{r}
$$

is the SVD of $X$ with $\Sigma_{r}=\operatorname{diag}\left(\sigma_{i}\right), r=\operatorname{rank}(X)$. Substituting this into the equation $T_{11} X=X T_{22}$ gives

$$
\left[\begin{array}{ll}
A_{11} & A_{12} \\
A_{21} & A_{22}
\end{array}\right]\left[\begin{array}{rr}
\Sigma_{r} & 0 \\
0 & 0
\end{array}\right]=\left[\begin{array}{rr}
\Sigma_{r} & 0 \\
0 & 0
\end{array}\right]\left[\begin{array}{ll}
B_{11} & B_{12} \\
B_{21} & B_{22}
\end{array}\right]
$$

where $U^{H} T_{11} U=\left(A_{i j}\right)$ and $V^{H} T_{22} V=\left(B_{i j}\right)$. By comparing blocks in this equation it is clear that $A_{21}=0, B_{12}=0$, and $\lambda\left(A_{11}\right)=\lambda\left(B_{11}\right)$. Consequently, $A_{11}$ and $B_{11}$ have an eigenvalue in common and that eigenvalue is in $\lambda\left(T_{11}\right) \cap \lambda\left(T_{22}\right)$. Thus, if $\phi$ is singular, then $T_{11}$ and $T_{22}$ have an eigenvalue in common. On the other hand, if $\lambda \in \lambda\left(T_{11}\right) \cap \lambda\left(T_{22}\right)$, then we have eigenvector equations $T_{11} x=\lambda x$ and $y^{H} T_{22}=\lambda y^{H}$. A calculation shows that $\phi\left(x y^{H}\right)=0$ confirming that $\phi$ is singular.

Finally, if $\phi$ is nonsingular, then $\phi(Z)=-T_{12}$ has a solution and

$$
Y^{-1} T Y=\left[\begin{array}{cc}
I_{p} & -Z \\
0 & I_{q}
\end{array}\right]\left[\begin{array}{cc}
T_{11} & T_{12} \\
0 & T_{22}
\end{array}\right]\left[\begin{array}{cc}
I_{p} & Z \\
0 & I_{q}
\end{array}\right]=\left[\begin{array}{cc}
T_{11} & T_{11} Z-Z T_{22}+T_{12} \\
0 & T_{22}
\end{array}\right]
$$

has the required block diagonal form.
By repcatedly applying this lemma, we can establish the following more general result.

Theorem 7.1.6 (Block Diagonal Decomposition). Suppose

$$
Q^{H} A Q=T=\left[\begin{array}{cccc}
T_{11} & T_{12} & \cdots & T_{1 q} \\
0 & T_{22} & \cdots & T_{2 q} \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & T_{q q}
\end{array}\right]
$$

is a Schur decomposition of $A \in \mathbb{C}^{n \times n}$ and that the $T_{i i}$ are square. If $\lambda\left(T_{i i}\right) \cap \lambda\left(T_{j j}\right)=\emptyset$ whenever $i \neq j$, then there exists a nonsingular matrix $Y \in \mathbb{C}^{n \times n}$ such that

$$
(Q Y)^{-1} A(Q Y)=\operatorname{diag}\left(T_{11}, \ldots, T_{q q}\right)
$$

Proof. See P7.1.2.
If each diagonal block $T_{i i}$ is associated with a distinct eigenvalue, then we obtain
Corollary 7.1.7. If $A \in \mathbb{C}^{n \times n}$, then there exists a nonsingular $X$ such that

$$
X^{-1} A X=\operatorname{diag}\left(\lambda_{1} I+N_{1}, \ldots, \lambda_{q} I+N_{q}\right) \quad N_{i} \in \mathbb{C}^{n_{i} \times n_{i}}
$$

where $\lambda_{1}, \ldots, \lambda_{q}$ are distinct, the integers $n_{1}, \ldots, n_{q}$ satisfy $n_{1}+\cdots+n_{q}=n$, and each $N_{i}$ is strictly upper triangular.

A number of important terms are connected with decomposition (7.1.11). The integer $n_{i}$ is referred to as the algebraic multiplicity of $\lambda_{i}$. If $n_{i}=1$, then $\lambda_{i}$ is said to be simple. The geometric multiplicity of $\lambda_{i}$ equals the dimensions of null $\left(N_{i}\right)$, i.e., the number of linearly independent eigenvectors associated with $\lambda_{i}$. If the algebraic multiplicity of $\lambda_{i}$ exceeds its geometric multiplicity, then $\lambda_{i}$ is said to be a defective eigenvalue. A matrix with a defective eigenvalue is referred to as a defective matrix. Nondefective matrices are also said to be diagonalizable.

Corollary 7.1.8 (Diagonal Form). $A \in \mathbb{C}^{n \times n}$ is nondefective if and only if there exists a nonsingular $X \in \mathbb{C}^{n \times n}$ such that

$$
X^{-1} A X=\operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{n}\right)
$$

Proof. $A$ is nondefective if and only if there exist independent vectors $x_{1} \ldots x_{n} \in \mathbb{C}^{n}$ and scalars $\lambda_{1}, \ldots, \lambda_{n}$ such that $A x_{i}=\lambda_{i} x_{i}$ for $i=1: n$. This is equivalent to the existence of a nonsingular $X=\left[x_{1}|\cdots| x_{n}\right] \in \mathbb{C}^{n \times n}$ such that $A X=X D$ where $D=\operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{n}\right)$.

Note that if $y_{i}^{H}$ is the $i$ th row of $X^{-1}$, then $y_{i}^{H} A=\lambda_{i} y_{i}^{H}$. Thus, the columns of $X^{-H}$ are left eigenvectors and the columns of $X$ are right eigenvectors.

If we partition the matrix $X$ in (7.1.11),

$$
X=\left[\underset{n_{1}}{X_{1}|\cdots|} \underset{n_{q}}{X_{q}}\right]
$$

then $\mathbb{C}^{n}=\operatorname{ran}\left(X_{1}\right) \oplus \ldots \oplus \operatorname{ran}\left(X_{q}\right)$, a direct sum of invariant subspaces. If the bases for these subspaces are chosen in a special way, then it is possible to introduce even more zeroes into the upper triangular portion of $X^{-1} A X$.

Theorem 7.1.9 (Jordan Decomposition). If $A \in \mathbb{C}^{n \times n}$, then there exists a nonsingular $X \in \mathbb{C}^{n \times n}$ such that $X^{-1} A X=\operatorname{diag}\left(J_{1}, \ldots, J_{q}\right)$ where

$$
J_{i}=\left[\begin{array}{ccccc}
\lambda_{i} & 1 & & \cdots & 0 \\
0 & \lambda_{i} & \ddots & & \vdots \\
& \ddots & \ddots & \ddots & \\
\vdots & & \ddots & \ddots & 1 \\
0 & \cdots & & 0 & \lambda_{i}
\end{array}\right] \in \mathbb{C}^{n_{i} \times n_{i}}
$$

and $n_{1}+\cdots+n_{q}=n$.
Proof. See Horn and Johnson (MA, p. 330)
The $J_{i}$ are referred to as Jordan blocks. The number and dimensions of the Jordan blocks associated with each distinct eigenvalue are unique, although their ordering along the diagonal is not.

### 7.1.5 Some Comments on Nonunitary Similarity

The Jordan block structure of a defective matrix is difficult to determine numerically. The set of $n$-by- $n$ diagonalizable matrices is dense in $\mathbb{C}^{n \times n}$, and thus, small changes in a defective matrix can radically alter its Jordan form. We have more to say about this in §7.6.5.

A related difficulty that arises in the eigenvalue problem is that a nearly defective matrix can have a poorly conditioned matrix of eigenvectors. For example, any matrix $X$ that diagonalizes

$$
A=\left[\begin{array}{cc}
1+\epsilon & 1 \\
0 & 1-\epsilon
\end{array}\right], \quad 0<\epsilon \ll 1,
$$

has a 2 -norm condition of order $1 / \epsilon$.
These observations serve to highlight the difficulties associated with ill-conditioned similarity transformations. Since

$$
\mathrm{fl}\left(X^{-1} A X\right)=X^{-1} A X+E
$$

where

$$
\|E\|_{2} \approx \mathbf{u} \cdot \kappa_{2}(X)\|A\|_{2},
$$

it is clear that large errors can be introduced into an eigenvalue calculation when we depart from unitary similarity.

### 7.1.6 Singular Values and Eigenvalues

Since the singular values of $A$ and its Schur decomposition $Q^{H} A Q=\operatorname{diag}\left(\lambda_{i}\right)+N$ are the same, it follows that

$$
\sigma_{\min }(A) \leq \min _{1 \leq i \leq n}\left|\lambda_{i}\right| \leq \max _{1 \leq i \leq n}\left|\lambda_{i}\right| \leq \sigma_{\max }(A) .
$$

From what we know about the condition of triangular matrices, it may be the case that

$$
\max _{1 \leq i, j \leq n} \frac{\left|\lambda_{i}\right|}{\left|\lambda_{j}\right|} \ll \kappa_{2}(A) .
$$

See §5.4.3. This is a reminder that for nonnormal matrices, eigenvalues do not have the "predictive power" of singular values when it comes to $A x=b$ sensitivity matters. Eigenvalues of nonnormal matrices have other shortcomings, a topic that is the focus of §7.9.

## Problems

P7.1.1 (a) Show that if $T \in \mathbb{C}^{n \times n}$ is upper triangular and normal, then $T$ is diagonal. (b) Show that if $A$ is normal and $Q^{H} A Q=T$ is a Schur decomposition, then $T$ is diagonal. (c) Use (a) and (b) to complete the proof of Corollary 7.1.4.
P7.1.2 Prove Theorem 7.1.6 by using induction and Lemma 7.1.5.
P7.1.3 Suppose $A \in \mathbb{C}^{n \times n}$ has distinct eigenvalues. Show that if $Q^{H} A Q=T$ is its Schur decomposition and $A B=B A$, then $Q^{H} B Q$ is upper triangular.
P7.1.4 Show that if $A$ and $B^{H}$ are in $\mathbb{C}^{m \times n}$ with $m \geq n$, then

$$
\lambda(A B)=\lambda(B A) \cup\{\underbrace{0, \ldots, 0}_{m-n}\} .
$$

P7.1.5 Given $A \in \mathbb{C}^{n \times n}$, use the Schur decomposition to show that for every $\epsilon>0$, there exists a diagonalizable matrix $B$ such that $\|A-B\|_{2} \leq \epsilon$. This shows that the set of diagonalizable matrices is dense in $\mathbb{C}^{n \times n}$ and that the Jordan decomposition is not a continuous matrix decomposition.
P7.1.6 Suppose $A_{k} \rightarrow A$ and that $Q_{k}^{H} A_{k} Q_{k}=T_{k}$ is a Schur decomposition of $A_{k}$. Show that $\left\{Q_{k}\right\}$ has a converging subsequence $\left\{Q_{k_{i}}\right\}$ with the property that

$$
\lim _{i \rightarrow \infty} Q_{k_{i}}=Q
$$

where $Q^{H} A Q=T$ is upper triangular. This shows that the eigenvalues of a matrix are continuous functions of its entries.
P7.1.7 Justify (7.1.14) and (7.1.15).
P7.1.8 Show how to compute the eigenvalues of

$$
M=\underset{k}{\left[\begin{array}{cc}
A & C \\
B & D \\
k & j
\end{array}\right]_{j}^{k}}
$$

where $A, B, C$, and $D$ are given real diagonal matrices.
P7.1.9 Use the Jordan decomposition to show that if all the eigenvalues of a matrix $A$ are strictly less than unity, then $\lim _{k \rightarrow \infty} A^{k}=0$.
P7.1.10 The initial value problem

$$
\begin{array}{ll}
\dot{x}(t)=y(t), & x(0)=1 \\
\dot{y}(t)=-x(t), & y(0)=0
\end{array}
$$

has solution $x(t)=\cos (t)$ and $y(t)=\sin (t)$. Let $h>0$. Here are three reasonable iterations that can be used to compute approximations $x_{k} \approx x(k h)$ and $y_{k} \approx y(k h)$ assuming that $x_{0}=1$ and $y_{k}=0$ :

Method 1:

$$
\begin{aligned}
& x_{k+1}=x_{k}+h y_{k}, \\
& y_{k+1}=y_{k}-h x_{k},
\end{aligned}
$$

Method 2:

$$
\begin{aligned}
& x_{k+1}=x_{k}+h y_{k}, \\
& y_{k+1}=y_{k}-h x_{k+1},
\end{aligned}
$$

Method 3:

$$
\begin{aligned}
& x_{k+1}=x_{k}+h y_{k+1}, \\
& y_{k+1}=y_{k}-h x_{k+1} .
\end{aligned}
$$

Express each method in the form

$$
\left[\begin{array}{l}
x_{k+1} \\
y_{k+1}
\end{array}\right]=A_{h}\left[\begin{array}{l}
x_{k} \\
y_{k}
\end{array}\right]
$$

where $A_{h}$ is a 2-by-2 matrix. For each case, compute $\lambda\left(A_{h}\right)$ and use the previous problem to discuss $\lim x_{k}$ and $\lim y_{k}$ as $k \rightarrow \infty$.
P7.1.11 If $J \in \mathbf{R}^{d \times d}$ is a Jordan block, what is $\kappa_{\infty}(J)$ ?
P7.1.12 Suppose $A, B \in \mathbb{C}^{n \times n}$. Show that the $2 n$-by- $2 n$ matrices

$$
M_{1}=\left[\begin{array}{cc}
A B & 0 \\
B & 0
\end{array}\right] \quad \text { and } \quad M_{2}=\left[\begin{array}{cc}
0 & 0 \\
B & B A
\end{array}\right]
$$

are similar thereby showing that $\lambda(A B)=\lambda(B A)$.
P7.1.13 Suppose $A \in \mathbf{R}^{n \times n}$. We say that $B \in \mathbf{R}^{n \times n}$ is the Drazin inverse of $A$ if (i) $A B=B A$, (ii) $B A B=B$, and (iii) the spectral radius of $A-A B A$ is zero. Give a formula for $B$ in terms of the Jordan decomposition of $A$ paying particular attention to the blocks associated with $A$ 's zero eigenvalues.
P7.1.14 Show that if $A \in \mathbf{R}^{n \times n}$, then $\rho(A) \geq\left(\sigma_{1} \cdots \sigma_{n}\right)^{1 / n}$ where $\sigma_{1}, \ldots, \sigma_{n}$ are the singular values of $A$.
P7.1.15 Consider the polynomial $q(x)=\operatorname{det}\left(I_{n}+x A\right)$ where $A \in \mathbf{R}^{n \times n}$. We wish to compute the coefficient of $x^{2}$. (a) Specify the coefficient in terms of the eigenvalues $\lambda_{1}, \ldots, \lambda_{n}$ of $A$. (b) Give a simple formula for the coefficient in terms of $\operatorname{tr}(A)$ and $\operatorname{tr}\left(A^{2}\right)$.
P7.1.16 Given $A \in \mathbf{R}^{2 \times 2}$, show that there exists a nonsingular $X \in \mathbf{R}^{2 \times 2}$ so $X^{-1} A X=A^{T}$. See Dubrulle and Parlett (2007).

## Notes and References for §7.1

For additional discussion about the linear algebra behind the eigenvalue problem, see Horn and Johnson (MA) and:
L. Mirsky (1963). An Introduction to Linear Algebra, Oxford University Press, Oxford, U.K.
M. Marcus and H. Minc (1964). A Survey of Matrix Theory and Matrix Inequalities, Allyn and Bacon, Boston.
R. Bellman (1970). Introduction to Matrix Analysis, second edition, McGraw-Hill, New York.
I. Gohberg, P. Lancaster, and L. Rodman (2006). Invariant Subspaces of Matrices with Applications, SIAM Publications, Philadelphia, PA.

For a general discussion about the similarity connection between a matrix and its transpose, see:
A.A. Dubrulle and B.N. Parlett (2010). "Revelations of a Transposition Matrix," J. Comp. and Appl. Math. 233, 1217-1219.

The Schur decomposition originally appeared in:
I. Schur (1909). "On the Characteristic Roots of a Linear Substitution with an Application to the Theory of Integral Equations." Math. Ann. 66, 488-510 (German).

A proof very similar to ours is given in:
H.W. Turnbull and A.C. Aitken (1961). An Introduction to the Theory of Canonical Forms, Dover, New York, 105.

### 7.2 Perturbation Theory

The act of computing eigenvalues is the act of computing zeros of the characteristic polynomial. Galois theory tells us that such a process has to be iterative if $n>4$ and so errors arise because of finite termination. In order to develop intelligent stopping criteria we need an informative perturbation theory that tells us how to think about approximate eigenvalues and invariant subspaces.

### 7.2.1 Eigenvalue Sensitivity

An important framework for eigenvalue computation is to produce a sequence of similarity transformations $\left\{X_{k}\right\}$ with the property that the matrices $X_{k}^{-1} A X_{k}$ are progressively "more diagonal." The question naturally arises, how well do the diagonal elements of a matrix approximate its eigenvalues?

Theorem 7.2.1 (Gershgorin Circle Theorem). If $X^{-1} A X=D+F$ where $D= \operatorname{diag}\left(d_{1}, \ldots, d_{n}\right)$ and $F$ has zero diagonal entries, then

$$
\lambda(A) \subseteq \bigcup_{i=1}^{n} D_{i}
$$

where $D_{i}=\left\{z \in \mathbb{C}:\left|z-d_{i}\right| \leq \sum_{j=1}^{n}\left|f_{i j}\right|\right\}$.
Proof. Suppose $\lambda \in \lambda(A)$ and assume without loss of generality that $\lambda \neq d_{i}$ for $i=1: n$. Since $(D-\lambda I)+F$ is singular, it follows from Lemma 2.3.3 that

$$
1 \leq\left\|(D-\lambda I)^{-1} F\right\|_{\infty}=\sum_{j=1}^{n} \frac{\left|f_{k j}\right|}{\left|d_{k}-\lambda\right|}
$$

for some $k, 1 \leq k \leq n$. But this implies that $\lambda \in D_{k}$.
It can also be shown that if the Gershgorin disk $D_{i}$ is isolated from the other disks, then it contains precisely one eigenvalue of $A$. Scc Wilkinson (AEP, pp. 71ff.).

For some methods it is possible to show that the computed eigenvalues are the exact eigenvalues of a matrix $A+E$ where $E$ is small in norm. Consequently, we should understand how the eigenvalues of a matrix can be affected by small perturbations.

Theorem 7.2.2 (Bauer-Fike). If $\mu$ is an eigenvalue of $A+E \in \mathbb{C}^{n \times n}$ and $X^{-1} A X= D=\operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{n}\right)$, then

$$
\min _{\lambda \in \lambda(A)}|\lambda-\mu| \leq \kappa_{p}(X)\|E\|_{p}
$$

where $\|\cdot\|_{p}$ denotes any of the $p$-norms.
Proof. If $\mu \in \lambda(A)$, then the theorem is obviously true. Otherwise if the matrix $X^{-1}(A+E-\mu I) X$ is singular, then so is $I+(D-\mu I)^{-1}\left(X^{-1} E X\right)$. Thus, from

Lemma 2.3.3 we obtain

$$
1 \leq\left\|(D-\mu I)^{-1}\left(X^{-1} E X\right)\right\|_{p} \leq\left\|(D-\mu I)^{-1}\right\|_{p}\|X\|_{p}\|E\|_{p}\left\|X^{-1}\right\|_{p}
$$

Since $(D-\mu I)^{-1}$ is diagonal and the $p$-norm of a diagonal matrix is the absolute value of the largest diagonal entry, it follows that

$$
\left\|(D-\mu I)^{-1}\right\|_{p}=\max _{\lambda \in \lambda(A)} \frac{1}{|\lambda-\mu|}
$$

completing the proof.
An analogous result can be obtained via the Schur decomposition:
Theorem 7.2.3. Let $Q^{H} A Q=D+N$ be a Schur decomposition of $A \in \mathbb{C}^{n \times n}$ as in (7.1.7). If $\mu \in \lambda(A+E)$ and $p$ is the smallest positive integer such that $|N|^{p}=0$, then

$$
\min _{\lambda \in \lambda(A)}|\lambda-\mu| \leq \max \left\{\theta, \theta^{1 / p}\right\}
$$

where

$$
\theta=\|E\|_{2} \sum_{k=0}^{p-1}\|N\|_{2}^{k}
$$

Proof. Define

$$
\delta=\min _{\lambda \in \lambda(A)}|\lambda-\mu|=\frac{1}{\left\|(\mu I-D)^{-1}\right\|_{2}}
$$

The theorem is clearly true if $\delta=0$. If $\delta>0$, then $I-(\mu I-A)^{-1} E$ is singular and by Lemma 2.3.3 we have

$$
\begin{aligned}
1 & \leq\left\|(\mu I-A)^{-1} E\right\|_{2} \leq\left\|(\mu I-A)^{-1}\right\|_{2}\|E\|_{2} \\
& =\left\|((\mu I-D)-N)^{-1}\right\|_{2}\|E\|_{2}
\end{aligned}
$$

Since $(\mu I-D)^{-1}$ is diagonal and $|N|^{p}=0$, it follows that $\left((\mu I-D)^{-1} N\right)^{p}=0$. Thus,

$$
((\mu I-D)-N)^{-1}=\sum_{k=0}^{p-1}\left((\mu I-D)^{-1} N\right)^{k}(\mu I-D)^{-1}
$$

and so

$$
\left\|((\mu I-D)-N)^{-1}\right\|_{2} \leq \frac{1}{\delta} \sum_{k=0}^{p-1}\left(\frac{\|N\|_{2}}{\delta}\right)^{k}
$$

If $\delta>1$, then

$$
\|(\mu I-D)-N)^{-1}\left\|_{2} \leq \frac{1}{\delta} \sum_{k=0}^{p-1}\right\| N \|_{2}^{k}
$$

and so from (7.2.1), $\delta \leq \theta$. If $\delta \leq 1$, then

$$
\|(\mu I-D)-N)^{-1}\left\|_{2} \leq \frac{1}{\delta^{p}} \sum_{k=0}^{p-1}\right\| N \|_{2}^{k}
$$

By using (7.2.1) again we have $\delta^{p} \leq \theta$ and so $\delta \leq \max \left\{\theta, \theta^{1 / p}\right\}$.
Theorems 7.2.2 and 7.2.3 suggest that the eigenvalues of a nonnormal matrix may be sensitive to perturbations. In particular, if $\kappa_{2}(X)$ or $\|N\|_{2}^{p-1}$ is large, then small changes in $A$ can induce large changes in the eigenvalues.

### 7.2.2 The Condition of a Simple Eigenvalue

Extreme eigenvalue sensitivity for a matrix $A$ cannot occur if $A$ is normal. On the other hand, nonnormality does not necessarily imply eigenvalue sensitivity. Indeed, a nonnormal matrix can have a mixture of well-conditioned and ill-conditioned eigenvalues. For this reason, it is beneficial to refine our perturbation theory so that it is applicable to individual eigenvalues and not the spectrum as a whole.

To this end, suppose that $\lambda$ is a simple eigenvalue of $A \in \mathbb{C}^{n \times n}$ and that $x$ and $y$ satisfy $A x=\lambda x$ and $y^{H} A=\lambda y^{H}$ with $\|x\|_{2}=\|y\|_{2}=1$. If $Y^{H} A X=J$ is the Jordan decomposition with $Y^{H}=X^{-1}$, then $y$ and $x$ are nonzero multiples of $X(:, i)$ and $Y(:, i)$ for some $i$. It follows from $1=Y(:, i)^{H} X(:, i)$ that $y^{H} x \neq 0$, a fact that we shall use shortly.

Using classical results from function theory, it can be shown that in a neighborhood of the origin there exist differentiable $x(\epsilon)$ and $\lambda(\epsilon)$ such that

$$
(A+\epsilon F) x(\epsilon)=\lambda(\epsilon) x(\epsilon), \quad\|F\|_{2}=1
$$

where $\lambda(0)=\lambda$ and $x(0)=x$. By differentiating this equation with respect to $\epsilon$ and setting $\epsilon=0$ in the result, we obtain

$$
A \dot{x}(0)+F x=\dot{\lambda}(0) x+\lambda \dot{x}(0) .
$$

Applying $y^{H}$ to both sides of this equation, dividing by $y^{H} x$, and taking absolute values gives

$$
|\dot{\lambda}(0)|=\left|\frac{y^{H} F x}{y^{H} x}\right| \leq \frac{1}{\left|y^{H} x\right|}
$$

The upper bound is attained if $F=y x^{H}$. For this reason we refer to the reciprocal of

$$
s(\lambda)=\left|y^{H} x\right|
$$

as the condition of the eigenvalue $\lambda$.
Roughly speaking, the above analysis shows that $O(\epsilon)$ perturbations in $A$ can induce $\epsilon / s(\lambda)$ changes in an eigenvalue. Thus, if $s(\lambda)$ is small, then $\lambda$ is appropriately regarded as ill-conditioned. Note that $s(\lambda)$ is the cosine of the angle between the left and right eigenvectors associated with $\lambda$ and is unique only if $\lambda$ is simple.

A small $s(\lambda)$ implies that $A$ is near a matrix having a multiple eigenvalue. In particular, if $\lambda$ is distinct and $s(\lambda)<1$, then there exists an $E$ such that $\lambda$ is a repeated eigenvalue of $A+E$ and

$$
\frac{\|E\|_{2}}{\|A\|_{2}} \leq \frac{s(\lambda)}{\sqrt{1-s(\lambda)^{2}}}
$$

This result is proved by Wilkinson (1972).

### 7.2.3 Sensitivity of Repeated Eigenvalues

If $\lambda$ is a repeated eigenvalue, then the eigenvalue sensitivity question is more complicated. For example, if

$$
A=\left[\begin{array}{ll}
1 & a \\
0 & 1
\end{array}\right] \quad \text { and } \quad F=\left[\begin{array}{ll}
0 & 0 \\
1 & 0
\end{array}\right]
$$

then $\lambda(A+\epsilon F)=\{1 \pm \sqrt{\epsilon a}\}$. Note that if $a \neq 0$, then it follows that the eigenvalues of $A+\epsilon F$ are not differentiable at zero; their rate of change at the origin is infinite. In general, if $\lambda$ is a defective eigenvalue of $A$, then $O(\epsilon)$ perturbations in $A$ can result in $O\left(\epsilon^{1 / p}\right)$ perturbations in $\lambda$ if $\lambda$ is associated with a $p$-dimensional Jordan block. See Wilkinson (AEP, pp. 77ff.) for a more detailed discussion.

### 7.2.4 Invariant Subspace Sensitivity

A collection of sensitive eigenvectors can define an insensitive invariant subspace provided the corresponding cluster of eigenvalues is isolated. To be precise, suppose

$$
Q^{H} A Q=\left[\begin{array}{cc}
T_{11} & T_{12} \\
0 & T_{22} \\
r & n-r
\end{array}\right]_{n-r}^{r}
$$

is a Schur decomposition of $A$ with

$$
Q=\left[\underset{r}{Q_{1}} \mid \underset{n-r}{Q_{2}}\right]
$$

It is clear from our discussion of eigenvector perturbation that the sensitivity of the invariant subspace ran $\left(Q_{1}\right)$ depends on the distance between $\lambda\left(T_{11}\right)$ and $\lambda\left(T_{22}\right)$. The proper measure of this distance turns out to be the smallest singular value of the linear transformation $X \rightarrow T_{11} X-X T_{22}$. (Recall that this transformation figures in Lemma 7.1.5.) In particular, if we define the separation between the matrices $T_{11}$ and $T_{22}$ by

$$
\operatorname{sep}\left(T_{11}, T_{22}\right)=\min _{X \neq 0} \frac{\left\|T_{11} X-X T_{22}\right\|_{F}}{\|X\|_{F}},
$$

then we have the following general result:

Theorem 7.2.4. Suppose that (7.2.3) and (7.2.4) hold and that for any matrix $E \in \mathbb{C}^{n \times n}$ we partition $Q^{H} E Q$ as follows:

$$
Q^{H} E Q=\left[\begin{array}{cc}
E_{11} & E_{12} \\
E_{21} & E_{22} \\
r & n-r
\end{array}\right]_{n-r}^{r}
$$

If $\operatorname{sep}\left(T_{11}, T_{22}\right)>0$ and

$$
\|E\|_{F}\left(1+\frac{5\left\|T_{12}\right\|_{F}}{\operatorname{sep}\left(T_{11}, T_{22}\right)}\right) \leq \frac{\operatorname{sep}\left(T_{11}, T_{22}\right)}{5}
$$

then there exists a $P \in \mathbb{C}^{(n-r) \times r}$ with

$$
\|P\|_{F} \leq 4 \frac{\left\|E_{21}\right\|_{F}}{\operatorname{sep}\left(T_{11}, T_{22}\right)}
$$

such that the columns of $\widetilde{Q}_{1}=\left(Q_{1}+Q_{2} P\right)\left(I+P^{H} P\right)^{-1 / 2}$ are an orthonormal basis for a subspace invariant for $A+E$.

Proof. This result is a slight recasting of Theorem 4.11 in Stewart (1973) which should be consulted for proof details. See also Stewart and Sun (MPA, p. 230). The matrix $\left(I+P^{H} P\right)^{-1 / 2}$ is the inverse of the square root of the symmetric positive definite matrix $I+P^{H} P$. Scc §4.2.4.

Corollary 7.2.5. If the assumptions in Theorem 7.2.4 hold, then

$$
\operatorname{dist}\left(\operatorname{ran}\left(Q_{1}\right), \operatorname{ran}\left(\tilde{Q}_{1}\right)\right) \leq 4 \frac{\left\|E_{21}\right\|_{F}}{\operatorname{sep}\left(T_{11}, T_{22}\right)}
$$

Proof. Using the SVD of $P$, it can be shown that

$$
\left\|P\left(I+P^{H} P\right)^{-1 / 2}\right\|_{2} \leq\|P\|_{2} \leq\|P\|_{F}
$$

Since the required distance is the 2 -norm of $Q_{2}^{H} \widetilde{Q}_{1}=P\left(I+P^{H} P\right)^{-1 / 2}$, the proof is complete.

Thus, the reciprocal of $\operatorname{sep}\left(T_{11}, T_{22}\right)$ can be thought of as a condition number that measures the sensitivity of $\operatorname{ran}\left(Q_{1}\right)$ as an invariant subspace.

### 7.2.5 Eigenvector Sensitivity

If we set $r=1$ in the preceding subsection, then the analysis addresses the issue of eigenvector sensitivity.

Corollary 7.2.6. Suppose $A, E \in \mathbb{C}^{n \times n}$ and that $Q=\left[q_{1} \mid Q_{2}\right] \in \mathbb{C}^{n \times n}$ is unitary with $q_{1} \in \mathbb{C}^{n}$. Assume

$$
Q^{H} A Q=\left[\begin{array}{cc}
\lambda & v^{H} \\
0 & T_{22} \\
1 & n-1
\end{array}\right]_{n-1}^{1}, \quad Q^{H} E Q=\left[\begin{array}{cc}
\epsilon & \gamma^{H} \\
\delta & E_{22} \\
1 & n-1
\end{array}\right]_{n-1}^{1}
$$

(Thus, $q_{1}$ is an eigenvector.) If $\sigma=\sigma_{\min }\left(T_{22}-\lambda I\right)>0$ and

$$
\|E\|_{F}\left(1+\frac{5\|v\|_{2}}{\sigma}\right) \leq \frac{\sigma}{5}
$$

then there exists $p \in \mathbb{C}^{n-1}$ with

$$
\|p\|_{2} \leq 4 \frac{\|\delta\|_{2}}{\sigma}
$$

such that $\tilde{q}_{1}=\left(q_{1}+Q_{2} p\right) / \sqrt{1+p^{H} p}$ is a unit 2-norm eigenvector for $A+E$. Moreover,

$$
\operatorname{dist}\left(\operatorname{span}\left\{q_{1}\right\}, \operatorname{span}\left\{\tilde{q}_{1}\right\}\right) \leq 4 \frac{\|\delta\|_{2}}{\sigma} .
$$

Proof. The result follows from Theorem 7.2.4, Corollary 7.2.5, and the observation that if $T_{11}=\lambda$, then $\operatorname{sep}\left(T_{11}, T_{22}\right)=\sigma_{\text {min }}\left(T_{22}-\lambda I\right)$.

Note that $\sigma_{\min }\left(T_{22}-\lambda I\right)$ roughly measures the separation of $\lambda$ from the eigenvalues of $T_{22}$. We have to say "roughly" because

$$
\operatorname{sep}\left(\lambda, T_{22}\right)=\sigma_{\min }\left(T_{22}-\lambda I\right) \leq \min _{\mu \in \lambda\left(T_{22}\right)}|\mu-\lambda|
$$

and the upper bound can be a gross overestimate.
That the separation of the eigenvalues should have a bearing upon eigenvector sensitivity should come as no surprise. Indeed, if $\lambda$ is a nondefective, repeated eigenvalue, then there are an infinite number of possible eigenvector bases for the associated invariant subspace. The preceding analysis merely indicates that this indeterminancy begins to be felt as the eigenvalues coalesce. In other words, the eigenvectors associated with nearby eigenvalues are "wobbly."

## Problems

P7.2.1 Suppose $Q^{H} A Q=\operatorname{diag}\left(\lambda_{1}\right)+N$ is a Schur decomposition of $A \in \mathbb{C}^{n \times n}$ and define $\nu(A)= \left\|A^{H} A-A A^{H}\right\|_{F}$. The upper and lower bounds in

$$
\frac{\nu(A)^{2}}{6\|A\|_{F}^{2}} \leq\|N\|_{F}^{2} \leq \sqrt{\frac{n^{3}-n}{12}} \nu(A)
$$

are established by Henrici (1962) and Eberlein (1965), respectively. Verify these results for the case $n=2$.
P7.2.2 Suppose $A \in \mathbb{C}^{n \times n}$ and $X^{-1} A X=\operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{n}\right)$ with distinct $\lambda_{i}$. Show that if the columns of $X$ have unit 2 -norm, then $\kappa_{F}(X)^{2}=n\left(1 / s\left(\lambda_{1}\right)^{2}+\cdots+1 / s\left(\lambda_{n}\right)^{2}\right.$ ).
P7.2.3 Suppose $Q^{H} A Q=\operatorname{diag}\left(\lambda_{i}\right)+N$ is a Schur decomposition of $A$ and that $X^{-1} A X=\operatorname{diag}\left(\lambda_{i}\right)$. Show $\kappa_{2}(X)^{2} \geq 1+\left(\|N\|_{F} /\|A\|_{F}\right)^{2}$. See Loizou (1969).
P7.2.4 If $X^{-1} A X=\operatorname{diag}\left(\lambda_{i}\right)$ and $\left|\lambda_{1}\right| \geq \cdots \geq\left|\lambda_{n}\right|$, then

$$
\frac{\sigma_{i}(A)}{\kappa_{2}(X)} \leq\left|\lambda_{i}\right| \leq \kappa_{2}(X) \sigma_{i}(A)
$$

Prove this result for the $n=2$ case. See Ruhe (1975).
P7.2.5 Show that if $A=\left[\begin{array}{ll}a & c \\ 0 & b\end{array}\right]$ and $a \neq b$, then $s(a)=s(b)=\left(1+|c /(a-b)|^{2}\right)^{-1 / 2}$.

P7.2.6 Suppose

$$
A=\left[\begin{array}{cc}
\lambda & v^{T} \\
0 & T_{22}
\end{array}\right]
$$

and that $\lambda \notin \lambda\left(T_{22}\right)$. Show that if $\sigma=\operatorname{sep}\left(\lambda, T_{22}\right)$, then

$$
s(\lambda)=\frac{1}{\sqrt{1+\left\|\left(T_{22}-\lambda I\right)^{-1} v\right\|_{2}^{2}}} \leq \frac{\sigma}{\sqrt{\sigma^{2}+\|v\|_{2}^{2}}}
$$

where $s(\lambda)$ is defined in (7.2.2).
P7.2.7 Show that the condition of a simple eigenvalue is preserved under unitary similarity transformations.
P7.2.8 With the same hypothesis as in the Bauer-Fike theorem (Theorem 7.2.2), show that

$$
\min _{\lambda \in \lambda(A)}|\lambda-\mu| \leq\left\|\left|X^{-1}\right||E||X|\right\|_{p} .
$$

P7.2.9 Verify (7.2.6).
P7.2.10 Show that if $B \in \mathbb{C}^{m \times m}$ and $C \in \mathbb{C}^{n \times n}$, then $\operatorname{sep}(B, C)$ is less than or equal to $|\lambda-\mu|$ for all $\lambda \in \lambda(B)$ and $\mu \in \lambda(C)$.

## Notes and References for §7.2

Many of the results presented in this section may be found in Wilkinson (AEP), Stewart and Sun (MPA) as well as:
F.L. Bauer and C.T. Fike (1960). "Norms and Exclusion Theorems," Numer. Math. 2, 123-44.
A.S. Householder (1964). The Theory of Matrices in Numerical Analysis. Blaisdell, New York.
R. Bhatia (2007). Perturbation Bounds for Matrix Eigenvalues, SIAM Publications, Philadelphia, PA.

Early papers concerned with the effect of perturbations on the eigenvalues of a general matrix include:
A. Ruhe (1970). "Perturbation Bounds for Means of Eigenvalues and Invariant Subspaces," BIT 10, 343-54.
A. Ruhe (1970). "Properties of a Matrix with a Very Ill-Conditioned Eigenproblem," Numer. Math. 15, 57-60.
J.H. Wilkinson (1972). "Note on Matrices with a Very Ill-Conditioned Eigenproblem," Numer. Math. 19, 176-78.
W. Kahan, B.N. Parlett, and E. Jiang (1982). "Residual Bounds on Approximate Eigensystems of Nonnormal Matrices," SIAM J. Numer. Anal. 19, 470-484.
J.H. Wilkinson (1984). "On Neighboring Matrices with Quadratic Elementary Divisors," Numer. Math. 44, 1-21.

Wilkinson's work on nearest defective matrices is typical of a growing body of literature that is concerned with "nearness" problems, see:
A. Ruhe (1987). "Closest Normal Matrix Found!," BIT 27, 585-598.
J.W. Demmel (1987). "On the Distance to the Nearest Ill-Posed Problem," Numer. Math. 51, 251-289.
J.W. Demmel (1988). "The Probability that a Numerical Analysis Problem is Difficult," Math. Comput. 50, 449-480.
N.J. Higham (1989). "Matrix Nearness Problems and Applications," in Applications of Matrix Theory, M.J.C. Gover and S. Barnett (eds.), Oxford University Press, Oxford, 1-27.
A.N. Malyshev (1999). "A Formula for the 2-norm Distance from a Matrix to the Set of Matrices with Multiple Eigenvalues," Numer. Math. 83, 443-454.
J.-M. Gracia (2005). "Nearest Matrix with Two Prescribed Eigenvalues," Lin. Alg. Applic. 401, 277-294.

An important subset of this literature is concerned with nearness to the set of unstable matrices. A matrix is unstable if it has an eigenvalue with nonnegative real part. Controllability is a related notion, see:
C. Van Loan (1985). "How Near is a Stable Matrix to an Unstable Matrix?," Contemp. Math. 47, 465-477.
J.W. Demmel (1987). "A Counterexample for two Conjectures About Stability," IEEE Trans. Autom. Contr. AC-32, 340-342.
R. Byers (1988). "A Bisection Method for Measuring the distance of a Stable Matrix to the Unstable Matrices," J. Sci. Stat. Comput. 9, 875-881.
J.V. Burke and M.L. Overton (1992). "Stable Perturbations of Nonsymmetric Matrices," Lin. Alg. Applic. 171, 249-273.
C. He and G.A. Watson (1998). "An Algorithm for Computing the Distance to Instability," SIAM J. Matrix Anal. Applic. 20, 101-116.
M. Gu, E. Mengi, M.L. Overton, J. Xia, and J. Zhu (2006). "Fast Methods for Estimating the Distance to Uncontrollability," SIAM J. Matrix Anal. Applic. 28, 477-502.

Aspects of eigenvalue condition are discussed in:
C. Van Loan (1987). "On Estimating the Condition of Eigenvalues and Eigenvectors," Lin. Alg. Applic. 88/89, 715-732.
C.D. Meyer and G.W. Stewart (1988). "Derivatives and Perturbations of Eigenvectors," SIAM J. Numer. Anal. 25, 679-691.
G.W. Stewart and G. Zhang (1991). "Eigenvalues of Graded Matrices and the Condition Numbers of Multiple Eigenvalues," Numer. Math. 58, 703-712.
J.-G. Sun (1992). "On Condition Numbers of a Nondefective Multiple Eigenvalue," Numer. Math. 61, 265-276.
S.M. Rump (2001). "Computational Error Bounds for Multiple or Nearly Multiple Eigenvalues," Lin. Alg. Applic. 324, 209-226.

The relationship between the eigenvalue condition number, the departure from normality, and the condition of the eigenvector matrix is discussed in:
P. Henrici (1962). "Bounds for Iterates, Inverses, Spectral Variation and Fields of Values of Nonnormal Matrices," Numer. Math. 4, 24-40.
P. Eberlein (1965). "On Measures of Non-Normality for Matrices," AMS Monthly 72, 995-996.
R.A. Smith (1967). "The Condition Numbers of the Matrix Eigenvalue Problem," Numer. Math. 10 232-240.
G. Loizou (1969). "Nonnormality and Jordan Condition Numbers of Matrices," J. ACM 16, 580-640.
A. van der Sluis (1975). "Perturbations of Eigenvalues of Non-normal Matrices," Commun. ACM 18, 30-36.
S.L. Lee (1995). "A Practical Upper Bound for Departure from Normality," SIAM J. Matrix Anal. Applic. 16, 462-468.

Gershgorin's theorem can be used to derive a comprehensive perturbation theory. The theorem itself can be generalized and extended in various ways, see:
R.S. Varga (1970). "Minimal Gershgorin Sets for Partitioned Matrices," SIAM J. Numer. Anal. 7, 493-507.
R.J. Johnston (1971). "Gershgorin Theorems for Partitioned Matrices," Lin. Alg. Applic. 4, 205-20.
R.S. Varga and A. Krautstengl (1999). "On Gergorin-type Problems and Ovals of Cassini," ETNA 8, 15-20.
R.S. Varga (2001). "Gergorin-type Eigenvalue Inclusion Theorems and Their Sharpness," ETNA 12, 113-133.
C. Beattie and I.C.F. Ipsen (2003). "Inclusion Regions for Matrix Eigenvalues," Lin. Alg. Applic. 358, 281-291.

In our discussion, the perturbations to the $A$-matrix are general. More can be said when the perturbations are structured, see:
G.W. Stewart (2001). "On the Eigensystems of Graded Matrices," Numer. Math. 90, 349-370.
J. Moro and F.M. Dopico (2003). "Low Rank Perturbation of Jordan Structure," SIAM J. Matrix Anal. Applic. 25, 495-506.
R. Byers and D. Kressner (2004). "On the Condition of a Complex Eigenvalue under Real Perturbations," BIT 44, 209-214.
R. Byers and D. Kressner (2006). "Structured Condition Numbers for Invariant Subspaces," SIAM J. Matrix Anal. Applic. 28, 326-347.

An absolute perturbation bound comments on the difference between an eigenvalue $\lambda$ and its perturbation $\bar{\lambda}$. A relative perturbation bound examines the quotient $|\lambda-\bar{\lambda}| /|\lambda|$, something that can be very important when there is a concern about a small eigenvalue. For results in this direction consult:
R.-C. Li (1997). "Relative Perturbation Theory. III. More Bounds on Eigenvalue Variation," Lin. Alg. Applic. 266, 337-345.
S.C. Eisenstat and I.C.F. Ipsen (1998). "Three Absolute Perturbation Bounds for Matrix Eigenvalues Imply Relative Bounds," SIAM J. Matrix Anal. Applic. 20, 149-158.
S.C. Eisenstat and I.C.F. Ipsen (1998). "Relative Perturbation Results for Eigenvalues and Eigenvectors of Diagonalisable Matrices," BIT 38, 502-509.
I.C.F. Ipsen (1998). "Relative Perturbation Results for Matrix Eigenvalues and Singular Values," Acta Numerica, 7, 151-201.
I.C.F. Ipsen (2000). "Absolute and Relative Perturbation Bounds for Invariant Subspaces of Matrices," Lin. Alg. Applic. 309, 45-56.
I.C.F. Ipsen (2003). "A Note on Unifying Absolute and Relative Perturbation Bounds," Lin. Alg. Applic. 358, 239-253.
Y. Wei, X. Li, F. Bu, and F. Zhang (2006). "Relative Perturbation Bounds for the Eigenvalues of Diagonalizable and Singular Matrices-Application to Perturbation Theory for Simple Invariant Subspaces," Lin. Alg. Applic. 419, 765-771.

The eigenvectors and invariant subspaces of a matrix also "move" when there are perturbations.
Tracking these changes is typically more challenging than tracking changes in the eigenvalues, see:
T. Kato (1966). Perturbation Theory for Linear Operators, Springer-Verlag, New York.
C. Davis and W.M. Kahan (1970). "The Rotation of Eigenvectors by a Perturbation, III," SIAM J. Numer. Anal. 7, 1-46.
G.W. Stewart (1971). "Error Bounds for Approximate Invariant Subspaces of Closed Linear Operators," SIAM. J. Numer. Anal. 8, 796-808.
G.W. Stewart (1973). "Error and Perturbation Bounds for Subspaces Associated with Certain Eigenvalue Problems," SIAM Review 15, 727-764.
J. Xie (1997). "A Note on the Davis-Kahan sin(20) Theorem," Lin. Alg. Applic. 258, 129-135.
S.M. Rump and J.-P.M. Zemke (2003). "On Eigenvector Bounds," BIT 43, 823-837.

Detailed analyses of the function sep(.,.) and the map $X \rightarrow A X+X A^{T}$ are given in:
J. Varah (1979). "On the Separation of Two Matrices," SIAM J. Numer. Anal. 16, 216-22.
R. Byers and S.G. Nash (1987). "On the Singular Vectors of the Lyapunov Operator," SIAM J. Alg. Disc. Methods 8, 59-66.

### 7.3 Power Iterations

Suppose that we are given $A \in \mathbb{C}^{n \times n}$ and a unitary $U_{0} \in \mathbb{C}^{n \times n}$. Recall from §5.2.10 that the Householder QR factorization can be extended to complex matrices and consider the following iteration:

$$
\begin{aligned}
& T_{0}=U_{0}^{H} A U_{0} \\
& \text { for } k=1,2, \ldots \\
& \qquad T_{k-1}=U_{k} R_{k} \quad \text { (QR factorization) } \\
& \quad T_{k}=R_{k} U_{k} \\
& \text { end }
\end{aligned}
$$

Since $T_{k}=R_{k} U_{k}=U_{k}^{H}\left(U_{k} R_{k}\right) U_{k}=U_{k}^{H} T_{k-1} U_{k}$ it follows by induction that

$$
T_{k}=\left(U_{0} U_{1} \cdots U_{k}\right)^{H} A\left(U_{0} U_{1} \cdots U_{k}\right) .
$$

Thus, each $T_{k}$ is unitarily similar to $A$. Not so obvious, and what is a central theme of this section, is that the $T_{k}$ almost always converge to upper triangular form, i.e., (7.3.2) almost always "converges" to a Schur decomposition of $A$.

Iteration (7.3.1) is called the $Q R$ iteration, and it forms the backbone of the most effective algorithm for computing a complete Schur decomposition of a dense general matrix. In order to motivate the method and to derive its convergence properties, two other eigenvalue iterations that are important in their own right are presented first: the power method and the method of orthogonal iteration.

### 7.3.1 The Power Method

Suppose $A \in \mathbb{C}^{n \times n}$ and $X^{-1} A X=\operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{n}\right)$ with $X=\left[x_{1}|\cdots| x_{n}\right]$. Assume that

$$
\left|\lambda_{1}\right|>\left|\lambda_{2}\right| \geq \cdots \geq\left|\lambda_{n}\right|
$$

Given a unit 2-norm $q^{(0)} \in \mathbb{C}^{n}$, the power method produces a sequence of vectors $q^{(k)}$ as follows:

$$
\begin{aligned}
& \text { for } k=1,2, \ldots \\
& \qquad \begin{aligned}
z^{(k)} & =A q^{(k-1)} \\
q^{(k)} & =z^{(k)} /\left\|z^{(k)}\right\|_{2} \\
\lambda^{(k)} & =\left[q^{(k)}\right]^{H} A q^{(k)}
\end{aligned}
\end{aligned}
$$

There is nothing special about using the 2 -norm for normalization except that it imparts a greater unity on the overall discussion in this section.

Let us examine the convergence properties of the power iteration. If

$$
q^{(0)}=a_{1} x_{1}+a_{2} x_{2}+\cdots+a_{n} x_{n}
$$

and $a_{1} \neq 0$, then

$$
A^{k} q^{(0)}=a_{1} \lambda_{1}^{k}\left(x_{1}+\sum_{j=2}^{n} \frac{a_{j}}{a_{1}}\left(\frac{\lambda_{j}}{\lambda_{1}}\right)^{k} x_{j}\right) .
$$

Since $q^{(k)} \in \operatorname{span}\left\{A^{k} q^{(0)}\right\}$ we conclude that

$$
\operatorname{dist}\left(\operatorname{span}\left\{q^{(k)}\right\}, \operatorname{span}\left\{x_{1}\right\}\right)=O\left(\left|\frac{\lambda_{2}}{\lambda_{1}}\right|^{k}\right)
$$

It is also easy to verify that

$$
\left|\lambda_{1}-\lambda^{(k)}\right|=O\left(\left|\frac{\lambda_{2}}{\lambda_{1}}\right|^{k}\right)
$$

Since $\lambda_{1}$ is larger than all the other eigenvalues in modulus, it is referred to as a dominant eigenvalue. Thus, the power method converges if $\lambda_{1}$ is dominant and if $q^{(0)}$ has a component in the direction of the corresponding dominant eigenvector $x_{1}$. The behavior of the iteration without these assumptions is discussed in Wilkinson (AEP, p. 570) and Parlett and Poole (1973).

In practice, the usefulness of the power method depends upon the ratio $\left|\lambda_{2}\right| /\left|\lambda_{1}\right|$, since it dictates the rate of convergence. The danger that $q^{(0)}$ is deficient in $x_{1}$ is less worrisome because rounding errors sustained during the iteration typically ensure that subsequent iterates have a component in this direction. Moreover, it is typically the case in applications that one has a reasonably good guess as to the direction of $x_{1}$. This guards against having a pathologically small coefficient $a_{1}$ in (7.3.4).

Note that the only thing required to implement the power method is a procedure for matrix-vector products. It is not necessary to store $A$ in an $n$-by- $n$ array. For this reason, the algorithm is of interest when the dominant eigenpair for a large sparse matrix is required. We have much more to say about large sparse eigenvalue problems in Chapter 10.

Estimates for the error $\left|\lambda^{(k)}-\lambda_{1}\right|$ can be obtained by applying the perturbation theory developed in §7.2.2. Define the vector

$$
r^{(k)}=A q^{(k)}-\lambda^{(k)} q^{(k)}
$$

and observe that $\left(A+E^{(k)}\right) q^{(k)}=\lambda^{(k)} q^{(k)}$ where $E^{(k)}=-r^{(k)}\left[q^{(k)}\right]^{H}$. Thus $\lambda^{(k)}$ is an eigenvalue of $A+E^{(k)}$ and

$$
\left|\lambda^{(k)}-\lambda_{1}\right| \approx \frac{\left\|E^{(k)}\right\|_{2}}{s\left(\lambda_{1}\right)}=\frac{\left\|r^{(k)}\right\|_{2}}{s\left(\lambda_{1}\right)} .
$$

If we use the power method to generate approximate right and left dominant eigenvectors, then it is possible to obtain an estimate of $s\left(\lambda_{1}\right)$. In particular, if $w^{(k)}$ is a unit 2 -norm vector in the direction of $\left(A^{H}\right)^{k} w^{(0)}$, then we can use the approximation $s\left(\lambda_{1}\right) \approx\left|w^{(k)^{H}} q^{(k)}\right|$.

### 7.3.2 Orthogonal Iteration

A straightforward generalization of the power method can be used to compute higherdimensional invariant subspaces. Let $r$ be a chosen integer satisfying $1 \leq r \leq n$. Given $A \in \mathbb{C}^{n \times n}$ and an $n$-by- $r$ matrix $Q_{0}$ with orthonormal columns, the method of orthogonal iteration generates a sequence of matrices $\left\{Q_{k}\right\} \subseteq \mathbb{C}^{n \times r}$ and a sequence of eigenvalue estimates $\left\{\lambda_{1}^{(k)}, \ldots, \lambda_{r}^{(k)}\right\}$ as follows:

$$
\begin{aligned}
& \text { for } k=1,2, \ldots \\
& \qquad Z_{k}=A Q_{k-1} \\
& \qquad Q_{k} R_{k}=Z_{k} \quad(\text { QR factorization }) \\
& \text { end } \quad \lambda\left(Q_{k}^{H} A Q_{k}\right)=\left\{\lambda_{1}^{(k)}, \ldots, \lambda_{r}^{(k)}\right\}
\end{aligned}
$$

Note that if $r=1$, then this is just the power method (7.3.3). Moreover, the sequence $\left\{Q_{k} e_{1}\right\}$ is precisely the sequence of vectors produced by the power iteration with starting vector $q^{(0)}=Q_{0} e_{1}$.

In order to analyze the behavior of this iteration, suppose that

$$
Q^{H} A Q=T=\operatorname{diag}\left(\lambda_{i}\right)+N, \quad\left|\lambda_{1}\right| \geq\left|\lambda_{2}\right| \geq \cdots \geq\left|\lambda_{n}\right|
$$

is a Schur decomposition of $A \in \mathbb{C}^{n \times n}$. Assume that $1 \leq r<n$ and partition $Q$ and $T$ as follows:

$$
\left.Q=\underset{r}{\left[Q_{\alpha} \mid\right.} \underset{n-r}{Q_{\beta}}\right], \quad T=\left[\begin{array}{cc}
T_{11} & T_{12} \\
0 & T_{22} \\
r & n-r
\end{array}\right]_{n-r}^{r}
$$

If $\left|\lambda_{r}\right|>\left|\lambda_{r+1}\right|$, then the subspace $D_{r}(A)=\operatorname{ran}\left(Q_{\alpha}\right)$ is rcferred to as a dominant invariant subspace. It is the unique invariant subspace associated with the eigenvalues $\lambda_{1}, \ldots, \lambda_{r}$. The following theorem shows that with reasonable assumptions, the subspaces $\operatorname{ran}\left(Q_{k}\right)$ generated by (7.3.6) converge to $D_{r}(A)$ at a rate proportional to $\left|\lambda_{r+1} / \lambda_{r}\right|^{k}$.

Theorem 7.3.1. Let the Schur decomposition of $A \in \mathbb{C}^{n \times n}$ be given by (7.3.7) and (7.3.8) with $n \geq 2$. Assume that $\left|\lambda_{r}\right|>\left|\lambda_{r+1}\right|$ and that $\mu \geq 0$ satisfies

$$
(1+\mu)\left|\lambda_{r}\right|>\|N\|_{F}
$$

Suppose $Q_{0} \in \mathbb{C}^{n \times r}$ has orthonormal columns and that $d_{k}$ is defined by

$$
d_{k}=\operatorname{dist}\left(D_{r}(A), \operatorname{ran}\left(Q_{k}\right)\right), \quad k \geq 0 .
$$

If

$$
d_{0}<1
$$

then the matrices $Q_{k}$ generated by (7.3.6) satisfy

$$
d_{k} \leq(1+\mu)^{n-2} \cdot\left(1+\frac{\left\|T_{12}\right\|_{F}}{\operatorname{sep}\left(T_{11}, T_{22}\right)}\right) \cdot\left[\frac{\left|\lambda_{r+1}\right|+\frac{\|N\|_{E}}{1+\mu}}{\left|\lambda_{r}\right|-\frac{\|N\|_{F}}{1+\mu}}\right]^{k} \cdot \frac{d_{0}}{\sqrt{1-d_{0}^{2}}}
$$

Proof. The proof is given in an appendix at the end of this section.
The condition (7.3.9) ensures that the initial matrix $Q_{0}$ is not deficient in certain eigendirections. In particular, no vector in the span of $Q_{0}$ 's columns is orthogonal to $D_{r}\left(A^{H}\right)$. The theorem essentially says that if this condition holds and if $\mu$ is chosen large enough, then

$$
\operatorname{dist}\left(D_{r}(A), \operatorname{ran}\left(Q_{k}\right)\right) \approx c\left|\frac{\lambda_{r+1}}{\lambda_{r}}\right|^{k}
$$

where $c$ depends on $\operatorname{sep}\left(T_{11}, T_{22}\right)$ and $A$ 's departure from normality.
It is possible to accelerate the convergence in orthogonal iteration using a technique described in Stewart (1976). In the accelerated scheme, the approximate eigenvalue $\lambda_{i}^{(k)}$ satisfies

$$
\left|\lambda_{i}^{(k)}-\lambda_{i}\right| \approx\left|\frac{\lambda_{r+1}}{\lambda_{i}}\right|^{k}, \quad i=1: r
$$

(Without the acceleration, the right-hand side is $\left|\lambda_{i+1} / \lambda_{i}\right|^{k}$.) Stewart's algorithm involves computing the Schur decomposition of the matrices $Q_{k}^{T} A Q_{k}$ every so often. The method can be very useful in situations where $A$ is large and sparse and a few of its largest eigenvalues are required.

### 7.3.3 The QR Iteration

We now derive the QR iteration (7.3.1) and examine its convergence. Suppose $r=n$ in (7.3.6) and the eigenvalues of $A$ satisfy

$$
\left|\lambda_{1}\right|>\left|\lambda_{2}\right|>\cdots>\left|\lambda_{n}\right| .
$$

Partition the matrix $Q$ in (7.3.7) and $Q_{k}$ in (7.3.6) as follows:

$$
Q=\left[q_{1}|\cdots| q_{n}\right], \quad Q_{k}=\left[q_{1}^{(k)}|\cdots| q_{n}^{(k)}\right] .
$$

If

$$
\operatorname{dist}\left(D_{i}\left(A^{H}\right), \operatorname{span}\left\{q_{1}^{(0)}, \ldots, q_{i}^{(0)}\right\}\right)<1, \quad i=1: n,
$$

then it follows from Theorem 7.3.1 that

$$
\operatorname{dist}\left(\operatorname{span}\left\{q_{1}^{(k)}, \ldots, q_{i}^{(k)}\right\}, \operatorname{span}\left\{q_{1}, \ldots, q_{i}\right\}\right) \rightarrow 0
$$

for $i=1: n$. This implies that the matrices $T_{k}$ defined by

$$
T_{k}=Q_{k}^{H} A Q_{k}
$$

are converging to upper triangular form. Thus, it can be said that the method of orthogonal iteration computes a Schur decomposition provided the original iterate $Q_{0} \in \mathbb{C}^{n \times n}$ is not deficient in the sense of (7.3.11).

The QR iteration arises naturally by considering how to compute the matrix $T_{k}$ directly from its predecessor $T_{k-1}$. On the one hand, we have from (7.3.6) and the definition of $T_{k-1}$ that

$$
T_{k-1}=Q_{k-1}^{H} A Q_{k-1}=Q_{k-1}^{H}\left(A Q_{k-1}\right)=\left(Q_{k-1}^{H} Q_{k}\right) R_{k} .
$$

On the other hand,

$$
T_{k}=Q_{k}^{H} A Q_{k}=\left(Q_{k}^{H} A Q_{k-1}\right)\left(Q_{k-1}^{H} Q_{k}\right)=R_{k}\left(Q_{k-1}^{H} Q_{k}\right)
$$

Thus, $T_{k}$ is determined by computing the QR factorization of $T_{k-1}$ and then multiplying the factors together in reverse order, precisely what is done in (7.3.1).

Note that a single QR iteration is an $O\left(n^{3}\right)$ calculation. Moreover, since convergence is only linear (when it exists), it is clear that the method is a prohibitively expensive way to compute Schur decompositions. Fortunately these practical difficulties can be overcome as we show in §7.4 and §7.5.

### 7.3.4 LR Iterations

We conclude with some remarks about power iterations that rely on the LU factorization rather than the QR factorizaton. Let $G_{0} \in \mathbb{C}^{n \times r}$ have rank $r$. Corresponding to (7.3.1) we have the following iteration:

$$
\begin{aligned}
& \text { for } k=1,2, \ldots \\
& \qquad \begin{aligned}
Z_{k} & =A G_{k-1} \\
Z_{k} & =G_{k} R_{k} \quad \text { (LU factorization) }
\end{aligned} \\
& \text { end }
\end{aligned}
$$

Suppose $r=n$ and that we define the matrices $T_{k}$ by

$$
T_{k}=G_{k}^{-1} A G_{k} .
$$

It can be shown that if we set $L_{0}=G_{0}$, then the $T_{k}$ can be generated as follows:

$$
\begin{aligned}
& T_{0}=L_{0}^{-1} A L_{0} \\
& \text { for } k=1,2, \ldots \\
& \qquad T_{k-1}=L_{k} R_{k} \quad \text { (LU factorization) } \\
& T_{k}=R_{k} L_{k} \\
& \text { end }
\end{aligned}
$$

Iterations (7.3.12) and (7.3.14) are known as treppeniteration and the $L R$ iteration, respectively. Under reasonable assumptions, the $T_{k}$ converge to upper triangular form. To successfully implement either method, it is necessary to pivot. See Wilkinson (AEP, p. 602).

## Appendix

In order to establish Theorem 7.3.1 we need the following lemma that bounds powers of a matrix and powers of its inverse.

Lemma 7.3.2. Let $Q^{H} A Q=T=D+N$ be a Schur decomposition of $A \in \mathbb{C}^{n \times n}$ where $D$ is diagonal and $N$ strictly upper triangular. Let $\lambda_{\max }$ and $\lambda_{\min }$ denote the largest and smallest eigenvalues of $A$ in absolute value. If $\mu \geq 0$, then for all $k \geq 0$ we have

$$
\left\|A^{k}\right\|_{2} \leq(1+\mu)^{n-1}\left(\left|\lambda_{\max }\right|+\frac{\|N\|_{F}}{1+\mu}\right)^{k} .
$$

If $A$ is nonsingular and $\mu \geq 0$ satisfies $(1+\mu)\left|\lambda_{\min }\right|>\|N\|_{F}$, then for all $k \geq 0$ we also have

$$
\left\|A^{-k}\right\|_{2} \leq(1+\mu)^{n-1}\left(\frac{1}{\left|\lambda_{\min }\right|-\|N\|_{F} /(1+\mu)}\right)^{k} .
$$

Proof. For $\mu \geq 0$, define the diagonal matrix $\Delta$ by

$$
\Delta=\operatorname{diag}\left(1,(1+\mu),(1+\mu)^{2}, \ldots,(1+\mu)^{n-1}\right)
$$

and note that $\kappa_{2}(\Delta)=(1+\mu)^{n-1}$. Since $N$ is strictly upper triangular, it is easy to verify that

$$
\left\|\Delta N \Delta^{-1}\right\|_{F} \leq \frac{\|N\|_{F}}{1+\mu}
$$

and thus

$$
\begin{aligned}
\left\|A^{k}\right\|_{2} & =\left\|T^{k}\right\|_{2}=\left\|\Delta^{-1}\left(D+\Delta N \Delta^{-1}\right)^{k} \Delta\right\|_{2} \\
& \leq \kappa_{2}(\Delta)\left(\|D\|_{2}+\left\|\Delta N \Delta^{-1}\right\|_{2}\right)^{k} \leq(1+\mu)^{n-1}\left(\left|\lambda_{\max }\right|+\frac{\|N\|_{F}}{1+\mu}\right)^{k} .
\end{aligned}
$$

On the other hand, if $A$ is nonsingular and $(1+\mu)\left|\lambda_{\min }\right|>\|N\|_{F}$, then

$$
\left\|\Delta D^{-1} N \Delta^{-1}\right\|_{2}=\left\|D^{-1}\left(\Delta N \Delta^{-1}\right)\right\|_{2} \leq \frac{1}{\left|\lambda_{\min }\right|}\left\|\Delta N \Delta^{-1}\right\|_{F}<1 .
$$

Using Lemma 2.3.3 we obtain

$$
\begin{aligned}
\left\|A^{-k}\right\|_{2} & =\left\|T^{-k}\right\|_{2}=\left\|\Delta^{-1}\left[\left(I+\Delta D^{-1} N \Delta^{-1}\right)^{-1} D^{-1}\right]^{k} \Delta\right\|_{2} \\
& \leq \kappa_{2}(\Delta)\left(\frac{\left\|D^{-1}\right\|_{2}}{1-\left\|\Delta D^{-1} N \Delta^{-1}\right\|_{2}}\right)^{k} \leq(1+\mu)^{n-1}\left(\frac{1}{|\mu|-\|N\|_{F} /(1+\mu)}\right)^{k}
\end{aligned}
$$

completing the proof of the lemma.
Proof of Theorem 7.3.1. By induction it is easy to show that the matrix $Q_{k}$ in (7.3.6) satisfies

$$
A^{k} Q_{0}=Q_{k}\left(R_{k} \cdots R_{1}\right)
$$

a QR factorization of $A^{k} Q_{0}$. By substituting the Schur decomposition (7.3.7)-(7.3.8) into this equation we obtain

$$
T^{k}\left[\begin{array}{c}
V_{0} \\
W_{0}
\end{array}\right]=\left[\begin{array}{c}
V_{k} \\
W_{k}
\end{array}\right]\left(R_{k} \cdots R_{1}\right)
$$

where

$$
V_{k}=Q_{\alpha}^{H} Q_{k}, \quad W_{k}=Q_{\beta}^{H} Q_{k}
$$

Our goal is to bound $\left\|W_{k}\right\|_{2}$ since by the definition of subspace distance given in §2.5.3 we have

$$
\left\|W_{k}\right\|_{2}=\operatorname{dist}\left(D_{r}(A), \operatorname{ran}\left(Q_{k}\right)\right) .
$$

Note from the thin CS decomposition (Theorem 2.5.2) that

$$
1=d_{k}^{2}+\sigma_{\min }\left(V_{k}\right)^{2}
$$

Since $T_{11}$ and $T_{22}$ have no eigenvalues in common, Lemma 7.1.5 tells us that the Sylvester equation $T_{11} X-X T_{22}=-T_{12}$ has a solution $X \in \mathbb{C}^{r \times(n-r)}$ and that

$$
\|X\|_{F} \leq \frac{\left\|T_{12}\right\|_{F}}{\operatorname{sep}\left(T_{11}, T_{22}\right)}
$$

It follows that

$$
\left[\begin{array}{cc}
I_{r} & X \\
0 & I_{n-r}
\end{array}\right]^{-1}\left[\begin{array}{cc}
T_{11} & T_{12} \\
0 & T_{22}
\end{array}\right]\left[\begin{array}{cc}
I_{r} & X \\
0 & I_{n-r}
\end{array}\right]=\left[\begin{array}{cc}
T_{11} & 0 \\
0 & T_{22}
\end{array}\right]
$$

By substituting this into (7.3.17) we obtain

$$
\left[\begin{array}{cc}
T_{11}^{k} & 0 \\
0 & T_{22}^{k}
\end{array}\right]\left[\begin{array}{c}
V_{0}-X W_{0} \\
W_{0}
\end{array}\right]--\left[\begin{array}{c}
V_{k}-X W_{k} \\
W_{k}
\end{array}\right]\left(R_{k} \cdots R_{1}\right)
$$

i.e.,

$$
\begin{aligned}
T_{11}^{k}\left(V_{0}-X W_{0}\right) & =\left(V_{k}-X W_{k}\right)\left(R_{k} \cdots R_{1}\right) \\
T_{22}^{k} W_{0} & =W_{k}\left(R_{k} \cdots R_{1}\right)
\end{aligned}
$$

The matrix $I+X X^{H}$ is Hermitian positive definite and so it has a Cholesky factorization

$$
I+X X^{H}=G G^{H} .
$$

It is clear that

$$
\sigma_{\min }(G) \geq 1
$$

If the matrix $Z \in \mathbb{C}^{n \times(n-r)}$ is defined by

$$
Z=Q\left[\begin{array}{c}
I_{r} \\
-X^{H}
\end{array}\right] G^{-H}=\left[Q_{\alpha} Q_{\beta}\right]\left[\begin{array}{c}
I_{r} \\
-X^{H}
\end{array}\right] G^{-H}=\left(Q_{\alpha}-Q_{\beta} X^{H}\right) G^{-H},
$$

then it follows from the equation $A^{H} Q=Q T^{H}$ that

$$
A^{H}\left(Q_{\alpha}-Q_{\beta} X^{H}\right)=\left(Q_{\alpha}-Q_{\beta} X^{H}\right) T_{11}^{H} .
$$

Since $Z^{H} Z=I_{r}$ and $\operatorname{ran}(Z)=\operatorname{ran}\left(Q_{\alpha}-Q_{\beta} X^{H}\right)$, it follows that the columns of $Z$ are an orthonormal basis for $D_{r}\left(A^{H}\right)$. Using the CS decomposition, (7.3.19), and the fact that $\operatorname{ran}\left(Q_{\beta}\right)=D_{r}\left(A^{H}\right)^{\perp}$, we have

$$
\begin{aligned}
\sigma_{\min }\left(Z^{T} Q_{0}\right)^{2} & =1-\operatorname{dist}\left(D_{r}\left(A^{H}\right), Q_{0}\right)^{2}=1-\left\|Q_{\beta}^{H} Q_{0}\right\| \\
& =\sigma_{\min }\left(Q_{\alpha}^{T} Q_{0}\right)^{2}=\sigma_{\min }\left(V_{0}\right)^{2}=1-d_{0}^{2}>0
\end{aligned}
$$

This shows that

$$
V_{0}-X W_{0}=\left[I_{r} \mid-X\right]\left[\begin{array}{c}
Q_{\alpha}^{H} Q_{0} \\
Q_{\beta}^{H} Q_{0}
\end{array}\right]=\left(Z G^{H}\right)^{H} Q_{0}=G\left(Z^{H} Q_{0}\right)
$$

is nonsingular and together with (7.3.24) we obtain

$$
\left\|\left(V_{0}-X W_{0}\right)^{-1}\right\|_{2} \leq\left\|G^{-1}\right\|_{2}\left\|\left(Z^{H} Q_{0}\right)^{-1}\right\|_{2} \leq \frac{1}{\sqrt{1-d_{0}^{2}}}
$$

Manipulation of (7.3.19) and (7.3.20) yields

$$
W_{k}=T_{22}^{k} W_{0}\left(R_{k} \cdots R_{1}\right)^{-1}=T_{22}^{k} W_{0}\left(V_{0}-X W_{0}\right)^{-1} T_{11}^{-k}\left(V_{k}-X W_{k}\right) .
$$

The verification of (7.3.10) is completed by taking norms in this equation and using (7.3.18), (7.3.19), (7.3.20), (7.3.26), and the following facts:

$$
\begin{aligned}
\left\|T_{22}^{k}\right\|_{2} & \leq(1+\mu)^{n-r-1}\left(\left|\lambda_{r+1}\right|+\|N\|_{F} /(1+\mu)\right)^{k}, \\
\left\|T_{11}^{-k}\right\|_{2} & \leq(1+\mu)^{r-1} /\left(\left|\lambda_{r}\right|-\|N\|_{F} /(1+\mu)\right)^{k}, \\
\left\|V_{k}-X W_{k}\right\|_{2} & \leq\left\|V_{k}\right\|_{2}+\|X\|_{2}\left\|W_{k}\right\|_{2} \leq 1+\left\|T_{12}\right\|_{F} / \operatorname{sep}\left(T_{11}, T_{22}\right) .
\end{aligned}
$$

The bounds for $\left\|T_{22}^{k}\right\|_{2}$ and $\left\|T_{11}^{-k}\right\|_{2}$ follow from Lemma 7.3.2.

## Problems

P7.3.1 Verify Equation (7.3.5).
P7.3.2 Suppose the eigenvalues of $A \in \mathbf{R}^{n \times n}$ satisfy $\left|\lambda_{1}\right|=\left|\lambda_{2}\right|>\left|\lambda_{3}\right| \geq \cdots \geq\left|\lambda_{n}\right|$ and that $\lambda_{1}$ and $\lambda_{2}$ are complex conjugates of one another. Let $S=\operatorname{span}\{y, z\}$ where $y, z \in \mathbf{R}^{n}$ satisfy $A(y+i z)= \lambda_{1}(y+i z)$. Show how the power method with a real starting vector can be used to compute an approximate basis for $S$.
P7.3.3 Assume $A \in \mathbf{R}^{n \times n}$ has eigenvalues $\lambda_{1}, \ldots, \lambda_{n}$ that satisfy

$$
\lambda=\lambda_{1}=\lambda_{2}=\lambda_{3}=\lambda_{4}>\left|\lambda_{5}\right| \geq \cdots \geq\left|\lambda_{n}\right|
$$

where $\lambda$ is positive. Assume that $A$ has two Jordan blocks of the form.

$$
\left[\begin{array}{ll}
\lambda & 1 \\
0 & \lambda
\end{array}\right]
$$

Discuss the convergence properties of the power method when applied to this matrix and how the convergence might be accelerated.
P7.3.4 A matrix $A$ is a positive matrix if $a_{i j}>0$ for all $i$ and $j$. A vector $v \in \mathbf{R}^{n}$ is a positive vector if $v_{i}>0$ for all $i$. Perron's theorem states that if $A$ is a positive square matrix, then it has a unique dominant eigenvalue equal to its spectral radius $\rho(A)$ and there is a positive vector $x$ so that $A x=\rho(A) \cdot x$. In this context, $x$ is called the Perron vector and $\rho(A)$ is called the Perron root. Assume that $A \in \mathbf{R}^{n \times n}$ is positive and $q \in \mathbf{R}^{n}$ is positive with unit 2 -norm. Consider the following implementation of the power method (7.3.3):

$$
\begin{aligned}
& z=A q, \lambda=q^{T} z \\
& \text { while }\|z-\lambda q\|_{2}>\delta \\
& \text { end } \quad q=z, q=q /\|q\|_{2}, z=A q, \lambda=q^{T} z
\end{aligned}
$$

(a) Adjust the termination criteria to guarantee (in principle) that the final $\lambda$ and $q$ satisfy $\tilde{A} q=\lambda q$, where $\|\tilde{A}-A\|_{2} \leq \delta$ and $\tilde{A}$ is positive. (b) Applied to a positive matrix $A \in \mathbf{R}^{n \times n}$, the CollatzWielandt formula states that $\rho(A)$ is the maximum value of the function $f$ defined by

$$
f(x)=\min _{1 \leq i \leq n} \frac{y_{i}}{x_{i}}
$$

where $x \in \mathbf{R}^{n}$ is positive and $y=A x$. Does it follow that $f(A q) \geq f(q)$ ? In other words, do the iterates $\left\{q^{(k)}\right\}$ in the power method have the property that $f\left(q^{(k)}\right)$ increases monotonically to the Perron root, assuming that $q^{(0)}$ is positive?
P7.3.5 (Read the previous problem for background.) A matrix $A$ is a nonnegative matrix if $a_{i j} \geq 0$ for all $i$ and $j$. A matrix $A \in \mathbf{R}^{n \times n}$ is reducible if there is a permutation $P$ so that $P^{T} A P$ is block triangular with two or more square diagonal blocks. A matrix that is not reducible is irreducible. The Perron-Frobenius theorem states that if $A$ is a square, nonnegative, and irreducible, then $\rho(A)$, the Perron root, is an eigenvalue for $A$ and there is a positive vector $x$, the Perron vector, so that $A x=\rho(A) \cdot x$. Assume that $A_{1}, A_{2}, A_{3} \in \mathbf{R}^{n \times n}$ are each positive and let the nonnegative matrix $A$ be defined by

$$
A=\left[\begin{array}{ccc}
0 & A_{1} & 0 \\
0 & 0 & A_{2} \\
A_{3} & 0 & 0
\end{array}\right]
$$

(a) Show that $A$ is irreducible. (b) Let $B=A_{1} A_{2} A_{3}$. Show how to compute the Perron root and vector for $A$ from the Perron root and vector for $B$. (c) Show that $A$ has other eigenvalues with absolute value equal to the Perron root. How could those eigenvalues and the associated eigenvectors be computed?
P7.3.6 (Read the previous two problems for background.) A nonnegative matrix $P \in \mathbf{R}^{n \times n}$ is stochastic if the entries in each column sum to 1 . A vector $v \in \mathbf{R}^{n}$ is a probability vector if its entries are nonnegative and sum to 1 . (a) Show that if $P \in \mathbf{R}^{n \times n}$ is stochastic and $v \in \mathbf{R}^{n}$ is a probability vector, then $w=P v$ is also a probability vector. (b) The entries in a stochastic matrix $P \in \mathbf{R}^{n \times n}$ can
be regarded as the transition probabilities associated with an $n$-state Markov Chain. Let $v_{j}$ be the probability of being in state $j$ at time $t=t_{\text {current }}$. In the Markov model, the probability of being in state $i$ at time $t=t_{\text {next }}$ is given by

$$
w_{i}=\sum_{j=1}^{n} p_{i j} v_{j} \quad i=1: n,
$$

i.e., $w=P v$. With the help of a biased coin, a surfer on the World Wide Web randomly jumps from page to page. Assume that the surfer is currently viewing web page $j$ and that the coin comes up heads with probability $\alpha$. Here is how the surfer determines the next page to visit:

Step 1. A coin is tossed.
Step 2. If it comes up heads and web page $j$ has at least one outlink, then the next page to visit is randomly selected from the list of outlink pages.

Step 3. Otherwise, the next page to visit is randomly selected from the list of all possible pages.
Let $P \in \mathbf{R}^{n \times n}$ be the matrix of transition probabilities that define this random process. Specify $P$ in terms of $\alpha$, the vector of ones $e$, and the link matrix $H \in \mathbf{R}^{n \times n}$ defined by

$$
h_{i j}= \begin{cases}1 & \text { if there is a link on web page } j \text { to web page } i \\ 0 & \text { otherwise }\end{cases}
$$

Hints: The number of nonzero components in $H(:, j)$ is the number of outlinks on web page $j . P$ is a convex combination of a very sparse sparse matrix and a very dense rank-1 matrix. (c) Detail how the power method can be used to determine a probability vector $x$ so that $P x=x$. Strive to get as much computation "outside the loop" as possible. Note that in the limit we can expect to find the random surfer viewing web page $i$ with probability $x_{i}$. Thus, a case can be made that more important pages are associated with the larger components of $x$. This is the basis of Google PageRank. If

$$
x_{i_{1}} \geq x_{i_{2}} \geq \cdots \geq x_{i_{n}}
$$

then web page $i_{k}$ has page rank $k$.
P7.3.7 (a) Show that if $X \in \mathbb{C}^{n \times n}$ is nonsingular, then

$$
\|A\|_{X}=\left\|X^{-1} A X\right\|_{2}
$$

defines a matrix norm with the property that

$$
\|A B\|_{X} \leq\|A\|_{X}\|B\|_{X}
$$

(b) Show that for any $\epsilon>0$ there exists a nonsingular $X \in \mathbb{C}^{n \times n}$ such that

$$
\|A\|_{X}=\left\|X^{-1} A X\right\|_{2} \leq \rho(A)+\epsilon
$$

where $\rho(A)$ is $A$ 's spectral radius. Conclude that there is a constant $M$ such that

$$
\left\|A^{k}\right\|_{2} \leq M(\rho(A)+\epsilon)^{k}
$$

for all non-negative integers $k$. (Hint: Set $X=Q \operatorname{diag}\left(1, a, \ldots, a^{n-1}\right)$ where $Q^{H} A Q=D+N$ is $A$ 's Schur decomposition.)
P7.3.8 Verify that (7.3.14) calculates the matrices $T_{k}$ defined by (7.3.13).
P7.3.9 Suppose $A \in \mathbb{C}^{n \times n}$ is nonsingular and that $Q_{0} \in \mathbb{C}^{n \times p}$ has orthonormal columns. The following iteration is referred to as inverse orthogonal iteration.
for $k=1,2, \ldots$
Solve $A Z_{k}=Q_{k-1}$ for $Z_{k} \in \mathbb{C}^{n \times p}$

$$
Z_{k}=Q_{k} R_{k} \quad \text { (QR factorization) }
$$

end
Explain why this iteration can usually be used to compute the $p$ smallest eigenvalues of $A$ in absolute value. Note that to implement this iteration it is necessary to be able to solve linear systems that involve $A$. If $p=1$, the method is referred to as the inverse power method.

## Notes and References for §7.3

For an excellent overview of the QR iteration and related procedures, see Watkins (MEP), Stewart (MAE), and Kressner (NMSE). A detailed, practical discussion of the power method is given in Wilkinson (AEP, Chap. 10). Methods are discussed for accelerating the basic iteration, for calculating nondominant eigenvalues, and for handling complex conjugate eigenvalue pairs. The connections among the various power iterations are discussed in:
B.N. Parlett and W.G. Poole (1973). "A Geometric Theory for the QR, LU, and Power Iterations," SIAM J. Numer. Anal. 10, 389-412.

The QR iteration was concurrently developed in:
J.G.F. Francis (1961). "The QR Transformation: A Unitary Analogue to the LR Transformation," Comput. J. 4, 265-71, 332-334.
V.N. Kublanovskaya (1961). "On Some Algorithms for the Solution of the Complete Eigenvalue Problem," USSR Comput. Math. Phys. 3, 637-657.

As can be deduced from the title of the first paper by Francis, the LR iteration predates the QR iteration. The former very fundamental algorithm was proposed by:
H. Rutishauser (1958). "Solution of Eigenvalue Problems with the LR Transformation," Nat. Bur. Stand. Appl. Math. Ser. 49, 47-81.

More recent, related work includes:
B.N. Parlett (1995). "The New qd Algorithms," Acta Numerica 5, 459-491.
C. Ferreira and B.N. Parlett (2009). "Convergence of the LR Algorithm for a One-Point Spectrum Tridiagonal Matrix," Numer. Math. 113, 417-431.

Numerous papers on the convergence and behavior of the QR iteration have appeared, see:
J.H. Wilkinson (1965). "Convergence of the LR, QR, and Related Algorithms," Comput. J. 8, 77-84.
B.N. Parlett (1965). "Convergence of the Q-R Algorithm," Numer. Math. 7, 187-93. (Correction in Numer. Math. 10, 163-164.)
B.N. Parlett (1966). "Singular and Invariant Matrices Under the QR Algorithm," Math. Comput. 20, 611-615.
B.N. Parlett (1968). "Global Convergence of the Basic QR Algorithm on Hessenberg Matrices," Math. Comput. 22, 803-817.
D.S. Watkins (1982). "Understanding the QR Algorithm," SIAM Review 24, 427-440.
T. Nanda (1985). "Differential Equations and the QR Algorithm," SIAM J. Numer. Anal. 22, 310-321.
D.S. Watkins (1993). "Some Perspectives on the Eigenvalue Problem," SIAM Review 35, 430-471.
D.S. Watkins (2008). "The QR Algorithm Revisited," SIAM Review 50, 133-145.
D.S. Watkins (2011). "Francis's Algorithm," AMS Monthly 118, 387-403.

A block analog of the QR iteration is discussed in:
M. Robbè and M. Sadkane (2005). "Convergence Analysis of the Block Householder Block Diagonalization Algorithm," BIT 45, 181-195.

The following references are concerned with various practical and theoretical aspects of simultaneous iteration:
H. Rutishauser (1970). "Simultaneous Iteration Method for Symmetric Matrices," Numer. Math. 16, 205-223.
M. Clint and A. Jennings (1971). "A Simultaneous Iteration Method for the Unsymmetric Eigenvalue Problem," J. Inst. Math. Applic. 8, 111-121.
G.W. Stewart (1976). "Simultaneous Iteration for Computing Invariant Subspaces of Non-Hermitian Matrices," Numer. Math. 25, 123-136.
A. Jennings (1977). Matrix Computation for Engineers and Scientists, John Wiley and Sons, New York.
Z. Bai and G.W. Stewart (1997). "Algorithm 776: SRRIT: a Fortran Subroutine to Calculate the Dominant Invariant Subspace of a Nonsymmetric Matrix," ACM Trans. Math. Softw. 23, 494513.

Problems P7.3.4-P7.3.6 explore the relevance of the power method to the problem of computing the Perron root and vector of a nonnegative matrix. For further background and insight, see:
A. Berman and R.J. Plemmons (1994). Nonnegative Matrices in the Mathematical Sciences, SIAM Publications,Philadelphia, PA.
A.N. Langville and C.D. Meyer (2006). Google's PageRank and Beyond, Princeton University Press, Princeton and Oxford. .

The latter volume is outstanding in how it connects the tools of numerical linear algebra to the design and analysis of Web browsers. See also:
W.J. Stewart (1994). Introduction to the Numerical Solution of Markov Chains, Princeton University Press, Princeton, NJ.
M.W. Berry, Z. Drmač, and E.R. Jessup (1999). "Matrices, Vector Spaces, and Information Retrieval," SIAM Review 41, 335-362.
A.N. Langville and C.D. Meyer (2005). "A Survey of Eigenvector Methods for Web Information Retrieval," SIAM Review 47, 135-161.
A.N. Langville and C.D. Meyer (2006). "A Reordering for the PageRank Problem", SIAM J. Sci. Comput. 27, 2112-2120.
A.N. Langville and C.D. Meyer (2006). "Updating Markov Chains with an Eye on Google's PageRank," SIAM J. Matrix Anal. Applic. 27, 968-987.

### 7.4 The Hessenberg and Real Schur Forms

In this and the next section we show how to make the QR iteration (7.3.1) a fast, effective method for computing Schur decompositions. Because the majority of eigenvalue/invariant subspace problems involve real data, we concentrate on developing the real analogue of (7.3.1) which we write as follows:

$$
\begin{aligned}
& H_{0}=U_{0}^{T} A U_{0} \\
& \text { for } k=1,2, \ldots \\
& \qquad \begin{array}{l}
H_{k-1}=U_{k} R_{k} \quad \text { (QR factorization) } \\
H_{k}=R_{k} U_{k}
\end{array} \\
& \text { end }
\end{aligned}
$$

Here, $A \in \mathbb{R}^{n \times n}$, each $U_{k} \in \mathbb{R}^{n \times n}$ is orthogonal, and each $R_{k} \in \mathbb{R}^{n \times n}$ is upper triangular. A difficulty associated with this real iteration is that the $H_{k}$ can never converge to triangular form in the event that $A$ has complex eigenvalues. For this reason, we must lower our expectations and be content with the calculation of an alternative decomposition known as the real Schur decomposition.

In order to compute the real Schur decomposition efficiently we must carefully choose the initial orthogonal similarity transformation $U_{0}$ in (7.4.1). In particular, if we choose $U_{0}$ so that $H_{0}$ is upper Hessenberg, then the amount of work per iteration is reduced from $O\left(n^{3}\right)$ to $O\left(n^{2}\right)$. The initial reduction to Hessenberg form (the $U_{0}$ computation) is a very important computation in its own right and can be realized by a sequence of Householder matrix operations.

### 7.4.1 The Real Schur Decomposition

A block upper triangular matrix with either 1-by-1 or 2-by-2 diagonal blocks is upper quasi-triangular. The real Schur decomposition amounts to a real reduction to upper quasi-triangular form.

Theorem 7.4.1 (Real Schur Decomposition). If $A \in \mathbb{R}^{n \times n}$, then there exists an orthogonal $Q \in \mathbb{R}^{n \times n}$ such that

$$
Q^{T} A Q=\left[\begin{array}{cccc}
R_{11} & R_{12} & \cdots & R_{1 m} \\
0 & R_{22} & \cdots & R_{2 m} \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & R_{m m}
\end{array}\right]
$$

where each $R_{i i}$ is either a 1-by-1 matrix or a 2-by-2 matrix having complex conjugate eigenvalues.

Proof. The complex eigenvalues of $A$ occur in conjugate pairs since the characteristic polynomial $\operatorname{det}(z I-A)$ has real coefficients. Let $k$ be the number of complex conjugate pairs in $\lambda(A)$. We prove the theorem by induction on $k$. Observe first that Lemma 7.1.2 and Theorem 7.1.3 have obvious real analogs. Thus, the theorem holds if $k=0$. Now suppose that $k \geq 1$. If $\lambda=\gamma+i \mu \in \lambda(A)$ and $\mu \neq 0$, then there exist vectors $y$ and $z$ in $\mathbb{R}^{n}(z \neq 0)$ such that $A(y+i z)=(\gamma+i \mu)(y+i z)$, i.e.,

$$
A[y \mid z]=[y \mid z]\left[\begin{array}{rr}
\gamma & \mu \\
-\mu & \gamma
\end{array}\right]
$$

The assumption that $\mu \neq 0$ implics that $y$ and $z$ span a 2 -dimensional, real invariant subspace for $A$. It then follows from Lemma 7.1.2 that an orthogonal $U \in \mathbb{R}^{n \times n}$ exists such that

$$
U^{T} A U=\left[\begin{array}{cc}
T_{11} & T_{12} \\
0 & T_{22} \\
2 & n-2
\end{array}\right]_{n-2}^{2}
$$

where $\lambda\left(T_{11}\right)=\{\lambda, \bar{\lambda}\}$. By induction, there exists an orthogonal $\tilde{U}$ so $\tilde{U}^{T} T_{22} \tilde{U}$ has the required structure. The theorem follows by setting $Q=U \cdot \operatorname{diag}\left(I_{2}, \tilde{U}\right)$.

The theorem shows that any real matrix is orthogonally similar to an upper quasitriangular matrix. It is clear that the real and imaginary parts of the complex eigenvalues can be easily obtained from the 2 -by- 2 diagonal blocks. Thus, it can be said that the real Schur decomposition is an eigenvalue-revealing decomposition.

### 7.4.2 A Hessenberg QR Step

We now turn our attention to the efficient execution of a single QR step in (7.4.1). In this regard, the most glaring shortcoming associated with (7.4.1) is that each step requires a full QR factorization costing $O\left(n^{3}\right)$ flops. Fortunately, the amount of work per iteration can be reduced by an order of magnitude if the orthogonal matrix $U_{0}$ is judiciously chosen. In particular, if $U_{0}^{T} A U_{0}=H_{0}=\left(h_{i j}\right)$ is upper Hessenberg ( $h_{i j}=0$, $i>j+1$ ), then each subsequent $H_{k}$ requires only $O\left(n^{2}\right)$ flops to calculate. To see this we look at the computations $H=Q R$ and $H_{+}=R Q$ when $H$ is upper Hessenberg. As described in §5.2.5, we can upper triangularize $H$ with a sequence of $n-1$ Givens rotations: $Q^{T} H \equiv G_{n-1}^{T} \cdots G_{1}^{T} H=R$. Here, $G_{i}=G\left(i, i+1, \theta_{i}\right)$. For the $n=4$ case there are three Givens premultiplications:

$$
\left[\begin{array}{cccc}
\times & \times & \times & \times \\
\times & \times & \times & \times \\
0 & \times & \times & \times \\
0 & 0 & \times & \times
\end{array}\right] \rightarrow\left[\begin{array}{llll}
\times & \times & \times & \times \\
0 & \times & \times & \times \\
0 & \times & \times & \times \\
0 & 0 & \times & \times
\end{array}\right] \rightarrow\left[\begin{array}{cccc}
\times & \times & \times & \times \\
0 & \times & \times & \times \\
0 & 0 & \times & \times \\
0 & 0 & \times & \times
\end{array}\right] \rightarrow\left[\begin{array}{cccc}
\times & \times & \times & \times \\
0 & \times & \times & \times \\
0 & 0 & \times & \times \\
0 & 0 & 0 & \times
\end{array}\right] .
$$

See Algorithm 5.2.5. The computation $R Q=R\left(G_{1} \cdots G_{n-1}\right)$ is equally easy to implement. In the $n=4$ case there are three Givens post-multiplications:

$$
\left[\begin{array}{cccc}
\times & \times & \times & \times \\
0 & \times & \times & \times \\
0 & 0 & \times & \times \\
0 & 0 & 0 & \times
\end{array}\right] \rightarrow\left[\begin{array}{cccc}
\times & \times & \times & \times \\
\times & \times & \times & \times \\
0 & 0 & \times & \times \\
0 & 0 & 0 & \times
\end{array}\right] \rightarrow\left[\begin{array}{cccc}
\times & \times & \times & \times \\
\times & \times & \times & \times \\
0 & \times & \times & \times \\
0 & 0 & 0 & \times
\end{array}\right] \rightarrow\left[\begin{array}{cccc}
\times & \times & \times & \times \\
\times & \times & \times & \times \\
0 & \times & \times & \times \\
0 & 0 & \times & \times
\end{array}\right] .
$$

Overall we obtain the following algorithm:
Algorithm 7.4.1 If $H$ is an $n$-by- $n$ upper Hessenberg matrix, then this algorithm overwrites $H$ with $H_{+}=R Q$ where $H=Q R$ is the QR factorization of $H$.
for $k=1: n-1$

$$
\begin{aligned}
& {\left[c_{k}, s_{k}\right]=\operatorname{givens}(H(k, k), H(k+1, k))} \\
& H(k: k+1, k: n)=\left[\begin{array}{rr}
c_{k} & s_{k} \\
-s_{k} & c_{k}
\end{array}\right]^{T} H(k: k+1, k: n)
\end{aligned}
$$

end
for $k=1: n-1$

$$
H(1: k+1, k: k+1)=H(1: k+1, k: k+1)\left[\begin{array}{rr}
c_{k} & s_{k} \\
-s_{k} & c_{k}
\end{array}\right]
$$

end

Let $G_{k}=G\left(k, k+1, \theta_{k}\right)$ be the $k$ th Givens rotation. It is easy to confirm that the matrix $Q=G_{1} \cdots G_{n-1}$ is upper Hessenberg. Thus, $R Q=H_{+}$is also upper Hessenberg. The algorithm requires about $6 n^{2}$ flops, an order of magnitude more efficient than a full matrix QR step (7.3.1).

### 7.4.3 The Hessenberg Reduction

It remains for us to show how the Hessenberg decomposition

$$
U_{0}^{T} A U_{0}=H, \quad U_{0}^{T} U_{0}=I
$$

can be computed. The transformation $U_{0}$ can be computed as a product of Householder matrices $P_{1}, \ldots, P_{n-2}$. The role of $P_{k}$ is to zero the $k$ th column below the subdiagonal. In the $n=6$ case, we have

$$
\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times
\end{array}\right] \xrightarrow{P_{1}}\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times
\end{array}\right] \xrightarrow{P_{2}}
$$

$$
\left[\begin{array}{llllll}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times
\end{array}\right] \xrightarrow{P_{3}}\left[\begin{array}{llllll}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times
\end{array}\right] \xrightarrow{P_{4}}\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times
\end{array}\right] .
$$

In general, after $k-1$ steps we have computed $k-1$ Householder matrices $P_{1}, \ldots, P_{k-1}$ such that

$$
\left(P_{1} \cdots P_{k-1}\right)^{T} A\left(P_{1} \cdots P_{k-1}\right)=\left[\begin{array}{ccc}
B_{11} & B_{12} & B_{13} \\
B_{21} & B_{22} & B_{23} \\
0 & B_{32} & B_{33}
\end{array}\right]_{n-k}{ }_{k-1}
$$

is upper Hessenberg through its first $k-1$ columns. Suppose $\tilde{P}_{k}$ is an order- $(n-k)$ Householder matrix such that $\tilde{P}_{k} B_{32}$ is a multiple of $e_{1}^{(n-k)}$. If $P_{k}=\operatorname{diag}\left(I_{k}, \tilde{P}_{k}\right)$, then

$$
\left(P_{1} \cdots P_{k}\right)^{T} A\left(P_{1} \cdots P_{k}\right)=\left[\begin{array}{ccc}
B_{11} & B_{12} & B_{13} \tilde{P}_{k} \\
B_{21} & B_{22} & B_{23} \tilde{P}_{k} \\
0 & \tilde{P}_{k} B_{32} & \tilde{P}_{k} B_{33} \tilde{P}_{k}
\end{array}\right]
$$

is upper Hessenberg through its first $k$ columns. Repeating this for $k=1: n-2$ we obtain

Algorithm 7.4.2 (Householder Reduction to Hessenberg Form) Given $A \in \mathbb{R}^{n \times n}$, the following algorithm overwrites $A$ with $H=U_{0}^{T} A U_{0}$ where $H$ is upper Hessenberg and $U_{0}$ is a product of Householder matrices.

```
for $k=1: n-2$
    $[v, \beta]=\operatorname{house}(A(k+1: n, k))$
    $A(k+1: n, k: n)=\left(I-\beta v v^{T}\right) A(k+1: n, k: n)$
    $A(1: n, k+1: n)=A(1: n, k+1: n)\left(I-\beta v v^{T}\right)$
end
```

This algorithm requires $10 n^{3} / 3$ flops. If $U_{0}$ is explicitly formed, an additional $4 n^{3} / 3$ flops are required. The $k$ th Householder matrix can be represented in $A(k+2: n, k)$. See Martin and Wilkinson (1968) for a detailed description.

The roundoff properties of this method for reducing $A$ to Hessenberg form are very desirable. Wilkinson (AEP, p. 351) states that the computed Hessenberg matrix $\hat{H}$ satisfies

$$
\hat{H}=Q^{T}(A+E) Q,
$$

where $Q$ is orthogonal and $\|E\|_{F} \leq c n^{2} \mathbf{u}\|A\|_{F}$ with $c$ a small constant.

### 7.4.4 Level-3 Aspects

The Hessenberg reduction (Algorithm 7.4.2) is rich in level-2 operations: half gaxpys and half outer product updates. We briefly mention two ideas for introducing level-3 computations into the process.

The first involves a block reduction to block Hessenberg form and is quite straightforward. Suppose (for clarity) that $n=r N$ and write

$$
A=\left[\begin{array}{ll}
A_{11} & A_{12} \\
A_{21} & A_{22} \\
r & n-r
\end{array}\right]_{n-r}^{r}
$$

Suppose that we have computed the QR factorization $A_{21}=\tilde{Q}_{1} R_{1}$ and that $\tilde{Q}_{1}$ is in WY form. That is, we have $W_{1}, Y_{1} \in \mathbb{R}^{(n-r) \times r}$ such that $\tilde{Q}_{1}=I+W_{1} Y_{1}^{T}$. (See §5.2.2 for details.) If $Q_{1}=\operatorname{diag}\left(I_{r}, \tilde{Q}_{1}\right)$ then

$$
Q_{1}^{T} A Q_{1}=\left[\begin{array}{cc}
A_{11} & A_{12} \tilde{Q}_{1} \\
R_{1} & \tilde{Q}_{1}^{T} A_{22} \tilde{Q}_{1}
\end{array}\right]
$$

Notice that the updates of the ( 1,2 ) and ( 2,2 ) blocks are rich in level- 3 operations given that $\tilde{Q}_{1}$ is in WY form. This fully illustrates the overall process as $Q_{1}^{T} A Q_{1}$ is block upper Hessenberg through its first block column. We next repeat the computations on the first $r$ columns of $\tilde{Q}_{1}^{T} A_{22} \tilde{Q}_{1}$. After $N-1$ such steps we obtain

$$
H=U_{0}^{T} A U_{0}=\left[\begin{array}{ccccc}
H_{11} & H_{12} & \cdots & \cdots & H_{1 N} \\
H_{21} & H_{22} & \cdots & \cdots & H_{2 N} \\
0 & \ddots & \ddots & \cdots & \vdots \\
\vdots & \vdots & \ddots & \ddots & \vdots \\
0 & 0 & \cdots & H_{N, N-1} & H_{N N}
\end{array}\right]
$$

where each $H_{i j}$ is $r$-by- $r$ and $U_{0}=Q_{1} \cdots Q_{N-2}$ with each $Q_{i}$ in WY form. The overall algorithm has a level-3 fraction of the form $1-O(1 / N)$. Note that the subdiagonal blocks in $H$ are upper triangular and so the matrix has lower bandwidth $r$. It is possible to reduce $H$ to actual Hessenberg form by using Givens rotations to zero all but the first subdiagonal.

Dongarra, Hammarling, and Sorensen (1987) have shown how to proceed directly to Hessenberg form using a mixture of gaxpys and level-3 updates. Their idea involves minimal updating after each Householder transformation is generated. For example, suppose the first Householder $P_{1}$ has been computed. To generate $P_{2}$ we need just the second column of $P_{1} A P_{1}$, not the full outer product update. To generate $P_{3}$ we need just the thirrd column of $P_{2} P_{1} A P_{1} P_{2}$, etc. In this way, the Householder matrices can be determined using only gaxpy operations. No outer product updates are involved. Once a suitable number of Householder matrices are known they can be aggregated and applied in level-3 fashion.

For more about the challenges of organizing a high-performance Hessenberg reduction, see Karlsson (2011).

### 7.4.5 Important Hessenberg Matrix Properties

The Hessenberg decomposition is not unique. If $Z$ is any $n$-by- $n$ orthogonal matrix and we apply Algorithm 7.4.2 to $Z^{T} A Z$, then $Q^{T} A Q=H$ is upper Hessenberg where $Q=Z U_{0}$. However, $Q e_{1}=Z\left(U_{0} e_{1}\right)=Z e_{1}$ suggesting that $H$ is unique once the first column of $Q$ is specified. This is essentially the case provided $H$ has no zero subdiagonal entries. Hessenberg matrices with this property are said to be unreduced. Here is important theorem that clarifies these issues.

Theorem 7.4.2 ( Implicit $\mathbf{Q}$ Theorem ). Suppose $Q=\left[q_{1}|\cdots| q_{n}\right]$ and $V= \left[v_{1}|\cdots| v_{n}\right]$ are orthogonal matrices with the property that the matrices $Q^{T} A Q=H$ and $V^{T} A V=G$ are each upper Hessenberg where $A \in \mathbb{R}^{n \times n}$. Let $k$ denote the smallest positive integer for which $h_{k+1, k}=0$, with the convention that $k=n$ if $H$ is unreduced. If $q_{1}=v_{1}$, then $q_{i}= \pm v_{i}$ and $\left|h_{i, i-1}\right|=\left|g_{i, i-1}\right|$ for $i=2: k$. Moreover, if $k<n$, then $g_{k+1, k}=0$.

Proof. Define the orthogonal matrix $W=\left[w_{1}|\cdots| w_{n}\right]=V^{T} Q$ and observe that $G W=W H$. By comparing column $i-1$ in this equation for $i=2: k$ we sce that

$$
h_{i, i-1} w_{i}=G w_{i-1}-\sum_{j=1}^{i-1} h_{j, i-1} w_{j} .
$$

Since $w_{1}=e_{1}$, it follows that $\left[w_{1}|\cdots| w_{k}\right]$ is upper triangular and so for $i=2: k$ we have $w_{i}= \pm I_{n}(:, i)= \pm e_{i}$. Since $w_{i}=V^{T} q_{i}$ and $h_{i, i-1}=w_{i}^{T} G w_{i-1}$ it follows that $v_{i}= \pm q_{i}$ and

$$
\left|h_{i, i-1}\right|=\left|q_{i}^{T} A q_{i-1}\right|=\left|v_{i}^{T} A v_{i-1}\right|=\left|g_{i, i-1}\right|
$$

for $i=2: k$. If $k<n$, then

$$
\begin{aligned}
g_{k+1, k} & =e_{k+1}^{T} G e_{k}= \pm e_{k+1}^{T} G W e_{k}= \pm e_{k+1}^{T} W H e_{k} \\
& = \pm e_{k+1}^{T} \sum_{i=1}^{k} h_{i k} W e_{i}= \pm \sum_{i=1}^{k} h_{i k} e_{k+1}^{T} e_{i}=0
\end{aligned}
$$

completing the proof of the theorem.
The gist of the implicit $Q$ theorem is that if $Q^{T} A Q=H$ and $Z^{T} A Z=G$ are each unreduced upper Hessenberg matrices and $Q$ and $Z$ have the same first column, then $G$ and $H$ are "essentially equal" in the sense that $G=D^{-1} H D$ where $D=\operatorname{diag}( \pm 1, \ldots, \pm 1)$.

Our next theorem involves a new type of matrix called a Krylov matrix. If $A \in \mathbb{R}^{n \times n}$ and $v \in \mathbb{R}^{n}$, then the Krylov matrix $K(A, v, j) \in \mathbb{R}^{n \times j}$ is defined by

$$
K(A, v, j)=\left[v|A v| \ldots \mid A^{j-1} v\right]
$$

It turns out that there is a connection between the Hessenberg reduction $Q^{T} A Q=H$ and the QR factorization of the Krylov matrix $K(A, Q(:, 1), n)$.

Theorem 7.4.3. Suppose $Q \in \mathbb{R}^{n \times n}$ is an orthogonal matrix and $A \in \mathbb{R}^{n \times n}$. Then $Q^{T} A Q=H$ is an unreduced upper Hessenberg matrix if and only if $Q^{T} K(A, Q(:, 1), n)= R$ is nonsingular and upper triangular.

Proof. Suppose $Q \in \mathbb{R}^{n \times n}$ is orthogonal and set $H=Q^{T} A Q$. Consider the identity

$$
Q^{T} K(A, Q(:, 1), n)=\left[e_{1}\left|H e_{1}\right| \ldots \mid H^{n-1} e_{1}\right] \equiv R
$$

If $H$ is an unreduced upper Hessenberg matrix, then it is clear that $R$ is upper triangular with $r_{i i}=h_{21} h_{32} \cdots h_{i, i-1}$ for $i=2: n$. Since $r_{11}=1$ it follows that $R$ is nonsingular.

To prove the converse, suppose $R$ is upper triangular and nonsingular. Since $R(:, k+1)=H R(:, k)$ it follows that $H(:, k) \in \operatorname{span}\left\{e_{1}, \ldots, e_{k+1}\right\}$. This implies that $H$ is upper Hessenberg. Since $r_{n n}=h_{21} h_{32} \cdots h_{n, n-1} \neq 0$ it follows that $H$ is also unreduced.

Thus, there is more or less a correspondence between nonsingular Krylov matrices and orthogonal similarity reductions to unreduced Hessenberg form.

Our last result is about the geometric multiplicity of an eigenvalue of an unreduced upper Hessenberg matrix.

Theorem 7.4.4. If $\lambda$ is an eigenvalue of an unreduced upper Hessenberg matrix $H \in \mathbb{R}^{n \times n}$, then its geometric multiplicity is 1 .

Proof. For any $\lambda \in \mathbb{C}$ we have $\operatorname{rank}(A-\lambda I) \geq n-1$ because the first $n-1$ columns of $H-\lambda I$ are independent.

### 7.4.6 Companion Matrix Form

Just as the Schur decomposition has a nonunitary analogue in the Jordan decomposition, so does the Hessenberg decomposition have a nonunitary analog in the companion matrix decomposition. Let $x \in \mathbb{R}^{n}$ and suppose that the Krylov matrix $K=K(A, x, n)$ is nonsingular. If $c=c(0: n-1)$ solves the linear system $K c=-A^{n} x$, then it follows that $A K=K C$ where $C$ has the form

$$
C=\left[\begin{array}{ccccl}
0 & 0 & \cdots & 0 & -c_{0} \\
1 & 0 & \cdots & 0 & -c_{1} \\
0 & 1 & \cdots & 0 & -c_{2} \\
\vdots & \vdots & \vdots & \vdots & \vdots \\
0 & 0 & \cdots & 1 & -c_{n-1}
\end{array}\right] .
$$

The matrix $C$ is said to be a companion matrix. Since

$$
\operatorname{det}(z I-C)=c_{0}+c_{1} z+\cdots+c_{n-1} z^{n-1}+z^{n}
$$

it follows that if $K$ is nonsingular, then the decomposition $K^{-1} A K=C$ displays $A$ 's characteristic polynomial. This, coupled with the sparseness of $C$, leads to "companion matrix methods" in various application areas. These techniques typically involve:

Step 1. Compute the Hessenberg decomposition $U_{0}^{T} A U_{0}=H$.
Step 2. Hope $H$ is unreduced and set $Y=\left[e_{1}\left|H e_{1}\right| \ldots \mid H^{n-1} e_{1}\right]$.
Step 3. Solve $Y C=H Y$ for $C$.

Unfortunately, this calculation can be highly unstable. $A$ is similar to an unreduced Hessenberg matrix only if each eigenvalue has unit geometric multiplicity. Matrices that have this property are called nonderogatory. It follows that the matrix $Y$ above can be very poorly conditioned if $A$ is close to a derogatory matrix.

A full discussion of the dangers associated with companion matrix computation can be found in Wilkinson (AEP, pp. 405ff.).

## Problems

P7.4.1 Suppose $A \in \mathbf{R}^{n \times n}$ and $z \in \mathbf{R}^{n}$. Give a detailed algorithm for computing an orthogonal $Q$ such that $Q^{T} A Q$ is upper Hessenberg and $Q^{T} z$ is a multiple of $e_{1}$. Hint: Reduce $z$ first and then apply Algorithm 7.4.2.
P7.4.2 Develop a similarity reduction to Hessenberg form using Gauss transforms with pivoting. How many flops are required. See Businger (1969).
P7.4.3 In some situations, it is necessary to solve the linear system $(A+z I) x=b$ for many different values of $z \in \mathbf{R}$ and $b \in \mathbf{R}^{n}$. Show how this problem can be efficiently and stably solved using the Hessenberg decomposition.
P7.4.4 Suppose $H \in \mathbf{R}^{n \times n}$ is an unreduced upper Hessenberg matrix. Show that there exists a diagonal matrix $D$ such that each subdiagonal element of $D^{-1} H D$ is equal to 1 . What is $\kappa_{2}(D)$ ?
P7.4.5 Suppose $W, Y \in \mathbf{R}^{n \times n}$ and define the matrices $C$ and $B$ by

$$
C=W+i Y, \quad B=\left[\begin{array}{cc}
W & -Y \\
Y & W
\end{array}\right]
$$

Show that if $\lambda \in \lambda(C)$ is real, then $\lambda \in \lambda(B)$. Relate the corresponding eigenvectors.
P7.4.6 Suppose

$$
A=\left[\begin{array}{ll}
w & x \\
y & z
\end{array}\right]
$$

is a real matrix having eigenvalues $\lambda \pm i \mu$, where $\mu$ is nonzero. Give an algorithm that stably determines $c=\cos (\theta)$ and $s=\sin (\theta)$ such that

$$
\left[\begin{array}{rr}
c & s \\
-s & c
\end{array}\right]^{T}\left[\begin{array}{ll}
w & x \\
y & z
\end{array}\right]\left[\begin{array}{rr}
c & s \\
-s & c
\end{array}\right]=\left[\begin{array}{ll}
\lambda & \beta \\
\alpha & \lambda
\end{array}\right]
$$

where $\boldsymbol{\alpha} \boldsymbol{\beta}=-\boldsymbol{\mu}^{2}$.
P7.4.7 Suppose ( $\lambda, x$ ) is a known eigenvalue-eigenvector pair for the upper Hessenberg matrix $H \in \mathbf{R}^{n \times n}$. Give an algorithm for computing an orthogonal matrix $P$ such that

$$
P^{T} H P=\left[\begin{array}{cc}
\lambda & w^{T} \\
0 & H_{1}
\end{array}\right]
$$

where $H_{1} \in \mathbf{R}^{(n-1) \times(n-1)}$ is upper Hessenberg. Compute $P$ as a product of Givens rotations.
P7.4.8 Suppose $H \in \mathbf{R}^{n \times n}$ has lower bandwidth $p$. Show how to compute $Q \in \mathbf{R}^{n \times n}$, a product of Givens rotations, such that $Q^{T} H Q$ is upper Hessenberg. How many flops are required?
P7.4.9 Show that if $C$ is a companion matrix with distinct eigenvalues $\lambda_{1}, \ldots, \lambda_{n}$, then $V C V^{-1}= \operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{n}\right)$ where

$$
V=\left[\begin{array}{cccc}
1 & \lambda_{1} & \cdots & \lambda_{1}^{n-1} \\
1 & \lambda_{2} & \cdots & \lambda_{2}^{n-1} \\
\vdots & \vdots & \ddots & \vdots \\
1 & \lambda_{n} & \cdots & \lambda_{n}^{n-1}
\end{array}\right]
$$

## Notes and References for §7.4

The real Schur decomposition was originally presented in:
F.D. Murnaghan and A. Wintner (1931). "A Canonical Form for Real Matrices Under Orthogonal Transformations," Proc. Nat. Acad. Sci. 17, 417-420.

A thorough treatment of the reduction to Hessenberg form is given in Wilkinson (AEP, Chap. 6), and Algol procedures appear in:
R.S. Martin and J.H. Wilkinson (1968). "Similarity Reduction of a General Matrix to Hessenberg Form," Numer. Math. 12, 349-368.

Givens rotations can also be used to compute the Hessenberg decomposition, see:
W. Rath (1982). "Fast Givens Rotations for Orthogonal Similarity," Numer. Math. 40, 47-56.

The high-performance computation of the Hessenberg reduction is a major challenge because it is a two-sided factorization, see:
J.J. Dongarra, L. Kaufman, and S. Hammarling (1986). "Squeezing the Most Out of Eigenvalue Solvers on High Performance Computers," Lin. Alg. Applic. 77, 113-136.
J.J. Dongarra, S. Hammarling, and D.C. Sorensen (1989). "Block Reduction of Matrices to Condensed Forms for Eigenvalue Computations," J. ACM 27, 215-227.
M.W. Berry, J.J. Dongarra, and Y. Kim (1995). "A Parallel Algorithm for the Reduction of a Nonsymmetric Matrix to Block Upper Hessenberg Form," Parallel Comput. 21, 1189-1211.
G. Quintana-Orti and R. Van Dc Gcijn (2006). "Improving the Performance of Reduction to Hessenberg Form," ACM Trans. Math. Softw. 32, 180-194.
S. Tomov, R. Nath, and J. Dongarra (2010). "Accelerating the Reduction to Upper Hessenberg, Tridiagonal, and Bidiagonal Forms Through Hybrid GPU-Based Computing," Parallel Comput. 36, 645-654.
L. Karlsson (2011). "Scheduling of Parallel Matrix Computations and Data Layout Conversion for HPC and Multicore Architectures," PhD Thesis, University of Umeå.

Reaching the Hessenberg form via Gauss transforms is discussed in:
P. Businger (1969). "Reducing a Matrix to Hessenberg Form," Math. Comput. 23, 819-821.
G.W. Howell and N. Diaa (2005). "Algorithm 841: BHESS: Gaussian Reduction to a Similar Banded Hessenberg Form," ACM Trans. Math. Softw. 31, 166-185.

Some interesting mathematical properties of the Hessenberg form may be found in:
B.N. Parlett (1967). "Canonical Decomposition of Hessenberg Matrices," Math. Comput. 21, 223227.

Although the Hessenberg decomposition is largely appreciated as a "front end" decomposition for the QR iteration, it is increasingly popular as a cheap alternative to the more expensive Schur decomposition in certain problems. For a sampling of applications where it has proven to be very useful, consult:
W. Enright (1979). "On the Efficient and Reliable Numerical Solution of Large Linear Systems of O.D.E.'s," IEEE Trans. Autom. Contr. AC-24, 905-908.
G.H. Golub, S. Nash and C. Van Loan (1979). "A Hessenberg-Schur Method for the Problem $A X+ X B=C$," IEEE Trans. Autom. Contr. AC-24, 909-913.
A. Laub (1981). "Efficient Multivariable Frequency Response Computations," IEEE Trans. Autom. Contr. AC-26, 407-408.
C.C. Paige (1981). "Properties of Numerical Algorithms Related to Computing Controllability," IEEE Trans. Auto. Contr. AC-26, 130-138.
G. Miminis and C.C. Paige (1982). "An Algorithm for Pole Assignment of Time Invariant Linear Systems," Int. J. Contr. 35, 341-354.
C. Van Loan (1982). "Using the Hessenberg Decomposition in Control Theory," in Algorithms and Theory in Filtering and Control , D.C. Sorensen and R.J. Wets (eds.), Mathematical Programming Study No. 18, North Holland, Amsterdam, 102-111.
C.D. Martin and C.F. Van Loan (2006). "Solving Real Linear Systems with the Complex Schur Decomposition," SIAM J. Matrix Anal. Applic. 29, 177-183.

The advisability of posing polynomial root problems as companion matrix eigenvalue problem is discussed in:
A. Edelman and H. Murakami (1995). "Polynomial Roots from Companion Matrix Eigenvalues," Math. Comput. 64, 763-776.

### 7.5 The Practical QR Algorithm

We return to the Hessenberg QR iteration, which we write as follows:

$$
\begin{aligned}
& H=U_{0}^{T} A U_{0} \quad \text { (Hessenberg reduction) } \\
& \text { for } k=1,2, \ldots \\
& \qquad \begin{aligned}
H & =U R \\
H & =R U
\end{aligned} \\
& \text { end }
\end{aligned}
$$

Our aim in this section is to describe how the $H$ 's converge to upper quasi-triangular form and to show how the convergence rate can be accelerated by incorporating shifts.

### 7.5.1 Deflation

Without loss of generality we may assume that each Hessenberg matrix $H$ in (7.5.1) is unreduced. If not, then at some stage we have

$$
H=\left[\begin{array}{cc}
H_{11} & H_{12} \\
0 & H_{22} \\
p & n-p
\end{array}\right]_{n-p}^{p}
$$

where $1 \leq p<n$ and the problem decouples into two smaller problems involving $H_{11}$ and $H_{22}$. The term deflation is also used in this context, usually when $p=n-1$ or $n-2$.

In practice, decoupling occurs whenever a subdiagonal entry in $H$ is suitably small. For example, if

$$
\left|h_{p+1, p}\right| \leq c \mathbf{u}\left(\left|h_{p p}\right|+\left|h_{p+1, p+1}\right|\right)
$$

for a small constant $c$, then $h_{p+1, p}$ can justifiably be set to zero because rounding errors of order $\mathbf{u}\|H\|$ are typically present throughout the matrix anyway.

### 7.5.2 The Shifted QR Iteration

Let $\mu \in \mathbb{R}$ and consider the iteration:

$$
\begin{aligned}
& H=U_{0}^{T} A U_{0} \quad \text { (Hessenberg reduction) } \\
& \text { for } k=1,2, \ldots \\
& \qquad \begin{array}{l}
\text { Determine a scalar } \mu \text {. } \\
H-\mu I=U R \quad \text { (QR factorization) } \\
H=R U+\mu I
\end{array} \\
& \text { end }
\end{aligned}
$$

The scalar $\mu$ is referred to as a shift . Each matrix $H$ generated in (7.5.3) is similar to $A$, since

$$
R U+\mu I=U^{T}(U R+\mu I) U=U^{T} H U
$$

If we order the eigenvalues $\lambda_{i}$ of $A$ so that

$$
\left|\lambda_{1}-\mu\right| \geq \cdots \geq\left|\lambda_{n}-\mu\right|
$$

and $\mu$ is fixed from iteration to iteration, then the theory of §7.3 says that the $p$ th subdiagonal entry in $H$ converges to zero with rate

$$
\left|\frac{\lambda_{p+1}-\mu}{\lambda_{p}-\mu}\right|^{k}
$$

Of course, if $\lambda_{p}=\lambda_{p+1}$, then there is no convergence at all. But if, for example, $\mu$ is much closer to $\lambda_{n}$ than to the other eigenvalues, then the zeroing of the $(n, n-1)$ entry is rapid. In the extreme case we have the following:

Theorem 7.5.1. Let $\mu$ be an eigenvalue of an $n$-by- $n$ unreduced Hessenberg matrix H. If

$$
\tilde{H}=R U+\mu I,
$$

where $H-\mu I=U R$ is the $Q R$ factorization of $H-\mu I$, then $\tilde{h}_{n, n-1}=0$ and $\tilde{h}_{n n}=\mu$.
Proof. Since $H$ is an unreduced Hessenberg matrix the first $n-1$ columns of $H-\mu I$ are independent, regardless of $\mu$. Thus, if $U R=(H-\mu I)$ is the QR factorization then $r_{i i} \neq 0$ for $i=1: n-1$. But if $H-\mu I$ is singular, then $r_{11} \cdots r_{n n}=0$. Thus, $r_{n n}=0$ and $\tilde{H}(n,:)=[0, \ldots, 0, \mu]$.

The theorem says that if we shift by an exact eigenvalue, then in exact arithmetic deflation occurs in one step.

### 7.5.3 The Single-Shift Strategy

Now let us consider varying $\mu$ from iteration to iteration incorporating new information about $\lambda(A)$ as the subdiagonal entries converge to zero. A good heuristic is to regard $h_{n n}$ as the best approximate eigenvalue along the diagonal. If we shift by this quantity during each iteration, we obtain the single-shift $Q R$ iteration:

$$
\begin{aligned}
& \text { for } k=1,2, \ldots \\
& \qquad \begin{aligned}
\mu & =H(n, n) \\
H & -\mu I=U R \quad \text { (QR factorization) } \\
H & =R U+\mu I
\end{aligned} \\
& \text { end }
\end{aligned}
$$

If the $(n, n-1)$ entry converges to zero, it is likely to do so at a quadratic rate. To see this, we borrow an example from Stewart (IMC, p. 366). Suppose $H$ is an unreduced upper Hessenberg matrix of the form

$$
H=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \epsilon & h_{n n}
\end{array}\right]
$$

and that we perform one step of the single-shift QR algorithm, i.e.,

$$
\begin{aligned}
& U R=H-h_{n n} \\
& \tilde{H}=R U+h_{n n} I .
\end{aligned}
$$

After $n-2$ steps in the orthogonal reduction of $H-h_{n n} I$ to upper triangular form we obtain a matrix with the following structure:

$$
H=\left[\begin{array}{lllll}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & a & b \\
0 & 0 & 0 & \epsilon & 0
\end{array}\right]
$$

It is not hard to show that

$$
\tilde{h}_{n, n-1}=-\frac{\epsilon^{2} b}{a^{2}+\epsilon^{2}}
$$

If we assume that $\epsilon \ll a$, then it is clear that the new $(n, n-1)$ entry has order $\epsilon^{2}$, precisely what we would expect of a quadratically converging algorithm.

### 7.5.4 The Double-Shift Strategy

Unfortunately, difficulties with (7.5.4) can be expected if at some stage the eigenvalues $a_{1}$ and $a_{2}$ of

$$
G=\left[\begin{array}{cc}
h_{m m} & h_{m n} \\
h_{n m} & h_{n n}
\end{array}\right], \quad m=n-1
$$

are complex for then $h_{n n}$ would tend to be a poor approximate eigenvalue.
A way around this difficulty is to perform two single-shift QR steps in succession using $a_{1}$ and $a_{2}$ as shifts:

$$
\begin{aligned}
H-a_{1} I & =U_{1} R_{1} \\
H_{1} & =R_{1} U_{1}+a_{1} I \\
H_{1}-a_{2} I & =U_{2} R_{2} \\
H_{2} & =R_{2} U_{2}+a_{2} I
\end{aligned}
$$

These equations can be manipulated to show that

$$
\left(U_{1} U_{2}\right)\left(R_{2} R_{1}\right)=M
$$

where $M$ is defined by

$$
M=\left(H-a_{1} I\right)\left(H-a_{2} I\right)
$$

Note that $M$ is a real matrix even if $G$ 's eigenvalues are complex since

$$
M=H^{2}-s H+t I
$$

where

$$
s=a_{1}+a_{2}=h_{m m}+h_{n n}=\operatorname{tr}(G) \in \mathbb{R}
$$

and

$$
t=a_{1} a_{2}=h_{m m} h_{n n}-h_{m n} h_{n m}=\operatorname{det}(G) \in \mathbb{R}
$$

Thus, (7.5.7) is the QR factorization of a real matrix and we may choose $U_{1}$ and $U_{2}$ so that $Z=U_{1} U_{2}$ is real orthogonal. It then follows that

$$
H_{2}=U_{2}^{H} H_{1} U_{2}=U_{2}^{H}\left(U_{1}^{H} H U_{1}\right) U_{2}=\left(U_{1} U_{2}\right)^{H} H\left(U_{1} U_{2}\right)=Z^{T} H Z
$$

is real.
Unfortunately, roundoff error almost always prevents an exact return to the real field. A real $H_{2}$ could be guaranteed if we

- explicitly form the real matrix $M=H^{2}-s H+t I$,
- compute the real QR factorization $M=Z R$, and
- set $H_{2}=Z^{T} H Z$.

But since the first of these steps requires $O\left(n^{3}\right)$ flops, this is not a practical course of action.

### 7.5.5 The Double-Implicit-Shift Strategy

Fortunately, it turns out that we can implement the double-shift step with $O\left(n^{2}\right)$ flops by appealing to the implicit Q theorem of §7.4.5. In particular we can effect the transition from $H$ to $H_{2}$ in $O\left(n^{2}\right)$ flops if we

- compute $M e_{1}$, the first column of $M$;
- determine a Householder matrix $P_{0}$ such that $P_{0}\left(M e_{1}\right)$ is a multiple of $e_{1}$;
- compute Householder matrices $P_{1}, \ldots, P_{n-2}$ such that if

$$
Z_{1}=P_{0} P_{1} \cdots P_{n-2}
$$

then $Z_{1}^{T} H Z_{1}$ is upper Hessenberg and the first columns of $Z$ and $Z_{1}$ are the same.
Under these circumstances, the implicit Q theorem permits us to conclude that, if $Z^{T} H Z$ and $Z_{1}^{T} H Z_{1}$ are both unreduced upper Hessenberg matrices, then they are essentially equal. Note that if these Hessenberg matrices are not unreduced, then we can effect a decoupling and proceed with smaller unreduced subproblems.

Let us work out the details. Observe first that $P_{0}$ can be determined in $O(1)$ flops since $M e_{1}=[x, y, z, 0, \ldots, 0]^{T}$ where

$$
\begin{aligned}
& x=h_{11}^{2}+h_{12} h_{21}-s h_{11}+t \\
& y=h_{21}\left(h_{11}+h_{22}-s\right) \\
& z=h_{21} h_{32}
\end{aligned}
$$

Since a similarity transformation with $P_{0}$ only changes rows and columns 1, 2, and 3, we see that

$$
P_{0} H P_{0}=\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times
\end{array}\right] .
$$

Now the mission of the Householder matrices $P_{1}, \ldots, P_{n-2}$ is to restore this matrix to upper Hessenberg form. The calculation proceeds as follows:

$$
\begin{gathered}
{\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times
\end{array}\right] \xrightarrow{P_{1}}\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times
\end{array}\right] \xrightarrow{P_{2}}} \\
{\left[\begin{array}{llllll}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times
\end{array}\right] \xrightarrow{P_{3}}\left[\begin{array}{llllll}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times
\end{array}\right] \xrightarrow{P_{3}}\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times
\end{array}\right] .}
\end{gathered}
$$

Each $P_{k}$ is the identity with a 3-by-3 or 2-by-2 Householder somewhere along its diagonal, e.g.,

$$
\begin{aligned}
& P_{1}=\left[\begin{array}{cccccc}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & \times & \times & \times & 0 & 0 \\
0 & \times & \times & \times & 0 & 0 \\
0 & \times & \times & \times & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 1
\end{array}\right], \quad P_{2}=\left[\begin{array}{cccccc}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & \times & \times & \times & 0 \\
0 & 0 & \times & \times & \times & 0 \\
0 & 0 & \times & \times & \times & 0 \\
0 & 0 & 0 & 0 & 0 & 1
\end{array}\right], \\
& P_{3}=\left[\begin{array}{llllll}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times
\end{array}\right], \quad P_{4}=\left[\begin{array}{cccccc}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & \times & 0 & 0 \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times
\end{array}\right] .
\end{aligned}
$$

The applicability of Theorem 7.4.3 (the implicit Q theorem) follows from the observation that $P_{k} e_{1}=e_{1}$ for $k=1: n-2$ and that $P_{0}$ and $Z$ have the same first column. Hence, $Z_{1} e_{1}=Z e_{1}$, and we can assert that $Z_{1}$ essentially equals $Z$ provided that the upper Hessenberg matrices $Z^{T} H Z$ and $Z_{1}^{T} H Z_{1}$ are each unreduced.

The implicit determination of $H_{2}$ from $H$ outlined above was first described by Francis (1961) and we refer to it as a Francis $Q R$ step. The complete Francis step is summarized as follows:

Algorithm 7.5.1 (Francis QR step) Given the unreduced upper Hessenberg matrix $H \in \mathbb{R}^{n \times n}$ whose trailing 2-by-2 principal submatrix has eigenvalues $a_{1}$ and $a_{2}$, this algorithm overwrites $H$ with $Z^{T} H Z$, where $Z$ is a product of Householder matrices and $Z^{T}\left(H-a_{1} I\right)\left(H-a_{2} I\right)$ is upper triangular.

```
$m=n-1$
\{Compute first column of $\left(H-a_{1} I\right)\left(H-a_{2} I\right)$ \}
$s=H(m, m)+H(n, n)$
$t=H(m, m) \cdot H(n, n)-H(m, n) \cdot H(n, m)$
$x=H(1,1) \cdot H(1,1)+H(1,2) \cdot H(2,1)-s \cdot H(1,1)+t$
$y=H(2,1) \cdot(H(1,1)+H(2,2)-s)$
$z=H(2,1) \cdot H(3,2)$
for $k=0: n-3$
    $[v, \beta]=$ house $\left(\left[\begin{array}{lll}x & y & z\end{array}\right]^{T}\right)$
    $q=\max \{1, k\}$.
    $H(k+1: k+3, q: n)=\left(I-\beta v v^{T}\right) \cdot H(k+1: k+3, q: n)$
    $r=\min \{k+4, n\}$
    $H(1: r, k+1: k+3)=H(1: r, k+1: k+3) \cdot\left(I-\beta v v^{T}\right)$
    $x=H(k+2, k+1)$
    $y=H(k+3, k+1)$
    if $k<n-3$
        $z=H(k+4, k+1)$
    end
end
$[v, \beta]=$ house $\left([x y]^{T}\right)$
$H(n-1: n, n-2: n)=\left(I-\beta v v^{T}\right) \cdot H(n-1: n, n-2: n)$
$H(1: n, n-1: n)=H(1: n, n-1: n) \cdot\left(I-\beta v v^{T}\right)$
```

This algorithm requires $10 n^{2}$ flops. If $Z$ is accumulated into a given orthogonal matrix, an additional $10 n^{2}$ flops are necessary.

### 7.5.6 The Overall Process

Reduction of $A$ to Hessenberg form using Algorithm 7.4.2 and then iteration with Algorithm 7.5.1 to produce the real Schur form is the standard means by which the dense unsymmetric eigenproblem is solved. During the iteration it is necessary to monitor the subdiagonal elements in $H$ in order to spot any possible decoupling. How this is done is illustrated in the following algorithm:

Algorithm 7.5.2 (QR Algorithm) Given $A \in \mathbb{R}^{n \times n}$ and a tolerance tol greater than the unit roundoff, this algorithm computes the real Schur canonical form $Q^{T} A Q=T$. If $Q$ and $T$ are desired, then $T$ is stored in $H$. If only the eigenvalues are desired, then diagonal blocks in $T$ are stored in the corresponding positions in $H$.

Use Algorithm 7.4.2 to compute the Hessenberg reduction

$$
H=U_{0}^{T} A U_{0} \text { where } U_{0}=P_{1} \cdots P_{n-2} .
$$

If $Q$ is desired form $Q=P_{1} \cdots P_{n-2}$. (See §5.1.6.)
until $q=n$
Set to zero all subdiagonal elements that satisfy:

$$
\left|h_{i, i-1}\right| \leq \text { tol } \cdot\left(\left|h_{i i}\right|+\left|h_{i-1, i-1}\right|\right) .
$$

Find the largest nonnegative $q$ and the smallest non-negative $p$ such that

$$
H=\begin{array}{ccc}
{\left[\begin{array}{ccc}
H_{11} & H_{12} & H_{13} \\
0 \\
0
\end{array}\right.} & \begin{array}{c}
H_{22} \\
0
\end{array} & \begin{array}{c}
H_{23} \\
H_{33}
\end{array} \\
p & n-p-q & q
\end{array} \begin{gathered}
p \\
n-p-q \\
q
\end{gathered}
$$

where $H_{33}$ is upper quasi-triangular and $H_{22}$ is unreduced.
if $q<n$
Perform a Francis QR step on $H_{22}: H_{22}=Z^{T} H_{22} Z$.
if $Q$ is required
$Q=Q \cdot \operatorname{diag}\left(I_{p}, Z, I_{q}\right)$
$H_{12}=H_{12} Z$
$H_{23}=Z^{T} H_{23}$
end.
end
end
Upper triangularize all 2-by-2 diagonal blocks in $H$ that have real eigenvalues and accumulate the transformations (if necessary).

This algorithm requires $25 n^{3}$ flops if $Q$ and $T$ are computed. If only the eigenvalues are desired, then $10 n^{3}$ flops are necessary. These flops counts are very approximate and are based on the empirical observation that on average only two Francis iterations are required before the lower 1-by-1 or 2-by-2 decouples.

The roundoff properties of the QR algorithm are what one would expect of any orthogonal matrix technique. The computed real Schur form $\hat{T}$ is orthogonally similar to a matrix near to $A$, i.e.,

$$
Q^{T}(A+E) Q=\hat{T}
$$

where $Q^{T} Q=I$ and $\|E\|_{2} \approx \mathbf{u}\|A\|_{2}$. The computed $\hat{Q}$ is almost orthogonal in the sense that $\hat{Q}^{T} \hat{Q}=I+F$ where $\|F\|_{2} \approx \mathbf{u}$.

The order of the eigenvalues along $\hat{T}$ is somewhat arbitrary. But as we discuss in §7.6, any ordering can be achieved by using a simple procedure for swapping two adjacent diagonal entries.

### 7.5.7 Balancing

Finally, we mention that if the elements of $A$ have widely varying magnitudes, then $A$ should be balanced before applying the QR algorithm. This is an $O\left(n^{2}\right)$ calculation in which a diagonal matrix $D$ is computed so that if

$$
D^{-1} A D=\left[c_{1}|\cdots| c_{n}\right]=\left[\begin{array}{c}
r_{1}^{T} \\
\vdots \\
r_{n}^{T}
\end{array}\right]
$$

then $\left\|r_{i}\right\|_{\infty} \approx\left\|c_{i}\right\|_{\infty}$ for $i=1: n$. The diagonal matrix $D$ is chosen to have the form

$$
D=\operatorname{diag}\left(\beta^{i_{1}}, \ldots, \beta^{i_{n}}\right)
$$

where $\beta$ is the floating point base. Note that $D^{-1} A D$ can be calculated without roundoff. When $A$ is balanced, the computed eigenvalues are usually more accurate although there are exceptions. See Parlett and Reinsch (1969) and Watkins(2006).

## Problems

P7.5.1 Show that if $\bar{H}=Q^{T} H Q$ is obtained by performing a single-shift QR step with

$$
H=\left[\begin{array}{ll}
w & x \\
y & z
\end{array}\right]
$$

then $\left|\bar{h}_{21}\right| \leq\left|y^{2} x\right| /\left[(w-z)^{2}+y^{2}\right]$.
P7.5.2 Given $A \in \mathbf{R}^{2 \times 2}$, show how to compute a diagonal $D \in \mathbf{R}^{2 \times 2}$ so that $\left\|D^{-1} A D\right\|_{F}$ is minimized.
P7.5.3 Explain how the single-shift QR step $H-\mu I=U R, \tilde{H}=R U+\mu I$ can be carried out implicitly. That is, show how the transition from $H$ to $\dot{H}$ can be carried out without subtracting the shift $\mu$ from the diagonal of $H$.
P7.5.4 Suppose $H$ is upper Hessenberg and that we compute the factorization $P H=L U$ via Gaussian elimination with partial pivoting. (See Algorithm 4.3.4.) Show that $H_{1}=U\left(P^{T} L\right)$ is upper Hessenberg and similar to $H$. (This is the basis of the modified LR algorithm.)
P7.5.5 Show that if $H=H_{0}$ is given and we generate the matrices $H_{k}$ via $H_{k}-\mu_{k} I=U_{k} R_{k}, H_{k+1} =R_{k} U_{k}+\mu_{k} I$, then $\left(U_{1} \cdots U_{j}\right)\left(R_{j} \cdots R_{1}\right)=\left(H-\mu_{1} I\right) \cdots\left(H-\mu_{j} I\right)$.

## Notes and References for §7.5

Historically important papers associated with the QR iteration include:
H. Rutishauser (1958). "Solution of Eigenvalue Problems with the LR Transformation," Nat. Bur. Stand. App. Math. Ser. 49, 47-81.
J.G.F. Francis (1961). "The QR Transformation: A Unitary Analogue to the LR Transformation, Parts I and II" Comput. J. 4, 265-72, 332-345.
V.N. Kublanovskaya (1961). "On Some Algorithms for the Solution of the Complete Eigenvalue Problem," Vychisl. Mat. Mat. Fiz 1(4), 555-570.
R.S. Martin and J.H. Wilkinson (1968). "The Modified LR Algorithm for Complex Hessenberg Matrices," Numer. Math. 12, 369-376.
R.S. Martin, G. Peters, and J.H. Wilkinson (1970). "The QR Algorithm for Real Hessenberg Matrices," Numer. Math. 14, 219-231.

For a general insight, we recommend:
D.S. Watkins (1982). "Understanding the QR Algorithm," SIAM Review 24, 427-440.
D.S. Watkins (1993). "Some Perspectives on the Eigenvalue Problem," SIAM Review 35, 430-471.
D.S. Watkins (2008). "The QR Algorithm Revisited," SIAM Review 50, 133-145.
D.S. Watkins (2011). "Francis's Algorithm," Amer. Math. Monthly 118, 387-403.

Papers concerned with the convergence of the method, shifting, deflation, and related matters include:
P.A. Businger (1971). "Numerically Stable Deflation of Hessenberg and Symmetric Tridiagonal Matrices, BIT 11, 262-270.
D.S. Watkins and L. Elsner (1991). "Chasing Algorithms for the Eigenvalue Problem," SIAM J. Matrix Anal. Applic. 12, 374-384.
D.S. Watkins and L. Elsner (1991). "Convergence of Algorithms of Decomposition Type for the Eigenvalue Problem," Lin. Alg. Applic. 143, 19-47.
J. Erxiong (1992). "A Note on the Double-Shift QL Algorithm," Lin. Alg. Applic. 171, 121-132.
A.A. Dubrulle and G.H. Golub (1994). "A Multishift QR Iteration Without Computation of the Shifts," Numer. Algorithms 7, 173-181.
D.S. Watkins (1996). "Forward Stability and Transmission of Shifts in the QR Algorithm," SIAM J. Matrix Anal. Applic. 16, 469-487.
D.S. Watkins (1996). "The Transmission of Shifts and Shift Blurring in the QR algorithm," Lin. Alg. Applic. 241-3, 877-896.
D.S. Watkins (1998). "Bulge Exchanges in Algorithms of QR Type," SIAM J. Matrix Anal. Applic. 19, 1074-1096.
R. Vandebril (2011). "Chasing Bulges or Rotations? A Metamorphosis of the QR-Algorithm" SIAM. J. Matrix Anal. Applic. 32, 217-247.

Aspects of the balancing problem are discussed in:
E.E. Osborne (1960). "On Preconditioning of Matrices," J. ACM 7, 338-345.
B.N. Parlett and C. Reinsch (1969). "Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors," Numer. Math. 13, 292-304.
D.S. Watkins (2006). "A Case Where Balancing is Harmful," ETNA 23, 1-4.

Versions of the algorithm that are suitable for companion matrices are discussed in:
D.A. Bini, F. Daddi, and L. Gemignani (2004). "On the Shifted QR iteration Applied to Companion Matrices," ETNA 18, 137-152.
M. Van Barel, R. Vandebril, P. Van Dooren, and K. Frederix (2010). "Implicit Double Shift QRAlgorithm for Companion Matrices," Numer. Math. 116, 177-212.

Papers that are concerned with the high-performance implementation of the QR iteration include:
Z. Bai and J.W. Demmel (1989). "On a Block Implementation of Hessenberg Multishift QR Iteration," Int. J. High Speed Comput. 1, 97-112.
R.A. Van De Geijn (1993). "Deferred Shifting Schemes for Parallel QR Methods," SIAM J. Matrix Anal. Applic. 14, 180-194.
D.S. Watkins (1994). "Shifting Strategies for the Parallel QR Algorithm," SIAM J. Sci. Comput. 15, 953-958.
G. Henry and R. van de Geijn (1996). "Parallelizing the QR Algorithm for the Unsymmetric Algebraic Eigenvalue Problem: Myths and Reality," SIAM J. Sci. Comput. 17, 870-883.
Z. Bai, J. Demmel, J. Dongarra, A. Petitet, H. Robinson, and K. Stanley (1997). "The Spectral Decomposition of Nonsymmetric Matrices on Distributed Memory Parallel Computers," SIAM J. Sci. Comput. 18, 1446-1461.
G. Henry, D.S. Watkins, and J. Dongarra (2002). "A Parallel Implementation of the Nonsymmetric QR Algorithm for Distributed Memory Architectures," SIAM J. Sci. Comput. 24, 284-311.
K. Braman, R. Byers, and R. Mathias (2002). "The Multishift QR Algorithm. Part I: Maintaining Well-Focused Shifts and Level 3 Performance," SIAM J. Matrix Anal. Applic. 23, 929-947.
K. Braman, R. Byers, and R. Mathias (2002). "The Multishift QR Algorithm. Part II: Aggressive Early Deflation," SIAM J. Matrix Anal. Applic. 23, 948-973.
M.R. Fahey (2003). "Algorithm 826: A Parallel Eigenvalue Routine for Complex Hessenberg Matrices," ACM Trans. Math. Softw. 29, 326-336.
D. Kressner (2005). "On the Use of Larger Bulges in the QR Algorithm," ETNA 20, 50-63.
D. Kressner (2008). "The Effect of Aggressive Early Deflation on the Convergence of the QR Algorithm," SIAM J. Matrix Anal. Applic. 30, 805-821.

### 7.6 Invariant Subspace Computations

Several important invariant subspace problems can be solved once the real Schur decomposition $Q^{T} A Q=T$ has been computed. In this section we discuss how to

- compute the eigenvectors associated with some subset of $\lambda(A)$,
- compute an orthonormal basis for a given invariant subspace,
- block-diagonalize $A$ using well-conditioned similarity transformations,
- compute a basis of eigenvectors regardless of their condition, and
- compute an approximate Jordan canonical form of $A$.

Eigenvector/invariant subspace computation for sparse matrices is discussed in §7.3.1 and §7.3.2 as well as portions of Chapters 8 and 10.

### 7.6.1 Selected Eigenvectors via Inverse Iteration

Let $q^{(0)} \in \mathbb{R}^{n}$ be a given unit 2 -norm vector and assume that $A-\mu I \in \mathbb{R}^{n \times n}$ is nonsingular. The following is referred to as inverse iteration:
for $k=1,2, \ldots$

$$
\begin{aligned}
& \text { Solve }(A-\mu I) z^{(k)}=q^{(k-1)} \text {. } \\
& q^{(k)}=z^{(k)} /\left\|z^{(k)}\right\|_{2} \\
& \lambda^{(k)}=q^{(k)^{T}} A q^{(k)}
\end{aligned}
$$

Inverse iteration is just the power method applied to $(A-\mu I)^{-1}$.
To analyze the behavior of (7.6.1), assume that $A$ has a basis of eigenvectors $\left\{x_{1}, \ldots, x_{n}\right\}$ and that $A x_{i}=\lambda_{i} x_{i}$ for $i=1: n$. If

$$
q^{(0)}=\sum_{i=1}^{n} \beta_{i} x_{i}
$$

then $q^{(k)}$ is a unit vector in the direction of

$$
(A-\mu I)^{-k} q^{(0)}=\sum_{i=1}^{n} \frac{\beta_{i}}{\left(\lambda_{i}-\mu\right)^{k}} x_{i} .
$$

Clearly, if $\mu$ is much closer to an eigenvalue $\lambda_{j}$ than to the other eigenvalues, then $q^{(k)}$ is rich in the direction of $x_{j}$ provided $\beta_{j} \neq 0$.

A sample stopping criterion for (7.6.1) might be to quit as soon as the residual

$$
r^{(k)}=(A-\mu I) q^{(k)}
$$

satisfies

$$
\left\|r^{(k)}\right\|_{\infty} \leq c \mathbf{u}\|A\|_{\infty}
$$

where $c$ is a constant of order unity. Since

$$
\left(A+E_{k}\right) q^{(k)}=\mu q^{(k)}
$$

with $E_{k}=-r^{(k)} q^{(k)^{T}}$, it follows that (7.6.2) forces $\mu$ and $q^{(k)}$ to be an exact eigenpair for a nearby matrix.

Inverse iteration can be used in conjunction with Hessenberg reduction and the QR algorithm as follows:

Step 1. Compute the Hessenberg decomposition $U_{0}^{T} A U_{0}=H$.
Step 2. Apply the double-implicit-shift Francis iteration to $H$ without accumulating transformations.

Step 3. For each computed eigenvalue $\lambda$ whose corresponding eigenvector $x$ is sought, apply (7.6.1) with $A=H$ and $\mu=\lambda$ to produce a vector $z$ such that $H z \approx \mu z$.

Step 4. Set $x=U_{0} z$.
Inverse iteration with $H$ is very economical because we do not have to accumulate transformations during the double Francis iteration. Moreover, we can factor matrices of the form $H-\lambda I$ in $O\left(n^{2}\right)$ flops, and (3) only one iteration is typically required to produce an adequate approximate eigenvector.

This last point is perhaps the most interesting aspect of inverse iteration and requires some justification since $\lambda$ can be comparatively inaccurate if it is ill-conditioned. Assume for simplicity that $\lambda$ is real and let

$$
H-\lambda I=\sum_{i=1}^{n} \sigma_{i} u_{i} v_{i}^{T}=U \Sigma V^{T}
$$

be the SVD of $H-\lambda I$. From what we said about the roundoff properties of the QR algorithm in §7.5.6, there exists a matrix $E \in \mathbb{R}^{n \times n}$ such that $H+E-\lambda I$ is singular and $\|E\|_{2} \approx \mathbf{u}\|H\|_{2}$. It follows that $\sigma_{n} \approx \mathbf{u} \sigma_{1}$ and

$$
\left\|(H-\hat{\lambda} I) v_{n}\right\|_{2} \approx \mathbf{u} \sigma_{1},
$$

i.e., $v_{n}$ is a good approximate eigenvector. Clearly if the starting vector $q^{(0)}$ has the expansion

$$
q^{(0)}=\sum_{i=1}^{n} \gamma_{i} u_{i}
$$

then

$$
z^{(1)}=\sum_{i=1}^{n} \frac{\gamma_{i}}{\sigma_{i}} v_{i}
$$

is "rich" in the direction $v_{n}$. Note that if $s(\lambda) \approx\left|u_{n}^{T} v_{n}\right|$ is small, then $z^{(1)}$ is rather deficient in the direction $u_{n}$. This explains (heuristically) why another step of inverse iteration is not likely to produce an improved eigenvector approximate, especially if $\lambda$ is ill-conditioned. For more details, see Peters and Wilkinson (1979).

### 7.6.2 Ordering Eigenvalues in the Real Schur Form

Recall that the real Schur decomposition provides information about invariant subspaces. If

$$
Q^{T} A Q=T=\left[\begin{array}{cc}
T_{11} & T_{12} \\
0 & T_{22} \\
p & q
\end{array}\right]_{q}^{p}
$$

and

$$
\lambda\left(T_{11}\right) \cap \lambda\left(T_{22}\right)=\emptyset
$$

then the first $p$ columns of $Q$ span the unique invariant subspace associated with $\lambda\left(T_{11}\right)$. (See §7.1.4.) Unfortunately, the Francis iteration supplies us with a real Schur decomposition $Q_{F}^{T} A Q_{F}=T_{F}$ in which the eigenvalues appear somewhat randomly along the diagonal of $T_{F}$. This poses a problem if we want an orthonormal basis for an invariant subspace whose associated eigenvalues are not at the top of $T_{F}$ 's diagonal. Clearly, we need a method for computing an orthogonal matrix $Q_{D}$ such that $Q_{D}^{T} T_{F} Q_{D}$ is upper quasi-triangular with appropriate eigenvalue ordering.

A look at the 2-by-2 case suggests how this can be accomplished. Suppose

$$
Q_{F}^{T} A Q_{F}=T_{F}=\left[\begin{array}{rr}
\lambda_{1} & t_{12} \\
0 & \lambda_{2}
\end{array}\right], \quad \lambda_{1} \neq \lambda_{2}
$$

and that we wish to reverse the order of the eigenvalues. Note that

$$
T_{F} x=\lambda_{2} x
$$

where

$$
x=\left[\begin{array}{c}
t_{12} \\
\lambda_{2}-\lambda_{1}
\end{array}\right]
$$

Let $Q_{D}$ be a Givens rotation such that the second component of $Q_{D}^{T} x$ is zero. If

$$
Q=Q_{F} Q_{D}
$$

then

$$
\left(Q^{T} A Q\right) e_{1}=Q_{D}^{T} T_{F}\left(Q_{D} e_{1}\right)=\lambda_{2} Q_{D}^{T}\left(Q_{D} e_{1}\right)=\lambda_{2} e_{1}
$$

The matrices $A$ and $Q^{T} A Q$ have the same Frobenius norm and so it follows that the latter must have the following form:

$$
Q^{T} A Q=\left[\begin{array}{rr}
\lambda_{2} & \pm t_{12} \\
0 & \lambda_{1}
\end{array}\right]
$$

The swapping gets a little more complicated if $T$ has 2-by-2 blocks along its diagonal. See Ruhe (1970) and Stewart (1976) for details.

By systematically interchanging adjacent pairs of eigenvalues (or 2-by-2 blocks), we can move any subset of $\lambda(A)$ to the top of $T$ 's diagonal. Here is the overall procedure for the case when there are no 2 -by- 2 bumps:

Algorithm 7.6.1 Given an orthogonal matrix $Q \in \mathbb{R}^{n \times n}$, an upper triangular matrix $T=Q^{T} A Q$, and a subset $\Delta=\left\{\lambda_{1}, \ldots, \lambda_{p}\right\}$ of $\lambda(A)$, the following algorithm computes an orthogonal matrix $Q_{D}$ such that $Q_{D}^{T} T Q_{D}=S$ is upper triangular and $\left\{s_{11}, \ldots, s_{p p}\right\} =\Delta$. The matrices $Q$ and $T$ are overwritten by $Q Q_{D}$ and $S$, respectively.

```
while $\left\{t_{11}, \ldots, t_{p p}\right\} \neq \Delta$
    for $k=1: n-1$
        if $t_{k k} \notin \Delta$ and $t_{k+1, k+1} \in \Delta$
            $[c, s]=\operatorname{givens}(T(k, k+1), T(k+1, k+1)-T(k, k))$
            $T(k: k+1, k: n)=\left[\begin{array}{rr}c & s \\ -s & c\end{array}\right]^{T} T(k: k+1, k: n)$
            $T(1: k+1, k: k+1)=T(1: k+1, k: k+1)\left[\begin{array}{rl}c & s \\ -s & c\end{array}\right]$
            $Q(1: n, k: k+1)=Q(1: n, k: k+1)\left[\begin{array}{rr}c & s \\ -s & c\end{array}\right]$
        end
    end
end
```

This algorithm requires $k(12 n)$ flops, where $k$ is the total number of required swaps. The integer $k$ is never greater than $(n-p) p$.

Computation of invariant subspaces by manipulating the real Schur decomposition is extremely stable. If $\hat{Q}=\left[\hat{q}_{1}|\cdots| \hat{q}_{n}\right]$ denotes the computed orthogonal matrix $Q$, then $\left\|\hat{Q}^{T} \hat{Q}-I\right\|_{2} \approx \mathbf{u}$ and there exists a matrix $E$ satisfying $\|E\|_{2} \approx \mathbf{u}\|A\|_{2}$ such that $(A+E) \hat{q}_{i} \in \operatorname{span}\left\{\hat{q}_{1}, \ldots, \hat{q}_{p}\right\}$ for $i=1: p$.

### 7.6.3 Block Diagonalization

Let

$$
T=\begin{array}{cccc}
{\left[\begin{array}{cccc}
T_{11} & T_{12} & \cdots & T_{1 q} \\
0 & T_{22} & \cdots & T_{2 q} \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & T_{q q} \\
n_{1} & n_{2} & & n_{q}
\end{array}\right] \begin{array}{c}
n_{1} \\
n_{2} \\
n_{q}
\end{array} \text { }}
\end{array}
$$

be a partitioning of some real Schur canonical form $Q^{T} A Q=T \in \mathbb{R}^{n \times n}$ such that $\lambda\left(T_{11}\right), \ldots, \lambda\left(T_{q q}\right)$ are disjoint. By Theorem 7.1.6 there exists a matrix $Y$ such that

$$
Y^{-1} T Y=\operatorname{diag}\left(T_{11}, \ldots, T_{q q}\right) .
$$

A practical procedure for determining $Y$ is now given together with an analysis of $Y$ 's sensitivity as a function of the above partitioning.

Partition $I_{n}=\left[E_{1}|\cdots| E_{q}\right]$ conformably with $T$ and define the matrix $Y_{i j} \in \mathbb{R}^{n \times n}$ as follows:

$$
Y_{i j}=I_{n}+E_{i} Z_{i j} E_{j}^{T}, \quad i<j, Z_{i j} \in \mathbb{R}^{n_{i} \times n_{j}} .
$$

In other words, $Y_{i j}$ looks just like the identity except that $Z_{i j}$ occupies the $(i, j)$ block position. It follows that if $Y_{i j}^{-1} T Y_{i j}=\bar{T}=\left(\bar{T}_{i j}\right)$, then $T$ and $\bar{T}$ are identical except that

$$
\begin{array}{ll}
\bar{T}_{i j}=T_{i i} Z_{i j}-Z_{i j} T_{j j}+T_{i j}, & \\
\bar{T}_{i k}=T_{i k}-Z_{i j} T_{j k}, & (k=j+1: q) \\
\bar{T}_{k j}=T_{k i} Z_{i j}+T_{k j}, & (k=1: i-1)
\end{array}
$$

Thus, $T_{i j}$ can be zeroed provided we have an algorithm for solving the Sylvester equation

$$
F Z-Z G=C
$$

where $F \in \mathbb{R}^{p \times p}$ and $G \in \mathbb{R}^{r \times r}$ are given upper quasi-triangular matrices and $C \in \mathbb{R}^{p \times r}$.
Bartels and Stewart (1972) have devised a method for doing this. Let $C= \left[c_{1}|\cdots| c_{r}\right]$ and $Z=\left[z_{1}|\cdots| z_{r}\right]$ be column partitionings. If $g_{k+1, k}=0$, then by comparing columns in (7.6.4) we find

$$
F z_{k}-\sum_{i=1}^{k} g_{i k} z_{i}=c_{k}
$$

Thus, once we know $z_{1}, \ldots, z_{k-1}$, then we can solve the quasi-triangular system

$$
\left(F-g_{k k} I\right) z_{k}=c_{k}+\sum_{i=1}^{k-1} g_{i k} z_{i}
$$

for $z_{k}$. If $g_{k+1, k} \neq 0$, then $z_{k}$ and $z_{k+1}$ can be simultaneously found by solving the $2 p$-by- $2 p$ system

$$
\left[\begin{array}{cc}
F-g_{k k} I & -g_{m k} I \\
-g_{k m} I & F-g_{m m} I
\end{array}\right]\left[\begin{array}{c}
z_{k} \\
z_{m}
\end{array}\right]=\left[\begin{array}{c}
c_{k} \\
c_{m}
\end{array}\right]+\sum_{i=1}^{k-1}\left[\begin{array}{c}
g_{i k} z_{i} \\
g_{i m} z_{i}
\end{array}\right]
$$

where $m=k+1$. By reordering the equations according to the perfect shuffle permutation $(1, p+1,2, p+2, \ldots, p, 2 p)$, a banded system is obtained that can be solved in $O\left(p^{2}\right)$ flops. The details may be found in Bartels and Stewart (1972). Here is the overall process for the case when $F$ and $G$ are each triangular.
Algorithm 7.6.2 (Bartels-Stewart Algorithm) Given $C \in \mathbb{R}^{p \times r}$ and upper triangular matrices $F \in \mathbb{R}^{p \times p}$ and $G \in \mathbb{R}^{r \times r}$ that satisfy $\lambda(F) \cap \lambda(G)=\emptyset$, the following algorithm overwrites $C$ with the solution to the equation $F Z-Z G=C$.
for $k=1: r$

$$
\begin{aligned}
& C(1: p, k)=C(1: p, k)+C(1: p, 1: k-1) \cdot G(1: k-1, k) \\
& \text { Solve }(F-G(k, k) I) z=C(1: p, k) \text { for } z \\
& C(1: p, k)=z
\end{aligned}
$$

This algorithm requires $p r(p+r)$ flops. By zeroing the superdiagonal blocks in $T$ in the appropriate order, the entire matrix can be reduced to block diagonal form.

Algorithm 7.6.3 Given an orthogonal matrix $Q \in \mathbb{R}^{n \times n}$, an upper quasi-triangular matrix $T=Q^{T} A Q$, and the partitioning (7.6.3), the following algorithm overwrites $Q$ with $Q Y$ where $Y^{-1} T Y=\operatorname{diag}\left(T_{11}, \ldots, T_{q q}\right)$.

```
for $j=2: q$
    for $i=1: j-1$
        for $k=j+1: q$
            $T_{i k}=T_{i k}-Z T_{j k}$
        end
        for $k=1: q$
            $Q_{k j}=Q_{k i} Z+Q_{k j}$
        end
    end
end
```

        Solve $T_{i i} Z-Z T_{j j}=-T_{i j}$ for $Z$ using the Bartels-Stewart algorithm.
    The number of flops required by this algorithm is a complicated function of the block sizes in (7.6.3).

The choice of the real Schur form $T$ and its partitioning in (7.6.3) determines the sensitivity of the Sylvester equations that must be solved in Algorithm 7.6.3. This in turn affects the condition of the matrix $Y$ and the overall usefulness of the block diagonalization. The reason for these dependencies is that the relative error of the computed solution $\hat{Z}$ to

$$
T_{i i} Z-Z T_{j j}=-T_{i j}
$$

satisfies

$$
\frac{\|\hat{Z}-Z\|_{F}}{\|Z\|_{F}} \approx \mathbf{u} \frac{\|T\|_{F}}{\operatorname{sep}\left(T_{i i}, T_{j j}\right)} .
$$

For details, see Golub, Nash, and Van Loan (1979). Since

$$
\operatorname{sep}\left(T_{i i}, T_{j j}\right)=\min _{X \neq 0} \frac{\left\|T_{i i} X-X T_{j j}\right\|_{F}}{\|X\|_{F}} \leq \min _{\substack{\lambda \in \lambda\left(T_{i i}\right) \\ \mu \in \lambda\left(T_{j j}\right)}}|\lambda-\mu|
$$

there can be a substantial loss of accuracy whenever the subsets $\lambda\left(T_{i i}\right)$ are insufficiently separated. Moreover, if $Z$ satisfies (7.6.6) then

$$
\|Z\|_{F} \leq \frac{\left\|T_{i j}\right\|_{F}}{\operatorname{sep}\left(T_{i i}, T_{j j}\right)}
$$

Thus, large norm solutions can be expected if $\operatorname{sep}\left(T_{i i}, T_{j j}\right)$ is small. This tends to make the matrix $Y$ in Algorithm 7.6.3 ill-conditioned since it is the product of the matrices

$$
Y_{i j}=\left[\begin{array}{rr}
I_{n_{i}} & Z \\
0 & I_{n_{j}}
\end{array}\right]
$$

Note that $\kappa_{F}\left(Y_{i j}\right)=n_{i}^{2}+n_{j}^{2}+\|Z\|_{F}^{2}$.

Confronted with these difficulties, Bavely and Stewart (1979) develop an algorithm for block diagonalizing that dynamically determines the eigenvalue ordering and partitioning in (7.6.3) so that all the $Z$ matrices in Algorithm 7.6.3 are bounded in norm by some user-supplied tolerance. Their research suggests that the condition of $Y$ can be controlled by controlling the condition of the $Y_{i j}$.

### 7.6.4 Eigenvector Bases

If the blocks in the partitioning (7.6.3) are all 1-by-1, then Algorithm 7.6.3 produces a basis of eigenvectors. As with the method of inverse iteration, the computed eigenvalueeigenvector pairs are exact for some "nearby" matrix. A widely followed rule of thumb for deciding upon a suitable eigenvector method is to use inverse iteration whenever fewer than $25 \%$ of the eigenvectors are desired.

We point out, however, that the real Schur form can be used to determine selected eigenvectors. Suppose

$$
Q^{T} A Q=\underset{k-1}{\left[\begin{array}{ccc}
T_{11} & u & T_{13} \\
0 & \lambda & v^{T} \\
0 & 0 & T_{33}
\end{array}\right]} \underset{n-k}{\begin{array}{c}
k-1 \\
1
\end{array}}
$$

is upper quasi-triangular and that $\lambda \notin \lambda\left(T_{11}\right) \cup \lambda\left(T_{33}\right)$. It follows that if we solve the linear systems $\left(T_{11}-\lambda I\right) w=-u$ and $\left(T_{33}-\lambda I\right)^{T} z=-v$ then

$$
x=Q\left[\begin{array}{c}
w \\
1 \\
0
\end{array}\right] \quad \text { and } \quad y=Q\left[\begin{array}{c}
0 \\
1 \\
z
\end{array}\right]
$$

are the associated right and left eigenvectors, respectively. Note that the condition of $\lambda$ is prescribed by

$$
1 / s(\lambda)=\sqrt{\left(1+w^{T} w\right)\left(1+z^{T} z\right)} .
$$

### 7.6.5 Ascertaining Jordan Block Structures

Suppose that we have computed the real Schur decomposition $A=Q T Q^{T}$, identified clusters of "equal" eigenvalues, and calculated the corresponding block diagonalization $T=Y \cdot \operatorname{diag}\left(T_{11}, \ldots, T_{q q}\right) Y^{-1}$. As we have seen, this can be a formidable task. However, even greater numerical problems confront us if we attempt to ascertain the Jordan block structure of each $T_{i i}$. A brief examination of these difficulties will serve to highlight the limitations of the Jordan decomposition.

Assume for clarity that $\lambda\left(T_{i i}\right)$ is real. The reduction of $T_{i i}$ to Jordan form begins by replacing it with a matrix of the form $C=\lambda I+N$, where $N$ is the strictly upper triangular portion of $T_{i i}$ and where $\lambda$, say, is the mean of its eigenvalues.

Recall that the dimension of a Jordan block $J(\lambda)$ is the smallest nonnegative integer $k$ for which $[J(\lambda)-\lambda I]^{k}=0$. Thus, if $p_{i}=\operatorname{dim}\left[\right.$ null $\left.\left(N^{i}\right)\right]$, for $i=0: n$, then $p_{i}-p_{i-1}$ equals the number of blocks in $C$ 's Jordan form that have dimension $i$ or
greater. A concrete example helps to make this assertion clear and to illustrate the role of the SVD in Jordan form computations.

Assume that $C$ is 7-by-7. Suppose we compute the SVD $U_{1}^{T} N V_{1}=\Sigma_{1}$ and "discover" that $N$ has rank 3. If we order the singular values from small to large then it follows that the matrix $N_{1}=V_{1}^{T} N V_{1}$ has the form

$$
N_{1}=\underset{4}{\left[\begin{array}{cc}
0 & K \\
0 & L
\end{array}\right]_{3}^{4}}
$$

At this point, we know that the geometric multiplicity of $\lambda$ is 4-i.e, $C$ 's Jordan form has four blocks ( $p_{1}-p_{0}=4-0=4$ ).

Now suppose $\tilde{U}_{2}^{T} L \tilde{V}_{2}=\Sigma_{2}$ is the SVD of $L$ and that we find that $L$ has unit rank. If we again order the singular values from small to large, then $L_{2}=\tilde{V}_{2}^{T} L \tilde{V}_{2}$ clearly has the following structure:

$$
L_{2}=\left[\begin{array}{lll}
0 & 0 & a \\
0 & 0 & b \\
0 & 0 & c
\end{array}\right] .
$$

However, $\lambda\left(L_{2}\right)=\lambda(L)=\{0,0,0\}$ and so $c=0$. Thus, if

$$
V_{2}=\operatorname{diag}\left(I_{4}, \tilde{V}_{2}\right)
$$

then $N_{2}=V_{2}^{T} N_{1} V_{2}$ has the following form:

$$
N_{2}=\left[\begin{array}{lllllll}
0 & 0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & 0 & 0 & a \\
0 & 0 & 0 & 0 & 0 & 0 & b \\
0 & 0 & 0 & 0 & 0 & 0 & 0
\end{array}\right] .
$$

Besides allowing us to introduce more zeros into the upper triangle, the SVD of $L$ also enables us to deduce the dimension of the nullspace of $N^{2}$. Since

$$
\text { . } \quad N_{1}^{2}=\left[\begin{array}{cc}
0 & K L \\
0 & L^{2}
\end{array}\right]=\left[\begin{array}{cc}
0 & K \\
0 & L
\end{array}\right]\left[\begin{array}{cc}
0 & K \\
0 & L
\end{array}\right]
$$

and $\left[\begin{array}{c}K \\ L\end{array}\right]$ has full column rank,

$$
p_{2}=\operatorname{dim}\left(\operatorname{null}\left(N^{2}\right)\right)=\operatorname{dim}\left(\operatorname{null}\left(N_{1}^{2}\right)\right)=4+\operatorname{dim}(\operatorname{null}(L))=p_{1}+2 .
$$

Hence, we can conclude at this stage that the Jordan form of $C$ has at least two blocks of dimension 2 or greater.

Finally, it is easy to see that $N_{1}^{3}=0$, from which we conclude that there is $p_{3}-p_{2} =7-6=1$ block of dimension 3 or larger. If we define $V=V_{1} V_{2}$ then it follows that
the decomposition

$$
V^{T} C V=\left[\begin{array}{ccccccc}
\lambda & 0 & 0 & 0 & \times & \times & \times \\
0 & \lambda & 0 & 0 & \times & \times & \times \\
0 & 0 & \lambda & 0 & \times & \times & \times \\
0 & 0 & 0 & \lambda & \times & \times & \times \\
0 & 0 & 0 & 0 & \lambda & \times & a \\
0 & 0 & 0 & 0 & 0 & \lambda & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & \lambda
\end{array}\right] \begin{aligned}
& \text { two blocks of order } 2 \text { or larger } \\
& \text { four blocks of order } 1 \text { or larger } \\
& \text { one block of order } 3 \text { or larger }
\end{aligned}
$$

displays $C$ 's Jordan block structure: two blocks of order 1, one block of order 2, and one block of order 3.

To compute the Jordan decomposition it is necessary to resort to nonorthogonal transformations. We refer the reader to Golub and Wilkinson (1976), Kågström and Ruhe (1980a, 1980b), and Demmel (1983) for more details. The above calculations with the SVD amply illustrate that difficult rank decisions must be made at each stage and that the final computed block structure depends critically on those decisions.

## Problems

P7.6.1 Give a complete algorithm for solving a real, $n$-by- $n$, upper quasi-triangular system $T x=b$.
P7.6.2 Suppose $U^{-1} A U=\operatorname{diag}\left(\alpha_{1}, \ldots, \alpha_{m}\right)$ and $V^{-1} B V=\operatorname{diag}\left(\beta_{1}, \ldots, \beta_{n}\right)$. Show that if

$$
\phi(X)=A X-X B
$$

then

$$
\lambda(\phi)=\left\{\alpha_{i}-\beta_{j}: i=1: m, j=1: n\right\} .
$$

What are the corresponding eigenvectors? How can these facts be used to solve $A X-X B=C$ ?
P7.6.3 Show that if $Z \in \mathbb{C}^{p \times q}$ and

$$
Y=\left[\begin{array}{rc}
I_{p} & Z \\
0 & I_{q}
\end{array}\right]
$$

then $\kappa_{2}(Y)=\left[2+\sigma^{2}+\sqrt{4 \sigma^{2}+\sigma^{4}}\right] / 2$ where $\sigma=\|Z\|_{2}$.
P7.6.4 Derive the system (7.6.5).
P7.6.5 Assume that $T \in \mathbf{R}^{n \times n}$ is block upper triangular and partitioned as follows:

$$
T=\left[\begin{array}{rrr}
T_{11} & T_{12} & T_{13} \\
0 & T_{22} & T_{23} \\
0 & 0 & T_{33}
\end{array}\right], \quad T \in \mathbf{R}^{n \times n}
$$

Suppose that the diagonal block $T_{22}$ is 2-by-2 with complex eigenvalues that are disjoint from $\lambda\left(T_{11}\right)$ and $\lambda\left(T_{33}\right)$. Give an algorithm for computing the 2 -dimensional real invariant subspace associated with $T_{22}$ 's eigenvalues.

P7.6.6 Suppose $H \in \mathbf{R}^{n \times n}$ is upper Hessenberg with a complex eigenvalue $\lambda+i \cdot \mu$. How could inverse iteration be used to compute $x, y \in \mathbf{R}^{n}$ so that $H(x+i y)=(\lambda+i \mu)(x+i y)$ ? Hint: Compare real and imaginary parts in this equation and obtain a $2 n$-by- $2 n$ real system.

## Notes and References for §7.6

Much of the material discussed in this section may be found in the following survey paper:
G.H. Golub and J.H. Wilkinson (1976). "Ill-Conditioned Eigensystems and the Computation of the Jordan Canonical Form," SIAM Review 18, 578-619.

The problem of ordering the eigenvalues in the real Schur form is the subject of:
A. Ruhe (1970). "An Algorithm for Numerical Determination of the Structure of a General Matrix," BIT 10, 196-216.
G.W. Stewart (1976). "Algorithm 406: HQR3 and EXCHNG: Fortran Subroutines for Calculating and Ordering the Eigenvalues of a Real Upper Hessenberg Matrix," ACM Trans. Math. Softw. 2, 275-280.
J.J. Dongarra, S. Hammarling, and J.H. Wilkinson (1992). "Numerical Considerations in Computing Invariant Subspaces," SIAM J. Matrix Anal. Applic. 13, 145-161.
Z. Bai and J.W. Demmel (1993). "On Swapping Diagonal Blocks in Real Schur Form," Lin. Alg. Applic. 186, 73-95

Procedures for block diagonalization including the Jordan form are described in:
C. Bavely and G.W. Stewart (1979). "An Algorithm for Computing Reducing Subspaces by Block Diagonalization," SIAM J. Numer. Anal. 16, 359-367.
B. Kågström and A. Ruhe (1980a). "An Algorithm for Numerical Computation of the Jordan Normal Form of a Complex Matrix," ACM Trans. Math. Softw. 6, 398-419.
B. Kågström and A. Ruhe (1980b). "Algorithm 560 JNF: An Algorithm for Numerical Computation of the Jordan Normal Form of a Complex Matrix," ACM Trans. Math. Softw. 6, 437-443.
J.W. Demmel (1983). "A Numerical Analyst's Jordan Canonical Form," PhD Thesis, Berkeley.
N. Ghosh, W.W. Hager, and P. Sarmah (1997). "The Application of Eigenpair Stability to Block Diagonalization," SIAM J. Numer. Anal. 34, 1255-1268.
S. Serra-Capizzano, D. Bertaccini, and G.H. Golub (2005). "How to Deduce a Proper Eigenvalue Cluster from a Proper Singular Value Cluster in the Nonnormal Case," SIAM J. Matrix Anal. Applic. 27, 82-86.

Before we offer pointers to the literature associated with invariant subspace computation, we remind the reader that in §7.3 we discussed the power method for computing the dominant eigenpair and the method of orthogonal iteration that can be used to compute dominant invariant subspaces. Inverse iteration is a related idea and is the concern of the following papers:
J. Varah (1968). "The Calculation of the Eigenvectors of a General Complex Matrix by Inverse Iteration," Math. Comput. 22, 785-791.
J. Varah (1970). "Computing Invariant Subspaces of a General Matrix When the Eigensystem is Poorly Determined," Math. Comput. 24, 137-149.
G. Peters and J.H. Wilkinson (1979). "Inverse Iteration, Ill-Conditioned Equations, and Newton's Method," SIAM Review 21, 339-360.
I.C.F. Ipsen (1997). "Computing an Eigenvector with Inverse Iteration," SIAM Review 39, 254-291.

In certain applications it is necessary to track an invariant subspace as the matrix changes, see:
L. Dieci and M.J. Friedman (2001). "Continuation of Invariant Subspaces," Num. Lin. Alg. 8, 317-327.
D. Bindel, J.W. Demmel, and M. Friedman (2008). "Continuation of Invariant Subsapces in Large Bifurcation Problems," SIAM J. Sci. Comput. 30, 637-656.

Papers concerned with estimating the error in a computed eigenvalue and/or eigenvector include:
S.P. Chan and B.N. Parlett (1977). "Algorithm 517: A Program for Computing the Condition Numbers of Matrix Eigenvalues Without Computing Eigenvectors," ACM Trans. Math. Softw. 3, 186-203.
H.J. Symm and J.H. Wilkinson (1980). "Realistic Error Bounds for a Simple Eigenvalue and Its Associated Eigenvector," Numer. Math. 35, 113-126.
C. Van Loan (1987). "On Estimating the Condition of Eigenvalues and Eigenvectors," Lin. Alg. Applic. 88/89, 715-732.
Z. Bai, J. Demmel, and A. McKenney (1993). "On Computing Condition Numbers for the Nonsymmetric Eigenproblem," ACM Trans. Math. Softw. 19, 202-223.

Some ideas about improving computed eigenvalues, eigenvectors, and invariant subspaces may be found in:
J. Varah (1968). "Rigorous Machine Bounds for the Eigensystem of a General Complex Matrix," Math. Comp. 22, 793-801.
J.J. Dongarra, C.B. Moler, and J.H. Wilkinson (1983). "Improving the Accuracy of Computed Eigenvalues and Eigenvectors," SIAM J. Numer. Anal. 20, 23-46.
J.W. Demmel (1987). "Three Methods for Refining Estimates of Invariant Subspaces," Comput. 38, 43-57.

As we have seen, the sep(.,.) function is of great importance in the assessment of a computed invariant subspace. Aspects of this quantity and the associated Sylvester equation are discussed in:
J. Varah (1979). "On the Separation of Two Matrices," SIAM J. Numer. Anal. 16, 212-222.
R. Byers (1984). "A Linpack-Style Condition Estimator for the Equation $A X-X B^{T}=C$," $I E E E$ Trans. Autom. Contr. AC-29, 926-928.
M. Gu and M.L. Overton (2006). "An Algorithm to Compute Sep ${ }_{\lambda}$," SIAM J. Matrix Anal. Applic. 28, 348-359.
N.J. Higham (1993). "Perturbation Theory and Backward Error for AX - XB = C," BIT 33, 124-136.

Sylvester equations arise in many settings, and there are many solution frameworks, see:
R.H. Bartels and G.W. Stewart (1972). "Solution of the Equation $A X+X B=C$," Commun. ACM 15, 820-826.
G.H. Golub, S. Nash, and C. Van Loan (1979). "A Hessenberg-Schur Method for the Matrix Problem $A X+X B=C$," IEEE Trans. Autom. Contr. AC-24, 909-913.
K. Datta (1988). "The Matrix Equation $X A-B X=R$ and Its Applications," Lin. Alg. Applic. 109, 91-105.
B. Kågström and P. Poromaa (1992). "Distributed and Shared Memory Block Algorithms for the Triangular Sylvester Equation with sep ${ }^{-1}$ Estimators," SIAM J. Matrx Anal. Applic. 13, 90101.
J. Gardiner, M.R. Wette, A.J. Laub, J.J. Amato, and C.B. Moler (1992). "Algorithm 705: A FORTRAN-77 Software Package for Solving the Sylvester Matrix Equation $A X B^{T}+C X D^{T}=E$," ACM Trans. Math. Softw. 18, 232-238.
V. Simoncini (1996). "On the Numerical Solution of AX -XB =C," BIT 36, 814-830.
C.H. Bischof, B.N Datta, and A. Purkayastha (1996). "A Parallel Algorithm for the Sylvester Observer Equation," SIAM J. Sci. Comput. 17, 686-698.
D. Calvetti, B. Lewis, L. Reichel (2001). "On the Solution of Large Sylvester-Observer Equations," Num. Lin. Alg. 8, 435-451.

The constrained Sylvester equation problem is considered in:
J.B. Barlow, M.M. Monahemi, and D.P. O'Leary (1992). "Constrained Matrix Sylvester Equations," SIAM J. Matrix Anal. Applic. 13, 1-9.
A.R. Ghavimi and A.J. Laub (1996). "Numerical Methods for Nearly Singular Constrained Matrix Sylvester Equations." SIAM J. Matrix Anal. Applic. 17, 212-221.

The Lyapunov problem $F X+X F^{T}=-C$ where $C$ is non-negative definite has a very important role to play in control theory, see:
G. Hewer and C. Kenney (1988). "The Sensitivity of the Stable Lyapunov Equation," SIAM J. Control Optim 26, 321-344.
A.R. Ghavimi and A.J. Laub (1995). "Residual Bounds for Discrete-Time Lyapunov Equations," IEEE Trans. Autom. Contr. 40, 1244-1249.
J.-R. Li and J. White (2004). "Low-Rank Solution of Lyapunov Equations," SIAM Review 46, 693713.

Several authors have considered generalizations of the Sylvester equation, i.e., $\Sigma F_{i} X G_{i}=C$. These include:
P. Lancaster (1970). "Explicit Solution of Linear Matrix Equations," SIAM Review 12, 544-566.
H. Wimmer and A.D. Ziebur (1972). "Solving the Matrix Equations $\Sigma f_{p}(A) g_{p}(A)=C$," SIAM Review 14, 318-323.
W.J. Vetter (1975). "Vector Structures and Solutions of Linear Matrix Equations," Lin. Alg. Applic. 10, 181-188.

### 7.7 The Generalized Eigenvalue Problem

If $A, B \in \mathbb{C}^{n \times n}$, then the set of all matrices of the form $A-\lambda B$ with $\lambda \in \mathbb{C}$ is a pencil. The generalized eigenvalues of $A-\lambda B$ are elements of the set $\lambda(A, B)$ defined by

$$
\lambda(A, B)=\{z \in \mathbb{C}: \operatorname{det}(A-z B)=0\}
$$

If $\lambda \in \lambda(A, B)$ and $0 \neq x \in \mathbb{C}^{n}$ satisfies

$$
A x=\lambda B x
$$

then $x$ is an eigenvector of $A-\lambda B$. The problem of finding nontrivial solutions to (7.7.1) is the generalized eigenvalue problem and in this section we survey some of its mathematical properties and derive a stable method for its solution. We briefly discuss how a polynomial eigenvalue problem can be converted into an equivalent generalized eigenvalue problem through a linearization process.

### 7.7.1 Background

The first thing to observe about the generalized eigenvalue problem is that there are $n$ eigenvalues if and only if $\operatorname{rank}(B)=n$. If $B$ is rank deficient then $\lambda(A, B)$ may be finite, empty, or infinite:

$$
\begin{aligned}
& A=\left[\begin{array}{ll}
1 & 2 \\
0 & 3
\end{array}\right], \quad B=\left[\begin{array}{ll}
1 & 0 \\
0 & 0
\end{array}\right] \Rightarrow \lambda(A, B)=\{1\}, \\
& A=\left[\begin{array}{ll}
1 & 2 \\
0 & 3
\end{array}\right], \quad B=\left[\begin{array}{ll}
0 & 1 \\
0 & 0
\end{array}\right] \Rightarrow \lambda(A, B)=\emptyset, \\
& A=\left[\begin{array}{ll}
1 & 2 \\
0 & 0
\end{array}\right], \quad B=\left[\begin{array}{ll}
1 & 0 \\
0 & 0
\end{array}\right] \Rightarrow \lambda(A, B)=\mathbb{C} .
\end{aligned}
$$

Note that if $0 \neq \lambda \in \lambda(A, B)$, then $(1 / \lambda) \in \lambda(B, A)$. Moreover, if $B$ is nonsingular, then $\lambda(A, B)=\lambda\left(B^{-1} A, I\right)=\lambda\left(B^{-1} A\right)$. This last observation suggests one method for solving the $A-\lambda B$ problem if $B$ is nonsingular:

Step 1. Solve $B C=A$ for $C$ using (say) Gaussian elimination with pivoting.
Step 2. Use the QR algorithm to compute the eigenvalues of $C$.
In this framework, $C$ is affected by roundoff errors of order $\mathbf{u}\|A\|_{2}\left\|B^{-1}\right\|_{2}$. If $B$ is illconditioned, then this precludes the possibility of computing any generalized eigenvalue accurately-even those eigenvalues that may be regarded as well-conditioned. For example, if

$$
A=\left[\begin{array}{rr}
1.746 & .940 \\
1.246 & 1.898
\end{array}\right] \quad \text { and } \quad B=\left[\begin{array}{ll}
.780 & .563 \\
.913 & .659
\end{array}\right]
$$

then $\lambda(A, B)=\left\{2,1.07 \times 10^{6}\right\}$. With 7-digit floating point arithmetic, we find $\lambda\left(\mathrm{f}\left(A B^{-1}\right)\right)=\left\{1.562539,1.01 \times 10^{6}\right\}$. The poor quality of the small eigenvalue is because $\kappa_{2}(B) \approx 2 \times 10^{6}$. On the other hand, we find that

$$
\lambda\left(I, \mathrm{fl}\left(A^{-1} B\right)\right) \approx\left\{2.000001,1.06 \times 10^{6}\right\} .
$$

The accuracy of the small eigenvalue is improved because $\kappa_{2}(A) \approx 4$.
The example suggests that we seek an alternative approach to the generalized eigenvalue problem. One idea is to compute well-conditioned $Q$ and $Z$ such that the matrices

$$
A_{1}=Q^{-1} A Z, \quad B_{1}=Q^{-1} B Z
$$

are each in canonical form. Note that $\lambda(A, B)=\lambda\left(A_{1}, B_{1}\right)$ since

$$
A x=\lambda B x \quad \Leftrightarrow \quad A_{1} y=\lambda B_{1} y, \quad x=Z y .
$$

We say that the pencils $A-\lambda B$ and $A_{1}-\lambda B_{1}$ are equivalent if (7.7.2) holds with nonsingular $Q$ and $Z$.

As in the standard eigenproblem $A-\lambda I$ there is a choice between canonical forms. Corresponding to the Jordan form is a decomposition of Kronecker in which both $A_{1}$ and $B_{1}$ are block diagonal with blocks that are similar in structure to Jordan blocks. The Kronecker canonical form poses the same numerical challenges as the Jordan form, but it provides insight into the mathematical properties of the pencil $A-\lambda B$. See Wilkinson (1978) and Demmel and Kågström (1987) for details.

### 7.7.2 The Generalized Schur Decomposition

From the numerical point of view, it makes to insist that the transformation matrices $Q$ and $Z$ be unitary. This leads to the following decomposition described in Moler and Stewart (1973).

Theorem 7.7.1 (Generalized Schur Decomposition). If $A$ and $B$ are in $\mathbb{C}^{n \times n}$, then there exist unitary $Q$ and $Z$ such that $Q^{H} A Z=T$ and $Q^{H} B Z=S$ are upper triangular. If for some $k, t_{k k}$ and $s_{k k}$ are both zero, then $\lambda(A, B)=\mathbb{C}$. Otherwise

$$
\lambda(A, B)=\left\{t_{i i} / s_{i i}: s_{i i} \neq 0\right\} .
$$

Proof. Let $\left\{B_{k}\right\}$ be a sequence of nonsingular matrices that converge to $B$. For each $k$, let

$$
Q_{k}^{H}\left(A B_{k}^{-1}\right) Q_{k}=R_{k}
$$

be a Schur decomposition of $A B_{k}^{-1}$. Let $Z_{k}$ be unitary such that

$$
Z_{k}^{H}\left(B_{k}^{-1} Q_{k}\right)=S_{k}^{-1}
$$

is upper triangular. It follows that $Q_{k}^{H} A Z_{k}=R_{k} S_{k}$ and $Q_{k}^{H} B_{k} Z_{k}=S_{k}$ are also upper triangular. Using the Bolzano-Weierstrass theorem, we know that the bounded sequence $\left\{\left(Q_{k}, Z_{k}\right)\right\}$ has a converging subsequence,

$$
\lim _{i \rightarrow \infty}\left(Q_{k_{i}}, Z_{k_{i}}\right)=(Q, Z) .
$$

It is easy to show that $Q$ and $Z$ are unitary and that $Q^{H} A Z$ and $Q^{H} B Z$ are upper triangular. The assertions about $\lambda(A, B)$ follow from the identity

$$
\operatorname{det}(A-\lambda B)=\operatorname{det}\left(Q Z^{H}\right) \prod_{i=1}^{n}\left(t_{i i}-\lambda s_{i i}\right)
$$

and that completes the proof of the theorem.
If $A$ and $B$ are real then the following decomposition, which corresponds to the real Schur decomposition (Theorem 7.4.1), is of interest.

Theorem 7.7.2 (Generalized Real Schur Decomposition). If $A$ and $B$ are in $\mathbb{R}^{n \times n}$ then there exist orthogonal matrices $Q$ and $Z$ such that $Q^{T} A Z$ is upper quasitriangular and $Q^{T} B Z$ is upper triangular.

Proof. See Stewart (1972).
In the remainder of this section we are concerned with the computation of this decomposition and the mathematical insight that it provides.

### 7.7.3 Sensitivity Issues

The generalized Schur decomposition sheds light on the issue of eigenvalue sensitivity for the $A-\lambda B$ problem. Clearly, small changes in $A$ and $B$ can induce large changes in the eigenvalue $\lambda_{i}=t_{i i} / s_{i i}$ if $s_{i i}$ is small. However, as Stewart (1978) argues, it may not be appropriate to regard such an eigenvalue as "ill-conditioned." The reason is that the reciprocal $\mu_{i}=s_{i i} / t_{i i}$ might be a very well-behaved eigenvalue for the pencil $\mu A-B$. In the Stewart analysis, $A$ and $B$ are treated symmetrically and the eigenvalues are regarded more as ordered pairs $\left(t_{i i}, s_{i i}\right)$ than as quotients. With this point of view it becomes appropriate to measure eigenvalue perturbations in the chordal metric chord( $a, b$ ) defined by

$$
\operatorname{chord}(a, b)=\frac{|a-b|}{\sqrt{1+a^{2}} \sqrt{1+b^{2}}}
$$

Stewart shows that if $\lambda$ is a distinct eigenvalue of $A-\lambda B$ and $\lambda_{\epsilon}$ is the corresponding eigenvalue of the perturbed pencil $\tilde{A}-\lambda \tilde{B}$ with $\|A-\tilde{A}\|_{2} \approx\|B-\tilde{B}\|_{2} \approx \epsilon$, then

$$
\operatorname{chord}\left(\lambda, \lambda_{\epsilon}\right) \leq \frac{\epsilon}{\sqrt{\left(y^{H} A x\right)^{2}+\left(y^{H} B x\right)^{2}}}+O\left(\epsilon^{2}\right)
$$

where $x$ and $y$ have unit 2 -norm and satisfy $A x=\lambda B x$ and $y^{H} A=\lambda y^{H} B$. Note that the denominator in the upper bound is symmetric in $A$ and $B$. The "truly" ill-conditioned eigenvalues are those for which this denominator is small.

The extreme case when both $t_{k k}$ and $s_{k k}$ are zero for some $k$ has been studied by Wilkinson (1979). In this case, the remaining quotients $t_{i i} / s_{i i}$ can take on arbitrary values.

### 7.7.4 Hessenberg-Triangular Form

The first step in computing the generalized real Schur decomposition of the pair ( $A, B$ ) is to reduce $A$ to upper Hessenberg form and $B$ to upper triangular form via orthogonal transformations. We first determine an orthogonal $U$ such that $U^{T} B$ is upper triangular. Of course, to preserve eigenvalues, we must also update $A$ in exactly the same way. Let us trace what happens in the $n=5$ case.

$$
A \leftarrow U^{T} A=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times
\end{array}\right], \quad B \leftarrow U^{T} B=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & \times
\end{array}\right] .
$$

Next, we reduce $A$ to upper Hessenberg form while preserving $B$ 's upper triangular form. First, a Givens rotation $Q_{45}$ is determined to zero $a_{51}$ :

$$
A \leftarrow Q_{45}^{T} A=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times
\end{array}\right], \quad B \leftarrow Q_{45}^{T} B=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & \times & \times
\end{array}\right] .
$$

The nonzero entry arising in the $(5,4)$ position in $B$ can be zeroed by postmultiplying with an appropriate Givens rotation $Z_{45}$ :

$$
A \leftarrow A Z_{45}=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times
\end{array}\right], \quad B \leftarrow B Z_{45}=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & \times
\end{array}\right] .
$$

Zeros are similarly introduced into the $(4,1)$ and $(3,1)$ positions in $A$ :

$$
\begin{array}{ll}
A \leftarrow Q_{34}^{T} A=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & \times & \times & \times & \times
\end{array}\right], & B \leftarrow Q_{34}^{T} B=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times
\end{array}\right] \\
A \leftarrow A Z_{34}=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & \times & \times & \times & \times
\end{array}\right], & B \leftarrow B Z_{34}=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & \times
\end{array}\right] \\
A \leftarrow Q_{23}^{T} A=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & \times & \times & \times & \times
\end{array}\right], & B \leftarrow Q_{23}^{T} B=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & \times
\end{array}\right]
\end{array}
$$

$$
A \leftarrow A Z_{23}=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & \times & \times & \times & \times
\end{array}\right], \quad B \leftarrow B Z_{23}=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & \times
\end{array}\right] .
$$

$A$ is now upper Hessenberg through its first column. The reduction is completed by zeroing $a_{52}, a_{42}$, and $a_{53}$. Note that two orthogonal transformations are required for each $a_{i j}$ that is zeroed-one to do the zeroing and the other to restore $B$ 's triangularity. Either Givens rotations or 2-by-2 modified Householder transformations can be used. Overall we have:

Algorithm 7.7.1 (Hessenberg-Triangular Reduction) Given $A$ and $B$ in $\mathbb{R}^{n \times n}$, the following algorithm overwrites $A$ with an upper Hessenberg matrix $Q^{T} A Z$ and $B$ with an upper triangular matrix $Q^{T} B Z$ where both $Q$ and $Z$ are orthogonal.

Compute the factorization $B=Q R$ using Algorithm 5.2.1 and overwrite $A$ with $Q^{T} A$ and $B$ with $Q^{T} B$.
for $j=1: n-2$

$$
\begin{aligned}
& \text { for } i=n:-1: j+2 \\
& \qquad \begin{array}{l}
{[c, s]=\operatorname{givens}(A(i-1, j), A(i, j))} \\
A(i-1: i, j: n)=\left[\begin{array}{rr}
c & s \\
-s & c
\end{array}\right]^{T} A(i-1: i, j: n) \\
B(i-1: i, i-1: n)=\left[\begin{array}{rr}
c & s \\
-s & c
\end{array}\right]^{T} B(i-1: i, i-1: n) \\
{[c, s]=\operatorname{givens}(-B(i, i), B(i, i-1))} \\
B(1: i, i-1: i)=B(1: i, i-1: i)\left[\begin{array}{rr}
c & s \\
-s & c
\end{array}\right] \\
A(1: n, i-1: i)=A(1: n, i-1: i)\left[\begin{array}{rr}
c & s \\
-s & c
\end{array}\right] \\
\text { end }
\end{array} \text { l }
\end{aligned}
$$

This algorithm requires about $8 n^{3}$ flops. The accumulation of $Q$ and $Z$ requires about $4 n^{3}$ and $3 n^{3}$ flops, respectively.

The reduction of $A-\lambda B$ to Hessenberg-triangular form serves as a "front end" decomposition for a generalized QR iteration known as the QZ iteration which we describe next.

### 7.7.5 Deflation

In describing the QZ iteration we may assume without loss of generality that $A$ is an unreduced upper Hessenberg matrix and that $B$ is a nonsingular upper triangular
matrix. The first of these assertions is obvious, for if $a_{k+1, k}=0$ then

$$
A-\lambda B=\left[\begin{array}{cc}
A_{11}-\lambda B_{11} & A_{12}-\lambda B_{12} \\
0 & A_{22}-\lambda B_{22} \\
k & n-k
\end{array}\right]_{n-k}^{k},
$$

and we may proceed to solve the two smaller problems $A_{11}-\lambda B_{11}$ and $A_{22}-\lambda B_{22}$. On the other hand, if $b_{k k}=0$ for some $k$, then it is possible to introduce a zero in $A$ 's $(n, n-1)$ position and thereby deflate. Illustrating by example, suppose $n=5$ and $k=3$ :

$$
A=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times
\end{array}\right], \quad B=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & \times
\end{array}\right] .
$$

The zero on $B$ 's diagonal can be "pushed down" to the $(5,5)$ position as follows using Givens rotations:

$$
\begin{array}{ll}
A \leftarrow Q_{34}^{T} A=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times
\end{array}\right], & B \leftarrow Q_{34}^{T} B=\left[\begin{array}{ccccc}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & \times \\
0 & 0 & 0 & 0 & \times
\end{array}\right], \\
A \leftarrow A Z_{23}=\left[\begin{array}{lllll}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times
\end{array}\right], & B \leftarrow B Z_{23}=\left[\begin{array}{lllll}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & \times \\
0 & 0 & 0 & 0 & \times
\end{array}\right], \\
A \leftarrow Q_{45}^{T} A=\left[\begin{array}{lllll}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & \times & \times & \times
\end{array}\right], & B \leftarrow Q_{45}^{T} B=\left[\begin{array}{lllll}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & \times \\
0 & 0 & 0 & 0 & 0
\end{array}\right], \\
A \leftarrow A Z_{34}=\left[\begin{array}{lllll}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times
\end{array}\right], & B \leftarrow B Z_{34}^{T}=\left[\begin{array}{lllll}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times \\
0 & 0 & 0 & 0 & 0
\end{array}\right], \\
A \leftarrow A Z_{45}=\left[\begin{array}{lllll}
\times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times
\end{array}\right], & B \leftarrow B Z_{45}=\left[\begin{array}{lllll}
\times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times \\
0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & 0
\end{array}\right],
\end{array}
$$

This zero-chasing technique is perfectly general and can be used to zero $a_{n, n-1}$ regardless of where the zero appears along $B$ 's diagonal.

### 7.7.6 The QZ Step

We are now in a position to describe a QZ step. The basic idea is to update $A$ and $B$ as follows

$$
(\bar{A}-\lambda \bar{B})=\bar{Q}^{T}(A-\lambda B) \bar{Z},
$$

where $\bar{A}$ is upper Hessenberg, $\bar{B}$ is upper triangular, $\bar{Q}$ and $\bar{Z}$ are each orthogonal, and $\bar{A} \bar{B}^{-1}$ is essentially the same matrix that would result if a Francis QR step (Algorithm 7.5.1) were explicitly applied to $A B^{-1}$. This can be done with some clever zero-chasing and an appeal to the implicit $Q$ theorem.

Let $M=A B^{-1}$ (upper Hessenberg) and let $v$ be the first column of the matrix $(M-a I)(M-b I)$, where $a$ and $b$ are the eigenvalues of $M$ 's lower 2-by-2 submatrix. Note that $v$ can be calculated in $O(1)$ flops. If $P_{0}$ is a Householder matrix such that $P_{0} v$ is a multiple of $e_{1}$, then

$$
A \leftarrow P_{0} A=\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & 0 & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times
\end{array}\right], \quad B \leftarrow P_{0} B=\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & 0 & \times
\end{array}\right] .
$$

The idea now is to restore these matrices to Hessenberg-triangular form by chasing the unwanted nonzero elements down the diagonal.

To this end, we first determine a pair of Householder matrices $Z_{1}$ and $Z_{2}$ to zero $b_{31}, b_{32}$, and $b_{21}$ :
$A \leftarrow A Z_{1} Z_{2}=\left[\begin{array}{cccccc}\times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ \times & \times & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times\end{array}\right], \quad B \leftarrow B Z_{1} Z_{2}=\left[\begin{array}{cccccc}\times & \times & \times & \times & \times & \times \\ 0 & \times & \times & \times & \times & \times \\ 0 & 0 & \times & \times & \times & \times \\ 0 & 0 & 0 & \times & \times & \times \\ 0 & 0 & 0 & 0 & \times & \times \\ 0 & 0 & 0 & 0 & 0 & \times\end{array}\right]$.
Then a Householder matrix $P_{1}$ is used to zero $a_{31}$ and $a_{41}$ :

$$
A \leftarrow P_{1} A=\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
\times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times
\end{array}\right], \quad B \leftarrow P_{1} B=\left[\begin{array}{cccccc}
\times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times \\
0 & 0 & 0 & 0 & \times & \times \\
0 & 0 & 0 & 0 & 0 & \times
\end{array}\right] .
$$

Notice that with this step the unwanted nonzero elements have been shifted down and to the right from their original position. This illustrates a typical step in the QZ iteration. Notice that $Q=Q_{0} Q_{1} \cdots Q_{n-2}$ has the same first column as $Q_{0}$. By the way the initial Householder matrix was determined, we can apply the implicit $Q$ theorem and assert that $A B^{-1}=Q^{T}\left(A B^{-1}\right) Q$ is indeed essentially the same matrix that we would obtain by applying the Francis iteration to $M=A B^{-1}$ directly. Overall we have the following algorithm.

Algorithm 7.7.2 (The QZ Step) Given an unreduced upper Hessenberg matrix $A \in \mathbb{R}^{n \times n}$ and a nonsingular upper triangular matrix $B \in \mathbb{R}^{n \times n}$, the following algorithm overwrites $A$ with the upper Hessenberg matrix $Q^{T} A Z$ and $B$ with the upper triangular matrix $Q^{T} B Z$ where $Q$ and $Z$ are orthogonal and $Q$ has the same first column as the orthogonal similarity transformation in Algorithm 7.5.1 when it is applied to $A B^{-1}$.

Let $M=A B^{-1}$ and compute $(M-a I)(M-b I) e_{1}=[x, y, z, 0, \ldots, 0]^{T}$ where $a$ and $b$ are the eigenvalues of $M$ 's lower 2-by-2.
for $k=1: n-2$
Find Householder $Q_{k}$ so $Q_{k}\left[\begin{array}{l}x \\ y \\ z\end{array}\right]=\left[\begin{array}{l}* \\ 0 \\ 0\end{array}\right]$.
$A=\operatorname{diag}\left(I_{k-1}, Q_{k}, I_{n-k-2}\right) \cdot A$
$B=\operatorname{diag}\left(I_{k-1}, Q_{k}, I_{n-k-2}\right) \cdot B$
Find Householder $Z_{k 1}$ so $\left[b_{k+2, k}\left|b_{k+2, k+1}\right| b_{k+2, k+2}\right] Z_{k 1}=[0|0| *]$.
$A=A \cdot \operatorname{diag}\left(I_{k-1}, Z_{k 1}, I_{n-k-2}\right)$
$B=B \cdot \operatorname{diag}\left(I_{k-1}, Z_{k 1}, I_{n-k-2}\right)$
Find Householder $Z_{k 2}$ so $\left[b_{k+1, k} \mid b_{k+1, k+1}\right] Z_{k 2}=[0 \mid *]$.
$A=A \cdot \operatorname{diag}\left(I_{k-1}, Z_{k 2}, I_{n-k-1}\right)$
$B=B \cdot \operatorname{diag}\left(I_{k-1}, Z_{k 2}, I_{n-k-1}\right)$
$x=a_{k+1, k} ; y=a_{k+2, k}$
if $k<n-2$

$$
z=a_{k+3, k}
$$

end
end
Find Householder $Q_{n-1}$ so $Q_{n-1}\left[\begin{array}{l}x \\ y\end{array}\right]=\left[\begin{array}{c}* \\ 0\end{array}\right]$.
$A=\operatorname{diag}\left(I_{n-2}, Q_{n-1}\right) \cdot A$
$B=\operatorname{diag}\left(I_{n-2}, Q_{n-1}\right) \cdot B$.
Find Householder $Z_{n-1}$ so $\left[b_{n, n-1} \mid b_{n n}\right] Z_{n-1}=[0 \mid *]$.
$A=A \cdot \operatorname{diag}\left(I_{n-2}, Z_{n-1}\right)$
$B=B \cdot \operatorname{diag}\left(I_{n-2}, Z_{n-1}\right)$
This algorithm requires $22 n^{2}$ flops. $Q$ and $Z$ can be accumulated for an additional $8 n^{2}$ flops and $13 n^{2}$ flops, respectively.

### 7.7.7 The Overall QZ Process

By applying a sequence of QZ steps to the Hessenberg-triangular pencil $A-\lambda B$, it is possible to reduce $A$ to quasi-triangular form. In doing this it is necessary to monitor $A$ 's subdiagonal and $B$ 's diagonal in order to bring about decoupling whenever possible. The complete process, due to Moler and Stewart (1973), is as follows:

Algorithm 7.7.3 Given $A \in \mathbb{R}^{n \times n}$ and $B \in \mathbb{R}^{n \times n}$, the following algorithm computes orthogonal $Q$ and $Z$ such that $Q^{T} A Z=T$ is upper quasi-triangular and $Q^{T} B Z=S$ is upper triangular. $A$ is overwritten by $T$ and $B$ by $S$.

Using Algorithm 7.7.1, overwrite $A$ with $Q^{T} A Z$ (upper Hessenberg) and $B$ with $Q^{T} B Z$ (upper triangular).
until $q=n$
Set to zero subdiagonal entries that satisfy $\left|a_{i, i-1}\right| \leq \epsilon\left(\left|a_{i-1, i-1}\right|+\left|a_{i i}\right|\right)$.
Find the largest nonnegative $q$ and the smallest nonnegative $p$ such that if

$$
A=\begin{array}{ccc}
{\left[\begin{array}{ccc}
A_{11} & A_{12} & A_{13} \\
0 & A_{22} & A_{23} \\
0 & 0 & A_{33}
\end{array}\right]} & \begin{array}{c}
p \\
p \\
n-p-q \\
q
\end{array}
\end{array}
$$

then $A_{33}$ is upper quasi-triangular and $A_{22}$ is upper Hessenberg and unreduced.

Partition $B$ conformably:

$$
B={ }_{\left[\begin{array}{ccc}
B_{11} & B_{12} & B_{13} \\
0 & B_{22} & B_{23} \\
0 & 0 & B_{33}
\end{array}\right]}^{\underset{p}{n-p-q}}{ }_{q}{ }_{q}
$$

if $q<n$
if $B_{22}$ is singular
Zero $a_{n-q, n-q-1}$
else
Apply Algorithm 7.7.2 to $A_{22}$ and $B_{22}$ and update:

$$
\begin{aligned}
& A=\operatorname{diag}\left(I_{p}, Q, I_{q}\right)^{T} A \cdot \operatorname{diag}\left(I_{p}, Z, I_{q}\right) \\
& B=\operatorname{diag}\left(I_{p}, Q, I_{q}\right)^{T} B \cdot \operatorname{diag}\left(I_{p}, Z, I_{q}\right)
\end{aligned}
$$

end
end
end
This algorithm requires $30 n^{3}$ flops. If $Q$ is desired, an additional $16 n^{3}$ are necessary. If $Z$ is required, an additional $20 n^{3}$ are needed. These estimates of work are based on the experience that about two QZ iterations per eigenvalue arc necessary. Thus, the convergence properties of QZ are the same as for QR. The speed of the QZ algorithm is not affected by rank deficiency in $B$.

The computed $S$ and $T$ can be shown to satisfy

$$
Q_{0}^{T}(A+E) Z_{0}=T, \quad Q_{0}^{T}(B+F) Z_{0}=S,
$$

where $Q_{0}$ and $Z_{0}$ are exactly orthogonal and $\|E\|_{2} \approx \mathbf{u}\|A\|_{2}$ and $\|F\|_{2} \approx \mathbf{u}\|B\|_{2}$.

### 7.7.8 Generalized Invariant Subspace Computations

Many of the invariant subspace computations discussed in §7.6 carry over to the generalized eigenvalue problem. For example, approximate eigenvectors can be found via inverse iteration:

```
$q^{(0)} \in \mathbb{C}^{n \times n}$ given.
for $k=1,2, \ldots$
    Solve $(A-\mu B) z^{(k)}=B q^{(k-1)}$.
    Normalize: $q^{(k)}=z^{(k)} /\left\|z^{(k)}\right\|_{2}$.
    $\lambda^{(k)}=\left[q^{(k)}\right]^{H} A q^{(k)} /\left[q^{(k)}\right]^{H} A q^{(k)}$
end
```

If $B$ is nonsingular, then this is equivalent to applying (7.6.1) with the matrix $B^{-1} A$. Typically, only a single iteration is required if $\mu$ is an approximate eigenvalue computed by the QZ algorithm. By inverse iterating with the Hessenberg-triangular pencil, costly accumulation of the $Z$-transformations during the QZ iteration can be avoided.

Corresponding to the notion of an invariant subspace for a single matrix, we have the notion of a deflating subspace for the pencil $A-\lambda B$. In particular, we say that a $k$-dimensional subspace $S \subseteq \mathbb{C}^{n}$ is deflating for the pencil $A-\lambda B$ if the subspace $\{A x+B y: x, y \in S\}$ has dimension $k$ or less. Note that if

$$
Q^{H} A Z=T, \quad Q^{H} B Z=S
$$

is a generalized Schur decomposition of $A-\lambda B$, then the columns of $Z$ in the generalized Schur decomposition define a family of deflating subspaces. Indeed, if

$$
Q=\left[q_{1}|\cdots| q_{n}\right], \quad Z=\left[z_{1}|\cdots| z_{n}\right]
$$

are column partitionings, then

$$
\begin{aligned}
& \operatorname{span}\left\{A z_{1}, \ldots, A z_{k}\right\} \subseteq \operatorname{span}\left\{q_{1}, \ldots, q_{k}\right\}, \\
& \operatorname{span}\left\{B z_{1}, \ldots, B z_{k}\right\} \subseteq \operatorname{span}\left\{q_{1}, \ldots, q_{k}\right\},
\end{aligned}
$$

for $k=1: n$. Properties of deflating subspaces and their behavior under perturbation are described in Stewart (1972).

### 7.7.9 A Note on the Polynomial Eigenvalue Problem

More general than the generalized eigenvalue problem is the polynomial eigenvalue problem. Here we are given matrices $A_{0}, \ldots, A_{d} \in \mathbb{C}^{n \times n}$ and determine $\lambda \in \mathbb{C}$ and $0 \neq x \in \mathbb{C}^{n}$ so that

$$
P(\lambda) x=0
$$

where the $\lambda$-matrix $P(\lambda)$ is defined by

$$
P(\lambda)=A_{0}+\lambda A_{1}+\cdots+\lambda^{d} A_{d} .
$$

We assume $A_{d} \neq 0$ and regard $d$ as the degree of $P(\lambda)$. The theory behind the polynomial eigenvalue problem is nicely developed in Lancaster (1966).

It is possible to convert (7.7.3) into an equivalent linear eigenvalue problem with larger dimension. For example, suppose $d=3$ and

$$
L(\lambda)=\left[\begin{array}{rrr}
0 & 0 & A_{0} \\
-I & 0 & A_{1} \\
0 & -I & A_{2}
\end{array}\right]+\lambda\left[\begin{array}{ccc}
I & 0 & 0 \\
0 & I & 0 \\
0 & 0 & A_{3}
\end{array}\right] .
$$

If

$$
L(\lambda)\left[\begin{array}{c}
u_{1} \\
u_{2} \\
x
\end{array}\right]=\left[\begin{array}{l}
0 \\
0 \\
0
\end{array}\right],
$$

then

$$
0=A_{0} x+\lambda u_{1}=A_{0}+\lambda\left(A_{1} x+\lambda u_{2}\right)=A_{0}+\lambda\left(A_{1} x+\lambda\left(A_{2}+\lambda A_{3}\right)\right) x=P(\lambda) x .
$$

In general, we say that $L(\lambda)$ is a linearization of $P(\lambda)$ if there are $d n$-by-dn $\lambda$-matrices $S(\lambda)$ and $T(\lambda)$, each with constant nonzero determinants, so that

$$
S(\lambda)\left[\begin{array}{cc}
P(\lambda) & 0 \\
0 & I_{(d-1) n}
\end{array}\right] T(\lambda)=L(\lambda)
$$

has unit degree. With this conversion, the $A-\lambda B$ methods just discussed can be applied to find the required eigenvalues and eigenvectors.

Recent work has focused on how to choose the $\lambda$-transformations $S(\lambda)$ and $T(\lambda)$ so that special structure in $P(\lambda)$ is reflected in $L(\lambda)$. See Mackey, Mackey, Mehl, and Mehrmann (2006). The idea is to think of (7.7.6) as a factorization and to identify the transformations that produce a properly structured $L(\lambda)$. To appreciate this solution framework it is necessary to have a facility with $\lambda$-matrix manipulation and to that end we briefly examine the $\lambda$-matrix transformations behind the above linearization. If

$$
P_{1}(\lambda)=A_{1}+\lambda A_{2}+\cdots+\lambda^{d-1} A_{d}
$$

then

$$
P(\lambda)=A_{0}+\lambda P_{1}(\lambda)
$$

and it is easy to verify that

$$
\left[\begin{array}{rr}
I_{n} & -\lambda I_{n} \\
0 & I_{n}
\end{array}\right]\left[\begin{array}{cc}
A_{0}+\lambda P_{1}(\lambda) & 0 \\
0 & I_{n}
\end{array}\right]\left[\begin{array}{rr}
0 & I_{n} \\
-I_{n} & P_{1}(\lambda)
\end{array}\right]=\left[\begin{array}{rr}
\lambda I_{n} & A_{0} \\
-I_{n} & P_{1}(\lambda)
\end{array}\right] .
$$

Notice that the transformation matrices have unit determinant and that the $\lambda$-matrix on the right-hand side has degree $d-1$. The process can be repeated. If

$$
P_{2}(\lambda)=A_{2}+\lambda A_{3}+\cdots+\lambda^{d-2} A_{d}
$$

then

$$
P_{1}(\lambda)=A_{1}+\lambda P_{2}(\lambda)
$$

and

$$
\begin{aligned}
{\left[\begin{array}{c|cr}
I_{n} & 0 & 0 \\
\hline 0 & I_{n} & -\lambda I_{n} \\
0 & 0 & I_{n}
\end{array}\right]\left[\begin{array}{c|cc}
\lambda I_{n} & A_{0} & 0 \\
\hline-I_{n} & P_{1}(\lambda) & 0 \\
0 & 0 & I_{n}
\end{array}\right] } & {\left[\begin{array}{c|cc}
I_{n} & 0 & 0 \\
\hline 0 & 0 & I_{n} \\
0 & -I_{n} & P_{2}(\lambda)
\end{array}\right]=} \\
& {\left[\begin{array}{c|cc}
\lambda I_{n} & 0 & A_{0} \\
\hline-I_{n} & \lambda I_{n} & A_{1} \\
0 & -I_{n} & P_{2}(\lambda)
\end{array}\right] . }
\end{aligned}
$$

Note that the matrix on the right has degree $d-2$. A straightforward induction argument can be assembled to establish that if the $d n$-by-dn matrices $S(\lambda)$ and $T(\lambda)$ are defined by
$S(\lambda)=\left[\begin{array}{ccccc}I_{n} & -\lambda I_{n} & 0 & \cdots & 0 \\ 0 & I_{n} & -\lambda I_{n} & & \vdots \\ 0 & & \ddots & \ddots & \\ \vdots & & & I_{n} & -\lambda I_{n} \\ 0 & 0 & \cdots & 0 & I_{n}\end{array}\right], \quad T(\lambda)=\left[\begin{array}{ccccc}0 & 0 & 0 & \cdots & I \\ -I_{n} & 0 & & & P_{1}(\lambda) \\ 0 & -I_{n} & \ddots & & \vdots \\ \vdots & & \ddots & \ddots & P_{d-2}(\lambda) \\ 0 & 0 & \cdots & -I_{n} & P_{d-1}(\lambda)\end{array}\right]$
where

$$
P_{k}(\lambda)=A_{k}+\lambda A_{k+1}+\cdots+\lambda^{d-k} A_{d}
$$

then

$$
S(\lambda)\left[\begin{array}{cc}
P(\lambda) & 0 \\
0 & I_{(d-1) n}
\end{array}\right] T(\lambda)=\left[\begin{array}{ccccc}
\lambda I_{n} & 0 & 0 & \cdots & A_{0} \\
-I_{n} & \lambda I_{n} & & & A_{1} \\
0 & -I_{n} & \ddots & & \vdots \\
\vdots & & \ddots & \lambda I_{n} & A_{d-2} \\
0 & 0 & \cdots & -I_{n} & A_{d-1}+\lambda A_{d}
\end{array}\right] .
$$

Note that, if we solve the linearized problem using the QZ algorithm, then $O\left((d n)^{3}\right)$ flops are required.

## Problems

P7.7.1 Suppose $A$ and $B$ are in $\mathbf{R}^{n \times n}$ and that

$$
U^{T} B V=\left[\begin{array}{cc}
D & 0 \\
0 & 0 \\
r & n-r
\end{array}\right]_{n-r}^{r}, \quad U=\left[\underset{r}{U_{1} \mid} \underset{n-r}{U_{2}}\right], \quad V=\left[\underset{r}{V_{1} \mid} \underset{n-r}{V_{2}}\right],
$$

is the SVD of $B$, where $D$ is $r$-by-r and $r=\operatorname{rank}(B)$. Show that if $\lambda(A, B)=\mathbb{C}$ then $U_{2}^{T} A V_{2}$ is singular.
P7.7.2 Suppose $A$ and $B$ are in $\mathbf{R}^{n \times n}$. Give an algorithm for computing orthogonal $Q$ and $Z$ such that $Q^{T} A Z$ is upper Hessenberg and $Z^{T} B Q$ is upper triangular.

P7.7.3 Suppose

$$
A=\left[\begin{array}{cc}
A_{11} & A_{12} \\
0 & A_{22}
\end{array}\right] \quad \text { and } \quad B=\left[\begin{array}{cc}
B_{11} & B_{12} \\
0 & B_{22}
\end{array}\right]
$$

with $A_{11}, B_{11} \in \mathbf{R}^{k \times k}$ and $A_{22}, B_{22} \in \mathbf{R}^{j \times j}$. Under what circumstances do there exist

$$
X=\left[\begin{array}{cc}
I_{k} & X_{12} \\
0 & I_{j}
\end{array}\right] \quad \text { and } \quad Y=\left[\begin{array}{cc}
I_{k} & Y_{12} \\
0 & I_{j}
\end{array}\right]
$$

so that $Y^{-1} A X$ and $Y^{-1} B X$ are both block diagonal? This is the generalized Sylvester equation problem. Specify an algorithm for the case when $A_{11}, A_{22}, B_{11}$, and $B_{22}$ are upper triangular. See Kågström (1994).
P7.7.4 Suppose $\mu \notin \lambda(A, B)$. Relate the eigenvalues and eigenvectors of $A_{1}=(A-\mu B)^{-1} A$ and $B_{1}=(A-\mu B)^{-1} B$ to the generalized eigenvalues and eigenvectors of $A-\lambda B$.
P7.7.5 What does the generalized Schur decomposition say about the pencil $A-\lambda A^{T}$ ? Hint: If $T \in \mathbf{R}^{n \times n}$ is upper triangular, then $\mathcal{E}_{n} T \mathcal{E}_{n}$ is lower triangular where $\mathcal{E}_{n}$ is the exchange permutation defined in §1.2.11.
P7.7.6 Prove that

$$
L_{1}(\lambda)=\left[\begin{array}{cccc}
A_{3}+\lambda A_{4} & A_{2} & A_{1} & A_{0} \\
-I_{n} & 0 & 0 & 0 \\
0 & -I_{n} & 0 & 0 \\
0 & 0 & -I_{n} & 0
\end{array}\right], \quad L_{2}(\lambda)=\left[\begin{array}{cccc}
A_{3}+\lambda A_{4} & -I_{n} & 0 & 0 \\
A_{2} & 0 & -I_{n} & 0 \\
A_{1} & 0 & 0 & -I_{n} \\
A_{0} & 0 & 0 & 0
\end{array}\right]
$$

are linearizations of

$$
P(\lambda)=A_{0}+\lambda A_{1}+\lambda^{2} A_{2}+\lambda^{3} A_{3}+\lambda^{4} A_{4}
$$

Specify the $\lambda$-matrix transformations that relate $\operatorname{diag}\left(P(\lambda), I_{3 n}\right)$ to both $L_{1}(\lambda)$ and $L_{2}(\lambda)$.

## Notes and References for §7.7

For background to the generalized eigenvalue problem we recommend Stewart(IMC), Stewart and Sun (MPT), and Watkins (MEP) and:
B. Kågström and A. Ruhe (1983). Matrix Pencils, Proceedings Pite Havsbad, 1982, Lecture Notes in Mathematics Vol. 973, Springer-Verlag, New York.
$Q Z$-related papers include:
C.B. Moler and G.W. Stewart (1973). "An Algorithm for Generalized Matrix Eigenvalue Problems," SIAM J. Numer. Anal. 10, 241-256.
L. Kaufman (1974). "The LZ Algorithm to Solve the Generalized Eigenvalue Problem," SIAM J. Numer. Anal. 11, 997-1024.
R.C. Ward (1975). "The Combination Shift QZ Algorithm," SIAM J. Numer. Anal. 12, 835-853.
C.F. Van Loan (1975). "A General Matrix Eigenvalue Algorithm," SIAM J. Numer. Anal. 12, 819-834.
L. Kaufman (1977). "Some Thoughts on the QZ Algorithm for Solving the Generalized Eigenvalue Problem," ACM Trans. Math. Softw. 3, 65-75.
R.C. Ward (1981). "Balancing the Generalized Eigenvalue Problem," SIAM J. Sci. Stat. Comput. 2, 141-152.
P. Van Dooren (1982). "Algorithm 590: DSUBSP and EXCHQZ: Fortran Routines for Computing Deflating Subspaces with Specified Spectrum," ACM Trans. Math. Softw. 8, 376-382.
K. Dackland and B. Kågström (1999). "Blocked Algorithms and Software for Reduction of a Regular Matrix Pair to Generalized Schur Form," ACM Trans. Math. Softw. 25, 425-454.
D.S. Watkins (2000). "Performance of the QZ Algorithm in the Presence of Infinite Eigenvalues," SIAM J. Matrix Anal. Applic. 22, 364-375.
B. Kågström, D. Kressner, E.S. Quintana-Orti, and G. Quintana-Orti (2008). "Blocked Algorithms for the Reduction to Hessenberg-Triangular Form Revisited," BIT 48, 563-584.
Many algorithmic ideas associated with the $A-\lambda I$ problem extend to the $A-\lambda B$ problem:
A. Jennings and M.R. Osborne (1977). "Generalized Eigenvalue Problems for Certain Unsymmetric Band Matrices," Lin. Alg. Applic. 29, 139-150.
V.N. Kublanovskaya (1984). "AB Algorithm and Its Modifications for the Spectral Problem of Linear Pencils of Matrices," Numer. Math. 43, 329-342.
Z. Bai, J. Demmel, and M. Gu (1997). "An Inverse Free Parallel Spectral Divide and Conquer Algorithm for Nonsymmetric Eigenproblems," Numer. Math. 76, 279-308.
G.H. Golub and Q. Ye (2000). "Inexact Inverse Iteration for Generalized Eigenvalue Problems," BIT 40, 671-684.
F. Tisseur (2001). "Newton's Method in Floating Point Arithmetic and Iterative Refinement of Generalized Eigenvalue Problems," SIAM J. Matrix Anal. Applic. 22, 1038-1057.
D. Lemonnier and P. Van Dooren (2006). "Balancing Regular Matrix Pencils," SIAM J. Matrix Anal. Applic. 28, 253-263.
R. Granat, B. Kågström, and D. Kressner (2007). "Computing Periodic Deflating Subspaces Associated with a Specified Set of Eigenvalues," BIT 47, 763-791.

The perturbation theory for the generalized eigenvalue problem is treated in:
G.W. Stewart (1972). "On the Sensitivity of the Eigenvalue Problem $A x=\lambda B x$," SIAM J. Numer. Anal. 9, 669-686.
G.W. Stewart (1973). "Error and Perturbation Bounds for Subspaces Associated with Certain Eigenvalue Problems," SIAM Review 15, 727-764.
G.W. Stewart (1975). "Gershgorin Theory for the Generalized Eigenvalue Problem $A x=\lambda B x$," Math. Comput. 29, 600-606.
A. Pokrzywa (1986). "On Perturbations and the Equivalence Orbit of a Matrix Pencil," Lin. Alg. Applic. 82, 99-121.
J. Sun (1995). "Perturbation Bounds for the Generalized Schur Decomposition," SIAM J. Matrix Anal. Applic. 16, 1328-1340.
R. Bhatia and R.-C. Li (1996). "On Perturbations of Matrix Pencils with Real Spectra. II," Math. Comput. 65, 637-645.
J.-P. Dedieu (1997). "Condition Operators, Condition Numbers, and Condition Number Theorem for the Generalized Eigenvalue Problem," Lin. Alg. Applic. 263, 1-24.
D.J. Higham and N.J. Higham (1998). "Structured Backward Error and Condition of Generalized Eigenvalue Problems," SIAM J. Matrix Anal. Applic. 20, 493-512.
R. Byers, C. He, and V. Mehrmann (1998). "Where is the Nearest Non-Regular Pencil?," Lin. Alg. Applic. 285, 81-105.
V. Frayss and V. Toumazou (1998). "A Note on the Normwise Perturbation Theory for the Regular Generalized Eigenproblem," Numer. Lin. Alg. 5, 1-10.
R.-C. Li (2003). "On Perturbations of Matrix Pencils with Real Spectra, A Revisit," Math. Comput. 72, 715-728.
S. Bora and V. Mehrmann (2006). "Linear Perturbation Theory for Structured Matrix Pencils Arising in Control Theory," SIAM J. Matrix Anal. Applic. 28, 148-169.
X.S. Chen (2007). "On Perturbation Bounds of Generalized Eigenvalues for Diagonalizable Pairs," Numer. Math. 107, 79-86.

The Kronecker structure of the pencil $A-\lambda B$ is analogous to Jordan structure of $A-\lambda I$ and it can provide useful information about the underlying application. Papers concerned with this important decomposition include:
J.H. Wilkinson (1978). "Linear Differential Equations and Kronecker's Canonical Form," in Recent Advances in Numerical Analysis, C. de Boor and G.H. Golub (eds.), Academic Press, New York, 231--265.
J.H. Wilkinson (1979). "Kronecker's Canonical Form and the QZ Algorithm," Lin. Alg. Applic. 28, 285-303.
P. Van Dooren (1979). "The Computation of Kronecker's Canonical Form of a Singular Pencil," Lin. Alg. Applic. 27, 103-140.
J.W. Demmel (1983). "The Condition Number of Equivalence Transformations that Block Diagonalize Matrix Pencils," SIAM J. Numer. Anal. 20, 599-610.
J.W. Demmel and B. Kågström (1987). "Computing Stable Eigendecompositions of Matrix Pencils," Linear Alg. Applic. 88/89, 139-186.
B. Kågström (1985). "The Generalized Singular Value Decomposition and the General $A-\lambda B$ Problem," BIT 24, 568-583.
B. Kågström (1986). "RGSVD: An Algorithm for Computing the Kronecker Structure and Reducing Subspaces of Singular $A-\lambda B$ Pencils," SIAM J. Sci. Stat. Comput. 7, 185-211.
J. Demmel and B. Kågström (1986). "Stably Computing the Kronecker Structure and Reducing Subspaces of Singular Pencils $A-\lambda B$ for Uncertain Data," in Large Scale Eigenvalue Problems, J. Cullum and R.A. Willoughby (eds.), North-Holland, Amsterdam.
T. Beelen and P. Van Dooren (1988). "An Improved Algorithm for the Computation of Kronecker's Canonical Form of a Singular Pencil," Lin. Alg. Applic. 105, 9-65.
E. Elmroth and B. Kågström(1996). "The Set of 2-by-3 Matrix Pencils - Kronecker Structures and Their Transitions under Perturbations," SIAM J. Matrix Anal. Applic. 17, 1-34.
A. Edelman, E. Elmroth, and B. Kågström (1997). "A Geometric Approach to Perturbation Theory of Matrices and Matrix Pencils Part I: Versal Defformations," SIAM J. Matrix Anal. Applic. 18, 653-692.
E. Elmroth, P. Johansson, and B. Kågström (2001). "Computation and Presentation of Graphs Displaying Closure Hierarchies of Jordan and Kronecker Structures," Num. Lin. Alg. 8, 381-399.

Just as the Schur decomposition can be used to solve the Sylvester equation problem $A_{1} X-X A_{2}=B$, the generalized Schur decomposition can be used to solve the generalized Sylvester equation problem where matrices $X$ and $Y$ are sought so that $A_{1} X-Y A_{2}=B_{1}$ and $A_{3} X-Y A_{4}=B_{2}$, see:
W. Enright and S. Serbin (1978). "A Note on the Efficient Solution of Matrix Pencil Systems," BIT 18, 276-81.
B. Kågström and L. Westin (1989). "Generalized Schur Methods with Condition Estimators for Solving the Generalized Sylvester Equation," IEEE Trans. Autom. Contr. AC-34, 745-751.
B. Kågström (1994). "A Perturbation Analysis of the Generalized Sylvester Equation ( $A R-L B, D R- L E)=(C, F), "$ SIAM J. Matrix Anal. Applic. 15, 1045-1060.
J.-G. Sun (1996). "Perturbation Analysis of System Hessenberg and Hessenberg-Triangular Forms," Lin. Alg. Applic. 241-3, 811-849.
B. Kågström and P. Poromaa (1996). "LAPACK-style Algorithms and Software for Solving the Generalized Sylvester Equation and Estimating the Separation Between Regular Matrix Pairs," ACM Trans. Math. Softw. 22, 78-103.
I. Jonsson and B. Kågström (2002). "Recursive Blocked Algorithms for Solving Triangular SystemsPart II: Two-sided and Generalized Sylvester and Lyapunov Matrix Equations," ACM Trans. Math. Softw. 28, 416-435.
R. Granat and B. Kågström (2010). "Parallel Solvers for Sylvester-Type Matrix Equations with Applications in Condition Estimation, Part I: Theory and Algorithms," ACM Trans. Math. Softw. 37, Article 32.

Rectangular generalized eigenvalue problems also arise. In this setting the goal is to reduce the rank of $A-\lambda B$, see:
G.W. Stewart (1994). "Perturbation Theory for Rectangular Matrix Pencils," Lin. Alg. Applic. 208/209, 297-301.
G. Boutry, M. Elad, G.H. Golub, and P. Milanfar (2005). "The Generalized Eigenvalue Problem for Nonsquare Pencils Using a Minimal Perturbation Approach," SIAM J. Matrix Anal. Applic. 27, 582-601.
D. Chu and G.H. Golub (2006). "On a Generalized Eigenvalue Problem for Nonsquare Pencils," SIAM J. Matrix Anal. Applic. 28, 770-787.

References for the polynomial eigenvalue problem include:
P. Lancaster (1966). Lambda-Matrices and Vibrating Systems, Pergamon Press, Oxford, U.K.
I. Gohberg, P. Lancaster, and L. Rodman (1982). Matrix Polynomials, Academic Press, New York.
F. Tisseur (2000). "Backward Error and Condition of Polynomial Eigenvalue Problems," Lin. Alg. Applic. 309, 339-361.
J.-P. Dedieu and F. Tisseur (2003). "Perturbation Theory for Homogeneous Polynomial Eigenvalue Problems," Lin. Alg. Applic. 358, 71-94.
N.J. Higham, D.S. Mackey, and F. Tisseur (2006). "The Conditioning of Linearizations of Matrix Polynomials," SIAM J. Matrix Anal. Applic. 28, 1005-1028.
D.S. Mackey, N. Mackey, C. Mehl, V. Mehrmann (2006). "Vector Spaces of Linearizations for Matrix Polynomials," SIAM J. Matrix Anal. Applic. 28, 971-1004.

The structured quadratic eigenvalue problem is discussed briefly in §8.7.9.

### 7.8 Hamiltonian and Product Eigenvalue Problems

Two structured unsymmetric eigenvalue problems are considered. The Hamiltonian matrix eigenvalue problem comes with its own special Schur decomposition. Orthogonal symplectic similarity transformations are used to bring about the required reduction. The product eigenvalue problem involves computing the eigenvalues of a product like $A_{1} A_{2}^{-1} A_{3}$ without actually forming the product or the designated inverses. For detailed background to these problems, see Kressner (NMGS) and Watkins (MEP).

### 7.8.1 Hamiltonian Matrix Eigenproblems

Hamiltonian and symplectic matrices are introduced in §1.3.10. Their 2-by-2 block structure provide a nice framework for practicing block matrix manipulation, see P1.3.2 and P2.5.4. We now describe some interesting eigenvalue problems that involve these matrices. For a given $n$, we define the matrix $J \in \mathbb{R}^{2 n \times 2 n}$ by

$$
J=\left[\begin{array}{cc}
0 & I_{n} \\
-I_{n} & 0
\end{array}\right]
$$

and proceed to work with the families of 2-by-2 block structured matrices that are displayed in Figure 7.8.1. We mention four important facts concerning these matrices.

| Family | Definition | What They Look Like |  |
| :--- | :--- | :--- | :--- |
| Hamiltonian | $J M=(J M)^{T}$ | $M=\left[\begin{array}{cc}A & G \\ F & -A^{T}\end{array}\right]$ | $G$ symmetric $F$ symmetric |
| Skew Hamiltonian | $J N=-(J N)^{T}$ | $N=\left[\begin{array}{cc}A & G \\ F & A^{T}\end{array}\right]$ | $G$ skew-symmetric $F$ skew-symmetric |
| Symplectic | $J S=S^{-T} J$ | $S=\left[\begin{array}{ll}S_{11} & S_{12} \\ S_{21} & S_{22}\end{array}\right]$ | $S_{11}^{T} S_{21}$ symmetric $S_{22}^{T} S_{12}$ symmetric $S_{11}^{T} S_{22}=I+S_{21}^{T} S_{12}$ |
| Orthogonal Symplectic | $J Q=Q J$ | $Q=\left[\begin{array}{rr}Q_{1} & Q_{2} \\ -Q_{2} & Q_{1}\end{array}\right]$ | $Q_{1}^{T} Q_{2}$ symmetric $I=Q_{1}^{T} Q_{1}+Q_{2}^{T} Q_{2}$ |

Figure 7.8.1. Hamiltonian and symplectic structures

(1) Symplectic similarity transformations preserve Hamiltonian structure:

$$
J\left(S^{-1} M S\right)=\left(J S^{-1} J^{T}\right)\left(J M J^{T}\right)(J S)=-S^{T} M^{T} S^{-T} J=\left(J\left(S^{-1} M S\right)\right)^{T}
$$

(2) The square of a Hamiltonian matrix is skew-Hamiltonian:

$$
J M^{2}=\left(J M J^{T}\right)(J M)=-M^{T}(J M)^{T}=-M^{2 T} J^{T}=-\left(J M^{2}\right)^{T} .
$$

(3) If $M$ is a Hamiltonian matrix and $\lambda \in \lambda(M)$, then $-\lambda \in \lambda(M)$ :

$$
M\left[\begin{array}{l}
u \\
v
\end{array}\right]=\lambda\left[\begin{array}{l}
u \\
v
\end{array}\right] \quad \Rightarrow \quad M^{T}\left[\begin{array}{c}
v \\
-u
\end{array}\right]=-\lambda\left[\begin{array}{c}
v \\
-u
\end{array}\right] .
$$

(4) If $S$ is symplectic and $\lambda \in \lambda(S)$, then $1 / \lambda \in \lambda(S)$ :

$$
S\left[\begin{array}{l}
u \\
v
\end{array}\right]=\lambda\left[\begin{array}{l}
u \\
v
\end{array}\right] \quad \Rightarrow \quad S^{T}\left[\begin{array}{c}
v \\
-u
\end{array}\right]=\frac{1}{\lambda}\left[\begin{array}{c}
v \\
-u
\end{array}\right] .
$$

Symplectic versions of Householder and Givens transformations have a prominanent role to play in Hamiltonian matrix computations. If $P=I_{n}-2 v v^{T}$ is a Householder matrix, then $\operatorname{diag}(P, P)$ is a symplectic orthogonal matrix. Likewise, if $G \in \mathbb{R}^{2 n \times 2 n}$ is a Givens rotation that involves planes $i$ and $i+n$, then $G$ is a symplectic orthogonal matrix. Combinations of these transformations can be used to introduce zeros. For example, a Householder-Givens-Householder sequence can do this:

$$
\left[\begin{array}{c}
\times \\
\times \\
\times \\
\times \\
\hline \times \\
\times \\
\times \\
\times
\end{array}\right] \xrightarrow{\operatorname{diag}\left(P_{1}, P_{1}\right)}\left[\begin{array}{c}
\times \\
\times \\
\times \\
\times \\
\hline \times \\
0 \\
0 \\
0
\end{array}\right] \xrightarrow{G_{1.5}}\left[\begin{array}{c}
\times \\
\times \\
\times \\
\times \\
\hline 0 \\
0 \\
0 \\
0
\end{array}\right] \xrightarrow{\operatorname{diag}\left(P_{2}, P_{2}\right)}\left[\begin{array}{l}
\times \\
0 \\
0 \\
0 \\
\hline 0 \\
0 \\
0 \\
0
\end{array}\right] .
$$

This kind of vector reduction can be sequenced to produce a constructive proof of a structured Schur decomposition for Hamiltonian matrices. Suppose $\lambda$ is a real eigenvalue of a Hamiltonian matrix $M$ and that $x \in \mathbb{R}^{2 n}$ is a unit 2 -norm vector with $M x=\lambda x$. If $Q_{1} \in \mathbb{R}^{2 n \times 2 n}$ is an orthogonal symplectic matrix and $Q_{1}^{T} x=e_{1}$, then it follows from $\left(Q_{1}^{T} M Q_{1}\right)\left(Q_{1}^{T} x\right)=\lambda\left(Q_{1}^{T} x\right)$ that

$$
Q_{1}^{T} M Q_{1}=\left[\begin{array}{cccc|cccc}
\lambda & \times & \times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times & \times & \times \\
\hline 0 & 0 & 0 & 0 & -\lambda & 0 & 0 & 0 \\
0 & \times & \times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times & \times & \times \\
0 & \times & \times & \times & \times & \times & \times & \times
\end{array}\right] .
$$

The "extra" zeros follow from the Hamiltonian structure of $Q_{1}^{T} M Q_{1}$. The process can be repeated on the 6-by-6 Hamiltonian submatrix defined by rows and columns 2-3-4-6-7-8. Together with the assumption that $M$ has no purely imaginary eigenvalues, it is possible to show that an orthogonal symplectic matrix $Q$ exists so that

$$
Q^{T} M Q=\left[\begin{array}{rr}
Q_{1} & Q_{2} \\
-Q_{2} & Q_{1}
\end{array}\right]^{T}\left[\begin{array}{cc}
A & F \\
G & -A^{T}
\end{array}\right]\left[\begin{array}{rr}
Q_{1} & Q_{2} \\
-Q_{2} & Q_{1}
\end{array}\right]=\left[\begin{array}{cc}
T & R \\
0 & -T^{T}
\end{array}\right]
$$

where $T \in \mathbb{R}^{n \times n}$ is upper quasi-triangular. This is the real Hamiltonian-Schur decomposition. See Paige and Van Loan (1981) and, for a more general version, Lin, Mehrmann, and Xu (1999).

One reason that the Hamiltonian eigenvalue problem is so important is its connection to the algebraic Ricatti equation

$$
G+X A+A^{T} X-X F X=0 .
$$

This quadratic matrix problem arises in optimal control and a symmetric solution is sought so that the eigenvalues of $A-F X$ are in the open left half plane. Modest assumptions typically ensure that $M$ has no eigenvalues on the imaginary axis and that the matrix $Q_{1}$ in (7.8.1) is nonsingular. If we compare ( 2,1 ) blocks in (7.8.1), then

$$
Q_{2}^{T} A Q_{1}-Q_{2}^{T} F Q_{2}+Q_{1}^{T} G Q_{1}+Q_{1}^{T} A^{T} Q_{2}=0 .
$$

It follows from $I_{n}=Q_{1}^{T} Q_{1}+Q_{2}^{T} Q_{2}$ that $X=Q_{2} Q_{1}^{-1}$ is symmetric and that it satisfies (7.8.2). From (7.8.1) it is easy to show that $A-F X=Q_{1} T Q_{1}^{-1}$ and so the eigenvalues of $A-F X$ are the eigenvalues of $T$. It follows that the desired solution to the algebraic Ricatti equation can be obtained by computing the real Hamiltonian-Schur decomposition and ordering the eigenvalues so that $\lambda(T)$ is in the left half plane.

How might the real Hamiltonian-Schur form be computed? One idea is to reduce $M$ to some condensed Hamiltonian form and then devise a structure-preserving QRiteration. Regarding the former task, it is easy to compute an orthogonal symplectic $U_{0}$ so that

$$
U_{0}^{T} M U_{0}=\left[\begin{array}{cc}
H & R \\
D & -H^{T}
\end{array}\right]
$$

where $H \in \mathbb{R}^{n \times n}$ is upper Hessenberg and $D$ is diagonal. Unfortunately, a structurepreserving QR iteration that maintains this condensed form has yet to be devised. This impasse prompts consideration of methods that involve the skew-Hamiltonian matrix $N=M^{2}$. Because the ( 2,1 ) block of a skew-Hamiltonian matrix is skew-symmetric, it has a zero diagonal. Symplectic similarity transforms preserve skew-Hamiltonian structure, and it is straightforward to compute an orthogonal symplectic matrix $V_{0}$ such that

$$
V_{0}^{T} M^{2} V_{0}=\left[\begin{array}{cc}
H & R \\
0 & H^{T}
\end{array}\right],
$$

where $H$ is upper Hessenberg. If $U^{T} H U=T$ is the real Schur form of $H$ and and $Q=V_{0} \cdot \operatorname{diag}(U, U)$, then

$$
Q^{T} M^{2} Q=\left[\begin{array}{cc}
T & U^{T} R U \\
0 & T^{T}
\end{array}\right]
$$

is the real skew-Hamiltonian Schur form. See Van Loan (1984). It does not follow that $Q^{T} M Q$ is in Schur-Hamiltonian form. Moreover, the quality of the computed small eigenvalues is not good because of the explicit squaring of $M$. However, these shortfalls can be overcome in an efficient numerically sound way, see Chu, Lie, and Mehrmann (2007) and the references therein. Kressner (NMSE, p. 175-208) and Watkins (MEP, p. 319-341) have in-depth treatments of the Hamiltonian eigenvalue problem.

### 7.8.2 Product Eigenvalue Problems

Using SVD and QZ, we can compute the eigenvalues of $A^{T} A$ and $B^{-1} A$ without forming products or inverses. The intelligent computation of the Hamiltonian-Schur decomposition involves a correspondingly careful handling of the product $M$-times- $M$. In this subsection we further develop this theme by discussing various product decompositions. Here is an example that suggests how we might compute the Hessenberg decomposition of

$$
A=A_{3} A_{2} A_{1}
$$

where $A_{1}, A_{2}, A_{3} \in \mathbb{R}^{n \times n}$. Instead of forming this product explicitly, we compute orthogonal $U_{1}, U_{2}, U_{3} \in \mathbb{R}^{n \times n}$ such that

$$
\begin{array}{ll}
U_{1}^{T} A_{3} U_{3}=H_{3} & \text { (upper Hessenberg), } \\
U_{3}^{T} A_{2} U_{2}=T_{2} & \text { (upper triangular), } \\
U_{2}^{T} A_{1} U_{1}=T_{1} & \text { (upper triangular). }
\end{array}
$$

It follows that

$$
U_{1}^{T} A U_{1}=\left(U_{1}^{T} A_{3} U_{3}\right)\left(U_{3}^{T} A_{2} U_{2}\right)\left(U_{2}^{T} A_{1} U_{1}\right)=H_{3} T_{2} T_{1}
$$

is upper Hessenberg. A procedure for doing this would start by computing the QR factorizations

$$
Q_{2}^{T} A_{1}=R_{1}, \quad Q_{3}^{T}\left(A_{2} Q_{2}\right)=R_{2}
$$

If $\tilde{A}_{3}=A_{3} Q_{3}$, then $A=\tilde{A}_{3} R_{2} R_{1}$. The next phase involves reducing $\tilde{A}_{3}$ to Hessenberg form with Givens transformations coupled with "bulge chasing" to preserve the triangular structures already obtained. The process is similar to the reduction of $A-\lambda B$ to Hessenberg-triangular form; see §7.7.4.

Now suppose we want to compute the real Schur form of $A$

$$
\begin{array}{ll}
Q_{1}^{T} A_{3} Q_{3}=T_{3} & \text { (upper quasi-triangular) }, \\
Q_{3}^{T} A_{2} Q_{2}=T_{2} & \text { (upper triangular), } \\
Q_{2}^{T} A_{1} Q_{1}=T_{1} & \text { (upper triangular), }
\end{array}
$$

where $Q_{1}, Q_{2}, Q_{3} \in \mathbb{R}^{n \times n}$ are orthogonal. Without loss of generality we may assume that $\left\{A_{3}, A_{2}, A_{1}\right\}$ is in Hessenberg-triangular-triangular form. Analogous to the QZ iteration, the next phase is to produce a sequence of converging triplets

$$
\left\{A_{3}^{(k)}, A_{2}^{(k)}, A_{1}^{(k)}\right\} \rightarrow\left\{T_{3}, T_{2}, T_{1}\right\}
$$

with the property that all the iterates are in Hessenberg-triangular-triangular form.

Product decompositions (7.8.5) and (7.8.6) can be framed as structured decompositions of block-cyclic 3-by-3 matrices. For example, if

$$
U=\left[\begin{array}{ccc}
U_{1} & 0 & 0 \\
0 & U_{2} & 0 \\
0 & 0 & U_{3}
\end{array}\right]
$$

then we have the following restatement of (7.8.5):

$$
U^{T}\left[\begin{array}{ccc}
0 & 0 & A_{3} \\
A_{1} & 0 & 0 \\
0 & A_{2} & 0
\end{array}\right] U=\left[\begin{array}{ccc}
0 & 0 & H_{3} \\
T_{1} & 0 & 0 \\
0 & T_{2} & 0
\end{array}\right]=\tilde{H} .
$$

Consider the zero-nonzero structure of this matrix for the case $n=4$ :

$$
\tilde{H}=\left[\begin{array}{cccc|cccc|cccc}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times & \times \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times & \times \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times & \times \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & \times \\
\hline \times & \times & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & \times & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & \times & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & \times & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
\hline 0 & 0 & 0 & 0 & \times & \times & \times & \times & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & \times & \times & \times & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & \times & \times & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & 0 & 0
\end{array}\right] .
$$

Using the perfect shuffle $\mathcal{P}_{34}$ (see §1.2.11) we also have

$$
\mathcal{P}_{34} \tilde{H} \mathcal{P}_{34}=\left[\begin{array}{ccc|ccc|ccc|ccc}
0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times \\
\times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 \\
0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 \\
\hline 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times \\
0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 \\
0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times & 0 \\
\hline 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 & \times \\
0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times & 0 \\
\hline 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 & \times \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \times & 0
\end{array}\right] .
$$

Note that this is a highly structured 12-by-12 upper Hessenberg matrix. This connection makes it possible to regard the product-QR iteration as a structure-preserving

QR iteration. For a detailed discussion about this connection and its implications for both analysis and computation, sec Kressner (NMSE, pp. 146-174) and Watkins(MEP, pp. 293-303). We mention that with the "technology" that has been developed, it is possible to solve product eigenvalue problems where the factor matrices that define $A$ are rectangular. Square nonsingular factors can also participate through their inverses, e.g., $A=A_{3} A_{2}^{-1} A_{1}$.

## Problems

P7.8.1 What can you say about the eigenvalues and eigenvectors of a symplectic matrix?
P7.8.2 Suppose $S_{1}, S_{2} \in \mathbf{R}^{n \times n}$ are both skew-symmetric and let $A=S_{1} S_{2}$. Show that the nonzero eigenvalues of $A$ are not simple. How would you compute these eigenvalues?
P7.8.3 Relate the eigenvalues and eigenvectors of

$$
A=\left[\begin{array}{cccc}
0 & A_{1} & 0 & 0 \\
0 & 0 & A_{2} & 0 \\
0 & 0 & 0 & A_{3} \\
A_{4} & 0 & 0 & 0
\end{array}\right] .
$$

to the eigenvalues and eigenvectors of $\tilde{A}=A_{1} A_{2} A_{3} A_{4}$. Assume that the diagonal blocks are square.

## Notes and References for §7.8

The books by Kressner(NMSE) and Watkins (MEP) have chapters on product eigenvalue problems and Hamiltonian eigenvalue problems. The sometimes bewildering network of interconnections that exist among various structured classes of matrices is clarified in:
A. Bunse-Gerstner, R. Bycrs, and V. Mehrmann (1992). "A Chart of Numerical Methods for Structured Eigenvalue Problems," SIAM J. Matrix Anal. Applic. 13, 419-453.

Papers concerned with the Hamiltonian Schur decomposition include:
A.J. Laub and K. Meyer (1974). "Canonical Forms for Symplectic and Hamiltonian Matrices," J. Celestial Mechanics 9, 213-238.
C.C. Paige and C. Van Loan (1981). "A Schur Decomposition for Hamiltonian Matrices," Lin. Alg. Applic. 41, 11-32.
V. Mehrmann (1991). Autonomous Linear Quadratic Control Problems, Theory and Numerical Solution, Lecture Notes in Control and Information Sciences No. 163, Springer-Verlag, Heidelberg.
W.-W. Lin, V. Mehrmann, and H. Xu (1999). "Canonical Forms for Hamiltonian and Symplectic Matrices and Pencils," Lin. Alg. Applic. 302/303, 469-533.

Various methods for Hamiltonian eigenvalue problems have been devised that exploit the rich underlying structure, see:
C. Van Loan (1984). "A Symplectic Method for Approximating All the Eigenvalues of a Hamiltonian Matrix," Lin. Alg. Applic. 61, 233-252.
R. Byers (1986) "A Hamiltonian QR Algorithm," SIAM J. Sci. Stat. Comput. 7, 212-229.
P. Benner, R. Byers, and E. Barth (2000). "Algorithm 800: Fortran 77 Subroutines for Computing the Eigenvalues of Hamiltonian Matrices. I: the Square-Reduced Method," ACM Trans. Math. Softw. 26, 49-77.
H. Fassbender, D.S. Mackey and N. Mackey (2001). "Hamilton and Jacobi Come Full Circle: Jacobi Algorithms for Structured Hamiltonian Eigenproblems," Lin. Alg. Applic. 332-4, 37-80.
D.S. Watkins (2006). "On the Reduction of a Hamiltonian Matrix to Hamiltonian Schur Form," ETNA 23, 141-157.
D.S. Watkins (2004). "On Hamiltonian and Symplectic Lanczos Processes," Lin. Alg. Applic. 385, 23-45.
D. Chu, X. Liu, and V. Mehrmann (2007). "A Numerical Method for Computing the Hamiltonian Schur Form," Numer. Math. 105, 375-412.

Generalized eigenvalue problems that involve Hamiltonian matrices also arise:
P. Benner, V. Mehrmann, and H. Xu (1998). "A Numerically Stable, Structure Preserving Method for Computing the Eigenvalues of Real Hamiltonian or Symplectic Pencils," Numer. Math. 78, 329-358.
C. Mehl (2000). "Condensed Forms for Skew-Hamiltonian/Hamiltonian Pencils," SIAM J. Matrix Anal. Applic. 21, 454-476.
V. Mehrmann and D.S. Watkins (2001). "Structure-Preserving Methods for Computing Eigenpairs of Large Sparse Skew-Hamiltonian/Hamiltonian Pencils," SIAM J. Sci. Comput. 22, 1905-1925.
P. Benner and R. Byers, V. Mehrmann, and H. Xu (2002). "Numerical Computation of Deflating Subspaces of Skew-Hamiltonian/Hamiltonian Pencils," SIAM J. Matrix Anal. Applic. 24, 165190.

Methods for symplectic eigenvalue problems are discussed in:
P. Benner, H. Fassbender and D.S. Watkins (1999). "SR and SZ Algorithms for the Symplectic (Butterfly) Eigenproblem," Lin. Alg. Applic. 287, 41-76.

The Golub-Kahan SVD algorithm that we discuss in the next chapter does not form $A^{T} A$ or $A A^{T}$ despite the rich connection to the Schur decompositions of those matrices. From that point on there has been an appreciation for the numerical dangers associated with explicit products. Here is a sampling of the literature:
C. Van Loan (1975). "A General Matrix Eigenvalue Algorithm," SIAM J. Numer. Anal. 12, 819-834.
M.T. Heath, A.J. Laub, C.C. Paige, and R.C. Ward (1986). "Computing the SVD of a Product of Two Matrices," SIAM J. Sci. Stat. Comput. 7, 1147-1159.
R. Mathias (1998). "Analysis of Algorithms for Orthogonalizing Products of Unitary Matrices," Num. Lin. Alg. 3, 125-145.
G. Golub, K. Solna, and P. Van Dooren (2000). "Computing the SVD of a General Matrix Product/Quotient," SIAM J. Matrix Anal. Applic. 22, 1-19.
D.S. Watkins (2005). "Product Eigenvalue Problems," SIAM Review 47, 3-40.
R. Granat and B. Kgstrom (2006). "Direct Eigenvalue Reordering in a Product of Matrices in Periodic Schur Form," SIAM J. Matrix Anal. Applic. 28, 285-300.

Finally we mention that there is a substantial body of work concerned with structured error analysis and structured perturbation theory for structured matrix problems, see:
F. Tisseur (2003). "A Chart of Backward Errors for Singly and Doubly Structured Eigenvalue Problems," SIAM J. Matrix Anal. Applic. 24, 877-897.
R. Byers and D. Kressner (2006). "Structured Condition Numbers for Invariant Subspaces," SIAM J. Matrix Anal. Applic. 28, 326-347.
M. Karow, D. Kressner, and F. Tisseur (2006). "Structured Eigenvalue Condition Numbers," SIAM J. Matrix Anal. Applic. 28, 1052-1068.

### 7.9 Pseudospectra

If the purpose of computing is insight, then it is easy to see why the well-conditioned eigenvector basis is such a valued commodity, for in many matrix problems, replacement of $A$ with its diagonalization $X^{-1} A X$ leads to powerful, analytic simplifications. However, the insight-through-eigensystem paradigm has diminished impact in problems where the matrix of eigenvectors is ill-conditioned or nonexistent. Intelligent invariant subspace computation as discussed in §7.6 is one way to address the shortfall; pseudospectra are another. In this brief section we discuss the essential ideas behind the theory and computation of pseudospectra. The central message is simple: if you are working with a nonnormal matrix, then a graphical pseudospectral analysis effectively tells you just how much to trust the eigenvalue/eigenvector "story."

A slightly awkward feature of our presentation has to do with the positioning of this section in the text. As we will see, SVD calculations are an essential part of the pseudospectra scene and we do not detail dense matrix algorithms for that important decomposition until the next chapter. However, it makes sense to introduce the pseudospectra concept here at the end of Chapter 7 while the challenges of the
unsymmetric eigenvalue problem are fresh in mind. Moreover, with this "early" foundation we can subsequently present various pseudospectra insights that concern the behavior of the matrix exponential (§9.3), the Arnoldi method for sparse unsymmetric eigenvalue problems (§10.5), and the GMRES method for sparse unsymmetric linear systems (§11.4).

For maximum generality, we investigate the pseudospectra of complex, nonnormal matrices. The definitive pseudospectra reference is Trefethen and Embree (SAP). Virtually everything we discuss is presented in greater detail in that excellent volume.

### 7.9.1 Motivation

In many settings, the eigenvalues of a matrix "say something" about an underlying phenomenon. For example, if

$$
A=\left[\begin{array}{cc}
\lambda_{1} & M \\
0 & \lambda_{2}
\end{array}\right], \quad M>0,
$$

then

$$
\lim _{k \rightarrow \infty}\left\|A^{k}\right\|_{2}=0
$$

if and only if $\left|\lambda_{1}\right|<1$ and $\left|\lambda_{2}\right|<1$. This follows from Lemma 7.3.1, a result that we needed to establish the convergence of the QR iteration. Applied to our 2-by-2 example, the lemma can be used to show that

$$
\left\|A^{k}\right\|_{2} \leq \frac{M}{\epsilon}(\rho(A)+\epsilon)^{k}
$$

for any $\epsilon>0$ where $\rho(A)=\max \left\{\left|\lambda_{1}\right|,\left|\lambda_{2}\right|\right\}$ is the spectral radius. By making $\epsilon$ small enough in this inequality, we can draw a conclusion about the asymptotic behavior of $A^{k}$ :

$$
\text { If } \rho(A)<1, \text { then asymptotically } A^{k} \text { converges to zero as } \rho(A)^{k} .
$$

However, while the eigenvalues adequately predict the limiting behavior of $\left\|A^{k}\right\|_{2}$, they do not (by themselves) tell us much about what is happening if $k$ is small. Indeed, if $\lambda_{1} \neq \lambda_{2}$, then using the diagonalization

$$
A=\left[\begin{array}{cc}
1 & M /\left(\lambda_{2}-\lambda_{1}\right) \\
0 & 1
\end{array}\right]\left[\begin{array}{cc}
\lambda_{1} & 0 \\
0 & \lambda_{2}
\end{array}\right]\left[\begin{array}{cc}
1 & M /\left(\lambda_{2}-\lambda_{1}\right) \\
0 & 1
\end{array}\right]^{-1}
$$

we can show that

$$
A^{k}=\left[\begin{array}{c|c}
\lambda_{1}^{k} & M \sum_{j=0}^{k-1} \lambda_{1}^{k-1-j} \lambda_{2}^{j} \\
\hline 0 & \lambda_{2}^{k}
\end{array}\right] .
$$

Consideration of the $(1,2)$ entry suggests that $A^{k}$ may grow before decay sets in. This is affirmed in Figure 7.9.1 where the size of $\left\|A^{k}\right\|_{2}$ is tracked for the example

$$
A=\left[\begin{array}{cc}
0.999 & 1000 \\
0.0 & 0.998
\end{array}\right]
$$

![](https://cdn.mathpix.com/cropped/59af40f2-d7a0-4751-8bf2-c2837a3551d5-82.jpg?height=679&width=1028&top_left_y=189&top_left_x=219)
Figure 7.9.1. $\left\|A^{k}\right\|_{2}$ can grow even if $\rho(A)<1$

Thus, it is perhaps better to augment (7.9.1) as follows:

If $\rho(A)<1$, then aymptotically $A^{k}$ converges to zero like $\rho(A)^{k}$.
However, $A^{k}$ may grow substantially before exponential decay sets in.

This example with its ill-conditioned eigenvector matrix displayed in (7.9.2), points to just why classical eigenvalue analysis is not so informative for nonnormal matrices. Ill-conditioned eigenvector bases create a discrepancy between how A behaves and how its diagonalization $X A X^{-1}$ behaves. Pseudospectra analysis and computation narrow this gap.

### 7.9.2 Definitions

The pseudospectra idea is a generalization of the eigenvalue idea. Whereas the spectrum $\Lambda(A)$ is the set of all $z \in \mathbb{C}$ that make $\sigma_{\text {min }}(A-\lambda I)$ zero, the $\epsilon$-pseudospectrum of a matrix $A \in \mathbb{C}^{n \times n}$ is the subset of the complex plane defined by

$$
\Lambda_{\epsilon}(A)=\left\{z \in \mathbb{C}: \sigma_{\min }(A-\lambda I) \leq \epsilon\right\} .
$$

If $\lambda \in \Lambda_{\epsilon}(A)$, then $\lambda$ is an $\epsilon$-pseudoeigenvalue of $A$. A unit 2 -norm vector $v$ that satisfies $\|(A-\lambda I) v\|_{2}=\epsilon$ is a corresponding $\epsilon$-pseudoeigenvector. Note that if $\epsilon$ is zero, then $\Lambda_{\epsilon}(A)$ is just the set of $A$ 's eigenvalues, i.e., $\Lambda_{0}(A)=\Lambda(A)$.

We mention that because of their interest in what pseudospectra say about general linear operators, Trefethen and Embree (2005) use a strict inequality in the definition (7.9.5). The distinction has no impact in the matrix case.

Equivalent definitions of $\Lambda_{\epsilon}(\cdot)$ include

$$
\Lambda_{\epsilon}(A)=\left\{z \in \mathbb{C}:\left\|(z I-A)^{-1}\right\|_{2} \geq \frac{1}{\epsilon}\right\}
$$

which highlights the resolvent $(z I-A)^{-1}$ and

$$
\Lambda_{\epsilon}(A)=\left\{z \in \mathbb{C}: z \in \Lambda(A+E),\|E\|_{2} \leq \epsilon\right\}
$$

which characterize pseudspectra as (traditional) eigenvalues of nearby matrices. The equivalence of these three definitions is a straightforward verification that makes use of Chapter 2 facts about singular values, 2 -norms, and matrix inverses. We mention that greater generality can be achieved in (7.9.6) and (7.9.7) by replacing the 2 -norm with an arbitrary matrix norm.

### 7.9.3 Display

The pseudospectrum of a matrix is a visible subset of the complex plane so graphical display has a critical role to play in pseudospectra analysis. The Matlab-based Eigtool system developed by Wright(2002) can be used to produce pseudospectra plots that are as pleasing to the eye as they are informative. Eigtool's pseudospectra plots are contour plots where each contour displays the $z$-values associated with a specified value of $\epsilon$. Since

$$
\epsilon_{1} \leq \epsilon_{2} \quad \Rightarrow \quad \Lambda_{\epsilon_{1}} \subseteq \Lambda_{\epsilon_{2}}
$$

the typical pseudospectral plot is basically a topographical map that depicts the function $f(z)=\sigma_{\min }(z I-A)$ in the vicinity of the eigenvalues.

We present three Eigtool-produced plots that serve as illuminating examples. The first involves the $n$-by- $n$ Kahan matrix $\mathrm{Kah}_{n}(s)$, e.g.,

$$
\operatorname{Kah}_{5}(s)=\left[\begin{array}{ccccc}
1 & -c & -c & -c & -c \\
0 & s & -s c & -s c & -s c \\
0 & 0 & s^{2} & -s^{2} c & -s^{2} c \\
0 & 0 & 0 & s^{3} & -s^{3} c \\
0 & 0 & 0 & 0 & s^{4}
\end{array}\right], \quad c^{2}+s^{2}=1
$$

Recall that we used these matrices in §5.4.3 to show that QR with column pivoting can fail to detect rank deficiency. The eigenvalues $\left\{1, s, s^{2}, \ldots, s^{n-1}\right\}$ of $\operatorname{Kah}_{n}(s)$ are extremely sensitive to perturbation. This is revealed by considering the $\epsilon=10^{-6}$ contour that is displayed in Figure 7.9.2 together with $\Lambda\left(\mathrm{Kah}_{n}(s)\right)$.

The second example is the Demmel matrix $\operatorname{Dem}_{n}(\beta)$, e.g.,

$$
\operatorname{Dem}_{5}(\beta)=-\left[\begin{array}{ccccc}
1 & \beta & \beta^{2} & \beta^{3} & \beta^{4} \\
0 & 1 & \beta & \beta^{2} & \beta^{3} \\
0 & 0 & 1 & \beta & \beta^{2} \\
0 & 0 & 0 & 1 & \beta \\
0 & 0 & 0 & 0 & 1
\end{array}\right] .
$$

![](https://cdn.mathpix.com/cropped/59af40f2-d7a0-4751-8bf2-c2837a3551d5-84.jpg?height=618&width=588&top_left_y=230&top_left_x=454)
Figure 7.9.2. $\quad \Lambda_{\epsilon}\left(\operatorname{Kah}_{30}(s)\right)$ with $s^{29}=0.1$ and contours for $\epsilon=10^{-2}, \ldots, 10^{-6}$

The matrix $\operatorname{Dem}_{n}(\beta)$ is defective and has the property that very small perturbations can move an original eigenvalue to a position that are relatively far out on the imaginary axis. See Figure 7.9.3. The example is used to illuminate the nearness-to-instability problem presented in P7.9.13.

![](https://cdn.mathpix.com/cropped/59af40f2-d7a0-4751-8bf2-c2837a3551d5-84.jpg?height=620&width=635&top_left_y=1415&top_left_x=438)
Figure 7.9.3. $\quad \Lambda_{\epsilon}\left(\operatorname{Dem}_{50}(\beta)\right)$ with $\beta^{49}=10^{8}$ and contours for $\epsilon=10^{-2}, \ldots, 10^{-6}$

The last example concerns the pseudospectra of the Matlab "Gallery(5)" matrix:

$$
G_{5}=\left[\begin{array}{rrrrr}
-9 & 11 & -21 & 63 & -252 \\
70 & -69 & 141 & -421 & 1684 \\
-575 & 575 & -1149 & 3451 & -13801 \\
3891 & -3891 & 7782 & -23345 & 93365 \\
1024 & -1024 & 2048 & -6144 & 24572
\end{array}\right] .
$$

Notice in Figure 7.9.4 that $\Lambda_{10^{-13.5}}\left(G_{5}\right)$ has five components. In general, it can be

![](https://cdn.mathpix.com/cropped/59af40f2-d7a0-4751-8bf2-c2837a3551d5-85.jpg?height=604&width=595&top_left_y=643&top_left_x=439)
Figure 7.9.4. $\quad \Lambda_{\epsilon}\left(G_{5}\right)$ with contours for $\epsilon=10^{-11.5}, 10^{-12}, \ldots, 10^{-13.5}, 10^{-14}$

shown that each connected component of $\Lambda_{\epsilon}(A)$ contains at least one eigenvalue of $A$.

### 7.9.4 Some Elementary Properties

Pseudospectra are subsets of the complex plane so we start with a quick summary of notation. If $S_{1}$ and $S_{2}$ are subsets of the complex plane, then their sum $S_{1}+S_{2}$ is defined by

$$
S_{1}+S_{2}=\left\{s: s=s_{1}+s_{2}, s_{1} \in S_{1}, s_{2} \in S_{2}\right\}
$$

If $S_{1}$ consists of a single complex number $\alpha$, then we write $\alpha+S_{2}$. If $S$ is a subset of the complex plane and $\beta$ is a complex number, then $\beta \cdot S$ is defined by

$$
\beta \cdot S=\{\beta z: z \in S\} .
$$

The disk of radius $\epsilon$ centered at the origin is denoted by

$$
\Delta_{\epsilon}=\{z:|z| \leq \epsilon\}
$$

Finally, the distance from a complex number $z_{0}$ to a set of complex numbers $S$ is defined by

$$
\operatorname{dist}\left(z_{0}, S\right)=\min \left\{\left|z_{0}-z\right|: z \in S\right\} .
$$

Our first result is about the effect of translation and scaling. For eigenvalues we have

$$
\Lambda(\alpha I+\beta A)=\alpha+\beta \cdot \Lambda(A) .
$$

The following theorem establishes an analogous result for pseudospectra.
Theorem 7.9.1. If $\alpha, \beta \in \mathbb{C}$ and $A \in \mathbb{C}^{n \times n}$, then $\Lambda_{\epsilon|\beta|}(\alpha I+\beta A)=\alpha+\beta \cdot \Lambda_{\epsilon}(A)$.
Proof. Note that

$$
\begin{aligned}
\Lambda_{\epsilon}(\alpha I+A) & =\left\{z:\left\|(z I-(\alpha I+A))^{-1}\right\| \geq 1 / \epsilon\right\} \\
& =\left\{z:\left\|((z-\alpha) I-A)^{-1}\right\| \geq 1 / \epsilon\right\} \\
& =\alpha+\left\{z-\alpha:\left\|((z-\alpha) I-A)^{-1}\right\| \geq 1 / \epsilon\right\} \\
& =\alpha+\left\{z:\left\|(z I-A)^{-1}\right\| \geq 1 / \epsilon\right\}=\Lambda_{\epsilon}(A)
\end{aligned}
$$

and

$$
\begin{aligned}
\Lambda_{\epsilon|\beta|}(\beta \cdot A) & =\left\{z:\left\|(z I-\beta A)^{-1}\right\| \geq 1 /|\beta| \epsilon\right\} \\
& \left.=\{z: \|(z / \beta) I-A)^{-1} \| \geq 1 / \epsilon\right\} \\
& \left.=\beta \cdot\{z / \beta: \|(z / \beta) I-A)^{-1} \| \geq 1 / \epsilon\right\} \\
& \left.=\beta \cdot\{z: \| z I-A)^{-1} \| \geq 1 / \epsilon\right\}=\beta \cdot \Lambda_{\epsilon}(A) .
\end{aligned}
$$

The theorem readily follows by composing these two results.

General similarity transforms preserve eigenvalues but not $\epsilon$-pseudoeigenvalues. However, a simple inclusion property holds in the pseudospectra case.

Theorem 7.9.2. If $B=X^{-1} A X$, then $\Lambda_{\epsilon}(B) \subseteq \Lambda_{\epsilon \kappa_{2}(X)}(A)$.
Proof. If $z \in \Lambda_{\epsilon}(B)$, then

$$
\frac{1}{\epsilon} \leq\left\|(z I-B)^{-1}\right\|=\left\|X^{-1}(z I-A)^{-1} X^{-1}\right\| \leq \kappa_{2}(X)\left\|(z I-A)^{-1}\right\|,
$$

from which the theorem follows.
Corollary 7.9.3. If $X \in \mathbb{C}^{n \times n}$ is unitary and $A \in \mathbb{C}^{n \times n}$, then $\Lambda_{\epsilon}\left(X^{-1} A X\right)=\Lambda_{\epsilon}(A)$.
Proof. The proof is left as an exercise.

The $\epsilon$-pseudospectrum of a diagonal matrix is the union of $\epsilon$-disks.
Theorem 7.9.4. If $D=\operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{n}\right)$, then $\Lambda_{\epsilon}(D)=\left\{\lambda_{1}, \ldots, \lambda_{n}\right\}+\Delta_{\epsilon}$.
Proof. The proof is left as an exercise.

Corollary 7.9.5. If $A \in \mathbb{C}^{n \times n}$ is normal, then $\Lambda_{\epsilon}(A)=\Lambda(A)+\Delta_{\epsilon}$.
Proof. Since $A$ is normal, it has a diagonal Schur form $Q^{H} A Q=\operatorname{diag}\left(\lambda_{1}, \ldots, \lambda_{n}\right)=D$ with unitary $Q$. The proof follows from Theorem 7.9.4.

If $T=\left(T_{i j}\right)$ is a 2-by-2 block triangular matrix, then $\Lambda(T)=\Lambda\left(T_{11}\right) \cup \Lambda\left(T_{22}\right)$. Here is the pseudospectral analog:

Theorem 7.9.6. If

$$
T=\left[\begin{array}{cc}
T_{11} & T_{12} \\
0 & T_{22}
\end{array}\right]
$$

with square diagonal blocks, then $\Lambda_{\epsilon}\left(T_{11}\right) \cup \Lambda_{\epsilon}\left(T_{22}\right) \subseteq \Lambda_{\epsilon}(T)$.
Proof. The proof is left as an exercise.

Corollary 7.9.7. If

$$
T=\left[\begin{array}{cc}
T_{11} & 0 \\
0 & T_{22}
\end{array}\right]
$$

with square diagonal blocks, then $\Lambda_{\epsilon}(T)=\Lambda_{\epsilon}\left(T_{11}\right) \cup \Lambda_{\epsilon}\left(T_{22}\right)$.
Proof. The proof is left as an exercise.

The last property in our gallery of facts connects the resolvant $\left(z_{0} I-A\right)^{-1}$ to the distance that separates $z_{0}$ from $\Lambda_{\epsilon}(A)$.

Theorem 7.9.8. If $z_{0} \in \mathbb{C}$ and $A \in \mathbb{C}^{n \times n}$, then

$$
\operatorname{dist}\left(z_{0}, \Lambda_{\epsilon}(A)\right) \geq \frac{1}{\left\|\left(z_{0} I-A\right)^{-1}\right\|_{2}}-\epsilon .
$$

Proof. For any $z \in \Lambda_{\epsilon}(A)$ we have from Corollary 2.4.4 and (7.9.6) that

$$
\epsilon \geq \sigma_{\min }(z I-A)=\sigma_{\min }\left(\left(z_{0} I-A\right)-\left(z-z_{0}\right) I\right) \geq \sigma_{\min }\left(z_{0} I-A\right)-\left|z-z_{0}\right|
$$

and thus

$$
\left|z-z_{0}\right| \geq \frac{1}{\left\|\left(z_{0} I-A\right)^{-1}\right\|}-\epsilon
$$

The proof is completed by minimizing over all $z \in \Lambda_{\epsilon(A)}$.

### 7.9.5 Computing Pseudospectra

The production of a pseudospectral contour plot such as those displayed above requires sufficiently accurate approximations of $\sigma_{\text {min }}(z I-A)$ on a grid that consists of (perhaps)

1000 's of $z$-values. As we will see in §8.6, the computation of the complete SVD of an $n$-by- $n$ dense matrix is an $O\left(n^{3}\right)$ endeavor. Fortunately, steps can be taken to reduce each grid point calculation to $O\left(n^{2}\right)$ or less by exploiting the following ideas:

1. Avoid SVD-type computations in regions where $\sigma_{\min }(z I-A)$ is slowly varying. See Gallestey (1998).
2. Exploit Theorem 7.9.6 by ordering the eigenvalues so that the invariant subspace associated with $\Lambda\left(T_{11}\right)$ captures the essential behavior of $(z I-A)^{-1}$. See Reddy, Schmid, and Henningson (1993).
3. Precompute the Schur decomposition $Q^{H} A Q=T$ and apply a $\sigma_{\text {min }}$ algorithm that is efficient for triangular matrices. See Lui (1997).

We offer a few comments on the last strategy since it has much in common with the condition estimation problem that we discussed in §3.5.4. The starting point is to recognize that since $Q$ is unitary,

$$
\sigma_{\min }(z I-A)=\sigma_{\min }(z I-T) .
$$

The triangular structure of the transformed problem makes it possible to obtain a satisfactory estimate of $\sigma_{\text {min }}(z I-A)$ in $O\left(n^{2}\right)$ flops. If $d$ is a unit 2 -norm vector and $(z I-T) y=d$, then it follows from the SVD of $z I-T$ that

$$
\sigma_{\min }(z I-T) \leq \frac{1}{\|y\|_{2}} .
$$

Let $u_{\text {min }}$ be a left singular vector associated with $\sigma_{\text {min }}(z I-T)$. If $d$ is has a significant component in the direction of $u_{\text {min }}$, then

$$
\sigma_{\min }(z I-T) \approx \frac{1}{\|y\|_{2}}
$$

Recall that Algorithm 3.5.1 is a cheap heuristic procedure that dynamically determines the right hand side vector $d$ so that the solution to a given triangular system is large in norm. This is tantamount to choosing $d$ so that it is rich in the direction of $u_{\text {min }}$. A complex arithmetic, 2 -norm variant of Algorithm 3.5.1 is outlined in P7.9.13. It can be applied to $z I-T$. The resulting $d$-vector can be refined using inverse iteration ideas, see Toh and Trefethen (1996) and §8.2.2. Other approaches are discussed by Wright and Trefethen (2001).

### 7.9.6 Computing the $\epsilon$-Pseudospectral Abscissa and Radius

The $\epsilon$-pseudospectral abscissa of a matrix $A \in \mathbb{C}^{n \times n}$ is the rightmost point on the boundary of $\Lambda_{\epsilon}$ :

$$
\alpha_{\epsilon}(A)=\max _{z \in \Lambda_{\epsilon}(A)} \operatorname{Re}(z) .
$$

Likewise, the $\epsilon$-pseudospectral radius is the point of largest magnitude on the boundary of $\Lambda_{\epsilon}$ :

$$
\rho_{\epsilon}(A)=\max _{z \in \Lambda_{\epsilon}(A)}|z| .
$$

These quantities arise in the analysis of dynamical systems and effective iterative algorithms for their estimation have been proposed by Burke, Lewis, and Overton (2003) and Mengi and Overton (2005). A complete presentation and analysis of their very clever optimization procedures, which build on the work of Byers (1988), is beyond the scope of the text. However, at their core they involve interesting intersection problems that can be reformulated as structured eigenvalue problems. For example, if $i \cdot r$ is an eigenvalue of the matrix

$$
M=\left[\begin{array}{cc}
i e^{i \theta} A^{H} & -\epsilon I \\
\epsilon I & i e^{-i \theta} A
\end{array}\right],
$$

then $\epsilon$ is a singular value of $A-r e^{i \theta} I$. To see this, observe that if

$$
\left[\begin{array}{cc}
i e^{i \theta} A^{H} & -\epsilon I \\
\epsilon I & i e^{-i \theta} A
\end{array}\right]\left[\begin{array}{l}
f \\
g
\end{array}\right]=i \cdot r\left[\begin{array}{l}
f \\
g
\end{array}\right],
$$

then

$$
\left(A-r e^{i \theta} I\right)^{H}\left(A-r e^{i \theta} I\right) g=\epsilon^{2} g .
$$

The complex version of the SVD (§2.4.4) says that $\epsilon$ is a singular value of $A-r e^{1 \theta} I$. It can be shown that if $i r_{\text {max }}$ is the largest pure imaginary eigenvalue of $M$, then

$$
\epsilon=\sigma_{\min }\left(A-r_{\max } e^{1 \theta} I\right)
$$

This result can be used to compute the intersection of the ray $\left\{r e^{i \theta}: R \geq 0\right\}$ and the boundary of $\Lambda_{\epsilon}(A)$. This computation is at the heart of computing the $\epsilon$-pseudospectral radius. See Mengi and Overton (2005).

### 7.9.7 Matrix Powers and the $\boldsymbol{\epsilon}$-Pseudospectral Radius

At the start of this section we used the example

$$
A=\left[\begin{array}{ll}
0.999 & 1000 \\
0.000 & 0.998
\end{array}\right]
$$

to show that $\left\|A^{k}\right\|_{2}$ can grow even though $\rho(A)<1$. This kind of transient behavior can be anticipated by the pseudospectral radius. Indeed, it can be shown that for any $\epsilon>0$,

$$
\sup _{k \geq 0}\left\|A^{k}\right\|_{2} \geq \frac{\rho_{\epsilon}(A)-1}{\epsilon} .
$$

See Trefethen and Embree (SAP, pp. 160-161). This says that transient growth will occur if there is a contour $\left\{z: \|(\| z I-A)^{-1}=1 / \epsilon\right\}$ that extends beyond the unit disk. For the above 2 -by- 2 example, if $\epsilon=10^{-8}$, then $\rho_{\epsilon}(A) \approx 1.0017$ and the inequality (7.9.11) says that for some $k,\left\|A^{k}\right\|_{2} \geq 1.7 \times 10^{5}$. This is consistent with what is displayed in Figure 7.9.1.

## Problems

P7.9.1 Show that the definitions (7.9.5), (7.9.6), and (7.9.7) are equivalent.
P7.9.2 Prove Corollary 7.9.3.
P7.9.3 Prove Theorem 7.9.4.
P7.9.4 Prove Theorem 7.9.6.
P7.9.5 Prove Corollary 7.9.7.
P7.9.6 Show that if $A, E \in \mathbb{C}^{n \times n}$, then $\Lambda_{\epsilon}(A+E) \subseteq \Lambda_{\epsilon+\|E\|_{2}}(A)$.
P7.9.7 Suppose $\sigma_{\text {min }}\left(z_{1} I-A\right)=\epsilon_{1}$ and $\sigma_{\text {min }}\left(z_{2} I-A\right)=\epsilon_{2}$. Prove that there exists a real number $\mu$ so that if $z_{3}=(1-\mu) z_{1}+\mu z_{2}$, then $\sigma_{\text {min }}\left(z_{3} I-A\right)=\left(\epsilon_{1}+\epsilon_{2}\right) / 2$ ?
P7.9.8 Suppose $A \in \mathbb{C}^{n \times n}$ is normal and $E \in \mathbb{C}^{n \times n}$ is nonnormal. State and prove a theorem about $\Lambda_{\epsilon}(A+E)$.
P7.9.9 Explain the connection between Theorem 7.9.2 and the Bauer-Fike Theorem (Theorem 7.2.2).
P7.9.10 Define the matrix $J \in \mathbf{R}^{2 n \times 2 n}$ by

$$
J=\left[\begin{array}{cc}
0 & I_{n} \\
-I_{n} & 0
\end{array}\right]
$$

(a) The matrix $H \in \mathbf{R}^{2 n \times 2 n}$ is a Hamiltonian matrix if $J^{T} H J=-H^{T}$. It is easy to show that if $H$ is Hamiltonian and $\lambda \in \Lambda(H)$, then $-\lambda \in \Lambda(H)$. Does it follow that if $\lambda \in \Lambda_{\epsilon}(H)$, then $-\lambda \in \Lambda_{\epsilon}(H)$ ? (b) The matrix $S \in \mathbf{R}^{2 n \times 2 n}$ is a symplectic matrix if $J^{T} S J=S^{-T}$. It is easy to show that if $S$ is symplectic and $\lambda \in \Lambda(S)$, then $1 / \lambda \in \Lambda(S)$. Does it follow that if $\lambda \in \Lambda_{\epsilon}(S)$, then $1 / \lambda \in \Lambda_{\epsilon}(S)$ ?
P7.9.11 Unsymmetric Toeplitz matrices tend to have very ill-conditioned eigensystems and thus have interesting pseudospectral properties. Suppose

$$
A=\left[\begin{array}{cccc}
0 & 1 & \cdots & 0 \\
\alpha & 0 & \ddots & \vdots \\
\vdots & \ddots & \ddots & 1 \\
0 & \cdots & \alpha & 0
\end{array}\right]
$$

(a) Construct a diagonal matrix $S$ so that $S^{-1} A S=B$ is symmetric and tridiagonal with 1 's on its subdiagonal and superdiagonal. (b) What can you say about the condition of $A$ 's eigenvector matrix?
P7.9.12 A matrix $A \in \mathbb{C}^{n \times n}$ is stable if all of its eigenvalues have negative real parts. Consider the problem of minimizing $\|E\|_{2}$ subject to the constraint that $A+E$ has an eigenvalue on the imaginary axis. Explain why this optimization problem is equivalent to minimizing $\sigma_{\min }(i r I-A)$ over all $r \in \mathbf{R}$. If $E_{*}$ is a minimizing $E$, then $\|E\|_{2}$ can be regarded as measure of $A$ 's nearness to instability. What is the connection between $A$ 's nearness to instability and $\alpha_{\epsilon}(A)$ ?
P7.9.13 This problem is about the cheap estimation of the minimum singular value of a matrix, a critical computation that is performed over an over again during the course of displaying the pseudospectrum of a matrix. In light of the discussion in §7.9.5, the challenge is to estimate the smallest singular value of an upper triangular matrix $U=T-z I$ where $T$ is the Schur form of $A \in \mathbf{R}^{n \times n}$. The condition estimation ideas of §3.5.4 are relevant. We want to determine a unit 2 -norm vector $d \in \mathbb{C}^{n}$ such that the solution to $U y=d$ has a large 2 -norm for then $\sigma_{\text {min }}(U) \approx 1 /\|y\|_{2}$. (a) Suppose

$$
U=\left[\begin{array}{cc}
u_{11} & u^{H} \\
0 & U_{1}
\end{array}\right] \quad y=\left[\begin{array}{l}
\tau \\
z
\end{array}\right] \quad d=\left[\begin{array}{c}
c \\
s d_{1}
\end{array}\right]
$$

where $u_{11}, \tau \in \mathbb{C}, u, z, d_{1} \in \mathbb{C}^{n-1}, U_{1} \in \mathbb{C}^{(n-1) \times(n-1)},\left\|d_{1}\right\|_{2}=1, U_{1} y_{1}=d_{1}$, and $c^{2}+s^{2}=1$. Give an algorithm that determines $c$ and $s$ so that if $U y=d$, then $\|y\|_{2}$ is as large as possible. Hint: This is a 2-by-2 SVD problem. (b) Using part (a), develop a nonrecursive method for estimating $\sigma_{\text {min }}(U(k: n, k: n))$ for $k=n:-1: 1$.

## Notes and References for §7.7

Besides Trefethen and Embree (SAP), the following papers provide a nice introduction to the pseudospectra idea:
M. Embree and L.N. Trefethen (2001). "Generalizing Eigenvalue Theorems to Pseudospectra Theorems," SIAM J. Sci. Comput. 23, 583-590.
L.N. Trefethen (1997). "Pseudospectra of Linear Operators," SIAM Review 39, 383-406.

For more details concerning the computation and display of pseudoeigenvalues, see:
S.C. Reddy, P.J. Schmid, and D.S. Henningson (1993). "Pseudospectra of the Orr-Sommerfeld Operator," SIAM J. Applic. Math. 53, 15-47.
S.-H. Lui (1997). "Computation of Pseudospectra by Continuation," SIAM J. Sci. Comput. 18, 565-573.
E. Gallestey (1998). "Computing Spectral Value Sets Using the Subharmonicity of the Norm of Rational Matrices," BIT, 38, 22-33.
L.N. Trefethen (1999). "Computation of Pseudospectra," Acta Numerica 8, 247-295.
T.G. Wright (2002). Eigtool, http://www.comlab.ox.ac.uk/pseudospectra/eigtool/.

Interesting extensions/generalizations/applications of the pseudospectra idea include:
L. Reichel and L.N. Trefethen (1992). "Eigenvalues and Pseudo-Eigenvalues of Toeplitz Matrices," Lin. Alg. Applic. 164-164, 153-185.
K-C. Toh and L.N. Trefethen (1994). "Pseudozeros of Polynomials and Pseudospectra of Companion Matrices," Numer. Math. 68, 403-425.
F. Kittaneh (1995). "Singular Values of Companion Matrices and Bounds on Zeros of Polynomials," SIAM J. Matrix Anal. Applic. 16, 333-340.
N.J. Higham and F. Tisseur (2000). "A Block Algorithm for Matrix 1-Norm Estimation, with an Application to 1-Norm Pseudospectra," SIAM J. Matrix Anal. Applic. 21, 1185-1201.
T.G. Wright and L.N. Trefethen (2002). "Pseudospectra of Rectangular matrices," IMA J. Numer. Anal. 22, 501-519.
R. Alam and S. Bora (2005). "On Stable Eigendecompositions of Matrices," SIAM J. Matrix Anal. Applic. 26, 830-848.

Pseudospectra papers that relate to the notions of controllability and stability of linear systems include:
J.V. Burke and A.S. Lewis. and M.L. Overton (2003). "Optimization and Pseudospectra, with Applications to Robust Stability," SIAM J. Matrix Anal. Applic. 25, 80-104.
J.V. Burke, A.S. Lewis, and M.L. Overton (2003). "Robust Stability and a Criss-Cross Algorithm for Pseudospectra," IMA J. Numer. Anal. 23, 359-375.
J.V. Burke, A.S. Lewis and M.L. Overton (2004). "Pseudospectral Components and the Distance to Uncontrollability," SIAM J. Matrix Anal. Applic. 26, 350-361.

The following papers are concerned with the computation of the numerical radius, spectral radius, and field of values:
C. He and G.A. Watson (1997). "An Algorithm for Computing the Numerical Radius," IMA J. Numer. Anal. 17, 329-342.
G.A. Watson (1996). "Computing the Numerical Radius" Lin. Alg. Applic. 234, 163-172.
T. Braconnier and N.J. Higham (1996). "Computing the Field of Values and Pseudospectra Using the Lanczos Method with Continuation," BIT 36, 422-440.
E. Mengi and M.L. Overton (2005). "Algorithms for the Computation of the Pseudospectral Radius and the Numerical Radius of a Matrix," IMA J. Numer. Anal. 25, 648-669.
N. Guglielmi and M. Overton (2011). "Fast Algorithms for the Approximation of the Pseudospectral Abscissa and Pseudospectral Radius of a Matrix," SIAM J. Matrix Anal. Applic. 32, 1166-1192.

For more insight into the behavior of matrix powers, see:
P. Henrici (1962). "Bounds for Iterates, Inverses, Spectral Variation, and Fields of Values of Nonnormal Matrices," Numer. Math.4, 24-40.
J. Descloux (1963). "Bounds for the Spectral Norm of Functions of Matrices," Numer. Math. 5, 185-90.
T. Ransford (2007). "On Pseudospectra and Power Growth," SIAM J. Matrix Anal. Applic. 29, 699-711.

As an example of what pseudospectra can tell us about highly structured matrices, see:
L. Reichel and L.N. Trefethen (1992). "Eigenvalues and Pseudo-eigenvalues of Toeplitz Matrices," Lin. Alg. Applic. 162/163/164, 153-186.

