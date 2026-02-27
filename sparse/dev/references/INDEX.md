# Paper Index -- Phase 0.1 SSIDS Reference Library

This index catalogs the 14 academic papers compiled for the Phase 0.1 literature
review of the SSIDS (Sparse Symmetric Indefinite Direct Solver) project. These
papers provide the theoretical and algorithmic foundations for implementing a
clean-room SSIDS solver in Rust using faer infrastructure.

Papers are organized by category reflecting their role in the SSIDS
implementation. Each entry includes a full bibliographic citation, category
assignment, and a one-sentence relevance summary.

---

## Paper Table

### Core APTP

These papers define the A Posteriori Threshold Pivoting algorithm that is the
central innovation of the SSIDS solver.

| File | Citation | Relevance |
|------|----------|-----------|
| `duff2020` | I. S. Duff, J. D. Hogg, and F. Lopez. "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting." *SIAM J. Sci. Comput.*, 42(2):C23--C42, 2020. DOI: 10.1137/18M1225963 | Defines the APTP algorithm and its task-based multifrontal implementation in SSIDS v2.0, serving as the primary algorithmic reference for this project. |
| `hogg2016` | J. D. Hogg, E. Ovtchinnikov, and J. A. Scott. "A Sparse Symmetric Indefinite Direct Solver for GPU Architectures." *ACM Trans. Math. Softw.*, 42(1), Article 1, 2016. DOI: 10.1145/2756548 | Describes the original GPU-based SSIDS v1.0 with threshold partial pivoting, providing the architectural foundation and multifrontal design that APTP builds upon. |

### Multifrontal Foundations

These papers establish the multifrontal method that underlies the SSIDS solver
structure.

| File | Citation | Relevance |
|------|----------|-----------|
| `duff1984` | I. S. Duff and J. K. Reid. "The Multifrontal Solution of Unsymmetric Sets of Linear Equations." *SIAM J. Sci. Statist. Comput.*, 5(3):633--641, 1984. | Extends the multifrontal method from symmetric to unsymmetric systems, establishing the assembly tree and frontal matrix concepts used throughout modern sparse solvers. |
| `liu1992` | J. W. H. Liu. "The Multifrontal Method for Sparse Matrix Solution: Theory and Practice." *SIAM Review*, 34(1):82--109, 1992. | Provides the definitive survey and formal definitions of frontal matrices, update matrices, and assembly trees from a pure matrix perspective, independent of finite-element context. |
| `davis2016` | T. A. Davis, S. Rajamanickam, and W. M. Sid-Lakhdar. "A Survey of Direct Methods for Sparse Linear Systems." *Acta Numerica*, 25:165--421, 2016. | Comprehensive modern survey of all sparse direct methods including supernodal, multifrontal, ordering, and symbolic analysis -- the primary reference for understanding the full algorithmic landscape. |

### Pivoting Strategies

These papers address pivoting techniques for symmetric indefinite systems,
providing alternative approaches that inform APTP design decisions.

| File | Citation | Relevance |
|------|----------|-----------|
| `schenk2006` | O. Schenk and K. Gartner. "On Fast Factorization Pivoting Methods for Sparse Symmetric Indefinite Systems." *Electronic Trans. Numer. Anal. (ETNA)*, 23:158--179, 2006. | Introduces the Supernode-Bunch-Kaufman (SBK) pivoting with perturbation used in PARDISO, serving as the primary competing approach against which APTP demonstrates superior numerical robustness. |
| `duff2005` | I. S. Duff and S. Pralet. "Strategies for Scaling and Pivoting for Sparse Symmetric Indefinite Problems." *SIAM J. Matrix Anal. Appl.*, 27(1):312--340, 2005. | Develops symmetric weighted matching and constrained ordering strategies for preselecting 1x1 and 2x2 pivots, directly applicable to SSIDS preprocessing. |

### Ordering and Analysis

These papers cover fill-reducing orderings and symbolic analysis algorithms
essential to the analyse phase of SSIDS.

| File | Citation | Relevance |
|------|----------|-----------|
| `george1973` | A. George. "Nested Dissection of a Regular Finite Element Mesh." *SIAM J. Numer. Anal.*, 10(2):345--363, 1973. | Introduces nested dissection ordering, proving O(n^3) operation count optimality for 2D meshes and establishing the recursive separator approach used by METIS and other modern orderers. |
| `george1989` | A. George and J. W. H. Liu. "The Evolution of the Minimum Degree Ordering Algorithm." *SIAM Review*, 31(1):1--19, 1989. | Traces the development of minimum degree ordering and its enhancements, documenting the algorithm that remains a primary fill-reducing heuristic alongside nested dissection. |
| `gilbert1992` | J. R. Gilbert, C. Moler, and R. Schreiber. "Sparse Matrices in MATLAB: Design and Implementation." *SIAM J. Matrix Anal. Appl.*, 13(1):333--356, 1992. | Establishes the compressed sparse column (CSC) storage format and sparse matrix operation design principles that faer and most modern sparse libraries follow. |
| `gilbert1994` | J. R. Gilbert, E. G. Ng, and B. W. Peyton. "An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization." *SIAM J. Matrix Anal. Appl.*, 15(4):1075--1091, 1994. | Presents near-linear-time algorithms for predicting factor nonzero counts from the elimination tree, critical for efficient storage allocation in the symbolic analysis phase. |
| `ng1993` | E. G. Ng and B. W. Peyton. "Block Sparse Cholesky Algorithms on Advanced Uniprocessor Computers." *SIAM J. Sci. Comput.*, 14(5):1034--1056, 1993. | Compares multifrontal and left-looking block Cholesky methods with supernodal techniques, providing the blocked factorization framework that SSIDS extends to the indefinite case. |

### Infrastructure

These papers address matrix preprocessing (matching, scaling, permutation)
that supports the factorization phase.

| File | Citation | Relevance |
|------|----------|-----------|
| `duff1999` | I. S. Duff and J. Koster. "The Design and Use of Algorithms for Permuting Large Entries to the Diagonal of Sparse Matrices." *SIAM J. Matrix Anal. Appl.*, 20(4):889--901, 1999. | Introduces bottleneck transversal algorithms for diagonal maximization, providing the foundation for MC64-style preprocessing used to improve pivot quality in indefinite solvers. |
| `duff2001` | I. S. Duff and J. Koster. "On Algorithms for Permuting Large Entries to the Diagonal of a Sparse Matrix." *SIAM J. Matrix Anal. Appl.*, 22(4):973--996, 2001. | Extends diagonal maximization to weighted bipartite matching with scaling, completing the MC64 algorithm suite for matrix preprocessing before factorization. |

---

## RAL Technical Reports

The Duff, Hogg, and Lopez 2020 publication in *SIAM J. Sci. Comput.*
(`duff2020`) supersedes earlier RAL Technical Reports that described the SSIDS
v2.0 solver and APTP algorithm during development. The journal version
incorporates revisions from peer review and is the authoritative reference for
the APTP algorithm. When citing SSIDS v2.0 methodology, always use the SIAM
journal publication (DOI: 10.1137/18M1225963) rather than the technical reports.

---

## Markdown Conversion Quality

All 14 papers were converted from PDF to Markdown. The PDF is always the
authoritative source. Below are notes on conversion quality for each file.

**Good quality** (text, equations, and structure well-preserved):
- `duff2020` -- Clean conversion; Algorithm 3.1 pseudocode is readable.
- `liu1992` -- Well-converted; formal definitions intact.
- `george1989` -- Clean; algorithmic figures preserved as code blocks.
- `gilbert1994` -- Good; mathematical notation clear.
- `duff1999` -- Good; algorithm pseudocode readable.
- `duff2001` -- Good; bipartite matching algorithms well-preserved.
- `duff1984` -- Good; tables and algorithm description clear.
- `davis2016` -- Good overall for a very long paper (2966 lines).

**Minor issues** (figures rendered as external image links, but text is intact):
- `hogg2016` -- 11 external image references (mathpix CDN); equations in text are fine.
- `schenk2006` -- 17 external image references; algorithm pseudocode readable.
- `duff2005` -- 22 external image references (most are charts/tables); text clear.
- `george1973` -- 14 external image references for figures; mathematical proofs intact.
- `ng1993` -- Contains a SIAM download watermark line; otherwise clean.
- `gilbert1992` -- 7 external image references for figures (spy plots, etc.); text fine.

All papers have their core algorithmic content (definitions, theorems,
algorithms, proofs) well-preserved in Markdown. Figures referenced as external
images will not render offline, but the surrounding text descriptions are
sufficient for understanding. Consult the PDF originals for any figure-dependent
analysis.

---

## faer Reference

The faer linear algebra library paper is not in this directory but is available
at `references/faer-rs/paper.md`:

> S. El Kazdadi. "faer: A linear algebra library for the Rust programming
> language." *Journal of Open Source Software (JOSS)*, 2023.

faer provides the core infrastructure for the SSIDS implementation: CSC sparse
matrix storage, elimination tree computation, AMD and COLAMD fill-reducing
orderings, sparse Cholesky factorization routines, permutation utilities, and
workspace management patterns. The faer sparse module (`faer::sparse`) is the
primary dependency for symbolic analysis and will be leveraged extensively
throughout development.

---

## Category Summary

| Category | Papers | Key Concepts |
|----------|--------|--------------|
| Core APTP | 2 | A posteriori threshold pivoting, fail-in-place, speculative execution |
| Multifrontal Foundations | 3 | Assembly trees, frontal/update matrices, elimination trees |
| Pivoting Strategies | 2 | SBK, Bunch-Kaufman, weighted matching, constrained ordering |
| Ordering and Analysis | 5 | Nested dissection, minimum degree, CSC storage, row/column counts, supernodes |
| Infrastructure | 2 | Diagonal maximization, bipartite matching, matrix scaling |
