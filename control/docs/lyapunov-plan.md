

## Claude's algorithm analysis

Balanced truncation does not use general Sylvester equations. It uses Lyapunov equations — specifically the Cholesky factor
variant.

The call chain is: AB09AD → AB09AX → SB03OU (called twice)

SB03OU solves for the Cholesky factor U where X = U'U:

- Continuous: A'X + XA = -scale² · B'B  (where X = U'U)
- Discrete: A'XA - X = -scale² · B'B  (where X = U'U)

The two calls compute:

1. Controllability Gramian (LTRANS=.TRUE.): solves with (A, B) to get the Cholesky factor of the controllability Gramian
2. Observability Gramian (LTRANS=.FALSE.): solves with (A', C') to get the Cholesky factor of the observability Gramian

The Hankel singular values are then computed from the SVD of the product of these two Cholesky factors (via MB03UD), and truncation is performed based on the singular value magnitudes.

Relationship to Sylvester

The Lyapunov equation A'X + XA = C is technically a special case of the Sylvester equation AX + XB = C where B = A'. However, in practice you'd never want to use a general Sylvester solver for this because:

1. Symmetry: C and X are symmetric (positive semi-definite, in fact). A general Sylvester solver ignores this structure and does roughly twice the work.
2. Cholesky factor approach: SB03OU computes U directly where X = U'U, never forming the full Gramian X. This is numerically superior — it guarantees the result is positive semi-definite by construction, and uses half the storage.
3. Single Schur decomposition: Since B = A', you only need one Schur decomposition (of A), not two independent ones.

**Implications for rivrs-linalg**

This means the next algorithmic building block for model reduction would be a Lyapunov solver (equivalent to the SB03 family in SLICOT), specifically:

- A Cholesky factor Lyapunov solver (SB03OU equivalent) — the workhorse for balanced truncation
- Both continuous and discrete variants
- The algorithm takes A already in Schur form (SB03OU assumes this), so the Schur decomposition infrastructure you already have in utils.rs via nalgebra would be reused

The existing Sylvester solver code would share some infrastructure (Schur decomposition, back-transformation) but the triangular solve phase would be different — it exploits the symmetry to work with triangular factors rather than full matrices.

## Key references

- **Continuous case:** S.J. Hammarling, "Numerical Solution of the Stable, Non-negative Definite Lyapunov Equation," IMA Journal of Numerical Analysis, 2(3):303–323, 1982.
- **Discrete case:** S.J. Hammarling, "Numerical Solution of the Discrete-Time, Convergent, Non-negative Definite Lyapunov Equation," Systems & Control Letters, 17(2):137–139, 1991.
