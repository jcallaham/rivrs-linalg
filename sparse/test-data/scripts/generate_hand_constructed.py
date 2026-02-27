#!/usr/bin/env python3
"""Generate hand-constructed test matrices with exact LDL^T factorizations.

Each matrix is written as a Matrix Market .mtx file with a companion .json
file containing the exact factorization (L, D, P, inertia).

All factorizations are verified numerically: ||P^T A P - L D L^T|| / ||A|| < epsilon.
"""

import os
import sys
import numpy as np
from scipy.io import mmread

# Ensure this script's directory is on the path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from mtx_utils import write_symmetric_mtx, read_mtx_header
from metadata_utils import create_matrix_entry, merge_metadata, write_metadata, load_metadata
from factorization_utils import (
    write_factorization,
    verify_factorization,
    reconstruct_L,
    reconstruct_D,
)

OUTPUT_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "hand-constructed")
METADATA_PATH = os.path.join(os.path.dirname(SCRIPT_DIR), "metadata.json")


def _dense_to_sparse_lower(A):
    """Extract lower-triangle entries from a dense symmetric matrix.

    Returns rows, cols, vals (all 0-indexed).
    """
    n = A.shape[0]
    rows, cols, vals = [], [], []
    for i in range(n):
        for j in range(i + 1):
            if A[i, j] != 0.0:
                rows.append(i)
                cols.append(j)
                vals.append(float(A[i, j]))
    return rows, cols, vals


def _compute_ldlt_1x1(A):
    """Compute LDL^T factorization using 1x1 pivots only (no permutation).

    Works for matrices where all leading minors are nonzero.
    Returns L (unit lower triangular), D (diagonal matrix), P (identity permutation).
    """
    n = A.shape[0]
    A = A.astype(float).copy()
    L = np.eye(n)
    d = np.zeros(n)

    for j in range(n):
        d[j] = A[j, j] - sum(L[j, k] ** 2 * d[k] for k in range(j))
        for i in range(j + 1, n):
            L[i, j] = (A[i, j] - sum(L[i, k] * L[j, k] * d[k] for k in range(j))) / d[j]

    D = np.diag(d)
    P = list(range(n))
    return L, D, P


def _compute_ldlt_scipy(A):
    """Compute LDL^T factorization using scipy for matrices that need pivoting.

    Returns L, D, P (as permutation vector).
    """
    from scipy.linalg import ldl

    L_full, D_full, perm = ldl(A)
    return L_full, D_full, list(perm)


def _l_to_entries(L):
    """Convert a unit lower triangular matrix L to sparse entry list (0-indexed)."""
    n = L.shape[0]
    entries = []
    for i in range(n):
        for j in range(i):
            if abs(L[i, j]) > 1e-16:
                entries.append({"row": i, "col": j, "value": float(L[i, j])})
    return entries


def _d_to_blocks(D):
    """Convert a (block) diagonal matrix D to the d_blocks format."""
    n = D.shape[0]
    blocks = []
    i = 0
    while i < n:
        if i + 1 < n and abs(D[i, i + 1]) > 1e-16:
            # 2x2 block
            blocks.append({
                "size": 2,
                "values": [
                    [float(D[i, i]), float(D[i, i + 1])],
                    [float(D[i + 1, i]), float(D[i + 1, i + 1])],
                ],
            })
            i += 2
        else:
            blocks.append({"size": 1, "values": [float(D[i, i])]})
            i += 1
    return blocks


def _count_inertia(D):
    """Count the inertia (positive, negative, zero eigenvalues) of block diagonal D."""
    n = D.shape[0]
    pos, neg, zero = 0, 0, 0
    i = 0
    while i < n:
        if i + 1 < n and abs(D[i, i + 1]) > 1e-16:
            # 2x2 block — compute eigenvalues
            block = D[i:i+2, i:i+2]
            eigs = np.linalg.eigvalsh(block)
            for e in eigs:
                if e > 1e-14:
                    pos += 1
                elif e < -1e-14:
                    neg += 1
                else:
                    zero += 1
            i += 2
        else:
            v = D[i, i]
            if v > 1e-14:
                pos += 1
            elif v < -1e-14:
                neg += 1
            else:
                zero += 1
            i += 1
    return {"positive": pos, "negative": neg, "zero": zero}


def _write_matrix(name, A, L, D, P, comment="", notes="", extra_props=None):
    """Write a matrix and its factorization to .mtx + .json files.

    Returns (metadata_entry, passed, residual).
    """
    n = A.shape[0]
    rows, cols, vals = _dense_to_sparse_lower(A)
    nnz = len(rows)

    mtx_path = os.path.join(OUTPUT_DIR, f"{name}.mtx")
    json_path = os.path.join(OUTPUT_DIR, f"{name}.json")

    write_symmetric_mtx(mtx_path, n, rows, cols, vals, comment=comment)

    l_entries = _l_to_entries(L)
    d_blocks = _d_to_blocks(D)
    inertia = _count_inertia(D)

    write_factorization(json_path, name, P, l_entries, d_blocks, inertia, notes=notes)

    # Verify
    passed, residual = verify_factorization(A, P, L, D)

    # Determine properties
    is_pd = inertia["negative"] == 0 and inertia["zero"] == 0
    is_singular = inertia["zero"] > 0

    properties = {
        "symmetric": True,
        "positive_definite": is_pd,
        "indefinite": not is_pd and not is_singular,
        "difficulty": "trivial",
    }
    if is_singular:
        properties["singular"] = True
        properties["indefinite"] = True  # singular matrices are also indefinite

    if extra_props:
        properties.update(extra_props)

    entry = create_matrix_entry(
        name=name,
        source="hand-constructed",
        category="hand-constructed",
        path=f"hand-constructed/{name}.mtx",
        size=n,
        nnz=nnz,
        properties=properties,
        factorization_path=f"hand-constructed/{name}.json",
    )

    return entry, passed, residual


# ============================================================================
# T008: Arrow matrix generators
# ============================================================================

def generate_arrow(n, positive_definite=True):
    """Generate an n x n arrow matrix (dense first row/column + diagonal).

    Arrow matrix structure:
        A[0, :] = dense
        A[:, 0] = dense
        A[i, i] = diagonal (for i > 0)

    For PD: diagonal entries large enough to ensure positive definiteness.
    For indefinite: mixed-sign diagonal.
    """
    A = np.zeros((n, n))

    if positive_definite:
        # First row/column = 1.0, diagonal = n (ensures PD via diagonal dominance)
        for i in range(1, n):
            A[0, i] = 1.0
            A[i, 0] = 1.0
        A[0, 0] = float(n)
        for i in range(1, n):
            A[i, i] = float(n)
    else:
        # Indefinite: alternating diagonal signs
        for i in range(1, n):
            A[0, i] = 1.0
            A[i, 0] = 1.0
        A[0, 0] = float(n)
        for i in range(1, n):
            if i % 2 == 0:
                A[i, i] = float(n)
            else:
                A[i, i] = -float(n)

    # Use scipy for indefinite, simple LDL for PD
    if positive_definite:
        L, D, P = _compute_ldlt_1x1(A)
    else:
        L, D, P = _compute_ldlt_scipy(A)

    return A, L, D, P


def generate_all_arrow():
    """Generate all arrow matrices: 5-pd, 10-indef, 15-indef."""
    results = []

    A, L, D, P = generate_arrow(5, positive_definite=True)
    results.append(_write_matrix(
        "arrow-5-pd", A, L, D, P,
        comment="5x5 arrow matrix, positive definite",
        notes="Arrow matrix with diagonal dominance. LDL^T computed analytically.",
        extra_props={"structure": "arrow"},
    ))

    A, L, D, P = generate_arrow(10, positive_definite=False)
    results.append(_write_matrix(
        "arrow-10-indef", A, L, D, P,
        comment="10x10 arrow matrix, indefinite",
        notes="Arrow matrix with alternating diagonal signs. LDL^T via scipy.",
        extra_props={"structure": "arrow"},
    ))

    A, L, D, P = generate_arrow(15, positive_definite=False)
    results.append(_write_matrix(
        "arrow-15-indef", A, L, D, P,
        comment="15x15 arrow matrix, indefinite",
        notes="Arrow matrix with alternating diagonal signs. LDL^T via scipy.",
        extra_props={"structure": "arrow"},
    ))

    return results


# ============================================================================
# T009: Tridiagonal matrix generators
# ============================================================================

def generate_tridiagonal(n, positive_definite=True):
    """Generate an n x n tridiagonal matrix.

    For PD: A = tridiag(-1, 4, -1) — diagonally dominant.
    For indefinite: A = tridiag(1, 0.5*(-1)^i, 1) — alternating diagonal.
    """
    A = np.zeros((n, n))

    if positive_definite:
        for i in range(n):
            A[i, i] = 4.0
        for i in range(n - 1):
            A[i, i + 1] = -1.0
            A[i + 1, i] = -1.0
    else:
        for i in range(n):
            A[i, i] = 0.5 * ((-1) ** i)
        for i in range(n - 1):
            A[i, i + 1] = 1.0
            A[i + 1, i] = 1.0

    # PD tridiagonal has simple 1x1 LDL^T by recurrence
    if positive_definite:
        L, D, P = _compute_ldlt_1x1(A)
    else:
        L, D, P = _compute_ldlt_scipy(A)

    return A, L, D, P


def generate_all_tridiagonal():
    """Generate all tridiagonal matrices: 5-pd, 10-indef, 20-indef."""
    results = []

    A, L, D, P = generate_tridiagonal(5, positive_definite=True)
    results.append(_write_matrix(
        "tridiag-5-pd", A, L, D, P,
        comment="5x5 tridiagonal matrix, positive definite",
        notes="Tridiag(-1, 4, -1). LDL^T via recurrence d_i = a_i - b^2/d_{i-1}.",
        extra_props={"structure": "tridiagonal"},
    ))

    A, L, D, P = generate_tridiagonal(10, positive_definite=False)
    results.append(_write_matrix(
        "tridiag-10-indef", A, L, D, P,
        comment="10x10 tridiagonal matrix, indefinite",
        notes="Tridiag with alternating diagonal. LDL^T via scipy.",
        extra_props={"structure": "tridiagonal"},
    ))

    A, L, D, P = generate_tridiagonal(20, positive_definite=False)
    results.append(_write_matrix(
        "tridiag-20-indef", A, L, D, P,
        comment="20x20 tridiagonal matrix, indefinite",
        notes="Tridiag with alternating diagonal. LDL^T via scipy.",
        extra_props={"structure": "tridiagonal"},
    ))

    return results


# ============================================================================
# T010: Block diagonal matrix generator
# ============================================================================

def generate_block_diagonal(block_sizes=None, block_types=None):
    """Generate a block-diagonal matrix with independently factored blocks.

    Default: three 5x5 blocks — one PD, one indefinite, one with 2x2 pivots.
    """
    if block_sizes is None:
        block_sizes = [5, 5, 5]
    if block_types is None:
        block_types = ["pd", "indefinite", "2x2_pivots"]

    n = sum(block_sizes)
    A = np.zeros((n, n))
    L = np.eye(n)
    D = np.zeros((n, n))
    P = list(range(n))

    offset = 0
    for bs, bt in zip(block_sizes, block_types):
        block = _make_block(bs, bt)
        A[offset:offset + bs, offset:offset + bs] = block

        # Factor each block independently
        bl, bd, bp = _compute_ldlt_scipy(block)
        L[offset:offset + bs, offset:offset + bs] = bl
        D[offset:offset + bs, offset:offset + bs] = bd
        # Adjust permutation for offset
        for i, p in enumerate(bp):
            P[offset + i] = offset + p

        offset += bs

    return A, L, D, P


def _make_block(size, block_type):
    """Create a dense block of given size and type."""
    if block_type == "pd":
        # Diagonally dominant SPD
        B = np.eye(size) * (size + 1.0)
        for i in range(size - 1):
            B[i, i + 1] = 1.0
            B[i + 1, i] = 1.0
        return B
    elif block_type == "indefinite":
        # Indefinite with mixed eigenvalues
        B = np.zeros((size, size))
        for i in range(size):
            B[i, i] = 2.0 * ((-1) ** i)
        for i in range(size - 1):
            B[i, i + 1] = 0.5
            B[i + 1, i] = 0.5
        return B
    elif block_type == "2x2_pivots":
        # Matrix that requires 2x2 pivots (zero diagonal entries)
        B = np.zeros((size, size))
        for i in range(0, size - 1, 2):
            B[i, i + 1] = 3.0
            B[i + 1, i] = 3.0
        if size % 2 == 1:
            B[size - 1, size - 1] = 1.0
        return B
    else:
        raise ValueError(f"Unknown block type: {block_type}")


def generate_all_block_diagonal():
    """Generate block diagonal matrix: 15x15 with 3 blocks."""
    results = []

    A, L, D, P = generate_block_diagonal()
    results.append(_write_matrix(
        "block-diag-15", A, L, D, P,
        comment="15x15 block diagonal: 5x5 PD + 5x5 indefinite + 5x5 with 2x2 pivots",
        notes="Three independent blocks factored separately. Adapted from SPRAL test patterns (BSD-3).",
        extra_props={"structure": "block-diagonal"},
    ))

    return results


# ============================================================================
# T011: Bordered block diagonal generator
# ============================================================================

def generate_bordered_block(block_sizes=None, border_width=2):
    """Generate a bordered block diagonal matrix.

    Structure: block diagonal with a dense border (last border_width rows/columns
    coupled to all blocks), causing fill-in during factorization.
    Adapted from SPRAL gen_bordered_block_diag pattern (BSD-3).
    """
    if block_sizes is None:
        block_sizes = [6, 6, 6]  # 18 interior + 2 border = 20

    n_interior = sum(block_sizes)
    n = n_interior + border_width
    A = np.zeros((n, n))

    # Fill diagonal blocks
    offset = 0
    for bs in block_sizes:
        block = np.eye(bs) * 4.0
        for i in range(bs - 1):
            block[i, i + 1] = -1.0
            block[i + 1, i] = -1.0
        A[offset:offset + bs, offset:offset + bs] = block

        # Border coupling: connect each block to border rows
        for b in range(border_width):
            border_row = n_interior + b
            # Connect first and last element of block to border
            coupling = 0.5 / (b + 1)
            A[offset, border_row] = coupling
            A[border_row, offset] = coupling
            A[offset + bs - 1, border_row] = coupling
            A[border_row, offset + bs - 1] = coupling

        offset += bs

    # Border diagonal entries (large to maintain well-conditioning)
    for b in range(border_width):
        A[n_interior + b, n_interior + b] = 10.0

    # Factor with scipy (needs pivoting due to border)
    L, D, P = _compute_ldlt_scipy(A)

    return A, L, D, P


def generate_all_bordered_block():
    """Generate bordered block diagonal: 20x20."""
    results = []

    A, L, D, P = generate_bordered_block()
    results.append(_write_matrix(
        "bordered-block-20", A, L, D, P,
        comment="20x20 bordered block diagonal: 3x6 blocks + 2-wide border",
        notes="Adapted from SPRAL gen_bordered_block_diag (BSD-3). Border causes fill-in.",
        extra_props={"structure": "bordered-block"},
    ))

    return results


# ============================================================================
# T012: Stress-test matrix generators
# ============================================================================

def generate_delayed_pivot_stress(n=10):
    """Generate a matrix with zero diagonal, forcing maximum 2x2 pivots.

    The matrix has all zeros on the diagonal and nonzero off-diagonal entries,
    requiring every pivot to be a 2x2 Bunch-Kaufman pivot.
    """
    A = np.zeros((n, n))
    for i in range(0, n - 1, 2):
        A[i, i + 1] = float(i + 2)
        A[i + 1, i] = float(i + 2)
    if n % 2 == 1:
        # Odd size: last element needs a diagonal entry
        A[n - 1, n - 1] = 1.0

    L, D, P = _compute_ldlt_scipy(A)
    return A, L, D, P


def generate_fill_in_stress(n=10):
    """Generate a matrix with a dense column causing worst-case fill-in.

    Column 0 is dense (connected to all other rows). Remaining structure is diagonal.
    Eliminating column 0 creates fill-in in the entire trailing submatrix.
    """
    A = np.zeros((n, n))
    A[0, 0] = float(n) * 2.0
    for i in range(1, n):
        A[i, 0] = 1.0
        A[0, i] = 1.0
        A[i, i] = float(n)

    L, D, P = _compute_ldlt_1x1(A)
    return A, L, D, P


def generate_ill_conditioned(n=10, cond=1e12):
    """Generate a near-singular matrix with known exact factorization.

    Uses a diagonal matrix with geometrically spaced eigenvalues to achieve
    the target condition number. The smallest eigenvalue is 1/cond.
    """
    eigenvalues = np.logspace(0, -np.log10(cond), n)
    # Alternate signs for indefinite character
    for i in range(1, n, 2):
        eigenvalues[i] = -eigenvalues[i]

    A = np.diag(eigenvalues)
    L = np.eye(n)
    D = A.copy()
    P = list(range(n))

    return A, L, D, P


def generate_all_stress():
    """Generate all stress-test matrices."""
    results = []

    A, L, D, P = generate_delayed_pivot_stress(10)
    results.append(_write_matrix(
        "stress-delayed-pivots", A, L, D, P,
        comment="10x10 zero diagonal, forces maximum 2x2 pivots",
        notes="Zero diagonal stress test. All pivots must be 2x2.",
        extra_props={
            "structure": "general",
            "expected_delayed_pivots": "high",
        },
    ))

    A, L, D, P = generate_fill_in_stress(10)
    results.append(_write_matrix(
        "stress-fill-in", A, L, D, P,
        comment="10x10 dense first column, worst-case fill-in",
        notes="Dense column causes fill-in in trailing submatrix.",
        extra_props={"structure": "general"},
    ))

    A, L, D, P = generate_ill_conditioned(10, 1e12)
    results.append(_write_matrix(
        "stress-ill-cond", A, L, D, P,
        comment="10x10 ill-conditioned diagonal, cond ~ 1e12",
        notes="Diagonal with geometrically spaced eigenvalues, alternating signs.",
        extra_props={
            "structure": "general",
            "condition_info": "ill-conditioned, cond ~ 1e12",
        },
    ))

    return results


# ============================================================================
# T013: Degenerate case generators
# ============================================================================

def generate_zero_diagonal():
    """3x3 matrix with zero diagonal. Adapted from SPRAL simple_mat_zero_diag (BSD-3).

    A = [0  2  1]
        [2  0  0]
        [1  0  3]
    """
    A = np.array([
        [0.0, 2.0, 1.0],
        [2.0, 0.0, 0.0],
        [1.0, 0.0, 3.0],
    ])
    L, D, P = _compute_ldlt_scipy(A)
    return A, L, D, P


def generate_singular():
    """3x3 singular matrix. Adapted from SPRAL simple_sing_mat (BSD-3).

    A = [2  1  0]
        [1  2  1]
        [0  1  0]

    Third column is nearly linearly dependent, leading to near-singularity.
    """
    A = np.array([
        [2.0, 1.0, 0.0],
        [1.0, 2.0, 1.0],
        [0.0, 1.0, 0.0],
    ])
    L, D, P = _compute_ldlt_scipy(A)
    return A, L, D, P


def generate_trivial_1x1():
    """1x1 trivial matrix: A = [5.0]."""
    A = np.array([[5.0]])
    L = np.eye(1)
    D = np.array([[5.0]])
    P = [0]
    return A, L, D, P


def generate_trivial_2x2():
    """2x2 matrix requiring a 2x2 pivot.

    A = [0  3]
        [3  0]
    """
    A = np.array([
        [0.0, 3.0],
        [3.0, 0.0],
    ])
    L, D, P = _compute_ldlt_scipy(A)
    return A, L, D, P


def generate_all_degenerate():
    """Generate all degenerate case matrices."""
    results = []

    A, L, D, P = generate_zero_diagonal()
    results.append(_write_matrix(
        "zero-diagonal-3", A, L, D, P,
        comment="3x3 zero diagonal. From SPRAL simple_mat_zero_diag (BSD-3)",
        notes="Adapted from SPRAL simple_mat_zero_diag. Requires 2x2 pivoting.",
        extra_props={
            "structure": "general",
            "expected_delayed_pivots": "high",
        },
    ))

    A, L, D, P = generate_singular()
    results.append(_write_matrix(
        "singular-3", A, L, D, P,
        comment="3x3 singular matrix. From SPRAL simple_sing_mat (BSD-3)",
        notes="Adapted from SPRAL simple_sing_mat. Structurally rank-deficient.",
        extra_props={"structure": "general"},
    ))

    A, L, D, P = generate_trivial_1x1()
    results.append(_write_matrix(
        "trivial-1x1", A, L, D, P,
        comment="1x1 trivial matrix: A = [5]",
        notes="Trivial 1x1 case. Factorization is trivial: L=I, D=A.",
    ))

    A, L, D, P = generate_trivial_2x2()
    results.append(_write_matrix(
        "trivial-2x2", A, L, D, P,
        comment="2x2 matrix requiring 2x2 pivot: A = [0 3; 3 0]",
        notes="Requires 2x2 Bunch-Kaufman pivot. Zero diagonal.",
        extra_props={
            "structure": "general",
            "expected_delayed_pivots": "high",
        },
    ))

    return results


# ============================================================================
# T014: Main entry point
# ============================================================================

def main():
    """Generate all hand-constructed matrices, verify factorizations, update metadata."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Generating hand-constructed test matrices...")
    print("=" * 60)

    all_results = []
    all_entries = []
    failures = []

    generators = [
        ("Arrow matrices", generate_all_arrow),
        ("Tridiagonal matrices", generate_all_tridiagonal),
        ("Block diagonal", generate_all_block_diagonal),
        ("Bordered block diagonal", generate_all_bordered_block),
        ("Stress-test matrices", generate_all_stress),
        ("Degenerate cases", generate_all_degenerate),
    ]

    for group_name, generator in generators:
        print(f"\n{group_name}:")
        results = generator()
        for entry, passed, residual in results:
            status = "PASS" if passed else "FAIL"
            print(f"  {entry['name']}: {entry['size']}x{entry['size']}, "
                  f"nnz={entry['nnz']}, residual={residual:.2e} [{status}]")
            all_entries.append(entry)
            if not passed:
                failures.append((entry['name'], residual))

    # Update metadata.json
    print("\n" + "=" * 60)
    metadata = load_metadata(METADATA_PATH)
    metadata = merge_metadata(metadata, all_entries)
    write_metadata(metadata, METADATA_PATH)
    print(f"Updated {METADATA_PATH}")

    # Summary
    print(f"\nSummary:")
    print(f"  Total matrices generated: {len(all_entries)}")
    print(f"  Verification failures: {len(failures)}")
    if failures:
        print("  FAILED matrices:")
        for name, residual in failures:
            print(f"    {name}: residual = {residual:.2e}")
        return 1

    print("  All factorizations verified successfully!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
