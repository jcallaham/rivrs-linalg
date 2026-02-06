"""Exact factorization I/O and verification utilities."""

import json
import numpy as np


class _NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy types."""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def write_factorization(filepath, matrix_name, permutation, l_entries, d_blocks, inertia, notes=""):
    """Write an exact LDL^T factorization to a JSON companion file.

    Args:
        filepath: Output .json file path.
        matrix_name: Name of the parent matrix.
        permutation: Permutation vector P (list of ints, 0-indexed).
        l_entries: List of dicts {"row": int, "col": int, "value": float} (0-indexed, unit diagonal implicit).
        d_blocks: List of dicts {"size": 1|2, "values": [d] or [[d11,d21],[d21,d22]]}.
        inertia: Dict {"positive": int, "negative": int, "zero": int}.
        notes: Optional string describing derivation method.
    """
    data = {
        "matrix_name": matrix_name,
        "permutation": permutation,
        "l_entries": l_entries,
        "d_blocks": d_blocks,
        "inertia": dict(inertia),
        "notes": notes,
    }
    with open(filepath, "w") as f:
        json.dump(data, f, indent=2, cls=_NumpyEncoder)
        f.write("\n")


def load_factorization(filepath):
    """Load a factorization from a JSON file.

    Returns:
        Dict with keys: matrix_name, permutation, l_entries, d_blocks, inertia, notes.
    """
    with open(filepath) as f:
        return json.load(f)


def reconstruct_L(n, l_entries):
    """Reconstruct the dense lower triangular factor L from sparse entries.

    L has unit diagonal. l_entries are 0-indexed {row, col, value}.
    """
    L = np.eye(n)
    for entry in l_entries:
        L[entry["row"], entry["col"]] = entry["value"]
    return L


def reconstruct_D(n, d_blocks):
    """Reconstruct the dense block diagonal factor D.

    d_blocks: list of {"size": 1|2, "values": ...}
    """
    D = np.zeros((n, n))
    pos = 0
    for block in d_blocks:
        if block["size"] == 1:
            D[pos, pos] = block["values"][0]
            pos += 1
        elif block["size"] == 2:
            vals = block["values"]
            D[pos, pos] = vals[0][0]
            D[pos, pos + 1] = vals[0][1]
            D[pos + 1, pos] = vals[1][0]
            D[pos + 1, pos + 1] = vals[1][1]
            pos += 2
        else:
            raise ValueError(f"Unsupported block size: {block['size']}")
    if pos != n:
        raise ValueError(f"D blocks cover {pos} rows but matrix is {n}x{n}")
    return D


def verify_factorization(A_dense, P, L, D, tol=1e-10):
    """Verify that P^T A P = L D L^T holds numerically.

    Args:
        A_dense: Dense n x n symmetric matrix (numpy array).
        P: Permutation vector (list of ints, 0-indexed).
        L: Dense lower triangular factor (n x n, unit diagonal).
        D: Dense block diagonal factor (n x n).
        tol: Tolerance for residual norm check.

    Returns:
        Tuple of (passed: bool, residual_norm: float).
    """
    n = A_dense.shape[0]
    P_mat = np.zeros((n, n))
    for i, p in enumerate(P):
        P_mat[i, p] = 1.0

    # P^T A P
    PTAP = P_mat.T @ A_dense @ P_mat
    # L D L^T
    LDLT = L @ D @ L.T

    residual = np.linalg.norm(PTAP - LDLT)
    norm_A = max(np.linalg.norm(A_dense), 1e-15)
    relative = residual / norm_A

    return relative < tol, relative


def verify_factorization_from_file(A_dense, factorization_path, tol=1e-10):
    """Load a factorization file and verify against a dense matrix.

    Args:
        A_dense: Dense n x n symmetric matrix.
        factorization_path: Path to the .json factorization file.
        tol: Tolerance.

    Returns:
        Tuple of (passed: bool, residual_norm: float).
    """
    fact = load_factorization(factorization_path)
    n = A_dense.shape[0]

    P = fact["permutation"]
    L = reconstruct_L(n, fact["l_entries"])
    D = reconstruct_D(n, fact["d_blocks"])

    return verify_factorization(A_dense, P, L, D, tol)
