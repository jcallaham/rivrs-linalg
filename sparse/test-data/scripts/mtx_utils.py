"""Matrix Market file I/O utilities for the test matrix collection."""


def write_symmetric_mtx(filepath, n, rows, cols, vals, comment=""):
    """Write a symmetric Matrix Market coordinate file (lower triangle only, 1-indexed).

    Args:
        filepath: Output file path.
        n: Matrix dimension (n x n).
        rows: Row indices (0-indexed).
        cols: Column indices (0-indexed).
        vals: Nonzero values.
        comment: Optional comment string.
    """
    # Filter to lower triangle only (row >= col) and deduplicate
    entries = {}
    for r, c, v in zip(rows, cols, vals):
        if r < c:
            r, c = c, r
        entries[(r, c)] = v

    nnz = len(entries)

    with open(filepath, "w") as f:
        f.write("%%MatrixMarket matrix coordinate real symmetric\n")
        if comment:
            for line in comment.strip().split("\n"):
                f.write(f"% {line}\n")
        f.write(f"{n} {n} {nnz}\n")
        for (r, c), v in sorted(entries.items()):
            f.write(f"{r + 1} {c + 1} {v:.15e}\n")


def read_mtx_header(filepath):
    """Read a Matrix Market file header and return (n, nnz).

    Args:
        filepath: Path to .mtx file.

    Returns:
        Tuple of (n, nnz) where n is the matrix dimension and nnz is the
        number of stored entries (lower triangle for symmetric).
    """
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith("%"):
                continue
            parts = line.split()
            n = int(parts[0])
            nnz = int(parts[2])
            return n, nnz
    raise ValueError(f"Could not parse header from {filepath}")
