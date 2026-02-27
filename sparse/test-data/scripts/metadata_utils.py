"""Metadata index utilities for the test matrix collection."""

import json
from datetime import datetime, timezone


def create_matrix_entry(name, source, category, path, size, nnz, properties, **kwargs):
    """Create a metadata entry dict matching the MetadataIndex schema.

    Args:
        name: Unique matrix identifier.
        source: Origin ("hand-constructed", "suitesparse", "spral-inspired", "paper-referenced").
        category: Classification ("hand-constructed", "easy-indefinite", "hard-indefinite", "positive-definite").
        path: Relative path from test-data/ to the .mtx file.
        size: Matrix dimension n.
        nnz: Number of stored nonzeros (lower triangle for symmetric).
        properties: Dict of MatrixProperties fields.
        **kwargs: Additional fields (factorization_path, paper_references, suitesparse_id, etc.)

    Returns:
        Dict matching the MatrixEntry schema from data-model.md.
    """
    entry = {
        "name": name,
        "source": source,
        "category": category,
        "path": path,
        "size": size,
        "nnz": nnz,
        "in_repo": kwargs.get("in_repo", True),
        "properties": properties,
        "paper_references": kwargs.get("paper_references", []),
        "reference_results": kwargs.get("reference_results", {}),
    }

    if "factorization_path" in kwargs:
        entry["factorization_path"] = kwargs["factorization_path"]
    if "download_command" in kwargs:
        entry["download_command"] = kwargs["download_command"]
    if "suitesparse_id" in kwargs:
        entry["suitesparse_id"] = kwargs["suitesparse_id"]
    if "suitesparse_meta" in kwargs:
        entry["suitesparse_meta"] = kwargs["suitesparse_meta"]

    return entry


def merge_metadata(existing_index, new_entries):
    """Merge new matrix entries into an existing metadata index.

    Entries are matched by name — existing entries are updated, new entries appended.
    Updates total_count and generated timestamp.

    Args:
        existing_index: The current metadata index dict.
        new_entries: List of new MatrixEntry dicts.

    Returns:
        Updated metadata index dict.
    """
    existing_by_name = {m["name"]: i for i, m in enumerate(existing_index["matrices"])}

    for entry in new_entries:
        if entry["name"] in existing_by_name:
            idx = existing_by_name[entry["name"]]
            existing_index["matrices"][idx] = entry
        else:
            existing_index["matrices"].append(entry)
            existing_by_name[entry["name"]] = len(existing_index["matrices"]) - 1

    existing_index["total_count"] = len(existing_index["matrices"])
    existing_index["generated"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    return existing_index


def write_metadata(index, filepath):
    """Write the metadata index to a JSON file with indent=2.

    Args:
        index: Metadata index dict.
        filepath: Output file path.
    """
    with open(filepath, "w") as f:
        json.dump(index, f, indent=2)
        f.write("\n")


def load_metadata(filepath):
    """Load a metadata index from a JSON file.

    Args:
        filepath: Path to metadata.json.

    Returns:
        Metadata index dict.
    """
    with open(filepath) as f:
        return json.load(f)
