#!/usr/bin/env python3
"""Download curated SuiteSparse matrices for SSIDS solver testing.

Uses ssgetpy to query and download matrices from the SuiteSparse Matrix Collection
in Matrix Market format. Matrices are organized by difficulty category.

Curated lists sourced from:
- Duff, Hogg, Lopez (2020) Tables 1-2
- Hogg, Ovtchinnikov, Scott (2016) Table III
- Duff & Pralet (2005) Tables 3.1-3.2
"""

import argparse
import os
import sys
import glob
import tarfile
import shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from mtx_utils import read_mtx_header
from metadata_utils import create_matrix_entry, merge_metadata, write_metadata, load_metadata

BASE_DIR = os.path.join(os.path.dirname(SCRIPT_DIR))
METADATA_PATH = os.path.join(BASE_DIR, "metadata.json")

LARGE_THRESHOLD_ROWS = 200000
# LARGE_THRESHOLD_ROWS = 99999999999

# ============================================================================
# Curated matrix lists (T016)
# ============================================================================

# Each entry: (Group/Name, properties_dict)
# Properties include paper references, difficulty, domain info

EASY_INDEFINITE = [
    ("Oberwolfach/t2dal", {
        "kind": "model reduction problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("GHS_indef/dixmaanl", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("Oberwolfach/rail_79841", {
        "kind": "model reduction problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("GHS_indef/dawson5", {
        "kind": "structural problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("Boeing/bcsstk39", {
        "kind": "structural problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("Boeing/pct20stif", {
        "kind": "structural problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("GHS_indef/copter2", {
        "kind": "computational fluid dynamics problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("Boeing/crystk03", {
        "kind": "structural problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("Koutsovasilis/F2", {
        "kind": "structural problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("Cunningham/qa8fk", {
        "kind": "acoustics problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("Oberwolfach/gas_sensor", {
        "kind": "model reduction problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("Oberwolfach/t3dh", {
        "kind": "model reduction problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("PARSEC/H2O", {
        "kind": "quantum chemistry problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("GHS_indef/sparsine", {
        "kind": "synthetic problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("Newman/astro-ph", {
        "kind": "undirected graph",
        "expected_delayed_pivots": "low",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("Andrianov/mip1", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("PARSEC/SiNa", {
        "kind": "quantum chemistry problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("INPRO/msdoor", {
        "kind": "structural problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("Boeing/pwtk", {
        "kind": "structural problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("BenElechi/BenElechi1", {
        "kind": "structural problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("PARSEC/Si10H16", {
        "kind": "quantum chemistry problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("PARSEC/Si5H12", {
        "kind": "quantum chemistry problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("GHS_indef/helm3d01", {
        "kind": "computational fluid dynamics problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Pralet (2005) Table 3.1"],
    }),
    ("Cote/vibrobox", {
        "kind": "structural problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Pralet (2005) Table 3.1"],
    }),
    ("GHS_indef/bloweybq", {
        "kind": "subsequent optimization problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Pralet (2005) Table 3.1"],
    }),
    ("GHS_indef/linverse", {
        "kind": "inverse problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Pralet (2005) Table 3.1"],
    }),
    ("GHS_indef/spmsrtls", {
        "kind": "least squares problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Duff, Pralet (2005) Table 3.1"],
    }),
    ("Boeing/crystk02", {
        "kind": "structural problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
    ("ND/nd3k", {
        "kind": "2D/3D problem",
        "expected_delayed_pivots": "low",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("Oberwolfach/filter3D", {
        "kind": "model reduction problem",
        "expected_delayed_pivots": "none",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 1"],
    }),
]

HARD_INDEFINITE = [
    ("TSOPF/TSOPF_FS_b39_c7", {
        "kind": "power network problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("TSOPF/TSOPF_FS_b162_c1", {
        "kind": "power network problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/cont-201", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "high",
        "killer_case": True,
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/stokes128", {
        "kind": "computational fluid dynamics problem",
        "expected_delayed_pivots": "high",
        "killer_case": True,
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/ncvxqp1", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/cont-300", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "high",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/bratu3d", {
        "kind": "nonlinear PDE problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/cvxqp3", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/d_pretok", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/turon_m", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/ncvxqp5", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/ncvxqp3", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "high",
        "killer_case": True,
        "paper_references": [
            "Hogg, Ovtchinnikov, Scott (2016) Table III",
            "Duff, Hogg, Lopez (2020) Table 2",
        ],
    }),
    ("GHS_indef/ncvxqp7", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("Schenk_IBMNA/c-big", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "high",
        "killer_case": True,
        "paper_references": ["Duff, Hogg, Lopez (2020) Table 2"],
    }),
    ("GHS_indef/c-71", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "high",
        "killer_case": True,
        "paper_references": [
            "Hogg, Ovtchinnikov, Scott (2016) Table III",
            "Duff, Hogg, Lopez (2020) Table 2",
        ],
    }),
    ("GHS_indef/aug3dcqp", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Pralet (2005) Table 3.2"],
    }),
    ("GHS_indef/blockqp1", {
        "kind": "optimization problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Pralet (2005) Table 3.2"],
    }),
    ("GHS_indef/mario001", {
        "kind": "computational fluid dynamics problem",
        "expected_delayed_pivots": "medium",
        "paper_references": ["Duff, Pralet (2005) Table 3.2"],
    }),
]

POSITIVE_DEFINITE = [
    ("GHS_psdef/crankseg_1", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("Rothberg/cfd2", {
        "kind": "computational fluid dynamics problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("DNVS/thread", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("DNVS/shipsec1", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("DNVS/shipsec8", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("Oberwolfach/boneS01", {
        "kind": "model reduction problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("GHS_psdef/crankseg_2", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("DNVS/shipsec5", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("GHS_psdef/bmwcra_1", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("DNVS/ship_003", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("Um/offshore", {
        "kind": "electromagnetics problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("ND/nd6k", {
        "kind": "2D/3D problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("ND/nd12k", {
        "kind": "2D/3D problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("Schenk_AFE/af_shell7", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("Schenk_AFE/af_0_k101", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("GHS_psdef/ldoor", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("GHS_psdef/inline_1", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("GHS_psdef/apache2", {
        "kind": "structural problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
    ("AMD/G3_circuit", {
        "kind": "circuit simulation problem",
        "paper_references": ["Hogg, Ovtchinnikov, Scott (2016) Table III"],
    }),
]


# ============================================================================
# Download logic (T017)
# ============================================================================

def download_matrix(group_name, category, dest_base, max_rows=None, dry_run=False):
    """Download a single matrix from SuiteSparse using ssgetpy.

    Args:
        group_name: "Group/Name" identifier (e.g., "GHS_indef/ncvxqp3").
        category: Target subdirectory ("easy-indefinite", "hard-indefinite", "positive-definite").
        dest_base: Base directory for downloads.
        max_rows: Maximum number of rows to download (skip if larger).
        dry_run: If True, only report what would be downloaded.

    Returns:
        Tuple of (metadata_dict, status_string) or (None, error_string).
    """
    import ssgetpy

    group, name = group_name.split("/")

    try:
        results = ssgetpy.search(name, limit=10)
    except Exception as e:
        return None, f"search failed: {e}"

    if not results:
        return None, "not found in SuiteSparse"

    # Find exact match by group and name
    mat = None
    for r in results:
        if r.group == group and r.name == name:
            mat = r
            break
    if mat is None:
        return None, f"no exact match for {group_name} in search results"

    rows = mat.rows
    nnz = mat.nnz

    if max_rows and rows > max_rows:
        return None, f"skipped: {rows} rows > max_rows={max_rows}"

    # Check if matrix should be download-on-demand (too large for repo)
    in_repo = rows <= LARGE_THRESHOLD_ROWS

    if dry_run:
        status = "would download" if in_repo else "metadata-only (too large)"
        meta = {
            "name": group_name,
            "rows": rows,
            "nnz": nnz,
            "in_repo": in_repo,
            "status": status,
        }
        return meta, status

    # Build destination path
    dest_dir = os.path.join(dest_base, "suitesparse", category, name)
    os.makedirs(dest_dir, exist_ok=True)
    mtx_path = os.path.join(dest_dir, f"{name}.mtx")

    if in_repo:
        # Skip if already downloaded
        if os.path.exists(mtx_path) and os.path.getsize(mtx_path) > 100:
            pass  # Already have it
        else:
            # Download the matrix
            try:
                mat.download(format="MM", destpath=dest_dir, extract=False)
            except Exception as e:
                return None, f"download failed: {e}"

            # Find and extract the tar.gz
            targz_candidates = glob.glob(os.path.join(dest_dir, "*.tar.gz"))
            if targz_candidates:
                targz_path = targz_candidates[0]
                try:
                    with tarfile.open(targz_path, "r:gz") as tar:
                        tar.extractall(path=dest_dir)
                except Exception as e:
                    return None, f"extraction failed: {e}"
                os.remove(targz_path)

            # Find the correct .mtx file (may be in a subdirectory)
            # SuiteSparse tarballs contain multiple files; we want the main one (name/name.mtx)
            mtx_candidates = glob.glob(os.path.join(dest_dir, "**", "*.mtx"), recursive=True)
            if not mtx_candidates:
                return None, "download succeeded but no .mtx file found"

            # Prefer the file named exactly {name}.mtx
            actual_mtx = None
            for candidate in mtx_candidates:
                if os.path.basename(candidate) == f"{name}.mtx":
                    actual_mtx = candidate
                    break
            if actual_mtx is None:
                actual_mtx = mtx_candidates[0]
            # Move to expected location if needed
            if actual_mtx != mtx_path:
                shutil.move(actual_mtx, mtx_path)

            # Clean up any extra extracted directories and tar.gz files
            for item in os.listdir(dest_dir):
                item_path = os.path.join(dest_dir, item)
                if os.path.isdir(item_path):
                    shutil.rmtree(item_path)
                elif item_path.endswith(".tar.gz"):
                    os.remove(item_path)

    # Build metadata
    is_pd = category == "positive-definite"

    properties = {
        "symmetric": True,
        "positive_definite": is_pd,
        "indefinite": not is_pd,
        "difficulty": "hard" if category == "hard-indefinite" else "easy",
    }

    relative_path = f"suitesparse/{category}/{name}/{name}.mtx"

    meta = create_matrix_entry(
        name=group_name,
        source="suitesparse",
        category=category,
        path=relative_path,
        size=rows,
        nnz=nnz,
        properties=properties,
        in_repo=in_repo,
    )

    if not in_repo:
        meta["download_command"] = (
            f"python download_suitesparse.py --name \"{group_name}\""
        )

    return meta, "downloaded" if in_repo else "metadata-only"


# ============================================================================
# CLI interface (T018)
# ============================================================================

def get_curated_list(category=None):
    """Return the list of (group_name, properties) for the given category."""
    lists = {
        "easy-indefinite": EASY_INDEFINITE,
        "hard-indefinite": HARD_INDEFINITE,
        "positive-definite": POSITIVE_DEFINITE,
    }
    if category:
        return [(category, lists[category])]
    return list(lists.items())


def main():
    parser = argparse.ArgumentParser(
        description="Download curated SuiteSparse matrices for SSIDS testing"
    )
    parser.add_argument(
        "--category",
        choices=["easy-indefinite", "hard-indefinite", "positive-definite"],
        help="Download only this category",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        help="Skip matrices with more than N rows",
    )
    parser.add_argument(
        "--name",
        help='Download a specific matrix by "Group/Name"',
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Check existing files without downloading",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List what would be downloaded without downloading",
    )
    args = parser.parse_args()

    metadata = load_metadata(METADATA_PATH)
    new_entries = []

    categories = get_curated_list(args.category)

    total = 0
    downloaded = 0
    skipped = 0
    meta_only = 0
    errors = 0

    for cat_name, matrices in categories:
        print(f"\n{'=' * 60}")
        print(f"Category: {cat_name} ({len(matrices)} matrices)")
        print(f"{'=' * 60}")

        for group_name, props in matrices:
            if args.name and args.name != group_name:
                continue

            total += 1

            if args.verify_only:
                # Check if the matrix exists in metadata
                existing = [m for m in metadata["matrices"] if m["name"] == group_name]
                if existing:
                    entry = existing[0]
                    if entry.get("in_repo", True):
                        mtx_path = os.path.join(BASE_DIR, entry["path"])
                        exists = os.path.exists(mtx_path)
                        status = "OK" if exists else "MISSING"
                    else:
                        status = "metadata-only (OK)"
                    print(f"  {group_name}: {status}")
                else:
                    print(f"  {group_name}: NOT IN METADATA")
                continue

            # Download or dry-run
            meta, status = download_matrix(
                group_name, cat_name, BASE_DIR,
                max_rows=args.max_rows,
                dry_run=args.dry_run,
            )

            if meta is None:
                print(f"  {group_name}: {status}")
                if "skipped" in status:
                    skipped += 1
                else:
                    errors += 1
                continue

            print(f"  {group_name}: {status} (rows={meta.get('rows', meta.get('size', '?'))})")

            if not args.dry_run:
                # Inject curated properties into the metadata entry
                if "kind" in props:
                    meta["properties"]["kind"] = props["kind"]
                if "expected_delayed_pivots" in props:
                    meta["properties"]["expected_delayed_pivots"] = props["expected_delayed_pivots"]
                if props.get("killer_case"):
                    meta["properties"]["killer_case"] = True
                if "paper_references" in props:
                    meta["paper_references"] = props["paper_references"]

                new_entries.append(meta)

                if meta.get("in_repo", True):
                    downloaded += 1
                else:
                    meta_only += 1

    if not args.verify_only and not args.dry_run and new_entries:
        metadata = merge_metadata(metadata, new_entries)
        write_metadata(metadata, METADATA_PATH)
        print(f"\nUpdated {METADATA_PATH}")

    print(f"\nSummary:")
    print(f"  Total matrices: {total}")
    if not args.verify_only:
        print(f"  Downloaded: {downloaded}")
        print(f"  Metadata-only: {meta_only}")
        print(f"  Skipped (size): {skipped}")
        print(f"  Errors: {errors}")


if __name__ == "__main__":
    main()
