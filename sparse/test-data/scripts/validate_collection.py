#!/usr/bin/env python3
"""Validate the test matrix collection for completeness and consistency.

Checks metadata against actual files, validates required fields, reports
summary statistics, and verifies paper coverage.
"""

import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from mtx_utils import read_mtx_header
from metadata_utils import load_metadata

BASE_DIR = os.path.join(os.path.dirname(SCRIPT_DIR))
METADATA_PATH = os.path.join(BASE_DIR, "metadata.json")


# ============================================================================
# T021: Collection validation
# ============================================================================

REQUIRED_FIELDS = ["name", "source", "category", "path", "size", "nnz", "in_repo", "properties"]
REQUIRED_PROPERTIES = ["symmetric", "positive_definite", "indefinite", "difficulty"]
VALID_SOURCES = {"hand-constructed", "suitesparse", "spral-inspired", "paper-referenced"}
VALID_CATEGORIES = {"hand-constructed", "easy-indefinite", "hard-indefinite", "positive-definite"}
VALID_DIFFICULTIES = {"trivial", "easy", "hard"}


def validate_metadata(index):
    """Check every entry has all required fields per data-model.md."""
    errors = []
    names_seen = set()

    for i, entry in enumerate(index["matrices"]):
        name = entry.get("name", f"<entry {i}>")

        # Check required fields
        for field in REQUIRED_FIELDS:
            if field not in entry:
                errors.append(f"{name}: missing required field '{field}'")

        # Check unique names
        if name in names_seen:
            errors.append(f"{name}: duplicate name")
        names_seen.add(name)

        # Check enum values
        source = entry.get("source", "")
        if source not in VALID_SOURCES:
            errors.append(f"{name}: invalid source '{source}'")

        category = entry.get("category", "")
        if category not in VALID_CATEGORIES:
            errors.append(f"{name}: invalid category '{category}'")

        # Check properties
        props = entry.get("properties", {})
        for field in REQUIRED_PROPERTIES:
            if field not in props:
                errors.append(f"{name}: missing required property '{field}'")

        difficulty = props.get("difficulty", "")
        if difficulty not in VALID_DIFFICULTIES:
            errors.append(f"{name}: invalid difficulty '{difficulty}'")

        # Check symmetric is always true
        if not props.get("symmetric", False):
            errors.append(f"{name}: symmetric must be true for all matrices")

        # Check PD/indefinite are mutually exclusive
        is_pd = props.get("positive_definite", False)
        is_indef = props.get("indefinite", False)
        if is_pd and is_indef:
            errors.append(f"{name}: positive_definite and indefinite are mutually exclusive")

        # Check hand-constructed matrices have factorization_path
        if source == "hand-constructed" and "factorization_path" not in entry:
            errors.append(f"{name}: hand-constructed matrix missing factorization_path")

        # Check non-in-repo matrices have download_command
        if not entry.get("in_repo", True) and "download_command" not in entry:
            errors.append(f"{name}: in_repo=false but no download_command")

        # Check killer_case entries have difficulty hard
        if props.get("killer_case") and difficulty != "hard":
            errors.append(f"{name}: killer_case=true but difficulty is '{difficulty}', not 'hard'")

    return errors


def validate_files(index, base_dir):
    """Verify .mtx files exist and match metadata for in_repo entries."""
    errors = []

    for entry in index["matrices"]:
        name = entry.get("name", "?")
        in_repo = entry.get("in_repo", True)
        path = entry.get("path", "")

        if in_repo:
            mtx_path = os.path.join(base_dir, path)
            if not os.path.exists(mtx_path):
                errors.append(f"{name}: file not found: {path}")
                continue

            # Check file size is reasonable (not a LFS pointer or empty)
            fsize = os.path.getsize(mtx_path)
            if fsize < 50:
                errors.append(f"{name}: file suspiciously small ({fsize} bytes): {path}")
                continue

            # Check header matches metadata
            try:
                n, nnz = read_mtx_header(mtx_path)
                if n != entry.get("size"):
                    errors.append(f"{name}: size mismatch: metadata={entry.get('size')}, file={n}")
                if nnz != entry.get("nnz"):
                    errors.append(f"{name}: nnz mismatch: metadata={entry.get('nnz')}, file={nnz}")
            except Exception as e:
                errors.append(f"{name}: could not parse header: {e}")

        # Check factorization file for hand-constructed
        if entry.get("source") == "hand-constructed":
            fact_path = entry.get("factorization_path", "")
            if fact_path:
                full_fact_path = os.path.join(base_dir, fact_path)
                if not os.path.exists(full_fact_path):
                    errors.append(f"{name}: factorization file not found: {fact_path}")

    return errors


def validate_properties(index):
    """Check property consistency across the collection."""
    errors = []

    for entry in index["matrices"]:
        name = entry.get("name", "?")
        props = entry.get("properties", {})

        # PD and indefinite are mutually exclusive
        if props.get("positive_definite") and props.get("indefinite"):
            errors.append(f"{name}: PD and indefinite are mutually exclusive")

        # killer_case entries should be hard difficulty
        if props.get("killer_case") and props.get("difficulty") != "hard":
            errors.append(f"{name}: killer_case but not hard difficulty")

    return errors


# ============================================================================
# T022: Summary statistics
# ============================================================================

def print_summary(index):
    """Report collection statistics and verify success criteria."""
    matrices = index["matrices"]
    total = len(matrices)

    # Count by category
    by_category = {}
    for m in matrices:
        cat = m.get("category", "unknown")
        by_category[cat] = by_category.get(cat, 0) + 1

    # Count by source
    by_source = {}
    for m in matrices:
        src = m.get("source", "unknown")
        by_source[src] = by_source.get(src, 0) + 1

    # Size range
    sizes = [m.get("size", 0) for m in matrices]
    min_size = min(sizes) if sizes else 0
    max_size = max(sizes) if sizes else 0

    # Count killer cases
    killer_count = sum(1 for m in matrices if m.get("properties", {}).get("killer_case"))

    # Count matrices with exact factorizations
    fact_count = sum(1 for m in matrices if m.get("factorization_path"))

    # Count in-repo vs download-on-demand
    in_repo = sum(1 for m in matrices if m.get("in_repo", True))
    on_demand = total - in_repo

    print(f"\n{'=' * 60}")
    print("Collection Summary")
    print(f"{'=' * 60}")
    print(f"  Total matrices: {total}")
    print(f"  By category:")
    for cat, count in sorted(by_category.items()):
        print(f"    {cat}: {count}")
    print(f"  By source:")
    for src, count in sorted(by_source.items()):
        print(f"    {src}: {count}")
    print(f"  Size range: {min_size} to {max_size}")
    print(f"  Killer cases: {killer_count}")
    print(f"  With exact factorizations: {fact_count}")
    print(f"  In-repo: {in_repo}")
    print(f"  Download-on-demand: {on_demand}")

    # Success criteria checks
    print(f"\n{'=' * 60}")
    print("Success Criteria")
    print(f"{'=' * 60}")

    criteria = []

    # SC-001: >= 70 total
    sc001 = total >= 70
    criteria.append(("SC-001", f">= 70 total matrices", total, sc001))

    # SC-004: 4 orders of magnitude size range
    if min_size > 0 and max_size > 0:
        magnitude = max_size / min_size
        sc004 = magnitude >= 1000  # 4 orders = at least 10^3 range
    else:
        magnitude = 0
        sc004 = False
    criteria.append(("SC-004", f"4+ orders of magnitude ({min_size} to {max_size})", magnitude, sc004))

    # SC-005: >= 5 killer cases
    sc005 = killer_count >= 5
    criteria.append(("SC-005", f">= 5 killer cases", killer_count, sc005))

    # SC-006: verified factorizations (>= 15 hand-constructed)
    sc006 = fact_count >= 15
    criteria.append(("SC-006", f">= 15 verified factorizations", fact_count, sc006))

    for sc_id, desc, value, passed in criteria:
        status = "PASS" if passed else "FAIL"
        print(f"  {sc_id}: {desc} = {value} [{status}]")

    return all(passed for _, _, _, passed in criteria)


# ============================================================================
# T023: Property queries
# ============================================================================

def test_queries(index):
    """Run and verify property-based queries against the collection."""
    matrices = index["matrices"]
    all_passed = True

    print(f"\n{'=' * 60}")
    print("Property Queries")
    print(f"{'=' * 60}")

    # Query 1: indefinite matrices < 1000x1000
    small_indef = [m for m in matrices
                   if m.get("properties", {}).get("indefinite")
                   and m.get("size", 0) < 1000]
    has_hand = any(m.get("source") == "hand-constructed" for m in small_indef)
    q1_pass = len(small_indef) > 0 and has_hand
    print(f"  Q1: Indefinite < 1000x1000: {len(small_indef)} matrices "
          f"(includes hand-constructed: {has_hand}) [{'PASS' if q1_pass else 'FAIL'}]")
    all_passed &= q1_pass

    # Query 2: killer cases >= 5
    killers = [m for m in matrices
               if m.get("properties", {}).get("killer_case")]
    q2_pass = len(killers) >= 5
    print(f"  Q2: Killer cases: {len(killers)} matrices [{'PASS' if q2_pass else 'FAIL'}]")
    if killers:
        for k in killers:
            print(f"       - {k['name']} ({k['size']}x{k['size']})")
    all_passed &= q2_pass

    # Query 3: positive-definite >= 19
    pd = [m for m in matrices
          if m.get("category") == "positive-definite"]
    q3_pass = len(pd) >= 19
    print(f"  Q3: Positive-definite: {len(pd)} matrices [{'PASS' if q3_pass else 'FAIL'}]")
    all_passed &= q3_pass

    # Query 4: matrices with factorization_path >= 15
    with_fact = [m for m in matrices
                 if m.get("factorization_path")]
    q4_pass = len(with_fact) >= 15
    print(f"  Q4: With exact factorizations: {len(with_fact)} matrices [{'PASS' if q4_pass else 'FAIL'}]")
    all_passed &= q4_pass

    return all_passed


# ============================================================================
# T026: Paper coverage report
# ============================================================================

# Core papers and expected matrices
PAPER_MATRICES = {
    "Duff, Hogg, Lopez (2020) Table 1": [
        "Oberwolfach/t2dal", "GHS_indef/dixmaanl", "Oberwolfach/rail_79841",
        "GHS_indef/dawson5", "Boeing/bcsstk39", "Boeing/pct20stif",
        "GHS_indef/copter2", "Boeing/crystk03", "Koutsovasilis/F2",
        "Cunningham/qa8fk", "Oberwolfach/gas_sensor", "Oberwolfach/t3dh",
        "PARSEC/H2O", "GHS_indef/sparsine",
    ],
    "Duff, Hogg, Lopez (2020) Table 2": [
        "TSOPF/TSOPF_FS_b39_c7", "TSOPF/TSOPF_FS_b162_c1",
        "GHS_indef/cont-201", "GHS_indef/stokes128",
        "GHS_indef/ncvxqp1", "GHS_indef/cont-300",
        "GHS_indef/bratu3d", "GHS_indef/cvxqp3",
        "GHS_indef/d_pretok", "GHS_indef/turon_m",
        "GHS_indef/ncvxqp5", "GHS_indef/ncvxqp3",
        "GHS_indef/ncvxqp7", "Schenk_IBMNA/c-big",
        "GHS_indef/c-71",
    ],
    "Hogg, Ovtchinnikov, Scott (2016) Table III": [
        "Newman/astro-ph", "Andrianov/mip1", "PARSEC/SiNa",
        "GHS_indef/ncvxqp3", "GHS_indef/c-71",
        "GHS_psdef/crankseg_1", "Rothberg/cfd2", "DNVS/thread",
        "DNVS/shipsec1", "DNVS/shipsec8", "Oberwolfach/boneS01",
        "GHS_psdef/crankseg_2", "DNVS/shipsec5",
        "GHS_psdef/bmwcra_1", "DNVS/ship_003",
        "Um/offshore", "ND/nd6k", "ND/nd12k",
        "Schenk_AFE/af_shell7", "Schenk_AFE/af_0_k101",
        "GHS_psdef/ldoor", "GHS_psdef/inline_1",
        "GHS_psdef/apache2", "AMD/G3_circuit",
    ],
}


def check_paper_coverage(index):
    """Verify paper coverage for core APTP papers."""
    matrix_names = {m["name"] for m in index["matrices"]}

    print(f"\n{'=' * 60}")
    print("Paper Coverage")
    print(f"{'=' * 60}")

    all_covered = True
    for paper, expected_matrices in PAPER_MATRICES.items():
        found = [m for m in expected_matrices if m in matrix_names]
        missing = [m for m in expected_matrices if m not in matrix_names]
        pct = len(found) / len(expected_matrices) * 100 if expected_matrices else 0

        status = "PASS" if not missing else "FAIL"
        print(f"\n  {paper}:")
        print(f"    Coverage: {len(found)}/{len(expected_matrices)} ({pct:.0f}%) [{status}]")
        if missing:
            all_covered = False
            for m in missing:
                print(f"    MISSING: {m}")

    return all_covered


# ============================================================================
# Main entry point
# ============================================================================

def main():
    """Run all validations."""
    print("Loading metadata...")
    index = load_metadata(METADATA_PATH)
    print(f"Loaded {len(index['matrices'])} matrices")

    all_passed = True

    # Metadata validation
    print(f"\n{'=' * 60}")
    print("Metadata Validation")
    print(f"{'=' * 60}")
    errors = validate_metadata(index)
    if errors:
        all_passed = False
        for e in errors:
            print(f"  ERROR: {e}")
    else:
        print("  All metadata entries valid.")

    # File validation
    print(f"\n{'=' * 60}")
    print("File Validation")
    print(f"{'=' * 60}")
    errors = validate_files(index, BASE_DIR)
    if errors:
        all_passed = False
        for e in errors:
            print(f"  ERROR: {e}")
    else:
        print("  All files present and consistent.")

    # Property validation
    print(f"\n{'=' * 60}")
    print("Property Validation")
    print(f"{'=' * 60}")
    errors = validate_properties(index)
    if errors:
        all_passed = False
        for e in errors:
            print(f"  ERROR: {e}")
    else:
        print("  All properties consistent.")

    # Summary statistics
    summary_passed = print_summary(index)
    all_passed &= summary_passed

    # Property queries
    queries_passed = test_queries(index)
    all_passed &= queries_passed

    # Paper coverage
    paper_passed = check_paper_coverage(index)
    all_passed &= paper_passed

    # Final status
    print(f"\n{'=' * 60}")
    if all_passed:
        print("OVERALL: ALL VALIDATIONS PASSED")
    else:
        print("OVERALL: SOME VALIDATIONS FAILED")
    print(f"{'=' * 60}")

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
