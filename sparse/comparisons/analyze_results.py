#!/usr/bin/env python3
"""Analyze SPRAL vs rivrs benchmark results from JSON benchmark files.

Reads the structured JSON output produced by `cargo run --bin spral-comparison`.
Computes summary statistics, breakdowns by matrix category and application
domain, and parallel scaling when multiple thread counts are provided.

Usage:
    # Single file (latest T8 run):
    python comparisons/analyze_results.py target/benchmarks/spral/spral-benchmark-1772124268.json

    # Multiple files (T1 + T8 for parallel scaling comparison):
    python comparisons/analyze_results.py \\
        target/benchmarks/spral/spral-benchmark-1772118484.json \\
        target/benchmarks/spral/spral-benchmark-1772124268.json

    # Glob all runs:
    python comparisons/analyze_results.py target/benchmarks/spral/*.json

    # Latest file by default (newest in target/benchmarks/spral/):
    python comparisons/analyze_results.py
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Application domain mapping (SuiteSparse Matrix Collection "Kind" metadata)
# ---------------------------------------------------------------------------
DOMAIN_MAP = {
    # Structural FEM (SuiteSparse "Kind": structural problem)
    "BenElechi/BenElechi1": "structural",
    "Koutsovasilis/F2": "structural",
    "Boeing/bcsstk39": "structural",
    "Boeing/crystk02": "structural",
    "Boeing/crystk03": "structural",
    "Boeing/pwtk": "structural",
    "INPRO/msdoor": "structural",
    "GHS_indef/dawson5": "structural",
    "GHS_indef/sparsine": "structural",
    "Schenk_AFE/af_0_k101": "structural",
    "Schenk_AFE/af_shell7": "structural",
    "GHS_psdef/apache2": "structural",
    "GHS_psdef/bmwcra_1": "structural",
    "GHS_psdef/crankseg_1": "structural",
    "GHS_psdef/crankseg_2": "structural",
    "GHS_psdef/inline_1": "structural",
    "GHS_psdef/ldoor": "structural",
    "DNVS/ship_003": "structural",
    "DNVS/shipsec1": "structural",
    "DNVS/shipsec5": "structural",
    "DNVS/shipsec8": "structural",
    "DNVS/thread": "structural",
    # PDE discretization (SuiteSparse "Kind": 2D/3D problem)
    "GHS_indef/bratu3d": "pde",
    "GHS_indef/helm3d01": "pde",
    "GHS_indef/aug3dcqp": "pde",
    "GHS_indef/d_pretok": "pde",
    "GHS_indef/mario001": "pde",
    "GHS_indef/turon_m": "pde",
    "GHS_indef/bloweybq": "pde",
    "ND/nd3k": "pde",
    "ND/nd6k": "pde",
    "ND/nd12k": "pde",
    "Um/offshore": "pde",
    # Optimization (QP, KKT, CUTEr test set)
    "GHS_indef/c-71": "optimization",
    "Schenk_IBMNA/c-big": "optimization",
    "GHS_indef/cont-201": "optimization",
    "GHS_indef/cont-300": "optimization",
    "GHS_indef/dixmaanl": "optimization",
    "GHS_indef/blockqp1": "optimization",
    "GHS_indef/cvxqp3": "optimization",
    "GHS_indef/ncvxqp1": "optimization",
    "GHS_indef/ncvxqp3": "optimization",
    "GHS_indef/ncvxqp5": "optimization",
    "GHS_indef/ncvxqp7": "optimization",
    "GHS_indef/linverse": "optimization",
    "GHS_indef/spmsrtls": "optimization",
    # Quantum chemistry (PARSEC collection: DFT)
    "PARSEC/H2O": "quantum chemistry",
    "PARSEC/Si10H16": "quantum chemistry",
    "PARSEC/Si5H12": "quantum chemistry",
    "PARSEC/SiNa": "quantum chemistry",
    # Model reduction (Oberwolfach benchmark collection)
    "Oberwolfach/filter3D": "model reduction",
    "Oberwolfach/gas_sensor": "model reduction",
    "Oberwolfach/rail_79841": "model reduction",
    "Oberwolfach/t2dal": "model reduction",
    "Oberwolfach/t3dh": "model reduction",
    "Oberwolfach/boneS01": "model reduction",
    # CFD (Stokes, rotor, fluid dynamics)
    "GHS_indef/copter2": "cfd",
    "GHS_indef/stokes128": "cfd",
    "Rothberg/cfd2": "cfd",
    # Acoustics / vibration
    "Cote/vibrobox": "acoustics",
    "Cunningham/qa8fk": "acoustics",
    # Circuit simulation
    "AMD/G3_circuit": "circuit",
    # Power systems (transient optimal power flow)
    "TSOPF/TSOPF_FS_b162_c1": "power systems",
    "TSOPF/TSOPF_FS_b39_c7": "power systems",
    # Graph / network
    "Newman/astro-ph": "graph/network",
}

DOMAIN_DESCRIPTIONS = {
    "structural": "Structural mechanics FEM (stiffness/mass from engineering structures)",
    "pde": "PDE discretization (2D/3D FE/FD, Helmholtz, geophysics, electromagnetics)",
    "optimization": "Optimization (convex/nonconvex QP, KKT systems, CUTEr test set)",
    "quantum chemistry": "Quantum chemistry (real-space pseudopotential DFT, PARSEC)",
    "model reduction": "Model reduction benchmarks (Oberwolfach: thermal, sensor, bio-mechanics)",
    "cfd": "Computational fluid dynamics (Stokes, Navier-Stokes, rotor grids)",
    "acoustics": "Acoustics and vibration (vibroacoustic FEM, wave propagation)",
    "circuit": "Circuit simulation (VLSI, IC interconnect networks)",
    "power systems": "Power network systems (transient optimal power flow)",
    "graph/network": "Graph/network problems (collaboration networks, social graphs)",
}


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load_json(path: Path) -> pd.DataFrame:
    """Load a spral-benchmark JSON file into a comparison DataFrame.

    Returns a DataFrame with columns:
        matrix, category, n, nnz, spral_fac, rivrs_fac, ratio,
        spral_be, rivrs_be, sl_st, sl_nd, threads
    """
    with open(path) as f:
        data = json.load(f)

    threads = data.get("omp_num_threads") or 1
    rows = []
    for rec in data["results"]:
        spral = rec.get("spral")
        rivrs = rec.get("rivrs")
        if spral is None or rivrs is None:
            continue

        spral_fac = spral["factor_s"]
        rivrs_fac = rivrs["factor_s"]
        ratio = rivrs_fac / spral_fac if spral_fac > 0 else float("inf")

        rows.append({
            "matrix": rec["matrix_name"],
            "category": rec.get("category", "unknown"),
            "n": rec["n"],
            "nnz": rec.get("nnz", 0),
            "spral_fac": spral_fac,
            "rivrs_fac": rivrs_fac,
            "ratio": ratio,
            "spral_be": spral.get("backward_error", float("nan")),
            "rivrs_be": rivrs.get("backward_error", float("nan")),
            "sl_st": rivrs.get("small_leaf_subtrees", 0),
            "sl_nd": rivrs.get("small_leaf_nodes", 0),
            "threads": threads,
            # Extra fields available from JSON
            "spral_analyse": spral.get("analyse_s", float("nan")),
            "rivrs_analyse": rivrs.get("analyse_s", float("nan")),
            "spral_solve": spral.get("solve_s", float("nan")),
            "rivrs_solve": rivrs.get("solve_s", float("nan")),
            "spral_total": spral.get("total_s", float("nan")),
            "rivrs_total": rivrs.get("total_s", float("nan")),
        })

    df = pd.DataFrame(rows)
    df.attrs["path"] = str(path)
    df.attrs["timestamp"] = data.get("timestamp", "")
    df.attrs["rivrs_version"] = data.get("rivrs_version", "?")
    df.attrs["spral_version"] = data.get("spral_version", "?")
    df.attrs["platform"] = data.get("platform", "?")
    return df


def find_default_dir() -> Path:
    """Find the benchmark output directory."""
    candidates = [
        Path("target/benchmarks/spral"),
        Path("sparse/target/benchmarks/spral"),
    ]
    for c in candidates:
        if c.is_dir():
            return c
    return candidates[0]


def load_inputs(paths: list[Path]) -> dict[str, pd.DataFrame]:
    """Load one or more JSON files, keyed by 'T{threads}'.

    When multiple files have the same thread count, the newest
    (by timestamp in the JSON) is kept.
    """
    results: dict[str, pd.DataFrame] = {}
    timestamps: dict[str, str] = {}

    for p in paths:
        df = load_json(p)
        if df.empty:
            print(f"  warning: {p.name} has no comparison results, skipping", file=sys.stderr)
            continue

        threads = df["threads"].iloc[0]
        label = f"T{threads}"
        ts = df.attrs.get("timestamp", "")

        if label not in results or ts > timestamps.get(label, ""):
            results[label] = df
            timestamps[label] = ts

    return results


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
def summary_stats(df: pd.DataFrame, label: str) -> str:
    """Compute and format summary statistics for a comparison DataFrame."""
    ratios = df["ratio"]
    n = len(ratios)

    total_spral = df["spral_fac"].sum()
    total_rivrs = df["rivrs_fac"].sum()
    time_weighted = total_rivrs / total_spral

    faster = (ratios < 0.90).sum()
    comparable = ((ratios >= 0.90) & (ratios <= 1.10)).sum()
    slower = (ratios > 1.10).sum()

    lines = [
        f"**{label} (threads={df['threads'].iloc[0]}) — {n} matrices:**",
        "",
        "| Statistic | Value |",
        "|-----------|-------|",
        f"| Median ratio | {ratios.median():.2f}× |",
        f"| Mean ratio | {ratios.mean():.2f}× |",
        f"| 25th percentile | {ratios.quantile(0.25):.2f}× |",
        f"| 75th percentile | {ratios.quantile(0.75):.2f}× |",
        f"| Time-weighted ratio | {time_weighted:.2f}× |",
        f"| Rivrs faster (<0.90×) | {faster}/{n} ({100*faster//n}%) |",
        f"| Comparable (0.90–1.10×) | {comparable}/{n} ({100*comparable//n}%) |",
        f"| SPRAL faster (>1.10×) | {slower}/{n} ({100*slower//n}%) |",
    ]

    # Top 5 wins for each
    top_rivrs = df.nsmallest(5, "ratio")[["matrix", "ratio"]]
    top_spral = df.nlargest(5, "ratio")[["matrix", "ratio"]]

    rivrs_str = ", ".join(
        f"{row.matrix.split('/')[-1]} {row.ratio:.2f}×" for row in top_rivrs.itertuples()
    )
    spral_str = ", ".join(
        f"{row.matrix.split('/')[-1]} {row.ratio:.2f}×" for row in top_spral.itertuples()
    )
    lines += [
        "",
        f"Top rivrs wins: {rivrs_str}.",
        f"Top SPRAL wins: {spral_str}.",
    ]

    return "\n".join(lines)


def category_table(dfs: dict[str, pd.DataFrame]) -> str:
    """Tabulate speedup by matrix category (easy/hard indefinite, PD)."""
    rows = []
    for label, df in sorted(dfs.items()):
        for cat in ["easy-indefinite", "hard-indefinite", "positive-definite"]:
            sub = df[df["category"] == cat]
            if sub.empty:
                continue
            tw = sub["rivrs_fac"].sum() / sub["spral_fac"].sum()
            rows.append({
                "threads": label,
                "category": cat,
                "count": len(sub),
                "median": sub["ratio"].median(),
                "mean": sub["ratio"].mean(),
                "time_weighted": tw,
                "min": sub["ratio"].min(),
                "max": sub["ratio"].max(),
            })

    tbl = pd.DataFrame(rows)
    lines = ["**By category:**", ""]
    for threads_label in sorted(tbl["threads"].unique()):
        sub = tbl[tbl["threads"] == threads_label]
        lines.append(f"*{threads_label}:*")
        lines.append("")
        lines.append("| Category | Count | Median | Mean | Time-wtd | Min | Max |")
        lines.append("|----------|------:|-------:|-----:|---------:|----:|----:|")
        for _, row in sub.iterrows():
            lines.append(
                f"| {row['category']} | {row['count']} "
                f"| {row['median']:.2f}× | {row['mean']:.2f}× "
                f"| {row['time_weighted']:.2f}× "
                f"| {row['min']:.2f}× | {row['max']:.2f}× |"
            )
        lines.append("")
    return "\n".join(lines)


def domain_table(dfs: dict[str, pd.DataFrame]) -> str:
    """Tabulate speedup by application domain."""
    rows = []
    for label, df in sorted(dfs.items()):
        df = df.copy()
        df["domain"] = df["matrix"].map(lambda m: DOMAIN_MAP.get(m, "unknown"))
        for domain in sorted(df["domain"].unique()):
            sub = df[df["domain"] == domain]
            if sub.empty:
                continue
            tw = sub["rivrs_fac"].sum() / sub["spral_fac"].sum()
            rows.append({
                "threads": label,
                "domain": domain,
                "count": len(sub),
                "median": sub["ratio"].median(),
                "mean": sub["ratio"].mean(),
                "time_weighted": tw,
            })

    tbl = pd.DataFrame(rows)
    lines = ["**By application domain:**", ""]
    for threads_label in sorted(tbl["threads"].unique()):
        sub = tbl[tbl["threads"] == threads_label].sort_values("median")
        lines.append(f"*{threads_label}:*")
        lines.append("")
        lines.append("| Domain | Count | Median | Mean | Time-wtd |")
        lines.append("|--------|------:|-------:|-----:|---------:|")
        for _, row in sub.iterrows():
            lines.append(
                f"| {row['domain']} | {row['count']} "
                f"| {row['median']:.2f}× | {row['mean']:.2f}× "
                f"| {row['time_weighted']:.2f}× |"
            )
        lines.append("")
    return "\n".join(lines)


def parallel_scaling(dfs: dict[str, pd.DataFrame]) -> str:
    """Compute parallel scaling stats across the lowest and highest thread counts."""
    labels = sorted(dfs.keys(), key=lambda k: dfs[k]["threads"].iloc[0])
    if len(labels) < 2:
        return ""

    lo_label, hi_label = labels[0], labels[-1]
    lo = dfs[lo_label].set_index("matrix")
    hi = dfs[hi_label].set_index("matrix")
    lo_t = lo["threads"].iloc[0]
    hi_t = hi["threads"].iloc[0]
    common = lo.index.intersection(hi.index)

    speedup = lo.loc[common, "rivrs_fac"] / hi.loc[common, "rivrs_fac"]

    top = speedup.nlargest(5)
    top_str = ", ".join(f"{m.split('/')[-1]} {v:.2f}×" for m, v in top.items())

    return (
        f"**Parallel scaling (rivrs T{lo_t}→T{hi_t}):** "
        f"Mean {speedup.mean():.2f}×, median {speedup.median():.2f}×. "
        f"Best: {top_str}."
    )


def phase_breakdown(dfs: dict[str, pd.DataFrame]) -> str:
    """Show analyse/factor/solve time breakdown (available from JSON)."""
    lines = ["**Phase breakdown (median, seconds):**", ""]
    for label in sorted(dfs.keys()):
        df = dfs[label]
        # Only show if we have the phase columns
        if "spral_analyse" not in df.columns:
            continue
        lines.append(f"*{label}:*")
        lines.append("")
        lines.append("| Phase | SPRAL median | rivrs median | Ratio |")
        lines.append("|-------|------------:|-----------:|------:|")
        for phase_s, phase_r, name in [
            ("spral_analyse", "rivrs_analyse", "Analyse"),
            ("spral_fac", "rivrs_fac", "Factor"),
            ("spral_solve", "rivrs_solve", "Solve"),
            ("spral_total", "rivrs_total", "Total"),
        ]:
            s = df[phase_s].median()
            r = df[phase_r].median()
            ratio = r / s if s > 0 else float("inf")
            lines.append(f"| {name} | {s:.3f} | {r:.3f} | {ratio:.2f}× |")
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Analyze SPRAL vs rivrs benchmark JSON results.",
        epilog="If no files are given, loads the newest file from target/benchmarks/spral/.",
    )
    parser.add_argument(
        "files",
        nargs="*",
        type=Path,
        help="One or more spral-benchmark-*.json files",
    )
    args = parser.parse_args()

    # Resolve input files
    if args.files:
        json_paths = [p for p in args.files if p.suffix == ".json" and p.exists()]
        if not json_paths:
            print(f"Error: no valid .json files found in {args.files}", file=sys.stderr)
            sys.exit(1)
    else:
        bench_dir = find_default_dir()
        if not bench_dir.is_dir():
            print(f"Error: {bench_dir} not found. Pass JSON file(s) explicitly.", file=sys.stderr)
            sys.exit(1)
        json_paths = sorted(bench_dir.glob("spral-benchmark-*.json"))
        if not json_paths:
            print(f"Error: no benchmark files in {bench_dir}", file=sys.stderr)
            sys.exit(1)
        # Take the newest file only
        json_paths = [json_paths[-1]]
        print(f"(using {json_paths[0].name})", file=sys.stderr)

    dfs = load_inputs(json_paths)

    if not dfs:
        print("No comparison results found in input files", file=sys.stderr)
        sys.exit(1)

    # Print metadata from first DataFrame
    first_df = next(iter(dfs.values()))
    print("=" * 72)
    print("SPRAL vs rivrs — Benchmark Analysis")
    print("=" * 72)
    platform = first_df.attrs.get("platform", "")
    rivrs_v = first_df.attrs.get("rivrs_version", "")
    spral_v = first_df.attrs.get("spral_version", "")
    if platform:
        print(f"Platform: {platform}")
    if rivrs_v or spral_v:
        print(f"rivrs: {rivrs_v}  SPRAL: {spral_v}")
    labels = sorted(dfs.keys(), key=lambda k: dfs[k]["threads"].iloc[0])
    threads_str = ", ".join(
        f"{l} ({dfs[l]['threads'].iloc[0]} threads, {len(dfs[l])} matrices)"
        for l in labels
    )
    print(f"Runs: {threads_str}")
    print()

    for label in labels:
        print(summary_stats(dfs[label], label))
        print()

    scaling = parallel_scaling(dfs)
    if scaling:
        print(scaling)
        print()

    print(category_table(dfs))
    print(domain_table(dfs))
    print(phase_breakdown(dfs))


if __name__ == "__main__":
    main()
