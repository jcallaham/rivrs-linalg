#!/usr/bin/env python3
"""Analyze solver benchmark results from JSON benchmark files.

Supports SPRAL, MUMPS, MA27 comparison JSON files (from the corresponding
``cargo run --bin *-comparison`` binaries) and rivrs baseline JSON files
(from ``cargo run --example baseline_collection``).

Usage:
    # SPRAL comparison (unchanged from original):
    python comparisons/analyze_results.py target/benchmarks/spral/spral-benchmark-*.json

    # MUMPS standalone (absolute timings, no rivrs ratios):
    python comparisons/analyze_results.py target/benchmarks/mumps/mumps-benchmark-*.json

    # MA27 standalone:
    python comparisons/analyze_results.py target/benchmarks/ma27/ma27-benchmark-*.json

    # MUMPS vs rivrs (join with baseline to compute ratios):
    python comparisons/analyze_results.py target/benchmarks/mumps/*.json \\
        --baseline target/benchmarks/baselines/baseline-1771946638.json

    # Multi-solver side-by-side:
    python comparisons/analyze_results.py \\
        target/benchmarks/spral/spral-benchmark-1772124268.json \\
        target/benchmarks/mumps/mumps-benchmark-1772176331.json \\
        target/benchmarks/ma27/ma27-benchmark-1772158310.json

    # Auto-discover latest results from all solver directories:
    python comparisons/analyze_results.py

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
# JSON format detection and loading
# ---------------------------------------------------------------------------
def detect_solver_type(data: dict) -> str:
    """Detect JSON type: 'spral', 'mumps', 'ma27', or 'baseline'."""
    if "baselines" in data:
        return "baseline"
    if "solver_name" in data:
        return data["solver_name"].lower()  # "mumps" or "ma27"
    if "spral_version" in data or "omp_num_threads" in data:
        return "spral"
    return "unknown"


def load_spral_json(path: Path) -> pd.DataFrame:
    """Load a spral-benchmark JSON file into a comparison DataFrame.

    Returns a DataFrame with columns:
        matrix, category, n, nnz, solver_name, solver_fac, rivrs_fac, ratio,
        solver_be, rivrs_be, threads, and phase timing columns.
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

        solver_fac = spral["factor_s"]
        rivrs_fac = rivrs["factor_s"]
        ratio = rivrs_fac / solver_fac if solver_fac > 0 else float("inf")

        rows.append({
            "matrix": rec["matrix_name"],
            "category": rec.get("category", "unknown"),
            "n": rec["n"],
            "nnz": rec.get("nnz", 0),
            "solver_name": "SPRAL",
            "solver_fac": solver_fac,
            "rivrs_fac": rivrs_fac,
            "ratio": ratio,
            "solver_be": spral.get("backward_error", float("nan")),
            "rivrs_be": rivrs.get("backward_error", float("nan")),
            "threads": threads,
            "solver_analyse": spral.get("analyse_s", float("nan")),
            "rivrs_analyse": rivrs.get("analyse_s", float("nan")),
            "solver_solve": spral.get("solve_s", float("nan")),
            "rivrs_solve": rivrs.get("solve_s", float("nan")),
            "solver_total": spral.get("total_s", float("nan")),
            "rivrs_total": rivrs.get("total_s", float("nan")),
        })

    df = pd.DataFrame(rows)
    df.attrs["path"] = str(path)
    df.attrs["timestamp"] = data.get("timestamp", "")
    df.attrs["rivrs_version"] = data.get("rivrs_version", "?")
    df.attrs["solver_version"] = data.get("spral_version", "?")
    df.attrs["solver_name"] = "SPRAL"
    df.attrs["platform"] = data.get("platform", "?")
    return df


def load_generic_json(path: Path) -> pd.DataFrame:
    """Load a MUMPS or MA27 benchmark JSON file into a normalized DataFrame.

    Handles the common BenchmarkSuite schema (solver_name, solver field).
    Gracefully skips entries where solver is null (e.g. MUMPS failing on
    singular matrices like astro-ph).
    """
    with open(path) as f:
        data = json.load(f)

    solver_name = data.get("solver_name", "unknown").upper()
    threads = data.get("threads") or 1
    rows = []
    skipped = []
    for rec in data["results"]:
        solver = rec.get("solver")
        rivrs = rec.get("rivrs")

        # Skip entries where the solver failed (null solver result)
        if solver is None:
            skipped.append(rec["matrix_name"])
            continue

        solver_fac = solver["factor_s"]
        rivrs_fac = rivrs["factor_s"] if rivrs else float("nan")

        if rivrs and solver_fac > 0:
            ratio = rivrs_fac / solver_fac
        else:
            ratio = float("nan")

        row = {
            "matrix": rec["matrix_name"],
            "category": rec.get("category", "unknown"),
            "n": rec["n"],
            "nnz": rec.get("nnz", 0),
            "solver_name": solver_name,
            "solver_fac": solver_fac,
            "rivrs_fac": rivrs_fac,
            "ratio": ratio,
            "solver_be": solver.get("backward_error", float("nan")),
            "rivrs_be": rivrs.get("backward_error", float("nan")) if rivrs else float("nan"),
            "threads": threads,
            "solver_analyse": solver.get("analyse_s", float("nan")),
            "rivrs_analyse": rivrs.get("analyse_s", float("nan")) if rivrs else float("nan"),
            "solver_solve": solver.get("solve_s", float("nan")),
            "rivrs_solve": rivrs.get("solve_s", float("nan")) if rivrs else float("nan"),
            "solver_total": float("nan"),
            "rivrs_total": rivrs.get("total_s", float("nan")) if rivrs else float("nan"),
        }
        # Compute solver total from components
        ana = solver.get("analyse_s", 0)
        fac = solver.get("factor_s", 0)
        slv = solver.get("solve_s", 0)
        if ana is not None and fac is not None and slv is not None:
            row["solver_total"] = ana + fac + slv

        rows.append(row)

    if skipped:
        print(
            f"  note: {solver_name} skipped {len(skipped)} failed matrices: "
            + ", ".join(m.split("/")[-1] for m in skipped),
            file=sys.stderr,
        )

    df = pd.DataFrame(rows)
    df.attrs["path"] = str(path)
    df.attrs["timestamp"] = data.get("timestamp", "")
    df.attrs["rivrs_version"] = data.get("rivrs_version", "?")
    df.attrs["solver_version"] = data.get("solver_version", "?")
    df.attrs["solver_name"] = solver_name
    df.attrs["platform"] = data.get("platform", "?")
    return df


def load_baseline_json(path: Path) -> pd.DataFrame:
    """Load a rivrs baseline JSON file.

    Baseline times are in milliseconds; converts to seconds for consistency
    with comparison JSON files. Returns a DataFrame keyed by matrix name with
    rivrs timing columns.
    """
    with open(path) as f:
        data = json.load(f)

    rows = []
    for bl in data["baselines"]:
        rows.append({
            "matrix": bl["matrix_name"],
            "n": bl["matrix_dim"],
            "nnz": bl["matrix_nnz"],
            "rivrs_fac": bl["factor_ms"] / 1000.0,
            "rivrs_analyse": (bl.get("ordering_ms", 0) + bl.get("symbolic_ms", 0)) / 1000.0,
            "rivrs_solve": bl["solve_ms"] / 1000.0,
            "rivrs_total": bl["total_ms"] / 1000.0,
            "rivrs_be": bl["backward_error"],
        })

    df = pd.DataFrame(rows)
    df.attrs["path"] = str(path)
    df.attrs["timestamp"] = data.get("timestamp", "")
    df.attrs["solver_version"] = data.get("solver_version", "?")
    df.attrs["platform"] = data.get("platform", "?")
    return df


def load_any_json(path: Path) -> tuple[str, pd.DataFrame]:
    """Load any benchmark JSON, returning (solver_type, DataFrame)."""
    with open(path) as f:
        data = json.load(f)

    solver_type = detect_solver_type(data)
    if solver_type == "spral":
        return solver_type, load_spral_json(path)
    elif solver_type in ("mumps", "ma27"):
        return solver_type, load_generic_json(path)
    elif solver_type == "baseline":
        return solver_type, load_baseline_json(path)
    else:
        print(f"  warning: unrecognized JSON format in {path.name}", file=sys.stderr)
        return solver_type, pd.DataFrame()


# ---------------------------------------------------------------------------
# Auto-discovery
# ---------------------------------------------------------------------------
def find_benchmark_root() -> Path:
    """Find the benchmark output root directory."""
    candidates = [
        Path("target/benchmarks"),
        Path("sparse/target/benchmarks"),
    ]
    for c in candidates:
        if c.is_dir():
            return c
    return candidates[0]


def find_default_files() -> list[Path]:
    """Find the newest benchmark file from each solver directory.

    Scans spral/, mumps/, ma27/ subdirectories under target/benchmarks/
    and returns the newest file from each.
    """
    root = find_benchmark_root()
    if not root.is_dir():
        return []

    found = []
    for solver_dir in ["spral", "mumps", "ma27"]:
        d = root / solver_dir
        if not d.is_dir():
            continue
        json_files = sorted(d.glob("*.json"))
        if json_files:
            found.append(json_files[-1])  # newest by filename (timestamp suffix)

    return found


# ---------------------------------------------------------------------------
# Input loading (multi-solver aware)
# ---------------------------------------------------------------------------
def load_inputs(paths: list[Path], baseline_path: Path | None = None) -> dict[str, pd.DataFrame]:
    """Load one or more JSON files, keyed by solver label.

    For single-solver mode (all files are the same solver type), keys
    are 'T{threads}' for backward compatibility with parallel scaling.

    For multi-solver mode (files from different solvers), keys are
    solver names like 'SPRAL', 'MUMPS', 'MA27'.

    If baseline_path is provided, rivrs timing data from the baseline is
    joined into solver DataFrames that lack rivrs results.
    """
    # Load baseline if provided
    baseline_df = None
    if baseline_path:
        baseline_df = load_baseline_json(baseline_path)
        baseline_df = baseline_df.set_index("matrix")

    # Load all solver files
    solver_types = set()
    loaded: list[tuple[str, pd.DataFrame]] = []
    for p in paths:
        stype, df = load_any_json(p)
        if stype == "baseline":
            # If someone passes a baseline as a positional arg, use it as the baseline
            if baseline_df is None:
                baseline_df = df.set_index("matrix")
                print(f"  note: using {p.name} as rivrs baseline", file=sys.stderr)
            continue
        if df.empty:
            print(f"  warning: {p.name} has no results, skipping", file=sys.stderr)
            continue
        solver_types.add(stype)
        loaded.append((stype, df))

    # Join baseline data into solver DataFrames that lack rivrs results
    if baseline_df is not None:
        for i, (stype, df) in enumerate(loaded):
            if df.empty:
                continue
            has_rivrs = df["rivrs_fac"].notna().any()
            if not has_rivrs:
                df = _join_baseline(df, baseline_df)
                loaded[i] = (stype, df)

    # Determine keying strategy
    if len(solver_types) <= 1:
        # Single solver type: key by thread count (backward compat)
        results: dict[str, pd.DataFrame] = {}
        timestamps: dict[str, str] = {}
        for _stype, df in loaded:
            if df.empty:
                continue
            threads = int(df["threads"].iloc[0])
            label = f"T{threads}"
            ts = df.attrs.get("timestamp", "")
            if label not in results or ts > timestamps.get(label, ""):
                results[label] = df
                timestamps[label] = ts
        return results
    else:
        # Multi-solver: key by solver name, keep newest per solver
        results = {}
        timestamps = {}
        for _stype, df in loaded:
            if df.empty:
                continue
            solver_name = df.attrs.get("solver_name", _stype.upper())
            ts = df.attrs.get("timestamp", "")
            if solver_name not in results or ts > timestamps.get(solver_name, ""):
                results[solver_name] = df
                timestamps[solver_name] = ts
        return results


def _join_baseline(solver_df: pd.DataFrame, baseline_df: pd.DataFrame) -> pd.DataFrame:
    """Join rivrs baseline data into a solver DataFrame.

    Fills rivrs_fac, rivrs_analyse, rivrs_solve, rivrs_total, rivrs_be
    from the baseline, and computes the ratio column.
    """
    df = solver_df.copy()
    rivrs_cols = ["rivrs_fac", "rivrs_analyse", "rivrs_solve", "rivrs_total", "rivrs_be"]
    for col in rivrs_cols:
        if col in baseline_df.columns:
            mapping = baseline_df[col].to_dict()
            df[col] = df["matrix"].map(mapping)

    # Recompute ratio
    mask = df["rivrs_fac"].notna() & df["solver_fac"].notna() & (df["solver_fac"] > 0)
    df.loc[mask, "ratio"] = df.loc[mask, "rivrs_fac"] / df.loc[mask, "solver_fac"]
    return df


# ---------------------------------------------------------------------------
# Statistics — single-solver comparison (solver vs rivrs)
# ---------------------------------------------------------------------------
def summary_stats(df: pd.DataFrame, label: str) -> str:
    """Compute and format summary statistics for a comparison DataFrame."""
    solver_name = df.attrs.get("solver_name", "solver")

    # If no rivrs data, show standalone stats
    has_ratios = df["ratio"].notna().any() if "ratio" in df.columns else False
    if not has_ratios:
        return standalone_summary(df, label)

    ratios = df["ratio"].dropna()
    n = len(ratios)
    if n == 0:
        return standalone_summary(df, label)

    total_solver = df.loc[ratios.index, "solver_fac"].sum()
    total_rivrs = df.loc[ratios.index, "rivrs_fac"].sum()
    time_weighted = total_rivrs / total_solver if total_solver > 0 else float("inf")

    faster = (ratios < 0.90).sum()
    comparable = ((ratios >= 0.90) & (ratios <= 1.10)).sum()
    slower = (ratios > 1.10).sum()

    lines = [
        f"**{label} (threads={df['threads'].iloc[0]}) — {n} matrices:**",
        "",
        "| Statistic | Value |",
        "|-----------|-------|",
        f"| Median ratio | {ratios.median():.2f}x |",
        f"| Mean ratio | {ratios.mean():.2f}x |",
        f"| 25th percentile | {ratios.quantile(0.25):.2f}x |",
        f"| 75th percentile | {ratios.quantile(0.75):.2f}x |",
        f"| Time-weighted ratio | {time_weighted:.2f}x |",
        f"| Rivrs faster (<0.90x) | {faster}/{n} ({100*faster//n}%) |",
        f"| Comparable (0.90-1.10x) | {comparable}/{n} ({100*comparable//n}%) |",
        f"| {solver_name} faster (>1.10x) | {slower}/{n} ({100*slower//n}%) |",
    ]

    # Top 5 wins for each
    valid_df = df.loc[ratios.index]
    top_rivrs = valid_df.nsmallest(5, "ratio")[["matrix", "ratio"]]
    top_solver = valid_df.nlargest(5, "ratio")[["matrix", "ratio"]]

    rivrs_str = ", ".join(
        f"{row.matrix.split('/')[-1]} {row.ratio:.2f}x" for row in top_rivrs.itertuples()
    )
    solver_str = ", ".join(
        f"{row.matrix.split('/')[-1]} {row.ratio:.2f}x" for row in top_solver.itertuples()
    )
    lines += [
        "",
        f"Top rivrs wins: {rivrs_str}.",
        f"Top {solver_name} wins: {solver_str}.",
    ]

    return "\n".join(lines)


def standalone_summary(df: pd.DataFrame, label: str) -> str:
    """Format absolute timing stats when no rivrs comparison is available."""
    solver_name = df.attrs.get("solver_name", "solver")
    n = len(df)

    lines = [
        f"**{label} — {solver_name} — {n} matrices (standalone):**",
        "",
        "| Statistic | Value |",
        "|-----------|-------|",
        f"| Median factor time | {df['solver_fac'].median():.3f}s |",
        f"| Mean factor time | {df['solver_fac'].mean():.3f}s |",
        f"| Total factor time | {df['solver_fac'].sum():.3f}s |",
        f"| Min factor time | {df['solver_fac'].min():.3f}s |",
        f"| Max factor time | {df['solver_fac'].max():.3f}s |",
    ]

    be = df["solver_be"].dropna()
    if len(be) > 0:
        lines += [
            f"| Median backward error | {be.median():.2e} |",
            f"| Max backward error | {be.max():.2e} |",
        ]

    return "\n".join(lines)


def standalone_table(df: pd.DataFrame) -> str:
    """Print per-matrix absolute timings (no rivrs comparison)."""
    solver_name = df.attrs.get("solver_name", "solver")
    lines = [
        f"**Per-matrix timings ({solver_name}):**",
        "",
        f"| {'Matrix':<35} | {'n':>8} | {'nnz':>10} | {'analyse':>8} | {'factor':>9} | {'solve':>8} | {'BE':>10} |",
        f"|{'-'*37}|{'-'*10}|{'-'*12}|{'-'*10}|{'-'*11}|{'-'*10}|{'-'*12}|",
    ]
    for _, row in df.sort_values("solver_fac", ascending=False).iterrows():
        be_str = f"{row['solver_be']:.2e}" if pd.notna(row.get("solver_be")) else "N/A"
        ana_str = f"{row['solver_analyse']:.3f}" if pd.notna(row.get("solver_analyse")) else "N/A"
        slv_str = f"{row['solver_solve']:.3f}" if pd.notna(row.get("solver_solve")) else "N/A"
        lines.append(
            f"| {row['matrix']:<35} | {row['n']:>8} | {row['nnz']:>10} "
            f"| {ana_str:>8} | {row['solver_fac']:>9.3f} | {slv_str:>8} | {be_str:>10} |"
        )
    lines.append("")
    return "\n".join(lines)


def category_table(dfs: dict[str, pd.DataFrame]) -> str:
    """Tabulate speedup by matrix category (easy/hard indefinite, PD)."""
    rows = []
    for label, df in sorted(dfs.items()):
        solver_name = df.attrs.get("solver_name", "solver")
        has_ratios = "ratio" in df.columns and df["ratio"].notna().any()
        if not has_ratios:
            continue
        for cat in ["easy-indefinite", "hard-indefinite", "positive-definite"]:
            sub = df[df["category"] == cat]
            sub = sub[sub["ratio"].notna()]
            if sub.empty:
                continue
            tw = sub["rivrs_fac"].sum() / sub["solver_fac"].sum() if sub["solver_fac"].sum() > 0 else float("nan")
            rows.append({
                "threads": label,
                "solver": solver_name,
                "category": cat,
                "count": len(sub),
                "median": sub["ratio"].median(),
                "mean": sub["ratio"].mean(),
                "time_weighted": tw,
                "min": sub["ratio"].min(),
                "max": sub["ratio"].max(),
            })

    if not rows:
        return ""

    tbl = pd.DataFrame(rows)
    lines = ["**By category:**", ""]
    for threads_label in sorted(tbl["threads"].unique()):
        sub = tbl[tbl["threads"] == threads_label]
        solver = sub["solver"].iloc[0] if len(sub) > 0 else "?"
        lines.append(f"*{threads_label} ({solver}):*")
        lines.append("")
        lines.append("| Category | Count | Median | Mean | Time-wtd | Min | Max |")
        lines.append("|----------|------:|-------:|-----:|---------:|----:|----:|")
        for _, row in sub.iterrows():
            lines.append(
                f"| {row['category']} | {row['count']} "
                f"| {row['median']:.2f}x | {row['mean']:.2f}x "
                f"| {row['time_weighted']:.2f}x "
                f"| {row['min']:.2f}x | {row['max']:.2f}x |"
            )
        lines.append("")
    return "\n".join(lines)


def domain_table(dfs: dict[str, pd.DataFrame]) -> str:
    """Tabulate speedup by application domain."""
    rows = []
    for label, df in sorted(dfs.items()):
        solver_name = df.attrs.get("solver_name", "solver")
        has_ratios = "ratio" in df.columns and df["ratio"].notna().any()
        if not has_ratios:
            continue
        df = df.copy()
        df["domain"] = df["matrix"].map(lambda m: DOMAIN_MAP.get(m, "unknown"))
        for domain in sorted(df["domain"].unique()):
            sub = df[df["domain"] == domain]
            sub = sub[sub["ratio"].notna()]
            if sub.empty:
                continue
            tw = sub["rivrs_fac"].sum() / sub["solver_fac"].sum() if sub["solver_fac"].sum() > 0 else float("nan")
            rows.append({
                "threads": label,
                "solver": solver_name,
                "domain": domain,
                "count": len(sub),
                "median": sub["ratio"].median(),
                "mean": sub["ratio"].mean(),
                "time_weighted": tw,
            })

    if not rows:
        return ""

    tbl = pd.DataFrame(rows)
    lines = ["**By application domain:**", ""]
    for threads_label in sorted(tbl["threads"].unique()):
        sub = tbl[tbl["threads"] == threads_label].sort_values("median")
        solver = sub["solver"].iloc[0] if len(sub) > 0 else "?"
        lines.append(f"*{threads_label} ({solver}):*")
        lines.append("")
        lines.append("| Domain | Count | Median | Mean | Time-wtd |")
        lines.append("|--------|------:|-------:|-----:|---------:|")
        for _, row in sub.iterrows():
            lines.append(
                f"| {row['domain']} | {row['count']} "
                f"| {row['median']:.2f}x | {row['mean']:.2f}x "
                f"| {row['time_weighted']:.2f}x |"
            )
        lines.append("")
    return "\n".join(lines)


def parallel_scaling(dfs: dict[str, pd.DataFrame]) -> str:
    """Compute parallel scaling stats across the lowest and highest thread counts."""
    labels = sorted(dfs.keys(), key=lambda k: dfs[k]["threads"].iloc[0])
    if len(labels) < 2:
        return ""

    # Only makes sense for same-solver comparisons with different thread counts
    solver_names = set(df.attrs.get("solver_name", "?") for df in dfs.values())
    if len(solver_names) > 1:
        return ""

    lo_label, hi_label = labels[0], labels[-1]
    lo = dfs[lo_label]
    hi = dfs[hi_label]
    lo_t = lo["threads"].iloc[0]
    hi_t = hi["threads"].iloc[0]
    if lo_t == hi_t:
        return ""

    lo = lo.set_index("matrix")
    hi = hi.set_index("matrix")
    common = lo.index.intersection(hi.index)

    # Use rivrs_fac if available, otherwise solver_fac
    fac_col = "rivrs_fac" if lo["rivrs_fac"].notna().any() else "solver_fac"

    speedup = lo.loc[common, fac_col] / hi.loc[common, fac_col]
    speedup = speedup.dropna()
    if speedup.empty:
        return ""

    top = speedup.nlargest(5)
    top_str = ", ".join(f"{m.split('/')[-1]} {v:.2f}x" for m, v in top.items())

    solver_label = "rivrs" if fac_col == "rivrs_fac" else next(iter(solver_names))
    return (
        f"**Parallel scaling ({solver_label} T{lo_t}->T{hi_t}):** "
        f"Mean {speedup.mean():.2f}x, median {speedup.median():.2f}x. "
        f"Best: {top_str}."
    )


def phase_breakdown(dfs: dict[str, pd.DataFrame]) -> str:
    """Show analyse/factor/solve time breakdown (available from JSON)."""
    solver_name_first = None
    for df in dfs.values():
        solver_name_first = df.attrs.get("solver_name", "solver")
        break

    has_any = False
    for df in dfs.values():
        if "solver_analyse" in df.columns and df["solver_analyse"].notna().any():
            has_any = True
            break
    if not has_any:
        return ""

    lines = ["**Phase breakdown (median, seconds):**", ""]
    for label in sorted(dfs.keys()):
        df = dfs[label]
        solver_name = df.attrs.get("solver_name", "solver")
        if "solver_analyse" not in df.columns:
            continue

        has_rivrs = "rivrs_analyse" in df.columns and df["rivrs_fac"].notna().any()

        lines.append(f"*{label}:*")
        lines.append("")

        if has_rivrs:
            lines.append(f"| Phase | {solver_name} median | rivrs median | Ratio |")
            lines.append("|-------|------------:|-----------:|------:|")
            for phase_s, phase_r, name in [
                ("solver_analyse", "rivrs_analyse", "Analyse"),
                ("solver_fac", "rivrs_fac", "Factor"),
                ("solver_solve", "rivrs_solve", "Solve"),
                ("solver_total", "rivrs_total", "Total"),
            ]:
                s = df[phase_s].median()
                r = df[phase_r].median()
                ratio = r / s if s > 0 else float("inf")
                lines.append(f"| {name} | {s:.3f} | {r:.3f} | {ratio:.2f}x |")
        else:
            lines.append(f"| Phase | {solver_name} median |")
            lines.append("|-------|------------:|")
            for phase_s, name in [
                ("solver_analyse", "Analyse"),
                ("solver_fac", "Factor"),
                ("solver_solve", "Solve"),
                ("solver_total", "Total"),
            ]:
                s = df[phase_s].median()
                lines.append(f"| {name} | {s:.3f} |")
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Multi-solver comparison
# ---------------------------------------------------------------------------
def multi_solver_table(dfs: dict[str, pd.DataFrame]) -> str:
    """Side-by-side factor time comparison across multiple solvers.

    Joins DataFrames from different solvers on matrix name and prints
    a combined table with one column per solver.
    """
    solver_names = sorted(dfs.keys())
    if len(solver_names) < 2:
        return ""

    # Build a master matrix info table (n, nnz, category) from all sources
    matrix_info = {}
    for name in solver_names:
        df = dfs[name]
        for _, row in df.iterrows():
            m = row["matrix"]
            if m not in matrix_info:
                matrix_info[m] = {"matrix": m}
            # Fill in from whichever solver has the data
            for col in ["n", "nnz", "category"]:
                if col in row.index and pd.notna(row[col]) and col not in matrix_info[m]:
                    matrix_info[m][col] = row[col]
    info_df = pd.DataFrame(matrix_info.values())

    # Collect factor times and rivrs data per solver
    merged = info_df.copy()
    for name in solver_names:
        df = dfs[name]
        solver_cols = df[["matrix", "solver_fac"]].rename(columns={"solver_fac": f"{name}_fac"})
        # Include rivrs if available and not yet merged
        if "rivrs_fac" in df.columns and df["rivrs_fac"].notna().any() and "rivrs_fac" not in merged.columns:
            solver_cols = solver_cols.merge(df[["matrix", "rivrs_fac"]], on="matrix")
        merged = merged.merge(solver_cols, on="matrix", how="outer")

    # Backfill n/nnz for any outer-joined rows
    for col in ["n", "nnz"]:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0).astype(int)

    if merged is None or merged.empty:
        return ""

    # Check if we have rivrs data
    has_rivrs = "rivrs_fac" in merged.columns and merged["rivrs_fac"].notna().any()

    # Build header
    fac_headers = [f"{name:>10}" for name in solver_names]
    if has_rivrs:
        fac_headers.append(f"{'rivrs':>10}")
    header = f"| {'Matrix':<35} | {'n':>8} | " + " | ".join(fac_headers) + " |"
    sep = f"|{'-'*37}|{'-'*10}|" + "|".join("-" * 12 for _ in fac_headers) + "|"

    lines = ["**Multi-solver factor time comparison (seconds):**", "", header, sep]

    for _, row in merged.sort_values("n").iterrows():
        n_str = f"{int(row['n']):>8}" if pd.notna(row.get("n")) else f"{'?':>8}"
        vals = []
        for name in solver_names:
            v = row.get(f"{name}_fac")
            vals.append(f"{v:>10.3f}" if pd.notna(v) else f"{'FAIL':>10}")
        if has_rivrs:
            v = row.get("rivrs_fac")
            vals.append(f"{v:>10.3f}" if pd.notna(v) else f"{'N/A':>10}")
        lines.append(f"| {row['matrix']:<35} | {n_str} | " + " | ".join(vals) + " |")

    lines.append("")

    # Summary: ratio of each solver vs rivrs (if rivrs available)
    if has_rivrs:
        lines.append("**Factor time ratio vs rivrs (rivrs_fac / solver_fac):**")
        lines.append("")
        lines.append(f"| {'Solver':<10} | {'Median':>8} | {'Mean':>8} | {'Time-wtd':>10} | {'Count':>6} |")
        lines.append(f"|{'-'*12}|{'-'*10}|{'-'*10}|{'-'*12}|{'-'*8}|")
        for name in solver_names:
            valid = merged[[f"{name}_fac", "rivrs_fac"]].dropna()
            if valid.empty:
                continue
            ratios = valid["rivrs_fac"] / valid[f"{name}_fac"]
            ratios = ratios.replace([np.inf, -np.inf], np.nan).dropna()
            if ratios.empty:
                continue
            tw = valid["rivrs_fac"].sum() / valid[f"{name}_fac"].sum()
            # Use appropriate precision for time-weighted ratio
            tw_str = f"{tw:.2f}" if tw >= 0.01 else f"{tw:.3f}"
            lines.append(
                f"| {name:<10} | {ratios.median():>7.2f}x | {ratios.mean():>7.2f}x "
                f"| {tw_str:>9}x | {len(ratios):>6} |"
            )
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Analyze solver benchmark JSON results (SPRAL, MUMPS, MA27, rivrs baselines).",
        epilog=(
            "If no files are given, auto-discovers the newest file from each solver directory.\n\n"
            "Examples:\n"
            "  # SPRAL comparison (unchanged):\n"
            "  python comparisons/analyze_results.py target/benchmarks/spral/*.json\n\n"
            "  # MUMPS standalone:\n"
            "  python comparisons/analyze_results.py target/benchmarks/mumps/*.json\n\n"
            "  # MUMPS with rivrs baseline:\n"
            "  python comparisons/analyze_results.py target/benchmarks/mumps/*.json \\\n"
            "      --baseline target/benchmarks/baselines/baseline-1771946638.json\n\n"
            "  # Multi-solver side-by-side:\n"
            "  python comparisons/analyze_results.py \\\n"
            "      target/benchmarks/spral/spral-benchmark-1772124268.json \\\n"
            "      target/benchmarks/mumps/mumps-benchmark-1772176331.json \\\n"
            "      target/benchmarks/ma27/ma27-benchmark-1772158310.json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "files",
        nargs="*",
        type=Path,
        help="One or more benchmark JSON files (SPRAL, MUMPS, MA27, or baseline)",
    )
    parser.add_argument(
        "--baseline",
        type=Path,
        default=None,
        help="Rivrs baseline JSON file for computing ratios with MUMPS/MA27",
    )
    args = parser.parse_args()

    # Resolve input files
    if args.files:
        json_paths = [p for p in args.files if p.suffix == ".json" and p.exists()]
        if not json_paths:
            print(f"Error: no valid .json files found in {args.files}", file=sys.stderr)
            sys.exit(1)
    else:
        json_paths = find_default_files()
        if not json_paths:
            # Fall back to old behavior: look for spral dir
            bench_dir = find_benchmark_root() / "spral"
            if bench_dir.is_dir():
                json_paths = sorted(bench_dir.glob("spral-benchmark-*.json"))
                if json_paths:
                    json_paths = [json_paths[-1]]
        if not json_paths:
            print("Error: no benchmark files found. Pass JSON file(s) explicitly.", file=sys.stderr)
            sys.exit(1)
        names = [p.name for p in json_paths]
        print(f"(auto-discovered: {', '.join(names)})", file=sys.stderr)

    dfs = load_inputs(json_paths, baseline_path=args.baseline)

    if not dfs:
        print("No results found in input files", file=sys.stderr)
        sys.exit(1)

    # Determine if this is a multi-solver comparison
    solver_names = set()
    for df in dfs.values():
        sn = df.attrs.get("solver_name", "?")
        solver_names.add(sn)

    is_multi_solver = len(solver_names) > 1

    # Print header
    first_df = next(iter(dfs.values()))
    platform = first_df.attrs.get("platform", "")

    if is_multi_solver:
        print("=" * 72)
        print("Multi-Solver Benchmark Analysis")
        print("=" * 72)
        if platform:
            print(f"Platform: {platform}")
        solvers_str = ", ".join(sorted(solver_names))
        print(f"Solvers: {solvers_str}")
        for label, df in sorted(dfs.items()):
            print(f"  {label}: {len(df)} matrices, version {df.attrs.get('solver_version', '?')}")
        print()

        # Multi-solver table
        print(multi_solver_table(dfs))

        # Per-solver summaries (if any have ratios)
        for label in sorted(dfs.keys()):
            df = dfs[label]
            has_ratios = "ratio" in df.columns and df["ratio"].notna().any()
            if has_ratios:
                print(summary_stats(df, label))
                print()

    else:
        solver_name = next(iter(solver_names))
        has_any_ratios = any(
            "ratio" in df.columns and df["ratio"].notna().any()
            for df in dfs.values()
        )

        print("=" * 72)
        if has_any_ratios:
            print(f"{solver_name} vs rivrs — Benchmark Analysis")
        else:
            print(f"{solver_name} — Benchmark Analysis")
        print("=" * 72)

        if platform:
            print(f"Platform: {platform}")
        rivrs_v = first_df.attrs.get("rivrs_version", "")
        solver_v = first_df.attrs.get("solver_version", "?")
        if has_any_ratios and rivrs_v:
            print(f"rivrs: {rivrs_v}  {solver_name}: {solver_v}")
        elif solver_v:
            print(f"{solver_name}: {solver_v}")
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

        # Standalone table for solver-only data
        if not has_any_ratios:
            for label in labels:
                print(standalone_table(dfs[label]))

        scaling = parallel_scaling(dfs)
        if scaling:
            print(scaling)
            print()

        cat = category_table(dfs)
        if cat:
            print(cat)
        dom = domain_table(dfs)
        if dom:
            print(dom)
        phase = phase_breakdown(dfs)
        if phase:
            print(phase)


if __name__ == "__main__":
    main()
