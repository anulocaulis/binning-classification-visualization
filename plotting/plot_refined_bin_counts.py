#!/usr/bin/env python3
"""
plot_refined_bin_counts.py

Summarize MetaWRAP refined bins and VAMB bins in one table/plot and optionally
overlay CheckM2 HQ/MQ counts when quality_report.tsv files are available.

If CheckM2 output is missing for a refined set, plotting still succeeds and the
set is shown with total bins only.

Usage:
  python scripts/data_visualization/plot_refined_bin_counts.py
"""
# Standard library imports
import argparse
import glob
import os
import shlex
import subprocess
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _choose_column(columns: List[str], candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in columns:
            return c
    return None

# Assembler to assembly type mapping
ASSEMBLER_TO_ASSEMBLY_TYPE: Dict[str, str] = {
    "flye": "long_read",
    "metamdbg": "long_read",
    "metaspades_hybrid": "hybrid_read",
    "opera-ms": "hybrid_read",
    "opera_ms": "hybrid_read",
    "metaspades": "short_read",
    "megahit": "short_read",
    "idbaud": "short_read",
    "spades": "short_read",
    "final": "short_read",
    "unknown": "short_read",
}

ASSEMBLY_TYPE_ORDER = ["long_read", "short_read", "hybrid_read"]
SAMPLE_ORDER = ["S1", "S2", "S5"]

ASSEMBLER_LABEL_ALIASES: Dict[str, str] = {
    "flye": "Flye",
    "metamdbg": "metaMDBG",
    "idbaud": "IDBA-UD",
    "megahit": "MEGAHIT",
    "metaspades": "metaSPAdes",
    "metaspades_hybrid": "metaSPAdes Hybrid",
    "opera-ms": "Opera-MS",
    "opera_ms": "Opera-MS",
    "spades": "SPAdes",
    "final": "Final",
}

DEFAULT_COLORS: Dict[str, List[str]] = {
    "short_read": ["#A1D99B", "#74C476", "#238B45"],
    "long_read": ["#C6DBEF", "#41B6C4", "#08519C"],
    "hybrid_read": ["#D0D1E6", "#6A51A3"],
}

DEFAULT_CONFIG_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "configs", "plotting_pipeline_config.sh")
)


# Load shared config defaults for roots, sample order, and output locations.
def load_shell_config(config_path: str) -> Dict[str, object]:
    if not os.path.isfile(config_path):
        raise SystemExit(f"ERROR: config file not found: {config_path}")

    scalar_keys = [
        "PROJECT_DATA_DIR",
        "PROJECT_PLOTTING_ROOT",
        "PROJECT_REFINED_BIN_COUNTS_OUTDIR",
        "PROJECT_REFINED_BINS_ROOT",
        "PROJECT_VAMB_OUTPUT_ROOT",
        "PROJECT_CHECKM2_REFINED_ROOT",
        "PROJECT_CHECKM2_VAMB_ROOT",
    ]
    array_keys = ["PROJECT_SAMPLE_NAMES"]

    command_parts = [f"source {shlex.quote(os.path.abspath(config_path))}"]
    for key in scalar_keys:
        command_parts.append(f'printf "%s\t%s\\0" {shlex.quote(key)} "${{{key}:-}}"')
    for key in array_keys:
        command_parts.append(f'printf "%s\t%s\\0" {shlex.quote(key)} "${{{key}[*]:-}}"')

    command = "set -euo pipefail; " + "; ".join(command_parts)
    result = subprocess.run(["bash", "-lc", command], check=True, capture_output=True)

    config: Dict[str, object] = {}
    for record in result.stdout.split(b"\0"):
        if not record:
            continue
        key, value = record.decode("utf-8", errors="replace").split("\t", 1)
        if key in array_keys:
            config[key] = [item for item in value.split() if item]
        else:
            config[key] = value
    return config


def _first_nonempty(*values: Optional[str], fallback: str) -> str:
    for value in values:
        if value:
            return value
    return fallback

# Assembly type classification and display helpers.
def classify_assembly_type(assembler: str) -> str:
    return ASSEMBLER_TO_ASSEMBLY_TYPE.get(assembler.strip().lower(), "short_read")


def display_label(assembler: str) -> str:
    return ASSEMBLER_LABEL_ALIASES.get(assembler, assembler)


def assembler_sort_key(assembler: str) -> Tuple[int, str]:
    atype = classify_assembly_type(assembler)
    rank = ASSEMBLY_TYPE_ORDER.index(atype) if atype in ASSEMBLY_TYPE_ORDER else len(ASSEMBLY_TYPE_ORDER)
    return (rank, assembler)


def assign_assembler_colors(assemblers: List[str]) -> Dict[str, str]:
    color_map: Dict[str, str] = {}
    for atype in ASSEMBLY_TYPE_ORDER:
        members = sorted([a for a in assemblers if classify_assembly_type(a) == atype])
        palette = DEFAULT_COLORS.get(atype, ["#888888"])
        for idx, assembler in enumerate(members):
            color_map[assembler] = palette[idx % len(palette)]
    return color_map


def make_assembly_type_legend_handles() -> List[object]:
    from matplotlib.patches import Patch

    handles = []
    for assembly_type in ASSEMBLY_TYPE_ORDER:
        if assembly_type in DEFAULT_COLORS:
            palette = DEFAULT_COLORS[assembly_type]
            handles.append(
                Patch(facecolor=palette[min(1, len(palette) - 1)], label=assembly_type.replace("_", " "))
            )
    return handles


def make_workflow_legend_handles() -> List[object]:
    from matplotlib.patches import Patch

    return [
        Patch(facecolor="#cccccc", edgecolor="black", hatch="", label="metawrap"),
        Patch(facecolor="#cccccc", edgecolor="black", hatch="//", label="vamb"),
    ]


def parse_set_name(set_name: str) -> Tuple[str, str]:
    stem = set_name
    for suffix in ["_metawrap_bins_refined", "_metawrap_bins", "_vamb"]:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break

    parts = stem.rsplit("_", 1)
    if len(parts) == 2 and parts[1].lower().startswith("s"):
        return parts[0], parts[1].upper()
    return stem, ""


# Locate the newest CheckM2 report for a given set directory.
def find_latest_quality_report(checkm2_root: str, set_name: str) -> Optional[str]:
    patterns = [
        os.path.join(checkm2_root, set_name, "quality_report.tsv"),
        os.path.join(checkm2_root, set_name, "*", "quality_report.tsv"),
        os.path.join(checkm2_root, set_name, "**", "quality_report.tsv"),
    ]

    candidates: List[str] = []
    for pattern in patterns:
        candidates.extend(glob.glob(pattern, recursive=True))

    candidates = [c for c in candidates if os.path.isfile(c)]
    if not candidates:
        return None

    candidates = sorted(set(candidates), key=lambda p: os.path.getmtime(p), reverse=True)
    return candidates[0]


# Count HQ, MQ, and LQ bins from a CheckM2 quality report.
def count_checkm2_quality(report_tsv: str) -> Tuple[int, int, int, int]:
    df = pd.read_csv(report_tsv, sep="\t")

    # Try both CheckM2 v1 (Completeness) and v2 (Completeness_General) column names
    comp_col = None
    if "Completeness_General" in df.columns:
        comp_col = "Completeness_General"
    elif "Completeness" in df.columns:
        comp_col = "Completeness"
    
    cont_col = "Contamination" if "Contamination" in df.columns else None
    if comp_col is None or cont_col is None:
        return 0, 0, 0, 0

    c = pd.to_numeric(df[comp_col], errors="coerce").fillna(0)
    x = pd.to_numeric(df[cont_col], errors="coerce").fillna(0)

    total = int(len(df))
    hq = int(((c >= 90) & (x <= 5)).sum())
    mq = int(((c >= 50) & (c < 90) & (x <= 10)).sum())
    lq = int(total - hq - mq)
    return total, hq, mq, lq


# Load per-bin CheckM2 rows with origin metadata for downstream plots.
def load_checkm2_bins_with_origin(
    report_tsv: str,
    set_name: str,
    workflow: str,
    binner: str,
    assembler: str,
    sample: str,
) -> pd.DataFrame:
    raw = pd.read_csv(report_tsv, sep="\t")
    if raw.empty:
        return pd.DataFrame()

    comp_col = _choose_column(raw.columns.tolist(), ["Completeness_General", "Completeness"])
    cont_col = _choose_column(raw.columns.tolist(), ["Contamination"])
    name_col = _choose_column(raw.columns.tolist(), ["Name", "Bin Id", "Bin", "bin", "Genome"])
    size_col = _choose_column(raw.columns.tolist(), ["Genome_Size", "Genome size"])
    gc_col = _choose_column(raw.columns.tolist(), ["GC_Content", "GC"])
    contig_col = _choose_column(raw.columns.tolist(), ["Total_Contigs", "N_contigs"])
    n50_col = _choose_column(raw.columns.tolist(), ["Contig_N50", "N50"])
    max_contig_col = _choose_column(
        raw.columns.tolist(),
        ["Max_Contig_Length", "Longest_Contig", "Max Contig Length", "Longest Contig"],
    )

    if comp_col is None or cont_col is None:
        return pd.DataFrame()

    comp = pd.to_numeric(raw[comp_col], errors="coerce")
    cont = pd.to_numeric(raw[cont_col], errors="coerce")

    quality = np.where(
        (comp >= 90) & (cont <= 5),
        "HQ",
        np.where((comp >= 50) & (cont <= 10), "MQ", "LQ"),
    )

    out = pd.DataFrame(
        {
            "set_name": set_name,
            "workflow": workflow,
            "binner": binner,
            "assembler": assembler,
            "sample": sample,
            "bin_name": raw[name_col].astype(str) if name_col is not None else [f"bin_{i+1}" for i in range(len(raw))],
            "completeness": comp,
            "contamination": cont,
            "quality_tier": quality,
            "quality_report": report_tsv,
        }
    )

    if size_col is not None:
        out["genome_size"] = pd.to_numeric(raw[size_col], errors="coerce")
    if gc_col is not None:
        out["gc_content"] = pd.to_numeric(raw[gc_col], errors="coerce")
    if contig_col is not None:
        out["total_contigs"] = pd.to_numeric(raw[contig_col], errors="coerce")
    if n50_col is not None:
        out["contig_n50"] = pd.to_numeric(raw[n50_col], errors="coerce")
    if max_contig_col is not None:
        out["max_contig_length"] = pd.to_numeric(raw[max_contig_col], errors="coerce")

    return out


# Count VAMB bins from FASTA files or a clusters.tsv fallback.
def count_vamb_bins(vamb_bins_dir: str, extensions: List[str]) -> Tuple[int, str]:
    fasta_count = 0
    if os.path.isdir(vamb_bins_dir):
        for ext in extensions:
            fasta_count += len(glob.glob(os.path.join(vamb_bins_dir, f"*{ext}")))
    if fasta_count > 0:
        return int(fasta_count), "vamb_bins_fasta"

    clusters_tsv = os.path.join(vamb_bins_dir, "clusters.tsv")
    if not os.path.isfile(clusters_tsv):
        return 0, "missing"

    clusters = set()
    with open(clusters_tsv, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            first = line.split("\t", 1)[0].strip()
            if not first:
                continue
            if first.lower() in {"cluster", "clustername", "bin", "label"}:
                continue
            clusters.add(first)
    return len(clusters), "clusters_tsv"


# Build the set-level summary table used by the overview plots.
def collect_all_bins_and_checkm2(
    refined_root: str,
    vamb_root: str,
    checkm2_metawrap_root: str,
    checkm2_vamb_root: str,
    bins_subdir: str,
    extensions: List[str],
    include_vamb: bool,
) -> pd.DataFrame:
    rows = []

    for entry in sorted(os.scandir(refined_root), key=lambda e: e.name):
        if not entry.is_dir() or not entry.name.endswith("_metawrap_bins_refined"):
            continue

        set_name = entry.name
        assembler, sample = parse_set_name(set_name)

        bins_dir = os.path.join(entry.path, bins_subdir)
        total_bins = 0
        if os.path.isdir(bins_dir):
            for ext in extensions:
                total_bins += len(glob.glob(os.path.join(bins_dir, f"*{ext}")))

        report_path = find_latest_quality_report(checkm2_metawrap_root, set_name)
        if report_path is not None:
            q_total, hq, mq, lq = count_checkm2_quality(report_path)
        else:
            q_total, hq, mq, lq = (np.nan, np.nan, np.nan, np.nan)

        rows.append(
            {
                "set_name": set_name,
                "assembler": assembler,
                "sample": sample,
                "workflow": "metawrap",
                "bins_dir": bins_dir,
                "refined_total_bins": int(total_bins),
                "refined_subdir": bins_subdir,
                "count_method": "fasta",
                "checkm2_total_bins": q_total,
                "hq_bins": hq,
                "mq_bins": mq,
                "lq_bins": lq,
                "checkm2_report": report_path if report_path is not None else "",
                "checkm2_available": int(report_path is not None),
            }
        )

    if include_vamb and os.path.isdir(vamb_root):
        for entry in sorted(os.scandir(vamb_root), key=lambda e: e.name):
            if not entry.is_dir() or not entry.name.endswith("_vamb"):
                continue

            set_name = entry.name
            assembler, sample = parse_set_name(set_name)

            bins_dir = os.path.join(entry.path, "vamb_bins")
            total_bins, count_method = count_vamb_bins(bins_dir, extensions)

            # Use only the configured VAMB CheckM2 root (expected to be complete runs).
            report_path = None
            if os.path.isdir(checkm2_vamb_root):
                report_path = find_latest_quality_report(checkm2_vamb_root, set_name)
            if report_path is not None:
                q_total, hq, mq, lq = count_checkm2_quality(report_path)
                # For full extraction: use CheckM2 count as total (only evaluated bins)
                total_bins = int(q_total) if not pd.isna(q_total) else total_bins
                count_method = "checkm2_evaluated"
            else:
                q_total, hq, mq, lq = (np.nan, np.nan, np.nan, np.nan)

            rows.append(
                {
                    "set_name": set_name,
                    "assembler": assembler,
                    "sample": sample,
                    "workflow": "vamb",
                    "bins_dir": bins_dir,
                    "refined_total_bins": int(total_bins),
                    "refined_subdir": "vamb_bins",
                    "count_method": count_method,
                    "checkm2_total_bins": q_total,
                    "hq_bins": hq,
                    "mq_bins": mq,
                    "lq_bins": lq,
                    "checkm2_report": report_path if report_path is not None else "",
                    "checkm2_available": int(report_path is not None),
                }
            )

    return pd.DataFrame(rows)


# Build the bin-level CheckM2 table with workflow and sample metadata.
def collect_checkm2_bins_aggregated(
    refined_root: str,
    vamb_root: str,
    checkm2_metawrap_root: str,
    checkm2_vamb_root: str,
    include_vamb: bool,
) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []

    for entry in sorted(os.scandir(refined_root), key=lambda e: e.name):
        if not entry.is_dir() or not entry.name.endswith("_metawrap_bins_refined"):
            continue
        set_name = entry.name
        assembler, sample = parse_set_name(set_name)
        report_path = find_latest_quality_report(checkm2_metawrap_root, set_name)
        if report_path is None:
            continue
        df_set = load_checkm2_bins_with_origin(
            report_tsv=report_path,
            set_name=set_name,
            workflow="metawrap",
            binner="metawrap",
            assembler=assembler,
            sample=sample,
        )
        if not df_set.empty:
            rows.append(df_set)

    if include_vamb and os.path.isdir(vamb_root):
        for entry in sorted(os.scandir(vamb_root), key=lambda e: e.name):
            if not entry.is_dir() or not entry.name.endswith("_vamb"):
                continue
            set_name = entry.name
            assembler, sample = parse_set_name(set_name)

            report_path = None
            if os.path.isdir(checkm2_vamb_root):
                report_path = find_latest_quality_report(checkm2_vamb_root, set_name)
            if report_path is None:
                continue

            df_set = load_checkm2_bins_with_origin(
                report_tsv=report_path,
                set_name=set_name,
                workflow="vamb",
                binner="vamb",
                assembler=assembler,
                sample=sample,
            )
            if not df_set.empty:
                rows.append(df_set)

    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    return out


# Apply optional workflow, assembler, sample, and binner filters.
def apply_bin_filters(
    df: pd.DataFrame,
    workflows: List[str],
    assemblers: List[str],
    samples: List[str],
    binners: List[str],
) -> pd.DataFrame:
    out = df.copy()
    if workflows:
        s = {x.strip().lower() for x in workflows if x.strip()}
        out = out[out["workflow"].str.lower().isin(s)]
    if assemblers:
        s = {x.strip().lower() for x in assemblers if x.strip()}
        out = out[out["assembler"].str.lower().isin(s)]
    if samples:
        s = {x.strip().upper() for x in samples if x.strip()}
        out = out[out["sample"].str.upper().isin(s)]
    if binners:
        s = {x.strip().lower() for x in binners if x.strip()}
        out = out[out["binner"].str.lower().isin(s)]
    return out.reset_index(drop=True)


# Plot completeness and contamination means for the bin-level table.
def plot_average_completeness_contamination(df_bins: pd.DataFrame, outdir: str) -> None:
    if df_bins.empty:
        return
    grp = (
        df_bins.groupby(["workflow", "assembler", "sample"], dropna=False)
        .agg(
            mean_completeness=("completeness", "mean"),
            mean_contamination=("contamination", "mean"),
            n_bins=("bin_name", "count"),
        )
        .reset_index()
    )
    grp["label"] = grp.apply(lambda r: f"{display_label(r['assembler'])}\n({r['sample']})\n{r['workflow']}", axis=1)

    x = np.arange(len(grp))
    fig, ax1 = plt.subplots(figsize=(max(12, len(grp) * 0.8), 6))
    ax2 = ax1.twinx()

    ax1.bar(x - 0.2, grp["mean_completeness"], width=0.4, color="#2b8cbe", label="Mean completeness")
    ax2.bar(x + 0.2, grp["mean_contamination"], width=0.4, color="#de2d26", alpha=0.7, label="Mean contamination")

    ax1.set_xticks(x)
    ax1.set_xticklabels(grp["label"].tolist(), rotation=35, ha="right", fontsize=9)
    ax1.set_ylabel("Mean completeness (%)", color="#2b8cbe")
    ax2.set_ylabel("Mean contamination (%)", color="#de2d26")
    ax1.set_title("Average CheckM2 completeness and contamination by origin")
    ax1.grid(axis="y", linestyle="--", alpha=0.4)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper left")

    fig.tight_layout()
    out_path = os.path.join(outdir, "checkm2_avg_completeness_contamination.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.close(fig)


# Plot HQ/MQ/LQ fractions by origin.
def plot_quality_tier_fractions(df_bins: pd.DataFrame, outdir: str) -> None:
    if df_bins.empty:
        return
    cnt = (
        df_bins.groupby(["workflow", "assembler", "sample", "quality_tier"], dropna=False)
        .size()
        .reset_index(name="n")
    )
    totals = cnt.groupby(["workflow", "assembler", "sample"], dropna=False)["n"].sum().reset_index(name="total")
    cnt = cnt.merge(totals, on=["workflow", "assembler", "sample"], how="left")
    cnt["fraction"] = cnt["n"] / cnt["total"].replace(0, np.nan)

    groups = (
        cnt[["workflow", "assembler", "sample"]]
        .drop_duplicates()
        .sort_values(["workflow", "assembler", "sample"])
        .reset_index(drop=True)
    )
    labels = [f"{display_label(a)}\n({s})\n{w}" for w, a, s in groups[["workflow", "assembler", "sample"]].values]
    x = np.arange(len(groups))

    tiers = ["LQ", "MQ", "HQ"]
    colors = {"LQ": "#fdae6b", "MQ": "#9ecae1", "HQ": "#74c476"}
    hatches = {"LQ": "xx", "MQ": "//", "HQ": ""}

    fig, ax = plt.subplots(figsize=(max(12, len(groups) * 0.9), 6))
    bottom = np.zeros(len(groups))
    for tier in tiers:
        vals = []
        for _, row in groups.iterrows():
            m = cnt[
                (cnt["workflow"] == row["workflow"])
                & (cnt["assembler"] == row["assembler"])
                & (cnt["sample"] == row["sample"])
                & (cnt["quality_tier"] == tier)
            ]
            vals.append(float(m["fraction"].iloc[0]) if not m.empty else 0.0)
        ax.bar(x, vals, bottom=bottom, color=colors[tier], hatch=hatches[tier], edgecolor="white", linewidth=0.6, label=tier)
        bottom = bottom + np.array(vals)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of bins")
    ax.set_title("Quality tier composition by origin (fractional)")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.legend(title="Tier", loc="upper left")

    fig.tight_layout()
    out_path = os.path.join(outdir, "checkm2_quality_tier_fractions.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.close(fig)


# Plot filtered bin metric summaries for one threshold pair.
def _plot_filtered_bin_metrics(
    df_bins: pd.DataFrame,
    outdir: str,
    min_completeness: float = 50.0,
    max_contamination: float = 10.0,
    tag: str = "gt50_lt10",
) -> None:
    """
    Plot and summarize bin-level metrics after quality filtering.

    Filter criterion (strict, per request):
      - completeness > min_completeness
      - contamination < max_contamination
    """
    if df_bins.empty:
        return

    filt = df_bins[
        (pd.to_numeric(df_bins["completeness"], errors="coerce") > min_completeness)
        & (pd.to_numeric(df_bins["contamination"], errors="coerce") < max_contamination)
    ].copy()

    if filt.empty:
        print(
            f"WARNING: no bins passed filter completeness>{min_completeness} and contamination<{max_contamination}"
        )
        return

    group_cols = ["workflow", "assembler", "sample"]
    summary = (
        filt.groupby(group_cols, dropna=False)
        .agg(
            bin_count=("bin_name", "count"),
            mean_genome_size=("genome_size", "mean"),
            median_genome_size=("genome_size", "median"),
            mean_total_contigs_per_bin=("total_contigs", "mean"),
            median_total_contigs_per_bin=("total_contigs", "median"),
            mean_max_contig_length=("max_contig_length", "mean"),
            median_max_contig_length=("max_contig_length", "median"),
        )
        .reset_index()
    )

    sample_rank = {s: i for i, s in enumerate(SAMPLE_ORDER)}
    workflow_rank = {"metawrap": 0, "vamb": 1}
    summary["_type_rank"] = summary["assembler"].map(assembler_sort_key).map(lambda x: x[0])
    summary["_sample_rank"] = summary["sample"].map(lambda s: sample_rank.get(s, 99))
    summary["_workflow_rank"] = summary["workflow"].map(lambda w: workflow_rank.get(w, 99))
    summary = summary.sort_values(
        ["_type_rank", "assembler", "_sample_rank", "_workflow_rank"]
    ).reset_index(drop=True)

    out_summary = os.path.join(outdir, f"filtered_{tag}_bin_metrics_summary.tsv")
    summary.drop(columns=["_type_rank", "_sample_rank", "_workflow_rank"]).to_csv(
        out_summary,
        sep="\t",
        index=False,
    )
    print(f"Saved: {out_summary}")

    labels = [
        f"{display_label(a)}\n({s})\n{w}"
        for w, a, s in summary[["workflow", "assembler", "sample"]].values
    ]
    x = np.arange(len(summary))
    assemblers = summary["assembler"].unique().tolist()
    color_map = assign_assembler_colors(assemblers)
    group_colors = [color_map.get(a, "#888888") for a in summary["assembler"]]
    workflow_hatch = {"metawrap": "", "vamb": "//"}

    # 1) Counts of bins passing the filter.
    fig, ax = plt.subplots(figsize=(max(12, len(summary) * 0.9), 6))
    for idx, row in summary.iterrows():
        ax.bar(
            x[idx],
            row["bin_count"],
            color=group_colors[idx],
            edgecolor="white",
            linewidth=0.6,
            hatch=workflow_hatch.get(row["workflow"], ""),
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Bin count")
    ax.set_title(
        f"Bins passing filter: completeness>{min_completeness:g}, contamination<{max_contamination:g}"
    )
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.legend(
        handles=make_assembly_type_legend_handles() + make_workflow_legend_handles(),
        title="Assembly type / Workflow",
        loc="upper left",
        fontsize=8,
    )
    fig.tight_layout()
    out_path = os.path.join(outdir, f"filtered_{tag}_bin_counts.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.close(fig)

    # 2) Total contigs per bin distribution.
    if "total_contigs" in filt.columns and filt["total_contigs"].notna().any():
        data = []
        for _, row in summary.iterrows():
            vals = filt[
                (filt["workflow"] == row["workflow"])
                & (filt["assembler"] == row["assembler"])
                & (filt["sample"] == row["sample"])
            ]["total_contigs"].dropna()
            data.append(vals.values)

        fig, ax = plt.subplots(figsize=(max(12, len(summary) * 0.9), 6))
        bp = ax.boxplot(data, positions=x + 1, widths=0.6, showfliers=False, patch_artist=True)
        for idx, box in enumerate(bp["boxes"]):
            box.set_facecolor(group_colors[idx])
            box.set_alpha(0.75)
            box.set_edgecolor("black")
        ax.set_xticks(x + 1)
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
        ax.set_ylabel("Total contigs per bin")
        ax.set_title(
            f"Total contigs per bin for filtered bins (completeness>{min_completeness:g}, contamination<{max_contamination:g})"
        )
        ax.grid(axis="y", linestyle="--", alpha=0.4)
        ax.legend(
            handles=make_assembly_type_legend_handles(),
            title="Assembly type",
            loc="upper left",
            fontsize=8,
        )
        fig.tight_layout()
        out_path = os.path.join(outdir, f"filtered_{tag}_total_contigs_per_bin.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {out_path}")
        plt.close(fig)
    else:
        print("WARNING: total_contigs column missing or empty; skipping total-contigs-per-bin plot")

    # 3) Max contig length distribution.
    if "max_contig_length" in filt.columns and filt["max_contig_length"].notna().any():
        data = []
        for _, row in summary.iterrows():
            vals = filt[
                (filt["workflow"] == row["workflow"])
                & (filt["assembler"] == row["assembler"])
                & (filt["sample"] == row["sample"])
            ]["max_contig_length"].dropna()
            data.append(vals.values)

        fig, ax = plt.subplots(figsize=(max(12, len(summary) * 0.9), 6))
        bp = ax.boxplot(data, positions=x + 1, widths=0.6, showfliers=False, patch_artist=True)
        for idx, box in enumerate(bp["boxes"]):
            box.set_facecolor(group_colors[idx])
            box.set_alpha(0.75)
            box.set_edgecolor("black")
        ax.set_xticks(x + 1)
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
        ax.set_ylabel("Max contig length")
        ax.set_title(
            f"Max contig length for filtered bins (completeness>{min_completeness:g}, contamination<{max_contamination:g})"
        )
        ax.grid(axis="y", linestyle="--", alpha=0.4)
        ax.legend(
            handles=make_assembly_type_legend_handles(),
            title="Assembly type",
            loc="upper left",
            fontsize=8,
        )
        fig.tight_layout()
        out_path = os.path.join(outdir, f"filtered_{tag}_max_contig_length.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {out_path}")
        plt.close(fig)
    else:
        print("WARNING: max_contig_length column missing or empty; skipping max-contig-length plot")

    # 4) Mean and median genome size for filtered bins.
    if "genome_size" in filt.columns and filt["genome_size"].notna().any():
        fig, ax = plt.subplots(figsize=(max(12, len(summary) * 0.9), 6))
        width = 0.38
        for idx, row in summary.iterrows():
            base_color = group_colors[idx]
            ax.bar(
                x[idx] - width / 2,
                row["mean_genome_size"],
                width=width,
                color=base_color,
                alpha=0.6,
                edgecolor="white",
                linewidth=0.6,
                hatch="//",
                label="Mean genome size" if idx == 0 else "_nolegend_",
            )
            ax.bar(
                x[idx] + width / 2,
                row["median_genome_size"],
                width=width,
                color=base_color,
                edgecolor="white",
                linewidth=0.6,
                label="Median genome size" if idx == 0 else "_nolegend_",
            )
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
        ax.set_ylabel("Genome size")
        ax.set_title(
            f"Mean and median genome size for filtered bins (completeness>{min_completeness:g}, contamination<{max_contamination:g})"
        )
        ax.grid(axis="y", linestyle="--", alpha=0.4)
        stat_handles, stat_labels = ax.get_legend_handles_labels()
        ax.legend(
            handles=make_assembly_type_legend_handles() + make_workflow_legend_handles() + stat_handles,
            labels=[h.get_label() for h in make_assembly_type_legend_handles() + make_workflow_legend_handles()] + stat_labels,
            loc="upper left",
            fontsize=8,
            title="Assembly type / Workflow / Statistic",
        )
        fig.tight_layout()
        out_path = os.path.join(outdir, f"filtered_{tag}_genome_size_mean_median.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {out_path}")
        plt.close(fig)
    else:
        print("WARNING: genome_size column missing or empty; skipping genome-size mean/median plot")


# Plot medium-quality and high-quality filtered count comparisons.
def plot_stacked_filtered_bin_counts(df_bins: pd.DataFrame, outdir: str) -> None:
    if df_bins.empty:
        return

    medium = df_bins[
        (pd.to_numeric(df_bins["completeness"], errors="coerce") > 50.0)
        & (pd.to_numeric(df_bins["contamination"], errors="coerce") < 10.0)
    ].copy()
    high = df_bins[
        (pd.to_numeric(df_bins["completeness"], errors="coerce") > 90.0)
        & (pd.to_numeric(df_bins["contamination"], errors="coerce") < 5.0)
    ].copy()

    if medium.empty:
        print("WARNING: no bins passed medium filter; skipping stacked filtered count plot")
        return

    group_cols = ["workflow", "assembler", "sample"]
    medium_counts = medium.groupby(group_cols, dropna=False).size().reset_index(name="medium_count")
    high_counts = high.groupby(group_cols, dropna=False).size().reset_index(name="high_count")

    summary = medium_counts.merge(high_counts, on=group_cols, how="left")
    summary["high_count"] = summary["high_count"].fillna(0).astype(int)
    summary["medium_count"] = summary["medium_count"].astype(int)
    summary["medium_only_count"] = (summary["medium_count"] - summary["high_count"]).clip(lower=0).astype(int)

    sample_rank = {s: i for i, s in enumerate(SAMPLE_ORDER)}
    workflow_rank = {"metawrap": 0, "vamb": 1}
    summary["_type_rank"] = summary["assembler"].map(assembler_sort_key).map(lambda x: x[0])
    summary["_sample_rank"] = summary["sample"].map(lambda s: sample_rank.get(s, 99))
    summary["_workflow_rank"] = summary["workflow"].map(lambda w: workflow_rank.get(w, 99))
    summary = summary.sort_values(
        ["_type_rank", "assembler", "_sample_rank", "_workflow_rank"]
    ).reset_index(drop=True)

    out_summary = os.path.join(outdir, "filtered_stacked_gt50_lt10_gt90_lt5_bin_counts_summary.tsv")
    summary.drop(columns=["_type_rank", "_sample_rank", "_workflow_rank"]).to_csv(
        out_summary,
        sep="\t",
        index=False,
    )
    print(f"Saved: {out_summary}")

    labels = [
        f"{display_label(a)}\n({s})\n{w}"
        for w, a, s in summary[["workflow", "assembler", "sample"]].values
    ]
    x = np.arange(len(summary))
    assemblers = summary["assembler"].unique().tolist()
    color_map = assign_assembler_colors(assemblers)
    group_colors = [color_map.get(a, "#888888") for a in summary["assembler"]]
    workflow_hatch = {"metawrap": "", "vamb": "//"}

    fig, ax = plt.subplots(figsize=(max(12, len(summary) * 0.9), 6))
    for idx, row in summary.iterrows():
        base_color = group_colors[idx]
        hatch = workflow_hatch.get(row["workflow"], "")
        ax.bar(
            x[idx],
            row["medium_only_count"],
            color=base_color,
            alpha=0.45,
            edgecolor="white",
            linewidth=0.6,
            hatch=hatch,
        )
        ax.bar(
            x[idx],
            row["high_count"],
            bottom=row["medium_only_count"],
            color=base_color,
            edgecolor="white",
            linewidth=0.6,
            hatch=hatch,
        )

    from matplotlib.patches import Patch

    filter_handles = [
        Patch(facecolor="#bbbbbb", edgecolor="black", alpha=0.45, label=">50 completeness, <10 contamination only"),
        Patch(facecolor="#bbbbbb", edgecolor="black", label=">90 completeness, <5 contamination"),
    ]

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Bin count")
    ax.set_title("Filtered bin counts stacked: >50/<10 split by >90/<5 subset")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.legend(
        handles=make_assembly_type_legend_handles() + make_workflow_legend_handles() + filter_handles,
        title="Assembly type / Workflow / Filter",
        loc="upper left",
        fontsize=8,
    )

    fig.tight_layout()
    out_path = os.path.join(outdir, "filtered_stacked_gt50_lt10_gt90_lt5_bin_counts.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.close(fig)


# Plot filtered metrics for the >50/<10 threshold pair.
def plot_medium_filtered_bin_metrics(df_bins: pd.DataFrame, outdir: str) -> None:
    _plot_filtered_bin_metrics(
        df_bins=df_bins,
        outdir=outdir,
        min_completeness=50.0,
        max_contamination=10.0,
        tag="gt50_lt10",
    )


# Plot filtered metrics for the >90/<5 threshold pair.
def plot_high_filtered_bin_metrics(df_bins: pd.DataFrame, outdir: str) -> None:
    _plot_filtered_bin_metrics(
        df_bins=df_bins,
        outdir=outdir,
        min_completeness=90.0,
        max_contamination=5.0,
        tag="gt90_lt5",
    )


# Sort plot rows by assembly type, sample order, workflow, and set name.
def sort_plot_df(df: pd.DataFrame) -> pd.DataFrame:
    sample_rank = {s: i for i, s in enumerate(SAMPLE_ORDER)}
    workflow_rank = {"metawrap": 0, "vamb": 1}
    out = df.copy()
    out["_type_rank"] = out["assembler"].map(assembler_sort_key).map(lambda x: x[0])
    out["_sample_rank"] = out["sample"].map(lambda s: sample_rank.get(s, 99))
    out["_workflow_rank"] = out["workflow"].map(lambda w: workflow_rank.get(w, 99))
    out = out.sort_values(["_type_rank", "assembler", "_sample_rank", "_workflow_rank", "set_name"]).reset_index(drop=True)
    return out


# Plot total bin counts for every set.
def plot_refined_totals(df: pd.DataFrame, outdir: str) -> None:
    plot_df = sort_plot_df(df)
    assemblers = plot_df["assembler"].unique().tolist()
    color_map = assign_assembler_colors(assemblers)

    x = np.arange(len(plot_df))
    colors = [color_map.get(a, "#888888") for a in plot_df["assembler"]]

    workflow_hatch = {"metawrap": "", "vamb": "//"}

    fig, ax = plt.subplots(figsize=(max(12, len(plot_df) * 1.0), 6))
    for idx, row in plot_df.iterrows():
        ax.bar(
            x[idx],
            row["refined_total_bins"],
            color=colors[idx],
            edgecolor="white",
            linewidth=0.6,
            hatch=workflow_hatch.get(row["workflow"], ""),
            zorder=3,
        )

    tick_labels = [
        f"{display_label(a)}\n({s})\n{w}"
        for a, s, w in zip(plot_df["assembler"], plot_df["sample"], plot_df["workflow"])
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Number of bins")
    ax.set_title("MetaWRAP refined and VAMB bin totals by assembly and sample")
    ax.grid(axis="y", linestyle="--", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    from matplotlib.patches import Patch

    type_handles = [
        Patch(facecolor=DEFAULT_COLORS[t][1], label=t.replace("_", " "))
        for t in ASSEMBLY_TYPE_ORDER
        if t in DEFAULT_COLORS
    ]
    workflow_handles = [
        Patch(facecolor="#cccccc", edgecolor="black", hatch="", label="metawrap"),
        Patch(facecolor="#cccccc", edgecolor="black", hatch="//", label="vamb"),
    ]
    ax.legend(handles=type_handles + workflow_handles, title="Assembly type / Workflow", loc="upper left", fontsize=8)

    fig.tight_layout()
    out_path = os.path.join(outdir, "refined_total_bins.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.close(fig)


# Plot total bins with CheckM2 HQ/MQ markers overlaid.
def plot_stacked_quality_breakdown(df: pd.DataFrame, outdir: str) -> None:
    """
    Stacked bar chart showing quality breakdown:
    - HQ (solid, no hatch)
    - MQ (hatched //)
    - LQ (cross-hatched xx)
    - Unclassified remainder (dotted)
    
    Uses same color scheme as other plots per assembler.
    For entries without CheckM2 data, shows full bar as solid.
    """
    plot_df = sort_plot_df(df)
    assemblers = plot_df["assembler"].unique().tolist()
    color_map = assign_assembler_colors(assemblers)

    x = np.arange(len(plot_df))
    colors = [color_map.get(a, "#888888") for a in plot_df["assembler"]]

    fig, ax = plt.subplots(figsize=(max(14, len(plot_df) * 1.2), 6.5))

    # Build stacked segments for each row
    hq_vals = []
    mq_vals = []
    lq_vals = []
    unclass_vals = []

    for idx, row in plot_df.iterrows():
        if row["checkm2_available"] == 1:
            hq = row["hq_bins"] if not pd.isna(row["hq_bins"]) else 0
            mq = row["mq_bins"] if not pd.isna(row["mq_bins"]) else 0
            lq = row["lq_bins"] if not pd.isna(row["lq_bins"]) else 0
            total = row["refined_total_bins"]
            unclass = total - hq - mq - lq
            unclass = max(0, unclass)  # Ensure non-negative
        else:
            # No CheckM2 data: full bar is unclassified
            hq = 0
            mq = 0
            lq = 0
            unclass = row["refined_total_bins"]

        hq_vals.append(hq)
        mq_vals.append(mq)
        lq_vals.append(lq)
        unclass_vals.append(unclass)

    # Create stacked bars: LQ, MQ, HQ (bottom to top) with Unclassified on top
    ax.bar(x, lq_vals, color=colors, hatch="xx", label="LQ", edgecolor="white", linewidth=0.6, zorder=2)
    ax.bar(x, mq_vals, bottom=lq_vals, color=colors, hatch="//", label="MQ", edgecolor="white", linewidth=0.6, zorder=2)
    ax.bar(x, hq_vals, bottom=np.array(lq_vals) + np.array(mq_vals), color=colors, label="HQ", edgecolor="white", linewidth=0.6, zorder=2)

    unclass_bottom = np.array(lq_vals) + np.array(mq_vals) + np.array(hq_vals)
    ax.bar(x, unclass_vals, bottom=unclass_bottom, color="lightgray", hatch="...", label="Unclassified", edgecolor="white", linewidth=0.6, zorder=2)

    tick_labels = [
        f"{display_label(a)}\n({s})\n{w}"
        for a, s, w in zip(plot_df["assembler"], plot_df["sample"], plot_df["workflow"])
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Bin count")
    ax.set_title("Refined bins by CheckM2 quality tier (HQ/MQ/LQ/Unclassified)")
    ax.grid(axis="y", linestyle="--", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    from matplotlib.patches import Patch
    quality_handles = [
        Patch(facecolor="#cccccc", edgecolor="black", hatch="xx", label="LQ (Cross-hatched)"),
        Patch(facecolor="#cccccc", edgecolor="black", hatch="//", label="MQ (Hatched)"),
        Patch(facecolor="#cccccc", edgecolor="black", hatch="", label="HQ (Solid)"),
        Patch(facecolor="lightgray", edgecolor="black", hatch="...", label="Unclassified (Dotted)"),
    ]
    ax.legend(handles=quality_handles, title="Quality tier", loc="upper left", fontsize=9)

    fig.tight_layout()
    out_path = os.path.join(outdir, "refined_stacked_quality_breakdown.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.close(fig)


# Plot total bins with optional CheckM2 HQ/MQ scatter overlays.
def plot_with_checkm2_overlay(df: pd.DataFrame, outdir: str) -> None:
    plot_df = sort_plot_df(df)
    assemblers = plot_df["assembler"].unique().tolist()
    color_map = assign_assembler_colors(assemblers)

    x = np.arange(len(plot_df))
    colors = [color_map.get(a, "#888888") for a in plot_df["assembler"]]

    workflow_hatch = {"metawrap": "", "vamb": "//"}

    fig, ax = plt.subplots(figsize=(max(13, len(plot_df) * 1.1), 6))

    for idx, row in plot_df.iterrows():
        ax.bar(
            x[idx],
            row["refined_total_bins"],
            color=colors[idx],
            edgecolor="white",
            linewidth=0.6,
            alpha=0.75,
            hatch=workflow_hatch.get(row["workflow"], ""),
            zorder=2,
            label="_nolegend_",
        )

    # Only sets with CheckM2 reports get these markers.
    has_qc = plot_df["checkm2_available"] == 1
    ax.scatter(x[has_qc], plot_df.loc[has_qc, "hq_bins"], marker="o", s=70, color="#1b9e77", label="CheckM2 HQ", zorder=4)
    ax.scatter(x[has_qc], plot_df.loc[has_qc, "mq_bins"], marker="^", s=70, color="#d95f02", label="CheckM2 MQ", zorder=4)

    tick_labels = [
        f"{display_label(a)}\n({s})\n{w}"
        for a, s, w in zip(plot_df["assembler"], plot_df["sample"], plot_df["workflow"])
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Bin count")
    ax.set_title("MetaWRAP refined and VAMB totals with optional CheckM2 HQ/MQ overlay")
    ax.grid(axis="y", linestyle="--", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    from matplotlib.patches import Patch
    workflow_handles = [
        Patch(facecolor="#cccccc", edgecolor="black", hatch="", label="metawrap"),
        Patch(facecolor="#cccccc", edgecolor="black", hatch="//", label="vamb"),
    ]
    h, l = ax.get_legend_handles_labels()
    ax.legend(handles=workflow_handles + h, title="Workflow / CheckM2", loc="upper left", fontsize=9)

    fig.tight_layout()
    out_path = os.path.join(outdir, "refined_with_checkm2_overlay.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.close(fig)


# Compare VAMB and metaWRAP in separate panels for quick contrast.
def plot_triptych(df: pd.DataFrame, outdir: str) -> None:
    """
    Create a 2-panel comparison with independent scaling:
    - Left: VAMB bins (long-read)
    - Right: Hybrid + short-read assemblers (metawrap)
    Each panel has its own y-axis scale for within-category comparison.
    """
    # Categorize by assembly type
    def get_assembly_type(assembler: str) -> str:
        return ASSEMBLER_TO_ASSEMBLY_TYPE.get(assembler, "short_read")
    
    df_temp = df.copy()
    df_temp["assembly_type"] = df_temp["assembler"].apply(get_assembly_type)
    
    # Filter metawrap and vamb
    df_metawrap = df_temp[df_temp["workflow"] == "metawrap"].copy()
    df_vamb = df_temp[df_temp["workflow"] == "vamb"].copy()
    
    # Separate into panels
    short_hybrid = df_metawrap[df_metawrap["assembly_type"].isin(["short_read", "hybrid_read"])].copy()
    vamb = df_vamb.copy()
    
    # Sort each panel
    short_hybrid = sort_plot_df(short_hybrid) if not short_hybrid.empty else short_hybrid
    vamb = sort_plot_df(vamb) if not vamb.empty else vamb
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))  # Independent y-axes
    
    panels = [
        (vamb, axes[0], "VAMB (Long-read)"),
        (short_hybrid, axes[1], "Hybrid & Short-read (metaWRAP)"),
    ]
    
    for panel_df, ax, title in panels:
        if panel_df.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            ax.set_xticks([])
            continue
        
        assemblers = panel_df["assembler"].unique().tolist()
        color_map = assign_assembler_colors(assemblers)
        
        x = np.arange(len(panel_df))
        colors = [color_map.get(a, "#888888") for a in panel_df["assembler"]]
        
        # Build stacked bars with quality tiers
        hq_vals = []
        mq_vals = []
        lq_vals = []
        unclass_vals = []
        
        for _, row in panel_df.iterrows():
            if row["checkm2_available"] == 1:
                hq = row["hq_bins"] if not pd.isna(row["hq_bins"]) else 0
                mq = row["mq_bins"] if not pd.isna(row["mq_bins"]) else 0
                lq = row["lq_bins"] if not pd.isna(row["lq_bins"]) else 0
                total = row["refined_total_bins"]
                unclass = total - hq - mq - lq
                unclass = max(0, unclass)
            else:
                hq = 0
                mq = 0
                lq = 0
                unclass = row["refined_total_bins"]
            
            hq_vals.append(hq)
            mq_vals.append(mq)
            lq_vals.append(lq)
            unclass_vals.append(unclass)
        
        # Stacked bars: LQ, MQ, HQ (bottom to top) with Unclassified on top
        ax.bar(x, lq_vals, color=colors, hatch="xx", label="LQ", edgecolor="white", linewidth=0.5, zorder=2)
        ax.bar(x, mq_vals, bottom=lq_vals, color=colors, hatch="//", label="MQ", edgecolor="white", linewidth=0.5, zorder=2)
        ax.bar(x, hq_vals, bottom=np.array(lq_vals) + np.array(mq_vals), color=colors, label="HQ", edgecolor="white", linewidth=0.5, zorder=2)
        
        unclass_bottom = np.array(lq_vals) + np.array(mq_vals) + np.array(hq_vals)
        ax.bar(x, unclass_vals, bottom=unclass_bottom, color="lightgray", hatch="...", label="Unclassified", edgecolor="white", linewidth=0.5, zorder=2)
        
        tick_labels = [f"{display_label(a)}\n({s})" for a, s in zip(panel_df["assembler"], panel_df["sample"])]
        ax.set_xticks(x)
        ax.set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=8)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_ylabel("Bin count", fontsize=10)
        ax.grid(axis="y", linestyle="--", alpha=0.5, zorder=0)
        ax.set_axisbelow(True)
    
    from matplotlib.patches import Patch
    quality_handles = [
        Patch(facecolor="#cccccc", edgecolor="black", hatch="xx", label="LQ (Cross-hatched)"),
        Patch(facecolor="#cccccc", edgecolor="black", hatch="//", label="MQ (Hatched)"),
        Patch(facecolor="#cccccc", edgecolor="black", hatch="", label="HQ (Solid)"),
        Patch(facecolor="lightgray", edgecolor="black", hatch="...", label="Unclassified (Dotted)"),
    ]
    fig.legend(handles=quality_handles, title="Quality tier", loc="upper center", ncol=4, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    out_path = os.path.join(outdir, "triptych_comparison.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.close(fig)


# Parse config defaults, build summary tables, and write all outputs.
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--config",
        default=DEFAULT_CONFIG_PATH,
        help="Path to plotting_pipeline_config.sh (sources configs_master.conf)",
    )
    parser.add_argument(
        "--refined_root",
        default=None,
        help="Path to refined bins root",
    )
    parser.add_argument(
        "--vamb_root",
        default=None,
        help="Path to VAMB output root",
    )
    parser.add_argument(
        "--checkm2_root",
        default=None,
        help="Path to CheckM2 output root for metawrap sets",
    )
    parser.add_argument(
        "--checkm2_vamb_root",
        default=None,
        help="Path to complete CheckM2 output root for vamb sets",
    )
    parser.add_argument(
        "--bins_subdir",
        default="fast_mode_bins",
        help="Subdirectory inside each refined set used for counting bin FASTAs (default: fast_mode_bins)",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Directory for output plots and table",
    )
    parser.add_argument(
        "--extensions",
        default=".fa,.fasta,.fna",
        help="Comma-separated FASTA extensions to count in bins_subdir",
    )
    parser.add_argument(
        "--no_vamb",
        action="store_true",
        help="Disable inclusion of VAMB outputs",
    )
    parser.add_argument(
        "--filter_workflow",
        default="",
        help="Optional comma-separated workflows to keep in bin-level aggregation (e.g., metawrap,vamb)",
    )
    parser.add_argument(
        "--filter_assembler",
        default="",
        help="Optional comma-separated assemblers to keep in bin-level aggregation",
    )
    parser.add_argument(
        "--filter_sample",
        default="",
        help="Optional comma-separated samples to keep in bin-level aggregation (e.g., S1,S2,S5)",
    )
    parser.add_argument(
        "--filter_binner",
        default="",
        help="Optional comma-separated binners to keep in bin-level aggregation (e.g., metawrap,vamb)",
    )
    args = parser.parse_args()

    config = load_shell_config(args.config)
    sample_order = config.get("PROJECT_SAMPLE_NAMES", [])
    if isinstance(sample_order, list) and sample_order:
        global SAMPLE_ORDER
        SAMPLE_ORDER = [str(sample) for sample in sample_order]

    plotting_root = str(config.get("PROJECT_PLOTTING_ROOT") or os.path.join(os.path.dirname(__file__), "..", "assembly_stats"))
    args.outdir = _first_nonempty(
        args.outdir,
        config.get("PROJECT_REFINED_BIN_COUNTS_OUTDIR") if isinstance(config.get("PROJECT_REFINED_BIN_COUNTS_OUTDIR"), str) else None,
        os.path.join(plotting_root, "plots", "refined_bin_counts"),
    )
    args.refined_root = _first_nonempty(
        args.refined_root,
        config.get("PROJECT_REFINED_BINS_ROOT") if isinstance(config.get("PROJECT_REFINED_BINS_ROOT"), str) else None,
        os.path.join(str(config.get("PROJECT_DATA_DIR") or ""), "binning_outputs", "refined_bins"),
    )
    args.vamb_root = _first_nonempty(
        args.vamb_root,
        config.get("PROJECT_VAMB_OUTPUT_ROOT") if isinstance(config.get("PROJECT_VAMB_OUTPUT_ROOT"), str) else None,
        os.path.join(str(config.get("PROJECT_DATA_DIR") or ""), "binning_outputs", "long_mapped_binning"),
    )
    args.checkm2_root = _first_nonempty(
        args.checkm2_root,
        config.get("PROJECT_CHECKM2_REFINED_ROOT") if isinstance(config.get("PROJECT_CHECKM2_REFINED_ROOT"), str) else None,
        os.path.join(str(config.get("PROJECT_DATA_DIR") or ""), "binning_outputs", "checkm2_refined_bins"),
    )
    args.checkm2_vamb_root = _first_nonempty(
        args.checkm2_vamb_root,
        config.get("PROJECT_CHECKM2_VAMB_ROOT") if isinstance(config.get("PROJECT_CHECKM2_VAMB_ROOT"), str) else None,
        os.path.join(str(config.get("PROJECT_DATA_DIR") or ""), "binning_outputs", "checkm2_vamb_bins_full"),
    )

    if not os.path.isdir(args.refined_root):
        raise SystemExit(f"ERROR: refined root not found: {args.refined_root}")

    os.makedirs(args.outdir, exist_ok=True)

    extensions = [x.strip() for x in args.extensions.split(",") if x.strip()]
    df = collect_all_bins_and_checkm2(
        refined_root=args.refined_root,
        vamb_root=args.vamb_root,
        checkm2_metawrap_root=args.checkm2_root,
        checkm2_vamb_root=args.checkm2_vamb_root,
        bins_subdir=args.bins_subdir,
        extensions=extensions,
        include_vamb=(not args.no_vamb),
    )

    if df.empty:
        raise SystemExit("ERROR: no metawrap or vamb set directories found")

    df = sort_plot_df(df)

    out_tsv = os.path.join(args.outdir, "metawrap_vamb_checkm2_counts.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Saved: {out_tsv}")

    # Build per-bin aggregated CheckM2 table with origin metadata.
    bin_df = collect_checkm2_bins_aggregated(
        refined_root=args.refined_root,
        vamb_root=args.vamb_root,
        checkm2_metawrap_root=args.checkm2_root,
        checkm2_vamb_root=args.checkm2_vamb_root,
        include_vamb=(not args.no_vamb),
    )

    fw = [x.strip() for x in args.filter_workflow.split(",") if x.strip()]
    fa = [x.strip() for x in args.filter_assembler.split(",") if x.strip()]
    fs = [x.strip() for x in args.filter_sample.split(",") if x.strip()]
    fb = [x.strip() for x in args.filter_binner.split(",") if x.strip()]

    if not bin_df.empty:
        bin_df = apply_bin_filters(bin_df, workflows=fw, assemblers=fa, samples=fs, binners=fb)

    out_bins_tsv = os.path.join(args.outdir, "checkm2_bins_aggregated.tsv")
    if bin_df.empty:
        print("WARNING: bin-level CheckM2 aggregation is empty after collection/filtering")
    else:
        bin_df.to_csv(out_bins_tsv, sep="\t", index=False)
        print(f"Saved: {out_bins_tsv}")

    print("\nSummary preview:")
    preview_cols = [
        "set_name",
        "workflow",
        "assembler",
        "sample",
        "refined_subdir",
        "count_method",
        "refined_total_bins",
        "checkm2_total_bins",
        "hq_bins",
        "mq_bins",
        "checkm2_available",
    ]
    print(df[preview_cols].to_string(index=False))

    plot_refined_totals(df, args.outdir)
    plot_with_checkm2_overlay(df, args.outdir)
    plot_stacked_quality_breakdown(df, args.outdir)
    plot_triptych(df, args.outdir)
    plot_average_completeness_contamination(bin_df, args.outdir)
    plot_quality_tier_fractions(bin_df, args.outdir)
    plot_stacked_filtered_bin_counts(bin_df, args.outdir)
    plot_medium_filtered_bin_metrics(bin_df, args.outdir)
    plot_high_filtered_bin_metrics(bin_df, args.outdir)


if __name__ == "__main__":
    main()
