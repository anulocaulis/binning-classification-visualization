#!/usr/bin/env python3
"""
poster_figures.py

Poster-focused plots for samples S1, S2, and S5:
- Rarefaction curves per sample (raw cumulative contig length)
- N50 vs L50 scatter plot
- NG50 vs L50 scatter plot
- Yield vs complexity scatter plot

All plots use a consistent assembler color mapping grouped by assembly type:
long_read, short_read, hybrid_read.
"""

# imports
import argparse
import base64
import gzip
import os
import re
import sys
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# Ensure the script can find local modules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_summary_stats import (  # noqa: E402
    DEFAULT_NG50_TARGET_BP,
    add_recalculated_ng50,
    classify_assembly_type,
    parse_summary_stats_log,
    parse_threshold_stats_log,
)

# Set default orders and color palettes for assembly types and samples
ASSEMBLY_TYPE_ORDER = ["long_read", "short_read", "hybrid_read"]
SAMPLE_ORDER = ["S1", "S2", "S5"]
DEFAULT_THRESHOLD_SERIES_KBP = [1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000]

DEFAULT_COLORS = {
    "short_read": ["#A1D99B", "#74C476","#238B45"],
    "long_read": ["#C6DBEF", "#41B6C4", "#08519C"],
    "hybrid_read": ["#D0D1E6", "#6A51A3"],
}

ASSEMBLER_LABEL_ALIASES = {
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

ASSEMBLER_COLOR_OVERRIDES = {}

# Helper functions for assembler sorting, color validation, and parsing
def _assembler_sort_key(assembler: str) -> Tuple[int, str]:
    assembly_type = classify_assembly_type(assembler)
    type_rank = (
        ASSEMBLY_TYPE_ORDER.index(assembly_type)
        if assembly_type in ASSEMBLY_TYPE_ORDER
        else len(ASSEMBLY_TYPE_ORDER)
    )
    return (type_rank, str(assembler))


def _valid_hex(color: str) -> bool:
    return bool(re.fullmatch(r"#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})", color.strip()))


def _parse_hex_list(raw: str, label: str) -> List[str]:
    colors = [c.strip() for c in raw.split(",") if c.strip()]
    if not colors:
        raise ValueError(f"{label} palette is empty.")
    invalid = [c for c in colors if not _valid_hex(c)]
    if invalid:
        raise ValueError(
            f"Invalid hex color(s) for {label}: {', '.join(invalid)}. Use #RGB or #RRGGBB."
        )
    return colors


def _assign_assembler_colors(
    assemblers: Sequence[str], color_map_by_type: Dict[str, Sequence[str]]
) -> Dict[str, str]:
    assembler_color: Dict[str, str] = {}
    for assembly_type in ASSEMBLY_TYPE_ORDER:
        members = sorted([a for a in assemblers if classify_assembly_type(a) == assembly_type])
        if not members:
            continue
        palette = list(color_map_by_type[assembly_type])
        for idx, assembler in enumerate(members):
            assembler_color[assembler] = palette[idx % len(palette)]
    assembler_color.update(ASSEMBLER_COLOR_OVERRIDES)
    return assembler_color


def _display_assembler_name(assembler: str) -> str:
    label = ASSEMBLER_LABEL_ALIASES.get(assembler, str(assembler))
    assembly_type = classify_assembly_type(assembler)
    type_label = {
        "short_read": "short",
        "long_read": "long",
        "hybrid_read": "hybrid",
    }.get(assembly_type, str(assembly_type))
    return f"{label} ({type_label})"

# Helper functions for sample sorting, file handling, and rarefaction curve preparation
def _sample_key(sample: str):
    if isinstance(sample, str) and sample.startswith("S") and sample[1:].isdigit():
        return int(sample[1:])
    return str(sample)


def _open_text_maybe_gzip(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

# Helper function to extract contig lengths from a FASTA file, with optional minimum length filtering
def _contig_lengths_from_fasta(path: str, min_length_bp: int = 1000) -> List[int]:
    lengths: List[int] = []
    current = 0

    with _open_text_maybe_gzip(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current >= min_length_bp:
                    lengths.append(current)
                current = 0
                continue
            current += len(line)

    if current >= min_length_bp:
        lengths.append(current)
    return lengths

# Helper function to downsample rarefaction curves to a manageable number of points
def _downsample_curve(
    contig_count: Sequence[int],
    cumulative_mbp: Sequence[float],
    max_points_per_curve: int,
) -> Tuple[List[int], List[float]]:
    if max_points_per_curve <= 0 or len(contig_count) <= max_points_per_curve:
        return list(contig_count), list(cumulative_mbp)

    n = len(contig_count)
    step = max(1, (n - 1) // (max_points_per_curve - 1))
    idxs = list(range(0, n, step))
    if idxs[-1] != n - 1:
        idxs.append(n - 1)

    return [contig_count[i] for i in idxs], [cumulative_mbp[i] for i in idxs]

# Helper function to prepare rarefaction data from threshold tables
def _prepare_rarefaction(threshold_df: pd.DataFrame) -> pd.DataFrame:
    df = threshold_df.drop_duplicates(
        subset=["sample", "assembler", "threshold_label", "contig_length_bp", "num_contigs"],
        keep="first",
    ).copy()
    if "assembly_type" not in df.columns:
        df["assembly_type"] = df["assembler"].map(classify_assembly_type)
    df["assembly_type"] = df["assembly_type"].fillna("short_read")

    # Keep only poster samples.
    df = df[df["sample"].isin(SAMPLE_ORDER)].copy()
    df = df.dropna(subset=["assembler", "num_contigs", "contig_length_bp", "threshold_kbp"])
    df = df[(df["threshold_kbp"] >= 1.0) | (df["threshold_label"] == "All")].copy()

    # Prepare raw cumulative length in Mb.
    df["contig_length_mb"] = df["contig_length_bp"] / 1_000_000.0

    # Match plot_rarefaction.py ordering: longest contigs first.
    df = df.sort_values(
        ["assembly_type", "assembler", "sample", "threshold_kbp"],
        ascending=[True, True, True, False],
    ).reset_index(drop=True)
    return df

#### Dataframe setup
# Helper function to check if a DataFrame looks like a threshold table
def _looks_like_threshold_table(df: pd.DataFrame) -> bool:
    required = {"sample", "assembler", "threshold_label", "threshold_kbp", "contig_length_bp", "num_contigs"}
    return required.issubset(set(df.columns))


# Helper function to prepare summary data for plotting
def _prepare_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    df = summary_df.drop_duplicates(
        subset=["sample", "assembler", "sequence_total_mb", "scaffold_total", "n50_kbp", "l50"],
        keep="first",
    ).copy()

    if "assembly_type" not in df.columns:
        df["assembly_type"] = df["assembler"].map(classify_assembly_type)
    df["assembly_type"] = df["assembly_type"].fillna("short_read")

    df = df[df["sample"].isin(SAMPLE_ORDER)].copy()
    return df

##### Rarefaction
# Helper function to prepare true rarefaction curves from summary data
def _prepare_true_rarefaction_from_summary(
    summ_df: pd.DataFrame,
    max_points_per_curve: int,
    min_length_bp: int = 1000,
) -> pd.DataFrame:
    records: List[Dict[str, object]] = []

    source_rows = summ_df.dropna(subset=["sample", "assembler", "assembly_path"]).drop_duplicates(
        subset=["sample", "assembler"], keep="first"
    )
    for _, row in source_rows.iterrows():
        assembly_path = str(row["assembly_path"])
        if not os.path.isfile(assembly_path):
            print(f"Skipping rarefaction FASTA (missing): {assembly_path}")
            continue

        lengths = sorted(_contig_lengths_from_fasta(assembly_path, min_length_bp=min_length_bp), reverse=True)
        if not lengths:
            continue

        cumulative_mbp: List[float] = []
        running_bp = 0
        for length in lengths:
            running_bp += int(length)
            cumulative_mbp.append(running_bp / 1_000_000.0)

        sampled_count, sampled_cumulative = _downsample_curve(
            contig_count=list(range(1, len(lengths) + 1)),
            cumulative_mbp=cumulative_mbp,
            max_points_per_curve=max_points_per_curve,
        )

        for contig_count, contig_length_mb in zip(sampled_count, sampled_cumulative):
            records.append(
                {
                    "sample": str(row["sample"]),
                    "assembler": str(row["assembler"]),
                    "num_contigs": int(contig_count),
                    "contig_length_mb": float(contig_length_mb),
                }
            )

    if not records:
        return pd.DataFrame(columns=["sample", "assembler", "num_contigs", "contig_length_mb"])

    out = pd.DataFrame(records)
    out["assembly_type"] = out["assembler"].map(classify_assembly_type).fillna("short_read")
    out = out.sort_values(["sample", "assembly_type", "assembler", "num_contigs"]).reset_index(drop=True)
    return out

# Helper function to draw a rarefaction curve for a single sample
def _draw_rarefaction_sample(
    sample_df: pd.DataFrame,
    sample: str,
    assembler_order: Sequence[str],
    assembler_color: Dict[str, str],
    outpath: str,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))

    for assembler in assembler_order:
        asm_df = sample_df[sample_df["assembler"] == assembler].dropna(
            subset=["num_contigs", "contig_length_mb"]
        )
        if asm_df.empty:
            continue

        x_vals = pd.concat(
            [pd.Series([0.0]), asm_df["num_contigs"].astype(float)], ignore_index=True
        )
        y_vals = pd.concat(
            [pd.Series([0.0]), asm_df["contig_length_mb"].astype(float)], ignore_index=True
        )

        color = assembler_color[assembler]
        label = _display_assembler_name(assembler)
        ax.fill_between(x_vals, y_vals, alpha=0.12, color=color)
        ax.plot(x_vals, y_vals, color=color, linewidth=2, label=label)
        ax.scatter(asm_df["num_contigs"], asm_df["contig_length_mb"], color=color, s=2, zorder=4)

    ax.set_title(f"{sample} Rarefaction", fontsize=13)
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f"{int(v):,}" if v >= 1 else f"{v:.2g}")
    )
    ax.set_xlabel("Cumulative contig count (longest to shortest, log scale)")
    ax.set_ylabel("Cumulative Contig Length (Mb, raw)")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f} Mb"))
    ax.legend(title="Assembler", bbox_to_anchor=(1.01, 1), loc="upper left", borderaxespad=0)

    fig.tight_layout()
    fig.savefig(outpath, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")

##### Rarefaction plotting
# Helper function to draw side-by-side rarefaction curves for multiple samples
def _draw_rarefaction_triptych(
    rare_df: pd.DataFrame,
    samples: Sequence[str],
    assembler_order: Sequence[str],
    assembler_color: Dict[str, str],
    outpath: str,
) -> None:
    """Side-by-side rarefaction panels for the given samples, sharing a single legend."""
    n = len(samples)
    fig, axes = plt.subplots(1, n, figsize=(9 * n, 5.5), sharey=True)
    if n == 1:
        axes = [axes]

    legend_handles: List = []
    legend_labels: List[str] = []

    for ax, sample in zip(axes, samples):
        sample_df = rare_df[rare_df["sample"] == sample].copy()
        for assembler in assembler_order:
            asm_df = sample_df[sample_df["assembler"] == assembler].dropna(
                subset=["num_contigs", "contig_length_mb"]
            )
            if asm_df.empty:
                continue

            x_vals = pd.concat(
                [pd.Series([0.0]), asm_df["num_contigs"].astype(float)], ignore_index=True
            )
            y_vals = pd.concat(
                [pd.Series([0.0]), asm_df["contig_length_mb"].astype(float)], ignore_index=True
            )

            color = assembler_color[assembler]
            display = _display_assembler_name(assembler)
            line, = ax.plot(x_vals, y_vals, color=color, linewidth=2, label=display)
            ax.fill_between(x_vals, y_vals, alpha=0.12, color=color)
            ax.scatter(asm_df["num_contigs"], asm_df["contig_length_mb"], color=color, s=6, zorder=4)

            if display not in legend_labels:
                legend_handles.append(line)
                legend_labels.append(display)

        ax.set_title(f"{sample} Rarefaction", fontsize=13)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, _: f"{int(v):,}" if v >= 1 else f"{v:.2g}")
        )
        ax.set_xlabel("Cumulative contig count (longest to shortest, log scale)")
        if ax is axes[0]:
            ax.set_ylabel("Cumulative Contig Length (Mb, raw)")
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f} Mb"))

    fig.legend(
        legend_handles,
        legend_labels,
        title="Assembler",
        bbox_to_anchor=(1.01, 0.5),
        loc="center left",
        borderaxespad=0,
    )
    fig.tight_layout()
    fig.savefig(outpath, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


# Helper function to format threshold labels for display
def _format_threshold_label(threshold_kbp: float) -> str:
    if threshold_kbp == 0:
        return "all"
    if threshold_kbp < 1:
        return f"{int(threshold_kbp * 1000)}bp"
    if threshold_kbp < 1000:
        return f"{int(threshold_kbp)}kbp" if threshold_kbp == int(threshold_kbp) else f"{threshold_kbp}kbp"
    mbp = threshold_kbp / 1000.0
    return f"{int(mbp)}mbp" if mbp == int(mbp) else f"{mbp}mbp"


# Helper function to create a display label for a threshold value
def _threshold_display_label(threshold_kbp: float) -> str:
    if threshold_kbp < 1000:
        return f"{int(threshold_kbp)} Kbp" if threshold_kbp == int(threshold_kbp) else f"{threshold_kbp:g} Kbp"
    mbp = threshold_kbp / 1000.0
    return f"{int(mbp)} Mbp" if mbp == int(mbp) else f"{mbp:g} Mbp"

# Helper function to build threshold ticks and labels for plotting
def _build_threshold_ticks(values: Sequence[float]):
    ordered = sorted({float(v) for v in values if pd.notna(v)})
    labels = []
    for v in ordered:
        if v == 0:
            labels.append("All")
        elif v < 1:
            labels.append(f"{int(v * 1000)} bp")
        elif v < 1000:
            labels.append(f"{v:g} Kbp")
        else:
            labels.append(f"{(v / 1000.0):g} Mbp")
    return ordered, labels

##### Threshold curve plotting
# Helper function to draw threshold curves for a single sample
def _draw_threshold_curves_for_sample(
    sample_df: pd.DataFrame,
    sample: str,
    assembler_order: Sequence[str],
    assembler_color: Dict[str, str],
    outdir: str,
) -> None:
    threshold_vals, threshold_labels = _build_threshold_ticks(sample_df["threshold_kbp"].tolist())
    if not threshold_vals:
        return

    # Normalize each assembler to its "All" row for percent-retained plots.
    base = (
        sample_df[sample_df["threshold_kbp"] == 0.0]
        .drop_duplicates(subset=["assembler"], keep="first")[["assembler", "num_contigs", "contig_length_mb"]]
        .rename(columns={"num_contigs": "base_num_contigs", "contig_length_mb": "base_contig_length_mb"})
    )
    curve_df = sample_df.merge(base, on="assembler", how="left")
    curve_df["contig_count_pct"] = (curve_df["num_contigs"] / curve_df["base_num_contigs"]) * 100.0
    curve_df["contig_length_pct"] = (curve_df["contig_length_mb"] / curve_df["base_contig_length_mb"]) * 100.0

    plot_specs = [
        ("contig_count_pct", "Contigs Retained (%)", f"{sample}_threshold_contig_count_pct.png", (0, 105), mticker.FuncFormatter(lambda v, _: f"{v:.0f}%")),
        ("contig_length_pct", "Contig Length Retained (%)", f"{sample}_threshold_contig_length_pct.png", (0, 105), mticker.FuncFormatter(lambda v, _: f"{v:.0f}%")),
        ("num_contigs", "Contig Count", f"{sample}_threshold_contig_count_abs.png", None, mticker.FuncFormatter(lambda v, _: f"{int(v):,}")),
        ("contig_length_mb", "Contig Length (Mb)", f"{sample}_threshold_contig_length_abs.png", None, mticker.FuncFormatter(lambda v, _: f"{v:,.0f} Mb")),
    ]

    for y_col, y_label, filename, ylim, yformatter in plot_specs:
        fig, ax = plt.subplots(figsize=(9, 5.5))
        for assembler in assembler_order:
            asm_df = curve_df[curve_df["assembler"] == assembler].dropna(subset=["threshold_kbp", y_col]).copy()
            if asm_df.empty:
                continue
            asm_df = asm_df.sort_values("threshold_kbp")
            ax.plot(
                asm_df["threshold_kbp"],
                asm_df[y_col],
                marker="o",
                markersize=4,
                linewidth=2,
                color=assembler_color[assembler],
                label=_display_assembler_name(assembler),
            )

        ax.set_title(f"{sample} Threshold Curve")
        ax.set_xlabel("Minimum Contig Length")
        ax.set_ylabel(y_label)
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.yaxis.set_major_formatter(yformatter)
        ax.set_xticks(threshold_vals)
        ax.set_xticklabels(threshold_labels, rotation=45, ha="right")
        ax.legend(title="Assembler", bbox_to_anchor=(1.01, 1), loc="upper left", borderaxespad=0)
        fig.tight_layout()
        outpath = os.path.join(outdir, filename)
        fig.savefig(outpath, dpi=600, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {outpath}")


##### Threshold bar chart plotting

# Helper function to draw a threshold bar chart for a specific threshold
def _draw_threshold_bar_chart(
    df: pd.DataFrame,
    threshold_kbp: float,
    assembler_order: Sequence[str],
    assembler_color: Dict[str, str],
    outpath: str,
) -> None:
    subset = df[df["threshold_kbp"] == threshold_kbp].copy()
    if subset.empty:
        return

    summary = (
        subset.groupby("assembler", as_index=False)
        .agg(
            mean_num_contigs=("num_contigs", "mean"),
        )
    )
    # Ensure every assembler in assembler_order appears, even those with 0 contigs at this threshold.
    all_assemblers_df = pd.DataFrame({"assembler": list(assembler_order)})
    summary = all_assemblers_df.merge(summary, on="assembler", how="left")
    summary["mean_num_contigs"] = summary["mean_num_contigs"].fillna(0.0)
    summary["assembler"] = pd.Categorical(summary["assembler"], categories=list(assembler_order), ordered=True)
    summary = summary.sort_values("assembler")

    sample_order = [sample for sample in ["S1", "S2", "S5"] if sample in set(subset["sample"])]
    sample_order.extend(sorted(set(subset["sample"]) - set(sample_order)))
    sample_markers = {
        "S1": "o",
        "S2": "s",
        "S5": "^",
    }
    if len(sample_order) == 1:
        sample_offsets = {sample_order[0]: 0.0}
    else:
        offset_positions = np.linspace(-0.18, 0.18, len(sample_order))
        sample_offsets = dict(zip(sample_order, offset_positions))

    fig, ax = plt.subplots(figsize=(10, 6))
    x_positions = range(len(summary))
    colors = [assembler_color[str(assembler)] for assembler in summary["assembler"]]
    labels = [_display_assembler_name(str(assembler)) for assembler in summary["assembler"]]

    ax.bar(
        list(x_positions),
        summary["mean_num_contigs"],
        color=colors,
        edgecolor="white",
        linewidth=0.7,
    )

    position_lookup = {str(assembler): idx for idx, assembler in enumerate(summary["assembler"])}
    for sample in sample_order:
        sample_subset = subset[subset["sample"] == sample].copy()

        # Fill in 0 for assemblers missing at this threshold for this sample.
        present = set(sample_subset["assembler"].tolist())
        missing_rows = [
            {"assembler": a, "num_contigs": 0, "sample": sample}
            for a in assembler_order if a not in present
        ]
        if missing_rows:
            sample_subset = pd.concat(
                [sample_subset, pd.DataFrame(missing_rows)], ignore_index=True
            )

        sample_subset["assembler"] = pd.Categorical(
            sample_subset["assembler"],
            categories=list(summary["assembler"]),
            ordered=True,
        )
        sample_subset = sample_subset.sort_values("assembler")

        for _, row in sample_subset.iterrows():
            assembler = str(row["assembler"])
            if assembler not in position_lookup:
                continue

            ax.scatter(
                position_lookup[assembler] + sample_offsets[sample],
                row["num_contigs"],
                color=assembler_color[assembler],
                marker=sample_markers.get(sample, "D"),
                edgecolors="black",
                linewidths=0.35,
                s=34,
                zorder=3,
            )

    ax.set_title(f"Contigs >= {_threshold_display_label(threshold_kbp)}")
    ax.set_xlabel("Assembler")
    ax.set_ylabel("Contig Count")
    ax.set_xticks(list(x_positions))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.grid(True, axis="y", alpha=0.25)

    sample_legend_handles = [
        Line2D(
            [0],
            [0],
            marker=sample_markers.get(sample, "D"),
            linestyle="None",
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=0.9,
            markersize=7,
            label=sample,
        )
        for sample in sample_order
    ]
    ax.legend(handles=sample_legend_handles, title="Sample", loc="upper right", frameon=False)

    fig.tight_layout()
    fig.savefig(outpath, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")

##### N50 vs L50 scatter plotting
def _draw_n50_l50_scatter(
    df: pd.DataFrame,
    assembler_order: Sequence[str],
    assembler_color: Dict[str, str],
    outpath: str,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))

    plot_df = df.dropna(subset=["assembler", "sample", "n50_kbp", "l50"]).copy()

    for assembler in assembler_order:
        sub = plot_df[plot_df["assembler"] == assembler]
        if sub.empty:
            continue

        ax.scatter(
            sub["l50"],
            sub["n50_kbp"],
            label=_display_assembler_name(assembler),
            color=assembler_color[assembler],
            s=70,
            alpha=0.85,
            edgecolors="white",
            linewidths=0.4,
        )

        for _, row in sub.iterrows():
            ax.annotate(
                row["sample"],
                xy=(row["l50"], row["n50_kbp"]),
                xytext=(4, 3),
                textcoords="offset points",
                fontsize=6.5,
                color=assembler_color[assembler],
                alpha=0.9,
            )

    ax.set_title("N50 vs L50")
    ax.set_xlabel("L50 (number of contigs needed to reach 50% of assembly)")
    ax.set_ylabel("N50 (Kbp)")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.legend(title="Assembler", bbox_to_anchor=(1.01, 1), loc="upper left")

    fig.tight_layout()
    fig.savefig(outpath, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


##### NG50 vs L50 scatter plotting

def _draw_ng50_l50_scatter(
    df: pd.DataFrame,
    assembler_order: Sequence[str],
    assembler_color: Dict[str, str],
    target_genome_bp: int,
    outpath: str,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))

    plot_df = df.dropna(subset=["assembler", "sample", "ng50_kbp", "l50"]).copy()

    for assembler in assembler_order:
        sub = plot_df[plot_df["assembler"] == assembler]
        if sub.empty:
            continue

        ax.scatter(
            sub["l50"],
            sub["ng50_kbp"],
            label=_display_assembler_name(assembler),
            color=assembler_color[assembler],
            s=70,
            alpha=0.85,
            edgecolors="white",
            linewidths=0.4,
        )

        for _, row in sub.iterrows():
            ax.annotate(
                row["sample"],
                xy=(row["l50"], row["ng50_kbp"]),
                xytext=(4, 3),
                textcoords="offset points",
                fontsize=6.5,
                color=assembler_color[assembler],
                alpha=0.9,
            )

    ax.set_title(f"NG50 vs L50 (Target: {target_genome_bp:,} bp)")
    ax.set_xlabel("L50 (number of contigs needed to reach 50% of assembly)")
    ax.set_ylabel("NG50 (Kbp)")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.legend(title="Assembler", bbox_to_anchor=(1.01, 1), loc="upper left")

    fig.tight_layout()
    fig.savefig(outpath, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


##### Yield vs Complexity scatter plotting
# Helper function to draw a yield vs complexity scatter plot
def _draw_yield_complexity(
    df: pd.DataFrame,
    assembler_order: Sequence[str],
    assembler_color: Dict[str, str],
    outpath: str,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 7))

    plot_df = df.dropna(subset=["assembler", "sample", "sequence_total_mb", "scaffold_total"]).copy()

    for assembler in assembler_order:
        sub = plot_df[plot_df["assembler"] == assembler]
        if sub.empty:
            continue

        ax.scatter(
            sub["sequence_total_mb"],
            sub["scaffold_total"],
            label=_display_assembler_name(assembler),
            color=assembler_color[assembler],
            s=80,
            alpha=0.85,
            edgecolors="white",
            linewidths=0.5,
            zorder=3,
        )

        for _, row in sub.iterrows():
            ax.annotate(
                row["sample"],
                xy=(row["sequence_total_mb"], row["scaffold_total"]),
                xytext=(4, 3),
                textcoords="offset points",
                fontsize=6.5,
                color=assembler_color[assembler],
                alpha=0.9,
            )

    ax.set_title("Assembly Yield vs Contig Count by Assembler")
    ax.set_xlabel("Total Contig Sequence (Mb)")
    ax.set_ylabel("Total Contig Count")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.legend(title="Assembler", bbox_to_anchor=(1.01, 1), loc="upper left")

    fig.tight_layout()
    fig.savefig(outpath, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")

##### Poster report HTML generation
# Helper function to build an HTML report for the poster figures
def _build_poster_report_html(
    outdir: str,
    rare_df: pd.DataFrame,
    summ_df: pd.DataFrame,
    threshold_list_kbp: Sequence[float],
) -> str:
    rare_tbl = rare_df.copy()
    rare_tbl["assembler_label"] = rare_tbl["assembler"].map(_display_assembler_name)
    rare_tbl = rare_tbl[
        [
            "sample",
            "assembler_label",
            "threshold_kbp",
            "num_contigs",
            "contig_length_mb",
        ]
    ].rename(
        columns={
            "sample": "Sample",
            "assembler_label": "Assembler",
            "threshold_kbp": "Threshold (Kbp)",
            "num_contigs": "Contigs >= Threshold",
            "contig_length_mb": "Cumulative Length (Mb)",
        }
    )
    rare_tbl = rare_tbl.sort_values(["Sample", "Assembler", "Threshold (Kbp)"])

    summ_tbl = summ_df.copy()
    summ_tbl["assembler_label"] = summ_tbl["assembler"].map(_display_assembler_name)
    summary_columns = [
        "sample",
        "assembler_label",
        "sequence_total_mb",
        "scaffold_total",
        "n50_kbp",
        "ng50_kbp",
        "l50",
        "lg50",
    ]
    summary_columns = [col for col in summary_columns if col in summ_tbl.columns]
    summ_tbl = summ_tbl[summary_columns].rename(
        columns={
            "sample": "Sample",
            "assembler_label": "Assembler",
            "sequence_total_mb": "Total Sequence (Mb)",
            "scaffold_total": "Total Contigs",
            "n50_kbp": "N50 (Kbp)",
            "ng50_kbp": "NG50 (Kbp)",
            "l50": "L50",
            "lg50": "LG50",
        }
    )
    summ_tbl = summ_tbl.sort_values(["Sample", "Assembler"])

    threshold_summary = (
        rare_df[rare_df["threshold_kbp"].isin(list(threshold_list_kbp))]
        .groupby(["threshold_kbp", "assembler"], as_index=False)
        .agg(mean_contigs=("num_contigs", "mean"), sd_contigs=("num_contigs", "std"), n=("num_contigs", "count"))
    )
    threshold_summary["sd_contigs"] = threshold_summary["sd_contigs"].fillna(0.0)
    threshold_summary["assembler"] = threshold_summary["assembler"].map(_display_assembler_name)
    threshold_summary = threshold_summary.rename(
        columns={
            "threshold_kbp": "Threshold (Kbp)",
            "assembler": "Assembler",
            "mean_contigs": "Mean Contigs",
            "sd_contigs": "SD Contigs",
            "n": "N Samples",
        }
    ).sort_values(["Threshold (Kbp)", "Assembler"])

    figure_files = sorted(
        [f for f in os.listdir(outdir) if f.lower().endswith(".png")],
        key=lambda x: (
            0 if x.endswith("_rarefaction_mb.png") else 1,
            x,
        ),
    )
    figure_blocks_list = []
    for fname in figure_files:
        image_path = os.path.join(outdir, fname)
        with open(image_path, "rb") as image_handle:
            image_b64 = base64.b64encode(image_handle.read()).decode("ascii")
        data_uri = f"data:image/png;base64,{image_b64}"
        figure_blocks_list.append(
            f"<div class='card'><h3>{fname}</h3>"
            f"<a href='{data_uri}'><img src='{data_uri}' alt='{fname}' loading='lazy'></a></div>"
        )
    figure_blocks = "\n".join(figure_blocks_list)

    css = """
body { font-family: Arial, sans-serif; margin: 24px; color: #1a1a1a; }
h1, h2 { margin: 0 0 10px 0; }
h2 { margin-top: 26px; }
.muted { color: #555; margin-bottom: 14px; }
.grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(360px, 1fr)); gap: 16px; }
.card { border: 1px solid #d9d9d9; border-radius: 10px; padding: 12px; background: #fbfbfb; }
.card h3 { margin: 0 0 10px 0; font-size: 14px; font-weight: 600; }
.card img { width: 100%; height: auto; border-radius: 6px; border: 1px solid #ececec; }
table { border-collapse: collapse; width: 100%; font-size: 12px; margin-top: 8px; }
th, td { border: 1px solid #dcdcdc; padding: 6px 8px; text-align: right; }
th:first-child, td:first-child, th:nth-child(2), td:nth-child(2) { text-align: left; }
thead { background: #f1f1f1; position: sticky; top: 0; }
.table-wrap { max-height: 420px; overflow: auto; border: 1px solid #e2e2e2; border-radius: 8px; }
"""

    html = f"""<!doctype html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1'>
  <title>Poster Figures Report</title>
  <style>{css}</style>
</head>
<body>
    <h1>Poster Figure Report</h1>
  <p class='muted'>Generated from poster_figures.py outputs. Click any image to open full resolution.</p>

    <h2>Figures</h2>
  <div class='grid'>
    {figure_blocks}
  </div>

    <h2>Summary Statistics</h2>
  <div class='table-wrap'>
    {summ_tbl.to_html(index=False, float_format=lambda x: f"{x:,.2f}")}
  </div>

    <h2>Rarefaction Data</h2>
  <div class='table-wrap'>
    {rare_tbl.to_html(index=False, float_format=lambda x: f"{x:,.2f}")}
  </div>

    <h2>Threshold Contig Summary</h2>
  <div class='table-wrap'>
    {threshold_summary.to_html(index=False, float_format=lambda x: f"{x:,.2f}")}
  </div>
</body>
</html>
"""

    report_path = os.path.join(outdir, "report.html")
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write(html)
    print(f"Saved: {report_path}")
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Poster-ready S1/S2/S5 figures: per-sample rarefaction (raw Mb), "
            "N50-vs-L50 scatter, and yield-vs-complexity scatter."
        )
    )
    parser.add_argument("--log", default="ge1000_summary_stats.txt")
    parser.add_argument(
        "--input-csv",
        default="filtered_assembly/data/rarefaction_data.csv",
        help="Optional precomputed rarefaction CSV from plot_rarefaction.py.",
    )
    parser.add_argument(
        "--summary-csv",
        default="filtered_assembly/data/parsed_summary_stats.csv",
        help="Optional precomputed summary CSV for N50/L50 and yield/complexity.",
    )
    parser.add_argument(
        "--ng50-target-bp",
        type=int,
        default=DEFAULT_NG50_TARGET_BP,
        help="Expected genome size in bp for NG50 (default: 5,000,000).",
    )
    parser.add_argument(
        "--data-dir",
        default="filtered_assembly/data",
        help="Directory to write normalized/parsed CSV outputs.",
    )
    parser.add_argument("--outdir", default="plots/ge1000_poster")
    parser.add_argument(
        "--long-colors",
        default=",".join(DEFAULT_COLORS["long_read"]),
        help="Comma-separated hex colors for long_read assemblers.",
    )
    parser.add_argument(
        "--short-colors",
        default=",".join(DEFAULT_COLORS["short_read"]),
        help="Comma-separated hex colors for short_read assemblers.",
    )
    parser.add_argument(
        "--hybrid-colors",
        default=",".join(DEFAULT_COLORS["hybrid_read"]),
        help="Comma-separated hex colors for hybrid_read assemblers.",
    )
    parser.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=DEFAULT_THRESHOLD_SERIES_KBP,
        help=(
            "Thresholds in Kbp for contig-count bar charts "
            "(default: 1 2.5 5 10 25 50 100 250 500 1000 2500 5000)."
        ),
    )
    parser.add_argument(
        "--max-tsv-points-per-curve",
        type=int,
        default=5000,
        help="Maximum points written per true rarefaction curve series (0 writes all points).",
    )
    args = parser.parse_args()

    try:
        colors_by_type = {
            "long_read": _parse_hex_list(args.long_colors, "long_read"),
            "short_read": _parse_hex_list(args.short_colors, "short_read"),
            "hybrid_read": _parse_hex_list(args.hybrid_colors, "hybrid_read"),
        }
    except ValueError as exc:
        raise SystemExit(str(exc))

    if os.path.exists(args.input_csv):
        threshold_df = pd.read_csv(args.input_csv)
        if not _looks_like_threshold_table(threshold_df):
            print(
                f"Input CSV does not look like a threshold table ({args.input_csv}); "
                f"falling back to parsing threshold rows from log: {args.log}"
            )
            threshold_df = parse_threshold_stats_log(args.log)
    else:
        threshold_df = parse_threshold_stats_log(args.log)
    if threshold_df.empty:
        raise SystemExit("No threshold rows found. Check --input-csv or --log path.")

    if os.path.exists(args.summary_csv):
        summary_df = pd.read_csv(args.summary_csv)
    else:
        summary_df = parse_summary_stats_log(args.log)
    summary_df = add_recalculated_ng50(summary_df, args.ng50_target_bp)
    if summary_df.empty:
        raise SystemExit("No summary rows found. Check --summary-csv or --log path.")

    os.makedirs(args.data_dir, exist_ok=True)
    threshold_df.to_csv(os.path.join(args.data_dir, "parsed_threshold_stats.csv"), index=False)
    summary_df.to_csv(os.path.join(args.data_dir, "parsed_summary_stats.csv"), index=False)

    rare_df = _prepare_rarefaction(threshold_df)
    summ_df = _prepare_summary(summary_df)
    rare_true_df = _prepare_true_rarefaction_from_summary(
        summ_df,
        max_points_per_curve=args.max_tsv_points_per_curve,
        min_length_bp=1000,
    )

    if rare_df.empty:
        raise SystemExit("No rarefaction rows available after filtering to S1/S2/S5.")
    if summ_df.empty:
        raise SystemExit("No summary rows available after filtering to S1/S2/S5.")
    if rare_true_df.empty:
        raise SystemExit("No FASTA-derived rarefaction rows available after filtering to S1/S2/S5.")

    all_assemblers = sorted(
        set(rare_true_df["assembler"].dropna().unique().tolist())
        | set(summ_df["assembler"].dropna().unique().tolist()),
        key=_assembler_sort_key,
    )
    assembler_color = _assign_assembler_colors(all_assemblers, colors_by_type)

    os.makedirs(args.outdir, exist_ok=True)

    legacy_outputs = [
        "rarefaction_triptych_mb.png",
        "s1_s2_s5_n50.png",
        "s1_s2_s5_l50.png",
        "s1_s2_s5_rarefaction_raw_single_panel.png",
        "s1_s2_s5_rarefaction_raw_data.csv",
        "s1_s2_s5_fragmentation_data.csv",
    ]
    for filename in legacy_outputs:
        path = os.path.join(args.outdir, filename)
        if os.path.exists(path):
            os.remove(path)

    rare_df.to_csv(os.path.join(args.outdir, "s1_s2_s5_rarefaction_data.csv"), index=False)
    rare_true_df.to_csv(os.path.join(args.outdir, "s1_s2_s5_rarefaction_true_data.csv"), index=False)
    summ_df.to_csv(os.path.join(args.outdir, "s1_s2_s5_summary_data.csv"), index=False)
    rare_df.to_csv(os.path.join(args.data_dir, "rarefaction_data.csv"), index=False)
    rare_true_df.to_csv(os.path.join(args.data_dir, "rarefaction_true_data.csv"), index=False)
    rare_true_df.to_csv(os.path.join(args.data_dir, "s1_s2_s5_true_cumulative_contig_lengths.tsv"), sep="\t", index=False)
    summ_df.to_csv(os.path.join(args.data_dir, "yield_complexity_data.csv"), index=False)
    rare_df.to_csv(os.path.join(args.data_dir, "threshold_curves_data.csv"), index=False)

    for sample in sorted(rare_true_df["sample"].unique(), key=_sample_key):
        sample_df = rare_true_df[rare_true_df["sample"] == sample].copy()
        _draw_rarefaction_sample(
            sample_df=sample_df,
            sample=sample,
            assembler_order=all_assemblers,
            assembler_color=assembler_color,
            outpath=os.path.join(args.outdir, f"{sample}_rarefaction_mb.png"),
        )
    for sample in sorted(rare_df["sample"].unique(), key=_sample_key):
        sample_df = rare_df[rare_df["sample"] == sample].copy()
        _draw_threshold_curves_for_sample(
            sample_df=sample_df,
            sample=sample,
            assembler_order=all_assemblers,
            assembler_color=assembler_color,
            outdir=args.outdir,
        )

    triptych_samples = [s for s in ["S1", "S2"] if s in rare_true_df["sample"].values]
    if len(triptych_samples) >= 2:
        _draw_rarefaction_triptych(
            rare_df=rare_true_df,
            samples=triptych_samples,
            assembler_order=all_assemblers,
            assembler_color=assembler_color,
            outpath=os.path.join(args.outdir, "S1_S2_rarefaction_triptych.png"),
        )

    _draw_n50_l50_scatter(
        df=summ_df,
        assembler_order=all_assemblers,
        assembler_color=assembler_color,
        outpath=os.path.join(args.outdir, "n50_l50_scatter.png"),
    )

    _draw_ng50_l50_scatter(
        df=summ_df,
        assembler_order=all_assemblers,
        assembler_color=assembler_color,
        target_genome_bp=args.ng50_target_bp,
        outpath=os.path.join(args.outdir, "ng50_l50_scatter.png"),
    )

    _draw_yield_complexity(
        df=summ_df,
        assembler_order=all_assemblers,
        assembler_color=assembler_color,
        outpath=os.path.join(args.outdir, "yield_vs_complexity.png"),
    )

    for threshold_kbp in sorted(set(args.thresholds)):
        _draw_threshold_bar_chart(
            df=rare_df,
            threshold_kbp=threshold_kbp,
            assembler_order=all_assemblers,
            assembler_color=assembler_color,
            outpath=os.path.join(
                args.outdir,
                f"{_format_threshold_label(threshold_kbp)}_contig_count_bar.png",
            ),
        )

    _build_poster_report_html(
        outdir=args.outdir,
        rare_df=rare_df,
        summ_df=summ_df,
        threshold_list_kbp=sorted(set(args.thresholds)),
    )


if __name__ == "__main__":
    main()
