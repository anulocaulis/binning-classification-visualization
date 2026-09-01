#!/usr/bin/env python3
"""Contig count triptych: S1/S2/S5 panels with assembler bars by threshold."""

import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

from viz_config import THEMES, ASSEMBLY_TYPE_COLORS, ASSEMBLER_LABELS, ASSEMBLY_TYPE_ORDER
from data_loader import SAMPLE_ORDER

try:
    from plot_summary_stats import classify_assembly_type
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from plot_summary_stats import classify_assembly_type


def plot_contig_counts_triptych(
    threshold_df: pd.DataFrame,
    outpath: str = "contig_counts_triptych.png",
    theme: str = "default",
    samples: list = None,
    thresholds_kbp: list = None,
) -> None:
    """Create triptych of contig counts (S1, S2, S5) by assembler and threshold.
    
    Args:
        threshold_df: DataFrame with columns [sample, assembler, threshold_label, 
                      threshold_kbp, num_contigs, ...]
        outpath: Output PNG path
        theme: Theme name from THEMES
        samples: List of samples to plot (default: S1, S2, S5)
        thresholds_kbp: List of threshold values in Kbp to plot (default: [10, 100, 1000, 2500])
    """
    if samples is None:
        samples = ["S1", "S2", "S5"]
    if thresholds_kbp is None:
        thresholds_kbp = [10, 100, 1000, 2500]
    
    config = THEMES[theme]
    plt.style.use(config["style"])
    
    # Extract fontsize values
    fontsize = config["fontsize"]
    title_fontsize = fontsize + 2
    label_fontsize = fontsize
    tick_fontsize = fontsize - 1
    edge_width = config["linewidth"]
    
    # Filter data
    plot_df = threshold_df[threshold_df["sample"].isin(samples)].copy()
    plot_df = plot_df[plot_df["threshold_kbp"].isin(thresholds_kbp)]
    
    if plot_df.empty:
        print("ERROR: No data after filtering")
        return
    
    # Get unique assemblers and sort by assembly type
    assemblers = sorted(plot_df["assembler"].unique(), key=lambda a: (
        ASSEMBLY_TYPE_ORDER.index(classify_assembly_type(a) or "short_read"),
        a
    ))
    
    # Build color map: one color per threshold
    threshold_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    threshold_color_map = dict(zip(thresholds_kbp, threshold_colors))
    
    # Create figure with 3 subplots (one per sample)
    fig, axes = plt.subplots(
        1, 3,
        figsize=(16, 5),
        dpi=config["dpi"],
        sharey=True,
    )
    fig.suptitle("Contig Counts by Threshold (≥ threshold bp)", 
                 fontsize=title_fontsize, 
                 fontweight="bold")
    
    for idx, sample in enumerate(samples):
        ax = axes[idx]
        sample_df = plot_df[plot_df["sample"] == sample]
        
        if sample_df.empty:
            ax.text(0.5, 0.5, f"No data for {sample}", 
                   ha="center", va="center", transform=ax.transAxes)
            continue
        
        # Prepare data: assembler x threshold
        x_positions = np.arange(len(assemblers))
        bar_width = 0.2
        
        for t_idx, threshold_kbp in enumerate(thresholds_kbp):
            threshold_label = f"{threshold_kbp//1000}M" if threshold_kbp >= 1000 else f"{threshold_kbp}K"
            counts = []
            
            for asm in assemblers:
                row = sample_df[(sample_df["assembler"] == asm) & 
                               (sample_df["threshold_kbp"] == threshold_kbp)]
                if len(row) > 0:
                    counts.append(row.iloc[0]["num_contigs"])
                else:
                    counts.append(0)
            
            offset = bar_width * (t_idx - len(thresholds_kbp)/2 + 0.5)
            ax.bar(
                x_positions + offset,
                counts,
                bar_width,
                label=threshold_label,
                color=threshold_color_map[threshold_kbp],
                alpha=0.8,
                edgecolor="black",
                linewidth=edge_width,
            )
        
        # Format axes
        ax.set_xlabel("Assembler", fontsize=label_fontsize)
        if idx == 0:
            ax.set_ylabel("Contig Count", fontsize=label_fontsize)
        
        ax.set_title(sample, fontsize=title_fontsize, fontweight="bold")
        ax.set_xticks(x_positions)
        ax.set_xticklabels(
            [ASSEMBLER_LABELS.get(a, a) for a in assemblers],
            rotation=45,
            ha="right",
            fontsize=tick_fontsize,
        )
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f"{int(x/1000)}K" if x >= 1000 else f"{int(x)}"))
        ax.tick_params(axis="y", labelsize=tick_fontsize)
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        ax.set_axisbelow(True)
    
    # Legend
    axes[1].legend(
        title="Threshold",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=4,
        frameon=True,
        fontsize=tick_fontsize,
    )
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=config["dpi"], bbox_inches="tight")
    print(f"✓ Saved: {outpath}")
    plt.close()


if __name__ == "__main__":
    import argparse
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from plot_summary_stats import parse_threshold_stats_log
    
    parser = argparse.ArgumentParser()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_log = os.path.join(script_dir, "..", "ge1000_summary_stats.txt")
    parser.add_argument("--log-path", default=default_log)
    parser.add_argument("--outdir", default="sandbox")
    parser.add_argument("--theme", default="default", choices=list(THEMES.keys()))
    args = parser.parse_args()
    
    # Load data
    threshold_df = parse_threshold_stats_log(args.log_path)
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    # Generate plot
    outpath = os.path.join(args.outdir, f"contig_counts_triptych_{args.theme}.png")
    plot_contig_counts_triptych(threshold_df, outpath, theme=args.theme)
