#!/usr/bin/env python3
"""Contig counts split by threshold: 4 graphs, one per threshold."""

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


def plot_contig_counts_by_threshold(
    threshold_df: pd.DataFrame,
    outpath: str = "contig_counts_by_threshold.png",
    theme: str = "default",
    samples: list = None,
    thresholds_kbp: list = None,
) -> None:
    """Create 4 subplots: one per threshold, showing S1/S2/S5 samples as grouped bars.
    
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
    
    # Build color map by assembly type
    assembler_colors = {}
    for asm in assemblers:
        asm_type = classify_assembly_type(asm)
        if asm_type == "short_read":
            color = ASSEMBLY_TYPE_COLORS["short_read"][0]  # Green
        elif asm_type == "long_read":
            color = ASSEMBLY_TYPE_COLORS["long_read"][0]   # Blue
        elif asm_type == "hybrid_read":
            color = ASSEMBLY_TYPE_COLORS["hybrid_read"][0]  # Purple
        else:
            color = "#333333"
        assembler_colors[asm] = color
    
    # Extract fontsize values
    fontsize = config["fontsize"]
    title_fontsize = fontsize + 2
    label_fontsize = fontsize
    tick_fontsize = fontsize - 1
    
    # Create figure with 2x2 subplots (one per threshold)
    fig, axes = plt.subplots(
        2, 2,
        figsize=(14, 10),
        dpi=config["dpi"],
    )
    axes = axes.flatten()
    
    fig.suptitle("Contig Counts by Threshold", 
                 fontsize=title_fontsize + 1, 
                 fontweight="bold")
    
    # Sample colors and markers
    sample_colors = {"S1": "#1f77b4", "S2": "#ff7f0e", "S5": "#2ca02c"}
    
    for threshold_idx, threshold_kbp in enumerate(thresholds_kbp):
        ax = axes[threshold_idx]
        threshold_label = f"{threshold_kbp//1000}M" if threshold_kbp >= 1000 else f"{threshold_kbp}K"
        threshold_data = plot_df[plot_df["threshold_kbp"] == threshold_kbp]
        
        if threshold_data.empty:
            ax.text(0.5, 0.5, f"No data for {threshold_label}", 
                   ha="center", va="center", transform=ax.transAxes)
            continue
        
        # Prepare data: assembler x sample
        x_positions = np.arange(len(assemblers))
        bar_width = 0.25
        
        for sample_idx, sample in enumerate(samples):
            counts = []
            
            for asm in assemblers:
                row = threshold_data[(threshold_data["assembler"] == asm) & 
                                    (threshold_data["sample"] == sample)]
                if len(row) > 0:
                    counts.append(row.iloc[0]["num_contigs"])
                else:
                    counts.append(0)
            
            offset = bar_width * (sample_idx - len(samples)/2 + 0.5)
            ax.bar(
                x_positions + offset,
                counts,
                bar_width,
                label=sample,
                color=sample_colors[sample],
                alpha=0.8,
                edgecolor="black",
                linewidth=config["linewidth"] / 2,
            )
        
        # Format axes
        ax.set_xlabel("Assembler", fontsize=label_fontsize)
        ax.set_ylabel("Contig Count", fontsize=label_fontsize)
        ax.set_title(f"Threshold: ≥ {threshold_label} bp", 
                    fontsize=title_fontsize, fontweight="bold")
        
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
        
        # Add legend only on first subplot
        if threshold_idx == 0:
            ax.legend(
                title="Sample",
                loc="upper left",
                fontsize=tick_fontsize,
                title_fontsize=label_fontsize,
            )
        
        # Add assembly type color legend to right side on last plot
        if threshold_idx == 3:
            legend_elements = []
            for asm_type in ASSEMBLY_TYPE_ORDER:
                if asm_type == "short_read":
                    color = ASSEMBLY_TYPE_COLORS["short_read"][0]
                    label = "Short-read"
                elif asm_type == "long_read":
                    color = ASSEMBLY_TYPE_COLORS["long_read"][0]
                    label = "Long-read"
                elif asm_type == "hybrid_read":
                    color = ASSEMBLY_TYPE_COLORS["hybrid_read"][0]
                    label = "Hybrid-read"
                else:
                    continue
                
                legend_elements.append(
                    plt.Rectangle((0, 0), 1, 1, fc=color, alpha=0.8, 
                                 edgecolor="black", linewidth=0.5)
                )
            
            ax.legend(
                handles=legend_elements,
                labels=["Short-read", "Long-read", "Hybrid-read"],
                title="Assembly Type",
                loc="upper right",
                fontsize=tick_fontsize,
                title_fontsize=label_fontsize,
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
    outpath = os.path.join(args.outdir, f"contig_counts_by_threshold_{args.theme}.png")
    plot_contig_counts_by_threshold(threshold_df, outpath, theme=args.theme)
