#!/usr/bin/env python3
"""Four triptychs of contig counts: one PNG per threshold (10K/100K/1M/2.5M).
Each PNG shows S1/S2/S5 panels with bars colored by assembly type."""

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
    from plot_summary_stats import classify_assembly_type, parse_threshold_stats_log
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from plot_summary_stats import classify_assembly_type, parse_threshold_stats_log


def assign_colors(assemblers: list, color_map: dict) -> dict:
    """Assign colors to assemblers by type, matching rarefaction plot style."""
    assembler_color = {}
    for assembly_type in ASSEMBLY_TYPE_ORDER:
        members = sorted([a for a in assemblers if classify_assembly_type(a) == assembly_type])
        if not members:
            continue
        palette = list(color_map.get(assembly_type, []))
        for idx, assembler in enumerate(members):
            assembler_color[assembler] = palette[idx % len(palette)]
    return assembler_color


def plot_contig_counts_threshold_triptychs(
    threshold_df: pd.DataFrame,
    outdir: str = ".",
    theme: str = "default",
    samples: list = None,
    thresholds_kbp: list = None,
) -> None:
    """Generate 4 separate triptych PNGs: one per threshold.
    Each PNG has 3 panels (S1, S2, S5) with bars colored by assembly type.
    
    Args:
        threshold_df: DataFrame with columns [sample, assembler, threshold_label, 
                      threshold_kbp, num_contigs, ...]
        outdir: Output directory
        theme: Theme name from THEMES
        samples: List of samples to plot (default: S1, S2, S5)
        thresholds_kbp: List of threshold values in Kbp (default: [10, 100, 1000, 2500])
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
    
    # Filter data to samples and thresholds we'll use
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
    
    # Assign colors per assembler (matching rarefaction plot style)
    assembler_colors = assign_colors(assemblers, ASSEMBLY_TYPE_COLORS)
    
    # Create one triptych per threshold
    for threshold_kbp in thresholds_kbp:
        threshold_label_title = f"{threshold_kbp//1000}M bp" if threshold_kbp >= 1000 else f"{threshold_kbp}K bp"
        threshold_label_file = f"{threshold_kbp//1000}M" if threshold_kbp >= 1000 else f"{threshold_kbp}K"
        
        fig, axes = plt.subplots(
            1, 3,
            figsize=(18, 5),
            dpi=config["dpi"],
        )
        fig.suptitle(f"Contig Counts (≥ {threshold_label_title})", 
                     fontsize=title_fontsize, 
                     fontweight="bold")
        
        # Plot one panel per sample
        for idx, sample in enumerate(samples):
            ax = axes[idx]
            sample_data = plot_df[
                (plot_df["sample"] == sample) & 
                (plot_df["threshold_kbp"] == threshold_kbp)
            ]
            
            if sample_data.empty:
                ax.text(0.5, 0.5, f"No data for {sample}", 
                       ha="center", va="center", transform=ax.transAxes)
                continue
            
            # Prepare data for plotting
            x_positions = np.arange(len(assemblers))
            
            for asm_idx, asm in enumerate(assemblers):
                row = sample_data[sample_data["assembler"] == asm]
                if len(row) > 0:
                    count = row.iloc[0]["num_contigs"]
                else:
                    count = 0
                
                # Use pre-assigned color per assembler (matches rarefaction plot)
                color = assembler_colors[asm]
                
                ax.bar(
                    asm_idx,
                    count,
                    width=0.7,
                    color=color,
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
        
        # Add legend for assembler colors (positioned to right of all panels)
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=assembler_colors[asm], edgecolor="black", 
                  label=ASSEMBLER_LABELS.get(asm, asm), alpha=0.8)
            for asm in assemblers
        ]
        fig.legend(
            handles=legend_elements,
            title="Assembler",
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            ncol=1,
            frameon=True,
            fontsize=tick_fontsize,
        )
        
        plt.tight_layout(rect=[0, 0, 0.92, 1])
        
        # Save PNG
        outpath = os.path.join(outdir, f"contig_counts_triptych_{threshold_label_file}_{theme}.png")
        plt.savefig(outpath, dpi=config["dpi"], bbox_inches="tight")
        print(f"✓ Saved: {outpath}")
        plt.close()


if __name__ == "__main__":
    import argparse
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from plot_summary_stats import parse_threshold_stats_log
    
    parser = argparse.ArgumentParser(description=__doc__)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_log = os.path.join(script_dir, "..", "ge1000_summary_stats.txt")
    parser.add_argument("--log-path", default=default_log, help="Path to ge1000_summary_stats.txt")
    parser.add_argument("--outdir", default="plots/sandbox", help="Output directory")
    parser.add_argument("--theme", default="default", choices=list(THEMES.keys()), help="Theme name")
    args = parser.parse_args()
    
    # Load data
    threshold_df = parse_threshold_stats_log(args.log_path)
    
    os.makedirs(args.outdir, exist_ok=True)
    plot_contig_counts_threshold_triptychs(threshold_df=threshold_df, outdir=args.outdir, theme=args.theme)
