#!/usr/bin/env python3
"""N50 vs L50 scatter plots with custom genome size targets (LG50 metric)."""

import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import numpy as np

from viz_config import THEMES, ASSEMBLY_TYPE_COLORS, ASSEMBLER_LABELS, ASSEMBLY_TYPE_ORDER
from data_loader import SAMPLE_ORDER

try:
    from plot_summary_stats import (
        classify_assembly_type,
        parse_summary_stats_log,
        add_recalculated_ng50,
    )
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from plot_summary_stats import (
        classify_assembly_type,
        parse_summary_stats_log,
        add_recalculated_ng50,
    )


def assign_colors(assemblers: list, color_map: dict) -> dict:
    """Assign colors to assemblers by type."""
    assembler_color = {}
    for assembly_type in ASSEMBLY_TYPE_ORDER:
        members = sorted([a for a in assemblers if classify_assembly_type(a) == assembly_type])
        if not members:
            continue
        palette = list(color_map.get(assembly_type, []))
        for idx, assembler in enumerate(members):
            assembler_color[assembler] = palette[idx % len(palette)]
    return assembler_color


def plot_ng50_l50_custom(
    summary_df: pd.DataFrame,
    outpath: str = "ng50_l50_custom.png",
    theme: str = "default",
    target_genome_bp: int = 35_000_000,
    color_map: dict = None,
    samples: list = None,
) -> None:
    """Plot Ng50 vs L50 for custom target genome size (LG50 metric).
    
    Args:
        summary_df: DataFrame with columns [sample, assembler, n50_kbp, l50, ...]
        outpath: Output PNG path
        theme: Theme name from THEMES
        target_genome_bp: Target genome size in bp for Ng50 calculation
        color_map: Dict mapping assembly_type -> list of hex colors
        samples: List of samples to plot (default: S1, S2, S5)
    """
    if color_map is None:
        color_map = ASSEMBLY_TYPE_COLORS
    if samples is None:
        samples = ["S1", "S2", "S5"]
    
    config = THEMES[theme]
    plt.style.use(config["style"])
    
    # Add recalculated Ng50 for target genome size
    plot_df = add_recalculated_ng50(summary_df.copy(), target_genome_bp=target_genome_bp)
    
    # Filter to target samples
    plot_df = plot_df[plot_df["sample"].isin(samples)]
    
    if plot_df.empty or "ng50_kbp" not in plot_df.columns:
        print(f"ERROR: Could not calculate Ng50 for target {target_genome_bp} bp")
        return
    
    # Extract fontsize values
    fontsize = config["fontsize"]
    title_fontsize = fontsize + 2
    label_fontsize = fontsize
    tick_fontsize = fontsize - 1
    dpi = config["dpi"]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7), dpi=dpi)
    
    # Target genome size in Mbp for title
    target_mbp = target_genome_bp / 1_000_000
    fig.suptitle(f"Ng50 vs L50 (Target Genome: {target_mbp:.0f} Mbp)", 
                 fontsize=title_fontsize, fontweight="bold")
    
    # Assign colors
    assemblers = sorted(plot_df["assembler"].unique())
    color_map_assembled = assign_colors(assemblers, color_map)
    
    # Sample markers
    sample_markers = {"S1": "o", "S2": "s", "S5": "^"}
    sample_marker_size = 150
    
    # Plot points
    for assembler in assemblers:
        asm_df = plot_df[plot_df["assembler"] == assembler]
        
        for sample in samples:
            sample_asm_df = asm_df[asm_df["sample"] == sample]
            if sample_asm_df.empty:
                continue
            
            # Get Ng50 and L50 values
            ng50_kbp = sample_asm_df.iloc[0].get("ng50_kbp", None)
            l50 = sample_asm_df.iloc[0].get("l50", None)
            
            if ng50_kbp is None or l50 is None:
                continue
            
            ax.scatter(
                ng50_kbp,
                l50,
                s=sample_marker_size,
                c=color_map_assembled.get(assembler, "#333333"),
                marker=sample_markers.get(sample, "o"),
                alpha=0.7,
                edgecolors="black",
                linewidth=1,
                label=f"{ASSEMBLER_LABELS.get(assembler, assembler)} ({sample})" if False else "",
            )
    
    # Formatting
    ax.set_xlabel("Ng50 (kbp)", fontsize=label_fontsize)
    ax.set_ylabel("L50 (count)", fontsize=label_fontsize)
    ax.tick_params(axis="both", labelsize=tick_fontsize)
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.set_axisbelow(True)
    
    # Add legend for assemblers and samples
    legend_elements = []
    for assembler in assemblers:
        legend_elements.append(
            plt.Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=color_map_assembled.get(assembler, "#333333"),
                      markersize=8, label=ASSEMBLER_LABELS.get(assembler, assembler),
                      markeredgecolor="black", markeredgewidth=0.5)
        )
    
    legend_elements.append(plt.Line2D([0], [0], color="w", label=""))  # Spacer
    
    for sample, marker in sample_markers.items():
        legend_elements.append(
            plt.Line2D([0], [0], marker=marker, color="w", markerfacecolor="gray",
                      markersize=8, label=sample, markeredgecolor="black", markeredgewidth=0.5)
        )
    
    ax.legend(handles=legend_elements, loc="best", fontsize=tick_fontsize, 
             title="Assembler / Sample", title_fontsize=label_fontsize)
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches="tight")
    print(f"✓ Saved: {outpath}")
    plt.close()


if __name__ == "__main__":
    import argparse
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from plot_summary_stats import parse_summary_stats_log
    
    parser = argparse.ArgumentParser()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_log = os.path.join(script_dir, "..", "ge1000_summary_stats.txt")
    parser.add_argument("--log-path", default=default_log)
    parser.add_argument("--outdir", default="sandbox")
    parser.add_argument("--theme", default="default", choices=list(THEMES.keys()))
    parser.add_argument("--targets", nargs="+", type=int, default=[35_000_000, 280_000_000, 70_000_000],
                       help="Target genome sizes in bp")
    args = parser.parse_args()
    
    # Load data
    summary_df = parse_summary_stats_log(args.log_path)
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    # Generate plots for each target
    for target_bp in args.targets:
        target_mbp = target_bp / 1_000_000
        outpath = os.path.join(args.outdir, f"ng50_l50_lg50_{target_mbp:.0f}Mbp_{args.theme}.png")
        plot_ng50_l50_custom(summary_df, outpath, theme=args.theme, 
                            target_genome_bp=target_bp)
