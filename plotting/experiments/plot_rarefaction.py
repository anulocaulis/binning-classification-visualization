"""Rarefaction curve plotting - cumulative contig length by count."""
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

from viz_config import THEMES, ASSEMBLY_TYPE_COLORS, ASSEMBLER_LABELS, ASSEMBLY_TYPE_ORDER
from data_loader import SAMPLE_ORDER

try:
    from plot_summary_stats import classify_assembly_type
except ImportError:
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from plot_summary_stats import classify_assembly_type


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


def plot_rarefaction(
    rare_df: pd.DataFrame,
    samples: list = None,
    outpath: str = "rarefaction.png",
    theme: str = "default",
    color_map: dict = None,
) -> None:
    """Plot rarefaction curves for specified samples.
    
    Args:
        rare_df: DataFrame with columns [sample, assembler, num_contigs, contig_length_mb]
        samples: List of samples to plot (default: S1, S2, S5)
        outpath: Output PNG path
        theme: Theme name from THEMES
        color_map: Dict mapping assembly_type -> list of hex colors
    """
    if samples is None:
        samples = SAMPLE_ORDER
    if color_map is None:
        color_map = ASSEMBLY_TYPE_COLORS
    
    config = THEMES[theme]
    plt.style.use(config["style"])
    
    # Get assemblers and colors
    assemblers = sorted(rare_df["assembler"].unique())
    assembler_color = assign_colors(assemblers, color_map)
    
    # Create figure
    n_samples = len(samples)
    fig, axes = plt.subplots(1, n_samples, figsize=config["figsize_large"], sharey=True)
    if n_samples == 1:
        axes = [axes]
    
    for ax, sample in zip(axes, samples):
        sample_df = rare_df[rare_df["sample"] == sample].copy()
        
        for assembler in assemblers:
            asm_df = sample_df[sample_df["assembler"] == assembler].dropna(
                subset=["num_contigs", "contig_length_mb"]
            )
            if asm_df.empty:
                continue
            
            # Add origin point
            x_vals = pd.concat(
                [pd.Series([0.0]), asm_df["num_contigs"].astype(float)], ignore_index=True
            )
            y_vals = pd.concat(
                [pd.Series([0.0]), asm_df["contig_length_mb"].astype(float)], ignore_index=True
            )
            
            color = assembler_color[assembler]
            label = ASSEMBLER_LABELS.get(assembler, assembler)
            
            ax.fill_between(x_vals, y_vals, alpha=0.12, color=color)
            ax.plot(x_vals, y_vals, color=color, linewidth=config["linewidth"], label=label)
            ax.scatter(asm_df["num_contigs"], asm_df["contig_length_mb"], 
                      color=color, s=4, zorder=4)
        
        ax.set_title(f"{sample} Rarefaction", fontsize=config["fontsize"] + 1)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, _: f"{int(v):,}" if v >= 1 else f"{v:.2g}")
        )
        ax.set_xlabel("Cumulative contig count (longest → shortest, log scale)", fontsize=config["fontsize"])
        if ax is axes[0]:
            ax.set_ylabel("Cumulative Contig Length (Mb, raw)", fontsize=config["fontsize"])
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f} Mb"))
        ax.grid(True, alpha=0.3)
    
    fig.legend(
        [plt.Line2D([0], [0], color=assembler_color[a], linewidth=config["linewidth"]) 
         for a in assemblers],
        [ASSEMBLER_LABELS.get(a, a) for a in assemblers],
        title="Assembler",
        bbox_to_anchor=(1.01, 0.5),
        loc="center left",
        fontsize=config["fontsize"],
    )
    
    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=config["dpi"], bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved: {outpath}")
