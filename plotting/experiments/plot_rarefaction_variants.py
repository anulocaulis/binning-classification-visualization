"""Rarefaction curve variants - multiple axes configurations and sampling strategies."""
import os
import numpy as np
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


def format_x_axis_thousands(v, pos):
    """Format x-axis values in thousands with 'K' suffix."""
    if v >= 1000:
        return f"{v/1000:.0f}K"
    elif v >= 1:
        return f"{int(v)}"
    else:
        return f"{v:.1g}"


def subsample_contigs(x_vals: np.ndarray, y_vals: np.ndarray, 
                      dense_threshold: int = 1000) -> tuple:
    """Subsample contigs beyond dense_threshold while keeping all dense ones.
    
    Args:
        x_vals: Array of x values (contig count)
        y_vals: Array of y values (cumulative length)
        dense_threshold: Keep all contigs up to this count, subsample beyond
        
    Returns:
        (subsampled_x, subsampled_y) - Arrays with dense region intact, sparse region subsampled
    """
    if len(x_vals) <= dense_threshold:
        return x_vals, y_vals
    
    # Keep all up to dense_threshold
    dense_mask = x_vals <= dense_threshold
    dense_x = x_vals[dense_mask]
    dense_y = y_vals[dense_mask]
    
    # Subsample beyond dense_threshold (keep every 2nd point to reduce clutter)
    sparse_mask = x_vals > dense_threshold
    sparse_x = x_vals[sparse_mask][::2]  # Every 2nd point
    sparse_y = y_vals[sparse_mask][::2]
    
    # Combine
    combined_x = np.concatenate([dense_x, sparse_x])
    combined_y = np.concatenate([dense_y, sparse_y])
    
    return combined_x, combined_y


def plot_rarefaction_loglog(
    rare_df: pd.DataFrame,
    samples: list = None,
    outpath: str = "rarefaction_loglog.png",
    theme: str = "default",
    color_map: dict = None,
) -> None:
    """Plot rarefaction curves with log-log axes (one figure per sample).
    
    Args:
        rare_df: DataFrame with columns [sample, assembler, num_contigs, contig_length_mb]
        samples: List of samples to plot
        outpath: Output PNG path (template: will be replaced with {sample} for each)
        theme: Theme name from THEMES
        color_map: Dict mapping assembly_type -> list of hex colors
    """
    if samples is None:
        samples = SAMPLE_ORDER
    if color_map is None:
        color_map = ASSEMBLY_TYPE_COLORS
    
    config = THEMES[theme]
    plt.style.use(config["style"])
    
    assemblers = sorted(rare_df["assembler"].unique())
    assembler_color = assign_colors(assemblers, color_map)
    
    # Generate one figure per sample
    for sample in samples:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        sample_df = rare_df[rare_df["sample"] == sample].copy()
        
        for assembler in assemblers:
            asm_df = sample_df[sample_df["assembler"] == assembler].dropna(
                subset=["num_contigs", "contig_length_mb"]
            ).sort_values("num_contigs")
            
            if asm_df.empty:
                continue
            
            x_vals = asm_df["num_contigs"].values.astype(float)
            y_vals = asm_df["contig_length_mb"].values.astype(float)
            
            # For log scale, don't add origin point (log can't represent 0)
            # Just use the data points directly
            x_plot, y_plot = subsample_contigs(x_vals, y_vals)
            
            color = assembler_color[assembler]
            label = ASSEMBLER_LABELS.get(assembler, assembler)
            
            ax.fill_between(x_plot, y_plot, alpha=0.12, color=color)
            ax.plot(x_plot, y_plot, color=color, linewidth=config["linewidth"], label=label)
            ax.scatter(asm_df["num_contigs"], asm_df["contig_length_mb"], 
                      color=color, s=4, zorder=4, alpha=0.6)
        
        ax.set_title(f"{sample} Rarefaction (log-log)", fontsize=config["fontsize"] + 2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(format_x_axis_thousands))
        ax.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, _: f"{v:.0f}" if v >= 1 else f"{v:.1g}")
        )
        ax.set_xlabel("Contig count (log scale)", fontsize=config["fontsize"])
        ax.set_ylabel("Cumulative Length (Mb, log scale)", fontsize=config["fontsize"])
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, which="both")
        
        ax.legend(
            title="Assembler",
            bbox_to_anchor=(1.01, 1),
            loc="upper left",
            fontsize=config["fontsize"],
        )
        
        # Generate output path with sample name
        sample_outpath = outpath.replace(".png", f"_{sample}.png")
        
        fig.tight_layout()
        os.makedirs(os.path.dirname(sample_outpath), exist_ok=True)
        fig.savefig(sample_outpath, dpi=config["dpi"], bbox_inches="tight")
        plt.close(fig)
        print(f"✓ Saved: {sample_outpath}")


def plot_rarefaction_linear(
    rare_df: pd.DataFrame,
    samples: list = None,
    outpath: str = "rarefaction_linear.png",
    theme: str = "default",
    color_map: dict = None,
) -> None:
    """Plot rarefaction curves with linear-linear axes (one figure per sample).
    
    Args:
        rare_df: DataFrame with columns [sample, assembler, num_contigs, contig_length_mb]
        samples: List of samples to plot
        outpath: Output PNG path (template: will be replaced with {sample} for each)
        theme: Theme name from THEMES
        color_map: Dict mapping assembly_type -> list of hex colors
    """
    if samples is None:
        samples = SAMPLE_ORDER
    if color_map is None:
        color_map = ASSEMBLY_TYPE_COLORS
    
    config = THEMES[theme]
    plt.style.use(config["style"])
    
    assemblers = sorted(rare_df["assembler"].unique())
    assembler_color = assign_colors(assemblers, color_map)
    
    # Generate one figure per sample
    for sample in samples:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        sample_df = rare_df[rare_df["sample"] == sample].copy()
        
        for assembler in assemblers:
            asm_df = sample_df[sample_df["assembler"] == assembler].dropna(
                subset=["num_contigs", "contig_length_mb"]
            ).sort_values("num_contigs")
            
            if asm_df.empty:
                continue
            
            x_vals = asm_df["num_contigs"].values.astype(float)
            y_vals = asm_df["contig_length_mb"].values.astype(float)
            
            # Add origin
            x_vals = np.concatenate([[0.0], x_vals])
            y_vals = np.concatenate([[0.0], y_vals])
            
            # Subsample for cleaner curves
            x_plot, y_plot = subsample_contigs(x_vals, y_vals)
            
            color = assembler_color[assembler]
            label = ASSEMBLER_LABELS.get(assembler, assembler)
            
            ax.fill_between(x_plot, y_plot, alpha=0.12, color=color)
            ax.plot(x_plot, y_plot, color=color, linewidth=config["linewidth"], label=label)
            ax.scatter(asm_df["num_contigs"], asm_df["contig_length_mb"], 
                      color=color, s=4, zorder=4, alpha=0.6)
        
        ax.set_title(f"{sample} Rarefaction (linear)", fontsize=config["fontsize"] + 2)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(format_x_axis_thousands))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f} Mb"))
        ax.set_xlabel("Contig count", fontsize=config["fontsize"])
        ax.set_ylabel("Cumulative Length (Mb)", fontsize=config["fontsize"])
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3)
        
        ax.legend(
            title="Assembler",
            bbox_to_anchor=(1.01, 1),
            loc="upper left",
            fontsize=config["fontsize"],
        )
        
        # Generate output path with sample name
        sample_outpath = outpath.replace(".png", f"_{sample}.png")
        
        fig.tight_layout()
        os.makedirs(os.path.dirname(sample_outpath), exist_ok=True)
        fig.savefig(sample_outpath, dpi=config["dpi"], bbox_inches="tight")
        plt.close(fig)
        print(f"✓ Saved: {sample_outpath}")


def plot_rarefaction_logx_lineary(
    rare_df: pd.DataFrame,
    first_n_contigs: int = 100,
    samples: list = None,
    outpath: str = "rarefaction_logx_lineary.png",
    theme: str = "default",
    color_map: dict = None,
) -> None:
    """Plot rarefaction curves with log-x and linear-y axes for shoulder detection (one figure per sample).
    
    Args:
        rare_df: DataFrame with columns [sample, assembler, num_contigs, contig_length_mb]
        first_n_contigs: Limit to first N contigs for shoulder detection (100, 200, 500, 1000)
        samples: List of samples to plot
        outpath: Output PNG path (template: will be replaced with {sample} for each)
        theme: Theme name from THEMES
        color_map: Dict mapping assembly_type -> list of hex colors
    """
    if samples is None:
        samples = SAMPLE_ORDER
    if color_map is None:
        color_map = ASSEMBLY_TYPE_COLORS
    
    config = THEMES[theme]
    plt.style.use(config["style"])
    
    assemblers = sorted(rare_df["assembler"].unique())
    assembler_color = assign_colors(assemblers, color_map)
    
    # Generate one figure per sample
    for sample in samples:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        sample_df = rare_df[rare_df["sample"] == sample].copy()
        
        for assembler in assemblers:
            asm_df = sample_df[sample_df["assembler"] == assembler].dropna(
                subset=["num_contigs", "contig_length_mb"]
            ).sort_values("num_contigs")
            
            if asm_df.empty:
                continue
            
            # Keep only first N contigs
            asm_df = asm_df.head(first_n_contigs)
            
            x_vals = asm_df["num_contigs"].values.astype(float)
            y_vals = asm_df["contig_length_mb"].values.astype(float)
            
            # Add origin at (1, 0) since log scale needs positive x
            x_vals = np.concatenate([[1.0], x_vals])
            y_vals = np.concatenate([[0.0], y_vals])
            
            color = assembler_color[assembler]
            label = ASSEMBLER_LABELS.get(assembler, assembler)
            
            ax.fill_between(x_vals, y_vals, alpha=0.12, color=color)
            ax.plot(x_vals, y_vals, color=color, linewidth=config["linewidth"], label=label)
            # Scatter plot without the artificial origin
            scatter_x = asm_df["num_contigs"].values
            scatter_y = asm_df["contig_length_mb"].values
            ax.scatter(scatter_x, scatter_y, color=color, s=6, zorder=4, alpha=0.6)
        
        ax.set_title(f"{sample} Rarefaction (first {first_n_contigs}, log x)", fontsize=config["fontsize"] + 2)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(format_x_axis_thousands))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f} Mb"))
        ax.set_xlabel("Contig count (log scale)", fontsize=config["fontsize"])
        ax.set_ylabel("Cumulative Length (Mb)", fontsize=config["fontsize"])
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, which="both")
        
        ax.legend(
            title="Assembler",
            bbox_to_anchor=(1.01, 1),
            loc="upper left",
            fontsize=config["fontsize"],
        )
        
        # Generate output path with sample name
        sample_outpath = outpath.replace(".png", f"_{sample}.png")
        
        fig.tight_layout()
        os.makedirs(os.path.dirname(sample_outpath), exist_ok=True)
        fig.savefig(sample_outpath, dpi=config["dpi"], bbox_inches="tight")
        plt.close(fig)
        print(f"✓ Saved: {sample_outpath}")


def plot_rarefaction_logx_lineary_first100(
    rare_df: pd.DataFrame,
    samples: list = None,
    outpath: str = "rarefaction_logx_lineary.png",
    theme: str = "default",
    color_map: dict = None,
) -> None:
    """Convenience wrapper for first 100 contigs (backward compatibility)."""
    return plot_rarefaction_logx_lineary(
        rare_df, first_n_contigs=100, samples=samples, outpath=outpath, theme=theme, color_map=color_map
    )
