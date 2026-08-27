"""Threshold curve plots - contig retention at different size cutoffs."""
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

from viz_config import THEMES, ASSEMBLY_TYPE_COLORS, ASSEMBLER_LABELS, ASSEMBLY_TYPE_ORDER

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


def _format_threshold_label(threshold_kbp: float) -> str:
    """Format threshold value for axis label."""
    if threshold_kbp == 0:
        return "All"
    if threshold_kbp < 1:
        return f"{int(threshold_kbp * 1000)}bp"
    if threshold_kbp < 1000:
        return f"{int(threshold_kbp)}kbp" if threshold_kbp == int(threshold_kbp) else f"{threshold_kbp}kbp"
    mbp = threshold_kbp / 1000.0
    return f"{int(mbp)}mbp" if mbp == int(mbp) else f"{mbp}mbp"


def plot_threshold_curves(
    threshold_df: pd.DataFrame,
    sample: str,
    metric: str = "num_contigs",
    outpath: str = "threshold_curve.png",
    theme: str = "default",
    color_map: dict = None,
) -> None:
    """Plot contig retention vs minimum contig length threshold.
    
    Args:
        threshold_df: DataFrame with columns [sample, assembler, threshold_kbp, num_contigs, contig_length_mb]
        sample: Sample to plot
        metric: "num_contigs" or "contig_length_mb"
        outpath: Output PNG path
        theme: Theme name from THEMES
        color_map: Dict mapping assembly_type -> list of hex colors
    """
    if color_map is None:
        color_map = ASSEMBLY_TYPE_COLORS
    
    config = THEMES[theme]
    plt.style.use(config["style"])
    
    sample_df = threshold_df[threshold_df["sample"] == sample].copy()
    if sample_df.empty:
        raise ValueError(f"No data for sample {sample}")
    
    fig, ax = plt.subplots(figsize=config["figsize_medium"])
    
    assemblers = sorted(sample_df["assembler"].unique())
    assembler_color = assign_colors(assemblers, color_map)
    
    # Build threshold ticks
    threshold_vals = sorted(sample_df["threshold_kbp"].dropna().unique())
    threshold_labels = [_format_threshold_label(t) for t in threshold_vals]
    
    for assembler in assemblers:
        asm_df = sample_df[sample_df["assembler"] == assembler].dropna(
            subset=["threshold_kbp", metric]
        ).sort_values("threshold_kbp")
        
        if asm_df.empty:
            continue
        
        ax.plot(
            asm_df["threshold_kbp"],
            asm_df[metric],
            marker="o",
            markersize=5,
            linewidth=config["linewidth"],
            color=assembler_color[assembler],
            label=ASSEMBLER_LABELS.get(assembler, assembler),
        )
    
    metric_label = "Contig Count" if metric == "num_contigs" else "Contig Length (Mb)"
    ax.set_title(f"{sample} - {metric_label} vs Minimum Length Threshold", fontsize=config["fontsize"] + 1)
    ax.set_xlabel("Minimum Contig Length", fontsize=config["fontsize"])
    ax.set_ylabel(metric_label, fontsize=config["fontsize"])
    ax.set_xticks(threshold_vals)
    ax.set_xticklabels(threshold_labels, rotation=45, ha="right", fontsize=config["fontsize"] - 1)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.legend(title="Assembler", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=config["fontsize"])
    ax.grid(True, alpha=0.3)
    
    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=config["dpi"], bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved: {outpath}")
