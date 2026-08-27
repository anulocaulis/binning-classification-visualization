"""N50 vs L50 scatter plot."""
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

from viz_config import THEMES, ASSEMBLY_TYPE_COLORS, ASSEMBLER_LABELS, SAMPLE_MARKERS, ASSEMBLY_TYPE_ORDER
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


def plot_n50_l50_scatter(
    summary_df: pd.DataFrame,
    outpath: str = "n50_l50_scatter.png",
    theme: str = "default",
    color_map: dict = None,
    show_sample_labels: bool = True,
) -> None:
    """Plot N50 vs L50 with sample annotations.
    
    Args:
        summary_df: DataFrame with columns [sample, assembler, n50_kbp, l50]
        outpath: Output PNG path
        theme: Theme name from THEMES
        color_map: Dict mapping assembly_type -> list of hex colors
        show_sample_labels: Whether to annotate each point with sample ID
    """
    if color_map is None:
        color_map = ASSEMBLY_TYPE_COLORS
    
    config = THEMES[theme]
    plt.style.use(config["style"])
    
    fig, ax = plt.subplots(figsize=config["figsize_medium"])
    
    plot_df = summary_df.dropna(subset=["assembler", "sample", "n50_kbp", "l50"]).copy()
    
    assemblers = sorted(plot_df["assembler"].unique())
    assembler_color = assign_colors(assemblers, color_map)
    
    for assembler in assemblers:
        sub = plot_df[plot_df["assembler"] == assembler]
        if sub.empty:
            continue
        
        ax.scatter(
            sub["l50"],
            sub["n50_kbp"],
            label=ASSEMBLER_LABELS.get(assembler, assembler),
            color=assembler_color[assembler],
            s=100,
            alpha=0.75,
            edgecolors="white",
            linewidths=0.5,
        )
        
        if show_sample_labels:
            for _, row in sub.iterrows():
                ax.annotate(
                    row["sample"],
                    xy=(row["l50"], row["n50_kbp"]),
                    xytext=(4, 3),
                    textcoords="offset points",
                    fontsize=config["fontsize"] - 1,
                    color=assembler_color[assembler],
                    alpha=0.8,
                )
    
    ax.set_title("N50 vs L50", fontsize=config["fontsize"] + 2)
    ax.set_xlabel("L50 (contigs to reach 50% assembly)", fontsize=config["fontsize"])
    ax.set_ylabel("N50 (Kbp)", fontsize=config["fontsize"])
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.legend(title="Assembler", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=config["fontsize"])
    ax.grid(True, alpha=0.3)
    
    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=config["dpi"], bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved: {outpath}")
