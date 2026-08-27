"""Data loading utilities for experimental visualizations."""
import os
import sys
from typing import Tuple

import pandas as pd

# Add parent plotting directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plot_summary_stats import (
    parse_summary_stats_log,
    parse_threshold_stats_log,
    add_recalculated_ng50,
    classify_assembly_type,
)

DEFAULT_NG50_TARGET_BP = 5_000_000
ASSEMBLY_TYPE_ORDER = ["long_read", "short_read", "hybrid_read"]
SAMPLE_ORDER = ["S1", "S2", "S5"]


def load_summary_and_threshold_data(
    log_path: str = "ge1000_summary_stats.txt",
    ng50_target_bp: int = DEFAULT_NG50_TARGET_BP,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and prepare summary + threshold data for plotting.
    
    Returns:
        (summary_df, threshold_df) - Both filtered to S1/S2/S5
    """
    # Parse logs
    summary_df = parse_summary_stats_log(log_path)
    threshold_df = parse_threshold_stats_log(log_path)
    
    if summary_df.empty:
        raise ValueError(f"No summary data found in {log_path}")
    if threshold_df.empty:
        raise ValueError(f"No threshold data found in {log_path}")
    
    # Add NG50 calculations
    summary_df = add_recalculated_ng50(summary_df, ng50_target_bp)
    
    # Add assembly type if missing
    if "assembly_type" not in summary_df.columns:
        summary_df["assembly_type"] = summary_df["assembler"].map(classify_assembly_type).fillna("short_read")
    if "assembly_type" not in threshold_df.columns:
        threshold_df["assembly_type"] = threshold_df["assembler"].map(classify_assembly_type).fillna("short_read")
    
    # Filter to S1/S2/S5
    summary_df = summary_df[summary_df["sample"].isin(SAMPLE_ORDER)].copy()
    threshold_df = threshold_df[threshold_df["sample"].isin(SAMPLE_ORDER)].copy()
    
    return summary_df, threshold_df


def get_assembler_order(df: pd.DataFrame) -> list:
    """Get consistent assembler ordering by type."""
    def sort_key(assembler: str):
        assembly_type = classify_assembly_type(assembler)
        type_rank = (
            ASSEMBLY_TYPE_ORDER.index(assembly_type)
            if assembly_type in ASSEMBLY_TYPE_ORDER
            else len(ASSEMBLY_TYPE_ORDER)
        )
        return (type_rank, str(assembler))
    
    return sorted(df["assembler"].unique(), key=sort_key)
