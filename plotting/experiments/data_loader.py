"""Data loading utilities for experimental visualizations."""
import os
import sys
from typing import Tuple

import pandas as pd

# Add parent plotting directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plot_summary_stats import classify_assembly_type
DEFAULT_NG50_TARGET_BP = 5_000_000  # Kept for compatibility, but not used
ASSEMBLY_TYPE_ORDER = ["long_read", "short_read", "hybrid_read"]
SAMPLE_ORDER = ["S1", "S2", "S5"]


def load_summary_and_threshold_data(
    log_path: str = "rarefaction_curves.csv",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load rarefaction curves data for plotting.
    
    Uses rarefaction_curves.csv which contains one row per contig (built from FASTA files).
    
    Args:
        log_path: Path to rarefaction_curves.csv (one point per contig)
    
    Returns:
        (rarefaction_df, threshold_df) - rarefaction_df has all contig points
    """
    # Load rarefaction curves (one point per contig)
    if not os.path.exists(log_path):
        raise FileNotFoundError(f"Rarefaction curves file not found: {log_path}")
    
    rare_df = pd.read_csv(log_path)
    
    if rare_df.empty:
        raise ValueError(f"No rarefaction data found in {log_path}")
    
    # Filter to S1/S2/S5 only
    rare_df = rare_df[rare_df["sample"].isin(SAMPLE_ORDER)].copy()
    
    # For threshold_df compatibility, create a dummy dataframe
    threshold_df = rare_df[["sample", "assembler", "assembly_type"]].drop_duplicates()
    
    return rare_df, threshold_df


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
