"""Centralized configuration for experimental visualizations."""

# Color themes for different contexts
THEMES = {
    "default": {
        "style": "seaborn-v0_8-darkgrid",
        "fontsize": 11,
        "linewidth": 2,
        "dpi": 300,
        "figsize_small": (9, 5.5),
        "figsize_medium": (10, 6),
        "figsize_large": (12, 8),
    },
    "presentation": {
        "style": "seaborn-v0_8-whitegrid",
        "fontsize": 13,
        "linewidth": 2.5,
        "dpi": 600,
        "figsize_small": (10, 6),
        "figsize_medium": (12, 7),
        "figsize_large": (14, 9),
    },
    "compact": {
        "style": "fast",
        "fontsize": 10,
        "linewidth": 1.5,
        "dpi": 200,
        "figsize_small": (8, 5),
        "figsize_medium": (9, 5.5),
        "figsize_large": (10, 6),
    },
}

# Assembly type color mapping
ASSEMBLY_TYPE_COLORS = {
    "short_read": ["#A1D99B", "#74C476", "#238B45"],
    "long_read": ["#C6DBEF", "#41B6C4", "#08519C"],
    "hybrid_read": ["#D0D1E6", "#6A51A3"],
}

# Assembler label aliases
ASSEMBLER_LABELS = {
    "flye": "Flye",
    "metamdbg": "metaMDBG",
    "idbaud": "IDBA-UD",
    "megahit": "MEGAHIT",
    "metaspades": "metaSPAdes",
    "metaspades_hybrid": "metaSPAdes Hybrid",
    "opera-ms": "Opera-MS",
    "opera_ms": "Opera-MS",
    "spades": "SPAdes",
    "final": "Final",
}

# Default markers for samples
SAMPLE_MARKERS = {
    "S1": "o",
    "S2": "s",
    "S5": "^",
}

# Assembly type order for consistent coloring
ASSEMBLY_TYPE_ORDER = ["long_read", "short_read", "hybrid_read"]
