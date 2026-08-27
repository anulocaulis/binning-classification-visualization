# Plotting Experiments Sandbox

**Exploratory figures for data analysis & visualization iteration — NOT final thesis figures.**

This is a **modular playground** for testing different visualization approaches, comparing axis scales, and exploring assembly/binning quality metrics. Figures here are **supplemental/experimental** — useful for understanding data, making suggestions, and iterating quickly. Final publication figures belong in the parent `plotting/` directory.

## Quick Start

### Generate all plots (default theme) - submit to SLURM

```bash
# From project root
sbatch plotting/experiments/sbatch_run_experimental.sh
```

### Generate specific experimental plots with different theme

```bash
# Rarefaction variants (log-log, linear, log-x linear-y)
sbatch plotting/experiments/sbatch_run_experimental.sh --plots rarefaction_variants --theme presentation

# N50 scatter only
sbatch plotting/experiments/sbatch_run_experimental.sh --plots n50_scatter --theme compact

# Original rarefaction + variants
sbatch plotting/experiments/sbatch_run_experimental.sh --plots rarefaction rarefaction_variants --theme default
```

## Architecture

```
experiments/
├── __init__.py                         # Package marker
├── viz_config.py                       # Themes, colors, labels (centralized config)
├── data_loader.py                      # Load & prepare data (reusable)
├── plot_rarefaction.py                 # Original rarefaction (log-x log-y)
├── plot_rarefaction_variants.py        # Rarefaction variants (log-log, linear, log-x linear-y)
├── plot_n50_scatter.py                 # N50 vs L50 scatter
├── plot_threshold_curves.py            # Threshold cutoff impact (per sample)
├── run_experimental.py                 # CLI orchestrator
└── sbatch_run_experimental.sh          # SLURM wrapper (modest resources)
```

## Purpose: Data Exploration & Iteration

This sandbox is designed for:

- **Comparing visualization approaches** — same data, different axes/scales (e.g., rarefaction log-log vs. linear)
- **Rapid iteration** — change a theme parameter, resubmit, see results in ~1-2 min
- **Supplemental analysis** — figures that help *understand* data quality, not necessarily *present* it
- **Per-plot customization** — each plot is independent, easy to modify without refactoring the whole pipeline
- **Publication-quality output** — when you settle on a visualization, migrate it to the parent `plotting/` directory for final figures

**These figures are exploratory.** They're meant to be examined, commented on, iterated on, and compared. The final thesis figures come from `plotting/plot_summary_stats.py`, `plotting/poster_figures.py`, etc.

## Rarefaction Variants

The `rarefaction_variants` plot generates **three different axis configurations** per sample to help assess different aspects of assembly quality:

1. **Log-Log** (`rarefaction_loglog_*.png`)
   - x-axis: log(contig count) — useful for seeing early diversity
   - y-axis: log(cumulative length) — good for power-law relationships
   - Sampling: All first 1000 contigs, then subsample beyond

2. **Linear** (`rarefaction_linear_*.png`)
   - x-axis: contig count (linear)
   - y-axis: cumulative length (linear)
   - Sampling: All first 1000 contigs, then subsample beyond
   - Shows absolute numbers without compression

3. **Log-x Linear-y (First 100)** (`rarefaction_logx_lineary_*.png`)
   - x-axis: log(contig count), limited to first 100 contigs only
   - y-axis: cumulative length (linear)
   - Best for zooming in on the high-quality region
   - Could be used as an inset on log-log if desired

### Generate rarefaction variants

```bash
# All three rarefaction variants with default theme
sbatch plotting/experiments/sbatch_run_experimental.sh --plots rarefaction_variants --theme default

# With presentation theme for publication
sbatch plotting/experiments/sbatch_run_experimental.sh --plots rarefaction_variants --theme presentation
```

## Workflows

### Fast iteration cycle (5 min per idea)

1. **Examine output**: `ls -lh plots/sandbox/`
2. **Get feedback**: "Try smaller fontsize" or "Use viridis colormap"
3. **Edit config**: Update `viz_config.py` or individual plot file
4. **Resubmit**: `sbatch plotting/experiments/sbatch_run_experimental.sh --plots PLOT --theme THEME`
5. **Repeat**: Typically generates output in 1-2 min (modest SLURM resources)

### Exploratory comparison

```bash
# Generate multiple theme variants to compare aesthetics
sbatch plotting/experiments/sbatch_run_experimental.sh --plots rarefaction_variants --theme default
sbatch plotting/experiments/sbatch_run_experimental.sh --plots rarefaction_variants --theme presentation
sbatch plotting/experiments/sbatch_run_experimental.sh --plots rarefaction_variants --theme compact

# Compare side-by-side in plots/sandbox/
```

### Customizing plots for exploration

Edit `viz_config.py` in the THEMES dict to adjust styling:

```python
# experiments/viz_config.py
THEMES = {
    "default": {
        "fontsize": 11,      # <- Change this
        "linewidth": 2,      # <- Or this
        "dpi": 300,
        ...
    },
    ...
}
```

```bash
# Try your changes
sbatch plotting/experiments/sbatch_run_experimental.sh --theme default
```

### Change colors or assembler labels

Edit `ASSEMBLY_TYPE_COLORS` or `ASSEMBLER_LABELS` in `viz_config.py`:

```python
ASSEMBLY_TYPE_COLORS = {
    "short_read": ["#A1D99B", "#74C476", "#238B45"],  # Edit hex colors
    "long_read": ["#C6DBEF", "#41B6C4", "#08519C"],
    ...
}

ASSEMBLER_LABELS = {
    "flye": "Flye",          # <- Update label
    "metaspades": "SPAdes",  # <- Change display name
    ...
}
```

## SLURM Resources

- **CPUs**: 2 (modest, won't compete with 32 CPU binning jobs)
- **Memory**: 4G (won't compete with 360GB binning jobs)
- **Time**: 30 min (plenty for all plots)
- **Output**: `plots/sandbox/sbatch_JOB_ID.{log,err}`

## Output

All plots saved to `plots/sandbox/` by default:

```
plots/sandbox/
├── rarefaction_default.png
├── rarefaction_presentation.png
├── n50_l50_scatter_default.png
├── n50_l50_scatter_presentation.png
├── threshold_S1_num_contigs_default.png
├── threshold_S1_contig_length_mb_default.png
├── threshold_S2_num_contigs_default.png
├── threshold_S2_contig_length_mb_default.png
├── threshold_S5_num_contigs_default.png
├── threshold_S5_contig_length_mb_default.png
└── sbatch_*.{log,err}
```

## Available Themes

| Theme | Use Case | DPI | Font Size | Figure Size |
|-------|----------|-----|-----------|-------------|
| **default** | Analysis/exploration | 300 | 11pt | medium |
| **presentation** | Poster/talk figures | 600 | 13pt | large |
| **compact** | Quick iterations | 200 | 10pt | small |

## From Experiments to Final Figures

When you find a visualization approach that works:

1. **Keep the exploratory figure** in `plots/sandbox/` for reference/documentation
2. **Promote the approach** to `plotting/` if it's going into your thesis:
   - Copy the working code from `experiments/plot_*.py` to a new function in `plotting/plot_summary_stats.py` or create a new module
   - Update configuration to match your final specs (DPI, colors, fonts)
   - Generate final output to `plots/ge1000_poster/` or your thesis figure directory
3. **Archive unsuccessful attempts** (optional): `mkdir plots/sandbox/archive_YYYYMMDD && mv plots/sandbox/*.png plots/sandbox/archive_YYYYMMDD/`

The sandbox is **disposable** — it's meant for exploration. Only code that makes it to `plotting/` or figures in thesis directories are "final."

## Adding New Plots

1. Create `plot_myplot.py` with a `plot_myplot(data_df, outpath, theme="default")` function
2. Import and call in `run_experimental.py`:
   ```python
   if "myplot" in args.plots:
       plot_myplot(summary_df, outpath=outpath, theme=args.theme)
   ```
3. Add to `sbatch_run_experimental.sh` choices if needed
4. Submit!

## Debugging

### View SLURM job output
```bash
tail -f plots/sandbox/sbatch_*.log
```

### Check SLURM queue
```bash
squeue -u $USER
```

### Run locally (not recommended, but possible)
```bash
cd plotting
conda activate visualization
python experiments/run_experimental.py --plots rarefaction --theme default
```

## Notes

- Data is loaded once per run (efficient)
- Individual plots are independent (can skip failing plots, continue with others)
- Themes centralized in `viz_config.py` for consistent styling
- Assembly type classification via `classify_assembly_type()` from `plot_summary_stats.py`
