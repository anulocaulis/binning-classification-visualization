#!/usr/bin/env python3
"""Orchestrate experimental plot generation with configurable themes - submit via sbatch."""
import argparse
import os
import sys
from pathlib import Path

# Add plotting module to path
plotting_dir = Path(__file__).parent.parent
sys.path.insert(0, str(plotting_dir))

from data_loader import load_summary_and_threshold_data, SAMPLE_ORDER
from plot_rarefaction import plot_rarefaction
from plot_rarefaction_variants import (
    plot_rarefaction_loglog,
    plot_rarefaction_linear,
    plot_rarefaction_logx_lineary,
    plot_rarefaction_logx_lineary_first100,
)
from plot_n50_scatter import plot_n50_l50_scatter
from plot_threshold_curves import plot_threshold_curves


def main():
    parser = argparse.ArgumentParser(
        description="Generate experimental plots with configurable themes (submit via sbatch)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Submit to sbatch with default theme (all plots)
  sbatch experiments/sbatch_run_experimental.sh
  
  # Submit rarefaction only with presentation theme
  sbatch experiments/sbatch_run_experimental.sh --plots rarefaction --theme presentation
  
  # Submit N50 scatter with compact theme
  sbatch experiments/sbatch_run_experimental.sh --plots n50_scatter --theme compact
  
  # Run locally (not recommended, use sbatch instead):
  python experiments/run_experimental.py --plots rarefaction --theme default
        """,
    )
    
    parser.add_argument(
        "--plots",
        nargs="+",
        choices=["rarefaction", "rarefaction_variants", "n50_scatter", "threshold"],
        default=["rarefaction", "n50_scatter", "threshold"],
        help="Which plots to generate (default: all)",
    )
    parser.add_argument(
        "--theme",
        choices=["default", "presentation", "compact"],
        default="default",
        help="Visualization theme (default: default)",
    )
    parser.add_argument(
        "--outdir",
        default="plots/sandbox",
        help="Output directory relative to project root (default: plots/sandbox)",
    )
    parser.add_argument(
        "--log-path",
        default="rarefaction_curves.csv",
        help="Path to rarefaction curves CSV file (relative to plotting dir)",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        default=SAMPLE_ORDER,
        help="Samples to include (default: S1 S2 S5)",
    )
    
    args = parser.parse_args()
    
    # Ensure we're in the plotting directory
    if not os.path.exists(args.log_path):
        plotting_path = os.path.join(plotting_dir, args.log_path)
        if os.path.exists(plotting_path):
            os.chdir(plotting_dir)
        else:
            print(f"Error: Could not find {args.log_path}", file=sys.stderr)
            print(f"  Checked: {args.log_path}", file=sys.stderr)
            print(f"  Checked: {plotting_path}", file=sys.stderr)
            sys.exit(1)
    
    print(f"\n{'='*60}")
    print(f"Loading data from: {os.path.abspath(args.log_path)}")
    print(f"{'='*60}\n")
    
    try:
        summary_df, threshold_df = load_summary_and_threshold_data(args.log_path)
        print(f"✓ Loaded {len(summary_df)} rarefaction records")
        print(f"  (from {len(summary_df['assembler'].unique())} assembler×sample combinations)")
    except ValueError as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    print(f"✓ Output directory: {os.path.abspath(args.outdir)}\n")
    
    # Generate requested plots
    print(f"Generating plots with theme='{args.theme}':\n")
    
    plot_count = 0
    
    if "rarefaction" in args.plots:
        outpath = os.path.join(args.outdir, f"rarefaction_{args.theme}.png")
        try:
            plot_rarefaction(summary_df, samples=args.samples, outpath=outpath, theme=args.theme)
            plot_count += 1
        except Exception as e:
            print(f"✗ Rarefaction plot failed: {e}", file=sys.stderr)
    
    if "rarefaction_variants" in args.plots:
        # Generate log-log and linear variants
        for variant_name, plot_func in [
            ("loglog", plot_rarefaction_loglog),
            ("linear", plot_rarefaction_linear),
        ]:
            outpath = os.path.join(args.outdir, f"rarefaction_{variant_name}_{args.theme}.png")
            try:
                plot_func(summary_df, samples=args.samples, outpath=outpath, theme=args.theme)
                plot_count += 1
            except Exception as e:
                print(f"✗ Rarefaction {variant_name} plot failed: {e}", file=sys.stderr)
        
        # Generate shoulder-detection log-x linear-y variants at different contig thresholds
        for contig_limit in [100, 200, 500, 1000]:
            variant_name = f"logx_lineary_first{contig_limit}"
            outpath = os.path.join(args.outdir, f"rarefaction_{variant_name}_{args.theme}.png")
            try:
                plot_rarefaction_logx_lineary(
                    summary_df,
                    first_n_contigs=contig_limit,
                    samples=args.samples,
                    outpath=outpath,
                    theme=args.theme,
                )
                plot_count += 1
            except Exception as e:
                print(f"✗ Rarefaction {variant_name} plot failed: {e}", file=sys.stderr)
    
    if "n50_scatter" in args.plots:
        outpath = os.path.join(args.outdir, f"n50_l50_scatter_{args.theme}.png")
        try:
            plot_n50_l50_scatter(summary_df, outpath=outpath, theme=args.theme)
            plot_count += 1
        except Exception as e:
            print(f"✗ N50/L50 scatter plot failed: {e}", file=sys.stderr)
    
    if "threshold" in args.plots:
        for sample in args.samples:
            for metric in ["num_contigs", "contig_length_mb"]:
                outpath = os.path.join(
                    args.outdir,
                    f"threshold_{sample}_{metric}_{args.theme}.png"
                )
                try:
                    plot_threshold_curves(
                        threshold_df,
                        sample=sample,
                        metric=metric,
                        outpath=outpath,
                        theme=args.theme,
                    )
                    plot_count += 1
                except Exception as e:
                    print(f"✗ Threshold plot for {sample}/{metric} failed: {e}", file=sys.stderr)
    
    print(f"\n{'='*60}")
    print(f"✓ Generated {plot_count} plots")
    print(f"✓ Outputs in: {os.path.abspath(args.outdir)}/")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
