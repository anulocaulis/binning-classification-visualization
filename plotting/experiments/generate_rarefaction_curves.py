#!/usr/bin/env python3
"""Generate rarefaction curves from FASTA files.

Reads all contigs from FASTA files and generates complete rarefaction curves
(one point per contig, sorted by length descending).

Output: rarefaction_curves.csv with columns:
  sample, assembler, assembly_type, num_contigs, contig_length_mb

Usage:
    python generate_rarefaction_curves.py [ASSEMBLIES_DIR] [OUTFILE]

Example:
    python generate_rarefaction_curves.py \
        /storage/biology/projects/miller-lowry/beitner/data/assemblies \
        rarefaction_curves.csv
"""
import argparse
import os
import pandas as pd
from pathlib import Path

# Add parent plotting directory to path
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plot_summary_stats import classify_assembly_type


def read_fasta_lengths(fasta_path: str) -> list:
    """Read contig lengths from FASTA file.
    
    Args:
        fasta_path: Path to FASTA file
    
    Returns:
        List of contig lengths in base pairs (sorted descending)
    """
    lengths = []
    current_length = 0
    
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if current_length > 0:
                        lengths.append(current_length)
                    current_length = 0
                else:
                    current_length += len(line)
            
            # Don't forget last contig
            if current_length > 0:
                lengths.append(current_length)
    except Exception as e:
        print(f"Error reading {fasta_path}: {e}")
        return []
    
    # Sort descending (longest first) for rarefaction curve
    return sorted(lengths, reverse=True)


def build_rarefaction_curve(contig_lengths: list) -> dict:
    """Build rarefaction curve from contig lengths.
    
    Args:
        contig_lengths: List of contig lengths (already sorted descending)
    
    Returns:
        Dict with 'num_contigs' and 'contig_length_mb' keys (lists of cumulative values)
    """
    cumsum = 0
    num_contigs_list = []
    length_mb_list = []
    
    for i, length_bp in enumerate(contig_lengths, 1):
        cumsum += length_bp
        num_contigs_list.append(i)
        length_mb_list.append(cumsum / 1_000_000.0)
    
    return {
        'num_contigs': num_contigs_list,
        'contig_length_mb': length_mb_list,
    }


def generate_rarefaction_curves(assemblies_dir: str, outfile: str) -> None:
    """Generate rarefaction curves for all assemblies.
    
    Args:
        assemblies_dir: Root directory containing samples (S1, S2, S5, etc.)
        outfile: Output CSV file path
    """
    records = []
    assemblies_path = Path(assemblies_dir)
    
    # Find all FASTA files
    for fasta_file in sorted(assemblies_path.glob("*/assembly.*/contigs.ge1000.fa")):
        parts = fasta_file.parts
        # Extract sample and assembler from path
        sample_idx = -3  # S1, S2, S5, etc.
        assembler_idx = -2  # assembly.XXX
        
        sample = parts[sample_idx]
        assembler_dir = parts[assembler_idx]
        
        # Extract assembler name from directory like "assembly.metaspades"
        if assembler_dir.startswith("assembly."):
            assembler = assembler_dir.replace("assembly.", "")
        else:
            assembler = assembler_dir
        
        assembly_type = classify_assembly_type(assembler)
        
        print(f"Processing: {sample}/{assembler_dir} ({fasta_file})")
        
        # Read contigs
        lengths = read_fasta_lengths(str(fasta_file))
        if not lengths:
            print(f"  ⚠ No contigs found, skipping")
            continue
        
        # Build rarefaction curve
        curve = build_rarefaction_curve(lengths)
        
        # Create records for each point on the curve
        for num_contigs, contig_length_mb in zip(curve['num_contigs'], curve['contig_length_mb']):
            records.append({
                'sample': sample,
                'assembler': assembler,
                'assembly_type': assembly_type,
                'num_contigs': num_contigs,
                'contig_length_mb': contig_length_mb,
            })
        
        print(f"  ✓ {len(lengths)} contigs, {len(records)} curve points")
    
    # Save to CSV
    df = pd.DataFrame(records)
    os.makedirs(os.path.dirname(outfile) if os.path.dirname(outfile) else '.', exist_ok=True)
    df.to_csv(outfile, index=False)
    print(f"\n✓ Saved {len(df)} records to {outfile}")
    print(f"  {df['sample'].nunique()} samples × {df['assembler'].nunique()} assemblers")


def main():
    parser = argparse.ArgumentParser(
        description="Generate rarefaction curves from FASTA files"
    )
    parser.add_argument(
        "assemblies_dir",
        default="/storage/biology/projects/miller-lowry/beitner/data/assemblies",
        nargs="?",
        help="Root directory containing assembly FASTA files",
    )
    parser.add_argument(
        "outfile",
        default="rarefaction_curves.csv",
        nargs="?",
        help="Output CSV file path",
    )
    
    args = parser.parse_args()
    generate_rarefaction_curves(args.assemblies_dir, args.outfile)


if __name__ == "__main__":
    main()
