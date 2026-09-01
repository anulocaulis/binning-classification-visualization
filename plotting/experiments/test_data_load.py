#!/usr/bin/env python3
"""Test data loading from ge1000_summary_stats.txt"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from plot_summary_stats import parse_threshold_stats_log
import pandas as pd

# Use absolute path
stats_file = os.path.join(os.path.dirname(__file__), '..', 'ge1000_summary_stats.txt')
print(f"Loading from: {stats_file}")
print(f"File exists: {os.path.exists(stats_file)}")

# Test loading the data
df = parse_threshold_stats_log(stats_file)
print("DataFrame shape:", df.shape)
print("\nColumns:", df.columns.tolist())
print("\nFirst few rows:")
print(df.head(10))
print("\nUnique samples:", sorted(df['sample'].unique()))
print("\nUnique assemblers:", sorted(df['assembler'].unique()))
print("\nSample counts for S1/S2/S5:")
for sample in sorted(set(df['sample'].unique()) & {'S1', 'S2', 'S5'}):
    count = len(df[df['sample']==sample])
    print(f"  {sample}: {count} rows")
