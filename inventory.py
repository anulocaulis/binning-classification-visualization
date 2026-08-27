import os
import subprocess

base_dir = "/storage/biology/projects/miller-lowry/beitner/data/binning_outputs"

# 1) Get/Print immediate subdirectories of data/binning_outputs
subdirs = sorted([d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))])
print("1) Immediate subdirectories of data/binning_outputs:")
for d in subdirs:
    print(f"  {d}")
print("\n" + "="*80 + "\n")

# 2) For each subdirectory, count directories named exactly bin_fasta within maxdepth 6
# 3) Print up to first 8 example bin_fasta paths per subdirectory
print("2 & 3) Bin Fasta Inventory (maxdepth 6):")
for d in subdirs:
    subdir_path = os.path.join(base_dir, d)
    # Run find command
    cmd = ["find", subdir_path, "-maxdepth", "6", "-type", "d", "-name", "bin_fasta"]
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        paths = [line.strip() for line in res.stdout.splitlines() if line.strip()]
    except Exception as e:
        paths = []
    
    count = len(paths)
    print(f"Subdirectory: {d}")
    print(f"  Total bin_fasta count: {count}")
    if count > 0:
        print("  Example paths (up to 8):")
        for p in paths[:8]:
            print(f"    - {p}")
    print()

print("="*80 + "\n")

# 4) Also count directories named fast_mode_bins and print up to first 8 examples globally.
print("4) Fast Mode Bins Inventory (Global search in binning_outputs):")
cmd_global = ["find", base_dir, "-type", "d", "-name", "fast_mode_bins"]
try:
    res_global = subprocess.run(cmd_global, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    global_paths = [line.strip() for line in res_global.stdout.splitlines() if line.strip()]
except Exception as e:
    global_paths = []

print(f"  Total fast_mode_bins count: {len(global_paths)}")
if global_paths:
    print("  Example paths (up to 8):")
    for p in global_paths[:8]:
        print(f"    - {p}")
