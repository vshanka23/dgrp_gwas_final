#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import argparse
import numpy as np

# ======================
# Parse command line arguments
# ======================
parser = argparse.ArgumentParser(description="Generate LD decay plot from PLINK LD file")
parser.add_argument("-i", "--input", required=True, help="PLINK LD gz file (e.g., chr1_LD.ld.gz)")
parser.add_argument("-o", "--output", default="LD_decay_plot.png", help="Output plot filename")
parser.add_argument("--max_dist", type=int, default=2000, help="Maximum distance to display (bp)")
parser.add_argument("--export", default="LD_decay_values.tsv", help="TSV file to export LD decay values")
args = parser.parse_args()

# ======================
# Load LD data
# ======================
print("Loading LD data...")
with gzip.open(args.input, 'rt') as f:
    ld_df = pd.read_csv(f, delim_whitespace=True)

# PLINK LD columns:
# CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
required_cols = ['CHR_A','BP_A','SNP_A','CHR_B','BP_B','SNP_B','R2']
if not all(col in ld_df.columns for col in required_cols):
    raise ValueError(f"Input file missing required columns. Found columns: {ld_df.columns.tolist()}")

# ======================
# Calculate distance between SNPs
# ======================
ld_df['dist'] = np.abs(ld_df['BP_A'] - ld_df['BP_B'])

# Filter to distances of interest (Drosophila LD decays fast)
ld_df = ld_df[ld_df['dist'] <= args.max_dist]

# ======================
# Bin distances for smoother decay curve
# ======================
bin_size = 10  # 10 bp bins
ld_df['dist_bin'] = (ld_df['dist'] // bin_size) * bin_size

# Compute average R2 per bin
ld_summary = ld_df.groupby('dist_bin')['R2'].mean().reset_index()

# ======================
# Export binned LD decay values
# ======================
ld_summary.to_csv(args.export, sep='\t', index=False)
print(f"LD decay values exported to {args.export}")

# ======================
# Plot LD decay
# ======================
sns.set(style="whitegrid")
plt.figure(figsize=(8,5))
sns.lineplot(data=ld_summary, x='dist_bin', y='R2', marker='o')
plt.xlabel("Distance between SNPs (bp)")
plt.ylabel("Average r²")
plt.title("LD Decay Plot")
plt.ylim(0, 1)
plt.xlim(0, args.max_dist)
plt.tight_layout()
plt.savefig(args.output, dpi=300)
print(f"Plot saved to {args.output}")