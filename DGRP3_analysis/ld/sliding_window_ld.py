#!/usr/bin/env python3

import gzip
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

# -------------------------
# PARAMETERS
# -------------------------
WINDOW_SIZE = 1_000_000  # 1 Mb
STEP_SIZE   = 100_000    # 100 kb
MIN_DIST    = 50
MAX_DIST    = 150
MIN_PAIRS   = 5          # skip windows with fewer pairs to avoid noise
PLOT_SVG    = True
PLOT_EPS    = True

# -------------------------
# ARGUMENTS
# -------------------------
if len(sys.argv) != 3:
    sys.stderr.write("Usage: sliding_window_ld.py <input.ld.gz> <output.tsv>\n")
    sys.exit(1)

ld_file = sys.argv[1]
out_file = sys.argv[2]

# -------------------------
# READ PLINK LD OUTPUT
# -------------------------
# Columns: CHR BP_A ... BP_B ... R2
df = pd.read_csv(ld_file, compression='gzip', sep=r'\s+', engine='python', header=0)
df = df[['CHR_A','BP_A','BP_B','R2']]
df = df.sort_values(['CHR_A','BP_A','BP_B']).reset_index(drop=True)

# -------------------------
# WINDOW COORDINATES
# -------------------------
chrom = df['CHR_A'].iloc[0]
chr_start = df[['BP_A','BP_B']].min().min()
chr_end   = df[['BP_A','BP_B']].max().max()
window_starts = np.arange(chr_start, chr_end, STEP_SIZE)

# -------------------------
# COMPUTE MEAN R2 PER WINDOW
# -------------------------
results = []

for w_start in window_starts:
    w_end = w_start + WINDOW_SIZE

    # Select pairs inside window
    window_pairs = df[(df['BP_A'] >= w_start) & (df['BP_B'] <= w_end)]
    # Distance filter
    window_pairs = window_pairs[(window_pairs['BP_B'] - window_pairs['BP_A']).between(MIN_DIST, MAX_DIST)]

    if len(window_pairs) < MIN_PAIRS:
        mean_r2 = np.nan
    else:
        mean_r2 = window_pairs['R2'].mean()

    results.append((chrom, w_start, w_end, mean_r2, len(window_pairs)))

# -------------------------
# SAVE TABLE
# -------------------------
out_df = pd.DataFrame(results, columns=['chr','window_start','window_end','mean_r2','n_pairs'])
out_df.to_csv(out_file, sep='\t', index=False)
print(f"Output written to {out_file}")

# -------------------------
# PLOT
# -------------------------
plt.figure(figsize=(12,4))
plt.plot(out_df['window_start'], out_df['mean_r2'], marker='.', linestyle='-')
plt.xlabel("Genomic position (bp)")
plt.ylabel("Mean r²")
plt.title(f"Sliding-window LD ({MIN_DIST}-{MAX_DIST} bp pairs, {WINDOW_SIZE//1_000_000}Mb windows, {STEP_SIZE//1000}kb steps)")
plt.tight_layout()

# Save SVG
if PLOT_SVG:
    svg_file = out_file.rsplit('.',1)[0] + ".svg"
    plt.savefig(svg_file, format='svg')
    print(f"Plot saved as SVG: {svg_file}")

# Save EPS
if PLOT_EPS:
    eps_file = out_file.rsplit('.',1)[0] + ".eps"
    plt.savefig(eps_file, format='eps')
    print(f"Plot saved as EPS: {eps_file}")

plt.close()