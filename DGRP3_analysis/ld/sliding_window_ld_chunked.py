#!/usr/bin/env python3

import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

# -------------------------
# PARAMETERS
# -------------------------
WINDOW_SIZE = 1_000_000
STEP_SIZE   = 100_000
MIN_DIST    = 50
MAX_DIST    = 150
MIN_PAIRS   = 5
CHUNK_SIZE  = 1_000_000

# -------------------------
# ARGUMENTS
# -------------------------
if len(sys.argv) != 3:
    sys.stderr.write("Usage: sliding_window_ld_fast.py <input.ld.gz> <output.tsv>\n")
    sys.exit(1)

ld_file = sys.argv[1]
out_file = sys.argv[2]

# -------------------------
# FIRST PASS: get coordinate range
# -------------------------
with gzip.open(ld_file, 'rt') as f:
    header = f.readline()
    coords = []
    for i, line in enumerate(f):
        flds = line.split()
        coords.append((int(flds[1]), int(flds[4])))
        if i == 100000:
            break

coords = np.array(coords)
chr_start = coords.min()
chr_end   = coords.max()

window_starts = np.arange(chr_start, chr_end + STEP_SIZE, STEP_SIZE)
n_windows = len(window_starts)

window_sum   = np.zeros(n_windows, dtype=np.float64)
window_count = np.zeros(n_windows, dtype=np.int64)

# -------------------------
# PROCESS IN CHUNKS
# -------------------------
chunks = pd.read_csv(
    ld_file,
    sep=r'\s+',
    engine='python',
    compression='gzip',
    chunksize=CHUNK_SIZE
)

for chunk in chunks:
    chunk = chunk[['CHR_A','BP_A','BP_B','R2']]

    bp1 = chunk['BP_A'].to_numpy()
    bp2 = chunk['BP_B'].to_numpy()
    r2  = chunk['R2'].to_numpy()

    dist = np.abs(bp2 - bp1)
    mask = (dist >= MIN_DIST) & (dist <= MAX_DIST)

    bp1 = bp1[mask]
    bp2 = bp2[mask]
    r2  = r2[mask]

    if len(bp1) == 0:
        continue

    # Compute window index range per pair
    first_win = (bp2 - WINDOW_SIZE) // STEP_SIZE
    last_win  = bp1 // STEP_SIZE

    first_win = np.clip(first_win, 0, n_windows - 1)
    last_win  = np.clip(last_win, 0, n_windows - 1)

    for i in range(len(bp1)):
        for w in range(first_win[i], last_win[i] + 1):
            window_sum[w]   += r2[i]
            window_count[w] += 1

# -------------------------
# FINAL OUTPUT
# -------------------------
mean_r2 = np.full(n_windows, np.nan)
valid = window_count >= MIN_PAIRS
mean_r2[valid] = window_sum[valid] / window_count[valid]

chrom = chunk['CHR_A'].iloc[0]

out_df = pd.DataFrame({
    'chr': chrom,
    'window_start': window_starts,
    'window_end': window_starts + WINDOW_SIZE,
    'mean_r2': mean_r2,
    'n_pairs': window_count
})

out_df.to_csv(out_file, sep='\t', index=False)

# -------------------------
# PLOT
# -------------------------
plt.figure(figsize=(12,4))
plt.plot(out_df['window_start'], out_df['mean_r2'], '.', lw=1)
plt.xlabel("Genomic position (bp)")
plt.ylabel("Mean r²")
plt.title("Sliding-window LD (50–150 bp pairs)")
plt.tight_layout()

plt.savefig(out_file.replace('.tsv','.svg'))
plt.savefig(out_file.replace('.tsv','.eps'))
plt.close()