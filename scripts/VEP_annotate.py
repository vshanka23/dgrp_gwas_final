#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import sys

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 5:
    print("Usage: VEP_annotate.py <path to DGRP-3 VEP file> <path to gwas file> <path to gene_symbol file> <path to write out the annot2 file>")
    sys.exit(1)

# Parse command-line arguments
raw_file_path = sys.argv[1]
gwas_file_path = sys.argv[2]
gene_symbol_file_path = sys.argv[3]
output_file_path = sys.argv[4]

# Import VEP data for all positions in DGRP-3
raw = pd.read_csv(raw_file_path, sep="\t")

# Import output from GWAS
gwas = pd.read_csv(gwas_file_path, sep="\t")

# Replace / with _ to use Uploaded_variation column for filtering and merging
raw["#Uploaded_variation"] = raw["#Uploaded_variation"].str.replace("/", "_")

# Filter VEP results from the entire DGRP-3 annotation to those found in the GWAS output
filt = raw[raw["#Uploaded_variation"].isin(gwas["rs"].tolist())]
filt = filt.rename(columns={"#Uploaded_variation": "rs"})

# Merge dataframes
annot = pd.merge(filt, gwas, how='outer', on="rs")

# Import gene_symbol data
gene_symbol = pd.read_csv(gene_symbol_file_path, sep="\t")
gene_symbol = gene_symbol.rename(columns={"primary_FBid": "Gene"})

# Merge with annot dataframe
annot2 = pd.merge(annot, gene_symbol[['Gene', 'current_symbol']], on="Gene", how='left')

# Replace NaNs with "-"
annot2 = annot2.replace(np.NaN, '-')

# Select columns for the final output
annot2 = annot2[['rs', 'Location', 'ps', 'n_miss', 'allele1', 'allele0', 'af', 'beta', 'se', 'logl_H1', 'l_remle',
                 'p_wald', 'Allele', 'Gene', 'current_symbol', 'Feature', 'Feature_type', 'Consequence',
                 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons',
                 'Existing_variation', 'Extra']]

# Rename column
annot2 = annot2.rename(columns={"Allele": "Allele_used_for_VEP"})

#Sort by p_wald
annot2 = annot2.sort_values(by='p_wald', ascending=True)

# Write to the specified output file
annot2.to_csv(output_file_path, header=True, index=False)