#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
import re

def reformat_data(input_file, output_file):
    # Read data from CSV
    data = pd.read_csv(input_file)

    # Format by delimiting and appending
    data[["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]] = data["Unnamed: 0"].str.split("|", expand=True)

    # Remove extraneous details
    data["Kingdom"] = data.Kingdom.str.replace("k__", "")
    data["Phylum"] = data.Phylum.str.replace("p__", "")
    data["Class"] = data.Class.str.replace("c__", "")
    data["Order"] = data.Order.str.replace("o__", "")
    data["Family"] = data.Family.str.replace("f__", "")
    data["Genus"] = data.Genus.str.replace("g__", "")
    data["Species"] = data.Species.str.replace("s__", "")
    
    pattern = r'(DGRP-[2-3]-[0-9]+_S\d+)'

    sample_columns = [col for col in data.columns if re.match(pattern, col)]
        
    # Check new data structure
    data_reform = data[["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"] + sample_columns]

    # Save reformatted data to CSV
    data_reform.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat Metaphlan data.")
    parser.add_argument("-i", "--input", help="Input CSV file", required=True)
    parser.add_argument("-o", "--output", help="Output CSV file", required=True)
    args = parser.parse_args()

    reformat_data(args.input, args.output)
