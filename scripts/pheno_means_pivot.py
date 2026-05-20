#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Flatten PLINK .qassoc.means into one row per SNP and export as TSV")
    parser.add_argument("--input", required=True, help="Input .qassoc.means file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    # Read file
    df = pd.read_csv(args.input, delim_whitespace=True)

    # Keep only SNP, VALUE, G11, G22 (ignore CHR and G12)
    df = df[["SNP", "VALUE", "G11", "G22"]]

    # Pivot the table: VALUE types become columns
    wide = df.pivot(index="SNP", columns="VALUE", values=["G11", "G22"])

    # Flatten MultiIndex columns
    wide.columns = [f"{allele}_{val}" for allele, val in wide.columns]

    # Build output DataFrame
    out = pd.DataFrame(index=wide.index)
    out["SNP"] = wide.index
    out["G11"] = wide["G11_GENO"]
    out["G22"] = wide["G22_GENO"]
    out["G11_Count"] = wide["G11_COUNTS"].astype(float)
    out["G22_Count"] = wide["G22_COUNTS"].astype(float)
    out["G11_Mean"] = wide["G11_MEAN"].astype(float)
    out["G22_Mean"] = wide["G22_MEAN"].astype(float)

    # SE = SD / sqrt(Count), round to 3 decimals
    out["G11_SE"] = (wide["G11_SD"].astype(float) / out["G11_Count"]**0.5).round(3)
    out["G22_SE"] = (wide["G22_SD"].astype(float) / out["G22_Count"]**0.5).round(3)

    # Reorder columns
    out = out[["SNP", "G11", "G22", "G11_Count", "G22_Count", "G11_Mean", "G22_Mean", "G11_SE", "G22_SE"]]

    # Write TSV
    out.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
