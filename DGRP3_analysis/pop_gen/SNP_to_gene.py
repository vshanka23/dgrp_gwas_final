#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys, getopt

def main(argv):
    inputfile1 = ''
    inputfile2 = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:u:o:",["in1=","in2=","out="])
    except getopt.GetoptError:
        print ('SNP_to_gene.py -i <gene info> -u <snv list> -o <output file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('SNP_to_gene.py -i <gene info> -u <snv list> -o <output file>')
            sys.exit()
        elif opt in ("-i", "--in1"):
            gene_info = arg
        elif opt in ("-u", "--in2"):
            snv_list = arg
        elif opt in ("-o", "--out"):
            output = arg
            gene_info_df = pd.read_csv(gene_info,sep="\t")
            snv_list_df = pd.read_csv(snv_list,sep="\t")
            
            gene_coord = pd.concat([gene_info_df['GENE'],gene_info_df['BODY']],axis=1)
            gene_coord[["CHROM","COORD"]] = gene_coord.BODY.str.split(":",expand=True)
            gene_coord[["MIN","MAX"]] = gene_coord.COORD.str.split("-",expand=True)
            gene_coord['MIN'] = pd.to_numeric(gene_coord['MIN'])
            gene_coord['MAX'] = pd.to_numeric(gene_coord['MAX'])

            gene_coord['MIN'] = gene_coord['MIN'] - 1000
            gene_coord['MAX'] = gene_coord['MAX'] + 1000
            snv_list_df.insert(0,"VARIANT",snv_list_df['SNP'])
            snv_list_df[['CHROM','LOC','SNP']] = snv_list_df.SNP.str.split("_",expand=True)

            snv_list_df['GENE']=""
            
            for l in range(0,len(snv_list_df)):
                snv_list_df.iloc[l,4] = (gene_coord[gene_coord['CHROM'].str.contains(snv_list_df.iloc[l,2]) & (pd.to_numeric(snv_list_df.iloc[l,3]) > gene_coord['MIN']) & (pd.to_numeric(snv_list_df.iloc[l,3]) < gene_coord['MAX'])].GENE).str.cat(sep=",")

            snv_list_df.to_csv(output,sep="\t")

if __name__ == "__main__":
    main(sys.argv[1:])