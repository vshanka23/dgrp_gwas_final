#!/bin/bash

# Genome lenghts are from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001215.4/

# Filter the final biallelic VCF for snps and indels. Then get position (BED format) data
ml bcftools
bcftools view -v snps ../DGRP-3_final_biallelic.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > snps.tsv
bcftools view -v indels ../DGRP-3_final_biallelic.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > indels.tsv

# Extract gene coordinates from GTF
cut -f 1,4,5 /data/databases/d_mel/d_mel_6.22_ensembl/Drosophila_melanogaster.BDGP6.22.98.gtf > genes.tsv

# Filter BED files for canonical chromosome (arms)
grep -e $'^X\t' -e $'^2L\t' -e $'^2R\t' -e $'^3L\t' -e $'^3R\t' -e $'^4\t' genes.tsv > genes_filtered.tsv
grep -e $'^X\t' -e $'^2L\t' -e $'^2R\t' -e $'^3L\t' -e $'^3R\t' -e $'^4\t' snps.tsv > snps_filtered.tsv
grep -e $'^X\t' -e $'^2L\t' -e $'^2R\t' -e $'^3L\t' -e $'^3R\t' -e $'^4\t' indels.tsv > indels_filtered.tsv

# Run the Circos (circlize) script
ml R/4.1.2
Rscript scripts/circos.R --snp snps_filtered.tsv --indel indels_filtered.tsv --gene genes_filtered.tsv --bin 50000 --out test.svg

# Testing with inversions
Rscript scripts/circos_test.R --snp snps_filtered.tsv --indel indels_filtered.tsv --gene genes_filtered.tsv --inv inversions.tsv --bin 50000 --out test_2.svg

# MT commands
# Filter BED files for canonical chromosome (arms)
cut -f 1,4,5 /data/databases/d_mel/d_mel_6.13_ncbi/GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.sorted.gtf | grep -e $'^KJ947872.2\t' > genes_MTfiltered.tsv
grep -e $'^MT\t' snps.tsv > snps_MTfiltered.tsv
grep -e $'^MT\t' indels.tsv > indels_MTfiltered.tsv

# Run the Circos (circlize) script
ml R/4.1.2
Rscript scripts/circos_MT.R --snp snps_MTfiltered.tsv --indel indels_MTfiltered.tsv --gene genes_MTfiltered.tsv --bin 50 --out circos_MT.svg