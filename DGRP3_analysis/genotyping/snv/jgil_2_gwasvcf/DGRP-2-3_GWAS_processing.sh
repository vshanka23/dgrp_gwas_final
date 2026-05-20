#!/bin/bash
#
#SBATCH --job-name=DGRP
#SBATCH --cpus-per-task=8   
#SBATCH --partition=bigmem
#SBATCH --time=10-00:00:00
#SBATCH --mem=320gb
#SBATCH --output=/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/DGRP-2-3_out_%j.txt
#SBATCH --error=/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/DGRP-2-3_err_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

cd /data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall

ml bcftools
ml gatk
ml java_jdk
ml htslib


list=$(bcftools query -l combined_filtered_jjremoved.vcf | grep "DGRP-*-*" | tr '\n' ',' | sed 's/,$//')

bcftools view \
    -Oz \
    --threads 8 \
    -s $list \
    combined_filtered_jjremoved.vcf.gz > \
    combined_filtered_DGRP-2-3.vcf.gz

bcftools index \
    combined_filtered_DGRP-2-3.vcf.gz

#Filter for specific chromosomes

bcftools filter \
    --threads 8 \
    combined_filtered_DGRP-2-3.vcf.gz \
    -r AE014134.6,AE013599.5,AE014296.5,AE014297.3,AE014135.4,AE014298.5,CP007106.1,KJ947872.2 \
    -Ozv > \
    combined_filtered_DGRP-2-3_final.vcf.gz

#Rename chromsomomes (I created the rename_chr_map.txt from reference genome page on NCBI)

bcftools annotate \
    -Ozv \
    --threads 8 \
    --rename-chrs rename_chr_map.txt \
    combined_filtered_DGRP-2-3_final.vcf.gz > \
    combined_filtered_DGRP-2-3_rename_final.vcf.gz

#Remove monomorphic sites

bcftools view \
    --threads 8 \
    -Ozv \
    -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
    combined_filtered_DGRP-2-3_rename_final.vcf.gz > \
    DGRP-2-3_rename_polymorphic.vcf.gz

#Filter for biallelic sites

bcftools view \
    -Oz \
    --threads 8 \
    -m2 \
    -M2 \
    DGRP-2-3_rename_polymorphic.vcf.gz > \
    DGRP-2-3_final_biallelic.vcf.gz

#Generates variant ID for each variant based on chromosome, position, reference and alt alleles

bcftools annotate \
    -Ozv \
    --threads 8 \
    --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
    DGRP-2-3_final_biallelic.vcf.gz > \
    DGRP-2-3_final_biallelic_ID.vcf.gz

#Renaming chromosomes to use PLINK/GCTA correctly (they expect human reference chromosome designations)

bcftools annotate \
    -Ozv \
    --threads 8 \
    --rename-chrs rename_chr_gcta.txt \
    DGRP-2-3_final_biallelic_ID.vcf.gz > \
    DGRP-2-3_final_biallelic_ID_chr-rename.vcf.gz

#Rename samples to use PLINK correctly (the current naming convention was dumping _S# into sample ID and keeps the line number as family ID because of PLINK's default behavior)

bcftools query \
    -l \
    combined_filtered_jjremoved.vcf | \
    grep "DGRP-*-*" > \
    DGRP-2-3_sample_rename_map.txt

bcftools query \
    -l combined_filtered_jjremoved.vcf | \
    grep "DGRP-*-*" | \
    sed 's/_S.*//' | \
    paste DGRP-2-3_sample_rename_map.txt - > \
    DGRP-2-3_sample_rename_map_1.txt

rm DGRP-2-3_sample_rename_map.txt

mv DGRP-2-3_sample_rename_map_1.txt DGRP-2-3_sample_rename_map.txt

bcftools reheader \
    DGRP-2-3_final_biallelic_ID_chr-rename.vcf.gz \
    --threads 8 \
    --samples DGRP-2-3_sample_rename_map.txt \
    --output DGRP-2-3_final_biallelic_ID_sample-chr-rename.vcf.gz

gatk IndexFeatureFile \
    -I DGRP-2-3_final_biallelic_ID_sample-chr-rename.vcf.gz

gatk VariantFiltration \
    -V DGRP-2-3_final_biallelic_ID_sample-chr-rename.vcf.gz \
    -O DGRP-2-3_final_biallelic_ID_sample-chr-rename_nohetflag.vcf.gz \
    --genotype-filter-expression "isHet == 1" \
    --genotype-filter-name "isHetFilter"

gatk SelectVariants \
    -V DGRP-2-3_final_biallelic_ID_sample-chr-rename_nohetflag.vcf.gz \
    --set-filtered-gt-to-nocall \
    -O DGRP-2-3_final_biallelic_ID_sample-chr-rename_nohet.vcf.gz