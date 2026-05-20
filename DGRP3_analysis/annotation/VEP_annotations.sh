#!/bin/bash
#
#SBATCH --job-name=DGRP3_VEP
#SBATCH --cpus-per-task=40   
#SBATCH --partition=gen-mk-compute-1
#SBATCH --time=10-00:00:00
#SBATCH --mem=200gb
#SBATCH --output=/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/5_annotation/VEP_v5_output_%j.txt
#SBATCH --error=/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/5_annotation/VEP_v5_error_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

unset PERL5LIB

#STRIPPED=/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/5_annotation/DGRP-3_var_info.vcf
#ORIGINAL=/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/DGRP-3_final_biallelic_ID.vcf.gz
#THREADS=40

#ml bcftools

#bcftools view -G -Ov -o $STRIPPED --threads $THREADS $ORIGINAL

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh
mamba activate vep
INPUT=/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/5_annotation/DGRP-3_var_info.vcf
OUTPUT=/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/5_annotation/DGRP-3_final_biallelic.VEP_v3.vcf
SPECIES=drosophila_melanogaster
DIR_CACHE=/data/databases/d_mel/d_mel_6.22_ensembl/VEP_cache/
CV=111
BF=50000
THREADS=40

vep \
-i $INPUT \
-o $OUTPUT \
--vcf \
--species $SPECIES \
--offline \
--cache \
--dir_cache $DIR_CACHE \
--cache_version $CV \
--buffer_size $BF \
--fork $THREADS \
--appris \
--biotype \
--buffer_size 5000 \
--check_existing \
--distance 1000 \
--mane \
--show_ref_allele \
--symbol \
--transcript_version \
--tsl \
--uploaded_allele
