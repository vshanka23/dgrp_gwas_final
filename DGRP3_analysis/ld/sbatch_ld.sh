#!/bin/bash
#
#SBATCH --job-name=DGRP3_LD
#SBATCH --cpus-per-task=4   
#SBATCH --partition=bigmem,zenyatta,dgx
#SBATCH --time=30-00:00:00
#SBATCH --mem=500gb
#SBATCH --output=/scratch2/vshanka/Projects/vshanka_dgrp3_dna/4_varcall/9_ld/log/test_output_%j.txt
#SBATCH --error=/scratch2/vshanka/Projects/vshanka_dgrp3_dna/4_varcall/9_ld/log/test_error_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

cd /scratch2/vshanka/Projects/vshanka_dgrp3_dna/4_varcall/9_ld

VCF=/scratch2/vshanka/Projects/vshanka_dgrp3_dna/4_varcall/DGRP-3_final_biallelic_ID_sample-chr-rename.vcf.gz
BFILE=DGRP-3_final_biallelic_ID
CHR=7

ml plink
plink \
    --allow-extra-chr \
    --vcf $VCF \
    --make-bed \
    --double-id \
    --geno 0.1 \
    --maf 0.01 \
    --chr 1, 2, 3, 4, 5, 6, 7, 8 \
    --out $BFILE

ml plink
plink \
    -bfile $BFILE \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --maf 0.01 \
    --geno 0.1 \
    --mind 0.5 \
    --chr $CHR \
    -r2 gz \
    --ld-window 1000 \
    --ld-window-kb 1000 \
    --ld-window-r2 0 \
    --make-bed \
    --out $BFILE"_LD_"$CHR


python scripts/ld_decay_calc.py -i $BFILE"_LD_"$CHR.ld.gz -o $BFILE"_LD_"$CHR"_MEAN"

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh
mamba activate pandas

python3 scripts/sliding_window_ld.py $BFILE"_LD_"$CHR.ld.gz $BFILE"_LD_"$CHR"_genome_wide.tsv"