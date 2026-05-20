#!/bin/bash
#
#SBATCH --job-name=DGRP3_ND
#SBATCH --cpus-per-task=16   
#SBATCH --partition=gen-mk-compute-1
#SBATCH --time=30-00:00:00
#SBATCH --mem=100gb
#SBATCH --output=/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/10_nuc_div/log/test_output_%j.txt
#SBATCH --error=/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/10_nuc_div/log/test_error_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

set -euo pipefail

############################################
# USER INPUTS
############################################
VCF="/scratch2/vshanka/Projects/vshanka_dgrp3_dna/4_varcall/DGRP-3_final_biallelic.vcf.gz"
GTF="/data/databases/d_mel/d_mel_6.22_ensembl/Drosophila_melanogaster.BDGP6.22.98.gtf"
OUTDIR="/data/Palmetto_sync/Projects/vshanka_dgrp3_dna/4_varcall/10_nuc_div"

WINDOW_SIZE=${4:-10000}
STEP_SIZE=${5:-10000}
GENE_FLANK=${6:-1000}
THREADS=${7:-16}

mkdir -p ${OUTDIR}

############################################
# Load tool modules
############################################

ml bcftools
ml vcftools
ml bedtools2

############################################
# STEP 0: FILTER VCF - MAF filter removed!
############################################

echo "Filtering VCF..."

vcftools --gzvcf ${VCF} \
    --max-missing 0.9 \
    --recode --stdout | bgzip -@ ${THREADS} > ${OUTDIR}/filtered.vcf.gz

tabix -p vcf ${OUTDIR}/filtered.vcf.gz


############################################
# STEP 1: Split SNPs and INDELs
############################################

echo "Splitting SNPs and INDELs..."

bcftools view --threads ${THREADS} -v snps \
    ${OUTDIR}/filtered.vcf.gz -Oz -o ${OUTDIR}/snps.vcf.gz

bcftools view --threads ${THREADS} -v indels \
    ${OUTDIR}/filtered.vcf.gz -Oz -o ${OUTDIR}/indels.vcf.gz

bcftools index --threads ${THREADS} ${OUTDIR}/snps.vcf.gz
bcftools index --threads ${THREADS} ${OUTDIR}/indels.vcf.gz


############################################
# STEP 2: Tajima's D
############################################

echo "Calculating Tajima's D..."

vcftools --gzvcf ${OUTDIR}/snps.vcf.gz \
    --TajimaD ${WINDOW_SIZE} \
    --out ${OUTDIR}/snps_tajimaD

vcftools --gzvcf ${OUTDIR}/indels.vcf.gz \
    --TajimaD ${WINDOW_SIZE} \
    --out ${OUTDIR}/indels_tajimaD

echo "Calculating Tajima's D (1kb windows for gene-level)..."

# SNPs
vcftools --gzvcf ${OUTDIR}/snps.vcf.gz \
    --TajimaD 1000 \
    --out ${OUTDIR}/snps_tajimaD_1kb

# INDELs
vcftools --gzvcf ${OUTDIR}/indels.vcf.gz \
    --TajimaD 1000 \
    --out ${OUTDIR}/indels_tajimaD_1kb


############################################
# STEP 3: Genome-wide Tajima's D (FIXED)
############################################

echo "Calculating genome-wide Tajima's D..."

# SNPs
awk 'NR>1 && $4!="nan" && $4!="NaN" && $4!="." {
    sum+=$4; n++
}
END {
    if(n>0) print sum/n;
    else print "NA"
}' ${OUTDIR}/snps_tajimaD.Tajima.D \
> ${OUTDIR}/snps_tajimaD.genome_mean.txt

# INDELs
awk 'NR>1 && $4!="nan" && $4!="NaN" && $4!="." {
    sum+=$4; n++
}
END {
    if(n>0) print sum/n;
    else print "NA"
}' ${OUTDIR}/indels_tajimaD.Tajima.D \
> ${OUTDIR}/indels_tajimaD.genome_mean.txt


############################################
# STEP 4: Chromosome-level Tajima's D (FIXED)
############################################

echo "Calculating chromosome-level Tajima's D..."

# SNPs
awk 'NR>1 && $4!="nan" && $4!="NaN" && $4!="." {
    chr=$1;
    sum[chr]+=$4;
    count[chr]++;
}
END {
    for (c in sum) {
        print c "\t" sum[c]/count[c];
    }
}' ${OUTDIR}/snps_tajimaD.Tajima.D \
> ${OUTDIR}/snps_tajimaD.per_chrom.txt

# INDELs
awk 'NR>1 && $4!="nan" && $4!="NaN" && $4!="." {
    chr=$1;
    sum[chr]+=$4;
    count[chr]++;
}
END {
    for (c in sum) {
        print c "\t" sum[c]/count[c];
    }
}' ${OUTDIR}/indels_tajimaD.Tajima.D \
> ${OUTDIR}/indels_tajimaD.per_chrom.txt


############################################
# STEP 5: Prepare gene regions (FIXED for GTF parsing)
############################################

echo "Preparing gene regions with names..."

awk '$3=="gene"' ${GTF} | \
awk 'BEGIN{OFS="\t"} {

    start=$4-'"${GENE_FLANK}"';
    if(start<0) start=0;
    end=$5+'"${GENE_FLANK}"';

    # Reconstruct full attribute field (columns 9+)
    attr="";
    for(i=9;i<=NF;i++){
        attr = attr $i " ";
    }

    gene="NA";

    if (match(attr, /gene_name "([^"]+)"/, a)) {
        gene=a[1];
    } else if (match(attr, /gene_id "([^"]+)"/, a)) {
        gene=a[1];
    }

    print $1, start, end, gene

}' > ${OUTDIR}/genes_1kb_named.bed


############################################
# STEP 7: Gene-level Tajima's D (1kb windows → aggregated per gene)
############################################

echo "Calculating gene-level Tajima's D..."

############################################
# 7A: Convert Tajima's D (1kb) to BED
############################################

# SNPs
awk 'NR>1 && $4!="nan" && $4!="NaN" && $4!="." {
    OFS="\t";
    print $1, $2, $2+1000, $4
}' ${OUTDIR}/snps_tajimaD_1kb.Tajima.D \
> ${OUTDIR}/snps_tajimaD_1kb.bed

# INDELs
awk 'NR>1 && $4!="nan" && $4!="NaN" && $4!="." {
    OFS="\t";
    print $1, $2, $2+1000, $4
}' ${OUTDIR}/indels_tajimaD_1kb.Tajima.D \
> ${OUTDIR}/indels_tajimaD_1kb.bed


############################################
# 7B: Intersect windows with genes
############################################

# SNPs
bedtools intersect \
    -a ${OUTDIR}/genes_1kb_named.bed \
    -b ${OUTDIR}/snps_tajimaD_1kb.bed \
    -wa -wb \
> ${OUTDIR}/snps_gene_window_overlap.txt

# INDELs
bedtools intersect \
    -a ${OUTDIR}/genes_1kb_named.bed \
    -b ${OUTDIR}/indels_tajimaD_1kb.bed \
    -wa -wb \
> ${OUTDIR}/indels_gene_window_overlap.txt


############################################
# 7C: Aggregate per gene
############################################

# SNPs
awk '{
    gene=$4;
    val=$8;

    sum[gene]+=val;
    count[gene]++;
}
END {
    for (g in sum) {
        print g "\t" sum[g]/count[g] "\t" count[g];
    }
}' ${OUTDIR}/snps_gene_window_overlap.txt \
> ${OUTDIR}/snps_genes_tajimaD_per_gene.txt


# INDELs
awk '{
    gene=$4;
    val=$8;

    sum[gene]+=val;
    count[gene]++;
}
END {
    for (g in sum) {
        print g "\t" sum[g]/count[g] "\t" count[g];
    }
}' ${OUTDIR}/indels_gene_window_overlap.txt \
> ${OUTDIR}/indels_genes_tajimaD_per_gene.txt


############################################
# 7D: Create tidy combined table
############################################

echo "Creating tidy gene-level table..."

# Add headers
awk 'BEGIN{print "GENE\tTajimaD_SNP\tN_WINDOWS_SNP"} {print}' \
${OUTDIR}/snps_genes_tajimaD_per_gene.txt \
> ${OUTDIR}/snps_tmp.txt

awk 'BEGIN{print "GENE\tTajimaD_INDEL\tN_WINDOWS_INDEL"} {print}' \
${OUTDIR}/indels_genes_tajimaD_per_gene.txt \
> ${OUTDIR}/indels_tmp.txt

# Join on gene name
join -t $'\t' -a1 -a2 -e NA -o auto \
    <(sort -k1,1 ${OUTDIR}/snps_tmp.txt) \
    <(sort -k1,1 ${OUTDIR}/indels_tmp.txt) \
> ${OUTDIR}/gene_tajimaD_combined.txt

# Clean header (fix duplicate column names)
sed -i '1s/.*/GENE\tTajimaD_SNP\tN_WINDOWS_SNP\tTajimaD_INDEL\tN_WINDOWS_INDEL/' \
${OUTDIR}/gene_tajimaD_combined.txt

# Cleanup temp files
rm ${OUTDIR}/snps_tmp.txt ${OUTDIR}/indels_tmp.txt


############################################
# STEP 8: pi
############################################

vcftools --gzvcf ${OUTDIR}/snps.vcf.gz \
    --window-pi ${WINDOW_SIZE} \
    --window-pi-step ${STEP_SIZE} \
    --out ${OUTDIR}/snps_pi

vcftools --gzvcf ${OUTDIR}/indels.vcf.gz \
    --window-pi ${WINDOW_SIZE} \
    --window-pi-step ${STEP_SIZE} \
    --out ${OUTDIR}/indels_pi


############################################
# STEP 9: Genome-wide pi (FIXED)
############################################

echo "Calculating genome-wide pi..."

# SNPs
awk 'NR>1 && $5!="nan" && $5!="NaN" && $5!="." {
    sum+=$5; n++
}
END {
    if(n>0) print sum/n;
    else print "NA"
}' ${OUTDIR}/snps_pi.windowed.pi \
> ${OUTDIR}/snps_pi.genome_mean.txt

# INDELs
awk 'NR>1 && $5!="nan" && $5!="NaN" && $5!="." {
    sum+=$5; n++
}
END {
    if(n>0) print sum/n;
    else print "NA"
}' ${OUTDIR}/indels_pi.windowed.pi \
> ${OUTDIR}/indels_pi.genome_mean.txt


############################################
# STEP 10: Chromosome-level pi (FIXED)
############################################

echo "Calculating chromosome-level pi..."

# SNPs
awk 'NR>1 && $5!="nan" && $5!="NaN" && $5!="." {
    chr=$1;
    sum[chr]+=$5;
    count[chr]++;
}
END {
    for (c in sum) {
        print c "\t" sum[c]/count[c];
    }
}' ${OUTDIR}/snps_pi.windowed.pi \
> ${OUTDIR}/snps_pi.per_chrom.txt

# INDELs
awk 'NR>1 && $5!="nan" && $5!="NaN" && $5!="." {
    chr=$1;
    sum[chr]+=$5;
    count[chr]++;
}
END {
    for (c in sum) {
        print c "\t" sum[c]/count[c];
    }
}' ${OUTDIR}/indels_pi.windowed.pi \
> ${OUTDIR}/indels_pi.per_chrom.txt


############################################
# STEP 11: Combine SNP + INDEL windows and export to Excel
############################################

echo "Combining SNP and INDEL windows and exporting to Excel..."

python3 << EOF
import pandas as pd

outdir = "${OUTDIR}"

############################################
# Load Tajima's D (10kb)
############################################

snps_td = pd.read_csv(f"{outdir}/snps_tajimaD.Tajima.D", sep="\t")
indels_td = pd.read_csv(f"{outdir}/indels_tajimaD.Tajima.D", sep="\t")

# vcftools TajimaD format: CHROM BIN_START N_SNPS TajimaD
snps_td.columns = ["CHROM","BIN_START","N_SNPS","TajimaD_SNP"]
indels_td.columns = ["CHROM","BIN_START","N_INDEL","TajimaD_INDEL"]

############################################
# Merge SNP + INDEL (Tajima's D)
############################################

td = pd.merge(snps_td, indels_td,
              on=["CHROM","BIN_START"],
              how="outer")

td["N_SNPS"] = td["N_SNPS"].fillna(0)
td["N_INDEL"] = td["N_INDEL"].fillna(0)
td["TajimaD_SNP"] = td["TajimaD_SNP"].fillna(0)
td["TajimaD_INDEL"] = td["TajimaD_INDEL"].fillna(0)

# Sort properly
td = td.sort_values(["CHROM","BIN_START"])


############################################
# Load pi (10kb)
############################################

snps_pi = pd.read_csv(f"{outdir}/snps_pi.windowed.pi", sep="\t")
indels_pi = pd.read_csv(f"{outdir}/indels_pi.windowed.pi", sep="\t")

# vcftools pi format: CHROM BIN_START BIN_END N_VARIANTS PI
snps_pi.columns = ["CHROM","BIN_START","BIN_END","N_SNPS","PI_SNP"]
indels_pi.columns = ["CHROM","BIN_START","BIN_END","N_INDEL","PI_INDEL"]

############################################
# Merge SNP + INDEL (pi)
############################################

pi = pd.merge(snps_pi, indels_pi,
              on=["CHROM","BIN_START","BIN_END"],
              how="outer")

pi["N_SNPS"] = pi["N_SNPS"].fillna(0)
pi["N_INDEL"] = pi["N_INDEL"].fillna(0)
pi["PI_SNP"] = pi["PI_SNP"].fillna(0)
pi["PI_INDEL"] = pi["PI_INDEL"].fillna(0)

pi = pi.sort_values(["CHROM","BIN_START"])


############################################
# Excel writer with formatting
############################################

def write_excel(df, filename, include_bin_end=False):

    writer = pd.ExcelWriter(filename, engine="xlsxwriter")
    workbook = writer.book

    header_format = workbook.add_format({
        'bold': True,
        'font_name': 'Arial',
        'font_size': 12,
        'bg_color': 'black',
        'font_color': 'white'
    })

    cell_format = workbook.add_format({
        'font_name': 'Arial',
        'font_size': 12
    })

    # Enforce chromosome order
    chrom_order = ["X", "2L", "2R", "3L", "3R", "4"]
    ordered_chroms = [c for c in chrom_order if c in df["CHROM"].unique()]

    for chrom in ordered_chroms:
        sub = df[df["CHROM"] == chrom].copy()

        # Sort bins within chromosome
        sub = sub.sort_values("BIN_START")

        if include_bin_end:
            sub = sub[["CHROM","BIN_START","BIN_END",
                       "N_SNPS","PI_SNP","N_INDEL","PI_INDEL"]]
        else:
            sub = sub[["CHROM","BIN_START",
                       "N_SNPS","TajimaD_SNP","N_INDEL","TajimaD_INDEL"]]

        sub.to_excel(writer, sheet_name=str(chrom), index=False)

        worksheet = writer.sheets[str(chrom)]

        # Header formatting
        for col_num, value in enumerate(sub.columns.values):
            worksheet.write(0, col_num, value, header_format)

        # Apply Arial 12 to all cells
        worksheet.set_column(0, len(sub.columns)-1, None, cell_format)

    writer.close()


############################################
# Write Excel outputs
############################################

write_excel(td,
            f"{outdir}/TajimaD_10kb_by_chrom.xlsx",
            include_bin_end=False)

write_excel(pi,
            f"{outdir}/Pi_10kb_by_chrom.xlsx",
            include_bin_end=True)

EOF

############################################
echo "Pipeline complete!"