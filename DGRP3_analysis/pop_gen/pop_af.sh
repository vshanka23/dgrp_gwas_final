# Calculating population allele frequencies

# Global
ml plink
plink   \
  --bfile ../9_ld/DGRP-3_final_biallelic_ID \
  --extract RF_imp_var_chr-pos.txt \
  --freq \
  --out RF_var_imp_AF.txt

# Specific populations

ml plink
plink   \
  --bfile ../9_ld/DGRP-3_final_biallelic_ID \
  --keep PCA_pop_1.txt \
  --extract RF_imp_var_chr-pos.txt \
  --freq \
  --out RF_var_imp_AF_pop1.txt

plink   \
  --bfile ../9_ld/DGRP-3_final_biallelic_ID \
  --keep PCA_pop_2.txt \
  --extract RF_imp_var_chr-pos.txt \
  --freq \
  --out RF_var_imp_AF_pop2.txt

plink   \
  --bfile ../9_ld/DGRP-3_final_biallelic_ID \
  --keep PCA_pop_3.txt \
  --extract RF_imp_var_chr-pos.txt \
  --freq \
  --out RF_var_imp_AF_pop3.txt

  # Merge

awk 'NR>1{print $2, $5}' RF_var_imp_AF_pop1.txt | sort -k1,1 > pop1.maf
awk 'NR>1{print $2, $5}' RF_var_imp_AF_pop2.txt | sort -k1,1 > pop2.maf
awk 'NR>1{print $2, $5}' RF_var_imp_AF_pop3.txt | sort -k1,1 > pop3.maf

join pop1.maf pop2.maf > tmp12.maf
join tmp12.maf pop3.maf > merged.maf

(echo -e "SNP\tPop1\tPop2\tPop3"; cat merged.maf) > final.tsv