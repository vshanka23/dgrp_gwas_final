#Run this first
ml bcftools
ml htslib
bcftools index ../4_varcall/combined_filtered_DGRP-3_rename_final.vcf.gz
bcftools filter -r MT -Ozv ../4_varcall/combined_filtered_DGRP-3_rename_final.vcf.gz > combined_filtered_DGRP-3_rename_final_MT.vcf.gz

#LD_pruned_region is from SNPRelate. Check commands.R in scripts/
bgzip -d combined_filtered_DGRP-3_rename_final_MT.vcf.gz
bgzip -c combined_filtered_DGRP-3_rename_final_MT.vcf > combined_filtered_DGRP-3_rename_final_MT.vcf.gz
bcftools index combined_filtered_DGRP-3_rename_final_MT.vcf.gz
cut -f 3,4 2_phylotree/LD_pruned.txt | tail -n +2 > LD_pruned_region.txt
bcftools view -R LD_pruned_region.txt combined_filtered_DGRP-3_rename_final_MT.vcf.gz -Ovz > MT_LD_pruned.vcf
bgzip -c MT_LD_pruned.vcf > MT_LD_pruned.vcf.gz
bcftools sort -Ovz MT_LD_pruned.vcf.gz > MT_LD_pruned_sorted.vcf
bgzip -c MT_LD_pruned_sorted.vcf > MT_LD_pruned_sorted.vcf.gz
tabix -p vcf MT_LD_pruned_sorted.vcf.gz

python3 scripts/vcf2phylip/vcf2phylip.py --input MT_LD_pruned_sorted.vcf.gz --nexus --fasta