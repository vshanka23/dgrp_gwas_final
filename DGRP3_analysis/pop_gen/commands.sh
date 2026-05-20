# Bash code

ml plink
plink \
    --allow-extra-chr \
    --vcf ../DGRP-3_final_biallelic_ID_sample-chr-rename.vcf.gz \
    --make-bed \
    --double-id \
    --geno 0.1 \
    --maf 0.01 \
    --chr 1, 2, 3, 4, 5, 6, 7, 8 \
    --out DGRP-3_plink

plink \
    --allow-extra-chr \
    --bfile DGRP-3_plink \
    --indep-pairwise 50 5 0.2 \
    --out DGRP-3_plink_LD

plink \
    --allow-extra-chr \
    --bfile DGRP-3_plink \
    --extract DGRP-3_plink_LD.prune.in \
    --make-bed \
    --out DGRP-3_plink_LD_pruned

# R code

library("genio")
library("rrBLUP")
#library("sommer")
data <- read_plink("DGRP-3_plink_LD_pruned",verbose=TRUE)
geno <- as.data.frame(data$X)
write.table(geno,"DGRP-3_plink_genotype_matrix.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

# Run PCA
pca_result <- prcomp(t(geno),center = TRUE, scale.=TRUE)
write.table(pca_result$x[,1:3],"DGRP-3_plink_pca_x.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
write.table(pca_result$rotation[,1:3],"DGRP-3_plink_pca_rotation.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)

# Create a grouping variable using the PCA x matrix (PCs) based on Line separation

# Run RF
library(randomForest)
group <- read.table("DGRP-3_RF_grouping.txt",header=TRUE,row.names=1)
group <- as.factor(group$Group)
RF <- randomForest(t(geno),group=group,importance=TRUE,proximity=TRUE)
write.table(RF$importance,"DGRP-3_RF_3pop_importance_scores.txt", sep="\t", quote=FALSE,row.names=TRUE,col.names=NA)

# Construct GRM for visualization
geno[geno == 0] <- -1
geno[geno == 1] <- 0
geno[geno == 2] <- 1
X <- as.matrix(t(geno))
data_a <- A.mat(X,return.imputed=FALSE)
diag_vals <- diag(data_a)
# Dividing by the diagonal leads to warping of the GRM. Instead, I have divided by 2.
#scaled_GRM <- data_a / diag_vals
scaled_GRM <- data_a / 2
write.table(scaled_GRM,"DGRP-3_plink_LD_pruned_GRM.txt",sep="\t")

# Bash code

Rscript standard_GRM/process_grm.R -i scaled_GRM/DGRP-3_plink_LD_pruned_GRM2.txt -o scaled_GRM/DGRP-3_plink_LD_pruned_GRM2_linear.txt