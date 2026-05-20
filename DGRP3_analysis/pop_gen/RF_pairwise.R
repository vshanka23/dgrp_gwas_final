# Pairwise RF

# Load Libraries
library("randomForest")

# Load data
setwd("/scratch2/vshanka/Projects/vshanka_dgrp3_dna/4_varcall/6_rel_lines/")
load("RF.RData")

# Subset genotype matrix by comparison
geno_G1_G2 <- t(geno)[group %in% c("G1", "G2"), ]
geno_G1_G3 <- t(geno)[group %in% c("G1", "G3"), ]
geno_G2_G3 <- t(geno)[group %in% c("G2", "G3"), ]

# Subset group by comparison
group_1_2 <- group[group %in% c("G1", "G2")]
group_1_3 <- group[group %in% c("G1", "G3")]
group_2_3 <- group[group %in% c("G2", "G3")]

# Run individual RF models
RF_1_2 <- randomForest(geno_G1_G2,group=group_1_2,importance=TRUE,proximity=TRUE)
RF_1_3 <- randomForest(geno_G1_G3,group=group_1_3,importance=TRUE,proximity=TRUE)
RF_2_3 <- randomForest(geno_G2_G3,group=group_2_3,importance=TRUE,proximity=TRUE)
save.image("RF_pairwise.RData")