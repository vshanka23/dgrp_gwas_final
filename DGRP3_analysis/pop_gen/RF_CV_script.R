#R/4.3.0 (requires a lot of memory >32GB)

library(caret)
library(doParallel)
setwd("/scratch2/vshanka/Projects/vshanka_dgrp3_dna/4_varcall/6_rel_lines/")
load("RF.RData")
#cores <- detectCores(cl)
cores <- 40
cl <- makeCluster(cores)
registerDoParallel(cl)

# Full model
training_parameters <- trainControl(method="cv", p = 1.0, number = 20, classProbs = TRUE, savePredictions = TRUE, allowParalle=TRUE)
discriminant_CV <- train(x = t(geno), y = group, method = 'rf', tuneLength = 15, trControl = training_parameters)
save("RF_CV_full_model.RData")

# Partial model with 216 variables
import <- read.table("important_SNV.txt",header=TRUE)
geno_subset <- geno[rownames(geno) %in% import$Subset, ]

training_parameters <- trainControl(method="cv", p = 1.0, number = 20, classProbs = TRUE, savePredictions = TRUE, allowParalle=TRUE)
discriminant_CV <- train(x = t(geno_subset), y = group, method = 'rf', tuneLength = 15, trControl = training_parameters)


library(pROC)
library(dplyr)

# 1. Filter for the best mtry
plot_data <- discriminant_CV$pred %>% filter(mtry == 2) # <-chosen based on best performance 

# 2. Create the ROC objects
roc_G1 <- roc(response = (plot_data$obs == "G1"), predictor = plot_data$G1)
roc_G2 <- roc(response = (plot_data$obs == "G2"), predictor = plot_data$G2)
roc_G3 <- roc(response = (plot_data$obs == "G3"), predictor = plot_data$G3)

# 3. Calculate Confidence Intervals for each (this may take a few seconds)
ci_G1 <- ci.se(roc_G1, specificities = seq(0, 1, 0.01))
ci_G2 <- ci.se(roc_G2, specificities = seq(0, 1, 0.01))
ci_G3 <- ci.se(roc_G3, specificities = seq(0, 1, 0.01))


# 1. Open the SVG device
# You can change the filename and dimensions (in inches)
pdf("ROC_Multiclass_Results.pdf", width = 8, height = 7)

# 2. Run your plotting code (ensure plot_data and ROC objects exist)
# --- START OF PLOTTING ---
plot(roc_G1, col="red", lwd=2, main="Full model: One-vs-All ROC (mtry=7714)")
#plot(ci_G1, type="shape", col=rgb(1, 0, 0, 0.15))

plot(roc_G2, add=TRUE, col="blue", lwd=2)
#plot(ci_G2, type="shape", col=rgb(0, 0, 1, 0.15))

plot(roc_G3, add=TRUE, col="green", lwd=2)
#plot(ci_G3, type="shape", col=rgb(0, 1, 0, 0.15))

legend("bottomright", 
       legend=c(paste("G1 AUC:", round(auc(roc_G1), 2)),
                paste("G2 AUC:", round(auc(roc_G2), 2)),
                paste("G3 AUC:", round(auc(roc_G3), 2))),
       col=c("red", "blue", "green"), lwd=2, bty="n")
# --- END OF PLOTTING ---

# 3. Close the device to save the file
dev.off()

# Multipanel figure

library(pROC)
library(dplyr)

# 1. Filter for the best mtry
plot_data <- discriminant_CV$pred %>% 
  filter(mtry == 2)   # chosen based on best performance

# 2. Create the ROC objects
roc_G1 <- roc(response = (plot_data$obs == "G1"), predictor = plot_data$G1)
roc_G2 <- roc(response = (plot_data$obs == "G2"), predictor = plot_data$G2)
roc_G3 <- roc(response = (plot_data$obs == "G3"), predictor = plot_data$G3)

# 3. Calculate confidence intervals
ci_G1 <- ci.se(roc_G1, specificities = seq(0, 1, 0.01))
ci_G2 <- ci.se(roc_G2, specificities = seq(0, 1, 0.01))
ci_G3 <- ci.se(roc_G3, specificities = seq(0, 1, 0.01))

# 4. Open PDF device
pdf("ROC_subsetmodel_Multiclass_Results_panels.pdf", width = 15, height = 5)

# 5. Set multi-panel layout: 1 row, 3 columns
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# Panel 1
plot(roc_G1,
     col = "red",
     lwd = 2,
     main = paste0("G1 vs All\nAUC = ", round(auc(roc_G1), 5)),
     xlab = "1 - Specificity",
     ylab = "Sensitivity")
# plot(ci_G1, type = "shape", col = rgb(1, 0, 0, 0.15))

# Panel 2
plot(roc_G2,
     col = "blue",
     lwd = 2,
     main = paste0("G2 vs All\nAUC = ", round(auc(roc_G2), 5)),
     xlab = "1 - Specificity",
     ylab = "Sensitivity")
# plot(ci_G2, type = "shape", col = rgb(0, 0, 1, 0.15))

# Panel 3
plot(roc_G3,
     col = "green3",
     lwd = 2,
     main = paste0("G3 vs All\nAUC = ", round(auc(roc_G3), 5)),
     xlab = "1 - Specificity",
     ylab = "Sensitivity")
# plot(ci_G3, type = "shape", col = rgb(0, 1, 0, 0.15))

# 6. Close device
dev.off()