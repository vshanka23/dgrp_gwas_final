# Function to filter and format covariates to use with GEMMA
filter_input <- function(inv_file, wolb_file, pheno_file, filtered_pheno_name, filtered_covar_name, sig_covar_name, indiv_IDs_name, threshold) {
  # Read into R and format
  inv <- read.table(inv_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1)
  wolb <- read.table(wolb_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1)
  pheno <- read.table(pheno_file, header=FALSE, sep="\t", stringsAsFactors=FALSE, row.names=1)
  
  inv[,1] <- NULL
  wolb[,1] <- NULL
  pheno[,1] <- NULL
  
  # Combine covariates. Wolbachia is dummy coded to [0,1] to represent the binary variable.
  covar <- cbind(wolb, inv[row.names(wolb),])
  
  # Filter combined covariate file for rows (line IDs) present in the phenotype file
  filtered_covar_temp <- na.omit(covar[row.names(pheno),])
  filtered_pheno <- na.omit(pheno[row.names(filtered_covar_temp), , drop = FALSE])
  indiv_IDs <- rownames(filtered_pheno)
  filtered_covar <- cbind(c(rep(1, length(rownames(filtered_pheno)))), filtered_covar_temp)

  # For loop to loop through each covariate, fit a GLM against the phenotype and test with a Wald test (Chi-square test or Log-ratio test) and save p-values
  output <- numeric()
  
  # Iterate over columns of covar
  for (i in 1:ncol(filtered_covar_temp)) {
    # Run the ANOVA and extract the p-value
    p_value <- anova(glm(filtered_pheno[, 1] ~ filtered_covar_temp[, i], family = gaussian()), test = "Chisq")$`Pr(>Chi)`[2]
    
    # Append the p-value to the output vector
    output[i] <- p_value
  }
  
  output <- as.data.frame(output)
  row.names(output) <- colnames(filtered_covar_temp)
  
  sig_covar <- cbind(c(rep(1, length(rownames(filtered_pheno)))),filtered_covar_temp[, row.names(output)[output[,1] < threshold], drop=FALSE])
  
  # Save significant covariates to the output file
  write.table(sig_covar, file=sig_covar_name, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Write the merged data to the output file
  write.table(cbind(cbind(rownames(filtered_pheno)),filtered_pheno), file = filtered_pheno_name, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
  write.table(filtered_covar, file = filtered_covar_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(indiv_IDs, file = indiv_IDs_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 8) {
  stop("Usage: Rscript processing.R inv_file wolb_file pheno_file filtered_pheno filtered_covar sig_covar indiv_IDs threshold")
}

inv_file <- args[1]
wolb_file <- args[2]
pheno_file <- args[3]
filtered_pheno_name <- args[4]
filtered_covar_name <- args[5]
sig_covar_name <- args[6]
indiv_IDs_name <- args[7]
threshold <- as.numeric(args[8])

filter_input(inv_file, wolb_file, pheno_file, filtered_pheno_name, filtered_covar_name, sig_covar_name, indiv_IDs_name, threshold)
