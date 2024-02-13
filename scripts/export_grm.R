calculate_kinship <- function(input_file, output_file) {
  if (missing(input_file) || missing(output_file)) {
    stop("Usage: calculate_kinship(input_file, output_file)")
  }
  
  library(genio)
  grm <- read_grm(input_file)
  write.table(grm$kinship, output_file, row.names = FALSE, col.names = FALSE, sep = "\t")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript export.R DGRP-3_plink_filtered_grm DGRP-3_gcta_kinship_inbred.txt")
}

input_file <- args[1]
output_file <- args[2]

calculate_kinship(input_file, output_file)

