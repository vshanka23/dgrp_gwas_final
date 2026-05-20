#!/usr/bin/env Rscript

# Load required libraries
library(tidyverse)
library(Hmisc)
library(optparse)
library(reshape2)

# Define function to compute correlations and export results
run_correlation_analysis <- function(input_file) {
  # Load the input data
  data <- read.csv(input_file, stringsAsFactors = FALSE)

  # Identify sample columns dynamically (all columns after "Species")
  sample_cols <- names(data)[-(1:9)]

  # Subset data for species-level analysis
  species_data_only <- data %>%
    select(Species, all_of(sample_cols))

  # Transpose the data for correlation (species as columns)
  species_transposed <- as.data.frame(t(species_data_only[-1]))
  colnames(species_transposed) <- species_data_only$Species

  # Perform Spearman correlation for species across samples
  species_correlation <- rcorr(as.matrix(species_transposed), type = "spearman")
  # Linearize the correlation matrix
  long_format_r <- melt(species_correlation$r, varnames = c("rowname", "columnname"), value.name = "correlation_value")
  long_format_p <- melt(species_correlation$P, varnames = c("rowname", "columnname"), value.name = "p_val")
  long_format_p$FDR <- p.adjust(long_format_p$p_val,method = "fdr")
  species_network_input <- merge(long_format_r,long_format_p, by = c("rowname","columnname"))
  
  write.table(species_network_input, "microbiome_species_netinput.txt", sep = "\t", row.names = FALSE, quote = FALSE)

  # Now, repeat the process for genus-level analysis

  # Consolidate genus-level data
  genus_data_only <- data %>%
    group_by(Genus) %>%
    summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = 'drop')

  # Transpose genus data for correlation (genus as rows)
  genus_transposed <- as.data.frame(t(genus_data_only[-1]))
  colnames(genus_transposed) <- genus_data_only$Genus

  # Perform Spearman correlation for genus across samples
  genus_correlation <- rcorr(as.matrix(genus_transposed), type = "spearman")

  # Export network input (three-column format) for significant genus correlations
  glong_format_r <- melt(genus_correlation$r, varnames = c("rowname", "columnname"), value.name = "correlation_value")
  glong_format_p <- melt(genus_correlation$P, varnames = c("rowname", "columnname"), value.name = "p_val")
  glong_format_p$FDR <- p.adjust(glong_format_p$p_val,method = "fdr")
  genus_network_input <- merge(glong_format_r,glong_format_p, by = c("rowname","columnname"))

  write.table(genus_network_input, "microbiome_genus_netinput.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript microbiome_correlations.R -i input_file.csv\n")
  quit(status = 1)
}

# Get the input file path
input_file <- args[which(args == "-i") + 1]

# Run the correlation analysis
run_correlation_analysis(input_file)