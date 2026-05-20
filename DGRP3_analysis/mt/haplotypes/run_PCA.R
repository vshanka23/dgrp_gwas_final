library(SNPRelate)
library(svglite)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: run_PCA.R <vcf_file> <output_directory>")
}
vcf_file <- args[1]
output_dir <- args[2]

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set output file names
gds_file <- file.path(output_dir, "genotype_data.gds")
pca_plot_png <- file.path(output_dir, "PCA_plot.png")
pca_plot_svg <- file.path(output_dir, "PCA_plot.svg")
pca_data_csv <- file.path(output_dir, "PCA_results.csv")
pca_eigenvalues_csv <- file.path(output_dir, "PCA_eigenvalues.csv")
tree_plot_png <- file.path(output_dir, "Phylogenetic_tree.png")
tree_plot_svg <- file.path(output_dir, "Phylogenetic_tree.svg")

# Convert VCF to GDS
cat("Converting VCF to GDS...\n")
snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic.only")

# Open the GDS file
cat("Opening GDS file...\n")
gds <- snpgdsOpen(gds_file)

# Perform LD pruning to reduce redundancy
cat("Performing LD pruning...\n")
snpset <- snpgdsLDpruning(gds, ld.threshold = 0.2,autosome.only=FALSE)
pruned_snps <- unlist(snpset)

# Perform PCA
cat("Running PCA analysis...\n")
pca_result <- snpgdsPCA(gds, snp.id = pruned_snps, autosome.only = FALSE)

# Extract PCA results
pca_variance <- pca_result$varprop * 100
pca_coordinates <- data.frame(SampleID = pca_result$sample.id,
                               PC1 = pca_result$eigenvect[, 1],
                               PC2 = pca_result$eigenvect[, 2],
                               PC3 = pca_result$eigenvect[, 3])

# Save PCA results to CSV
cat("Saving PCA results to CSV...\n")
write.csv(pca_coordinates, pca_data_csv, row.names = FALSE)

# Save PCA eigenvalues (variance) to CSV
cat("Saving PCA eigenvalues to CSV...\n")
pca_eigenvalues <- data.frame(PC = paste0("PC", 1:length(pca_variance)),
                               Variance = pca_variance)
write.csv(pca_eigenvalues, pca_eigenvalues_csv, row.names = FALSE)

# Plot PCA and save as PNG and SVG
cat("Creating PCA plots...\n")

# PNG plot (300 DPI)
png(pca_plot_png, width = 8, height = 6, units = "in", res = 300)
plot(pca_coordinates$PC1, pca_coordinates$PC2,
     xlab = paste0("PC1 (", round(pca_variance[1], 2), "%)"),
     ylab = paste0("PC2 (", round(pca_variance[2], 2), "%)"),
     main = "PCA of Genotype Data", col = "blue", pch = 19)
text(pca_coordinates$PC1, pca_coordinates$PC2,
     labels = pca_coordinates$SampleID, cex = 0.6, pos = 4)
dev.off()

# SVG plot
svglite(pca_plot_svg, width = 8, height = 6)
plot(pca_coordinates$PC1, pca_coordinates$PC2,
     xlab = paste0("PC1 (", round(pca_variance[1], 2), "%)"),
     ylab = paste0("PC2 (", round(pca_variance[2], 2), "%)"),
     main = "PCA of Genotype Data", col = "blue", pch = 19)
text(pca_coordinates$PC1, pca_coordinates$PC2,
     labels = pca_coordinates$SampleID, cex = 0.6, pos = 4)
dev.off()

# Construct phylogenetic tree
cat("Constructing phylogenetic tree...\n")
dist_matrix <- snpgdsDiss(gds, snp.id = pruned_snps,autosome.only=FALSE)
hc <- hclust(as.dist(dist_matrix$diss), method = "average")

# Plot and save phylogenetic tree
cat("Creating phylogenetic tree plots...\n")
png(tree_plot_png, width = 8, height = 6, units = "in", res = 300)
plot(hc, main = "Phylogenetic Tree of Samples", xlab = "Samples", sub = "", cex = 0.8)
dev.off()

svglite(tree_plot_svg, width = 8, height = 6)
plot(hc, main = "Phylogenetic Tree of Samples", xlab = "Samples", sub = "", cex = 0.8)
dev.off()

# Close GDS file
cat("Closing GDS file...\n")
snpgdsClose(gds)

cat("PCA analysis complete. Results saved in:", output_dir, "\n")