#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# mt_haplotype_pipeline_full.R
# End-to-end pipeline from a multi-sample mitochondrial VCF to:
#  - Per-sample consensus FASTA (using bcftools)
#  - Optional extraction of MT-only reference from a whole-genome reference
#  - Masking of AT-rich control region (default: 14915-19524 for D. melanogaster)
#  - Full-length masked alignment (FASTA)
#  - Variant-only alignment (FASTA) for PopART
#  - PCA + model-based clustering (mclust with BIC)
#  - fastbaps sequence clustering (if installed)
#  - rhierBAPS hierarchical clustering (if installed)
#  - IQ-TREE maximum-likelihood tree (if iqtree is installed)
#
# Notes:
#  * Requires command-line tools: bcftools, tabix, samtools, mafft; iqtree optional.
#  * R packages: ape, Biostrings, seqinr, mclust, cluster, fastbaps (optional), rhierbaps (optional)
#  * This script is defensive: it will subset VCF to MT and bgzip+index it.
#  * Mask coordinates are configurable; defaults suit D. melanogaster (see README for refs).
# ----------------------------------------------------------------------------
.libPaths("/opt/ohpc/pub/libs/gnu9/R/4.1.2/lib64/R/library")
suppressPackageStartupMessages({
  library(parallel)
  library(ape)
  library(Biostrings)
  library(seqinr)
  library(mclust)
  library(cluster)
})

# ----------------------------- USER SETTINGS ------------------------------- #
VCF          <- "/scratch2/vshanka/Projects/vshanka_dgrp3_dna/8_mt_hap/combined_filtered_DGRP-3_rename_final_MT.vcf.gz"   # Multi-sample VCF (bgzipped; .tbi index present or creatable)
REF          <- "GCA_000001215.4_Release_6_plus_ISO1_MT-only.fna"      # Reference FASTA (whole genome or MT-only). Must be indexed by samtools faidx.
OUTDIR       <- "mt_pipeline_results"# Output directory
CHR          <- "MT"                 # Name of mitochondrial contig in VCF/REF (e.g., "MT", "chrM", "M"...)
QUAL_MIN     <- 30                   # QUAL filter for variants in VCF subsetting step

# If REF is whole-genome, set this TRUE to create an MT-only reference on-the-fly (recommended)
EXTRACT_MT_REF <- TRUE

# Mask regions: comma-separated ranges "start-end" (1-based, inclusive). Default covers D. mel control region.
MASK_REGIONS_STR <- "5960-6060,14915-19524"

# mclust settings
PC_USE <- 6      # number of PCs (or MDS axes) used by mclust
MAX_G  <- 20     # test up to this many mixture components

# Toggle optional analyses
RUN_FASTBAPS   <- TRUE   # requires fastbaps R package
RUN_RHIERBAPS  <- TRUE   # requires rhierbaps R package
RUN_IQTREE     <- TRUE   # requires iqtree available on PATH
IQTREE_THREADS <- "AUTO" # or a number like "8"

# Parallelization: number of cores for consensus generation (set to 1 for sequential)
# Default: detectCores() - 1 (but not less than 1)
N_CORES <- max(1, parallel::detectCores(logical = TRUE) - 1)

# --------------------------- END USER SETTINGS ----------------------------- #

# Helpers ------------------------------------------------------------------- #
run <- function(cmd) {
  cat(">>", cmd, "\n");
  st <- system(cmd);
  if (st != 0) stop("Command failed: ", cmd)
}

ensure_tool <- function(name) {
  ok <- suppressWarnings(system2("bash", c("-lc", sprintf("command -v %s >/dev/null 2>&1", name))))
  if (ok != 0) stop(sprintf("Required tool '%s' not found in PATH.", name))
}

# Check external tools we depend on
invisible(ensure_tool("bcftools"))
invisible(ensure_tool("tabix"))
invisible(ensure_tool("samtools"))
if (RUN_IQTREE) {
  ok <- suppressWarnings(system2("bash", c("-lc", "command -v iqtree >/dev/null 2>&1")))
  if (ok != 0) { warning("iqtree not found; disabling RUN_IQTREE."); RUN_IQTREE <<- FALSE }
}

# Create output dir
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Index FASTA if needed
if (!file.exists(paste0(REF, ".fai"))) run(sprintf("samtools faidx %s", shQuote(REF)))

# Optionally extract MT-only reference (recommended if REF is whole-genome)
MT_REF <- REF
if (EXTRACT_MT_REF) {
  MT_REF <- file.path(OUTDIR, sprintf("%s_only.fa", CHR))
  run(sprintf("samtools faidx %s %s > %s", shQuote(REF), shQuote(CHR), shQuote(MT_REF)))
  if (!file.exists(MT_REF) || file.info(MT_REF)$size == 0) stop("Failed to extract MT contig. Check CHR name vs REF.")
  run(sprintf("samtools faidx %s", shQuote(MT_REF)))
}

# Subset VCF to MT + QUAL filter; bgzip + index
vcf_mt <- file.path(OUTDIR, "mt.filtered.vcf.gz")
run(sprintf("bcftools view -r %s -i 'QUAL>=%d' -Oz -o %s %s", shQuote(CHR), QUAL_MIN, shQuote(vcf_mt), shQuote(VCF)))
run(sprintf("tabix -p vcf %s", shQuote(vcf_mt)))

# Normalize VCF to avoid overlapping variants (left-align indels, merge MNPs)
vcf_mt_norm <- file.path(OUTDIR, "mt.norm.vcf.gz")
run(sprintf("bcftools norm -f %s -m +any -Oz -o %s %s", shQuote(MT_REF), shQuote(vcf_mt_norm), shQuote(vcf_mt)))
run(sprintf("tabix -p vcf %s", shQuote(vcf_mt_norm)))

# Get sample list
samples_file <- file.path(OUTDIR, "samples.txt")
run(sprintf("bcftools query -l %s > %s", shQuote(vcf_mt_norm), shQuote(samples_file)))
samples <- readLines(samples_file)
cat("Found", length(samples), "samples in VCF.\n")
if (length(samples) == 0) stop("No samples found in VCF after filtering.")

# Write mask BED file for bcftools consensus
mask_bed <- file.path(OUTDIR, "mask.bed")
mask_ranges <- strsplit(MASK_REGIONS_STR, ",")[[1]]
mask_lines <- vapply(mask_ranges, function(r) {
  se <- as.integer(strsplit(r, "-")[[1]])
  sprintf("%s\t%d\t%d", CHR, se[1]-1, se[2])  # BED is 0-based, end-exclusive
}, character(1))
writeLines(mask_lines, con = mask_bed)

# Generate per-sample consensus FASTA for MT contig (using MT_REF + mask)
cat("Generating per-sample consensus sequences (parallel)...\n")

fa_files <- mclapply(samples, function(s) {
  fa_out <- file.path(OUTDIR, paste0(s, ".fa"))
  cmd <- sprintf("bcftools consensus -s %s -f %s -m %s %s > %s",
              shQuote(s), shQuote(MT_REF), shQuote(mask_bed), shQuote(vcf_mt_norm), shQuote(fa_out))
  ret <- system(cmd)
  if (ret != 0) stop(sprintf("bcftools consensus failed for sample %s (exit %d)", s, ret))
  sed_ret <- system(sprintf("sed -i '1s/>.*/>%s/' %s", s, shQuote(fa_out)))
  if (sed_ret != 0) warning(sprintf("sed header replace returned %d for %s", sed_ret, s))
  return(fa_out)
}, mc.cores = N_CORES)

fa_files <- unlist(fa_files)

# Combine per-sample FASTAs into one multi-FASTA
combined_fa <- file.path(OUTDIR, "mt_consensus_all.fasta")
run(sprintf("cat %s > %s", paste(shQuote(fa_files), collapse = " "), shQuote(combined_fa)))

# Align with MAFFT ------------------------------------------------------------
aligned_fa <- file.path(OUTDIR, "mt_consensus_all_aligned.fasta")
cat("Running MAFFT alignment...\n")
run(sprintf("mafft --auto --thread %d %s > %s", N_CORES, shQuote(combined_fa), shQuote(aligned_fa)))

# Masking function: applies gaps ('-') over specified ranges (1-based inclusive)
mask_alignment <- function(dna_set, mask_str) {
  if (is.null(mask_str) || mask_str == "") return(dna_set)
  ranges <- strsplit(mask_str, ",")[[1]]
  seq_len_all <- unique(width(dna_set))
  if (length(seq_len_all) != 1) stop("Sequences have differing lengths; ensure consensus used the same reference.")
  L <- seq_len_all
  masked <- dna_set
  for (r in ranges) {
    se <- as.integer(strsplit(trimws(r), "-")[[1]])
    if (length(se) != 2) stop("Bad mask range: ", r)
    start <- max(1, se[1]); end <- min(L, se[2])
    if (start > end) next
    for (i in seq_along(masked)) {
      subseq(masked[[i]], start = start, end = end) <- DNAString(paste(rep("-", end - start + 1), collapse = ""))
    }
  }
  masked
}

# Read aligned FASTA
cat("Reading aligned FASTA...\n")
dna <- readDNAStringSet(aligned_fa)
if (length(dna) != length(samples)) warning("Sample count mismatch after FASTA read.")

# Masking function (unchanged)
cat("Masking regions: ", MASK_REGIONS_STR, "\n")
dna_masked <- mask_alignment(dna, MASK_REGIONS_STR)

# Write full-length masked FASTA
masked_fasta <- file.path(OUTDIR, "mt_masked_full.fasta")
writeXStringSet(dna_masked, filepath = masked_fasta, format = "fasta")
cat("Full-length masked FASTA written: ", masked_fasta, "\n")

# Variant-only alignment for PopART / networks --------------------------------
cat("Building variant-only alignment...\n")
# Convert to matrix of characters
aln_char <- as.character(as.matrix(as.DNAbin(dna_masked)))
# Identify variable columns (exclude gaps-only columns)
var_cols <- which(apply(aln_char, 2, function(col) {
  u <- unique(col)
  u <- u[u != 'N' & u != '-' ]
  length(u) > 1
}))
if (length(var_cols) == 0) stop("No variable sites found after masking.")
variant_mat <- aln_char[, var_cols, drop = FALSE]
variant_fasta <- file.path(OUTDIR, "mt_variants_only.fasta")
seqs_list <- apply(variant_mat, 1, function(x) paste(x, collapse = ""))
write.fasta(as.list(seqs_list), names = rownames(variant_mat), file.out = variant_fasta)
cat("Variant-only FASTA written: ", variant_fasta, "\n")

library(ape)

# Path to your variant-only FASTA
variant_fasta <- file.path(OUTDIR,"mt_variants_only.fasta")

# Read FASTA
dna_var <- read.dna(variant_fasta, format = "fasta")

# Optional: check dimensions
cat("Number of sequences:", nrow(dna_var), "\n")
cat("Number of sites:", ncol(dna_var), "\n")

# Write as Nexus
nexus_file <- file.path(OUTDIR,"mt_variants_only.nex")
write.nexus.data(dna_var, file = nexus_file, format = "dna")
cat("Nexus file written:", nexus_file, "\n")

# Distance, ordination (cmdscale), and mclust ---------------------------------
cat("Computing pairwise distances and MDS (PC-like axes)...\n")
dna_masked <- readDNAStringSet(aligned_fa) # <- skipping masking 
dmat <- dist.dna(as.DNAbin(dna_masked), model = "raw", as.matrix = TRUE)
mds  <- cmdscale(as.dist(dmat), k = max(PC_USE, 10), eig = TRUE)
pc   <- mds$points
rownames(pc) <- names(dna_masked)

cat("Running mclust (BIC selects K)...\n")
pc_use <- pc[, 1:PC_USE, drop = FALSE]
mc <- Mclust(pc_use, G = 1:MAX_G)
# Save BIC plot
pdf(file.path(OUTDIR, "mclust_BIC_plot.pdf"))
plot(mc, what = "BIC")
dev.off()
# Choose the correct model from the plot
mc <- Mclust(pc_use, G = 16, modelNames = "EEV") # <- EEV model has the second highest BIC at G = 12
# Cluster assignments and diagnostics
mclust_assign <- mc$classification
sil <- silhouette(mclust_assign, dist(pc_use))
write.csv(data.frame(sample = rownames(pc_use),
                     mclust_cluster = mclust_assign,
                     uncertainty = mc$uncertainty),
          file = file.path(OUTDIR, "mclust_clusters.csv"), row.names = FALSE)
cat(sprintf("mclust chose G=%d; avg silhouette=%.3f\n", mc$G, mean(sil[,3])))


library(fastbaps)
library(rhierbaps)

# Use your unmasked multiple sequence alignment from MAFFT
msa_file <- file.path(OUTDIR, "mt_consensus_all_aligned.fasta")

snp_sparse <- import_fasta_sparse_nt(msa_file)
snp_sparse <- optimize_prior(snp_sparse, type = "optimized.baps")

# 2. Perform hierarchical clustering (this gives an hclust tree)
snp_hclust <- fast_baps(snp_sparse)

# 3. Optimise partitioning of the tree
partitions <- best_baps_partition(snp_sparse, snp_hclust)

# 4. The result is a vector of cluster IDs, named by sample
cluster_df <- data.frame(sample = colnames(snp_sparse$snp.matrix), cluster = partitions)

head(cluster_df)
write.csv(cluster_df, "fastbaps_clusters.csv", row.names = FALSE)

# Load the data
snp_matrix <- load_fasta(file.path(OUTDIR, "mt_consensus_all_aligned.fasta"))
# Perform clustering
hb_result <- hierBAPS(snp_matrix, max.depth = 2, n.pops = floor(nrow(snp_matrix) / 5))
# Extract partitions
partitions <- hb_result$partitions
# Create a data frame with sample names and their corresponding clusters
cluster_df <- data.frame(sample = rownames(snp_matrix), cluster = partitions)

# IQ-TREE ML tree (optional) --------------------------------------------------
if (RUN_IQTREE) {
  cat("Inferring ML tree with IQ-TREE (GTR+G, 1000 ultrafast bootstraps)...\n")
  iq_input <- file.path(OUTDIR, "mt_masked_full.fasta")
  iq_prefix <- file.path(OUTDIR, "iqtree_run")
  cmd <- sprintf("iqtree -s %s -m GTR+G -bb 1000 -nt %s -pre %s",
                 shQuote(iq_input), IQTREE_THREADS, shQuote(iq_prefix))
  run(cmd)
  treefile <- paste0(iq_prefix, ".treefile")
  if (file.exists(treefile)) {
    file.copy(treefile, file.path(OUTDIR, "mt_ml_tree.nwk"), overwrite = TRUE)
    cat("ML tree written: ", file.path(OUTDIR, "mt_ml_tree.nwk"), "\n")
  } else {
    warning("IQ-TREE finished without treefile; check logs in OUTDIR.")
  }
}

# Using PEGAS

library(ape)
library(pegas)

# Read aligned sequences
dna <- read.dna(file.path(OUTDIR, "mt_consensus_all_aligned.fasta"), format = "fasta")

# Identify haplotypes
haps <- haplotype(dna)
haps
hap_freq <- haploFreq(haps)

hap_assign <- with(haps, setNames(hap, rownames(dna)))
hap_assign

net <- haploNet(haps)
plot(net, size=attr(net, "freq"), show.mutation=TRUE)

# Consolidate memberships -----------------------------------------------------
cat("Writing consolidated cluster membership table...\n")
members <- data.frame(sample = names(dna_masked), stringsAsFactors = FALSE)
# mclust
members$mclust <- NA_integer_
members$mclust[match(rownames(pc_use), members$sample)] <- mclust_assign
# fastbaps
fb_file <- file.path(OUTDIR, "fastbaps_clusters.csv")
if (file.exists(fb_file)) {
  fb <- read.csv(fb_file)
  members$fastbaps <- fb$fastbaps_cluster[match(members$sample, fb$sample)]
}
# rhierBAPS
hb_file <- file.path(OUTDIR, "hierBAPS_clusters.csv")
if (file.exists(hb_file)) {
  hb <- read.csv(hb_file)
  members$hierBAPS <- hb$hierBAPS_cluster[match(members$sample, hb$sample)]
}
write.csv(members, file = file.path(OUTDIR, "cluster_memberships_allmethods.csv"), row.names = FALSE)

cat("\nDone. Key outputs in:", normalizePath(OUTDIR), "\n")
cat("- Full-length masked FASTA: mt_masked_full.fasta\n")
cat("- Variant-only FASTA:       mt_variants_only.fasta\n")
cat("- mclust BIC plot:          mclust_BIC_plot.pdf\n")
cat("- Cluster tables:           mclust/fastbaps/hierBAPS CSVs + consolidated file\n")
cat("- (Optional) ML tree:       mt_ml_tree.nwk\n")