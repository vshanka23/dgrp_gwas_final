#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(circlize)
})

# ---- Arguments ----
option_list <- list(
  make_option(c("--snv"),  type = "character", help = "SNV TSV: chr pos"),
  make_option(c("--gene"), type = "character", help = "Gene TSV: chr start end"),
  make_option(c("--inv"),  type = "character", help = "Inversion TSV: chr start end"),
  make_option(c("--bin"),  type = "integer", default = 50000,
              help = "Bin size for gene density track [default %default]"),
  make_option(c("--out"),  type = "character", default = "circos_variants.svg",
              help = "Output SVG file [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$snv) || is.null(opt$gene) || is.null(opt$inv)) {
  stop("Provide --snv, --gene, and --inv files")
}

if (!grepl("\\.svg$", opt$out, ignore.case = TRUE)) {
  opt$out <- paste0(opt$out, ".svg")
}

bin_size <- opt$bin

# ---- Load data ----
snvs <- read.table(
  opt$snv,
  header = FALSE,
  col.names = c("chr", "pos"),
  stringsAsFactors = FALSE
)

genes <- read.table(
  opt$gene,
  header = FALSE,
  col.names = c("chr", "start", "end"),
  stringsAsFactors = FALSE
)

inversions <- read.table(
  opt$inv,
  header = FALSE,
  col.names = c("chr", "start", "end"),
  stringsAsFactors = FALSE
)

# ---- Chromosome lengths ----
chr_lengths <- data.frame(
  chr = c("X", "2L", "2R", "3L", "3R", "4"),
  length = c(23542271, 23513712, 25286936, 28110227, 32079331, 1348131),
  stringsAsFactors = FALSE
)
chr_order <- chr_lengths$chr

# ---- Keep only recognized chromosomes ----
snvs <- snvs %>%
  filter(chr %in% chr_order) %>%
  mutate(pos = as.numeric(pos))

genes <- genes %>%
  filter(chr %in% chr_order)

inversions <- inversions %>%
  filter(chr %in% chr_order)

# ---- Drop out-of-range SNVs ----
snvs <- snvs %>%
  left_join(chr_lengths, by = "chr") %>%
  filter(pos >= 0 & pos <= length) %>%
  select(chr, pos)

# ---- Function to bin positions (genes only) ----
make_bins <- function(chr, chr_len, bin_size) {
  starts <- seq(0, chr_len - 1, by = bin_size)
  ends <- pmin(starts + bin_size, chr_len)
  data.frame(chr = chr, bin_start = starts, bin_end = ends)
}

# ---- Bin genes per chromosome ----
message("Making GENE bins")
gene_bins_list <- lapply(chr_order, function(chr) {
  chr_genes <- genes %>% filter(chr == !!chr)
  chr_len <- chr_lengths$length[chr_lengths$chr == chr]
  bins <- make_bins(chr, chr_len, bin_size)
  bins$count <- vapply(seq_len(nrow(bins)), function(i) {
    sum(chr_genes$start < bins$bin_end[i] & chr_genes$end > bins$bin_start[i])
  }, integer(1))
  bins
})
gene_bins <- do.call(rbind, gene_bins_list)
message("Finished making GENE bins")

# ---- Circos plot ----
svg(opt$out, width = 8, height = 8)
circos.clear()
circos.par(start.degree = 90, gap.degree = 4)

bed_data <- data.frame(
  chr   = chr_lengths$chr,
  start = rep(0, nrow(chr_lengths)),
  end   = chr_lengths$length
)

circos.genomicInitialize(
  bed_data,
  sector.names = chr_order,
  plotType = c("axis", "labels"),
  tickLabelsStartFromZero = TRUE,
  axis.labels.cex = 0.4 * par("cex"),
  labels.cex = 0.8 * par("cex")
)

# ---- Inversion track ----
circos.track(
  ylim = c(0, 5),
  track.height = 0.10,
  bg.border = NA,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    df <- inversions %>% filter(chr == !!chr)
    if (nrow(df) > 0) {
      df$row <- (seq_len(nrow(df)) - 1) %% 5
      circos.genomicRect(
        df[, c("chr", "start", "end")],
        ybottom = df$row + 0.1,
        ytop    = df$row + 0.8,
        col     = rgb(0.6, 0.1, 0.7, 0.5),
        border  = NA
      )
    }
  }
)

# ---- Gene density track ----
gene_max <- max(gene_bins$count)
if (!is.finite(gene_max) || gene_max < 1) gene_max <- 1

circos.track(
  ylim = c(0, gene_max),
  track.height = 0.10,
  bg.border = NA,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    df <- gene_bins %>% filter(chr == !!chr & count > 0)
    if (nrow(df) > 0) {
      circos.genomicRect(
        df[, c("chr", "bin_start", "bin_end")],
        ybottom = 0,
        ytop    = df$count,
        col     = "darkgreen",
        border  = NA
      )
    }
  }
)

# ---- Single SNV marker track ----
circos.track(
  ylim = c(0, 1),
  track.height = 0.06,
  bg.border = NA,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    df <- snvs %>% filter(chr == !!chr)
    if (nrow(df) > 0) {
      circos.points(
        x   = df$pos,
        y   = rep(0.5, nrow(df)),
        pch = 16,
        cex = 0.45,
        col = "red"
      )
    }
  }
)

dev.off()
circos.clear()
message("Circos plot saved to: ", opt$out)
