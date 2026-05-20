#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(circlize)
})

# ---- Arguments ----
option_list <- list(
  make_option(c("--snp"), type="character", help="SNP TSV: chr pos"),
  make_option(c("--indel"), type="character", help="INDEL TSV: chr pos"),
  make_option(c("--gene"), type="character", help="Gene TSV: chr start end"),
  make_option(c("--bin"), type="integer", default=50000, help="Bin size"),
  make_option(c("--out"), type="character", default="circos_density.pdf")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$snp) || is.null(opt$gene)) stop("Provide --snp, --indel and --gene files")
bin_size <- opt$bin

# ---- Load data ----
snps  <- read.table(opt$snp, header=FALSE, col.names=c("chr","pos"), stringsAsFactors=FALSE)
indels  <- read.table(opt$indel, header=FALSE, col.names=c("chr","pos"), stringsAsFactors=FALSE)
genes <- read.table(opt$gene, header=FALSE, col.names=c("chr","start","end"), stringsAsFactors=FALSE)

# ---- Chromosome lengths ----
chr_lengths <- data.frame(
  chr=c("X", "2L","2R","3L","3R", "4"),
  length=c(23542271, 23513712, 25286936, 28110227, 32079331, 1348131),
  stringsAsFactors = FALSE
)
chr_order <- chr_lengths$chr

# ---- Function to bin positions ----
make_bins <- function(chr, chr_len, bin_size) {
  starts <- seq(0, chr_len-1, by=bin_size)
  ends <- pmin(starts + bin_size, chr_len)
  data.frame(chr=chr, bin_start=starts, bin_end=ends)
}

print("Making SNP bins")
# ---- Bin SNPs and Genes per chromosome ----
snp_bins_list <- lapply(chr_order, function(chr) {
  chr_snps <- snps %>% filter(chr==!!chr)
  chr_len <- chr_lengths$length[chr_lengths$chr==chr]
  bins <- make_bins(chr, chr_len, bin_size)
  bins$count <- sapply(1:nrow(bins), function(i) {
    sum(chr_snps$pos >= bins$bin_start[i] & chr_snps$pos < bins$bin_end[i])
  })
  bins
})
snp_bins <- do.call(rbind, snp_bins_list)
print("Finished making SNP bins")

print("Making INDEL bins")
indel_bins_list <- lapply(chr_order, function(chr) {
  chr_indels <- indels %>% filter(chr==!!chr)
  chr_len <- chr_lengths$length[chr_lengths$chr==chr]
  bins <- make_bins(chr, chr_len, bin_size)
  bins$count <- sapply(1:nrow(bins), function(i) {
    sum(chr_indels$pos >= bins$bin_start[i] & chr_indels$pos < bins$bin_end[i])
  })
  bins
})
indel_bins <- do.call(rbind, indel_bins_list)
print("Finished making INDEL bins")

print("Making GENE bins")
gene_bins_list <- lapply(chr_order, function(chr) {
  chr_genes <- genes %>% filter(chr==!!chr)
  chr_len <- chr_lengths$length[chr_lengths$chr==chr]
  bins <- make_bins(chr, chr_len, bin_size)
  bins$count <- sapply(1:nrow(bins), function(i) {
    sum(chr_genes$start < bins$bin_end[i] & chr_genes$end > bins$bin_start[i])
  })
  bins
})
gene_bins <- do.call(rbind, gene_bins_list)
print("Finished making GENE bins")


# ---- Circos plot ----
svg(opt$out, width=8, height=8)
circos.clear()
circos.par(start.degree=90, gap.degree=4)

# Convert chr_lengths (chr + length) into BED-like format (chr, start, end)
bed_data <- data.frame(
  chr   = chr_lengths$chr,
  start = rep(0, nrow(chr_lengths)),        # start at 0
  end   = chr_lengths$length                # end at chromosome length
)

# Now initialize
circos.genomicInitialize(
  bed_data,
  sector.names = chr_order,
  plotType = c("axis", "labels"),
  tickLabelsStartFromZero = TRUE,
  axis.labels.cex = 0.4*par("cex"),
  labels.cex = 0.8*par("cex")
)


# Outer track: Gene density
gene_max <- max(gene_bins$count)
circos.track(ylim=c(0,gene_max), track.height=0.1, bg.border=NA,
             panel.fun=function(x,y){
               chr <- CELL_META$sector.index
               df <- gene_bins %>% filter(chr==!!chr & count>0)
               if(nrow(df)>0){
                 circos.genomicRect(df[,c("chr","bin_start","bin_end")],
                                    ybottom=0, ytop=df$count,
                                    col="darkgreen", border=NA)
               }
             })

# First inner track: SNP density
snp_max <- max(snp_bins$count)
circos.track(ylim=c(0,snp_max), track.height=0.1, bg.border=NA,
             panel.fun=function(x,y){
               chr <- CELL_META$sector.index
               df <- snp_bins %>% filter(chr==!!chr & count>0)
               if(nrow(df)>0){
                 circos.genomicRect(df[,c("chr","bin_start","bin_end")],
                                    ybottom=0, ytop=df$count,
                                    col="steelblue", border=NA)
               }
             })

# Second inner track: SNP density
indel_max <- max(indel_bins$count)
circos.track(ylim=c(0,indel_max), track.height=0.1, bg.border=NA,
             panel.fun=function(x,y){
               chr <- CELL_META$sector.index
               df <- indel_bins %>% filter(chr==!!chr & count>0)
               if(nrow(df)>0){
                 circos.genomicRect(df[,c("chr","bin_start","bin_end")],
                                    ybottom=0, ytop=df$count,
                                    col="indianred", border=NA)
               }
             })

dev.off()
message("Circos density plot saved to: ", opt$out)