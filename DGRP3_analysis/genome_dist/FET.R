#!/usr/bin/env Rscript

## =========================
## Parameters
## =========================

MT_LENGTH <- 19524
AT_START  <- 14900
AT_END    <- 19524

L_AT   <- AT_END - AT_START + 1
L_REST <- MT_LENGTH - L_AT

## =========================
## Helper function
## =========================

run_fet <- function(tsv_file, label) {

  df <- read.table(
    tsv_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    col.names = c("CHR", "POS")
  )

  ## Keep MT only (defensive)
  df <- df[df$CHR == "MT", ]

  in_AT <- df$POS >= AT_START & df$POS <= AT_END

  a <- sum(in_AT)      # variants in AT-rich
  b <- sum(!in_AT)     # variants outside

  ## FIXED contingency table
  contingency <- matrix(
    c(a,
      L_AT   - a,
      b,
      L_REST - b),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Region = c("AT_rich", "Rest"),
      Status = c("Variant", "No_Variant")
    )
  )

  test <- fisher.test(contingency)

  list(
    label = label,
    counts = contingency,
    odds_ratio = test$estimate,
    p_value = test$p.value
  )
}

## =========================
## Main
## =========================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: fet_mt_at_enrichment.R <snps.tsv> <indels.tsv>")
}

snps_file   <- args[1]
indels_file <- args[2]

snps_res   <- run_fet(snps_file,   "SNPs")
indels_res <- run_fet(indels_file, "INDELs")

## =========================
## Output
## =========================

print_result <- function(res) {
  cat("\n=============================\n")
  cat(res$label, "\n")
  cat("=============================\n")
  print(res$counts)
  cat("\nOdds ratio:", unname(res$odds_ratio), "\n")
  cat("Fisher's exact p-value:", res$p_value, "\n")
}

print_result(snps_res)
print_result(indels_res)