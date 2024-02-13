adjust_fdr <- function(input_file, output_file, sorted, sorted_filt) {
	if (missing(input_file) || missing(output_file)) {
    stop("Usage: adjust_fdr(input_file, output_file, sorted, sorted_filt)")
}
data <- read.table(input_file,header=TRUE)
wald_FDR <- p.adjust(data$p_wald,method="fdr")
res <- cbind(data,wald_FDR)
write.table(res, file=output_file,row.names=FALSE,col.names=TRUE,sep="\t")
data_sorted <- data[order(data$p_wald),]
data_sorted_filt <- data_sorted[data_sorted$p_wald < 0.001, ]
write.table(data_sorted, file=sorted,row.names=FALSE,col.names=TRUE,sep="\t")
write.table(data_sorted_filt, file=sorted_filt,row.names=FALSE,col.names=TRUE,sep="\t")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript adjust_fdr.R <association file from WALD test> <association file with FDR adjustment> <sorted data> <sorted filtered data>")
}

input_file <- args[1]
output_file <- args[2]
sorted <- args[3]
sorted_filt <- args[4]

adjust_fdr(input_file, output_file, sorted, sorted_filt)