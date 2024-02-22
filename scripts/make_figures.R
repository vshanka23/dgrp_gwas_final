make_figures <- function(filename, output_qq_png, output_man_png){

library(qqman)
res3 <- read.table(filename, header=T)
new3<-data.frame(res3$rs, res3$chr, res3$ps, res3$p_wald)
names(new3)<-c("SNP", "CHR", "BP", "P")
png(output_qq_png)
qq(new3$P)
dev.off()
png(output_man_png)
#Spencer Hatfield identified some associations where some p-values were < 1E-10. These were being omitted because of the hard coded ylim in the plot. His solution is below:
max_p <- max(-log10(new3$P))
ylim_max <- ifelse(max_p > 10, max_p, 10)
manhattan(new3, main = "Manhattan Plot", ylim = c(0, ylim_max), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c("2L","2R","3L","3R","4","MT","X","Y"))
dev.off()
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
    stop("Usage: Rscript <Input_WALD_association_test_file.txt from PLINK/GCTA/GEMMA> <QQ_plot_name.png> <Manhattan_plot_name.png>")
}
filename <- args[1]
output_qq_png <- args[2]
output_man_png <- args[3]

make_figures(filename, output_qq_png, output_man_png)