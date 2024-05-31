make_figures <- function(filename, output_qq_png, output_man_png){

library(qqman)
res3 <- read.table(filename, header=T)
new3<-data.frame(res3$rs, res3$chr, res3$ps, res3$p_wald)
names(new3)<-c("SNP", "CHR", "BP", "P")
#Added this line after EH and KC found a particular scenario where NAs in p_wald column affects the Max P math operation from Spencer's bug fix.
new3 <- new3[!is.na(new3$P), ]
#Moving chromosome X to the front as chromosome 1
new3$CHR <- new3$CHR+1
new3$CHR[new3$CHR==8]<-1
#Removing chromsome Y and MT
new3 <- new3[new3$CHR<7,]
#png(output_qq_png)
#qq(new3$P)
#dev.off()

tiff(output_qq_png, width = 6, height = 6, units = 'in', res = 300)
qq(new3$P)
dev.off()

# png(output_man_png)
# #Spencer Hatfield identified some associations where some p-values were < 1E-10. These were being omitted because of the hard coded ylim in the plot. His solution is below:
# max_p <- max(-log10(new3$P))
# ylim_max <- ifelse(max_p > 10, max_p, 10)
# manhattan(new3, main = "Manhattan Plot", ylim = c(0, ylim_max), cex = 0.6, cex.axis = 0.9, col = c("purple3", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c("X","2L","2R","3L","3R","4"))
# dev.off()

#setEPS()
#postscript(paste(output_man_png,".eps",sep=""))
#max_p <- max(-log10(new3$P))
#ylim_max <- ifelse(max_p > 10, max_p, 10)
#manhattan(new3, main = "Manhattan Plot", ylim = c(0, ylim_max), cex = 0.6, cex.axis = 0.9, col = c("purple3", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c("X","2L","2R","3L","3R","4"))
#dev.off()

tiff(output_man_png, width = 6, height = 6, units = 'in', res = 300)
max_p <- max(-log10(new3$P))
ylim_max <- ifelse(max_p > 10, max_p, 10)
manhattan(new3, main = "Manhattan Plot", ylim = c(0, ylim_max), cex = 0.6, cex.axis = 0.9, col = c("purple3", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c("X","2L","2R","3L","3R","4"))
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