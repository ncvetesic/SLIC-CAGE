#/*==========================================================================#*/
#' ## Supplemental Figure S2A-C - correlation matrices, CTSS level
#'
#/*==========================================================================#*/

## Plot all correlations for supplementary - corrplot
library(corrplot)

# extract TPM values normalized for all samples
tag.count <- LICAGE@normalizedTpmMatrix

# filter so at least 1 sample has >=1 TPM
tag.count.filtr <- tag.count[rowSums(tag.count) >= 1, ]

#calculate correlations
corr_mat <- cor(as.matrix(tag.count.filtr), method = "pearson")

names <- c("SLIC 1 ng", "SLIC 2 ng", "SLIC 5 ng", 
           "SLIC 10 ng r1", "SLIC 10 ng r2",
           "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", 
           "nAnTi PCR", 
           "nAnTi 5 \u03bcg r2")
rownames(corr_mat) <- names
colnames(corr_mat) <- names

# visualize correlations
library(corrplot)
quartz(type = "pdf", file = "final_figures/SLIC_CAGE_corrPlot_All.pdf", width = 9, height = 7, family = "Helvetica", dpi = 300)
corrplot(corr_mat, method = "color", 
         outline = T, 
         addgrid.col = "darkgray",
         tl.col = "black", 
         tl.cex = 1.5,
         tl.srt = 45,
         tl.pos = "ld",
         addCoef.col = "white",
         number.digits = 2,
         number.cex = 1.25,
         cl.lim = c(0, 1),
         cl.cex = 1.25,
         is.corr = TRUE,
         type = "lower")
dev.off()
