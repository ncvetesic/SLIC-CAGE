#/*==========================================================================#*/
#' ## Supplemental Figure S14 - correlation CTSS and consensus cluster level
#'
#/*==========================================================================#*/

#------------ Correlation plots - selected samples at CTSS level -----------#

tag.count <- K562_CAGE@normalizedTpmMatrix

# I should filter for each compared pair if both do not have more than 1 TPM in that sample - make 1 plot and then create a plot from several samples

samples <- "K562_nanoCAGE_XL"
control <- "K562_CAGE"

filtr.data.l <- list()
filtr.data.log.l <- list()
cor.l <- list()

#plot log scale but calculate correlation non non-log scale
for (i in 1: length(samples)) {
  filtr.data.l[[i]] <- tag.count[ , c(samples[[i]], control)]
  filtr.data.l[[i]] <- filtr.data.l[[i]][(filtr.data.l[[i]][, 1] >= 1 | filtr.data.l[[i]][, 2] >= 1), ]
  cor.l[[i]] <- round(cor(filtr.data.l[[i]]), digits = 3)
  filtr.data.log.l[[i]] <- log10(filtr.data.l[[i]] + 1)
}
names(filtr.data.l) <- samples
names(filtr.data.log.l) <- samples

# plots, for manual layout in inkscape
library(ggplot2)
corrPlot <- function(data, cor, xlab, ylab) {
  p <- ggplot(data, aes(data[,2], y = data[,1])) +
    geom_hex(bins = 100, show.legend = FALSE) +
    scale_fill_gradient(low = "gray24", high = "gray56") +
    theme(text = element_text(size = 20,
                              family = "Helvetica"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.subtitle = element_text(size = 24, hjust = -0.125, vjust = 1, face="bold", color="black"),
          axis.text = element_text(size = 18, family = "Helvetica", color = "black"),
          legend.title = element_text(size = 16, family = "Helvetica", colour = "black"),
          legend.text = element_text(size = 14, family = "Helvetica", colour = "black")) +
    xlab(xlab) +
    ylab(ylab) +
    xlim(c(0, 4.25)) +
    ylim(c(0, 4.25)) +
    annotate("text", x = 0, y = 4, size = 6, hjust = 0, label = paste("R = ", cor))  +
    coord_fixed(ratio = 1)
  
  return(p)
}

labels_y <- "K562_nanoCAGE"
labels_x <- "K562_CAGE"
samples <- "K562_nanoCAGE"

p.l <- list()

for(i in 1:length(samples)) {
  pdf(file = paste0("../final_figures/", samples[[i]], ".pdf"), width = 7, height = 7)
  p <- corrPlot(data = filtr.data.log.l[[i]], xlab = labels_x[[i]], ylab = labels_y[[i]], 
                cor = cor.l[[i]][1,2])
  print(p)
  p.l[[i]] <- p
  dev.off()
}


pdf("results/nanoCAGE_XL.pdf", height = 6, width = 6)
print(p.l[[1]])
dev.off()


#----------- Correlation at consensus cluster level -----------# 
cons_tpm.df <- as.data.frame(K562_CAGE@consensusClustersTpmMatrix)

cor <- round(cor(cons_tpm.df), digits = 3)

# convert to log scale
cons_tpm_log.df <-  log10(cons_tpm.df + 1)

# plots, for manual layout in inkscape
library(ggplot2)
p <- ggplot(cons_tpm_log.df, aes(x = cons_tpm_log.df$K562_CAGE, y = cons_tpm_log.df$K562_nanoCAGE_XL)) +
  geom_hex(bins = 100, show.legend = FALSE) +
  scale_fill_gradient(low = "gray24", high = "gray56") +
  theme(text = element_text(size = 20,
                            family = "Helvetica"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 18, family = "Helvetica", color = "black"),
        legend.title = element_text(size = 16, family = "Helvetica", colour = "black"),
        legend.text = element_text(size = 14, family = "Helvetica", colour = "black")) +
  xlim(c(0, 4.25)) +
  ylim(c(0, 4.25)) +
  xlab("K562 CAGE") +
  ylab("K562 nanoCAGE XL") +
  annotate("text", x = 0, y = 4, size = 6, hjust = 0, label = paste("R = ",  cor[2]))  +
  coord_fixed(ratio = 1)

pdf("results/nanoCAGE_XL_consensus.pdf", height = 6, width = 6)
print(p)
dev.off()
