#/*==========================================================================#*/
#' ## Supplemental Figure S18 - correlation CTSS and consensus clusters
#'
#/*==========================================================================#*/
## Plot correlation of PGC E11.5 replicates
tag.count <- CAGE_dedup@normalizedTpmMatrix

# Filter each pair to have >=1 TPM per CTSS 
samples <- c("E11_5_rep1_dedup")
control <- "E11_5_rep1"

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
    xlim(c(0, 4.5)) +
    ylim(c(0, 4.5)) +
    annotate("text", x = 0, y = 4.5, size = 6, hjust = 0, label = paste("R = ", cor))  +
    coord_fixed(ratio = 1)
  
  return(p)
}

labels_y <- "E11.5 rep 2"
labels_x <- "E11.5 rep 2 deduplicated"
samples <- "PGC E11.5 dedup"


for(i in 1:length(samples)) {
  p <- corrPlot(data = filtr.data.log.l[[i]], xlab = labels_x[[i]], ylab = labels_y[[i]], 
                cor = cor.l[[i]][2])
  pdf(file = paste0("results/replicate/dedup/E11_5_vs_rep_", samples[[i]], ".pdf"), width = 6, height = 6)
  print(p)
  dev.off()
}


## Consensus clusters 
# aggregate the clusters across the samples (raw tag number): remember these are raw reads in Tpm slot
aggregateTagClusters(CAGE_dedup, 
                     tpmThreshold = 5, 
                     qLow = 0.1, 
                     qUp = 0.9, 
                     maxDist = 100)


## Correlation at consensus cluster level PGC E11.5 and replicate
cons_tpm.df <- as.data.frame(CAGE_dedup@consensusClustersTpmMatrix)

cor <- round(cor(cons_tpm.df ), digits = 3)

# convert to log scale
cons_tpm_log.df <-  log10(cons_tpm.df + 1)

# plots, for manual layout in inkscape
library(ggplot2)
p <- ggplot(cons_tpm_log.df, aes(x = cons_tpm_log.df$E11_5_rep1, y = cons_tpm_log.df$E11_5_rep1_dedup)) +
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
  xlab("PGC E11.5 SLIC rep1") +
  ylab("PGC E11.5 SLIC rep2") +
  annotate("text", x = 0, y = 4, size = 6, hjust = 0, label = paste("R = ",  cor[2]))  +
  coord_fixed(ratio = 1)

pdf("results/replicate/dedup/E11_5_replicate_dedup_cons_clusters.pdf", height = 6, width = 6)
print(p)
dev.off()
