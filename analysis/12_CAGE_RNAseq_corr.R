#/*==========================================================================#*/
#' ## Figure 5A - correlation of CAGE and RNA-seq data
#'
#/*==========================================================================#*/

# -------- Plot correlation of mESC and PGC E11.5 CAGE data -------- # 

tag.count <- myCAGEset@normalizedTpmMatrix

# filter for each compared pair to have more than 1TPM

samples <- c("E11_5_l1")
control <- "mESC_E14"

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

labels_y <- "E11.5 SLIC-CAGE"
labels_x <- "E14 nAnTi"
samples <- "PGC E11.5"


for(i in 1:length(samples)) {
  p <- corrPlot(data = filtr.data.log.l[[i]], xlab = labels_x[[i]], ylab = labels_y[[i]], 
                cor = cor.l[[i]][2])
  pdf(file = paste0("results/E11_5_vs_E14_nanti", samples[[i]], ".pdf"), width = 6, height = 6)
  print(p)
  dev.off()
}


# -------- Correlation at consensus cluster level mESC and PGC E11.5 (CAGE data) -------- #
cons_tpm.df <- as.data.frame(myCAGEset@consensusClustersTpmMatrix)

cor <- round(cor(cons_tpm.df), digits = 3)

# convert to log scale
cons_tpm_log.df <-  log10(cons_tpm.df + 1)

# plots, for manual layout in inkscape
library(ggplot2)
p <- ggplot(cons_tpm_log.df, aes(x = cons_tpm_log.df$mESC_E14, y = cons_tpm_log.df$E11_5_l1)) +
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
  xlab("mESC E14 nAnTi") +
  ylab("PGC E11.5 SLIC") +
  annotate("text", x = 0, y = 4, size = 6, hjust = 0, label = paste("R = ",  cor[2]))  +
  coord_fixed(ratio = 1)

pdf("results/E11_5_vs_E14_nanti_cons_clusters.pdf", height = 6, width = 6)
print(p)
dev.off()

# -------- Combine count matrices E11_5 and RNA-seq, normalize using Deseq and plot correlation -------- #
library(DESeq2)
cm_all <- data.frame("gene_id" = cm_E11_5_all_cage$gene_id, "CAGE" = cm_E11_5_all_cage$counts, "RNA_seq" = cm_E11_5_all_rnaseq$counts)

rownames(cm_all) <- cm_all$gene_id
cm_all<- cm_all[, -1]

countMatrix <- as.matrix.data.frame(cm_all)

samples <- colnames(countMatrix)
info.df <- data.frame(condition = factor(x = c("CAGE", "RNA_seq"),
                                         levels = c("CAGE", "RNA_seq")),
                      row.names = samples)

# create DESeqDataSet object - stores count matrix and design
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = info.df,
                              design = ~ condition)

# differential expression analysis
## analysis
dds <- DESeq(dds)

## results
res <- results(dds)

## results summary
summary(res)

## extract normalized data
dds <- estimateSizeFactors(dds)
normCounts <- counts(dds, normalized=TRUE)

# filter to have at least 1 tpm per sample
normCounts_filtr <- normCounts[rowSums(normCounts) >= 1, ]

# plot using gpplot2
cor.l <- list()
cor.l[[1]] <- round(cor(normCounts_filtr, method = "pearson"), digits = 3)[[2]]
cor.l[[2]] <- round(cor(normCounts_filtr, method = "spearman"), digits = 3)[[2]]


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
    xlim(c(0, 6)) +
    ylim(c(0, 6)) +
    annotate("text", x = 0, y = 6, size = 6, hjust = 0, label = paste("R = ", cor))  +
    coord_fixed(ratio = 1)
  
  return(p)
}

labels_y <- "E11.5 SLIC-CAGE"
labels_x <- "E11.5 RNA-seq"
samples <- "PGC E11.5"

p.l <- list()


# convert to log scale
normCounts_filtr <- log10(normCounts_filtr + 1)
normCounts_filtr.df <- as.data.frame(normCounts_filtr)

for(i in 1:length(samples)) {
  pdf(file = paste0("results/E11_5_RNAseq_vs_CAGE_", samples[[i]], ".pdf"), width = 6, height = 6)
  p <- corrPlot(data = normCounts_filtr.df, xlab = labels_x[[i]], ylab = labels_y[[i]], 
                cor = cor.l[[2]])
  print(p)
  p.l[[i]] <- p
  dev.off()
}

