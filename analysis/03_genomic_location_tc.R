#/*==========================================================================#*/
#' ## Figures 3A
#' 
#/*==========================================================================#*/

# add gene annotation using ChIPseeker
library(ChIPseeker)
library(biomaRt)
library(GenomicFeatures)

txdb <- makeTxDbFromBiomart(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl")

# - rename chromosome names in txdb to match NCBI naming scheme - #
seqlevels(txdb) <- c("NC_001133.9","NC_001134.8","NC_001135.5",
                     "NC_001136.10", "NC_001137.3", "NC_001138.5", 
                     "NC_001139.9", "NC_001140.6", "NC_001141.2", 
                     "NC_001142.9", "NC_001143.9", "NC_001144.5", 
                     "NC_001145.3", "NC_001146.8", "NC_001147.6", 
                     "NC_001148.4", "NC_001224.1")
seqlevels(txdb)
txdb


# select samples for annotation
samples <- c("Y1ng_carrier", "Y5ng_carrier", "Y10ng_carrier1", "BY4741")
tc_selected <- tc.grl[samples]

# change names
names <- c("SLIC 1 ng", "SLIC 5 ng", "SLIC 10 ng", "nAnTi 5 ug")
names(tc_selected) <- names

# use peak anno to annotate tag clusters
peakAnno_list <- lapply(tc_selected, function(x) annotatePeak(x, TxDb = txdb,  tssRegion = c(-500, 500), annoDb = "org.Sc.eg.db", sameStrand = TRUE, verbose = FALSE))

# set names
names(peakAnno_list) <- names

# plotting is limited within the package so I extract features in a list
feats.l <- lapply(peakAnno_list, function(x) x@annoStat)
names(feats.l) <- names

# add names to as a column to each sample
for (i in 1:length(names)) {
  feats.l[[i]]$sample <- rep(names[[i]], nrow(feats.l[[i]]))
}

# connect to 1 dataframe
feats.df <- do.call("rbind", feats.l)

# set levels for plotting (sample order and feature order)
feats.df$sample <- factor(feats.df$sample, levels = names[4:1])
features <- c("Promoter", "1st Exon", "Other Exon", "1st Intron", "Downstream (<=3kb)", "Distal Intergenic")
feats.df$Feature <- factor(feats.df$Feature, levels = features)

# set colors
library(RColorBrewer)
col1 <- brewer.pal(5, "Greys")[3:5]
col2 <- brewer.pal(8, "Blues")[8:6]
col <- c(col2, col1)

# new colour scheme
col <- c("#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
         "#0B877D", "#126872", "#031727")

# set feature factors
feature_factors <- c("Promoter", "Promoter (<=1kb)", "Promoter (1-3kb)", "5' UTR", "1st Exon", "Other Exon", 
                     "1st Intron", "Other Intron", "3' UTR", "Downstream (<=3kb)", "Distal Intergenic")
names(col) <- feature_factors
col_sel <- col[names(col) %in% unique(feats.df$Feature)]

# plot genomic features
library(ggplot2)
p <- ggplot(feats.df, aes(x = sample, y = Frequency, fill = Feature), alpha = 0.7) +
  geom_bar(stat = "identity", width = 0.75, col = "black", lwd = 0.125) +
  coord_flip() +
  scale_fill_manual("Features", values = col_sel) +
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Percentage", x = NULL)

pdf(file = "final_figures/genomicFeatures_sc_sel.pdf", height = 3, width = 5)
print(p)
dev.off()

#----- plot for all samples

# add alphabet letters to get the proper plotting order
names <- c("SLIC 1 ng", "SLIC 2 ng", "SLIC 5 ng", "SLIC 10 ng r1", "SLIC 10 ng r2", 
           "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi PCR", "nAnTi 5 ug")

# use peak anno to annotate tag clusters
peakAnno_list <- lapply(tc.grl, function(x) annotatePeak(x, TxDb = txdb,  tssRegion = c(-500, 500), annoDb = "org.Sc.eg.db", sameStrand = TRUE, verbose = FALSE))

# set names
names(peakAnno_list) <- names

# plotting is limited within the package so I extract features in a list
feats.l <- lapply(peakAnno_list, function(x) x@annoStat)
names(feats.l) <- names

# add names to as a column to each sample
for (i in 1:length(names)) {
  feats.l[[i]]$sample <- rep(names[[i]], nrow(feats.l[[i]]))
}

# connect to 1 dataframe
feats.df <- do.call("rbind", feats.l)

# set levels for plotting (sample order and feature order)
feats.df$sample <- factor(feats.df$sample, levels = names[length(names):1])
features <- c("Promoter", "1st Exon", "Other Exon", "1st Intron", "Downstream (<=3kb)", "Distal Intergenic")
feats.df$Feature <- factor(feats.df$Feature, levels = features)

# set colors
library(RColorBrewer)
col1 <- brewer.pal(5, "Greys")[3:5]
col2 <- brewer.pal(8, "Blues")[8:6]
col <- c(col2, col1)

# new colour scheme
col <- c("#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
         "#0B877D", "#126872", "#031727")

# set feature factors
feature_factors <- c("Promoter", "Promoter (<=1kb)", "Promoter (1-3kb)", "5' UTR", "1st Exon", "Other Exon", 
                     "1st Intron", "Other Intron", "3' UTR", "Downstream (<=3kb)", "Distal Intergenic")
names(col) <- feature_factors
col_sel <- col[names(col) %in% unique(feats.df$Feature)]

# plot genomic features
library(ggplot2)
p <- ggplot(feats.df, aes(x = sample, y = Frequency, fill = Feature), alpha = 0.7) +
  geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
  coord_flip() +
  scale_fill_manual("Features", values = col_sel) +
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Percentage", x = NULL)

pdf(file = "final_figures/genomicFeatures_sc_all.pdf", height = 3, width = 5)
print(p)
dev.off()