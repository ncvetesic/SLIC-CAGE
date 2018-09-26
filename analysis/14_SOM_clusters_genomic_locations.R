#/*==========================================================================#*/
#' ## Figure 5H - SOM clusters, genomic locations
#'
#/*==========================================================================#*/

# load som_cons_clusters
som_cons_clusters_4_1.grl <- readRDS("intermediate_data/SOM_consensus_clusters_4_1_grl.RDS")

# annotate each class with gene names and genomic locations
# create gene annotation using ChIPseeker
library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


# use peak anno from ChipSeeker to annotate tag clusters
peakAnno.l <- lapply(som_cons_clusters_4_1.grl, function(x) annotatePeak(x, TxDb = txdb,  annoDb = NULL, sameStrand = TRUE, verbose = FALSE))

# convert to GRanges 
som_cons_clusters_4_1_anno.grl <- lapply(peakAnno.l, as.GRanges)

# save annotated GRanges object 
saveRDS(som_cons_clusters_4_1_anno.grl, "intermediate_data/som_cons_clusters_4_1_anno_grl.RDS")

# plot annotated feature
png(file = "results/annoFeat_soms_4_1.png", height = 400, width = 800)
plotAnnoBar(peakAnno.l, title = NULL)
dev.off()

# plotting is limited within the package so I extract genomic locations
feats.l <- lapply(peakAnno.l, function(x) return(x@annoStat))

# set class
for (i in 1:length(feats.l)) {
  feats.l[[i]]$Frequency <- round(as.numeric(feats.l[[i]]$Frequency), digits = 2)
}

# colour scheme
col <- c("#FF9F27","#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
         "#0B877D", "#126872", "#031727")

# set feature factors
feature_factors <- c("Promoter", "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR", "1st Exon", "Other Exon", "1st Intron", "Other Intron", "3' UTR", "Downstream (<=3kb)", "Distal Intergenic")

names(col) <- feature_factors
col_sel <- col[names(col) %in% unique(unlist(lapply(feats.l, function(x) x$Feature)))]

# convert the list to a dataframe
feats.df <- do.call(rbind, feats.l)

# add column of sample names
feats.df$sample <- gsub("\\..*", "", rownames(feats.df))

# order by % of each feature
order_lev <- unique(feats.df$Feature[order(feats.df$Frequency, decreasing = FALSE)])
feats.df$Feature <- factor(feats.df$Feature, levels = order_lev)

# set levels for samples
feats.df$sample <- factor(feats.df$sample, levels = rev(names(som_cons_clusters_4_1_anno.grl)))

# plot genomic features
library(ggplot2)
p <- ggplot(feats.df, aes(x = sample, y = Frequency, fill = Feature), alpha = 0.7) +
  geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
  coord_flip() +
  scale_fill_manual("Features", values = col_sel) +
  theme_bw() +
  theme(text = element_text(size = 12, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Percentage", x = NULL) +
  guides(fill = guide_legend(reverse = TRUE))

pdf(file = "results/genomicFeatures_E11_5_mESC_soms_4_1.pdf", height = 3, width = 6)
print(p)
dev.off()