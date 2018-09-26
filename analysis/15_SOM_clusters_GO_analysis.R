#/*==========================================================================#*/
#' ## Figure 5I - GO analysis
#'
#/*==========================================================================#*/

# load SOM annotated GRanges
som_cons_clusters_anno.grl <- readRDS("intermediate_data/som_cons_clusters_anno_grl.RDS")

# # load libraries  
library(org.Mm.eg.db)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# genes 
gene <- genes(txdb)


# GO enrichment using cluster profiler
library(clusterProfiler)

# create gene universe - all genes from SOM classes

geneUniverse <- unlist(lapply(som_cons_clusters_anno.grl, function(x) return(x$geneId)))
samples <- names(som_cons_clusters_anno.grl)

geneSample.l <- lapply(som_cons_clusters_anno.grl, function(x) (x$geneId))
names(geneSample.l) <- samples


ego_MF.l <- lapply(geneSample.l, function(x) 
  enrichGO(gene = x,
           universe      = geneUniverse,
           OrgDb         = org.Mm.eg.db,
           ont           = "MF",
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = FALSE,
           keyType = "ENTREZID"))

ego_BP.l <- lapply(geneSample.l, function(x) 
  enrichGO(gene = x,
           universe      = geneUniverse,
           OrgDb         = org.Mm.eg.db,
           ont           = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = FALSE,
           keyType = "ENTREZID"))

ego_CC.l <- lapply(geneSample.l, function(x) 
  enrichGO(gene = x,
           universe      = geneUniverse,
           OrgDb         = org.Mm.eg.db,
           ont           = "CC",
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = FALSE,
           keyType = "ENTREZID"))

# print output in pdf
for (i in names(ego_BP.l)) {
  
  p <-  dotplot(ego_MF.l[[i]])
  pdf(paste0("results/SOM_4_2/GO_enrichment/MF/SOM_4_2_GO_MF_", i, ".pdf"), height = 4, width = 12)
  print(p)
  dev.off()
  
  p <-  dotplot(ego_BP.l[[i]])
  pdf(paste0("results/SOM_4_2/GO_enrichment/BP/SOM_4_2_GO_BP_", i, ".pdf"), height = 3.5, width = 7)
  print(p)
  dev.off()
  
  
  p <-  dotplot(ego_CC.l[[i]])
  pdf(paste0("results/SOM_4_2/GO_enrichment/CC/SOM_4_2_GO_CC_", i, ".pdf"), height = 4, width = 12)
  print(p)
  dev.off()
}

# save results in txt files
# CC mESC add info on NPLB class
for (i in 1:length(ego_CC.l)) {
  ego_CC.l[[i]] <- as.data.frame(ego_CC.l[[i]])
  ego_CC.l[[i]]$NPLB_class = rep(names(ego_CC.l[i]), times = nrow(ego_CC.l[[i]]))
}

write.csv(do.call(rbind, ego_CC.l), file = "results/SOM_4_2/GO_enrichment/BP/SOM_4_2_GO_BP.csv")

# BP add info on NPLB class
for (i in 1:length(ego_BP.l)) {
  ego_BP.l[[i]] <- as.data.frame(ego_BP.l[[i]])
  ego_BP.l[[i]]$NPLB_class = rep(names(ego_BP.l[i]), times = nrow(ego_BP.l[[i]]))
}

write.csv(do.call(rbind, ego_BP.l), file = "results/SOM_4_2/GO_enrichment/MF/SOM_4_2_GO_MF.csv")    
# MF add info on NPLB class
for (i in 1:length(ego_MF.l)) {
  ego_MF.l[[i]] <- as.data.frame(ego_MF.l[[i]])
  ego_MF.l[[i]]$NPLB_class = rep(names(ego_MF.l[i]), times = nrow(ego_MF.l[[i]]))
}

write.csv(do.call(rbind, ego_MF.l), file = "results/SOM_4_2/GO_enrichment/MF/SOM_4_2_GO_CC.csv")  