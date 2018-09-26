#/*==========================================================================#*/
#' ## Figure 5G - SOM expression profiling
#'
#/*==========================================================================#*/

# SOM - expression profiling
## 4_2 cluster
# create consensus clusters
aggregateTagClusters(myCAGEset, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
# 4_2 cluster
getExpressionProfiles(myCAGEset, what = "consensusClusters", tpmThreshold = 10, nrPassThreshold = 1, method = "som", xDim = 4, yDim = 2)
plotExpressionProfiles(myCAGEset, what = "consensusClusters")

# extract clusters
clusters <- c(paste0(0:3, "_", 0), paste0(0:3, "_", 1))
names(clusters) <- clusters

som_cons_clusters.l <- lapply(clusters, function(x) extractExpressionClass(myCAGEset, what = "consensusClusters", which = x))

# save SOM 4_2 classes
saveRDS(som_cons_clusters.l, "intermediate_data/SOM_consensus_clusters_4_2_l.RDS")

# convert from dataframe to GRanges object
som_cons_clusters.grl <- lapply(som_cons_clusters.l, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = TRUE, 
                                                                                          ignore.strand = FALSE,
                                                                                          seqnames.field = "chr",
                                                                                          start.field = "start",
                                                                                          end.field = "end",
                                                                                          strand.field = "strand"))

# save SOM classes as GRanges object 
saveRDS(som_cons_clusters.grl, "intermediate_data/SOM_consensus_clusters_4_2_grl.RDS")