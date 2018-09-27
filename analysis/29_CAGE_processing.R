#/*==========================================================================#*/
#' ## CAGE data processing - starting from bam files
#'
#/*==========================================================================#*/

# Load S. cerevisiae genome - rename the chromosomes 
library(BSgenome.Scerevisiae.UCSC.sacCer3)
seqnames(BSgenome.Scerevisiae.UCSC.sacCer3) <- c("NC_001133.9","NC_001134.8","NC_001135.5",
                                                 "NC_001136.10", "NC_001137.3", "NC_001138.5", 
                                                 "NC_001139.9", "NC_001140.6", "NC_001141.2", 
                                                 "NC_001142.9", "NC_001143.9", "NC_001144.5", 
                                                 "NC_001145.3", "NC_001146.8", "NC_001147.6", 
                                                 "NC_001148.4", "NC_001224.1")
seqnames(BSgenome.Scerevisiae.UCSC.sacCer3)

# Importing all samples
# loading in data from low input test 4 - 10 ng, 25 ng, 50 ng and 100 ng (2 replicates)
inputDir1 <- "/Users/ncvetesic/Documents/projects/carrierCAGE/data/mapped.AR1AK"
paths1 <- list.files(inputDir1, full.names = T)
pathsToInputFiles1 <- paths1[grep(pattern = "*sorted.bam$", paths1)]

# loading in data from 5ug yeast CAGE (2 replicates for comparison) 
inputDir2 <- "/Users/ncvetesic/Documents/projects/carrierCAGE/data/mapped_8techRepYeast"
paths2 <- list.files(inputDir2, full.names = T)
pathsToInputFiles2 <- paths2[grep(pattern = "*sorted.bam$", paths2)][1:2]

# loading in data from LICAGE5 - 2 rep 1ng, 2 rep 2ng, 2 rep 5ng, 2 rep 10ng
inputDir3 <- "/Users/ncvetesic/Documents/projects/carrierCAGE/data/mapped.AVEE8"
paths3 <- list.files(inputDir3, full.names = T)
pathsToInputFiles3 <- paths3[grep(pattern = "*sorted.bam$", paths3)]

# loading in data from low input test 2 - PCR amplification test on nanti library (sample8)
inputDir4 <- "/Users/ncvetesic/Documents/projects/carrierCAGE/data/mapped.AN8VR"
paths4 <- list.files(inputDir4, full.names = T)
pathsToInputFiles4 <- paths4[grep(pattern = "*8.sorted.bam$", paths4)]

# unified path for all input files
pathsToInputFiles <- c(pathsToInputFiles1, pathsToInputFiles2, pathsToInputFiles3, pathsToInputFiles4)


# Create a new CAGE object for all samples
library(CAGEr)

LICAGE <- new("CAGEset", genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3", 
              inputFiles = pathsToInputFiles, inputFilesType = "bam", 
              sampleLabels = c("Y10ng_carrier1_t1", "Y10ng_carrier2_t1", 
                               "Y25ng_carrier1", "Y25ng_carrier2", 
                               "Y50ng_carrier1", "Y50ng_carrier2",
                               "Y100ng_carrier1", "Y100ng_carrier2",
                               "BY4741_rep1", "BY4741_rep2", 
                               "Y1ng_carrier1", "Y1ng_carrier2",
                               "Y2ng_carrier1", "Y2ng_carrier2",
                               "Y5ng_carrier1", "Y5ng_carrier2",
                               "Y10ng_carrier1_t2", "Y10ng_carrier2_t2",
                               "Ynanti_PCR"))
getCTSS(LICAGE)

# Quality control
## Calculate correlation prior to normalization

# Pearson correlation
setwd("results/")
corr.m <- plotCorrelation(LICAGE, samples = "all", method = "pearson")


## Normalization and correlation
setwd("results/")
librarySizes(LICAGE)
plotReverseCumulatives(LICAGE, fitInRange = c(100, 1000), onePlot = TRUE)

# powerLaw normalization
normalizeTagCount(LICAGE, method = "powerLaw", fitInRange = c(100, 1000), alpha = 1.14, T = 1000000)
plotReverseCumulatives(LICAGE, fitInRange = c(100, 1000), onePlot = TRUE, values = "normalized")

# correlation of normalized data
corr.norm <- plotCorrelation(LICAGE, samples = "all", method = "pearson", values = "normalized")

[!Pearson correlation, power-law data](CTSS_raw_values_pairwise_correlation_pearson.png)

## Merge replicates and reorder samples
mergeSamples(LICAGE, mergeIndex = c(4,4,6,6,7,7,8,8,10,10,1,1,2,2,3,3,5,5,9), mergedSampleLabels = c("Y1ng_carrier", "Y2ng_carrier", "Y5ng_carrier", "Y10ng_carrier1", "Y10ng_carrier2", "Y25ng_carrier", "Y50ng_carrier", "Y100ng_carrier", "Ynanti_PCR", "BY4741"))

# Pearson correlation
setwd("results/")
corr.m <- plotCorrelation(LICAGE, samples = "all", method = "pearson")

[!Pearson correlation, raw data](CTSS_raw_values_pairwise_correlation_pearson.png)

## Normalization of merged replicates and correlation
setwd("results/")
librarySizes(LICAGE)
plotReverseCumulatives(LICAGE, fitInRange = c(100, 1000), onePlot = TRUE)

# powerLaw normalization
normalizeTagCount(LICAGE, method = "powerLaw", fitInRange = c(100, 1000), alpha = 1.14, T = 1000000)
plotReverseCumulatives(LICAGE, fitInRange = c(100, 1000), onePlot = TRUE, values = "normalized")

# correlation of normalized data - Pearson
corr.norm <- plotCorrelation(LICAGE, samples = "all", method = "pearson", values = "normalized")

# correlation of normalized data - spearman
setwd("results/")
corr.norm <- plotCorrelation(LICAGE, samples = "all", 
                             method = "spearman", values = "normalized",
                             tagCountThreshold = 1, applyThresholdBoth = TRUE)

# correlation of normalized data - spearman
setwd("results/")
corr.norm <- plotCorrelation(LICAGE, samples = "all", 
                             method = "spearman", values = "normalized",
                             tagCountThreshold = 1, applyThresholdBoth = FALSE)

# Export signal to bedgraph
setwd("tracks/")
exportCTSStoBedGraph(LICAGE, values = "normalized", format = "bedGraph", oneFile = FALSE)


# CTSS clustering
clusterCTSS(object = LICAGE, threshold = 1, thresholdIsTpm = TRUE, nrPassThreshold = 1, 
            method = "distclu", maxDist = 20, removeSingletons = TRUE, keepSingletonsAbove = 5)


# Calculate and plot promoter width

# - calculating promoter width - #
setwd("results/")
cumulativeCTSSdistribution(LICAGE, clusters = "tagClusters")
quantilePositions(LICAGE, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(LICAGE, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)

# Save CAGE object for later use
saveRDS(LICAGE, file = "data/LICAGE.RDS")

# Extract tag clusters 
```{r, message = FALSE, warning = FALSE, eval = FALSE}
sampleLabels <- unname(sampleLabels(LICAGE))
tc.l <- lapply(sampleLabels, 
               function (x) tagClusters(LICAGE, sample = x, 
                                        returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9))
names(tc.l) <- sampleLabels

## Convert tag clusters list into a GRanges list
tc.grl <- GRangesList()
tc.grl <- lapply(tc.l, function(x) GRanges(seqnames = x$chr, 
                                           ranges = IRanges(start = x$start, end = x$end),
                                           strand = x$strand,
                                           nr_ctss = x$nr_ctss,
                                           dominant_ctss = x$dominant_ctss,
                                           tpm = x$tpm,
                                           tpm.dominant_ctss = x$tpm.dominant_ctss,
                                           q_0.1 = x$q_0.1,
                                           q_0.9 = x$q_0.9,
                                           interquantile_width = x$interquantile_width,
                                           seqlengths = seqlengths(BSgenome.Scerevisiae.UCSC.sacCer3)))

# add seqinfo information                
for (i in 1:length(tc.grl)) {
  seqinfo(tc.grl[[i]]) <- seqinfo(BSgenome.Scerevisiae.UCSC.sacCer3) 
}

# save tc.grl as RDS object
saveRDS(object = tc.grl, file = "data/tagClusters_list.RDS")

## Dominant TSS position centered GRanges object
domTSS.grl <- list()

for (i in 1:length(tc.grl)) {
  domTSS.grl[[i]] <- GRanges(seqnames = seqnames(tc.grl[[i]]),
                             ranges = IRanges(start = tc.grl[[i]]$dominant_ctss, 
                                              end = tc.grl[[i]]$dominant_ctss),
                             strand = strand(tc.grl[[i]]),
                             nr_ctss = tc.grl[[i]]$nr_ctss,
                             dominant_ctss = tc.grl[[i]]$dominant_ctss,
                             tpm = tc.grl[[i]]$tpm, 
                             tpm.dominant_ctss = tc.grl[[i]]$tpm.dominant_ctss, 
                             q_0.1 = tc.grl[[i]]$q_0.1,
                             q_0.9 = tc.grl[[i]]$q_0.9,
                             interquantile_width = tc.grl[[i]]$interquantile_width,
                             tc_start = start(tc.grl[[i]]),
                             tc_end = end(tc.grl[[i]]))
  
  seqlevels(domTSS.grl[[i]]) <- seqlevels(tc.grl[[i]])
  #attach chromosome sequence lengths
  seqlengths(domTSS.grl[[i]]) <- seqlengths(tc.grl[[i]])
  
  #attach genome information
  genome(domTSS.grl[[i]]) <- genome(BSgenome.Scerevisiae.UCSC.sacCer3)
}

# attach sample names
names(domTSS.grl) <- names(tc.grl)

# save dominant TSS cluster GRanges object for later
saveRDS(domTSS.grl, file ="data/domTSS_grl.RDS")
