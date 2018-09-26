#/*===================================================================================#*/
#' ## Supplemental Figure S3 - simulation of low complexity and IQ-width distribution
#'
#/*===================================================================================#*/

# Using replicate one from nAnT-iCAGE yeast, downsampled using samtools view -s 2.2 -o BY4741_1_s2.2.sorted.bam BY4741_1.sorted.bam
# Fraction .2 is keep 20 % of the reads, while first integer number is seed

# Load S. cerevisiae genome - rename the chromosomes to NCBI
library(BSgenome.Scerevisiae.UCSC.sacCer3)
seqnames(BSgenome.Scerevisiae.UCSC.sacCer3) <- c("NC_001133.9","NC_001134.8","NC_001135.5",
                                                 "NC_001136.10", "NC_001137.3", "NC_001138.5", 
                                                 "NC_001139.9", "NC_001140.6", "NC_001141.2", 
                                                 "NC_001142.9", "NC_001143.9", "NC_001144.5", 
                                                 "NC_001145.3", "NC_001146.8", "NC_001147.6", 
                                                 "NC_001148.4", "NC_001224.1")
seqnames(BSgenome.Scerevisiae.UCSC.sacCer3)


# Import data - BAM files
# loading in data yeast CAGE replicate 1 and subsampled - 80 %, 60 %, 40 %, 20 % of the reads
inputDir <- "/Users/ncvetesic/Documents/projects/carrierCAGE/data/mapped_8techRepYeast"
paths <- list.files(inputDir, full.names = T)
pathsToInputFiles <- paths[grep(pattern = "BY4741_1.*sorted.bam$", paths)]

# reorder according to input amount / size
pathsToInputFiles <- pathsToInputFiles[c(6, 1:5)]

# Create a new CAGE object for all samples
library(CAGEr)

myCAGEset <- new("CAGEset", genomeName = "BSgenome.Scerevisiae.UCSC.sacCer3", 
                 inputFiles = pathsToInputFiles, inputFilesType = "bam", 
                 sampleLabels = c("nAnTi_100", "nAnTi_80", "nAnTi_60", "nAnTi_40", "nAnTi_20", "nAnTi_10"))
getCTSS(myCAGEset)

# Quality control
## Calculate correlation prior to normalization
# Pearson correlation
corr.m <- plotCorrelation(myCAGEset, samples = "all", method = "pearson")


## Normalization and correlation

librarySizes(myCAGEset)
plotReverseCumulatives(myCAGEset, fitInRange = c(5, 1000), onePlot = TRUE)

# powerLaw normalization
normalizeTagCount(myCAGEset, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.18, T = 1000000)
plotReverseCumulatives(myCAGEset, fitInRange = c(5, 1000), onePlot = TRUE, values = "normalized")

# correlation of normalized data
corr.norm <- plotCorrelation(myCAGEset, samples = "all", method = "pearson", values = "normalized")

## CTSS clustering
clusterCTSS(object = myCAGEset, threshold = 1, thresholdIsTpm = TRUE, nrPassThreshold = 1, 
            method = "distclu", maxDist = 20, removeSingletons = TRUE, keepSingletonsAbove = 5)

## Calculate and plot promoter width
# - calculating promoter width - #
cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters")
quantilePositions(myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(myCAGEset, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)

# Save CAGE object for later use
saveRDS(myCAGEset, file = "subsampled_nAnTi.RDS")

# Methods Paper figures
## Correlation plots - selected samples
tag.count <- myCAGEset@normalizedTpmMatrix

# Filter each pair to 1 or more than 1 TPM per CTSS in a pair

samples <- sampleLabels(myCAGEset)[2:6]
control <- "nAnTi_100"

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
corrPlot <- function(data, cor, xlab, ylab, sub) {
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
    annotate("text", x = 0, y = 4, size = 6, hjust = 0, label = paste("R = ", cor)) +
    labs(subtitle = sub) +
    coord_fixed(ratio = 1)
  
  return(p)
}

labels_y <- c("nAnT-iCAGE 80 %", "nAnT-iCAGE 60 %", "nAnT-iCAGE 40 %", "nAnT-iCAGE 20 %", "nAnT-iCAGE 10 %")
labels_x <- c(rep("nAnT-iCAGE 100 %", times = 5))
samples <- c("nAnT-iCAGE_80", "nAnT-iCAGE_60", "nAnT-iCAGE_40", "nAnT-iCAGE_20", "nAnT-iCAGE_10")
subtitles <- c("a\n", "b\n", "c\n", "d\n", "e\n")

p.l <- list()

for(i in 1:length(samples)) {
  pdf(file = paste0(samples[[i]], ".pdf"), width = 7, height = 7)
  p <- corrPlot(data = filtr.data.log.l[[i]], xlab = labels_x[[i]], ylab = labels_y[[i]], 
                cor = cor.l[[i]][1,2], sub = subtitles[[i]])
  print(p)
  p.l[[i]] <- p
  dev.off()
}


# plots, set up layout in R altogether 4by4 plot
library(cowplot)
pgrid <- plot_grid(plotlist = p.l, align = "h", ncol = 3)
save_plot("nAnTi_sub_corr_3by2.pdf", pgrid,
          ncol =3, 
          nrow =2, 
          base_aspect_ratio = 1)

# plots, set up layout in R altogether horizontal plot
pgrid <- plot_grid(plotlist = p.l, align = "h", ncol = 5)
save_plot("nAnTi_sub_corr_horiz.pdf", pgrid,
          ncol = 5, 
          nrow = 1, 
          base_aspect_ratio = 1)

## Extract tag clusters
sampleLabels <- unname(sampleLabels(myCAGEset))
tc.l <- lapply(sampleLabels, function (x) tagClusters(myCAGEset, sample = x, 
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
saveRDS(object = tc.grl, file = "tagClusters_list.RDS")

## Interquantile widths
# selected samples for main paper
samples <- sampleLabels(myCAGEset)

tc_selected <- tc.grl[samples]

# extract interquantile widths
iq.l <- lapply(tc_selected, function(x) data.frame(iq_width = x$interquantile_width))

# add names to as a column to each sample
for (i in 1:length(samples)) {
  iq.l[[i]]$sample <- rep(samples[[i]], nrow(iq.l[[i]]))
}

# combine to 1 dataframe
iq.df <- do.call(rbind, iq.l)

# set levels for plotting - sample levels
iq.df$sample <- factor(iq.df$sample, levels = samples)

# plotting
library(ggplot2)
p <- ggplot(iq.df, aes(iq_width)) +
  geom_histogram(aes(x = iq_width, y = ..ncount..* 100), binwidth = 1, 
                 fill = "gray60", col = "black", size = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 16, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black", fill = "gray87")) +
  labs(x = "Interquantile width", y = "Percentage") +
  coord_equal(ratio = 1) +
  xlim(0, 150)

pdf(file = "IQ_width_sel_nAnTi_sub_sc.pdf", height = 6, width = 12)
p + facet_wrap(~ sample, ncol = 3)
dev.off()