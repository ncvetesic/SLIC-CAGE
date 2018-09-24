#/*==========================================================================#*/
#' ## Figure 4F, G - H3K4me3 coverage, metaplot and heatmap
#'
#/*==========================================================================#*/

library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
# import bam file as GRanges
H3K4me3_bam <- import("../data/ENCODE_H3K4me3/merged2rep.bam")

# convert to GRanges
H3K4me3_bam.gr <- GRanges(seqnames = seqnames(H3K4me3_bam),
                          ranges = IRanges(start = start(H3K4me3_bam),
                                           end = end(H3K4me3_bam)),
                          strand = strand(H3K4me3_bam))

seqlevelsStyle(H3K4me3_bam.gr) <- "UCSC"

# remove bam object  
rm(H3K4me3_bam)

# add sequence information
seqinfo(H3K4me3_bam.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(H3K4me3_bam.gr)]

# keep only standard chromosomes
H3K4me3_bam.gr <- keepStandardChromosomes(H3K4me3_bam.gr, pruning.mode = "coarse")

# calculate coverage separately for + and - strand
H3K4me3_bam_plus.cov <- coverage(H3K4me3_bam.gr[strand(H3K4me3_bam.gr) == "+", ])
H3K4me3_bam_minus.cov <- coverage(H3K4me3_bam.gr[strand(H3K4me3_bam.gr) == "-", ])

# calculated subtracted coverage (plus-minus strand)
h3k4me3_subtracted.cov = H3K4me3_bam_plus.cov - H3K4me3_bam_minus.cov

# create coverage heatmaps around domTSS
# load RDS object centered on dominant TSS (merged lanes and replicates)
domTSS_E14slic_mergLandR.grl <- readRDS(file = "../intermediate_data/domTSS_E14slic_mergLandR_grl.RDS")

# plot coverage heatmaps for each sample, sorted from sharp to broad and make clusters on transition of sharp to broad
library(heatmaps)

# create windows for plotting
up <- 100
down <- 300
range <- c(-up, down)
win <- up + down

domTSS_E14_win.grl <- sapply(domTSS_E14slic_mergLandR.grl, function(x) promoters(x, up = up, down = down))
lapply(domTSS_E14_win.grl, length)

# filter out values with dominantTSS values below 5 tpm
domTSS_E14_win.grl <- lapply(domTSS_E14_win.grl, function(x) x[x$tpm.dominant_ctss >= 5])

# remove out of bound ranges
domTSS_E14_win.grl <- sapply(domTSS_E14_win.grl, function(x) x[width(trim(x)) == win])
lapply(domTSS_E14_win.grl, length)

# make index of sharp and broad promoters
sIdx.l <- lapply(domTSS_E14_win.grl, function(x) x$interquantile_width <= 3)
bIdx.l <- lapply(domTSS_E14_win.grl, function(x) x$interquantile_width > 3)

# extract sequence 
domTSS_E14_seq.l <- sapply(domTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# plot metaplots 
# produce heatmaps prior to plotting (no smoothing)
hm_broad.l <- list()
hm_sharp.l <- list()

# plot coverage using heatmaps - separately for broad and for sharp
for (i in 1:length(domTSS_E14_seq.l)) {
  hm_broad.l[[i]] <- CoverageHeatmap(domTSS_E14_win.grl[[i]][bIdx.l[[i]]],
                                     h3k4me3_subtracted.cov,
                                     coords = range)
  hm_sharp.l[[i]] <- CoverageHeatmap(domTSS_E14_win.grl[[i]][sIdx.l[[i]]],
                                     h3k4me3_subtracted.cov,
                                     coords = range)
}

#domTSS_E14_win.grl <- lapply(domTSS_E14_win.grl, function(x) dropSeqlevels(x, value = "chrM", pruning.mode = "coarse"))

# convert heatmaps to metaplot data - returns df 
meta_broad.l <- lapply(hm_broad.l, function(x) dataMetaplot(x, binsize = 3))
meta_sharp.l <- lapply(hm_sharp.l, function(x) dataMetaplot(x, binsize = 3))

# create a list of dataframes (with sharp and broad column)
H3K4me3_df.all <- list()

samples <- c("SLIC 5 ng", "SLIC 10 ng", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi 5 ug")

for (i in 1:length(samples)) {
  H3K4me3_df.all[[i]] <- data.frame("x_coord" = meta_broad.l[[i]]$x_coord, 
                                    "sharp_H3K4me3" =  meta_sharp.l[[i]]$occurrence,
                                    "broad_H3K4me3" = meta_broad.l[[i]]$occurrence,
                                    "sharp_binsums" = meta_sharp.l[[i]]$bin_sums,
                                    "broad_binsums" = meta_broad.l[[i]]$bin_sums,
                                    "sample" = rep(samples[i], times = length(meta_sharp.l[[1]]$occurrence)))
}
names(H3K4me3_df.all) <- samples

# collapse the list into a dataframe
H3K4me3_df <- do.call(rbind, H3K4me3_df.all)
rownames(H3K4me3_df) <- 1:nrow(H3K4me3_df)

# omit relative signal - occurence
H3K4me3_df <- H3K4me3_df[, c("x_coord", "sharp_binsums", "broad_binsums", "sample")]

# collapse df to ggplot compatible df
library(tidyr)
H3K4me3_df.gg <- gather(H3K4me3_df, 
                        ... = colnames(H3K4me3_df[2:3]),
                        key = "type",
                        value = "H3K4me3")

# set levels for plotting
H3K4me3_df.gg$type <- factor(H3K4me3_df.gg$type, levels = c("sharp_binsums", "broad_binsums"))

# plot metaplots using ggplot2 - but using only absolute signal
library(ggplot2)
p <- ggplot(H3K4me3_df.gg) +
  geom_line(lty = 1, aes(x = x_coord, y = H3K4me3, colour = type), size = 0.6) +
  scale_color_manual(values = c("firebrick", "royalblue4"), name = NULL, labels = c( "sharp", 
                                                                                     "broad")) +
  theme_bw() +
  theme(text = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(reverse=T)) +
  geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
  geom_hline(yintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
  scale_y_continuous(name = "H3K4me3 subtracted coverage") +
  scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))

p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/H3K4me3_positioning_metaplot_broad_sharp.pdf", height = 12, width = 12)
print(p)
dev.off()

# print selected 10 ng and 5 ug samples
library(ggplot2)
p <- ggplot(H3K4me3_df.gg[H3K4me3_df.gg$sample %in% c("SLIC 10 ng", "nAnTi 5 ug"), ]) +
  geom_line(lty = 1, aes(x = x_coord, y = H3K4me3, colour = type), size = 0.6) +
  scale_color_manual(values = c("firebrick", "royalblue4"), name = NULL, labels = c( "sharp", 
                                                                                     "broad")) +
  theme_bw() +
  theme(text = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(reverse=T)) +
  geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
  geom_hline(yintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
  scale_y_continuous(name = "H3K4me3 subtracted coverage") +
  scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))

p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/H3K4me3_positioning_metaplot_broad_sharp_selected.pdf", height = 4, width = 11)
print(p)
dev.off()


# create 1 heatmap sorted per IQ width - range -500 to 500
# create windows for plotting
up <- 500
down <- 500
range <- c(-up, down)
win <- up + down

domTSS_E14_win.grl <- sapply(domTSS_E14slic_mergLandR.grl, function(x) promoters(x, up = up, down = down))
lapply(domTSS_E14_win.grl, length)

# remove out of bound ranges
domTSS_E14_win.grl <- sapply(domTSS_E14_win.grl, function(x) x[width(trim(x)) == win])
lapply(domTSS_E14_win.grl, length)

# extract sequence 
domTSS_E14_seq.l <- sapply(domTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))


sort_idx.l <- lapply(domTSS_E14_win.grl, function(x) order(x$interquantile_width, decreasing = FALSE))
hm_all.l <- list()

# add cluster information - sharp and broad clusters
cluster.l <- lapply(domTSS_E14_win.grl, function(x) as.numeric(x$interquantile_width <= 3))


for (i in 1:length(domTSS_E14_win.grl)) {
  hm_all.l[[i]] <- CoverageHeatmap(domTSS_E14_win.grl[[i]][sort_idx.l[[i]]],
                                   track =  h3k4me3_subtracted.cov,
                                   coords = range,
                                   label = paste("H3K4me3, ", names[[i]], sep = ""))
  scale(hm_all.l[[i]]) <- c(-30, 30)
  pdf(paste("../final_figures/H3K4me3_coverage_",names(domTSS_E14_seq.l)[i], ".pdf", sep = ""), height = 8, width = 6)
  plotHeatmapList(hm_all.l[[i]],
                  cex.label=1.5,
                  color = c("red", "white", "blue"),
                  partition = c(sum(cluster.l[[i]] == 1), sum(cluster.l[[i]] == 0)),
                  partition.lines = TRUE,
                  legend = TRUE,
                  legend.pos = "l",
                  legend.width = 0.15)
  dev.off()
}

#'------------- Malcolm's metaplot function from heatmaps package with plotting removed -------------#

#' Plot a Meta-region plot from heatmaps
#'
#' @param hm_list A list of heatmaps
#' @param binsize Integer, size of bins to use in plot
#' @param colors Color to use for each heatmap
#' @param addReferenceLine Logical, add reference line at zero or not
#'
#' This function creates a meta-region plot from 1 or more heatmaps with the same
#' coordinates. A meta-region plot graphs the sum of the signal at each position in
#' each heatmap rather than visualising the signal in two dimensions. Often binning
#' is required to smooth noisy signal.
#'
#' @return invisible(0)
#'
#' @export
#' @importFrom graphics plot mtext legend axis lines
#' @examples
#' data(HeatmapExamples)
#' plotHeatmapMeta(hm, color="steelblue")

dataMetaplot = function(hm_list, binsize=1, colors=gg_col(length(hm_list)), addReferenceLine=FALSE) {
  if (class(hm_list) == "Heatmap") hm_list = list(hm_list) # allow single heatmap argument
  
  if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
    stop("heatmaps must have the same coordinates")
  
  
  if (!length(unique(lapply(hm_list, function(x) xm(x)))) == 1)
    stop("heatmaps must have the same xm values")
  
  n_seq = unique(sapply(hm_list, function(x) x@nseq))
  
  if (!length(n_seq) == 1)
    stop("heatmaps must have the same number of sequences")
  
  coords = hm_list[[1]]@coords
  if (binsize != 1) {
    if (!all(xm(hm_list[[1]]) == 1:width(hm_list[[1]]))) {
      stop("cannot set binsize for heatmaps which are already binned/smoothed")
    }
    breaks = seq(0, width(hm_list[[1]]), by=binsize)
    bin_sums = lapply(hm_list, bin_heatmap, breaks=breaks)
    x_coord = breaks[1:(length(breaks)-1)] + binsize/2 + coords[1]
  } else {
    breaks = xm(hm_list[[1]])
    x_coord = breaks + binsize/2 + coords[1]
    bin_sums = lapply(hm_list, function(x) colSums(image(x)))
  }
  scale_factor = n_seq*(width(hm_list[[1]])/length(x_coord))
  occurrence = lapply(bin_sums, function(x) x/scale_factor)
  max_value = max(vapply(occurrence, max, numeric(1)))
  output = data.frame("x_coord" = x_coord, "occurrence" = occurrence[[1]], "bin_sums" = bin_sums[[1]])
  return(output)
}

#' @importFrom stats aggregate
bin_heatmap = function(hm, breaks) {
  partition = data.frame(pos=xm(hm), value=colSums(image(hm)), bin=cut(xm(hm), breaks))
  aggregate(partition$value, sum, by=list(bin=partition$bin))$x
}

#' @importFrom grDevices hcl
gg_col = function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
```
