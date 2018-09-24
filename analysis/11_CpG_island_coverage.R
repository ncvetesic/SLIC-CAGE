#/*==========================================================================#*/
#' ## Figure 4H, I - CpG island coverage - heatmap and metaplot
#'
#/*==========================================================================#*/

library(BSgenome.Mmusculus.UCSC.mm10)
# import UCSC CpG islands track
CpG_UCSC <- read.table("../data/UCSC_CpG_islands_mm10.txt", header = FALSE, sep = "\t")
colnames(CpG_UCSC) <- c("bin", "chrom", "chromStart", "chromEnd", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")

# convert CpG_USCS to GRanges object
CpG_UCSC.gr <- GRanges(seqnames = CpG_UCSC$chrom,
                       ranges = IRanges(start = CpG_UCSC$chromStart, end = CpG_UCSC$chromEnd),
                       strand = "*",
                       length = CpG_UCSC$length,
                       cpgnum = CpG_UCSC$cpgNum,
                       gcNum = CpG_UCSC$gcNum,
                       perCpG = CpG_UCSC$perCpg,
                       perGc = CpG_UCSC$perGc,
                       obsExp = CpG_UCSC$obsExp)
# add sequence information
seqinfo(CpG_UCSC.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqnames(seqinfo(CpG_UCSC.gr))]

# keep only standard chromosome
CpG_UCSC.gr <- keepStandardChromosomes(CpG_UCSC.gr, pruning.mode = "coarse")

# create a coverage heatmap sorted by iqWidth of tag clusters
library(heatmaps)

# overlap CpG with TC 
CpG.df <- data.frame("nr_CpG" = sapply(tc_E14.grl, function(x) sum(countOverlaps(x, CpG_UCSC.gr))),
                     "nr of TAG clusters" = sapply(tc_E14.grl, length),
                     "CpG perc" = round(sapply(tc_E14.grl, function(x) sum(countOverlaps(x, CpG_UCSC.gr)))/sapply(tc_E14.grl, length)*100, 2))

# read dominantTSS data
domTSS_E14slic_mergLandR.grl <- readRDS("../intermediate_data/domTSS_E14slic_mergLandR_grl.RDS")

# create windows for plotting
up <- 500
down <- 500
range <- c(-up, down)
win <- up + down

domTSS_E14_win.grl <- lapply(domTSS_E14slic_mergLandR.grl, function(x) promoters(x, up = up, down = down))
lapply(domTSS_E14_win.grl, length)

# remove out of bound ranges
domTSS_E14_win.grl <- lapply(domTSS_E14_win.grl, function(x) x[width(trim(x)) == win])
lapply(domTSS_E14_win.grl, length)

# remove chrM (as its not in CpGs)
domTSS_E14_win.grl <- lapply(domTSS_E14_win.grl, function(x) dropSeqlevels(x, value = "chrM", pruning.mode = "coarse"))

# extract sequence 
domTSS_E14_seq.l <- lapply(domTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))


sort_idx.l <- lapply(domTSS_E14_win.grl, function(x) order(x$interquantile_width, decreasing = FALSE))
hm_all.l <- list()

# add cluster information - sharp and broad clusters
cluster.l <- lapply(domTSS_E14_win.grl, function(x) as.numeric(x$interquantile_width <= 3))

# add sample names
samples <- c("SLIC 5 ng", "SLIC 10 ng", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi 5 ug")


library(RColorBrewer)

for (i in 1:length(domTSS_E14_win.grl)) {
  hm_all.l[[i]] <- CoverageHeatmap(domTSS_E14_win.grl[[i]][sort_idx.l[[i]]],
                                   track =  CpG_UCSC.gr,
                                   coords = range,
                                   label = paste("CpG islands, ", samples[[i]], sep = ""))
  
  
  
  pdf(paste("../final_figures/CpGislands_", names(domTSS_E14_seq.l)[i], ".pdf", sep = ""), height = 8, width = 6)
  plotHeatmapList(hm_all.l[[i]],
                  cex.label=1.5,
                  color = "Purples",
                  partition = c(sum(cluster.l[[i]] == 1), sum(cluster.l[[i]] == 0)),
                  partition.lines = TRUE,
                  legend = TRUE,
                  legend.pos = "l",
                  legend.width = 0.15)
  dev.off()
}

# create metaplots separately for broad and sharp
# make index of sharp and broad promoters
sIdx.l <- lapply(domTSS_E14_win.grl, function(x) x$interquantile_width <= 3)
bIdx.l <- lapply(domTSS_E14_win.grl, function(x) x$interquantile_width > 3)

# plot metaplots 
# produce heatmaps prior to plotting (no smoothing)
hm_broad.l <- list()
hm_sharp.l <- list()

# plot coverage using heatmaps - separately for broad and for sharp
for (i in 1:length(domTSS_E14_seq.l)) {
  hm_broad.l[[i]] <- CoverageHeatmap(domTSS_E14_win.grl[[i]][bIdx.l[[i]]],
                                     CpG_UCSC.gr,
                                     coords = range)
  hm_sharp.l[[i]] <- CoverageHeatmap(domTSS_E14_win.grl[[i]][sIdx.l[[i]]],
                                     CpG_UCSC.gr,
                                     coords = range)
}

# convert heatmaps to metaplot data - returns df 
meta_broad.l <- lapply(hm_broad.l, function(x) dataMetaplot(x, binsize = 3))
meta_sharp.l <- lapply(hm_sharp.l, function(x) dataMetaplot(x, binsize = 3))

# create a list of dataframes (with sharp and broad column)
CpG_df.all <- list()

samples <- c("SLIC 5 ng", "SLIC 10 ng", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi 5 ug")

for (i in 1:length(samples)) {
  CpG_df.all[[i]] <- data.frame("x_coord" = meta_broad.l[[i]]$x_coord, 
                                "sharp_CpG" =  meta_sharp.l[[i]]$occurrence,
                                "broad_CpG" = meta_broad.l[[i]]$occurrence,
                                "sharp_binsums" = meta_sharp.l[[i]]$bin_sums,
                                "broad_binsums" = meta_broad.l[[i]]$bin_sums,
                                "sample" = rep(samples[i], times = length(meta_sharp.l[[1]]$occurrence)))
}
names(CpG_df.all) <- samples

# collapse the list into a dataframe
CpG_df <- do.call(rbind, CpG_df.all)
rownames(CpG_df) <- 1:nrow(CpG_df)

# omit relative plots - occurrence
CpG_df <- CpG_df[, c("x_coord", "sharp_binsums", "broad_binsums", "sample")]

# collapse df to ggplot compatible df
library(tidyr)
CpG_df.gg <- gather(CpG_df, 
                    ... = colnames(CpG_df[2:3]),
                    key = "type",
                    value = "CpG")

# set levels for plotting
CpG_df.gg$type <- factor(CpG_df.gg$type, levels = c("sharp_binsums", "broad_binsums"))

# plot metaplots using ggplot2
library(ggplot2)
p <- ggplot(CpG_df.gg) +
  geom_line(lty = 1, aes(x = x_coord, y = CpG, colour = type), size = 0.6) +
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
  scale_y_continuous(name = "CpG islands coverage") +
  scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-500, 500), seq(-500, 500, by = 200))

p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/cpg_positioning_metaplot_broad_sharp.pdf", height = 12, width = 12)
print(p)
dev.off()

# print selected 10 ng and 5 ug samples
library(ggplot2)
p <- ggplot(CpG_df.gg[CpG_df.gg$sample %in% c("SLIC 10 ng", "nAnTi 5 ug"), ]) +
  geom_line(lty = 1, aes(x = x_coord, y = CpG, colour = type), size = 0.6) +
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
  scale_y_continuous(name = "CpG islands coverage") +
  scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-500, 500), seq(-500, 500, by = 200))

p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/cpg_positioning_metaplot_broad_sharp_selected.pdf", height = 4, width = 11)
print(p)
dev.off()
