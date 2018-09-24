#/*==========================================================================#*/
#' ## Figure 4D - nucleosome positioning signal / WW periodicity metaplot
#'
#/*==========================================================================#*/

library(BSgenome.Mmusculus.UCSC.mm10)
library(seqPattern)

# load RDS object centered on dominant TSS (merged lanes and replicates)
domTSS_E14slic_mergLandR.grl <- readRDS(file = "../intermediate_data/domTSS_E14slic_mergLandR_grl.RDS")

# create windows centered on domTSS -100, +300
up <- 100
down <- 300
range <- c(-up, down)
win <- up + down

domTSS_E14_win.grl <- sapply(domTSS_E14slic_mergLandR.grl, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
domTSS_E14_win.grl <- sapply(domTSS_E14_win.grl, function(x) x[width(trim(x)) == win])

# make index of sharp and broad promoters
sIdx.l <- lapply(domTSS_E14_win.grl, function(x) x$interquantile_width <= 3)
bIdx.l <- lapply(domTSS_E14_win.grl, function(x) x$interquantile_width > 3)

# check the length of the each range in the list
sapply(domTSS_E14_win.grl, length)

# extract sequence 
domTSS_E14_seq.l <- sapply(domTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# use malcolms heatmaps to create smoothed pattern heatmaps and sum it into a metaplot..- plot using ggplot2
library(heatmaps)

# plot dinucleotide plots centered on dominant TSS
pattern <- "WW"
hm.l <- list()
hm_smoothed.l <- list()

# create sorting index for heatmaps
sort.sharp.l <- list()
sort.broad.l <- list()

for(i in 1:length(domTSS_E14_win.grl)) {
  sort.sharp.l[[i]] <- order(domTSS_E14_win.grl[[i]][sIdx.l[[i]]]$interquantile_width, decreasing = FALSE)
  sort.broad.l[[i]] <- order(domTSS_E14_win.grl[[i]][bIdx.l[[i]]]$interquantile_width, decreasing = FALSE)
}
names(sort.sharp.l) <- names(domTSS_E14_win.grl)
names(sort.broad.l) <- names(domTSS_E14_win.grl)

# plot metaplots 
# produce heatmaps prior to plotting (no smoothing)
hm_broad.l <- list()
hm_sharp.l <- list()

for(i in 1:length(domTSS_E14_seq.l)) {
  hm_broad.l[[i]] <- PatternHeatmap(domTSS_E14_seq.l[[i]][bIdx.l[[i]]][sort.broad.l[[i]]], pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
  hm_sharp.l[[i]] <- PatternHeatmap(domTSS_E14_seq.l[[i]][sIdx.l[[i]]][sort.sharp.l[[i]]], pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
}
names(hm_broad.l) <- names(domTSS_E14_seq.l)
names(hm_sharp.l) <- names(domTSS_E14_seq.l)

# convert heatmaps to metaplot data - returns df 
meta_broad.l <- lapply(hm_broad.l, function(x) dataMetaplot(x, binsize = 2))
meta_sharp.l <- lapply(hm_sharp.l, function(x) dataMetaplot(x, binsize = 2))

# create a list of dataframes (with sharp and broad column)
all_df.l <- list()

samples <- c("SLIC 5 ng", "SLIC 10 ng", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi 5 ug")

for (i in 1:length(samples)) {
  all_df.l[[i]] <- data.frame("x_coord" = meta_broad.l[[i]]$x_coord, 
                              "sharp_occurrence" =  meta_sharp.l[[i]]$occurrence,
                              "broad_occurrence" = meta_broad.l[[i]]$occurrence,
                              "sample" = rep(samples[i], times = length(meta_sharp.l[[1]]$occurrence)))
}
names(all_df.l) <- samples

# collapse the list into a dataframe
all_df <- do.call(rbind, all_df.l)
rownames(all_df) <- 1:nrow(all_df)

# collapse df to ggplot compatible df
library(tidyr)
all_df.gg <- gather(all_df, 
                    ... = colnames(all_df[2:3]),
                    key = "type",
                    value = "occurence")


# plot metaplots using ggplot2
library(ggplot2)
p <- ggplot(all_df.gg) +
  geom_line(lty = 1, aes(x = x_coord, y = occurence, colour = type), size = 0.6) +
  scale_color_manual(values = c("royalblue4", "firebrick"), name = NULL, labels = c( "broad", "sharp")) +
  theme_bw() +
  theme(text = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(reverse=T)) +
  geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
  scale_y_continuous(name = "Relative frequency", limits = c(0, 0.4), breaks = c(0, 0.2, 0.4)) +
  scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))

p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/WW_sharp_broad.pdf", height = 12, width = 11)
print(p)
dev.off()

# plot selected samples for Figure3
all_df_sub.gg <- all_df.gg[all_df.gg$sample %in% c("SLIC 10 ng", "nAnTi 5 ug"), ]
all_df_sub.gg$sample <- droplevels(all_df_sub.gg$sample)

library(ggplot2)
p <- ggplot(all_df_sub.gg) +
  geom_line(lty = 1, aes(x = x_coord, y = occurence, colour = type), size = 0.6) +
  scale_color_manual(values = c("royalblue4", "firebrick"), name = NULL, labels = c( "broad", "sharp")) +
  theme_bw() +
  theme(text = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(reverse=T)) +
  geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
  scale_y_continuous(name = "Relative frequency", limits = c(0, 0.4), breaks = c(0, 0.2, 0.4)) +
  scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))

p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/WW_sharp_broad_selected.pdf", height = 4, width = 11)
print(p)
dev.off()


# zoomed metaplots for the inset = only broad pattern occurence
library(ggplot2)
p <- ggplot(all_df.gg[all_df.gg$type %in% "broad_occurrence", ]) +
  geom_line(lty = 1, aes(x = x_coord, y = occurence, colour = type), size = 0.6, legend = FALSE) +
  scale_color_manual(values = "royalblue4", name = NULL) +
  theme_bw() +
  theme(text = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(55, 200), ylim = c(0.11, 0.14))
p <- p + facet_wrap(~ sample, ncol = 2)

# choose random TSS in broad promoters (not dominant) equal length as dominant and plot WW
# broad as defined as > 3, 
num_broad.l <- lapply(bIdx.l, sum)

lapply(domTSS_E14_win.grl, function(x) x$interquantile_width > 3)


p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/WW_sharp_broad_zoom.pdf", height = 12, width = 10)
print(p)
dev.off()

# check nucleosome positioning signal in windows centered on randomly selected true CTSS
#extract CTSS dataframe with normalized TPM per CTSS per sample
CTSSnorm <- CTSSnormalizedTpm(E14replicates_merge)

# filter to have at least 1TPM 
filtr_idx.list <- list()
samples_all <- colnames(CTSSnorm)[-c(1:3)]

for (i in 1:length(samples_all)) {
  filtr_idx.list[[i]] <- c(CTSSnorm[, samples_all[i]] >= 1)
}
names(filtr_idx.list) <- samples_all

# create GRanges object CTSSs (filter per sample to have >= 1 TPM)
filtr.ctss.grl <- lapply(filtr_idx.list, function(x) GRanges(seqnames = CTSSnorm[x, ]$chr, 
                                                             ranges = IRanges(start = CTSSnorm[x, ]$pos, width = 1),
                                                             strand = CTSSnorm[x, ]$strand))
# import tag cluster object
tc_E14_mergLandR.grl <- readRDS("../intermediate_data/tc_E14SLICmergedRep_grl.RDS")


# overlap tag clusters and CTSS - to connect which CTSSs are in which tagClusters, but only CTSSs with >= 1 TPM (not all CTSS)
tc_CTSS_idx.l <- list()

for (i in 1:length(tc_E14_mergLandR.grl)) {
  tc_CTSS_idx.l[[i]] <- findOverlaps(tc_E14_mergLandR.grl[[i]], filtr.ctss.grl[[i]], select = "all")
  # convert to a dataframe
  tc_CTSS_idx.l[[i]] <- as.data.frame(tc_CTSS_idx.l[[i]])
}
names(tc_CTSS_idx.l) <- names(tc_E14_mergLandR.grl)

# for each tag cluster select a random real CTSS to use as dominant - end object should be of equal length as each tag cluster
# collapse by tag cluster number
# create an index of tag clusters that have overlapping CTSSs with tpm >=1 
tc.idx.l <- lapply(tc_CTSS_idx.l, function(x) unique(x$queryHits))

# sample random CTSS per overlapping tag cluster 
# select random CTSS per tag cluster - take indexes of CTSSs that overlap tag clusters, and sample one of them
random.ctss.l <- list()
for (i in 1:length(tc.idx.l)) {
  random.ctss.l[[i]] <- sapply(tc.idx.l[[i]], function(x) sample(tc_CTSS_idx.l[[i]][tc_CTSS_idx.l[[i]]$queryHits == x, "subjectHits"], 1))
}
names(random.ctss.l) <- tc_E14_mergLandR.grl

# subset tag-clusters in each sample based on the ones that where overlapped with ctss's
tc_ctss_sub.l <- list()

for (i in 1:length(tc.idx.l)) {
  tc_ctss_sub.l[[i]] <- tc_E14_mergLandR.grl[[i]][tc.idx.l[[i]]]
}
names(tc_ctss_sub.l) <- names(tc_E14_mergLandR.grl)

# get random CTSS coordinates (I have indices)
random.ctss.coord.l <- list()

for(i in 1:length(random.ctss.l)) {
  random.ctss.coord.l[[i]] <- filtr.ctss.grl[[i]][random.ctss.l[[i]], ]
}

# attach CTSS information to corresponding tag cluster
for (i in 1:length(tc_ctss_sub.l)) {
  tc_ctss_sub.l[[i]]$random_ctss_coord <- start(random.ctss.coord.l[[i]])
  tc_ctss_sub.l[[i]]$random_ctss_chr <- chrom(random.ctss.coord.l[[i]])
  seqinfo(random.ctss.coord.l[[i]]) <- seqinfo(tc_ctss_sub.l[[i]])[seqlevels(random.ctss.coord.l[[i]])]
  random.ctss.coord.l[[i]]$interquantile_width <- tc_ctss_sub.l[[i]]$interquantile_width
}
names(random.ctss.coord.l) <- names(tc_ctss_sub.l)

# # create windows centered on random TSS from that tag cluster -100, +300
up <- 100
down <- 300
range <- c(-up, down)
win <- up + down

randomTSS_E14_win.grl <- sapply(random.ctss.coord.l, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
randomTSS_E14_win.grl <- sapply(randomTSS_E14_win.grl, function(x) x[width(trim(x)) == win])

# check the length of the each range in the list
sapply(randomTSS_E14_win.grl, length)

# extract sequence 
randomTSS_E14_seq.l <- sapply(randomTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# use malcolms heatmaps to create smoothed pattern heatmaps and sum it into a metaplot..- plot using ggplot2
library(heatmaps)

# plot dinucleotide plots centered on dominant TSS
pattern <- "WW"
hm.l <- list()
hm_smoothed.l <- list()

# make index of sharp and broad promoters
sIdx.l <- lapply(randomTSS_E14_win.grl, function(x) x$interquantile_width <= 3)
bIdx.l <- lapply(randomTSS_E14_win.grl, function(x) x$interquantile_width > 3)

# plot metaplots 
# produce heatmaps prior to plotting (no smoothing)
hm_broad.l <- list()
hm_sharp.l <- list()

for(i in 1:length(domTSS_E14_seq.l)) {
  hm_broad.l[[i]] <- PatternHeatmap(randomTSS_E14_seq.l[[i]][bIdx.l[[i]]], pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
  hm_sharp.l[[i]] <- PatternHeatmap(randomTSS_E14_seq.l[[i]][sIdx.l[[i]]], pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
}
names(hm_broad.l) <- names(randomTSS_E14_seq.l)
names(hm_sharp.l) <- names(randomTSS_E14_seq.l)

# convert heatmaps to metaplot data - returns df 
meta_broad.l <- lapply(hm_broad.l, function(x) dataMetaplot(x, binsize = 2))
meta_sharp.l <- lapply(hm_sharp.l, function(x) dataMetaplot(x, binsize = 2))

# create a list of dataframes (with sharp and broad column)
random_df.all <- list()

samples <- c("SLIC 5 ng", "SLIC 10 ng", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi 5 ug")

for (i in 1:length(samples)) {
  random_df.all[[i]] <- data.frame("x_coord" = meta_broad.l[[i]]$x_coord, 
                                   "random_sharp_occurrence" =  meta_sharp.l[[i]]$occurrence,
                                   "random_broad_occurrence" = meta_broad.l[[i]]$occurrence,
                                   "sample" = rep(samples[i], times = length(meta_sharp.l[[1]]$occurrence)))
}
names(random_df.all) <- samples

# collapse the list into a dataframe
random_df <- do.call(rbind, random_df.all)
rownames(random_df) <- 1:nrow(random_df)

# collapse df to ggplot compatible df
library(tidyr)
random_df.gg <- gather(random_df, 
                       ... = colnames(random_df[2:3]),
                       key = "type",
                       value = "occurence")

# combine random and normal TSS data
all_combined.gg <- rbind(all_df.gg, random_df.gg)

#set levels for plotting
all_combined.gg$type <- factor(all_combined.gg$type, levels = c("sharp_occurrence", "random_sharp_occurrence",
                                                                "broad_occurrence", "random_broad_occurrence"))

# plot metaplots using ggplot2
library(ggplot2)
p <- ggplot(all_combined.gg) +
  geom_line(lty = 1, aes(x = x_coord, y = occurence, colour = type), size = 0.6) +
  scale_color_manual(values = c("firebrick", "darksalmon", "royalblue4", "cornflowerblue"), name = NULL, labels = c( "sharp", "sharp random", "broad", "broad random")) +
  theme_bw() +
  theme(text = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(reverse=T)) +
  geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
  scale_y_continuous(name = "Relative frequency", limits = c(0, 0.4), breaks = c(0, 0.2, 0.4)) +
  scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))

p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/WW_sharp_broad_all_domTSS_random.pdf", height = 12, width = 12)
print(p)
dev.off()

# zoomed metaplots for the inset = only broad pattern occurence
library(ggplot2)
p <- ggplot(all_combined.gg[random_df.gg$type %in% c("random_broad_occurrence", "broad_occurrence"),  ]) +
  geom_line(lty = 1, aes(x = x_coord, y = occurence, colour = type), size = 0.6, legend = FALSE) +
  scale_color_manual(values = c("royalblue4", "goldenrod2"), name = NULL, labels = c( "broad", "broad random")) +
  theme_bw() +
  theme(text = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(reverse=T)) +
  scale_y_continuous(name = "Relative frequency") +
  scale_x_continuous(name = "Distance to dominant TSS (bp)") +
  coord_cartesian(xlim = c(55, 200), ylim = c(0.11, 0.14))
p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/WW_sharp_broad_random_zoom.pdf", height = 12, width = 12)
print(p)
dev.off()

# print only selected samples - 5 ng and nanti
all_combined_sub.gg <- all_combined.gg[all_combined.gg$sample %in% c("SLIC 10 ng", "nAnTi 5 ug"), ]
all_combined_sub.gg$sample <-droplevels(all_combined_sub.gg$sample)

# zoomed metaplots for the inset = only broad pattern occurence for selected samples
library(ggplot2)
p <- ggplot(all_combined_sub.gg[all_combined_sub.gg$type %in% c("random_broad_occurrence", "broad_occurrence"),  ]) +
  geom_line(lty = 1, aes(x = x_coord, y = occurence, colour = type), size = 0.6, legend = FALSE) +
  scale_color_manual(values = c("royalblue4", "goldenrod2"), name = NULL, labels = c( "broad", "broad random")) +
  theme_bw() +
  theme(text = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(reverse=T)) +
  scale_y_continuous(name = "Relative frequency") +
  scale_x_continuous(name = "Distance to dominant TSS (bp)") +
  coord_cartesian(xlim = c(55, 200), ylim = c(0.11, 0.14))
p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/WW_sharp_broad_random_selected_zoom.pdf", height = 4, width = 12)
print(p)
dev.off()