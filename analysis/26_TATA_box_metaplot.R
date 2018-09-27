#/*==========================================================================#*/
#' ## Supplemental Figure S17 - TATA-box metaplot
#'
#/*==========================================================================#*/

## TATA-box metaplots - supplementary
library(BSgenome.Mmusculus.UCSC.mm10)
library(heatmaps)

# load RDS object centered on dominant TSS (merged lanes and replicates)
domTSS_E14slic_mergLandR.grl <- readRDS(file = "../intermediate_data/domTSS_E14slic_mergLandR_grl.RDS")

# checked distribution of IQ-widths - will use <= 3 as definition of sharp
lapply(tc_E14_mergLandR.grl, function(x) sum(x$interquantile_width <= 3)/length(x$interquantile_width))

# create windows for plotting - rather zoomed in for TBP positioning
up <- 100
down <- 300
range <- c(-up, down)
win <- up + down

domTSS_E14_win.grl <- sapply(domTSS_E14slic_mergLandR.grl, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
domTSS_E14_win.grl <- sapply(domTSS_E14_win.grl, function(x) x[width(trim(x)) == win])

# check the length of the each range in the list
sapply(domTSS_E14_win.grl, length)

# extract sequence 
domTSS_E14_seq.l <- sapply(domTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# plot TBPpwm positioning
library(seqPattern)
data(TBPpwm)

# using seqpattern to get 80 % match to TBPpwm...
tbp_scores.l <- list()

for(i in 1:length(domTSS_E14_seq.l)) {
  tbp_scores.l[[i]] <- motifScanHits(domTSS_E14_seq.l[[i]], motifPWM = TBPpwm, minScore = "80%")
}
# add names to TBP scores
names(tbp_scores.l) <- names(domTSS_E14_seq.l)

# create matrices for each of the samples 
matrix.l <- list()
for (i in 1:length(domTSS_E14_seq.l)) {
  matrix.l[[i]] <- matrix(data = 0, nrow = length(domTSS_E14_seq.l[[i]]), ncol = unique(width(domTSS_E14_seq.l[[i]])))
}

names(matrix.l) <- names(domTSS_E14_seq.l)


# add values to empty matrices (1 where matching is >80%)
for (i in 1:length(matrix.l)) {
  matrix.l[[i]] <- replace(matrix.l[[i]], as.matrix(tbp_scores.l[[i]][, 1:2]), 1)
}
names(matrix.l) <- names(tbp_scores.l)

# create index for TBP matrix separation into broad and sharp
sIdx.l <- lapply(domTSS_E14_win.grl, function(x) x$interquantile_width <= 3)
bIdx.l <- lapply(domTSS_E14_win.grl, function(x) x$interquantile_width > 3)


# use heatmaps - need heatmaps object prior to making metaplots
library(heatmaps)
hm_broad.l <- list()
hm_sharp.l <- list()

# convert TBP scores to heatmaps object
names <- c("SLIC 5 ng", "SLIC 10 ng", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "SLIC 5 ug")    

for (i in 1:length(tbp_scores.l)) {
  hm_sharp.l[[i]] <- new("Heatmap",
                         image = matrix.l[[i]][sIdx.l[[i]], ],
                         scale = c(0, 1),
                         coords = as.integer(range),
                         nseq = nrow(matrix.l[[i]][sIdx.l[[i]], ]),
                         metadata = list())
  
  hm_broad.l[[i]] <- new("Heatmap",
                         image = matrix.l[[i]][bIdx.l[[i]], ],
                         scale = c(0, 1),
                         coords = as.integer(range),
                         nseq = nrow(matrix.l[[i]][bIdx.l[[i]], ]),
                         metadata = list())
}

#add names to heatmaps
names(hm_sharp.l) <- names(tbp_scores.l)
names(hm_broad.l) <- names(tbp_scores.l)

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
  scale_y_continuous(name = "TATA-box relative frequency", limits = c(0, 0.08)) +
  scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))

p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/TATA_metaplot_sharp_broad.pdf", height = 12, width = 11)
print(p)
dev.off()

# plot selected samples (if I want it in the main paper)
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
  scale_y_continuous(name = "Relative frequency", limits = c(0, 0.08)) +
  scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))

p <- p + facet_wrap(~ sample, ncol = 2)

pdf("../final_figures/TATA_metaplot_sharp_broad_selected.pdf", height = 4, width = 11)
print(p)
dev.off()

