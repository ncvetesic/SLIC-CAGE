#/*==========================================================================#*/
#' ## Figures 4A-C - heatmaps TA, TATA-box and GC density
#' 
#/*==========================================================================#*/

library(BSgenome.Mmusculus.UCSC.mm10)

# load RDS object centered on dominant TSS (merged lanes and replicates)
domTSS_E14slic_mergLandR.grl <- readRDS(file = "../intermediate_data/domTSS_E14slic_mergLandR_grl.RDS")

# checked distribution of IQ-widths - will use <= 3 as definition of sharp
lapply(tc_E14_mergLandR.grl, function(x) sum(x$interquantile_width <= 3)/length(x$interquantile_width))



# plot heatmaps for each sample, sorted from sharp to broad and make clusters on transition of sharp to broad
library(heatmaps)

# create windows for plotting
up <- 500
down <- 500
range <- c(-up, down)
win <- up + down

domTSS_E14_win.grl <- sapply(domTSS_E14slic_mergLandR.grl, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
domTSS_E14_win.grl <- sapply(domTSS_E14_win.grl, function(x) x[width(trim(x)) == win])

# check the length of the each range in the list
sapply(domTSS_E14_win.grl, length)

# extract sequence 
domTSS_E14_seq.l <- sapply(domTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# assign cluster numbers - 1 for iq.width <= 3, 0 for iq.width > 3
cluster.l <- lapply(domTSS_E14_win.grl, function(x) as.numeric(x$interquantile_width <= 3))

# sort index - no sequences are trimmed out so I can use the .gr object for sorting
sort_idx.l <- sapply(domTSS_E14_win.grl, function(x) order(x$interquantile_width, decreasing = FALSE))

# plot dinucleotide plots centered on dominant TSS
pattern <- c("TA", "GC")
hm.l <- list()
hm_smoothed.l <- list()
color.l <- c("Blues", "Reds")

# add scaling so its equal for each pattern
scale.l <- list(c(0, 0.12), c(0, 0.18))
names(scale.l) <- pattern
TA_images.l <- list()
GC_images.l <- list()

for (i in 1:length(pattern)) {
  for(j in 1:length(domTSS_E14_seq.l)) {
    hm.l[j] = PatternHeatmap(domTSS_E14_seq.l[[j]][sort_idx.l[[j]]], pattern = pattern[i], coords = range, label = paste(pattern[i], ", ", samples[j], sep = ""))
    
    hm_smoothed.l[j] = smoothHeatmap(hm.l[[j]], sigma = c(30, 3), output.size=c(10000, 1000))
    
    # change scale so its unique/equal for each pattern
    scale(hm_smoothed.l[[j]]) <- scale.l[[i]]
    
    # plot heatmaps
    pdf(paste("../final_figures/Figure3_analysis/domTSS_all_sharp_broad_",names(domTSS_E14_seq.l)[j],"_", pattern[[i]], ".pdf", sep = ""), height = 8, width = 6)
    plotHeatmapList(hm_smoothed.l[[j]],
                    cex.label = 1.5,
                    partition = c(sum(cluster.l[[j]] == 1), sum(cluster.l[[j]] == 0)),
                    partition.legend = FALSE,
                    partition.lines = TRUE, 
                    legend = TRUE,
                    legend.pos = "l",
                    legend.width = 0.15,
                    color = color.l[[i]])
    dev.off()
  }
  if(pattern[i] == "TA") {
    TA_images.l <- lapply(hm_smoothed.l, function(x) return(image(x)))
    names(TA_images.l) <- names(domTSS_E14_win.grl)
  } else {
    GC_images.l <- lapply(hm_smoothed.l, function(x) return(image(x)))
    names(GC_images.l) <- names(domTSS_E14_win.grl)
  }
}

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

# use heatmaps for plotting
library(heatmaps)
# convert TBP scores to heatmaps object
hm.l <- list()
hm_smoothed.l <- list()
TBP_images.l <- list()
names <- c("SLIC 5 ng", "SLIC 10 ng", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "SLIC 5 ug")    

for (i in 1:length(tbp_scores.l)) {
  hm.l[[i]] <- new("Heatmap",
                   image = matrix.l[[i]][sort_idx.l[[i]], ],
                   scale = c(0, 1),
                   coords = as.integer(range),
                   nseq = nrow(matrix.l[[i]]),
                   label = paste("TATA-box pwm \n", names[[i]], sep = ""),
                   metadata = list())
}
#add names to heatmaps
names(hm.l) <- names(tbp_scores.l)

# set smoothing and plotting
hm_smoothed.l <- list()

for (i in 1:length(hm.l)) {
  hm_smoothed.l[[i]] <- smoothHeatmap(hm.l[[i]],  algorithm = "kernel", sigma = c(30, 3), output.size=c(10000, 1000))
  scale(hm_smoothed.l[[i]]) <- c(0, 0.1)
  TBP_images.l[[i]] <- image(hm_smoothed.l[[i]])
  
  # plot heatmaps
  pdf(paste("../final_figures/Figure3_analysis/domTSS_all_sharp_broad_TBP_",names(domTSS_E14_seq.l)[i], ".pdf", sep = ""), height = 8, width = 6)
  plotHeatmapList(hm_smoothed.l[[i]],
                  cex.label = 1.5,
                  partition = c(sum(cluster.l[[j]] == 1), sum(cluster.l[[j]] == 0)),
                  partition.legend = FALSE,
                  partition.lines = TRUE, 
                  legend = TRUE,
                  legend.pos = "l",
                  legend.width = 0.15,
                  color = "Greens")
  dev.off()
}
