#/*==========================================================================#*/
#' ## Figure 4E - tag cluster coverage heatmap
#'
#/*==========================================================================#*/

# create coverage heatmaps around domTSS
# load RDS object centered on dominant TSS (merged lanes and replicates)
domTSS_E14slic_mergLandR.grl <- readRDS(file = "../intermediate_data/domTSS_E14slic_mergLandR_grl.RDS")

# plot coverage heatmaps for each sample, sorted from sharp to broad and make clusters on transition of sharp to broad
library(heatmaps)

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

# add sample names
samples <- c("SLIC 5 ng", "SLIC 10 ng", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi 5 ug")

library(RColorBrewer)

for (i in 1:length(domTSS_E14_win.grl)) {
  hm_all.l[[i]] <- CoverageHeatmap(domTSS_E14_win.grl[[i]][sort_idx.l[[i]]],
                                   track =  tc_E14_mergLandR.grl[[i]],
                                   coords = range,
                                   weight = log10(tc_E14_mergLandR.grl[[i]]$tpm + 1),
                                   label = paste("CAGE coverage, ", samples[[i]], sep = ""))
  scale(hm_all.l[[i]]) <- c(0,3)
  
  pdf(paste("../final_figures/tc_coverage_", names(domTSS_E14_seq.l)[i], ".pdf", sep = ""), height = 8, width = 6)
  plotHeatmapList(hm_all.l[[i]],
                  cex.label=1.5,
                  color = "Greys",
                  partition = c(sum(cluster.l[[i]] == 1), sum(cluster.l[[i]] == 0)),
                  partition.lines = TRUE,
                  legend = TRUE,
                  legend.pos = "l",
                  legend.width = 0.15)
  dev.off()
}

