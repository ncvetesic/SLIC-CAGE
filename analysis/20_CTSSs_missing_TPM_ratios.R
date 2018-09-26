#/*==========================================================================#*/
#' ## Supplemental Figure S9, S10 and S13 - missing CTSSs and TPM ratios
#'
#/*==========================================================================#*/

#CTSS ratio analysis (each sample vs control nanti)
# extract CTSS dataframe with normalized TPM per CTSS per sample
LICAGE_ctss <- CTSSnormalizedTpm(LICAGE)

# calculate ratios for each sample vs BY4741 in a list
all_samples <- colnames(LICAGE_ctss[, 4:12])
control <- "BY4741"

LICAGE_ctss_sub.l <- list()
for (i in 1:length(all_samples)) {
  LICAGE_ctss_sub.l[[i]] <- LICAGE_ctss[, c("chr", "pos", "strand", all_samples[[i]], "BY4741")]   
}

# extract control BY4741 CTSSs with TPM value for a control coverage plot
control_CTSStpm <- LICAGE_ctss[, c("chr", "pos", "strand", "BY4741")]
# filter out zeros
control_CTSStpm <- control_CTSStpm[control_CTSStpm$BY4741 != 0, ]
# convert to GRanges
control_CTSStpm.gr <- GRanges(seqnames = control_CTSStpm$chr,
                              strand = control_CTSStpm$strand,
                              ranges = IRanges(start = control_CTSStpm$pos, width = 1),
                              tpm = control_CTSStpm$BY4741,
                              log10tpm = log10(control_CTSStpm$BY4741 + 1))
# add seqinfo
seqinfo(control_CTSStpm.gr) <- seqinfo(BSgenome.Scerevisiae.UCSC.sacCer3)

### --- function to calculate TPM ratios ---###  
CTSStpmRatio <- function(LICAGE_ctss_sub) {
  # calculate tpm ratios per CTSS (nantiCAGE)/(slic-CAGE))
  LICAGE_ctss_sub$ratio <- LICAGE_ctss_sub[, 5] / LICAGE_ctss_sub[, 4]
  
  # filter NaN values - they are result of 0 vs 0 division
  LICAGE_ctss_sub <- LICAGE_ctss_sub[!is.na(LICAGE_ctss_sub$ratio), ]
  
  # filter Inf values
  LICAGE_ctss_sub <- LICAGE_ctss_sub[!is.infinite(LICAGE_ctss_sub$ratio), ]
  
  # filter 0 values
  LICAGE_ctss_sub <- LICAGE_ctss_sub[!(LICAGE_ctss_sub$ratio == 0), ]
  
  # convert to log10 ratio values
  log_ratio <- log10(LICAGE_ctss_sub$ratio) 
  
  # convert CTSSnorm to GRanges
  LICAGE_ctss_sub.gr <- GRanges(seqnames = LICAGE_ctss_sub$chr, 
                                strand = LICAGE_ctss_sub$strand,
                                ranges = IRanges(start = LICAGE_ctss_sub$pos, width = 1),
                                tpm_sample = LICAGE_ctss_sub[, 4],
                                tpm_control = LICAGE_ctss_sub[, 5],
                                ratio = LICAGE_ctss_sub$ratio,
                                log10ratio = log_ratio)
  seqinfo(LICAGE_ctss_sub.gr) <- seqinfo(domTSS.grl[[1]])
  
  return(LICAGE_ctss_sub.gr)
}
### ---------------------------------------###

# calculate tpm ratios for all samples vs BY4741 standard
LICAGE_ctss_sub.grl <- lapply(LICAGE_ctss_sub.l, CTSStpmRatio)
names(LICAGE_ctss_sub.grl)  <- all_samples

# plot ratios for each sample in a heatmap centered on dominant TSS
library(heatmaps)
coords <- c(-30, 30)
windows_domTSS <- promoters(domTSS.grl[["BY4741"]], -coords[1], coords[2])
heatmap.l <- list()
for (i in 1:length(LICAGE_ctss_sub.grl)) {
  heatmap.l[[i]] <- CoverageHeatmap(
    windows = windows_domTSS,
    track = LICAGE_ctss_sub.grl[[i]],
    coords = coords,
    weight = LICAGE_ctss_sub.grl[[i]]$log10ratio,
    label = paste0("Tpm ratio: BY4741 vs ", names(LICAGE_ctss_sub.grl[i]))
  )
}
# add list names
names(heatmap.l) <- names(LICAGE_ctss_sub.grl)

# create control CTSS coverage heatmap
control_heatmap <- CoverageHeatmap(
  windows = windows_domTSS,
  track = control_CTSStpm.gr,
  coords = coords,
  weight = control_CTSStpm.gr$log10tpm,
  label = "CTSS coverage in BY4741"
) 

# extract raw images for reordering
raw_image.l <- list()
control_raw_image.l <- list()
order.l <- list() # saving order for reordering heatmaps with missing CTSS

for (i in 1:length(heatmap.l)) {
  raw_image <- 0
  raw_image <- image(heatmap.l[[i]])
  order <- order(raw_image[, 31], decreasing = TRUE)
  raw_image <- raw_image[order, ]
  raw_image.l[[i]] <- raw_image
  control_raw_image <- image(control_heatmap)
  control_raw_image <- control_raw_image[order, ]
  control_raw_image.l[[i]] <- control_raw_image
  order.l[[i]] <- order
}
# add names to ordered raw images
names(raw_image.l) <- names(LICAGE_ctss_sub.grl)
names(control_raw_image.l) <- names(LICAGE_ctss_sub.grl)
names(order.l) <- names(LICAGE_ctss_sub.grl)

# create heatmap list object
heatmap_order.l <- list()
for (i in 1:length(raw_image.l)) {
  heatmap_order.l[[i]] = Heatmap(
    raw_image.l[[i]],
    coords = coords,
    label = paste0("CTSS tpm ratios: BY4741 vs ", names(LICAGE_ctss_sub.grl[i])),
    scale = scale(heatmap.l[[i]])
  )
}

# create control heatmap list object
heatmap_control_order.l <- list()
for (i in 1:length(control_raw_image.l)) {
  heatmap_control_order.l[[i]] = Heatmap(
    control_raw_image.l[[i]],
    coords = coords,
    label = paste0("BY4741 CTSS cov ordered by ", names(LICAGE_ctss_sub.grl[i])),
    scale = scale(control_heatmap)
  )
}

# add names to the heatmap list object
names(heatmap_order.l) <- names(LICAGE_ctss_sub.grl)
names(heatmap_control_order.l) <- names(LICAGE_ctss_sub.grl)

# plot the heatmap list - each heatmap separately
for (i in 1:length(heatmap_order.l)) {
  png(paste0("results/CTSS_ratio_analysis/tpm_ratio_BY4741_", names(heatmap_order.l)[[i]], ".png"), width = 600, height = 800)
  plotHeatmapList(heatmap_order.l[[i]], color = c("black", "white", "red"), legend = TRUE, cex.label = 1, legend.width = 0.1)
  dev.off()
  png(paste0("results/CTSS_ratio_analysis/control_cov_", names(heatmap_order.l)[[i]], ".png"), width = 600, height = 800)
  plotHeatmapList(heatmap_control_order.l[[i]], color = c("white", "red", "black"), legend = TRUE, cex.label = 1, legend.width = 0.1)
  dev.off()
  
}

# cluster by tpm - form 4 quantiles - all samples
# add bin/quartile information to tpm.dominant_ctss
library(dplyr)

# calculate quantiles
tpm.df <- data.frame(tpm = domTSS.grl[["BY4741"]]$tpm.dominant_ctss)
tpm.df <- tpm.df %>% mutate(tpm_quantile = ntile(tpm, 4))

# attach quantile information as mcols to GRanges object
mcols(domTSS.grl[["BY4741"]])$tpm_quantile <- tpm.df$tpm_quantile

# plot heatmaps clustered per TPM quantiles (they are preordered per ratio of dominant TSS)
raw_image_clust.l <- list()
# order tpm.df and create a separate ordered df per each sample
tpm.dfl <- list()
for (i in 1:length(order.l)) {
  tpm.dfl[[i]] <- tpm.df[order.l[[i]], ]
}
names(tpm.dfl) <- names(order.l)

for (i in 1:length(heatmap_order.l)) {
  raw_image_clust <- 0
  raw_image_clust <- image(heatmap_order.l[[i]])
  raw_image_clust <- raw_image_clust[order(tpm.dfl[[i]]$tpm_quantile), ]
  raw_image_clust.l[[i]]<- raw_image_clust
}
# add names to ordered raw images
names(raw_image_clust.l) <- names(LICAGE_ctss_sub.grl)

# create heatmap object list
heatmap_order_clust.l <- list()
for(i in 1:length(raw_image_clust.l)) {
  heatmap_order_clust.l[[i]] = Heatmap(
    raw_image_clust.l[[i]],
    coords = coords,
    label = paste0("CTSS tpm ratios: BY4741 vs ", names(raw_image_clust.l[i])),
    scale = c(-2.5, 2.5)
  )
}
# add names to heatmap list
names(heatmap_order_clust.l) <- names(heatmap_order.l)

# plot control heatmaps clustered per TPM quantiles
control_raw_image_clust.l <- list()
for(i in 1:length(heatmap_control_order.l)) {
  control_raw_image_clust <- 0
  control_raw_image_clust <- image(heatmap_control_order.l[[i]])
  control_raw_image_clust <- control_raw_image_clust[order(tpm.dfl[[i]]$tpm_quantile), ]
  control_raw_image_clust.l[[i]] <- control_raw_image_clust
}
# add names to ordered raw images
names(control_raw_image_clust.l) <- names(heatmap_control_order.l)

# create heatmap control cluster list
heatmap_control_clust.l <- list()
for (i in 1:length(control_raw_image_clust.l)){
  heatmap_control_clust.l[[i]] = Heatmap(
    control_raw_image_clust.l[[i]],
    coords = coords,
    label = paste0("BY4741 CTSS cov ordered by  ", names(raw_image_clust.l[i])),
    scale = scale(heatmap_control_order.l[[i]])
  )
}
# add names to control heatmap list
names(heatmap_control_clust.l) <- names(control_raw_image_clust.l)

# plot the heatmap cluster list
for (i in 1:length(heatmap_order_clust.l)) {
  pdf(paste0("results/CTSS_ratio_analysis/tpm_ratio_BY4741_quantile_", names(heatmap_order_clust.l[i]), ".pdf"), width = 6, height = 8)
  plotHeatmapList(heatmap_order_clust.l[[i]], 
                  color = c("black", "white", "red"), 
                  legend = TRUE, 
                  cex.label = 1, 
                  legend.width = 0.15,
                  legend.pos="l",
                  partition = c(sum(tpm.dfl[[i]]$tpm_quantile == 1), 
                                sum(tpm.dfl[[i]]$tpm_quantile == 2),
                                sum(tpm.dfl[[i]]$tpm_quantile == 3),
                                sum(tpm.dfl[[i]]$tpm_quantile == 4)),
                  partition.legend = FALSE,
                  partition.lines = TRUE)
  dev.off()
}

# plot the control heatmap cluster list
for (i in 1:length(heatmap_control_clust.l)) {
  pdf(paste0("results/CTSS_ratio_analysis/control_cov_quantile", names(heatmap_control_clust.l[i]), ".pdf"), width = 6, height = 8)
  plotHeatmapList(heatmap_control_clust.l[[i]], 
                  color = c("white", "red", "black"), 
                  legend = TRUE, 
                  cex.label = 1, 
                  legend.width = 0.15,
                  legend.pos="l",
                  partition = c(sum(tpm.dfl[[i]]$tpm_quantile == 1), 
                                sum(tpm.dfl[[i]]$tpm_quantile == 2),
                                sum(tpm.dfl[[i]]$tpm_quantile == 3),
                                sum(tpm.dfl[[i]]$tpm_quantile == 4)),
                  partition.legend = FALSE,
                  partition.lines = TRUE)
  dev.off()
}

# Positional analysis of missing CTSS
# extract CTSS dataframe with normalized TPM per CTSS per sample
LICAGE_ctss <- CTSSnormalizedTpm(LICAGE)

# extract names of samples and make a list of cross-comparisons with controls
all_samples <- colnames(LICAGE_ctss[, 4:12])
control <- "BY4741"

LICAGE_ctss_sub.l <- list()
for (i in 1:length(all_samples)) {
  LICAGE_ctss_sub.l[[i]] <- LICAGE_ctss[, c("chr", "pos", "strand", all_samples[[i]], "BY4741")]   
}

### --- function to extract missing CTSSS ---###  

missingCTSS <- function(LICAGE_ctss_sub) {
  
  # extract CTSS examples were there is 0 TPM in one sample
  LICAGE_ctss_sub <- LICAGE_ctss_sub[(LICAGE_ctss_sub[, 4] == 0 & LICAGE_ctss_sub[, 5] != 0) | (LICAGE_ctss_sub[, 4] != 0 & LICAGE_ctss_sub[, 5] == 0), ]
  
  # assign -1 TPM value for CTSS missing in sample vs control, and +TPM value for CTSS in sample missing in control
  # create a penalty column
  LICAGE_ctss_sub$penalty <- 0
  LICAGE_ctss_sub$log10penalty <- 0
  
  # add penalties
  LICAGE_ctss_sub$penalty[which(0 == LICAGE_ctss_sub[, 4])] <- -LICAGE_ctss_sub[which(0 == LICAGE_ctss_sub[, 4]), 5]
  LICAGE_ctss_sub$penalty[which(0 == LICAGE_ctss_sub[, 5])] <- LICAGE_ctss_sub[which(0 == LICAGE_ctss_sub[, 5]), 4]
  
  # add log10(TPM) penalties
  LICAGE_ctss_sub$log10penalty[which(0 == LICAGE_ctss_sub[, 4])] <- log10(LICAGE_ctss_sub[which(0 == LICAGE_ctss_sub[, 4]), 5] + 1)
  LICAGE_ctss_sub$log10penalty[which(0 == LICAGE_ctss_sub[, 5])] <- -log10(LICAGE_ctss_sub[which(0 == LICAGE_ctss_sub[, 5]), 4] + 1)
  
  # convert to GRanges
  LICAGE_ctss_sub.gr <- GRanges(seqnames = LICAGE_ctss_sub$chr, 
                                strand = LICAGE_ctss_sub$strand,
                                ranges = IRanges(start = LICAGE_ctss_sub$pos, width = 1),
                                tpm_sample = LICAGE_ctss_sub[, 4],
                                tpm_control = LICAGE_ctss_sub[, 5],
                                penalty = LICAGE_ctss_sub$penalty,
                                log10penalty = LICAGE_ctss_sub$log10penalty)
  seqinfo(LICAGE_ctss_sub.gr) <- seqinfo(domTSS.grl[[1]])
  
  return(LICAGE_ctss_sub.gr)
}

### ---------------------------------------###
# calculate penalties for all samples
missed_CTSS.grl <- lapply(LICAGE_ctss_sub.l, missingCTSS)
names(missed_CTSS.grl)  <- all_samples

# plot ratios for each sample in a heatmap centered on dominant TSS
library(heatmaps)
coords <- c(-30, 30)
windows_domTSS <- promoters(domTSS.grl[["BY4741"]], -coords[1], coords[2])
heatmap.l <- list()
for (i in 1:length(LICAGE_ctss_sub.grl)) {
  heatmap.l[[i]] <- CoverageHeatmap(
    windows = windows_domTSS,
    track = missed_CTSS.grl[[i]],
    coords = coords,
    weight = missed_CTSS.grl[[i]]$log10penalty,
    label = paste0("Tpm ratio: BY4741 vs ", names(missed_CTSS.grl[i]))
  )
}

# add list names
names(heatmap.l) <- names(LICAGE_ctss_sub.grl)

# extract raw image for reordering - order by ordering in the ratio analysis
raw_image.l <- list()
for (i in 1:length(heatmap.l)) {
  raw_image <- image(heatmap.l[[i]])
  raw_image <- raw_image[order.l[[i]], ]
  raw_image.l[[i]]<- raw_image
}
# add names to ordered raw images
names(raw_image.l) <- names(missed_CTSS.grl)

# create heatmap list object
heatmap_order.l <- list()
for (i in 1:length(raw_image.l)) {
  heatmap_order.l[[i]] = Heatmap(
    raw_image.l[[i]],
    coords = coords,
    label = paste0("Tpm ratio: BY4741 vs ", names(missed_CTSS.grl[i])),
    scale = c(-2, 2)
  )
}
# add names to the heatmap list object
names(heatmap_order.l) <- names(missed_CTSS.grl)

# plot the heatmap list - each heatmap separately
for (i in 1:length(heatmap_order.l)) {
  png(paste0("results/CTSS_ratio_analysis/missingCTSS_", names(heatmap_order.l)[[i]], ".png"), width = 600, height = 800)
  plotHeatmapList(heatmap_order.l[[i]], color = c("black", "white", "red"), legend = TRUE, cex.label = 1, legend.width = 0.1)
  dev.off()
}

# plot heatmaps clustered per TPM quantiles (they are preordered per ratio of dominant TSS)
raw_image_clust.l <- list()

for (i in 1:length(heatmap_order.l)) {
  raw_image_clust <- 0
  raw_image_clust <- image(heatmap_order.l[[i]])
  raw_image_clust <- raw_image_clust[order(tpm.dfl[[i]]$tpm_quantile), ]
  raw_image_clust.l[[i]]<- raw_image_clust
}
# add names to ordered raw images
names(raw_image_clust.l) <- names(LICAGE_ctss_sub.grl)

# create heatmap object list
heatmap_order_clust.l <- list()
for(i in 1:length(raw_image_clust.l)) {
  heatmap_order_clust.l[[i]] = Heatmap(
    raw_image_clust.l[[i]],
    coords = coords,
    label = paste0("Missing CTSSs: BY4741 vs ", names(raw_image_clust.l[i])),
    scale = scale(heatmap_order.l[[i]])
  )
}
# add names to heatmap list
names(heatmap_order_clust.l) <- names(heatmap_order.l)

# plot the heatmap cluster list
for (i in 1:length(heatmap_order_clust.l)) {
  pdf(paste0("results/CTSS_ratio_analysis/missingCTSS_quantile_", names(heatmap_order_clust.l[i]), ".pdf"), width = 6, height = 8)
  plotHeatmapList(heatmap_order_clust.l[[i]], 
                  color = c( "dodgerblue4", "white", "darkviolet"), 
                  legend = TRUE, 
                  cex.label = 1, 
                  legend.width = 0.15,
                  legend.pos="l",
                  partition = c(sum(tpm.df$tpm_quantile == 1), 
                                sum(tpm.df$tpm_quantile == 2),
                                sum(tpm.df$tpm_quantile == 3),
                                sum(tpm.df$tpm_quantile == 4)),
                  partition.legend = FALSE,
                  partition.lines = TRUE)
  dev.off()
}
