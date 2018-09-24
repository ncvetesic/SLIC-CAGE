#/*==========================================================================#*/
#' ## Figures 1E, I
#' 
#/*==========================================================================#*/

library(Gviz)
library(GenomicFeatures)

# prepare CTSS data for plotting - separate samples, filter everything <1TPM and convert to GRangesList
# extract CTSS dataframe with normalized TPM per CTSS per sample
LICAGE_ctss <- CTSSnormalizedTpm(LICAGE)

# change chr names to UCSC
UCSC_chr <- fetchExtendedChromInfoFromUCSC("sacCer3",
                                           goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
                                           quiet=FALSE)
UCSC_chr_names <- UCSC_chr$UCSC_seqlevel

# list of new chr names 
chr_naming.df <- data.frame(chr_NCBI = unique(LICAGE_ctss$chr), UCSC_chr_names = UCSC_chr_names)

# attach new chr names to CTSS dataframe
LICAGE_ctss_renam <- merge(LICAGE_ctss, chr_naming.df, by.x = "chr", by.y = "chr_NCBI")

# set strand to the same and if -strand change TPM values to -TPM
LICAGE_ctss_renam_noStrand <- LICAGE_ctss_renam

idx_min <- LICAGE_ctss_renam_noStrand$strand %in% "-"
LICAGE_ctss_renam_noStrand[idx_min, -c(1:3, ncol(LICAGE_ctss_renam_noStrand))] = - LICAGE_ctss_renam_noStrand[idx_min, -c(1:3, ncol(LICAGE_ctss_renam_noStrand))]

LICAGE_ctss_renam_noStrand$strand <- "*"

# select samples to correspond correlation plots
samples <- c("Y5ng_carrier", "Y10ng_carrier1", "Y10ng_carrier2", "BY4741")

CTSS_select_all.grl <- list()

for (i in 1:length(samples)) {
  CTSS_select_all.grl[[i]] <- GRanges(seqnames = LICAGE_ctss_renam_noStrand$UCSC_chr_names,
                                      ranges = IRanges(start = LICAGE_ctss_renam_noStrand$pos,
                                                       end = LICAGE_ctss_renam_noStrand$pos),
                                      strand = LICAGE_ctss_renam_noStrand$strand)
  
  TPM <- LICAGE_ctss_renam_noStrand[, samples[[i]]]
  strand <- LICAGE_ctss_renam$strand
  
  # add strand info as mcols - TPM values are sample specific
  mcols(CTSS_select_all.grl[[i]]) <- data.frame(TPM = TPM, strand = strand)
}

names(CTSS_select_all.grl) <- samples

# create selected range for tracks
range.gr <- GRanges(seqnames = "chrXII", 
                    ranges = IRanges(start = 1026799, end = 1032009),
                    strand = "*")

subset.grl <- lapply(CTSS_select_all.grl, function(x) subsetByOverlaps(x, range.gr))
subset.grl <- GenomicRangesList(subset.grl)



# create CAGE tracks as datatrack
SLIC_sub_track.l <- list()
names <- c("5 ng SLIC-CAGE", "10 ng SLIC-CAGE r1", "10 ng SLIC-CAGE r2", "nAnTi 5 \u03bcg")

# create track as histogram type in gviz
for (i in 1:length(subset.grl)) {
  SLIC_sub_track.l[[i]] <-  DataTrack(subset.grl[[i]], 
                                      genome = "sacCer3", 
                                      showTitle = FALSE,
                                      chr = "chrXII", 
                                      cex = 1,
                                      col.baseline = "black",
                                      col = "black",
                                      ylim = c(-60, 85),
                                      type = "histogram",
                                      baseline = 0,
                                      fill.histogram = "gray60",
                                      col = "black",
                                      col.histogram="black",
                                      fontsize = 16,
                                      background.title="white",
                                      col.axis="black",
                                      col.title="black",
                                      fontfamily="Helvetica",
                                      lwd.mountain = 0.25,
                                      cex.axis = 1,
                                      lwd = 1,
                                      lwd.mountain = 0.5)
}
names(SLIC_sub_track.l) <- names

# save SLIC_sub_tracl.l object
saveRDS(object = SLIC_sub_track.l, file = "intermediate/SLIC_sub_track.RDS")

# zoomed in CTSS region to put as inset - region 1027920-1027960
# create selected range for tracks
range_zoom.gr <- GRanges(seqnames = "chrXII", 
                         ranges = IRanges(start = 1027920, end = 1027975),
                         strand = "*")

subset_zoom.grl <- lapply(CTSS_select_all.grl, function(x) subsetByOverlaps(x, range_zoom.gr))
subset_zoom.grl <- GenomicRangesList(subset_zoom.grl)

# create zoomed track as histogram type in gviz
SLIC_sub_zoom_track.l <- list()
for (i in 1:length(subset_zoom.grl)) {
  SLIC_sub_zoom_track.l[[i]] <-  DataTrack(subset_zoom.grl[[i]], 
                                           genome = "sacCer3", 
                                           showTitle = FALSE,
                                           chr = "chrXII", 
                                           cex = 1,
                                           col.baseline = "black",
                                           ylim = c(-70, 0),
                                           type = "histogram",
                                           baseline = 0,
                                           fill.histogram = "gray60",
                                           col.histogram = "black",
                                           fontsize = 18,
                                           background.title="white",
                                           col.axis="black",
                                           col.title="black",
                                           fontfamily="Helvetica",
                                           cex.axis = 1,
                                           lwd = 1.5)
}

names(SLIC_sub_zoom_track.l) <- names

# save SLIC_sub_zoom_track.l object
saveRDS(object = SLIC_sub_zoom_track.l, file = "intermediate/SLIC_sub_zoom_track.RDS")


# new zoomed in region for comparison with  nanoCAGE

# zoomed in CTSS region to put as inset - region 1028750-1028860
# create selected range for tracks
range_zoom.gr <- GRanges(seqnames = "chrXII", 
                         ranges = IRanges(start = 1028750, end = 1028860),
                         strand = "*")

subset_4nc_zoom.grl <- lapply(CTSS_select_all.grl, function(x) subsetByOverlaps(x, range_zoom.gr))
subset_4nc_zoom.grl <- GenomicRangesList(subset_4nc_zoom.grl)

# create zoomed track as histogram type in gviz
sub_4nc_zoom_track.l <- list()
for (i in 1:length(subset_4nc_zoom.grl)) {
  sub_4nc_zoom_track.l[[i]] <-  DataTrack(subset_4nc_zoom.grl[[i]], 
                                          genome = "sacCer3", 
                                          showTitle = FALSE,
                                          chr = "chrXII", 
                                          cex = 1,
                                          ylim = c(0, 100),
                                          col.baseline = "black",
                                          type = "histogram",
                                          baseline = 0,
                                          fill.histogram = "gray60",
                                          col.histogram = "black",
                                          fontsize = 16,
                                          background.title="white",
                                          col.axis="black",
                                          col.title="black",
                                          fontfamily="Helvetica",
                                          cex.axis = 1,
                                          lwd = 1.5)
}

names(sub_4nc_zoom_track.l) <- names

# save nc_sub_zoom_track.l object
saveRDS(object = sub_4nc_zoom_track.l, file = "intermediate/sub_4nc_zoom_track.RDS")

# create genome axis
genomeAxis <- GenomeAxisTrack(range = range.gr, labelPos = "alternating", 
                              add53 = FALSE, add35 = FALSE, cex = 1,
                              col.range= "white",
                              fill.range="white",
                              fontcolor = "black",
                              fontsize = 16,
                              frame = FALSE,
                              font.family = "Helvetica",
                              stackHeight = 0.05)

# create ideogram track
itrack <- IdeogramTrack(genome = "sacCer3", chromosome = UCSC_chr_names)

# create gene track
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

gene_track <- GeneRegionTrack(txdb, chromosome = "chrXII", start = start(range.gr), end = end(range.gr),
                              showId = TRUE, geneSymbols = TRUE,
                              col.sampleNames = "black",
                              background.title = "white",
                              col.axis = "black",
                              col.title = "black",
                              col = "black",
                              fill = "darkgray",
                              fontfamily = "Helvetica",
                              fontsize = 16,
                              cex = 1,
                              shape = "arrow",
                              lwd = 1,
                              transcriptAnnotation = "symbol",
                              stackHeight = 0.5, thinBoxFeature = "UTR", 
                              arrowHeadWidth = 3,
                              arrowHeadMaxWidth = 5,
                              geneSymbol = TRUE,
                              fontsize.group = 22,
                              just.group="left",
                              min.distance = 5,
                              stacking = "squish",
                              fontcolor.group = "black")


pdf(file = "final_figures/sc_slic_bg.pdf", height = 7, width = 11)
plotTracks(trackList = c(genomeAxis,
                         gene_track,
                         SLIC_sub_track.l))
dev.off()

pdf(file = "final_figures/sc_slic_bg_zoom.pdf", height = 7, width = 4)
plotTracks(trackList = SLIC_sub_zoom_track.l)
dev.off()