#/*==========================================================================#*/
#' ## Supplemental Figure S15 - heatmap similarity - Jaccard index
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

# compare 10 ng and nanti heatmaps - extract correlation

# ---- compare TA heatmaps ---- #
library(philentropy)


jaccard_TA <- distance(rbind(as.vector(TA_images.l[["E14_10ng"]]), as.vector(TA_images.l[["E14_nanti"]])) , method = "jaccard")

# randomise TA images columnwise before vectorising and plot jaccard index
jaccard_random_TA <- vector()
for (i in 1:10000) {
  print(i)
  jaccard_random_TA[i] <- distance(rbind(as.vector(TA_images.l[["E14_10ng"]][, sample(ncol(TA_images.l[["E14_10ng"]]))]), as.vector(TA_images.l[["E14_nanti"]][, sample(ncol(TA_images.l[["E14_nanti"]]))])), method = "jaccard")
}
jaccard_random_TA <- as.data.frame(jaccard_random_TA)

# plot jaccard_distances as histogram overlaid with kernel density curve
library(ggplot2)
p <- ggplot(jaccard_random_TA, aes(x = jaccard_random_TA)) +
  geom_histogram(aes(y=..density.. ),
                 binwidth=.001,
                 colour="black", fill="white", size = 0.5) +
  geom_density(alpha=.5, fill="#FF6666", size = 0.7) +
  geom_segment(x = jaccard_TA, xend = jaccard_TA, y = 0, yend = 120, colour = "dodgerblue2", linetype = "longdash", size = 0.7) +
  geom_segment(x = quantile(jaccard_random_TA$jaccard_random_TA, 0.01), xend = quantile(jaccard_random_TA$jaccard_random_TA, 0.01),y = 0, yend = 120, colour = "firebrick2", linetype = "longdash", size = 0.7) +
  scale_x_continuous(limits = c(0.04, 0.2), name = "Jaccard distance") +
  scale_y_continuous(name = "Density") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  annotate("text", x = 0.07, y = 130, label = "TA, pvalue <<< 0.01", size = 6)
pdf("../final_figures/jaccard_TA_10ng_vs_nanti.pdf", height = 4, width = 6)
print(p)
dev.off()

# ---- compare GC heatmaps ---- #
library(philentropy)

jaccard_GC <- distance(rbind(as.vector(GC_images.l[["E14_10ng"]]), as.vector(GC_images.l[["E14_nanti"]])) , method = "jaccard")

# randomise GC images columnwise before vectorising and plot jaccard index
jaccard_random_GC <- vector()
for (i in 1:10000) {
  print(i)
  jaccard_random_GC[i] <- distance(rbind(as.vector(GC_images.l[["E14_10ng"]][, sample(ncol(GC_images.l[["E14_10ng"]]))]), as.vector(GC_images.l[["E14_nanti"]][, sample(ncol(GC_images.l[["E14_nanti"]]))])), method = "jaccard")
}

jaccard_random_GC <- as.data.frame(jaccard_random_GC)

# plot jaccard_distances as histogram overlaid with kernel density curve
library(ggplot2)
p <- ggplot(jaccard_random_GC, aes(x = jaccard_random_GC)) +
  geom_histogram(aes(y=..density.. ),
                 binwidth=.001,
                 colour="black", fill="white", size = 0.5) +
  geom_density(alpha=.5, fill="#FF6666", size = 0.7) +
  geom_segment(x = jaccard_GC, xend = jaccard_GC, y = 0, yend = 180, colour = "dodgerblue2", linetype = "longdash", size = 0.7) +
  geom_segment(x = quantile(jaccard_random_GC$jaccard_random_GC, 0.01), xend = quantile(jaccard_random_GC$jaccard_random_GC, 0.01),y = 0, yend = 180, colour = "firebrick2", linetype = "longdash", size = 0.7) +
  scale_x_continuous(limits = c(0.01, 0.1), name = "Jaccard distance") +
  scale_y_continuous(name = "Density") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  annotate("text", x = 0.03, y = 200, label = "GC, pvalue <<< 0.01", size = 6)
pdf("../final_figures/jaccard_GC_10ng_vs_nanti.pdf", height = 4, width = 6)
print(p)
dev.off()
