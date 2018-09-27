#/*==========================================================================#*/
#' ## Supplemental Figure S17 - IQ-width distribution
#'
#/*==========================================================================#*/

# load RDS object centered on dominant TSS (merged lanes and replicates)
domTSS_E14slic_mergLandR.grl <- readRDS(file = "../intermediate_data/domTSS_E14slic_mergLandR_grl.RDS")

# separate sharp and broad promoters
threshold <- 1:100
samples <- c("SLIC 5 ng", "SLIC 10 ng", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi 5 ug")

# get the number of sharp and broad promoters for each sample and threshold
prom_width_df.l <- list()
for(i in 1:length(domTSS_E14_win.grl)) {
  nr_sharp <- sapply(threshold, function(x) sum(domTSS_E14slic_mergLandR.grl[[i]]$interquantile_width <= x))
  nr_broad <- sapply(threshold, function(x) sum(domTSS_E14slic_mergLandR.grl[[i]]$interquantile_width > x))
  prom_width_df.l[[i]] <- data.frame(nr_sharp = nr_sharp, nr_broad = nr_broad, threshold = threshold, sample = samples[i])
}
names(prom_width_df.l) <- names(domTSS_E14slic_mergLandR.grl)

# convert dataframe list to concatenated dataframe
prom_width_df <- do.call(rbind, prom_width_df.l)
rownames(prom_width_df) <- 1:nrow(prom_width_df)

# prepare dataframe for plotting
library(tidyr)

prom_width_df.gg <- gather(prom_width_df, 
                           ...= colnames(prom_width_df[,1:2]), 
                           key = "type",
                           value = "number")

# visualization of promoter width/number of sharp vs broad
library(viridis)
library(ggplot2)
p <- ggplot(data = prom_width_df.gg, aes(x = threshold, y = number, fill = type), alpha = 0.8) +
  geom_bar(colour = "black", size = 0.1, stat = "identity") +
  scale_fill_manual(values = c("royalblue4", "firebrick"), name = NULL, labels = c("broad", "sharp")) +
  theme_bw() +
  theme(text = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  geom_vline(xintercept = 3, lty = 2, colour = "gray66", size = 0.8) +
  scale_y_continuous(name = "Number of tag clusters") +
  scale_x_continuous(name = "Threshold (tag cluster interquantile width)")

p <- p + facet_wrap(~ sample, ncol = 2, scales = "free")

pdf("../final_figures/sharp_broad_distribution.pdf", height = 12, width = 10)
print(p)
dev.off()