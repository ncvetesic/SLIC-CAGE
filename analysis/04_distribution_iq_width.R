#/*==========================================================================#*/
#' ## Figures 2B and F - distribution of tag cluster/promoter interquantile widths
#' 
#/*==========================================================================#*/

# selected samples for main paper
samples <- c("Y1ng_carrier", "Y5ng_carrier", "Y10ng_carrier1", "BY4741")
tc_selected <- tc_anno.grl[samples]

# add alphabet letters to get the proper plotting order
names <- c("SLIC 1 ng", "SLIC 5 ng", "SLIC 10 ng", "nAnTi 5 ug")

# extract interquantile widths
iq.l <- lapply(tc_selected, function(x) data.frame(iq_width = x$interquantile_width))

# add names to as a column to each sample
for (i in 1:length(names)) {
  iq.l[[i]]$sample <- rep(names[[i]], nrow(iq.l[[i]]))
}

# combine to 1 dataframe
iq.df <- do.call(rbind, iq.l)

# set levels for plotting - sample levels
iq.df$sample <- factor(iq.df$sample, levels = names)

# plotting
library(ggplot2)
p <- ggplot(iq.df, aes(iq_width)) +
  geom_histogram(aes(x = iq_width, y = ..ncount..* 100), binwidth = 1, 
                 fill = "gray60", col = "black", size = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 16, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black", fill = "gray87")) +
  labs(x = "Interquantile width", y = "Percentage") +
  coord_equal(ratio = 1) +
  xlim(0, 150)

pdf(file = "final_figures/IQ_width_sel_sc.pdf", height = 6, width = 6)
p + facet_wrap(~ sample, ncol = 2)
dev.off()


# plot all interquantile widths for supplementary
iq_all.l <- lapply(tc_anno.grl, function(x) data.frame(iq_width = x$interquantile_width))

# set names for plotting
names <- c("SLIC 1 ng", "SLIC 2 ng", "SLIC 5 ng", 
           "SLIC 10 ng r1", "SLIC 10 ng r2", "SLIC 25 ng", 
           "SLIC 50 ng", "SLIC 100 ng", "nAnTi PCR",
           "nAnTi 5 ug")

names(iq_all.l) <- names

# add names to as a column to each sample
for (i in 1:length(names)) {
  iq_all.l[[i]]$sample <- rep(names[[i]], nrow(iq_all.l[[i]]))
}

# combine to 1 dataframe
iq_all.df <- do.call(rbind, iq_all.l)

# set levels for plotting - sample levels
iq_all.df$sample <- factor(iq_all.df$sample, levels = names)

# plot iq-widths for all samples
# plotting
library(ggplot2)
p <- ggplot(iq_all.df, aes(iq_width)) +
  geom_histogram(aes(x = iq_width, y = ..ncount..* 100), binwidth = 1, 
                 fill = "gray60", col = "black", size = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 16, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black", fill = "gray87")) +
  labs(x = "Interquantile width", y = "Percentage") +
  coord_equal(ratio = 1) +
  xlim(0, 150)

pdf(file = "final_figures/IQ_width_all_slic.pdf", height = 4, width = 15)
p + facet_wrap(~ sample, ncol = 5)
dev.off()