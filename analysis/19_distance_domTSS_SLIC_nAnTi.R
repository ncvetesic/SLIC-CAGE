#/*==========================================================================#*/
#' ## Supplemental Figure S6 - distance of dominant TSS SLIC- vs nAnT-iCAGE
#'
#/*==========================================================================#*/

# Distance between dominant TSS in SLIC-CAGE vs nAnT-iCAGE
# load tc.grl object
tc.grl <- readRDS("data/tagClusters_list.RDS")

# separate samples and the control nAnTi sample
sample.grl <- tc.grl[-length(tc.grl)]

control.gr <- tc.grl[length(tc.grl)][[1]]

# find overlapping tag cluster -each sample vs control
overlap.l <- lapply(sample.grl, function(x) findOverlapPairs(x, control.gr))

# extract dominant TSS distances in the list
distance.l <- list()
for(i in 1:length(overlap.l)) {
  first.gr <- first(overlap.l[[i]])
  second.gr <- second(overlap.l[[i]])
  
  # extract positive strands
  positive_idx <- as.logical(strand(second.gr) == "+")
  
  distance_plus.l <- data.frame(distance = first.gr$dominant_ctss[positive_idx] - second.gr$dominant_ctss[positive_idx])
  distance_minus.l <- data.frame(distance = -(first.gr$dominant_ctss[!positive_idx] - second.gr$dominant_ctss[!positive_idx]))
  
  distance.l[[i]] <- rbind(distance_plus.l, distance_minus.l)
  distance.l[[i]] <- cbind(distance.l[[i]], sample = rep(names(overlap.l[i]), times = length(overlap.l[[i]])))
}

# convert to a dataframe
distance.df <- do.call(rbind, distance.l)

# ggplot distance
col = rep("dodgerblue", times = 9)

library(ggplot2)
p <- ggplot(distance.df, aes(x = distance,  fill = sample), alpha = 0.5) +
  geom_histogram (binwidth = 1, col = "black", lwd = 0.125) +
  scale_fill_manual("sample", values = col) +
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Frequency", x = NULL) + 
  coord_cartesian( xlim = c(-100, 100)) +
  guides(fill = FALSE)
p + facet_wrap(~ sample, ncol = 3)

pdf(file = "final_figures/distance_domTSS_slic.pdf", height = 6, width = 8)
print(p + facet_wrap(~ sample, ncol = 3))
dev.off()

# plot zoom-in
library(ggplot2)
p <- ggplot(distance.df, aes(x = distance,  fill = sample), alpha = 0.5) +
  geom_histogram (binwidth = 1, col = "black", lwd = 0.125) +
  scale_fill_manual("sample", values = col) +
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Frequency", x = NULL) + 
  coord_cartesian( xlim = c(-100, 100), ylim = c(0, 300)) +
  guides(fill = FALSE)
p + facet_wrap(~ sample, ncol = 3)

pdf(file = "final_figures/distance_domTSS_slic_zoom.pdf", height = 6, width = 8)
print(p + facet_wrap(~ sample, ncol = 3))
dev.off()


# extract number of 0 distance
dist_table <- lapply(distance.l, function(x) data.frame(table(x)))
dist_table_0 <- do.call(rbind, lapply(dist_table, function(x) x[x$distance == 0, ]))

# recalculate to proportion of matched domTSS
dist_table_0$prop <- dist_table_0$Freq/sapply(overlap.l, length)*100

# quantify how many are within +-10 bp
distance_10bp.l <- lapply(distance.l, function(x) x[abs(x$distance) <= 10, ])

distance_10bp.df <- sapply(distance_10bp.l, nrow)/sapply(overlap.l, length)*100
