#/*==========================================================================#*/
#' ## Figures 1B-D, F-H
#' ## Supplementary Figures S1C-I
#/*==========================================================================#*/

tag.count <- LICAGE@normalizedTpmMatrix

# example correlation

# filter  each compared pair to have 1 TPM in at least one sample; plot log scale but calculate correlation on non-log scale
samples <- c("Y5ng_carrier", "Y10ng_carrier1",  "Y10ng_carrier2")
control <- "BY4741"

filtr.data.l <- list()
filtr.data.log.l <- list()
cor.l <- list()

for (i in 1: length(samples)) {
  filtr.data.l[[i]] <- tag.count[ , c(samples[[i]], control)]
  filtr.data.l[[i]] <- filtr.data.l[[i]][(filtr.data.l[[i]][, 1] >= 1 | filtr.data.l[[i]][, 2] >= 1), ]
  cor.l[[i]] <- round(cor(filtr.data.l[[i]]), digits = 3)
  filtr.data.log.l[[i]] <- log10(filtr.data.l[[i]] + 1)
}
names(filtr.data.l) <- samples
names(filtr.data.log.l) <- samples

# add comparison of 10ng rep1 and 10ng rep2
filtr.data.l[[4]] <- tag.count[ , c("Y10ng_carrier1", "Y10ng_carrier2")]
filtr.data.l[[4]] <- filtr.data.l[[4]][(filtr.data.l[[4]][, 1] >= 1 | filtr.data.l[[4]][, 2] >= 1), ]
cor.l[[4]] <- round(cor(filtr.data.l[[4]]), digits = 3)
filtr.data.log.l[[4]] <- log10(filtr.data.l[[4]] + 1)


# plots, for manual layout in inkscape
library(ggplot2)
corrPlot <- function(data, cor, xlab, ylab, sub) {
  p <- ggplot(data, aes(data[,2], y = data[,1])) +
    geom_hex(bins = 100, show.legend = FALSE) +
    scale_fill_gradient(low = "gray24", high = "gray56") +
    theme(text = element_text(size = 20,
                              family = "Helvetica"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.subtitle = element_text(size = 24, hjust = -0.125, vjust = 1, face="bold", color="black"),
          axis.text = element_text(size = 18, family = "Helvetica", color = "black"),
          legend.title = element_text(size = 16, family = "Helvetica", colour = "black"),
          legend.text = element_text(size = 14, family = "Helvetica", colour = "black")) +
    xlab(xlab) +
    ylab(ylab) +
    xlim(c(0, 4.25)) +
    ylim(c(0, 4.25)) +
    annotate("text", x = 0, y = 4, size = 6, hjust = 0, label = paste("R = ", cor)) +
    labs(subtitle = sub) +
    coord_fixed(ratio = 1)
  
  return(p)
}

labels_y <- c("5 ng SLIC-CAGE", "10 ng SLIC-CAGE r1", "10 ng SLIC-CAGE r2", "10 ng SLIC-CAGE r1")
labels_x <- c(expression(paste("5 ", mu,"g nAnT-iCAGE")), 
              expression(paste("5 ", mu,"g nAnT-iCAGE")), 
              expression(paste("5 ", mu,"g nAnT-iCAGE")), 
              "10 ng SLIC-CAGE r2")
samples <- c("Y5ng_carrier", "Y10ng_carrier1",  "Y10ng_carrier2", "Y10ng_carrier1_2")
subtitles <- c("b\n", "c\n", "d\n", "e\n")

p.l <- list()

for(i in 1:length(samples)) {
  pdf(file = paste0("final_figures/", samples[[i]], ".pdf"), width = 7, height = 7)
  p <- corrPlot(data = filtr.data.log.l[[i]], xlab = labels_x[[i]], ylab = labels_y[[i]], 
                cor = cor.l[[i]][1,2], sub = subtitles[[i]])
  print(p)
  p.l[[i]] <- p
  dev.off()
}

# plots, set up layout in R 4by4 plot
library(cowplot)
pgrid <- plot_grid(plotlist = p.l, align = "h", ncol = 2)
save_plot("final_figures/slicCAGE_select_corr_2by2.pdf", pgrid,
          ncol =2, 
          nrow =2, 
          base_aspect_ratio = 1)

# plots, set up layout in R horizontal plot
pgrid <- plot_grid(plotlist = p.l, align = "h", ncol = 4)
save_plot("final_figures/slicCAGE_select_corr_horiz.pdf", pgrid,
          ncol = 4, 
          nrow = 1, 
          base_aspect_ratio = 1)