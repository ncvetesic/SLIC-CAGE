#/*==========================================================================#*/
#' ## Supplemental Figure S23 - consensus cluster level correlation
#'
#/*==========================================================================#*/

### Plot selected correlations for Supplementary - consensus cluster level
tag.count <- as.data.frame(myCAGEset@consensusClustersTpmMatrix)

# plot corr - just 5 ng, 10 ng1, 10 ng2 50 and 500 and nanti 


samples <- c("ncSC_5ng", "ncSC_10ng_r1",  "ncSC_10ng_r2", "ncSC_50ng_r1", "ncSC_500ng_r1")
control <- "sc_nanti_r1"

filtr.data.l <- list()
filtr.data.log.l <- list()
cor.l <- list()

#plot log scale but calculate correlation non non-log scale
for (i in 1: length(samples)) {
  filtr.data.l[[i]] <- tag.count[ , c(samples[[i]], control)]
  filtr.data.l[[i]] <- filtr.data.l[[i]][(filtr.data.l[[i]][, 1] >= 1 | filtr.data.l[[i]][, 2] >= 1), ]
  cor.l[[i]] <- round(cor(filtr.data.l[[i]]), digits = 3)
  filtr.data.log.l[[i]] <- log10(filtr.data.l[[i]] + 1)
}
names(filtr.data.l) <- samples
names(filtr.data.log.l) <- samples

# add comparison of 10ng rep1 and 10ng rep2
filtr.data.l[[6]] <- tag.count[ , c("ncSC_10ng_r1", "ncSC_10ng_r2")]
filtr.data.l[[6]] <- filtr.data.l[[6]][(filtr.data.l[[6]][, 1] >= 1 | filtr.data.l[[6]][, 2] >= 1), ]
cor.l[[6]] <- round(cor(filtr.data.l[[6]]), digits = 3)
filtr.data.log.l[[6]] <- log10(filtr.data.l[[6]] + 1)


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

labels_y <- c("5 ng nanoCAGE", "10 ng nanoCAGE r1", "10 ng nanoCAGE r2", 
              "50 ng nanoCAGE", "500 ng nanoCAGE", "10 ng nanoCAGE r1")
labels_x <- c(expression(paste("5 ", mu,"g nAnT-iCAGE")), 
              expression(paste("5 ", mu,"g nAnT-iCAGE")), 
              expression(paste("5 ", mu,"g nAnT-iCAGE")), 
              expression(paste("5 ", mu,"g nAnT-iCAGE")),
              expression(paste("5 ", mu,"g nAnT-iCAGE")),
              "10 ng nanoCAGE r2")
samples <- c("ncSC_5ng", "ncSC_10ng_r1",  "ncSC_10ng_r2", "ncSC_50ng", "ncSC_500ng", "ncSC_10ng_r1_r2")
subtitles <- c("a\n", "b\n", "c\n", "d\n", "e\n", "f\n")

p.l <- list()

for(i in 1:length(samples)) {
  pdf(file = paste0("../final_figures/", samples[[i]], "cons_cluster.pdf"), width = 7, height = 7)
  p <- corrPlot(data = filtr.data.log.l[[i]], xlab = labels_x[[i]], ylab = labels_y[[i]], 
                cor = cor.l[[i]][1,2], sub = subtitles[[i]])
  print(p)
  p.l[[i]] <- p
  dev.off()
}


# plots, set up layout in R altogether 3by2 plot
library(cowplot)
pgrid <- plot_grid(plotlist = p.l, align = "h", ncol = 3)
save_plot("../final_figures/nanoCAGE_select_corr_3by2_cons_cluster.pdf", pgrid,
          ncol =3, 
          nrow =2, 
          base_aspect_ratio = 1)

# plots, set up layout in R altogether horizontal plot
pgrid <- plot_grid(plotlist = p.l, align = "h", ncol = 6)
save_plot("../final_figures/nanoCAGE_select_corr_horiz_cons_cluster.pdf", pgrid,
          ncol = 6, 
          nrow = 1, 
          base_aspect_ratio = 1)

