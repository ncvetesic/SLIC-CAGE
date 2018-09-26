#/*==========================================================================#*/
#' ## Supplemental Figure S5 - ROC curves
#'
#/*==========================================================================#*/

##------------ ROC curves on individual CTSS level ------------## 
CTSSnorm <- CTSSnormalizedTpm(LICAGE)

# convert to dataframe and add columns TP(true positive) FP(false positive)
CTSSnorm <- as.data.frame(CTSSnorm)

# extract TPM values as pair sample vs nanti
samples <- sampleLabels(LICAGE)[-length(sampleLabels(LICAGE))]
tpm_pairs.dfl <- lapply(samples, function(x) {
  df <- CTSSnorm[, c(x, "BY4741")]
  df <- df[rowSums(df) > 0, ]
  return(df) })
names(tpm_pairs.dfl) <- samples

# calculate TP and FP per different TPM threshold
TP_FP_l <- list()
TP_FP.gg <- list()

for(j in 1:length(tpm_pairs.dfl)){
  threshold <- seq(0, 10, 0.05) 
  for (i in 1:length(threshold)) {
    tpm_filt <- tpm_pairs.dfl[[j]][tpm_pairs.dfl[[j]][, 1] >= threshold[i], ]
    TP_FP_l[[i]] <- data.frame("TP" = tpm_filt[,1] > 0 & tpm_filt[,2] > 0,
                               "FP" = tpm_filt[,1] > 0 & tpm_filt[,2] == 0)
    
  }
  TP_FP.gg[[j]] <- as.data.frame(do.call(rbind, lapply(TP_FP_l, colSums)))
  TP_FP.gg[[j]]$sample <- rep(names(tpm_pairs.dfl[j]), times = nrow(TP_FP.gg[[j]]))
}

# flatten for plotting
TP_FP.gg <- as.data.frame(do.call(rbind, TP_FP.gg))
TP_FP.gg$sample <- factor(TP_FP.gg$sample, levels = samples)

library(viridis)
col <- magma(10)

p <- ggplot(TP_FP.gg, aes(x = FP, y = TP, fill = sample)) +
  geom_line(aes(col = sample), size = 1, alpha = 0.7) +
  geom_point(shape = 21, size = 1, alpha = 0.7) +
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) +
  theme(text = element_text(size = 12,
                            family = "Helvetica"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12, family = "Helvetica", color = "black"),
        legend.title = element_text(size = 12, family = "Helvetica", colour = "black"),
        legend.text = element_text(size = 12, family = "Helvetica", colour = "black")) +
  xlab("False positive CTSSs") +
  ylab("True positive CTSSs") +
  xlim(c(0, 120000)) +
  ylim(c(0, 80000)) +
  
  pdf(file = "final_figures/ROC_yeast_slic.pdf", height = 4, width = 6)
print(p)
dev.off()

##------------ ROC curves on consensus cluster level ------------## 
# aggregate the clusters across the samples 
aggregateTagClusters(LICAGE, 
                     tpmThreshold = 5, 
                     qLow = 0.1, 
                     qUp = 0.9, 
                     maxDist = 100)

samples <- sampleLabels(LICAGE)[-length(sampleLabels(LICAGE))]

CTSSnorm <- consensusClustersTpm(LICAGE)

# convert to dataframe and add columns TP(true positive) FP(false positive)
CTSSnorm <- as.data.frame(CTSSnorm)

# extract TPM values as pair sample vs nanti
samples <- sampleLabels(LICAGE)[-length(sampleLabels(LICAGE))]
tpm_pairs.dfl <- lapply(samples, function(x) {
  df <- CTSSnorm[, c(x, "BY4741")]
  df <- df[rowSums(df) > 0, ]
  return(df) })
names(tpm_pairs.dfl) <- samples

# calculate TP and FP per different TPM threshold
TP_FP_l <- list()
TP_FP.gg <- list()

for(j in 1:length(tpm_pairs.dfl)){
  threshold <- seq(0, 500, 0.5) 
  for (i in 1:length(threshold)) {
    tpm_filt <- tpm_pairs.dfl[[j]][tpm_pairs.dfl[[j]][, 1] >= threshold[i], ]
    TP_FP_l[[i]] <- data.frame("TP" = tpm_filt[,1] > 0 & tpm_filt[,2] > 0,
                               "FP" = tpm_filt[,1] > 0 & tpm_filt[,2] == 0)
    
  }
  TP_FP.gg[[j]] <- as.data.frame(do.call(rbind, lapply(TP_FP_l, colSums)))
  TP_FP.gg[[j]]$sample <- rep(names(tpm_pairs.dfl[j]), times = nrow(TP_FP.gg[[j]]))
}

# flatten for plotting
TP_FP.gg <- as.data.frame(do.call(rbind, TP_FP.gg))
TP_FP.gg$sample <- factor(TP_FP.gg$sample, levels = samples)

library(viridis)
col <- magma(10)

p <- ggplot(TP_FP.gg, aes(x = FP, y = TP, fill = sample)) +
  geom_line(aes(col = sample), size = 1,alpha = 0.7) +
  geom_point(shape = 21, size = 1, alpha = 0.7) +
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) +
  theme(text = element_text(size = 12,
                            family = "Helvetica"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12, family = "Helvetica", color = "black"),
        legend.title = element_text(size = 12, family = "Helvetica", colour = "black"),
        legend.text = element_text(size = 12, family = "Helvetica", colour = "black")) +
  xlab("False positive TCs") +
  ylab("True positive TCs") +
  xlim(c(0, 1600)) +
  ylim(c(0, 5000))

pdf(file = "final_figures/ROC_yeast_slic_TCs.pdf", height = 4, width = 6)
print(p)
dev.off()
