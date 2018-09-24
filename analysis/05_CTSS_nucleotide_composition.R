#/*==========================================================================#*/
#' ## Figures 2C and G - nucleotide composition of identified CTSSs
#' 
#/*==========================================================================#*/

library(tidyr)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# change chromosome names to NCBI format
seqnames(BSgenome.Scerevisiae.UCSC.sacCer3) <- c("NC_001133.9","NC_001134.8","NC_001135.5",
                                                 "NC_001136.10", "NC_001137.3", "NC_001138.5", 
                                                 "NC_001139.9", "NC_001140.6", "NC_001141.2", 
                                                 "NC_001142.9", "NC_001143.9", "NC_001144.5", 
                                                 "NC_001145.3", "NC_001146.8", "NC_001147.6", 
                                                 "NC_001148.4", "NC_001224.1")
seqnames(BSgenome.Scerevisiae.UCSC.sacCer3)

#extract CTSS dataframe with normalized TPM per CTSS per sample
LICAGE_ctss <- CTSSnormalizedTpm(LICAGE)

#filter to have at least 1TPM 
filtr_idx.list <- list()
samples_all <- colnames(LICAGE_ctss)[-c(1:3)]

for (i in 1:length(samples_all)) {
  filtr_idx.list[[i]] <- c(LICAGE_ctss[, samples_all[i]] >= 1)
}
names(filtr_idx.list) <- samples_all

#create GRanges object CTSSs (filter per sample to have >= 1 TPM)
filtr.ctss.grl <- lapply(filtr_idx.list, function(x) GRanges(seqnames = LICAGE_ctss[x, ]$chr, 
                                                             ranges = IRanges(start = LICAGE_ctss[x, ]$pos, width = 1),
                                                             strand = LICAGE_ctss[x, ]$strand))
#get CTSS sequence
filtr.ctss.seq <- lapply(filtr.ctss.grl, function(x) getSeq(BSgenome.Scerevisiae.UCSC.sacCer3, x))


#change names
samples <- c("SLIC 1 ng", "SLIC 2 ng", "SLIC 5 ng", "SLIC 10 ng r1", 
             "SLIC 10 ng r2", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi PCR", "nAnTi 5 ug" )
names(filtr.ctss.seq) <- samples

#calculate nucleotide frequency detected CTSS (filtr to have >=1TPM per sample)
CTSS_nucl_freq.l <- lapply(filtr.ctss.seq, function(x) nucleotideFrequencyAt(x, at = 1) / length(x) * 100)
CTSS_nucl_freq <- data.frame(matrix(unlist(CTSS_nucl_freq.l), nrow = length(CTSS_nucl_freq.l), byrow = T))
colnames(CTSS_nucl_freq) <- names(CTSS_nucl_freq.l[[1]])
rownames(CTSS_nucl_freq) <- samples
CTSS_nucl_freq$samples <- rownames(CTSS_nucl_freq)

CTSS_nucl_freq.tidy <- gather(CTSS_nucl_freq, Feature, Percentage, (c("A", "C", "G", "T")))

# plot start nucleotide frequencies as a histogram for selected samples
# set levels
CTSS_nucl_freq.tidy$samples <- factor(CTSS_nucl_freq.tidy$samples, 
                                      levels = samples[length(samples):1])

# set levels for nucleotides
CTSS_nucl_freq.tidy$Feature <- factor(CTSS_nucl_freq.tidy$Feature, levels = c("T", "G", "C", "A"))

# select samples for ploting
selected <-  c("SLIC 5 ng", "SLIC 10 ng r1", "SLIC 10 ng r2", "nAnTi 5 ug" )
CTSS_nucl_sel <- CTSS_nucl_freq.tidy[CTSS_nucl_freq.tidy$samples %in% selected, ]
CTSS_nucl_sel$samples <- factor(CTSS_nucl_sel$samples, levels = selected[length(selected):1])

# select colours
library(viridis)
col <- magma(5, alpha = 0.7)[c(1,4,3,2)]

# select colours - standard nucleotide color scheme - A, C, G, T (green, blue, yellow, red)
col <- c("darkolivegreen3", "cornflowerblue", "gold", "firebrick3")[4:1]

library(ggplot2)
p <- ggplot(CTSS_nucl_sel, aes(x = samples, y = Percentage, fill = Feature)) +
  geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
  coord_flip() +
  scale_fill_manual("Nucleotide", values = col) +
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Percentage", x = NULL)

pdf(file = "final_figures/slicCAGE_CTSS_freq_sel_sc.pdf", height = 2, width = 5)
print(p)
dev.off()

# plot start nucleotide frequencies as a histogram for all samples
# select colours
library(viridis)
col <- magma(5, alpha = 0.7)[c(1,4,3,2)]

# select colours - standard nucleotide color scheme - A, C, G, T (green, blue, yellow, red)
col <- c("darkolivegreen3", "cornflowerblue", "gold", "firebrick3")[4:1]

library(ggplot2)
p <- ggplot(CTSS_nucl_freq.tidy, aes(x = samples, y = Percentage, fill = Feature)) +
  geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
  coord_flip() +
  scale_fill_manual("Nucleotide", values = col) +
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Percentage", x = NULL)

pdf(file = "final_figures/slicCAGE_CTSS_freq_all_sc.pdf", height = 3, width = 5)
print(p)
dev.off()