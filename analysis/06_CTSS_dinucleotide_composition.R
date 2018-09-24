#/*==========================================================================#*/
#' ## Figures 2D and H - dinucleotide composition of identified CTSSs
#' 
#/*==========================================================================#*/

#-------  ALL CTSSs -------# 
# extract CTSS dataframe with normalized TPM per CTSS per sample
LICAGE_ctss <- CTSSnormalizedTpm(LICAGE)

# filter to have at least 1TPM 
filtr_idx.list <- list()
samples_all <- colnames(LICAGE_ctss)[-c(1:3)]

for (i in 1:length(samples_all)) {
  filtr_idx.list[[i]] <- c(LICAGE_ctss[, samples_all[i]] >= 1)
}
names(filtr_idx.list) <- samples_all

# create GRanges object CTSSs (filter per sample to have >= 1 TPM)
filtr.ctss.grl <- lapply(filtr_idx.list, function(x) GRanges(seqnames = LICAGE_ctss[x, ]$chr, 
                                                             ranges = IRanges(start = LICAGE_ctss[x, ]$pos, width = 1),
                                                             strand = LICAGE_ctss[x, ]$strand))

# add genome info
for (i in 1:length(filtr.ctss.grl)) {
  seqinfo(filtr.ctss.grl[[i]]) <- seqinfo(BSgenome.Scerevisiae.UCSC.sacCer3)[seqlevels(filtr.ctss.grl[[i]])]
}

# expand to include upstream nucleotide (to get (-1, +1))
filtr.ctss_up.grl <- lapply(filtr.ctss.grl, function(x) promoters(x, upstream = 1, downstream = 1))

# get initiation dinucleotide sequence
filtr.ctss_up.seq <- lapply(filtr.ctss_up.grl, function(x) getSeq(BSgenome.Scerevisiae.UCSC.sacCer3, x))

# count dinucleotides in initiator CTSS (filtered to have >=1TPM per sample)
CTSS_nucl_count.l <- lapply(filtr.ctss_up.seq, table)

# convert dinucleotide count to percentages
CTSS_nucl_freq.l <- lapply(CTSS_nucl_count.l, function(x) x/sum(x) * 100)

# set levels of dinucleotide usage for plotting
dinucl.levels <- names(CTSS_nucl_freq.l[[10]])[order(CTSS_nucl_freq.l[[10]])]

# convert dataframe to a list
CTSS_nucl_freq.df <- as.data.frame(CTSS_nucl_freq.l)

# tidy the dataframe
CTSS_nucl_freq_tidy.df <- CTSS_nucl_freq.df[, grep(x = colnames(CTSS_nucl_freq.df), pattern = "Freq")]

# attach dinucleotide information
CTSS_nucl_freq_tidy.df <- cbind(CTSS_nucl_freq_tidy.df, CTSS_nucl_freq.df$Y1ng_carrier.Var1)

# rename columns
names <- c("SLIC 1 ng", "SLIC 2 ng", "SLIC 5 ng", "SLIC 10 ng r1", 
           "SLIC 10 ng r2", "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi PCR", "nAnTi 5 ug" )
colnames(CTSS_nucl_freq_tidy.df) <- c(names, "dinucleotide")
rownames(CTSS_nucl_freq_tidy.df) <- CTSS_nucl_freq_tidy.df$dinucleotide

# prepare dataframe for ggplot
CTSS_nucl_freq_tidy.gg <- gather(CTSS_nucl_freq_tidy.df, 
                                 ...= colnames(CTSS_nucl_freq_tidy.df[,1:10]),
                                 key = "samples",
                                 value = "percentage")

# plot initiators as a histogram - selected samples
# select samples
selected <- c("SLIC 5 ng", "SLIC 10 ng r1", "SLIC 10 ng r2", "nAnTi 5 ug")
initiator_sel <- CTSS_nucl_freq_tidy.gg[CTSS_nucl_freq_tidy.gg$samples %in% selected, ]

initiator_sel$samples <- factor(initiator_sel$samples, levels = selected[length(selected):1])
initiator_sel$dinucleotide <- factor(initiator_sel$dinucleotide, levels = dinucl.levels)

library(ggplot2)
library(viridis)
col = magma(4, alpha = 0.8)[4:1]

p <- ggplot(data = initiator_sel, aes(x = dinucleotide,
                                      y = percentage,
                                      fill = samples)) + 
  scale_fill_manual(values = col) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
  coord_flip() +
  xlab("Initiator dinucleotide") + ylab("Percentage") + 
  ggtitle(NULL) + 
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 30)) +
  labs(fill = NULL)

pdf(file = "final_figures/initiator_sel_sc.pdf", height = 5, width = 4)
print(p)
dev.off()

# plot initiators as a histogram - all samples
# set levels
CTSS_nucl_freq_tidy.gg$dinucleotide <- factor(CTSS_nucl_freq_tidy.gg$dinucleotide, levels = dinucl.levels)
CTSS_nucl_freq_tidy.gg$samples <- factor(CTSS_nucl_freq_tidy.gg$samples, levels = names[length(names):1])

library(ggplot2)
library(viridis)
col = magma(10, alpha = 0.8)[10:1]

p <- ggplot(data = CTSS_nucl_freq_tidy.gg, aes(x = dinucleotide,
                                               y = percentage,
                                               fill = samples)) + 
  scale_fill_manual(values = col) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
  coord_flip() +
  xlab("Initiator dinucleotide") + ylab("Percentage") + 
  ggtitle(NULL) + 
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 30)) +
  labs(fill = NULL)

pdf(file = "final_figures/initiator_all_sc.pdf", height = 5, width = 4)
print(p)
dev.off()

#-------  Only dominant CTSSs -------# 
# load domTSS centered GRanges object
domTSS.grl <- readRDS(file = "data/domTSS_grl.RDS")

# get initiatior region -> (-1, +1)
domTSS.grl <- lapply(domTSS.grl, function(x) promoters(x, upstream = 1, downstream  = 1))

# filter to have => 1 TPM
for (i in 1:length(domTSS.grl)) {
  domTSS.grl[[i]] <- domTSS.grl[[i]][domTSS.grl[[i]]$tpm.dominant_ctss >= 1, ]
}

# extract dominant TSS initiator sequence
domTSS_seq.l <- lapply(domTSS.grl, function(x) getSeq(BSgenome.Scerevisiae.UCSC.sacCer3, x))

# count dominant TSS initiator sequences
domTSS_init_count.l <- lapply(domTSS_seq.l, table)

# convert dominant TSS initiator counts to percentages
domTSS_init_count.l <- lapply(domTSS_init_count.l, function(x) x/sum(x) * 100)

# order of dinucleotides according to control
dinucl.levels <- names(domTSS_init_count.l[[10]])[order(domTSS_init_count.l[[10]], decreasing = FALSE)]

# convert list to a dataframe
domTSS_init_count.df <- as.data.frame(domTSS_init_count.l)

# tidy the dataframe
domTSS_init_freq_tidy.df <- domTSS_init_count.df[, grep(x = colnames(domTSS_init_count.df), pattern = "Freq")]

# attach dinucleotide information
domTSS_init_freq_tidy.df <- cbind(domTSS_init_freq_tidy.df, domTSS_init_count.df$Y1ng_carrier.Var1)

# rename columns
names <- c("SLIC 1 ng", "SLIC 2 ng", "SLIC 5 ng", "SLIC 10 ng r1", "SLIC 10 ng r2", 
           "SLIC 25 ng", "SLIC 50 ng", "SLIC 100 ng", "nAnTi PCR", "nAnTi 5 ug" )
colnames(domTSS_init_freq_tidy.df) <- c(names, "dinucleotide")
rownames(domTSS_init_freq_tidy.df) <- domTSS_init_freq_tidy.df$dinucleotide

# prepare dataframe for ggplot
ddomTSS_init_freq_tidy.gg <- gather(domTSS_init_freq_tidy.df, 
                                    ...= colnames(domTSS_init_freq_tidy.df[,1:10]),
                                    key = "samples",
                                    value = "percentage")

# plot dominant start nucleotide frequencies as a histogram - all samples
# set levels
ddomTSS_init_freq_tidy.gg$samples <- factor(ddomTSS_init_freq_tidy.gg$samples, levels = names[length(names):1])
ddomTSS_init_freq_tidy.gg$dinucleotide <- factor(ddomTSS_init_freq_tidy.gg$dinucleotide,
                                                 levels = dinucl.levels)
library(ggplot2)
library(viridis)

col = magma(10, alpha = 0.8)[10:1]

p <- ggplot(data = ddomTSS_init_freq_tidy.gg, aes(x = dinucleotide,
                                                  y = percentage,
                                                  fill = samples)) + 
  scale_fill_manual(values = col) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
  coord_flip() +
  xlab("Initiator dinucleotide") + ylab("Percentage") + 
  ggtitle(NULL) + 
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(fill = NULL)

pdf(file = "final_figures/domTSS_initiator_all_sc.pdf", height = 5, width = 4)
print(p)
dev.off()

# plot dominant start nucleotide frequencies as a histogram - selected samples
# set levels
selected <- c("SLIC 5 ng", "SLIC 10 ng r1", "SLIC 10 ng r2", "nAnTi 5 ug")
domTSS_select <- ddomTSS_init_freq_tidy.gg[ddomTSS_init_freq_tidy.gg$samples %in% selected, ]
domTSS_select$samples <- factor(domTSS_select$samples, levels = selected[length(selected):1])
library(ggplot2)
library(viridis)

col = magma(4, alpha = 0.8)[4:1]

p <- ggplot(data = domTSS_select, aes(x = dinucleotide,
                                      y = percentage,
                                      fill = samples)) + 
  scale_fill_manual(values = col) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
  coord_flip() +
  xlab("Initiator dinucleotide") + ylab("Percentage") + 
  ggtitle(NULL) + 
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(fill = NULL)

pdf(file = "final_figures/domTSS_initiator_sel_sc.pdf", height = 5, width = 4)
print(p)
dev.off()