#/*==========================================================================#*/
#' ## Supplemental Figure S12- dinucl. comp. of domTSS split by TPM value
#'
#/*==========================================================================#*/

# Stratify domTSS initiators per TPM value 
# load domTSS centered GRanges object
domTSS.grl <- readRDS(file = "../intermediate/domTSS_nanoCAGE_yeast_grl")

# get initiatior region -> (-1, +1)
domTSS_init.grl <- lapply(domTSS.grl, function(x) promoters(x, upstream = 1, downstream  = 1))

# filter domTSS with value less than 1 TPM
for (i in 1:length(domTSS_init.grl)) {
  domTSS_init.grl[[i]] <- domTSS_init.grl[[i]][domTSS_init.grl[[i]]$tpm.dominant_ctss >= 1, ]
}

# extract dominant TSS initiator sequence
domTSS_init_seq.grl <- lapply(domTSS_init.grl, function(x) getSeq(BSgenome.Scerevisiae.UCSC.sacCer3, x))

# add bin/quartile information to tpm.dominant_ctss
library(dplyr)

tpm.l <- list()
for(i in 1:length(domTSS_init.grl)) {
  tpm.l[[i]] <- data.frame(tpm.dominant_ctss = mcols(domTSS_init.grl[[i]])$tpm.dominant_ctss)
  tpm.l[[i]] <- tpm.l[[i]] %>% mutate(quartile_tpm.dominant_ctss = ntile(tpm.dominant_ctss, 4))
}
names(tpm.l) <- names(domTSS_init.grl)

# create idx for quartile GRanges retrieval
idx_q1.l <- lapply(tpm.l, function(x) x$quartile_tpm.dominant_ctss == 1)
idx_q2.l <- lapply(tpm.l, function(x) x$quartile_tpm.dominant_ctss == 2)
idx_q3.l <- lapply(tpm.l, function(x) x$quartile_tpm.dominant_ctss == 3)
idx_q4.l <- lapply(tpm.l, function(x) x$quartile_tpm.dominant_ctss == 4)

# combine all indices to a list
idx_all <- list(idx_q1.l, idx_q2.l, idx_q3.l, idx_q4.l)
names(idx_all) <- c("idx_q1", "idx_q2", "idx_q3", "idx_q4")

# extract domTSS initiators for each domTPM quartile
quart.l <- list()
quart_nest.l <- list()

for(i in 1:length(idx_all)) {
  for (j in 1:length(idx_all[[i]])) {
    quart_nest.l[j] <- domTSS_init_seq.grl[[j]][idx_all[[i]][[j]]]
  }
  names(quart_nest.l) <- names(idx_all[[i]])
  quart.l[[i]] <- quart_nest.l 
}
names(quart.l) <- c("q1", "q2", "q3", "q4") # add quartile info

# function to count dominant TSS initiatior, convert to percentages, format to a dataframe
# input is a DNAStringSet list of dinucleotides
TSS_initiators_count <- function(seq.l) {
  
  # count dominant TSS initiator sequences
  seq.count.l <- lapply(seq.l, table)
  
  # convert counts to percentages
  seq.perc.l <- lapply(seq.count.l, function(x) x/sum(x)*100)
  
  # convert list to a dataframe - need to add 0 to missing dinucleotides
  # create a dinucleotide vector - all 16
  dinucleotides <- expand.grid(x = as.vector(c("A", "C", "G", "T")), 
                               y = as.vector(c("A", "C", "G", "T")),
                               stringsAsFactors = FALSE) 
  
  # convert dataframe to a vector of dinucleotides
  dinucl.vector <- vector()
  for (i in 1:nrow(dinucleotides)) {
    dinucl.vector[i] <- paste(dinucleotides[i, ], collapse = "")
  }
  
  # add 0 value to missing dinucleotides
  missing <- list()
  names.l <- list()
  
  for (i in 1:length(seq.perc.l)) {
    missing[[i]] <- setdiff(dinucl.vector, names(seq.perc.l[[i]]))
    seq.perc.l[[i]][missing[[i]]] <- 0 # add zero value
    seq.perc.l[[i]] <- seq.perc.l[[i]][order(names(seq.perc.l[[i]]), decreasing = FALSE)] #order alphabetically
    names.l[[i]] <- names(seq.perc.l[[i]]) # store names
    seq.perc.l[[i]] <- as.numeric(seq.perc.l[[i]]) # convert elements to numeric (from table)
    names(seq.perc.l[[i]]) <- names.l[[i]] # restore names
  }
  
  # convert initiator list to a dataframe
  seq.perc.df <- as.data.frame(seq.perc.l)
  seq.perc.df$dinucleotide <- rownames(seq.perc.df) # attach dinucleotide column
  return(seq.perc.df)
}


# calculate initiation percentages in each TPM quartile
quartile_init.l <- lapply(quart.l, TSS_initiators_count)

library(tidyr)
# 4-highest TPM, 1-lowest TPM
# prepare dataframes for ggplot
quartile_init_gg.l <- lapply(quartile_init.l, 
                             function(x) {
                               gather(x, ...=colnames(x[, 1:10]),
                                      key = "samples",
                                      value = "percentage")
                             })


samples <- sampleLabels(myCAGEset)

# set levels
for (i in 1:length(quartile_init_gg.l)) {
  quartile_init_gg.l[[i]]$samples <- factor(quartile_init_gg.l[[i]]$samples,
                                            levels = samples)
}


# plot initiator frequency for each quantile
library(ggplot2)
library(viridis)
library(gridExtra)

p.l <- list()
p.l <- lapply(quartile_init_gg.l, 
              function(var) {
                ggplot(data = var, aes(x = dinucleotide,
                                       y = percentage,
                                       fill = samples)) + 
                  geom_bar(stat = "identity", position = position_dodge(), colour = "black") + 
                  scale_fill_manual(values = viridis(n = 10, alpha = 0.8, begin = 0, option = "viridis")) +
                  xlab("Dominant TSS initiator dinucleotide") + ylab("Percentage") + 
                  ggtitle("Initiator frequency") + 
                  theme_bw() +
                  theme(text = element_text(size = 16), 
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14)) +
                  scale_y_continuous(limits = c(0, 40)) +
                  labs(fill = "Samples") 
              })
# add plot annotation
for(i in 1:length(p.l)){
  p.l[[i]] <- p.l[[i]] + annotate("text", x = 13, y = 40, label = paste("dominant TPM quantile", i, sep = ""), size = 6)
}
pdf(file = "../results/nanoCAGE_domTSS_stratTPM_freq.pdf", width = 18, height = 12)
grid.arrange(grobs = p.l, nrow = 2, ncol = 2)
dev.off()
