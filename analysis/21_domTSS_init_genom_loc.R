#/*==========================================================================#*/
#' ## Supplemental Figure S11- dinucl. comp. of domTSS split in genomic loc
#'
#/*==========================================================================#*/

# Initiation dinucleotides stratified per genomic feature
# annotate initiation ranges with genomic features
# load domTSS centered GRanges object
domTSS.grl <- readRDS(file = "../intermediate/domTSS_nanoCAGE_yeast_grl")

# get initiatior region -> (-1, +1)
domTSS_init.grl <- lapply(domTSS.grl, function(x) promoters(x, upstream = 1, downstream  = 1))

for (i in 1:length(domTSS_init.grl)) {
  mcols(domTSS_init.grl[[i]]) <- mcols(tc_anno.grl[[i]])
}

# filter domTSS with value less than 1 TPM
for (i in 1:length(domTSS_init.grl)) {
  domTSS_init.grl[[i]] <- domTSS_init.grl[[i]][domTSS_init.grl[[i]]$tpm.dominant_ctss >= 1, ]
}

# extract sequence 
domTSS_init_seq.grl <- lapply(domTSS_init.grl, function(x) getSeq(BSgenome.Scerevisiae.UCSC.sacCer3, x))

# convert to list of dataframes
domTSS_init_seq.dfl <- lapply(domTSS_init_seq.grl, as.data.frame)

# add information on genomic features
for (i in 1:length(domTSS_init_seq.dfl)) {
  domTSS_init_seq.dfl[[i]] <- cbind(domTSS_init_seq.dfl[[i]], feature =  domTSS_init.grl[[i]]$annotation)
}

# add abbreviated feature
for (i in 1:length(domTSS_init_seq.dfl)) {
  
  # duplicate features and convert to character
  domTSS_init_seq.dfl[[i]] <- cbind(domTSS_init_seq.dfl[[i]], abbr_feature =  domTSS_init.grl[[i]]$annotation)
  domTSS_init_seq.dfl[[i]]$abbr_feature <- as.character(domTSS_init_seq.dfl[[i]]$abbr_feature)
  
  # change exonic features to exon
  domTSS_init_seq.dfl[[i]][grep("Exon", domTSS_init_seq.dfl[[i]]$feature), ]$abbr_feature <- "Exon"
  
  # change intronic features to intron - added if clause because not all sets have introns
  if (sum(grepl("Intron", domTSS_init_seq.dfl[[i]]$feature)) > 0) {
    domTSS_init_seq.dfl[[i]][grep("Intron", domTSS_init_seq.dfl[[i]]$feature), ]$abbr_feature <- "Intron"
  }
  
  # convert abbreviated features to factor
  domTSS_init_seq.dfl[[i]]$abbr_feature <- as.factor(domTSS_init_seq.dfl[[i]]$abbr_feature)
}

# summarize dinucleotide count for each sample per genomic feature
library(dplyr)
dinucl_count.l <- list()
inits_feats <- list()
for (i in 1:length(domTSS_init_seq.dfl)) {
  dinucl_count.l[[i]] <- domTSS_init_seq.dfl[[i]] %>% group_by(abbr_feature, x) %>% tally()
  # calculate the total amount of initiators per genomic feature
  inits_feats[[i]]<- dinucl_count.l[[i]] %>% group_by(abbr_feature) %>% summarize(sum = sum(n))
  # add the total_count column - total count per genomic feature
  dinucl_count.l[[i]] <- dinucl_count.l[[i]] %>% mutate(total_count = rep(inits_feats[[i]]$sum, table(abbr_feature)))
  # calculate percentage of each dinucleotide per genomic feature
  dinucl_count.l[[i]] <- dinucl_count.l[[i]] %>% mutate(percentage = round(n/total_count * 100, digits = 2))
}
# add names to the dinucleotide count
names(dinucl_count.l) <- names(domTSS_init_seq.dfl)

# add sample name column
for (i in 1:length(dinucl_count.l)) {
  dinucl_count.l[[i]] <- dinucl_count.l[[i]] %>% mutate(sample = names(dinucl_count.l[i]))
}                                                        

# combine dataframes from a list to 1 dataframe
dinucl_count.df <- do.call("rbind", dinucl_count.l)

# change column names
colnames(dinucl_count.df) <- c("abbr_feature", "dinucleotide", "count", "total_count", "percentage", "sample")

# change abbr_feature, dinucleotide and sample columns to factors
dinucl_count.df$abbr_feature <- as.factor(dinucl_count.df$abbr_feature)
dinucl_count.df$dinucleotide <- as.factor(dinucl_count.df$dinucleotide)
dinucl_count.df$sample <- as.factor(dinucl_count.df$sample)

# plot start nucleotide frequencies as a histogram
# set sample levels
dinucl_count.df$sample <- factor(dinucl_count.df$sample,levels = sampleLabels(myCAGEset))

library(ggplot2)
library(viridis)

# plot with feature faceting
p <- ggplot(data = dinucl_count.df, 
            aes(x = dinucleotide,
                y = percentage,
                fill = sample)) + 
  geom_bar(stat = "identity", position = position_dodge(), colour = "black") + 
  scale_fill_manual(values = viridis(n = 10, alpha = 0.8, begin = 0, option = "viridis")) +
  xlab("Dominant TSS initiator dinucleotide") + 
  ylab("Percentage") + 
  ggtitle("Dominant TSS initiator frequency") + 
  theme_bw() +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(fill = "Samples")

pdf(file = "../results/domTSS_dinucl_feats_feat.pdf", width = 10, height = 20)
p + facet_grid(abbr_feature ~ .)
dev.off()

# plot with sample faceting
p <- ggplot(data = dinucl_count.df, 
            aes(x = dinucleotide,
                y = percentage,
                fill = abbr_feature)) + 
  geom_bar(stat = "identity", position = position_dodge(), colour = "black") + 
  scale_fill_manual(values = viridis(n = 7, alpha = 0.8, begin = 0, option = "viridis")) +
  xlab("Dominant TSS initiator dinucleotide") + 
  ylab("Percentage") + 
  ggtitle("Dominant TSS initiator frequency") + 
  theme_bw() +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(fill = "Samples")

png(filename = "../results/domTSS_dinucl_feats_sample.png", width = 1000, height = 2000)
p + facet_grid(sample ~ .)
dev.off()