# SLIC-CAGE: high-resolution transcription start site mapping using nanogram-levels of total RNA

## Abstract
Cap analysis of gene expression (CAGE) is a methodology for genome-wide quantitative mapping of mRNA 5’ends to precisely capture transcription start sites at a single nucleotide resolution. In combination with high-throughput sequencing, CAGE has revolutionized our understanding of rules of transcription initiation, led to discovery of new core promoter sequence features and discovered transcription initiation at enhancers genome-wide. The biggest limitation of CAGE is that even the most recently improved version (nAnT-iCAGE) still requires large amounts of total cellular RNA (5 micrograms), preventing its application to scarce biological samples such as those from early embryonic development or rare cell types. Here, we present SLIC-CAGE, a Super-Low Input Carrier-CAGE approach to capture 5’ends of RNA polymerase II transcripts from as little as 5-10 ng of total RNA. The dramatic increase in sensitivity is achieved by specially designed, selectively degradable carrier RNA. We demonstrate the ability of SLIC-CAGE to generate data for genome-wide promoterome with 1000-fold less material than required by existing CAGE methods by generating a complex, high quality library from mouse embryonic day (E) 11.5 primordial germ cells.

## Requirements 

* Install following R packages from CRAN: 
dplyr, tidyr, ggplot2, RColorBrewer, viridis

* Install following R packages from bioconductor:
CAGEr, Seqpattern, Heatmaps, GenomicRanges, BSgenome.Scerevisiae.UCSC.sacCer3, BSgenome.Mmusculus.UCSC.mm10, GenomicFeatures, org.Mm.eg.db, ChIPseeker, TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, biomaRt

