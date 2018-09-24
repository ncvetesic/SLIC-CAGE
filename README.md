# SLIC-CAGE: high-resolution transcription start site mapping using nanogram-levels of total RNA

## Abstract
Cap analysis of gene expression (CAGE) is a methodology for genome-wide quantitative mapping of mRNA 5’ends to precisely capture transcription start sites at a single nucleotide resolution. In combination with high-throughput sequencing, CAGE has revolutionized our understanding of rules of transcription initiation, led to discovery of new core promoter sequence features and discovered transcription initiation at enhancers genome-wide. The biggest limitation of CAGE is that even the most recently improved version (nAnT-iCAGE) still requires large amounts of total cellular RNA (5 micrograms), preventing its application to scarce biological samples such as those from early embryonic development or rare cell types. Here, we present SLIC-CAGE, a Super-Low Input Carrier-CAGE approach to capture 5’ends of RNA polymerase II transcripts from as little as 5-10 ng of total RNA. The dramatic increase in sensitivity is achieved by specially designed, selectively degradable carrier RNA. We demonstrate the ability of SLIC-CAGE to generate data for genome-wide promoterome with 1000-fold less material than required by existing CAGE methods by generating a complex, high quality library from mouse embryonic day (E) 11.5 primordial germ cells.

## Requirements 

* Install following R packages from CRAN: 
dplyr, tidyr, ggplot2, RColorBrewer, viridis, reshape2, cowplot, corrplot, zoo

* Install following R packages from bioconductor:
CAGEr, SeqPattern, heatmaps, GenomicRanges, BSgenome.Scerevisiae.UCSC.sacCer3, BSgenome.Mmusculus.UCSC.mm10, TxDb.Mmusculus.UCSC.mm10.knownGene, GenomicFeatures, org.Mm.eg.db, ChIPseeker, TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, biomaRt, Gviz, DESeq2

## Figure to code map 

* Figure 1 [B-D and F-H](analysis/01_pairwise_ctss_corr.R), [E and I](analysis/02_gviz_gbrowser_views.R)  
* Figure 2 [A](analysis/03_genomic_location_tc.R)


## Additional information 

* Pre-processing pipelines are available upon request.

* Processed data (*.bam, *.bw, CTSS or tag cluster tables, etc.) and intermediate are available upon request.

## Reference 
**SLIC-CAGE: high-resolution transcription start site mapping using nanogram-levels of total RNA**

[biorxiv preprint](https://www.biorxiv.org/content/early/2018/07/19/368795)