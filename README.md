# SLIC-CAGE: high-resolution transcription start site mapping using nanogram-levels of total RNA

## Abstract
Cap analysis of gene expression (CAGE) is a methodology for genome-wide quantitative mapping of mRNA 5’ends to precisely capture transcription start sites at a single nucleotide resolution. In combination with high-throughput sequencing, CAGE has revolutionized our understanding of rules of transcription initiation, led to discovery of new core promoter sequence features and discovered transcription initiation at enhancers genome-wide. The biggest limitation of CAGE is that even the most recently improved version (nAnT-iCAGE) still requires large amounts of total cellular RNA (5 micrograms), preventing its application to scarce biological samples such as those from early embryonic development or rare cell types. Here, we present SLIC-CAGE, a Super-Low Input Carrier-CAGE approach to capture 5’ends of RNA polymerase II transcripts from as little as 5-10 ng of total RNA. The dramatic increase in sensitivity is achieved by specially designed, selectively degradable carrier RNA. We demonstrate the ability of SLIC-CAGE to generate data for genome-wide promoterome with 1000-fold less material than required by existing CAGE methods by generating a complex, high quality library from mouse embryonic day (E) 11.5 primordial germ cells.

## Requirements 

* Install following R packages from CRAN: 
dplyr, tidyr, ggplot2, RColorBrewer, viridis, reshape2, cowplot, corrplot, zoo, philentropy

* Install following R packages from bioconductor:
CAGEr, SeqPattern, heatmaps, GenomicRanges, BSgenome.Scerevisiae.UCSC.sacCer3, BSgenome.Mmusculus.UCSC.mm10, TxDb.Mmusculus.UCSC.mm10.knownGene, GenomicFeatures, org.Mm.eg.db, ChIPseeker, TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, biomaRt, Gviz, DESeq2

## Figure to code map 

* Figure 1 [B-D and F-H](analysis/01_pairwise_ctss_corr.R), [E and I](analysis/02_gviz_gbrowser_views.R)  
* Figure 2 [A](analysis/03_genomic_location_tc.R), [B and F](analysis/04_distribution_iq_width.R), [C and G](analysis/05_CTSS_nucleotide_composition.R), [D and H](analysis/06_CTSS_dinucleotide_composition.R)  
* Figure 3 [A-F](analysis/01_pairwise_ctss_corr.R), [G](analysis/02_gviz_gbrowser_views.R), [H](analysis/03_genomic_location_tc.R), [I](analysis/04_distribution_iq_width.R), [J](analysis/05_CTSS_nucleotide_composition.R), [K](analysis/06_CTSS_dinucleotide_composition.R)  
* Figure 4 [A-C](analysis/07_heatmaps_TA_TATA_GC.R), [D](analysis/08_WW_periodicity_metaplot.R), [E](09_tag_cluster_coverage.R), [F-G](analysis/10_H3K4me3_coverage.R), [H-I](analysis/11_CpG_island_coverage.R)  
* Figure 5 [A](analysis/12_CAGE_RNAseq_corr.R), [B](analysis/04_distribution_iq_width.R), [C](analysis/03_genomic_location_tc.R), [D](analysis/05_CTSS_nucleotide_composition.R), [E](analysis/06_CTSS_dinucleotide_composition.R), [F](analysis/07_heatmaps_TA_TATA_GC.R), [G](analysis/13_SOM_expression_profiling.R), [H](analysis/14_SOM_clusters_genomic_locations.R), [I](analysis/15_SOM_clusters_GO_analysis.R), [J-K](analysis/02_gviz_gbrowser_views.R)  
  
* Supplemental Figure S1 [C-I](analysis/01_pairwise_ctss_corr.R), [J](analysis/03_genomic_location_tc.R), [K](analysis/04_distribution_iq_width.R), [L](analysis/05_CTSS_nucleotide_composition.R), [M](analysis/06_CTSS_dinucleotide_composition.R)  
* Supplemental Figure S2 [A-C](6_correlation_matrices.R), [D-F](analysis/03_genomic_location_tc.R), [G-I](analysis/05_CTSS_nucleotide_composition.R), [J-O](analysis/06_CTSS_dinucleotide_composition.R)  
* [Supplemental Figure S3](analysis/17_low_complex_simulation.R)  
* Supplemental Figure S4 [A-C](analysis/04_distribution_iq_width.R)  
* Supplemental Figure S5 [A-C](analysis/18_ROC_curves.R)  
* Supplemental Figure [S6-S8](analysis/19_distance_domTSS_SLIC_nAnTi.R)  
* Supplemental Figure [S9-S10, S13](analysis/20_CTSSs_missing_TPM_ratios.R)  
* [Supplemental Figure S11](analysis/21_domTSS_init_genom_loc.R)  
* [Supplemental Figure S12](analysis/22_domTSS_init_TPM_strat.R)  
* Supplemental Figure S14 [A-B](analysis/23_nanoCAGE_XL_corr.R), [C](analysis/19_distance_domTSS_SLIC_nAnTi.R), [D](analysis/analysis/04_distribution_iq_width.R), [E](analysis/03_genomic_location_tc.R), [F](analysis/05_CTSS_nucleotide_composition.R), [G](06_CTSS_dinucleotide_composition.R)  
* [Supplemental Figure S15](24_heatmap_similarity.R)  
* Supplemental Figure S16 - see Figure 4  
* Supplemental Figure S17 [A](25_IQ_width_distrib_threshold.R), [B](analysis/08_WW_periodicity_metaplot.R), [C](analysis/07_heatmaps_TA_TATA_GC.R), [D](26_TATA_box_metaplot.R)  
* Supplemental Figure S18 [B](27_CTSS_cons_correlation.R), [C](analysis/03_genomic_location_tc.R), [D](analysis/04_distribution_iq_width.R), [E](analysis/05_CTSS_nucleotide_composition.R), [F](analysis/06_CTSS_dinucleotide_composition.R), [G](analysis/07_heatmaps_TA_TATA_GC.R)  
* Supplemental Figure S19 [A](27_CTSS_cons_correlation.R), [B](analysis/04_distribution_iq_width.R)  
* [Supplemental Figure S23](analysis/28_nanoCAGE_cons_corr.R)  
  
For general processing of CAGE libraries starting from bam files see the CAGE processing [script.](analysis/29_CAGE_processing.R)









## Additional information 

* Pre-processing pipelines are available upon request.

* Processed data (*.bam, *.bw, CTSS or tag cluster tables, etc.) and intermediate data is available upon request.

## Reference 
**SLIC-CAGE: high-resolution transcription start site mapping using nanogram-levels of total RNA**

[biorxiv preprint](https://www.biorxiv.org/content/early/2018/07/19/368795)