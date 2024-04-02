01_singlecell/00_scRNA_get_expression_profile.sh is used to get expression profile from scRNA-seq fastq file.
01_singlecell/01_seurat_filtering.R is used to create Seurat object and filter some contaminated cells.
01_singlecell/02_monocle_pseudotime.R is used to perform pseudotime analysis.
01_singlecell/03_function_analysis.R is used to perform enrichment analysis, Bcl11b activity analysis, pathway activity analysis and transcription factor enrichment analysis along pseudotime
01_singlecell/04_cellCycle_analysis.R is used to perform cell cycle analysis.
01_singlecell/05_TCGA_data_usage.R is used to perform analysis on TCGA data.
01_singlecell/06_scEntropy.m is used to perform entropy analysis with MATLAB.
02_ChIP-seq/01_ChIPseq_getBAM.sh is used to get mapped bam file from ChIPseq fastq file.
02_ChIP-seq/02_ChIPseq_getPeak.sh is used to get binding peaks from mapped bam file.
02_ChIP-seq/03_ChIPseq_annotation.R is used to annotate peaks.
