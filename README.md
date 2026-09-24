# RNASeq-lncRNA-Pipeline
command-line bash scripts and R workflows for differential expression and stress-responsive lncRNA profiling.
# Stress-Responsive lncRNA Analysis Pipeline

An end-to-end transcriptomic pipeline for processing raw RNA-Seq data, identifying long non-coding RNAs (lncRNAs), and modeling gene-regulatory networks under abiotic stress.

## Workflow Overview
1. **Quality Control & Preprocessing:** `FastQC` and `Trimmomatic` for adapter trimming and quality filtering.
2. **Read Alignment:** `HISAT2` indexing and spliced alignment to reference genome.
3. **Transcript Assembly & Feature Quantification:** `StringTie` assembly and `featureCounts` quantification.
4. **Differential Expression:** `DESeq2` analysis in R to identify significant log2 fold-changes.
5. **lncRNA Filtering & Target Prediction:** Coding potential calculation (`GffCompare`/`CPC2`), cis/trans target prediction, and pathway enrichment (`clusterProfiler`).
6. **Network Visualization:** Importing lncRNA-mRNA interactions into `Cytoscape`.

## Repository Contents
- `pipeline.sh`: Automated bash script for quality control, alignment, and transcript quantification.
- `DESeq2_analysis.R`: R script for differential gene expression analysis, volcano plots, and functional enrichment.
