# Chickpea Stress-Responsive lncRNA Analysis & AI Predictive Pipeline

An end-to-end hybrid bioinformatic and machine learning pipeline for processing raw RNA-Seq data, identifying long non-coding RNAs (lncRNAs), and predicting stress-responsive candidate regulators using AI classification models.

## Pipeline Architecture
1. **Quality Control & Read Processing (`01_ssp_download_qc.sh`):** Quality assessment with `FastQC` and adapter trimming using `Trimmomatic`.
2. **Alignment & Transcript Assembly (`02_lncRNA_identification_pipeline.sh`):** Reference mapping via `HISAT2`, transcript quantification via `StringTie` & `featureCounts`.
3. **Differential Expression & Functional Profiling (`DESeq2_analysis.R`):** Identification of significant log2 fold-changes and pathway enrichment using `clusterProfiler`.
4. **AI/ML Non-Coding Classifier (`03_ai_lncRNA_classifier.py`):** Supervised Machine Learning model (Random Forest / Feature Importance) trained on transcript sequence features (k-mers, ORF length, GC content) to classify novel lncRNAs and rank top biomarker candidates.

## Tech & AI Stack
- **Bioinformatics:** `FastQC`, `HISAT2`, `StringTie`, `featureCounts`, `DESeq2`
- **Machine Learning & Data Science:** Python (`scikit-learn`, `pandas`, `numpy`), R (`ggplot2`, `clusterProfiler`)
