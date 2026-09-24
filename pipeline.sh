#!/bin/bash
# RNA-Seq & lncRNA Processing Pipeline

# 1. Quality Control
fastqc -o qc_reports/ raw_data/*.fastq.gz

# 2. Adapter & Quality Trimming
trimmomatic PE -phred33 \
  raw_data/sample_R1.fastq.gz raw_data/sample_R2.fastq.gz \
  trimmed/sample_R1_paired.fq.gz trimmed/sample_R1_unpaired.fq.gz \
  trimmed/sample_R2_paired.fq.gz trimmed/sample_R2_unpaired.fq.gz \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# 3. Read Alignment (HISAT2)
hisat2 -x genome_index -1 trimmed/sample_R1_paired.fq.gz -2 trimmed/sample_R2_paired.fq.gz | \
  samtools view -bS - | samtools sort -o sorted_bam/sample.sorted.bam

# 4. Transcript Assembly & Quantification
stringtie sorted_bam/sample.sorted.bam -G reference.gff -o assembly/sample.gtf
featureCounts -p -t exon -g gene_id -a reference.gff -o counts/matrix.txt sorted_bam/*.bam
