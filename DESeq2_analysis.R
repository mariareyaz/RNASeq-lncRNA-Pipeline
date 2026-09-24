# Differential Expression Analysis using DESeq2
library(DESeq2)
library(ggplot2)
library(clusterProfiler)

# Load count matrix and metadata
counts <- read.table("counts/matrix.txt", header=TRUE, row.names=1)
colData <- read.csv("metadata.csv", row.names=1)

# Construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)

# Extract results
res <- results(dds, pAdjustMethod="BH", alpha=0.05)
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Generate Volcano Plot
ggplot(as.data.frame(res), aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=padj < 0.05 & abs(log2FoldChange) > 1), alpha=0.6) +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal() +
  labs(title="Volcano Plot of Differential Expression", x="log2 Fold Change", y="-log10 Adjusted P-value")
