# Install necessary libraries if not already installed
install.packages(c("ggplot2", "reshape2", "pheatmap"))
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot", "org.Mm.eg.db"))

# Load libraries
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(reshape2)

# Load count data
count_data <- read.csv("em.csv", header = TRUE)

# Load annotation data
annotation_data <- read.csv("annotation.csv", header = TRUE)

# Merge count data with annotation data
merged_data <- merge(count_data, annotation_data, by = "ID")

# Load column data (metadata)
col_data <- read.csv("col_data.csv", header = TRUE)

# Prepare DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = merged_data[, 2:ncol(merged_data)],
                              colData = col_data,
                              design = ~ condition)




# Run DESeq analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Variance stabilization for downstream analysis
vsd <- vst(dds)

ma_data <- data.frame(log2FoldChange = res$log2FoldChange,
                      meanExpression = log10(res$baseMean + 1),
                      significant = ifelse(res$padj < 0.05, "Yes", "No"))

ggplot(ma_data, aes(x = meanExpression, y = log2FoldChange, color = significant)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
  labs(title = "MA Plot", x = "Log10 Mean Expression", y = "Log2 Fold Change") +
  theme_minimal()

# Volcano plot
res_df <- as.data.frame(res)
res_df$neg_log10_pvalue <- -log10(res_df$pvalue)
res_df$category <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Up",
                          ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Down", "No change"))

ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_pvalue, color = category)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "No change" = "black")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal()


# Violin plot
melted_data <- melt(count_data[, -1], variable.name = "sample", value.name = "expression")
melted_data$condition <- col_data$condition[match(melted_data$sample, col_data$sampleName)]

ggplot(melted_data, aes(x = sample, y = expression, fill = condition)) +
  geom_violin(trim = TRUE) +
  labs(title = "Violin Plot of Gene Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))


# Heat map
normalized_data <- log2(count_data[, -1] + 1)

pheatmap(normalized_data,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         scale = "row",
         main = "Heatmap of Gene Expression")


# GO Enrichment Analysis
# Load gene list (replace with your actual gene list)
gene_list <- c("Gene1", "Gene2", "Gene3") # Replace with your significant gene list

# Perform GO enrichment analysis (Biological Process category)
ego <- enrichGO(gene          = gene_list,
                OrgDb         = org.Mm.eg.db, # Replace with appropriate organism database
                keyType       = "SYMBOL",    # Use "ENSEMBL" if your gene IDs are Ensembl IDs
                ont           = "BP",       # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

# View top results
head(ego)


# Create dot plot
dot_plot <- dotplot(ego, showCategory = 10) +
  labs(title = "GO Enrichment Dot Plot")

# Save the dot plot to your PC
ggsave("go_dot_plot.png", plot = dot_plot, width = 8, height = 6)

# Create enrichment map
emap_plot <- emapplot(ego, showCategory = 10, layout = "kk") +
  labs(title = "GO Enrichment Map")

# Save the enrichment map to your PC
ggsave("go_enrichment_map.png", plot = emap_plot, width = 10, height = 8)


# Create category network plot (cnetplot)
cnet_plot <- cnetplot(ego,
                      showCategory = 10,
                      foldChange = NULL) + # Add fold change data if available
  labs(title = "GO Category Network Plot")

# Save the network plot to your PC
ggsave("go_network_plot.png", plot = cnet_plot, width = 10, height = 8)


# Create bar plot
bar_plot <- barplot(ego, showCategory = 10, color = "p.adjust") +
  labs(title = "GO Enrichment Bar Plot")

# Save the bar plot to your PC
ggsave("go_bar_plot.png", plot = bar_plot, width = 8, height = 6)



