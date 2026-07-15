library(rtracklayer)
library(tidyverse)
library(tximport)
library(DESeq2)
library(matrixStats)
library(pheatmap)


PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"
pdf(file.path(PATH, "Dataset_0", "1_RNA-Seq", "RNA-Seq Graphs D0.pdf"), height = 10, width = 15)

### 1. File and reference download ###

## Downloading the count files
#Setting file directory
setwd(file.path(PATH, "Dataset_0", "1_RNA-Seq", "rnaseq"))

files_df <- read.delim(file.path(PATH, "Dataset_0", "Ground Truth", "samples.tsv"))
files <- paste0(files_df$sample, ".sf")
names(files) <- paste0("S", 1:100)

## Importing counts
tx_matrix <- tximport(files,
                      type  = "salmon",
                      txOut = TRUE)     #txOut=TRUE: no transcript->gene collapse (data is gene-level)


## Creating the metadata table
metadata <- data.frame(condition = factor(files_df$condition, levels = c("ctrl", "case")),
                       pair = factor(files_df$pair))
rownames(metadata) <- paste0("S", seq_len(nrow(files_df)))



### 2. DE analisys: DESeq2 ###

dds <- DESeqDataSetFromTximport(tx_matrix,
                                colData = metadata,
                                design = ~ pair + condition)

#Filtering: NO filtering for the Positive control dataset


##DESeq2:
dds <- DESeq(dds)


## Export: 
saveRDS(dds, file = file.path(PATH, "Dataset_0", "1_RNA-Seq", "dds.rds"))



### 3. Explorative analysis ###

## For this experiment there are 6 samples in a paired set-up. For initial visualization only PCA will show clustering and outliers (if any). 
## Heatmap will confirm clustering and add information on expression. 
## Then boxplots of overall expression for each sample will show whether any samples have a general over/under expression, and will verify that normalization was successful

vst_dds <- vst(dds, blind = FALSE)                      #Normalization: variance stabilizing transformation
counts_norm <- assay(vst_dds)                           #Extracting counts
counts_norm <- counts_norm[rowVars(counts_norm) > 0, ]  #Deleting rows with variance = 0
counts_norm_t <- t(counts_norm)                         #The matrix is needed transposed for pca

## PCA:
pca <- prcomp(counts_norm_t)
summary(pca)

#PCA graph:
#Creating the plot dataframe
pca_df <- as.data.frame(pca$x[, 1:2])
#Adding metadata information
pca_df <- merge(pca_df, metadata, by = "row.names")
pca_df$Row.names <- NULL

#Plotting
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +                                                                                    #Adds the data as points
  labs(title = "Explorative analysis - PCA",
       x = paste0("PCA1 (", round(100 * summary(pca)$importance[2,1], 1), "% of the variance)"),    #Title and labels
       y = paste0("PCA2 (", round(100 * summary(pca)$importance[2,2], 1), "% of the variance)")) +
  theme_bw() +      
  geom_line(aes(group = pair), color = "grey70", alpha = 0.4) + #Setting white background
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +                   #Formatting title
  scale_color_manual(values = c("case" = "red", "ctrl" ="darkblue"))                                     #Adding legend


## Heatmap
#Again using the counts after transposition: counts_norm_t
#Filtering: Heatmap is too heavy
g_vars <- rowVars(counts_norm) 
counts_norm_filt <- counts_norm[order(g_vars, decreasing=TRUE)[1:1000], ]

#Preparing data
counts_norm_t_filt <- t(counts_norm_filt)
counts_norm_matrix <- scale(counts_norm_t_filt)
heatmap_matrix <- t(counts_norm_matrix)

Conditions <- data.frame(
  Condition  = metadata$condition,
  Pair = metadata$pair
)
row.names(Conditions) <- colnames(heatmap_matrix)

#Plotting
pheat <- pheatmap(heatmap_matrix,
                  cluster_rows = FALSE,
                  cluster_cols = TRUE,
                  annotation_col = Conditions,
                  show_rownames = FALSE,
                  show_colnames = TRUE,
                  main = "Explorative analysis - Heatmap")


## Boxplot for Sample
#For the boxplot using the counts are needed in the "long" format: transforming
counts_norm_long <- counts_norm %>% 
  as.data.frame() %>% 
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "expr")
counts_norm_long$Sample <- factor(counts_norm_long$Sample, levels = paste0("S", 1:100))

#Plotting
ggplot(counts_norm_long, aes(x = Sample, y = expr)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Explorative analysis - Expression boxplots",
       x = "Sample",
       y = "VST gene expression")



### 4. Significant genes extraction and visualization ###

## Extracting Differentially expressed genes (by adj_pvalue)
#Considering the contrast on the condition: "case"
DE_res <- as.data.frame(results(dds, contrast = c("condition", "case", "ctrl")))
#Mapping the gene names as HUGO symbols
DE_res$'hugo_symbol' <- rownames(DE_res)


#Obtaining significant results
sign_DE_res <- DE_res %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)

## Exporting results
write.csv(sign_DE_res, file = file.path(PATH, "Dataset_0", "1_RNA-Seq", "DE_results.csv"), row.names = FALSE)

## Exporting the universe
RNAseq_universe <- DE_res[, c("hugo_symbol", "log2FoldChange", "padj")]
write.csv(RNAseq_universe, file = file.path(PATH, "Dataset_0", "1_RNA-Seq", "RNAseq_universe.csv"), row.names = FALSE)


## Visualization: volcano plot

ggplot(RNAseq_universe, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_point(data = RNAseq_universe[!is.na(RNAseq_universe$padj) & 
                                    RNAseq_universe$padj < 0.05 & 
                                    abs(RNAseq_universe$log2FoldChange) > 1, ], color = "red") +
  theme_minimal() +
  labs(title = "Volcano plot: DE analysis", x = "log2 Fold Change", y = "-log10(adj pvalue)")


## Output data for comparison with ground truth
Output_data <- data.frame(metric = c("M2 genes", "DE_genes"),
                          value = c(length(DE_res$padj), length(sign_DE_res$padj)))
write.csv(Output_data, file.path(PATH, "Dataset_0", "Ground Truth", "Output_data_tmp_1.csv"), row.names = FALSE)

dev.off()