library(rtracklayer)
library(tidyverse)
library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)


setwd("C:/Users/pierp/Desktop/Thesis PROJECT/Dataset_2/1_RNA-Seq")
pdf("C:/Users/pierp/Desktop/Thesis PROJECT/Dataset_2/1_RNA-Seq/RNA-Seq Graphs.pdf", height = 10, width = 15)

### 1. File and reference download ###

## Downloading the reference: hg38 version 49
gtf <- import("C:/Users/pierp/Desktop/Thesis PROJECT/references/gencode.v49.annotation.gtf.gz")
#mapping ensembl id and hugo symbols
gtf_genes <- gtf[gtf$type == "gene"]
gene_map <- data.frame(ensembl_id = gtf_genes$gene_id,
                       hugo_symbol = gtf_genes$gene_name)


## Creating the gene-transcript correspondence
#Selecting only "transcripts
gtf_tx <- gtf[gtf$type == "transcript"]
#Creating the table (Using id, not names, because they are unique)
tx2gene <- data.frame(transcript = gtf_tx$transcript_id,   
                      gene = gtf_tx$gene_id)
#Removing duplicates (if any)
tx2gene <- distinct(tx2gene)


## Downloading the count files
files <- paste0("sample_", setdiff(1988:2026, c(1989,2007,2009,2010,2011)), "_quant_D2/quant.sf")
names(files) <- paste0("S", setdiff(1988:2026, c(1989,2007,2009,2010,2011)))

#Importing counts
tx_matrix <- tximport(files,
                      type = "salmon",
                      tx2gene = tx2gene, # dataframe: txID, geneID
                      txOut = FALSE, 
                      ignoreAfterBar = TRUE)


## Creating the metadata table
#Adding metadata for covariates
metadata <- data.frame(
  sample = c(paste0("rep_", 1:16, "_schizo"),
             paste0("rep_", 1:17, "_ctrl")),
  group = c(rep("schizo", 16), rep("ctrl", 17)),
  sex = c(
    # 16 schizo: 
    "F","M","F","M","F","M","F","F","M","F","F","F","M","M","F","M",
    # 17 ctrl:
    "M","F","M","F","M","F","F","M","M","M","M","M","F","M","M","F","M"
  ),
  stringsAsFactors = FALSE
)

rownames(metadata) <- paste0("S", setdiff(1988:2026, c(1989,2007,2009,2010,2011)))

metadata$group     <- factor(metadata$group,     levels = c("ctrl", "schizo"))
metadata$sex       <- factor(metadata$sex)

#Output metadata
write.csv(metadata, "C:/Users/pierp/Desktop/Thesis PROJECT/Dataset_2/metadata.csv")


### 2. DE analisys: DESeq2 ###

## !The covariates are a bit unbalanced, but not extremely. So they can be inserted in the design
dds <- DESeqDataSetFromTximport(tx_matrix,
                                colData = metadata,
                                design = ~ sex + group)
#Filtering: keeping only genes with 10 counts in at least one of the conditions (half of the samples)
keep <- rowSums(counts(dds) >= 19) >= 10
dds <- dds[keep, ]


##DESeq2:
dds <- DESeq(dds)


## Export: 
#Before exporting the dds it's important to shift the names from ENSEMBL ID to SYMBOLS. This will be useful later.
dds_exp <- dds   #I'll work on a parallel dds, and modify only this for the export

#Clean mapping
#Maybe this is not needed, it's just to be sure there is a 1:1 mapping, and no NA.
gene_map_clean <- gene_map %>%
  distinct() %>%
  group_by(ensembl_id) %>%
  summarise(hugo_symbol = first(hugo_symbol), .groups = "drop") %>%
  mutate(hugo_symbol = ifelse(is.na(hugo_symbol) | hugo_symbol == "",
                              ensembl_id, hugo_symbol))

#Stripping version
rownames(dds_exp) <- sub("\\.\\d+$", "", rownames(dds_exp))
gene_map_clean$ensembl_id <- sub("\\.\\d+$", "", gene_map_clean$ensembl_id)

#Changing the dds row names with SYMBOLS
new_names <- gene_map_clean$hugo_symbol[
  match(rownames(dds_exp), gene_map_clean$ensembl_id)]
new_names[is.na(new_names)] <- rownames(dds_exp)[is.na(new_names)]
new_names <- make.unique(new_names)      #make.unique will add ".1", ".2", etc..
rownames(dds_exp) <- new_names

##Export
saveRDS(dds_exp, file = "C:/Users/pierp/Desktop/Thesis PROJECT/Dataset_2/1_RNA-Seq/dds.rds")



### 3. Explorative analysis ###

## For this experiment there are 6 samples in a paired set-up. For initial visualization only PCA will show clustering and outliers (if any). 
## t-SNE and UMAP need more sample to start being informative so they will not be implemented. 
## Heatmap will confirm clustering and add information on expression. 
## Then boxplots of overall expression for each sample will show whether any samples have a general over/under expression, and will verify that normalization was successful

vst_dds <- vst(dds, blind = FALSE)      #Normalization: variance stabilizing transformation
counts_norm <- assay(vst_dds)           #Extracting counts
counts_norm <- counts_norm[rowVars(counts_norm) > 0, ]  #Deleting rows with variance = 0
counts_norm_t <- t(counts_norm)         #The matrix is needed transposed for pca

##PCA:
pca <- prcomp(counts_norm_t)
summary(pca)

##PCA graph:
#Creating the plot dataframe
pca_df <- as.data.frame(pca$x[, 1:2])
#Adding metadata information
pca_df <- merge(pca_df, metadata, by = "row.names")
pca_df$Row.names <- NULL

#ggplot
ggplot(pca_df, aes(x = PC1, y = PC2, shape = sex, color = group)) +
  geom_point(size = 3) +                                                                                    #Adds the data as points
  labs(title = "Explorative analysis - PCA",
       x = paste0("PCA1 (", round(100 * summary(pca)$importance[2,1], 1), "% of the variance)"),    #Title and labels
       y = paste0("PCA2 (", round(100 * summary(pca)$importance[2,2], 1), "% of the variance)")) +
  theme_bw() +                                                                                      #Setting white background
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +                   #Formatting title
  scale_shape_manual(values= c(15:19, 8)) +                                                         #Setting points shapes. From 15 on shapes are full inside
  scale_color_manual(values = c("ctrl" = "red", "schizo" ="darkblue")) +                                #Setting colors
  guides(fill = guide_legend(override.aes = list(shape = 21)))                                      #Adding legend


## Heatmap
#Again using the counts after transposition: counts_norm_t
counts_norm_matrix <- scale(counts_norm_t)
heatmap_matrix <- t(counts_norm_matrix)

Conditions <- data.frame(
  group  = metadata$group,
  Sex = metadata$sex)
row.names(Conditions) <- colnames(heatmap_matrix)

pheat <- pheatmap(heatmap_matrix,
                  cluster_rows = FALSE,
                  cluster_cols = TRUE,
                  annotation_col = Conditions,
                  show_rownames = FALSE,
                  show_colnames = TRUE,
                  main = "Explorative analysis - Heatmap")


##Boxplot for Sample
#For the boxplot using the counts are needed in the "long" format: transforming
counts_norm_long <- counts_norm %>% 
  as.data.frame() %>% 
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "expr")
counts_norm_long$Sample <- factor(counts_norm_long$Sample, levels = paste0("S", setdiff(1988:2026, c(1989,2007,2009,2010,2011))))

#Boxplot
ggplot(counts_norm_long, aes(x = Sample, y = expr)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Explorative analysis - Expression boxplots",
       x = "Sample",
       y = "VST gene expression")

## Results
## PCA -> Strong Clustering for condition (group). There is also a slight clustering for sample (on PCA2): paired design is correct.
## Heatmap -> Confirms clustering seen for PCA
## boxplot -> From the boxplot profiles it appears that: 1) Normalization was correctly done. 2) Overall expression is coherent among different samples: no sample is an anomaly. 
## 3) All samples have an expected distribution of gene expression, (with a peak for housekeeping genes)  



### 4. Significant genes extraction and visualization ###

##Extracting Differentially expressed genes (by adj_pvalue)
#Considering the contrast on the condition: group
DE_res <- as.data.frame(results(dds, contrast = c("group", "schizo", "ctrl")))
#Mapping the gene names as HUGO symbols
DE_res$ENSEMBL <- row.names(DE_res)
DE_res <- merge(DE_res, gene_map,
                by.x = "ENSEMBL",
                by.y = "ensembl_id",
                all.x = TRUE)

#Removing duplicates (there are some, as shown by "sum(duplicated(sign_DE_res$hugo_symbol))")
DE_res <- DE_res %>%
  group_by(hugo_symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

sign_DE_res <- DE_res %>% filter(pvalue < 0.05 & abs(log2FoldChange) > 1)
#Exporting results
write.csv(sign_DE_res, "C:/Users/pierp/Desktop/Thesis PROJECT/Dataset_2/1_RNA-Seq/DE_results.csv", row.names = FALSE)


#Exporting the universe
RNAseq_universe <- DE_res[, c("hugo_symbol", "log2FoldChange", "padj")]
write.csv(RNAseq_universe, file = "C:/Users/pierp/Desktop/Thesis PROJECT/Dataset_2/1_RNA-Seq/RNAseq_universe.csv", row.names = FALSE)



##Visualization: volcano plot

ggplot(RNAseq_universe, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_point(data = RNAseq_universe[!is.na(RNAseq_universe$padj) & 
                                    RNAseq_universe$padj < 0.05 & 
                                    abs(RNAseq_universe$log2FoldChange) > 1, ], color = "red") +
  theme_minimal() +
  labs(title = "Volcano plot: group", x = "log2 Fold Change", y = "-log10(adj pvalue)")

dev.off()







#Deleted
#Adding metadata for covariates
# metadata <- data.frame(
#   sample = c(paste0("rep_", 1:21, "_schizo"),
#              paste0("rep_", 1:19, "_ctrl")),
#   group = c(rep("schizo", 21), rep("ctrl", 19)),
#   sex = c(
#     # 21 schizo
#     "F","F","M","F","M","M","F","F","M","F","F","F","F","F","M","M","F","M","M","F","M",
#     # 19 ctrl
#     "M","F","M","F","M","F","F","F","M","M","M","M","M","M","F","M","M","F","M"
#   ),
#   age = c(
#     # schizo
#     46,85,32,73,33,26,43,73,40,50,52,42,85,77,56,46,32,61,26,56,66,
#     # ctrl
#     40,65,36,68,41,55,53,49,80,76,75,66,64,49,75,66,52,36,43
#   ),
#   PMI = c(
#     # schizo
#     23,9,27,36,20,24,11,14,9,27.15,27.5,15.5,11.5,14.7,14.7,21.7,12.3,28,16,10.5,16.47,
#     # ctrl
#     10,11,23,19,24,25,23,15.3,14,16,11.5,13.2,17.5,15,12.1,16.97,13.12,18.08,14.68
#   ),
#   RIN = c(
#     # schizo: 1505,1506,1507,1508,1510,1511,1512,1513,1514,1515,1518,1523,4336,4361,4395,4448,4730,4804,AN09634,AN17799,AN18099
#     9.2,8.4,9.9,7.4,9.0,10.0,6.7,6.5,9.4,8.5,8.1,8.3,4.3,5.9,7.1,7.1,6.7,5.8,8.0,6.1,5.1,
#     # ctrl: 1524,1525,1527,1532,1536,1538,1539,1541,3545,3586,3590,3602,3611,4615,AN03398,AN05483,AN10090,AN15240,AN16799
#     8.4,8.3,9.8,8.8,9.3,7.6,7.6,8.6,4.8,6.4,5.5,6.6,2.7,7.2,5.1,6.5,6.6,7.9,5.8
#   ),
#   hemisphere = c(
#     # schizo
#     "L","R","R","R","L","L","R","L","R","R","L","R","L","L","L","L","L","L","R","L","L",
#     # ctrl
#     "R","L","L","R","R","R","R","L","L","L","L","L","R","L","L","L","R","L","R"
#   ),
#   brainbank = c(
#     # schizo: 12 Ctlab, 6 HSB, 3 HBTC
#     rep("Ctlab", 12), rep("HSB", 6), rep("HBTC", 3),
#     # ctrl: 8 Ctlab, 6 HSB, 5 HBTC
#     rep("Ctlab", 8),  rep("HSB", 6), rep("HBTC", 5)
#   ),
#   stringsAsFactors = FALSE
# )
# 
# rownames(metadata) <- paste0("S", setdiff(38:85, c(49, 55, 63, 65, 66, 67, 72, 84)))
# metadata$RIN <- scale(metadata$RIN, center = TRUE, scale = FALSE)[, 1]
# metadata$group     <- factor(metadata$group,     levels = c("ctrl", "schizo"))
# metadata$sex       <- factor(metadata$sex)
# metadata$brainbank <- factor(metadata$brainbank)