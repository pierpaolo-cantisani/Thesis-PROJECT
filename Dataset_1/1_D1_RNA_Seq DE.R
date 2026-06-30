library(rtracklayer)
library(tidyverse)
library(tximport)
library(DESeq2)
library(matrixStats)
library(pheatmap)


PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"
pdf(file.path(PATH, "Dataset_1", "1_RNA-Seq", "RNA-Seq Graphs D1.pdf"), height = 10, width = 15)

### 1. File and reference download ###

## Downloading the reference: hg38 version 49
gtf <- import(file.path(PATH, "references", "gencode.v49.annotation.gtf.gz"))
#mapping ensembl id and hugo symbols
gtf_genes <- gtf[gtf$type == "gene"]
gene_map <- data.frame(ensembl_id = gtf_genes$gene_id,
                       hugo_symbol = gtf_genes$gene_name)
#Clean mapping
#Maybe this is not needed, it's just to be sure there is a 1:1 mapping, and no NA.
gene_map <- gene_map %>%
  distinct() %>%
  group_by(ensembl_id) %>%
  summarise(hugo_symbol = first(hugo_symbol), .groups = "drop") %>%
  mutate(hugo_symbol = ifelse(is.na(hugo_symbol) | hugo_symbol == "",
                              ensembl_id, hugo_symbol))


## Creating the gene-transcript correspondence
#Selecting only "transcripts
gtf_tx <- gtf[gtf$type == "transcript"]
#Creating the table (Using id, not names, because they are unique)
tx2gene <- data.frame(transcript = gtf_tx$transcript_id,   
                      gene = gtf_tx$gene_id)
#Removing duplicates (if any)
tx2gene <- distinct(tx2gene)

#Setting file directory
setwd(file.path(PATH, "Dataset_1", "1_RNA-Seq"))

## Downloading the count files
files <- c("rep_1_TB_quant/quant.sf", "rep_1_NI_quant/quant.sf", "rep_2_TB_quant/quant.sf",
           "rep_2_NI_quant/quant.sf", "rep_3_TB_quant/quant.sf", "rep_3_NI_quant/quant.sf",
           "rep_4_TB_quant/quant.sf", "rep_4_NI_quant/quant.sf", "rep_5_TB_quant/quant.sf",
           "rep_5_NI_quant/quant.sf", "rep_6_TB_quant/quant.sf", "rep_6_NI_quant/quant.sf")
names(files) <- paste0("S", 1:12)

#Importing counts
tx_matrix <- tximport(files,
                      type = "salmon",
                      tx2gene = tx2gene, # dataframe: txID, geneID
                      txOut = FALSE, 
                      ignoreAfterBar = TRUE)


## Creating the metadata table
metadata <- data.frame(infection = c("TB", "CTRL", "TB", "CTRL", "TB", "CTRL", "TB", "CTRL", "TB", "CTRL", "TB", "CTRL"),
                       replicate = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5", "6", "6"))
rownames(metadata) <- paste0("S", 1:12)
metadata$infection <- factor(metadata$infection)   #CTRL is first
metadata$replicate <- factor(metadata$replicate)




### 2. DE analisys: DESeq2 ###
dds <- DESeqDataSetFromTximport(tx_matrix,
                                colData = metadata,
                                design = ~ replicate + infection)
#Filtering: keeping only genes with at least 10 counts in half of the samples
keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep, ]

##DESeq2:
dds <- DESeq(dds)


## Export: 
#Before exporting the dds it's important to shift the names from ENSEMBL ID to SYMBOLS. This will be useful later.
dds_exp <- dds   #I'll work on a parallel dds, and modify only this for the export



#Stripping version
rownames(dds_exp) <- sub("\\.\\d+$", "", rownames(dds_exp))
gene_map$ensembl_id <- sub("\\.\\d+$", "", gene_map$ensembl_id)

#Changing the dds row names with SYMBOLS
new_names <- gene_map$hugo_symbol[
  match(rownames(dds_exp), gene_map$ensembl_id)]
new_names[is.na(new_names)] <- rownames(dds_exp)[is.na(new_names)]
new_names <- make.unique(new_names)      #make.unique will add ".1", ".2", etc..
rownames(dds_exp) <- new_names

##Export
saveRDS(dds_exp, file = file.path(PATH, "Dataset_1", "1_RNA-Seq", "dds.rds"))



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
ggplot(pca_df, aes(x = PC1, y = PC2, shape = replicate, color = infection)) +
  geom_point(size = 3) +                                                                                    #Adds the data as points
  labs(title = "Explorative analysis - PCA",
       x = paste0("PCA1 (", round(100 * summary(pca)$importance[2,1], 1), "% of the variance)"),    #Title and labels
       y = paste0("PCA2 (", round(100 * summary(pca)$importance[2,2], 1), "% of the variance)")) +
  theme_bw() +      
  geom_line(aes(group = replicate), color = "grey70", alpha = 0.4) + #Setting white background
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +                   #Formatting title
  scale_shape_manual(values= c(15:19, 8)) +                                                         #Setting points shapes. From 15 on shapes are full inside
  scale_color_manual(values = c("CTRL" = "red", "TB" ="darkblue"))                                     #Adding legend


## Heatmap
#Again using the counts after transposition: counts_norm_t
#Filtering: Heatmap is too heavy
g_vars <- rowVars(counts_norm) 
counts_norm_filt <- counts_norm[order(g_vars, decreasing=TRUE)[1:1000], ]

#now preparing data
counts_norm_t_filt <- t(counts_norm_filt)
counts_norm_matrix <- scale(counts_norm_t_filt)
heatmap_matrix <- t(counts_norm_matrix)

Conditions <- data.frame(
  Infection  = metadata$infection,
  Replicate = metadata$replicate
)
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
counts_norm_long$Sample <- factor(counts_norm_long$Sample, levels = paste0("S", 1:12))

#Boxplot
ggplot(counts_norm_long, aes(x = Sample, y = expr)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Explorative analysis - Expression boxplots",
       x = "Sample",
       y = "VST gene expression")

## Results
## PCA -> Strong Clustering for condition (infection). There is also a slight clustering for sample (on PCA2): paired design is correct.
## Heatmap -> Confirms clustering seen for PCA
## boxplot -> From the boxplot profiles it appears that: 1) Normalization was correctly done. 2) Overall expression is coherent among different samples: no sample is an anomaly. 
## 3) All samples have an expected distribution of gene expression, (with a peak for housekeeping genes)  



### 4. Significant genes extraction and visualization ###

##Extracting Differentially expressed genes (by adj_pvalue)
#Considering the contrast on the condition: infection
DE_res <- as.data.frame(results(dds, contrast = c("infection", "TB", "CTRL")))
#Mapping the gene names as HUGO symbols
DE_res$ENSEMBL <- sub("\\.\\d+$", "", rownames(DE_res)) #Stripping .1, .2 etc..
DE_res <- merge(DE_res, gene_map,
                by.x = "ENSEMBL",
                by.y = "ensembl_id",
                all.x = TRUE)

#Removing duplicates (there are some, as shown by "sum(duplicated(sign_DE_res$hugo_symbol))")
DE_res <- DE_res %>%
  group_by(hugo_symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

sign_DE_res <- DE_res %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
#Exporting results
write.csv(sign_DE_res, file = file.path(PATH, "Dataset_1", "1_RNA-Seq", "DE_results.csv"), row.names = FALSE)


#Exporting the universe
RNAseq_universe <- DE_res[, c("hugo_symbol", "log2FoldChange", "padj")]
write.csv(RNAseq_universe, file = file.path(PATH, "Dataset_1", "1_RNA-Seq", "RNAseq_universe.csv"), row.names = FALSE)



##Visualization: volcano plot

ggplot(RNAseq_universe, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_point(data = RNAseq_universe[!is.na(RNAseq_universe$padj) & 
                                    RNAseq_universe$padj < 0.05 & 
                                    abs(RNAseq_universe$log2FoldChange) > 1, ], color = "red") +
  theme_minimal() +
  labs(title = "Volcano plot: infection", x = "log2 Fold Change", y = "-log10(adj pvalue)")

dev.off()