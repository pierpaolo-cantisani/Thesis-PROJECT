library(rtracklayer)
library(matrixStats)
library(tidyverse)
library(DESeq2)
library(pheatmap)


PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"

pdf(file.path(PATH, "Dataset_3", "1_RNA-Seq", "RNA-Seq Graphs D3.pdf"), height = 10, width = 15)

### 1. File and reference download ###

## Downloading the reference: hg38 version 49
gtf <- import(file.path(PATH, "references", "gencode.v49.annotation.gtf.gz"))
#Mapping ensembl id and hugo symbols
gtf_genes <- gtf[gtf$type == "gene"]
gene_map <- data.frame(
  ensembl_id  = sub("\\..*$", "", gtf_genes$gene_id),  # stripped version
  hugo_symbol = gtf_genes$gene_name,
  stringsAsFactors = FALSE) |> dplyr::distinct()


#Importing counts
counts <- read.delim(file.path(PATH, "Dataset_3", "1_RNA-Seq", "Additional_f6_ST4.txt"),
                     header = TRUE, row.names = 1,
                     check.names = FALSE)

#Filtering: keeping only those in common with BS-Seq:
keep_subj <- c("1","14","16","17","22","34","47","58","64","65",
               "149","151","155","161","166","169","170","171","177","179",
               "184","185","188","190","194","196","203","207","213","214",
               "215","219","220","221","224","227","238","244","249","260",
               "395","456")
sample_subject <- str_extract(colnames(counts), "^[^_]+")
keep_cols <- sample_subject %in% keep_subj
counts <- counts[, keep_cols]

#Sorting numerically (rather than alphabetically)
subj_num <- as.integer(str_extract(colnames(counts), "^[^_]+"))
tp       <- factor(str_extract(colnames(counts), "D0$|D28$"), levels = c("D0","D28"))

#Ordering: first numeric, then D0 < D28
ord <- order(subj_num, tp)
counts <- counts[, ord]

#As matrix for DESeq2
counts_mat <- as.matrix(counts)
#Checking:
stopifnot(typeof(counts_mat) == "integer")



## Creating the metadata table
sample_ids <- colnames(counts_mat)
metadata <- data.frame(
  condition = factor(str_extract(sample_ids, "D0$|D28$"), levels = c("D0","D28")),
  donor   = factor(str_extract(sample_ids, "^[^_]+"),
                     levels = as.character(keep_subj)),
  row.names = sample_ids
)



### 2. DE analisys: DESeq2 ###
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = metadata,
                              design = ~ donor + condition)

#Filtering: keeping only genes with at least 10 counts in half of the samples
keep <- rowSums(counts(dds) >= 10) >= 42
dds <- dds[keep, ]

## DESeq2:
dds <- DESeq(dds)


## Export: 

dds_exp <- dds   #Working on a parallel dds, and modify only this for the export

#Stripping version
rownames(dds_exp) <- sub("\\.\\d+$", "", rownames(dds_exp))

#Changing the dds row names with SYMBOLS
new_names <- gene_map$hugo_symbol[
  match(rownames(dds_exp), gene_map$ensembl_id)]
new_names[is.na(new_names)] <- rownames(dds_exp)[is.na(new_names)]
new_names <- make.unique(new_names)      #make.unique will add ".1", ".2", etc..
rownames(dds_exp) <- new_names

#Exporting
saveRDS(dds_exp, file = file.path(PATH, "Dataset_3", "1_RNA-Seq", "dds.rds"))



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
  geom_point(size = 3) +                                                                            
  labs(title = "Explorative analysis - PCA",
       x = paste0("PCA1 (", round(100 * summary(pca)$importance[2,1], 1), "% of the variance)"),    #Title and labels
       y = paste0("PCA2 (", round(100 * summary(pca)$importance[2,2], 1), "% of the variance)")) +
  theme_bw() +                                                                                      #Setting white background
  geom_line(aes(group = donor), color = "grey70", alpha = 0.4) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +                   #Formatting title
  scale_color_manual(values = c("D0" = "red", "D28" ="darkblue"))                                   #Setting colors


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
  Donor = metadata$donor
)
row.names(Conditions) <- colnames(heatmap_matrix)

#Plotting:
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
counts_norm_long$Sample <- factor(counts_norm_long$Sample, levels = colnames(counts_mat))

#Plotting:
ggplot(counts_norm_long, aes(x = Sample, y = expr)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Explorative analysis - Expression boxplots",
       x = "Sample",
       y = "VST gene expression")



### 4. Significant genes extraction and visualization ###

##Extracting Differentially expressed genes (by adj_pvalue)

#Considering the contrast on the condition: D28
DE_res <- as.data.frame(results(dds, contrast = c("condition", "D28", "D0")))

#Mapping the gene names as HUGO symbols
DE_res$ENSEMBL <- sub("\\.\\d+$", "", rownames(DE_res))  #Stripping .1, .2 etc..
DE_res <- merge(DE_res, gene_map,
                by.x = "ENSEMBL",
                by.y = "ensembl_id",
                all.x = TRUE)


#Obtaining significant results
sign_DE_res <- DE_res %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)

## Exporting results
write.csv(sign_DE_res, file.path(PATH, "Dataset_3", "1_RNA-Seq", "DE_results.csv"), row.names = FALSE)

## Exporting the universe
RNAseq_universe <- DE_res[, c("hugo_symbol", "log2FoldChange", "padj")]
write.csv(RNAseq_universe, file = file.path(PATH, "Dataset_3", "1_RNA-Seq", "RNAseq_universe.csv"), row.names = FALSE)


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

dev.off()