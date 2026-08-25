library(rtracklayer)
library(tidyverse)
library(tximport)
library(DESeq2)
library(matrixStats)
library(pheatmap)


PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"

pdf(file.path(PATH, "Dataset_2", "1_RNA-Seq", "RNA-Seq Graphs D2.pdf"), height = 10, width = 15)


### 1. File and reference download ###

## Downloading the reference: hg38 version 49
gtf <- import(file.path(PATH, "references", "gencode.v49.annotation.gtf.gz"))
#Mapping ensembl id and hugo symbols
gtf_genes <- gtf[gtf$type == "gene"]
gene_map <- data.frame(ensembl_id = gtf_genes$gene_id,
                       hugo_symbol = gtf_genes$gene_name) #This will be needed later

#Clean mapping
#Maybe this is not needed, it's just to be sure there is a 1:1 mapping, and no NA.
gene_map <- gene_map %>%
  distinct() %>%
  group_by(ensembl_id) %>%
  summarise(hugo_symbol = dplyr::first(hugo_symbol), .groups = "drop") %>%
  mutate(hugo_symbol = ifelse(is.na(hugo_symbol) | hugo_symbol == "",
                              ensembl_id, hugo_symbol))
gene_map$ensembl_id <- sub("\\.\\d+$", "", gene_map$ensembl_id)  #Stripping .1, .2 etc..

#Downloading intron reference
t2g <- read.table(file.path(PATH, "references", "t2g_spliced_intron.tsv"),
                  header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("transcript", "gene"))

#Removing duplicates (if any)
t2g <- distinct(t2g)


## Downloading the count files
files <- c(
  # === Olig2 Control (17) ===
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2005_quant_D2", "quant.sf"),   # 1524
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2006_quant_D2", "quant.sf"),   # 1525
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2008_quant_D2", "quant.sf"),   # 1527
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2012_quant_D2", "quant.sf"),   # 1532
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2013_quant_D2", "quant.sf"),   # 1536
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2014_quant_D2", "quant.sf"),   # 1539
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2015_quant_D2", "quant.sf"),   # 1541
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2016_quant_D2", "quant.sf"),   # 3545
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2017_quant_D2", "quant.sf"),   # 3586
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2018_quant_D2", "quant.sf"),   # 3590
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2019_quant_D2", "quant.sf"),   # 3602
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2020_quant_D2", "quant.sf"),   # 4615
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2021_quant_D2", "quant.sf"),   # AN03398
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2022_quant_D2", "quant.sf"),   # AN05483
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2023_quant_D2", "quant.sf"),   # AN10090
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2024_quant_D2", "quant.sf"),   # AN15240
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_olig", "sample_2025_quant_D2", "quant.sf"),   # AN16799
  
  # === NeuN Control (17) ===
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_61_quant_D2", "quant.sf"),     # 1524
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_62_quant_D2", "quant.sf"),     # 1525
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_64_quant_D2", "quant.sf"),     # 1527
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_68_quant_D2", "quant.sf"),     # 1532
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_69_quant_D2", "quant.sf"),     # 1536
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_71_quant_D2", "quant.sf"),     # 1539
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_73_quant_D2", "quant.sf"),     # 1541
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_74_quant_D2", "quant.sf"),     # 3545
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_75_quant_D2", "quant.sf"),     # 3586
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_76_quant_D2", "quant.sf"),     # 3590
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_77_quant_D2", "quant.sf"),     # 3602
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_79_quant_D2", "quant.sf"),     # 4615
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_80_quant_D2", "quant.sf"),     # AN03398
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_81_quant_D2", "quant.sf"),     # AN05483
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_82_quant_D2", "quant.sf"),     # AN10090
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_83_quant_D2", "quant.sf"),     # AN15240
  file.path(PATH, "Dataset_2", "1_RNA-Seq", "salmon_output_neun", "sample_85_quant_D2", "quant.sf")      # AN16799
)
names(files) <- c(paste0("rep_", 1:17, "_olig2"),
                  paste0("rep_", 1:17, "_neun"))

#Importing counts
tx_matrix <- tximport(files,
                      type = "salmon",
                      tx2gene = t2g, # dataframe: txID, geneID
                      txOut = FALSE, 
                      ignoreAfterBar = TRUE)

#Uniting spliced and intron rows:
counts_total <- rowsum(tx_matrix$counts,
                       group = sub("-I$", "", rownames(tx_matrix$counts)))



## Creating the metadata table
metadata <- data.frame(
  sample = c(paste0("rep_", 1:17, "_olig2"),
             paste0("rep_", 1:17, "_neun")),
  group = c(rep("olig2", 17), rep("neun", 17)),
  donor = c(
    # olig (17): 
    1:17,
    # neun (17):
    1:17
  ),
  stringsAsFactors = FALSE
)

metadata$group <- factor(metadata$group, levels = c("neun", "olig2"))
metadata$donor <- factor(metadata$donor, levels = 1:17)
rownames(metadata) <- c(paste0("rep_", 1:17, "_olig2"),
                        paste0("rep_", 1:17, "_neun"))




### 2. DE analisys: DESeq2 ###

dds <- DESeqDataSetFromMatrix(countData = round(counts_total),
                                colData = metadata,
                                design = ~ donor + group)

#Filtering: keeping only genes with at least 10 counts in half of the samples
keep <- rowSums(counts(dds) >= 10) >= 17
dds <- dds[keep, ]

## DESeq2:
dds <- DESeq(dds)


## dds export: 
#Working on a parallel dds, and modify only this for the export
dds_exp <- dds   

#Stripping version
rownames(dds_exp) <- sub("\\.\\d+$", "", rownames(dds_exp))

#Changing the dds row names with SYMBOLS
new_names <- gene_map$hugo_symbol[match(rownames(dds_exp), gene_map$ensembl_id)]
new_names[is.na(new_names)] <- rownames(dds_exp)[is.na(new_names)]
new_names <- make.unique(new_names)      #make.unique will add ".1", ".2", etc..
rownames(dds_exp) <- new_names

#exporting
saveRDS(dds_exp, file = file.path(PATH, "Dataset_2", "1_RNA-Seq", "dds.rds"))



### 3. Explorative analysis ###

## For this experiment there are 17 samples in a paired set-up. For initial visualization only PCA will show clustering and outliers (if any). 
## Heatmap will confirm clustering and add information on expression. 
## Then boxplots of overall expression for each sample will show whether any samples have a general over/under expression, and will verify that normalization was successful

vst_dds <- vst(dds, blind = FALSE)                      #Normalization: variance stabilizing transformation
counts_norm <- assay(vst_dds)                           #Extracting counts
counts_norm <- counts_norm[rowVars(counts_norm) > 0, ]  #Deleting rows with variance = 0
counts_norm_t <- t(counts_norm)                         #The matrix is needed transposed for pca

#PCA:
pca <- prcomp(counts_norm_t)
summary(pca)

## PCA graph:
#Creating the plot dataframe
pca_df <- as.data.frame(pca$x[, 1:2])
#Adding metadata information
pca_df <- merge(pca_df, metadata, by = "row.names")
pca_df$Row.names <- NULL

#Plotting:
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +                                                                                    #Adds the data as points
  labs(title = "Explorative analysis - PCA",
       x = paste0("PCA1 (", round(100 * summary(pca)$importance[2,1], 1), "% of the variance)"),    #Title and labels
       y = paste0("PCA2 (", round(100 * summary(pca)$importance[2,2], 1), "% of the variance)")) +
  theme_bw() +                                                                                      #Setting white background
  geom_line(aes(group = donor), color = "grey70", alpha = 0.4) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +                   #Formatting title                                                     #Setting points shapes. From 15 on shapes are full inside
  scale_color_manual(values = c("neun" = "red", "olig2" ="darkblue"))                                    #Adding legend


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
  Donor  = metadata$donor,
  Group = metadata$group
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


## Boxplot for sample
#For the boxplot the counts are needed in the "long" format: transforming
counts_norm_long <- counts_norm %>% 
  as.data.frame() %>% 
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "expr")

counts_norm_long$Sample <- factor(counts_norm_long$Sample, levels = c(paste0("rep_", 1:17, "_olig2"),
                                                                      paste0("rep_", 1:17, "_neun")))

#Plotting:
ggplot(counts_norm_long, aes(x = Sample, y = expr)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Explorative analysis - Expression boxplots",
       x = "Sample",
       y = "VST gene expression")



### 4. Significant genes extraction and visualization ###

#Extracting Differentially expressed genes (by adj_pvalue). Considering the contrast on the condition: group
DE_res <- as.data.frame(results(dds, contrast = c("group", "olig2", "neun")))

#Mapping the gene names as HUGO symbols
DE_res$ENSEMBL <- sub("\\.\\d+$", "", rownames(DE_res))  #Stripping .1, .2 etc..
DE_res <- merge(DE_res, gene_map,
                by.x = "ENSEMBL",
                by.y = "ensembl_id",
                all.x = TRUE)


#Obtaining significant results
sign_DE_res <- DE_res %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)

## Exporting results
write.csv(sign_DE_res, file.path(PATH, "Dataset_2", "1_RNA-Seq", "DE_results.csv"), row.names = FALSE)

## Exporting the universe
RNAseq_universe <- DE_res[, c("hugo_symbol", "log2FoldChange", "padj")]
write.csv(RNAseq_universe, file = file.path(PATH, "Dataset_2", "1_RNA-Seq", "RNAseq_universe.csv"), row.names = FALSE)


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