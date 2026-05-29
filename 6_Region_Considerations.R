library(dplyr)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyr)
library(ggplot2)
library(cowplot)

## Variabile confondente da considerare!! Questa annotazione è con CHIPseeker (methods 2 e 3). Quindi questa analisi
#  potrebbe avere risultati "biased"

##! Metodi 1, 2 e 5 sono uguali. E metodi 3 e 4 sono sottoinsiemi di 1/2/5. Poco interessante


#Importing DE genes
DE_results <- read.csv("C:/Users/pierp/Desktop/THESIS PROJECT/1_RNA-Seq/dataset_1/DE_results.csv")
DE_results <- DE_results %>% dplyr::rename(SYMBOL = hugo_symbol)  #renaming for coherence

#Creating list for graphs and table
graph_list <- list()
Stats_region <- list()


for(METHOD in 1:5){
  
  #Importing association df
  DM_sites <- read.csv(sprintf("C:/Users/pierp/Desktop/THESIS PROJECT/3_Benchmark/dataset_1/DM_sites_Met%d.csv", METHOD))
  
  #Converting to GRanges
  DM_GR <- makeGRangesFromDataFrame(DM_sites,
                                    seqnames.field   = "seqnames",
                                    start.field      = "start",
                                    end.field        = "start", 
                                    keep.extra.columns = FALSE,           
                                    starts.in.df.are.0based = FALSE)
  seqlevelsStyle(DM_GR) <- "UCSC"
  
  ## Annotation: CpG that were part of an association for each method are reannotated with the same tool (annotatr): in this way
  #              they will be assigned to a gene region, and they will be comparable.
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  peakAnno <- annotatePeak(DM_GR,
                           tssRegion = c(-2000, 200),   #This is the standard definition for Promoter. Can be arbitrarly changed
                           TxDb      = txdb,
                           annoDb    = "org.Hs.eg.db")
  DM_ann_df <- as.data.frame(peakAnno)
  
  
  
  ### Region Analysis ###
  
  #Intersection:
  intersect_genes <- intersect(DE_results$SYMBOL, DM_ann_df$SYMBOL)
  
  ## DMR on promoter/exon/intron/UTR/intergenic genes
  #! exon contains 5' and 3' UTR. "others" contains intergenic and 1to5kb upstream
  prom_genes <- unique(DM_ann_df$SYMBOL[grepl("^Promoter", DM_ann_df$annotation)])
  exon_genes <- unique(DM_ann_df$SYMBOL[grepl("^Exon", DM_ann_df$annotation)])
  intron_genes <- unique(DM_ann_df$SYMBOL[grepl("^Intron", DM_ann_df$annotation)])
  interg_genes <- unique(DM_ann_df$SYMBOL[grepl("Distal Intergenic", DM_ann_df$annotation)])
  UTR5_genes <- unique(DM_ann_df$SYMBOL[grepl("5' UTR", DM_ann_df$annotation)])
  UTR3_genes <- unique(DM_ann_df$SYMBOL[grepl("3' UTR", DM_ann_df$annotation)])
  
  #Obtaining only exon/intron/UTR intersecting genes
  prom_inter_genes <- intersect(prom_genes, DE_results$SYMBOL)
  exon_inter_genes <- intersect(exon_genes, DE_results$SYMBOL)
  intron_inter_genes <- intersect(intron_genes, DE_results$SYMBOL)
  interg_inter_genes <- intersect(interg_genes, DE_results$SYMBOL)
  UTR5_inter_genes <- intersect(UTR5_genes, DE_results$SYMBOL)
  UTR3_inter_genes <- intersect(UTR3_genes, DE_results$SYMBOL)
  
  ## Some of the cases have very low numbers: must be careful with the test interpretation. 
  
  ## Hypergeometric tests
  #Are intersecting genes associated to methylation specific to introns/exons/3UTR/promoters?
  
  #Promoter
  N <- as.numeric(length(unique(DM_ann_df$SYMBOL)))             # all DM genes
  m <- as.numeric(length(prom_genes))                           # prom DM genes
  k <- as.numeric(length(intersect_genes))                      # all intersecting DM genes
  q <- as.numeric(length(prom_inter_genes))                     # intersecting prom on DM genes
  hyp_prom_trend <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  #exon
  N <- as.numeric(length(unique(DM_ann_df$SYMBOL)))             # all DM genes
  m <- as.numeric(length(exon_genes))                           # exon DM genes
  k <- as.numeric(length(intersect_genes))                      # all intersecting DM genes
  q <- as.numeric(length(exon_inter_genes))                     # intersecting exon DM genes
  hyp_ex <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  #intron
  N <- as.numeric(length(unique(DM_ann_df$SYMBOL)))             # all DM genes
  m <- as.numeric(length(intron_genes))                         # intron DM genes
  k <- as.numeric(length(intersect_genes))                      # all intersecting DM genes
  q <- as.numeric(length(intron_inter_genes))                   # intersecting intron DM genes
  hyp_int <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  #Intergenic
  N <- as.numeric(length(unique(DM_ann_df$SYMBOL)))             # all DM genes
  m <- as.numeric(length(interg_genes))                         # intergenic DM genes
  k <- as.numeric(length(intersect_genes))                      # all intersecting DM genes
  q <- as.numeric(length(interg_inter_genes))                   # intersecting interg on DM genes
  hyp_interg <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  #5' UTR
  N <- as.numeric(length(unique(DM_ann_df$SYMBOL)))             # all DM genes
  m <- as.numeric(length(UTR5_genes))                           # 5'UTR DM genes
  k <- as.numeric(length(intersect_genes))                      # all intersecting DM genes
  q <- as.numeric(length(UTR5_inter_genes))                     # intersecting 5'UTR on DM genes
  hyp_5UTR <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  #Intergenic
  N <- as.numeric(length(unique(DM_ann_df$SYMBOL)))              # all DM genes
  m <- as.numeric(length(UTR3_genes))                            # 3'UTR DM genes
  k <- as.numeric(length(intersect_genes))                       # all intersecting DM genes
  q <- as.numeric(length(UTR3_inter_genes))                      # intersecting 3'UTR on DM genes
  hyp_3UTR <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  
  perc_prom <- as.numeric(length(prom_inter_genes))*100/as.numeric(length(prom_genes))     #as.numeric is probably useless here
  perc_ex = as.numeric(length(exon_inter_genes))*100/as.numeric(length(exon_genes))
  perc_int = as.numeric(length(intron_inter_genes))*100/as.numeric(length(intron_genes))
  perc_interg = as.numeric(length(interg_inter_genes))*100/as.numeric(length(interg_genes))
  perc_5UTR = as.numeric(length(UTR5_inter_genes))*100/as.numeric(length(UTR5_genes))
  perc_3UTR = as.numeric(length(UTR3_inter_genes))*100/as.numeric(length(UTR3_genes))
  
  Stats_region[[METHOD]] <- c(paste0(round(perc_prom, 1)), hyp_prom_trend, paste0(round(perc_ex, 1)), hyp_ex, 
                              paste0(round(perc_int, 1)), hyp_int, paste0(round(perc_interg, 1)), hyp_interg,
                              paste0(round(perc_5UTR, 1)), hyp_5UTR, paste0(round(perc_3UTR, 1)), hyp_3UTR)
  names(Stats_region[[METHOD]]) <- c("Perc ∩ prom", "Sign prom (p)", "Perc ∩ exon", "Sign exon (p)", "Perc ∩ intron", 
                                     "Sign intron (p)", "Perc ∩ interg", "Sign interg (p)", "Perc ∩ 5'UTR", "Sign 5'UTR (p)",
                                     "Perc ∩ 3'UTR", "Sign 3'UTR (p)")
  
  
  
  ## Graph: Region type ##
  
  # Building counts
  region_summary <- data.frame(
    region = c("Promoter", "5'UTR", "Exon", "Intron", "3'UTR", "Intergenic"),
    DM_all = c(length(prom_genes), length(UTR5_genes), length(exon_genes), length(intron_genes), length(UTR3_genes), length(interg_genes)),
    inters_DM = c(length(prom_inter_genes), length(UTR5_inter_genes), length(exon_inter_genes), length(intron_inter_genes), length(UTR3_inter_genes), length(interg_inter_genes))
  )
  
  # Obtaining long format for ggplot
  region_long <- region_summary %>%
    pivot_longer(cols = c(DM_all, inters_DM),
                 names_to = "group",
                 values_to = "n") %>%
    group_by(group) %>%
    mutate(pct = n / sum(n) * 100,
           group = factor(group,
                          levels = c("DM_all", "inters_DM"),
                          labels = c("All DM genes", "DM ∩ DE genes")),
           region = factor(region, levels = c("Promoter", "5'UTR", "Exon", "Intron", "3'UTR", "Intergenic")))
  
  region_colors <- c(
    "Promoter" = "#E41A1C",
    "5'UTR"    = "darkgreen",
    "Exon"     = "#FF7F00",
    "Intron"   = "#4393C3",
    "3'UTR"    = "lightgreen",
    "Intergenic"   = "#984EA3"
  )
  
  # Tot of annotation on bar
  totals <- region_long %>%
    group_by(group) %>%
    summarise(total = sum(n), .groups = "drop")
  
  if (METHOD == 1) {
    graph_list[[METHOD]] <- ggplot(region_long, aes(x = group, y = pct, fill = region)) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_text(
        data = totals,
        aes(x = group, y = 102, label = paste0("n = ", total)),
        inherit.aes = FALSE,
        size = 4, fontface = "bold"
      ) +
      scale_fill_manual(values = region_colors, name = "Genomic region") +
      scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 20),
                         labels = function(x) paste0(x, "%")) +
      labs(
        title = "Method 1&2&5: Genomic region distribution",
        subtitle = "All DM genes vs DM ∩ DE genes",
        x = NULL,
        y = "Percentage of genes"
      ) +
      theme_classic(base_size = 13) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12))
  } else if (METHOD == 2 || METHOD == 5) {
    next
  } else {
    graph_list[[METHOD-1]] <- ggplot(region_long, aes(x = group, y = pct, fill = region)) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_text(
        data = totals,
        aes(x = group, y = 102, label = paste0("n = ", total)),
        inherit.aes = FALSE,
        size = 4, fontface = "bold"
      ) +
      scale_fill_manual(values = region_colors, name = "Genomic region") +
      scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 20),
                         labels = function(x) paste0(x, "%")) +
      labs(
        title = sprintf("Method %d: Genomic region distribution", METHOD),
        subtitle = "All DM genes vs DM ∩ DE genes",
        x = NULL,
        y = "Percentage of genes"
      ) +
      theme_classic(base_size = 13) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12))
  }
}

#Graph Output
pdf("C:/Users/pierp/Desktop/THESIS PROJECT/4_Integration_results/dataset_1/Region_analysis_ChIPseeker.pdf", height = 10, width = 15)
plot_grid(plotlist = graph_list, ncol = 2)
dev.off()

#Stats output
Stats_df <- as.data.frame(do.call(cbind, Stats_region))
colnames(Stats_df) <- sprintf("M%d", seq_along(Stats_region))
write.csv(Stats_df, "C:/Users/pierp/Desktop/THESIS PROJECT/4_Integration_results/dataset_1/Stats_table_region.csv", row.names = FALSE)