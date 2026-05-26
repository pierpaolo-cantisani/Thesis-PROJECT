library(dplyr)
library(GenomicRanges)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyr)
library(ggplot2)
library(cowplot)

## Variabile confondente da considerare!! Questa annotazione è con il ref TxDb, uguale a CHIPseeker (meth 2 e 3). Quindi quei metodi 
#  potrebbero avere annotazioni influenzate da ciò

##!!! Da decidere cosa fare. Con annotatr le regioni "intergenic" non vengono associate a SYMBOL, quindi si perdono.
##    Due opzioni:
##    - Togliere "intergenic" e "1to5kb" anche da "All DM genes", cioè non considerarle nel confronto del grafico
##    - Sostituire annotatr con CHIPseeker: così anche intergenic viene considerato
## Da capire come il test statistic è influenzato da questo. Controlla..


#Importing DE genes
DE_results <- read.csv("C:/Users/pierp/Desktop/THESIS PROJECT/1_RNA-Seq/dataset_1/DE_results.csv")
DE_results <- DE_results %>% dplyr::rename(SYMBOL = hugo_symbol)  #renaming for coherence

#Creating list for graphs and table
graph_list <- list()
Stats_region <- list()

row.names(Stats_region) = Stats_region$metric

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
  annot <- build_annotations(genome = "hg38", annotations = c("hg38_genes_promoters", "hg38_genes_5UTRs", "hg38_genes_exons", "hg38_genes_introns",  
                                                              "hg38_genes_3UTRs", "hg38_genes_1to5kb", "hg38_genes_intergenic"))
  
  DM_ann_GR <- annotate_regions(regions = DM_GR,
                                annotations = annot,
                                minoverlap = 1L,
                                ignore.strand = TRUE,
                                quiet = FALSE)
  DM_ann <- as.data.frame(DM_ann_GR)
  DM_ann$coord_key <- paste(DM_ann$seqnames, DM_ann$start, sep="_")
  
  #annotatr finds many associations for each CpG. Collapsing into 1, according to following hierarchy: 
  hyerarchy <- c("hg38_genes_promoters"   = 1,
                 "hg38_genes_5UTRs"       = 2,
                 "hg38_genes_exons"       = 3,
                 "hg38_genes_introns"     = 4,
                 "hg38_genes_3UTRs"       = 5,
                 "hg38_genes_1to5kb"      = 6,
                 "hg38_genes_intergenic"  = 7)
  
  DM_ann_df <- DM_ann %>% mutate(rank = hyerarchy[annot.type]) %>%
                      group_by(coord_key) %>%
                      slice_min(rank, n = 1, with_ties = FALSE) %>%
                      ungroup() # keeps the feature with the highest priority
  
  DM_ann_df <- DM_ann_df %>% dplyr::rename(SYMBOL = "annot.symbol")
  
  
  
  ### Region Analysis ###
  
  #Intersection:
  intersect_genes <- intersect(DE_results$SYMBOL, DM_ann_df$SYMBOL)
  
  ## DMR on promoter/exon/intron/UTR/intergenic genes
  #! exon contains 5' and 3' UTR. "others" contains intergenic and 1to5kb upstream
  prom_genes <- unique(DM_ann_df$SYMBOL[grepl("hg38_genes_promoters", DM_ann_df$annot.type)])
  exon_genes <- unique(DM_ann_df$SYMBOL[grepl("hg38_genes_exons|hg38_genes_5UTRs|hg38_genes_3UTRs", DM_ann_df$annot.type)])
  intron_genes <- unique(DM_ann_df$SYMBOL[grepl("hg38_genes_introns", DM_ann_df$annot.type)])
  #!! Molti SYMBOL si perdono: sono mappati come NA su annotatr. Questo perchè cadono in regioni intergeniche/dubbie/non bene annotate
  
  #Obtaining only exon/intron/UTR intersecting genes
  prom_inter_genes <- intersect(prom_genes, DE_results$SYMBOL)
  exon_inter_genes <- intersect(exon_genes, DE_results$SYMBOL)
  intron_inter_genes <- intersect(intron_genes, DE_results$SYMBOL)
  
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
  # N <- as.numeric(length(unique(DM_ann_df$SYMBOL)))           # all DM genes
  # m <- as.numeric(length(other_genes))                        # intergenic DM genes
  # k <- as.numeric(length(intersect_genes))                    # all intersecting DM genes
  # q <- as.numeric(length(other_inter_genes))                  # intersecting prom on DM genes
  # hyp_other <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  
  perc_prom <- as.numeric(length(prom_inter_genes))*100/as.numeric(length(prom_genes))
  perc_ex = as.numeric(length(exon_inter_genes))*100/as.numeric(length(exon_genes))
  perc_int = as.numeric(length(intron_inter_genes))*100/as.numeric(length(intron_genes))
  #perc_other = as.numeric(length(other_inter_genes))*100/as.numeric(length(other_genes))
  
  Stats_region[[METHOD]] <- c(paste0(round(perc_prom, 1)), hyp_prom_trend, paste0(round(perc_ex, 1)), hyp_ex, 
                                    paste0(round(perc_int, 1)), hyp_int)
  names(Stats_region[[METHOD]]) <- c("Perc ∩ prom", "Sign prom (p)", "Perc ∩ exon", "Sign exon (p)", "Perc ∩ intron", 
                                         "Sign intron (p)")
  
  
  
  ## Graph: Region type ##
  
  # Recreating the region variables, but with coord_key (so not to exclude NA SYMBOL)
  prom_coord <- unique(DM_ann_df$coord_key[grepl("hg38_genes_promoters", DM_ann_df$annot.type)])
  UTR5_coord <- unique(DM_ann_df$coord_key[grepl("hg38_genes_5UTRs", DM_ann_df$annot.type)])
  exon_coord <- unique(DM_ann_df$coord_key[grepl("hg38_genes_exons", DM_ann_df$annot.type)])
  intron_coord <- unique(DM_ann_df$coord_key[grepl("hg38_genes_introns", DM_ann_df$annot.type)])
  UTR3_coord <- unique(DM_ann_df$coord_key[grepl("hg38_genes_3UTRs", DM_ann_df$annot.type)])
  to15kb_coord <- unique(DM_ann_df$coord_key[grepl("hg38_genes_1to5kb", DM_ann_df$annot.type)])
  interg_coord <- unique(DM_ann_df$coord_key[grepl("hg38_genes_intergenic", DM_ann_df$annot.type)])
  
  #Integrating missing before:
  UTR5_genes <- unique(DM_ann_df$SYMBOL[grepl("hg38_genes_5UTRs", DM_ann_df$annot.type)])
  UTR3_genes <- unique(DM_ann_df$SYMBOL[grepl("hg38_genes_3UTRs", DM_ann_df$annot.type)])
  to15kb_genes <- unique(DM_ann_df$SYMBOL[grepl("hg38_genes_1to5kb", DM_ann_df$annot.type)])
  interg_genes <- unique(DM_ann_df$SYMBOL[grepl("hg38_genes_intergenic", DM_ann_df$annot.type)])
  
  UTR5_inter_genes <- intersect(UTR5_genes, DE_results$SYMBOL)
  UTR3_inter_genes <- intersect(UTR3_genes, DE_results$SYMBOL)
  to15kb_inter_genes <- intersect(to15kb_genes, DE_results$SYMBOL)
  interg_inter_genes <- intersect(interg_genes, DE_results$SYMBOL)
  
  # Building counts
  region_summary <- data.frame(
    region = c("Promoter", "5'UTR", "Exon", "Intron", "3'UTR", "1to5kb up", "Intergenic"),
    DM_all = c(length(prom_coord), length(UTR5_coord), length(exon_coord), length(intron_coord), length(UTR3_coord), length(to15kb_coord), length(interg_coord)),
    inters_DM = c(length(prom_inter_genes), length(UTR5_inter_genes), length(exon_inter_genes), length(intron_inter_genes), length(UTR3_inter_genes), length(to15kb_inter_genes), length(interg_inter_genes))
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
           region = factor(region, levels = c("Promoter", "5'UTR", "Exon", "Intron", "3'UTR", "1to5kb up", "Intergenic")))
  
  region_colors <- c(
    "Promoter" = "#E41A1C",
    "5'UTR"    = "darkgreen",
    "Exon"     = "#FF7F00",
    "Intron"   = "#4393C3",
    "3'UTR"    = "lightgreen",
    "1to5kb up"    = "yellow",
    "Intergenic"   = "#984EA3"
  )
  
  # Tot of annotation on bar
  totals <- region_long %>%
    group_by(group) %>%
    summarise(total = sum(n), .groups = "drop")
  
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
      title = sprintf("Method %d: Genomic region distribution", METHOD),
      subtitle = "All DM genes vs DM ∩ DE genes",
      x = NULL,
      y = "Percentage of genes"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 12)
    )
}

#Graph Output
pdf("C:/Users/pierp/Desktop/THESIS PROJECT/4_Integration_results/dataset_1/Region_analysis.pdf", height = 10, width = 15)
plot_grid(plotlist = graph_list, ncol = 2)
dev.off()

#Stats output
Stats_df <- as.data.frame(do.call(cbind, Stats_region))
colnames(Stats_df) <- sprintf("M%d", seq_along(Stats_region))
write.csv(Stats_df, "C:/Users/pierp/Desktop/THESIS PROJECT/4_Integration_results/dataset_1/Stats_table_region.csv", row.names = FALSE)