library(dplyr)
library(ggplot2)
library(tidyr)
library(openxlsx)
library(cowplot)

#In sezione 2.3 (hypo/hyper vs up/down) costruisci N, m, k, q contando righe (siti), non geni, e lo dichiari esplicitamente 
#("analysis is for site, not for gene"). Concettualmente legittimo, ma l'ipergeometrica assume estrazioni senza rimpiazzo da 
#una popolazione di unità indipendenti. I siti sullo stesso gene non sono indipendenti (condividono lo stesso log2FC: il gene 
#uno solo!). Quindi un gene con 5 siti DM hypo e fortemente up conta 5 volte nello stesso "successo", gonfiando artificialmente 
#q e N. Questo viola l'indipendenza alla base del test e rende i p-value anti-conservativi (troppo significativi). 
#Per M4/M5, dove la molteplicità è alta, l'effetto è marcato. Non è un bug di codice, è un problema di validità del test: 
#in tesi è esattamente il punto su cui un biostatistico ti incalzerebbe. La via difendibile è fare questi test a livello di 
#(un gene = un'osservazione, classificato per il suo trend prevalente o per il sito a maggior |meth.diff|), oppure dichiarare apertamente il limite.
#Secondo punto statistico: stai facendo ~15 test ipergeometrici per metodo × 5 metodi ≈ 75 test, senza nessuna correzione per multiple testing su questa famiglia. 
#Le sezioni 2-4 producono p-value che poi interpreti come "YES se < 0.05". Con 75 test, diversi saranno < 0.05 per caso. 
#Anche solo un BH sulla famiglia di test concettualmente omogenei renderebbe le conclusioni molto più solide.
#Terzo: in sezione 2.1 hyp_DM_up e poi hyp_fin_up usano DM_up_genes calcolato come all_up$SYMBOL[all_up$SYMBOL %in% DM_genes] — ma DM_genes è unique(DM_sites$SYMBOL) 
#(tutti i DM, non intersecato con l'universo), mentre N usa length(universe). C'è un mismatch tra la popolazione (universe) e
#il sottoinsieme (DM_genes non ristretto a universe): alcuni DM_genes potrebbero non essere nell'universo, quindi k e q sono 
#contati su una base leggermente diversa da N. Andrebbe usato DM_genes_univ (che hai già!) ovunque per coerenza. Questo è un bug sottile di consistenza dei conteggi che falsa leggermente i p-value.



### !!! Summary:
#This analysis answer the following questions (p-value < 0.05 means YES):

#For the METHOD == X:
## Section 1
# (Significance of intersection): Are DM genes more Differentially Expressed compared to all genes?
## Section 2
# (i)All DM (p-value): Do DM genes have a trend towards up/down regulation?
# (ii)DM ∩ DE DM (p-value): Do DM ∩ DE genes have a trend towards up/down regulation?
# (iii)DM ∩ DE vs All DM (p-value): Do DM ∩ DE have a significantly higher trend towards up/down regulation compared to All DM?
##!!! Queste tre domande insieme vedono se c'è un trend di tipo hypo/up, hyper/down,e se questo è specifico per DM ∩ DE, o una regola generale dei DM. 
# Hypomethylated DM: Do hypomethylated DM genes have a trend towards up/down regulation?
# Hypermethylated DM: Do hypermethylated DM genes have a trend towards up/down regulation?
## Section 3
# (meth vs expr): Is the magnitude of the methylation (in the category) correlated to the genes up/down regulation?
## Section 4
# (multi DM genes) : Are genes with 2 or more (non intergenic) methylation sites more likely to be DE than those with 1?
## Section 5
# (sites vs log2FC): Do the number of DM sites in a gene correlate with its expression (making it more likely DE)?
# (sites vs |log2FC|): Do the number of DM sites in a DE gene correlate with the magnitude of its differential expression?

# All p-values within a method are BH-adjusted as a single family of 15 tests.
# This includes 12 hypergeometric tests + 3 Spearman correlation tests.
# Significance threshold: padj < 0.05.

PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"

### 0. Importing DE genes ###
#Importing RNA-Seq DE genes
DE_results <- read.csv(file.path(PATH, "Dataset_2", "1_RNA-Seq", "DE_results.csv"))
DE_results <- DE_results %>% dplyr::rename(SYMBOL = hugo_symbol)  #renaming for coherence

## Obtaining the universe N: 
#Importing RNA-seq universe (all genes considered for the DESeq2 analysis)
RNAseq_universe <- read.csv(file.path(PATH, "Dataset_2", "1_RNA-Seq", "RNAseq_universe.csv"), col.names = c("SYMBOL", "log2FoldChange", "padj"))

## Choosing "universe" as the RNAseq_universe. So that the METHOD's association can theorically connect to any of those genes.
universe <- RNAseq_universe$SYMBOL

# Final genes ready for comparison are:
DE_genes <- unique(DE_results$SYMBOL)

#And DM sites
Scalar_M <- readRDS(file.path(PATH, "Dataset_2", "4_Integration_results", "Method_final_df.rds"))


#Output data for each method will be temporarily inserted in this list. Same for Graph outputs.
Data_list <- list()
magnitude_list <- list()     #For section 3
num_barplot_list <- list()   #For section 4
scatter_list <- list()       #For section 5

### Iterating on METHOD:
for(METHOD in 1:5) {
  
  ##### CODE START: ##### 
  
  #Importing DM METHOD lists:
  DM_sites <- read.csv(sprintf(file.path(PATH, "Dataset_2", "3_Benchmark", "DM_sites_Met%d.csv"), METHOD))
  DM_genes <- unique(DM_sites$SYMBOL)
  
  DM_genes_univ <- unique(intersect(DM_sites$SYMBOL, universe))
  
  
  ### 1. Analysis of the intersection between DE genes and DM sites ###
  
  intersect_genes <- unique(DM_sites$SYMBOL[DM_sites$SYMBOL %in% DE_results$SYMBOL])
  intersect_df <- merge(DM_sites[, c("coord_key", "meth.diff", "SYMBOL")], DE_results[, c("log2FoldChange", "padj", "SYMBOL")], by = "SYMBOL")
  
  ## Hypergeometric test
  #Are DE genes more Differentially methylated compared to how differentially methylated all genes are?
  #!! THis is not so informative, since it's obvious that DM genes will have association (but NOT CAUSATION) to DE genes
  N <- as.numeric(length(universe))                         # universe: all genes from RNA-Seq
  m <- as.numeric(length(DM_genes_univ))                    # DM genes (intersected with universe)
  k <- as.numeric(length(DE_genes))                         # DE genes                        
  q <- as.numeric(length(intersect_genes))                  # intersecting genes
  
  hyp_intersect <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  
  
  ### Exploring trends ###
  
  ### 2: Do intersected genes have a trend towards up/down regulation? ###
  
  #2.1: UPREGULATION
  DE_up <- DE_results %>% filter(log2FoldChange > 0)
  
  #Obtaining up and down intersecting genes
  intersect_up <- intersect_df %>% filter(log2FoldChange > 0)
  intersect_up_genes <- unique(intersect_up$SYMBOL)     # upregolated intersecting genes
  intersect_down <- intersect_df %>% filter(log2FoldChange < 0)
  intersect_down_genes <- unique(intersect_down$SYMBOL) # downregolated intersecting genes
  
  # universe restricted to DE_genes only, because we are testing 
  # directionality of expression within the DE subset, not across all genes
  N <- as.numeric(length(DE_genes))                         # universe: all DE_genes
  m <- as.numeric(length(DE_up$SYMBOL))                     # all DE up genes
  k <- as.numeric(length(intersect_genes))                  # all intersected DM genes
  q <- as.numeric(length(intersect_up_genes))               # upregulated intersected DM genes
  
  hyp_up <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  
  ## But what if all DM genes are significantly enriched for upregulation? If that is the case this result is less interesting. 
  
  ##Let's see:
  #Obtaining all up genes (also non DE)
  all_up <- RNAseq_universe %>% filter(log2FoldChange > 0)
  DM_up_genes <- all_up$SYMBOL[all_up$SYMBOL %in% DM_genes]
  
  ## Hypergeometric test on all rel DM sites
  #How significant is the trend of all DM genes towards being upregulated
  N <- as.numeric(length(universe))                         # universe: all RNA-Seq genes
  m <- as.numeric(length(all_up$SYMBOL))                    # all upregulated genes (also non DE)
  k <- as.numeric(length(DM_genes_univ))                    # all DM genes (inters univ)
  q <- as.numeric(length(DM_up_genes))                      # upregulated DM genes
  
  hyp_DM_up <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  
  ## The next question is: is the DM ∩ DE case significantly more upregulated than all DM?
  
  ## Hypergeometric test on all rel DM sites
  #How significant is the trend of all DM genes towards being upregulated
  N <- as.numeric(length(DM_genes_univ))                    # all DM genes (inters univ)
  m <- as.numeric(length(unique(DM_up_genes)))              # all DM upregulated genes (also non DE)
  k <- as.numeric(length(intersect_genes))                  # all intersected DM genes
  q <- as.numeric(length(intersect_up_genes))               # upregulated intersected DM genes
  
  hyp_fin_up <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  
  
  ##2.2: DOWNREGULATION:
  
  DE_down <- DE_results %>% filter(log2FoldChange < 0)
  
  # All down genes in universe (also non DE)
  all_down <- RNAseq_universe %>% filter(log2FoldChange < 0)
  DM_down_genes <- all_down$SYMBOL[all_down$SYMBOL %in% DM_genes]
  
  
  # Hypergeometric test: are all DM genes enriched for downregulation?
  N <- as.numeric(length(universe))                 # universe
  m <- as.numeric(length(all_down$SYMBOL))          # all downregulated genes in universe
  k <- as.numeric(length(DM_genes_univ))            # all DM genes (inters univ)
  q <- as.numeric(length(DM_down_genes))            # downregulated DM genes
  hyp_DM_down <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  # Hypergeometric test: are DM ∩ DE genes enriched for downregulation?
  N <- as.numeric(length(DE_genes))                 # universe: all DE genes
  m <- as.numeric(length(DE_down$SYMBOL))           # all DE down genes
  k <- as.numeric(length(intersect_genes))          # all intersected DM genes
  q <- as.numeric(length(intersect_down_genes))     # downregulated intersected DM genes
  hyp_down <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  # Hypergeometric test: is DM ∩ DE more downregulated than all relevant DM?
  N <- as.numeric(length(DM_genes_univ))            # all DM genes (inters univ)
  m <- as.numeric(length(DM_down_genes))            # all DM downregulated genes (also non DE)
  k <- as.numeric(length(intersect_genes))          # all intersecting genes
  q <- as.numeric(length(intersect_down_genes))     # all downregulated intersected genes
  hyp_fin_down <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  
  
  ##2.3: Association: UP/DOWN DE vs Hyper/Hypo DM
  
  #Obtaining hypo/hyper lists
  hypo_DM <- DM_sites %>% filter(meth.diff < 0)
  hyper_DM <- DM_sites %>% filter(meth.diff > 0)
  hypo_int_DM <- hypo_DM %>% filter(SYMBOL %in% DE_results$SYMBOL)
  hyper_int_DM <- hyper_DM %>% filter(SYMBOL %in% DE_results$SYMBOL)
  #I'm not filtering for unique genes. So there will be duplicates: this analysis will be done for site, rather than gene
  hypo_int_up_DM <- hypo_DM %>% filter(SYMBOL %in% DE_up$SYMBOL)
  hyper_int_up_DM <- hyper_DM %>% filter(SYMBOL %in% DE_up$SYMBOL)
  hypo_int_down_DM <- hypo_DM %>% filter(SYMBOL %in% DE_down$SYMBOL)
  hyper_int_down_DM <- hyper_DM %>% filter(SYMBOL %in% DE_down$SYMBOL)
  
  # Hypergeometric test: Are hypo/hyper DM sites associated with strong up/downregulation?
  # hypo-up
  N <- as.numeric(length(hypo_int_DM$SYMBOL) + length(hyper_int_DM$SYMBOL))            # All intersecting genes
  m <- as.numeric(length(hypo_int_up_DM$SYMBOL) + length(hyper_int_up_DM$SYMBOL))      # All up intersecting genes
  k <- as.numeric(length(hypo_int_DM$SYMBOL))                                          # All hypo intersecting genes
  q <- as.numeric(length(hypo_int_up_DM$SYMBOL))                                       # hypo-up intersecting genes
  hyp_up_hypo <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  # hypo-down
  N <- as.numeric(length(hypo_int_DM$SYMBOL) + length(hyper_int_DM$SYMBOL))            # All intersecting genes
  m <- as.numeric(length(hypo_int_down_DM$SYMBOL) + length(hyper_int_down_DM$SYMBOL))  # All down intersecting genes
  k <- as.numeric(length(hypo_int_DM$SYMBOL))                                          # All hypo intersecting genes
  q <- as.numeric(length(hypo_int_down_DM$SYMBOL))                                     # hypo-down intersecting genes
  hyp_down_hypo <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  # hyper-up
  N <- as.numeric(length(hypo_int_DM$SYMBOL) + length(hyper_int_DM$SYMBOL))            # All intersecting genes
  m <- as.numeric(length(hypo_int_up_DM$SYMBOL) + length(hyper_int_up_DM$SYMBOL))      # All up intersecting genes
  k <- as.numeric(length(hyper_int_DM$SYMBOL))                                         # All hyper intersecting genes
  q <- as.numeric(length(hyper_int_up_DM$SYMBOL))                                      # hyper-up intersecting genes
  hyp_up_hyper <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  # hyper-down
  N <- as.numeric(length(hypo_int_DM$SYMBOL) + length(hyper_int_DM$SYMBOL))            # All intersecting genes
  m <- as.numeric(length(hypo_int_down_DM$SYMBOL) + length(hyper_int_down_DM$SYMBOL))  # All down intersecting genes
  k <- as.numeric(length(hyper_int_DM$SYMBOL))                                         # All hypo intersecting genes
  q <- as.numeric(length(hyper_int_down_DM$SYMBOL))                                    # hypo-down intersecting genes
  hyp_down_hyper <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  
  ##! Duplicates were not cleaned: this is wanted, as different sites on the same gene may have different trends.
  ##  In this way all information is kept. Therefore the analysis is for "site", not for "gene".
  
  
  
  
  ### 3: The more the sites are differentially methylated (in magnitude), the more the gene is up/down regulated? (Correlation between magnitude of methyl and up/down regulation) ###
  
  #Considering unique genes: meth will be the mean of the sites
  #Using all of the DM genes (in the universe), not the intersection (the intersection will be then highlighted in the graph)
  DM_univ_df <-  merge(DM_sites[, c("coord_key", "meth.diff", "SYMBOL")], RNAseq_universe[, c("log2FoldChange", "padj", "SYMBOL")], by = "SYMBOL")
  
  #Importing M-values
  gene_level_df <- Scalar_M[[METHOD]] %>%
    group_by(SYMBOL) %>% summarise(
      mean_M = mean(Mv),
      log2FC = unique(log2FC))
  
  #Intersected
  #all
  all_cor <- cor.test(gene_level_df$mean_M, gene_level_df$log2FC, method = "spearman")
  
  ## Graph: Scatter plot ##
  magnitude_list[[METHOD]] <- ggplot(gene_level_df, aes(x = mean_M, y = log2FC)) +
    geom_point() +
    geom_point(data = gene_level_df[gene_level_df$SYMBOL %in% intersect_df$SYMBOL, ], color = "red") +
    geom_smooth(method = "lm") +
    labs(title = sprintf("M%d: Body methylation vs Expression", METHOD),
         x        = "Mean meth.diff per gene",
         y        = "log2FC")
  
  
  
  
  
  ### 4: Effects of number of DM sites on Expression ###
  
  ##4.1: The more DM sites a gene has, the more probability of being DE it has?
  
  #Extracting the genes with more than 1 DM site
  DM_num <- as.data.frame(table(DM_sites$SYMBOL))
  colnames(DM_num) <- c("SYMBOL", "n_DM_sites")
  DM_num$SYMBOL <- as.character(DM_num$SYMBOL)       #Avoid it being "factor"
  
  multi_DM_genes  <- DM_num$SYMBOL[DM_num$n_DM_sites >= 2]   # ≥ 2 sites
  #Intersecting:
  multi_intersect_genes <- intersect(multi_DM_genes, intersect_genes)
  
  ## Hypergeometric test: Are genes with 2 or more (non intergenic) methyl sites more likely to be DE than those with 1?
  N <- as.numeric(length(DM_genes))                           # all DM genes
  m <- as.numeric(length(intersect_genes))                    # all intersecting DM DE genes
  k <- as.numeric(length(multi_DM_genes))                     # genes with 2 or more DM sites
  q <- as.numeric(length(multi_intersect_genes))              # intersecting DM DE genes with 2 or more DM sites
  hyp_multi <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
  
  
  ## Graph: Probability of being DE depending on number of DM sites ##
  # Summary dataframe
  
  single_DM <- DM_sites %>% filter(SYMBOL %in% DM_num$SYMBOL[DM_num$n_DM_sites < 2])
  single_DM_genes <- DM_num$SYMBOL[DM_num$n_DM_sites < 2]    # = 1 site

  single_intersect_genes <- setdiff(intersect_genes, multi_intersect_genes)
  dm_summary <- data.frame(
    group    = c("All DM genes", "DM%DE genes"),
    single   = c(length(single_DM_genes),  length(single_intersect_genes)),
    multi    = c(length(multi_DM_genes),   length(multi_intersect_genes))
  )
  
  # Long format
  dm_long <- dm_summary %>%
    pivot_longer(cols = c(single, multi),
                 names_to  = "dm_class",
                 values_to = "n") %>%
    group_by(group) %>%
    mutate(
      total    = sum(n),
      pct      = n / total * 100,
      group    = factor(group, levels = c("All DM genes", "DM%DE genes")),
      dm_class = factor(dm_class,
                        levels = c("single", "multi"),
                        labels = c("1 DM site", "≥2 DM sites"))
    )
  
  # Totals for bar annotation
  totals_dm <- dm_long %>%
    distinct(group, total)
  
  # Colors
  dm_colors <- c(
    "1 DM site"   = "#4393C3",
    "≥2 DM sites" = "#D6604D"
  )
  
  # Plot
  num_barplot_list[[METHOD]] <- ggplot(dm_long, aes(x = group, y = pct, fill = dm_class)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_text(
      data = totals_dm,
      aes(x = group, y = 102, label = paste0("n = ", total)),
      inherit.aes = FALSE,
      size = 4, fontface = "bold"
    ) +
    scale_fill_manual(values = dm_colors, name = "DM sites per gene") +
    scale_y_continuous(
      limits = c(0, 108),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    labs(
      title    = sprintf("M%d: DM site multiplicity", METHOD),
      subtitle = "All DM genes vs DM%DE genes",
      x        = NULL,
      y        = "Percentage of genes"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position  = "right",
      plot.title       = element_text(face = "bold"),
      axis.text.x      = element_text(size = 12)
    )
  
  
  
  
  
  ### 5: Correlation between number of DM sites and log2FC ###
  
  ## Qui si considerano tutti i geni, non l'intersezione. Da capire se mi convince.
  
  # Counting DM sites per gene
  #Using DM_num, called at the beginning of section 4
  
  # Joining expression values from the RNA-seq universe
  sites_fc_df <- DM_num %>%
    left_join(RNAseq_universe, by = "SYMBOL") %>%
    filter(!is.na(log2FoldChange))
  
  # Spearman correlation: number of DM sites vs magnitude of regulation
  sites_abs_fc_cor <- cor.test(sites_fc_df$n_DM_sites, abs(sites_fc_df$log2FoldChange), method = "spearman")
  
  # Spearman correlation: number of DM sites vs signed log2FC
  sites_fc_cor <- cor.test(sites_fc_df$n_DM_sites, sites_fc_df$log2FoldChange, method = "spearman")
  
  
  
  ## Graph: Scatter plot of number of DM sites vs log2FC ##
  jit <- position_jitter(width = 0.12, height = 0, seed = 42)
  #Magnitude
  p_sites_fc_abs <- ggplot(sites_fc_df, aes(x = n_DM_sites, y = abs(log2FoldChange))) +
    geom_point(alpha = 0.7, position = jit) +
    geom_point(data     = sites_fc_df[sites_fc_df$SYMBOL %in% intersect_df$SYMBOL, ],
               color    = "red", alpha = 0.8,
               position = jit) +
    geom_smooth(method = "lm") +
    labs(title = sprintf("M%d: Correlation # DM sites vs expr (magnitude)", METHOD),
         x = "Number of DM sites",
         y = "log2FoldChange") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  #Direction
  p_sites_fc <- ggplot(sites_fc_df, aes(x = n_DM_sites, y = log2FoldChange)) +
    geom_point(alpha = 0.7, position = jit) +
    geom_point(data     = sites_fc_df[sites_fc_df$SYMBOL %in% intersect_df$SYMBOL, ],
               color    = "red", alpha = 0.8,
               position = jit) +
    geom_smooth(method = "lm") +
    labs(title = sprintf("M%d: Correlation # DM sites vs expr (direction)", METHOD),
         x = "Number of DM sites",
         y = "log2FoldChange") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  scatter_list[[METHOD]] <- plot_grid(p_sites_fc_abs, p_sites_fc, nrow = 2, ncol = 1)
  
  
  
  #Multiple test correction: before output
  pvals_raw <- c(
    hyp_intersect    = hyp_intersect,
    hyp_up           = hyp_up,
    hyp_DM_up        = hyp_DM_up,
    hyp_fin_up       = hyp_fin_up,
    hyp_down         = hyp_down,
    hyp_DM_down      = hyp_DM_down,
    hyp_fin_down     = hyp_fin_down,
    hyp_up_hypo      = hyp_up_hypo,
    hyp_down_hypo    = hyp_down_hypo,
    hyp_up_hyper     = hyp_up_hyper,
    hyp_down_hyper   = hyp_down_hyper,
    hyp_multi        = hyp_multi,
    all_cor_p        = all_cor$p.value,
    sites_abs_fc_p   = sites_abs_fc_cor$p.value,
    sites_fc_p       = sites_fc_cor$p.value
  )
  
  # Applying BH
  pvals_adj <- p.adjust(pvals_raw, method = "BH")
  
  
  ##Final output data:
  Data_list[[sprintf("M%d", METHOD)]] <- c(
    pvals_adj["hyp_intersect"], pvals_adj["hyp_DM_up"],
    pvals_adj["hyp_up"],  pvals_adj["hyp_fin_up"],
    pvals_adj["hyp_DM_down"], pvals_adj["hyp_down"], pvals_adj["hyp_fin_down"],
    pvals_adj["hyp_up_hypo"], pvals_adj["hyp_down_hypo"],
    pvals_adj["hyp_up_hyper"], pvals_adj["hyp_down_hyper"],
    all_cor$estimate, pvals_adj["all_cor_p"],
    length(multi_DM_genes), pvals_adj["hyp_multi"],
    sites_abs_fc_cor$estimate, pvals_adj["sites_abs_fc_p"],
    sites_fc_cor$estimate, pvals_adj["sites_fc_p"]
  )
} 

#Before output, to improve table visualization:
options(scipen = 999)

#Creating the Stats table
final_df <- as.data.frame(Data_list)
row.names(final_df) <- c("Sign of inters (padj)", 
                         "All DM up (padj)", "DM ∩ DE DM up (padj)", "DM ∩ DE vs All DM up (padj)",
                         "All DM down (padj)", "DM ∩ DE DM down (padj)", "DM ∩ DE vs All DM down (padj)",
                         "Hypo-up (padj)", "Hypo-down (padj)", "Hyper-up (padj)", "Hyper-down (padj)",
                         "meth Spear (rho)", "meth Spear (padj)", 
                         "Multi DM genes", "Multi DM ∩ DE (padj)",
                         "sites vs |log2FC| (rho)", "sites vs |log2FC| (padj)",
                         "sites vs log2FC (rho)", "sites vs log2FC (padj)")

#Importing old Stats table:
Stats_table <- read.csv(file.path(PATH, "Dataset_2", "4_Integration_results", "Stats_table.csv"))
rownames(Stats_table) <- Stats_table$metric
Stats_table$metric <- NULL

#Sanity check:
stopifnot(setequal(colnames(Stats_table), colnames(final_df)))

#Merging: output
Stats_table_final <- bind_rows(Stats_table, final_df)
Stats_table_final$metric <- rownames(Stats_table_final)
Stats_table_final <- Stats_table_final[ , c("metric", setdiff(names(Stats_table_final), "metric"))]  # Puts "metric" as first column


#Writing on file
write.xlsx(Stats_table_final, file.path(PATH, "Dataset_2", "4_Integration_results", "Stats_table_final.xlsx"), rowNames = FALSE)


##Plotting graphs:
pdf(file.path(PATH, "Dataset_2", "4_Integration_results", "Additional_comparison D2.pdf"), height = 10, width = 15)

plot_grid(plotlist = magnitude_list, nrow = 2, ncol= 3)
plot_grid(plotlist = num_barplot_list, nrow = 2, ncol= 3)
grid_1 <- plot_grid(plotlist = scatter_list[1:2], nrow = 1, ncol= 2)    #Dividing for graph clarity
print(grid_1)
grid_2 <- plot_grid(plotlist = scatter_list[3:4], nrow = 1, ncol= 2)
print(grid_2)
grid_3 <- plot_grid(plotlist = scatter_list[5], nrow = 1, ncol = 2)
print(grid_3)

dev.off()
