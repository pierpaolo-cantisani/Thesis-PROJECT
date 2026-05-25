library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

### Da cambiare leggermente: Prendere per ogni metodo le coordinate genomiche delle associazioni (considerando intersezione DE DM!!), 
# invece dei SYMBOLS (considerando ogni coordinata 1 sola volta! Attenzione nei M4 e M5). Su queste
# rifare l'annotazione (magari CHIPseeker, o annotatr meglio) sovrascrivendola, così la hanno tutta uguale.
# Il risultato da ottenere alla fine è: quale regione genomica (metilata) tende ad avere più/meno effetto sull'espressione?

#Explaination of Stats_region:
# Perc ∩ prom: % of sites for that region that are intersected, compared to total sites of that region
# Sign prom (p): Are the genes with DM sites in this region more differentially expressed compared to other regions?
#               i.e.: Are intersected genes specifically associated to exons/introns/promoters/3UTR? ###


#!! The 2 methods have different "promoter" definition:
#        - Method 1: (-1000, 1000)b from TSS
#        - Method 2: (-2000, 200)b from TSS
# Difference in "promoter" DM genes might be due to this


### CODE START ###
### Importing DE and DM data ###

#Importing RNA-Seq DE genes
DE_results <- read.csv("C:/Users/pierp/Desktop/THESIS PROJECT/1_RNA-Seq/dataset_1/DE_results.csv")
DE_results <- DE_results %>% dplyr::rename(SYMBOL = hugo_symbol)  #renaming for coherence

#Considerations will be made for METHOD 1 (genomation) and METHOD 2 (CHIPseeker)

##### METHOD 1 (genomation) #####
#Importing DM sites
DM_sites_M1 <- read.csv("C:/Users/pierp/Desktop/THESIS PROJECT/3_Benchmark/dataset_1/DM_sites_Met1.csv")

#Intersection:
intersect_genes <- intersect(DE_results$SYMBOL, DM_sites_M1$SYMBOL)

## DMR on promoter/exon/intron/UTR/intergenic genes
prom_genes <- unique(DM_sites_M1$SYMBOL[grepl("promoter", DM_sites_M1$region_type)])
exon_genes <- unique(DM_sites_M1$SYMBOL[grepl("exon", DM_sites_M1$region_type)])
intron_genes <- unique(DM_sites_M1$SYMBOL[grepl("intron", DM_sites_M1$region_type)])
other_genes <- unique(DM_sites_M1$SYMBOL[grepl("intergenic", DM_sites_M1$region_type)])

#Obtaining only exon/intron/UTR intersecting genes
prom_inter_genes <- intersect(prom_genes, DE_results$SYMBOL)
exon_inter_genes <- intersect(exon_genes, DE_results$SYMBOL)
intron_inter_genes <- intersect(intron_genes, DE_results$SYMBOL)
other_inter_genes <- intersect(other_genes, DE_results$SYMBOL)

## Some of the cases have very low numbers: must be careful with the test interpretation. 

## Hypergeometric tests
#Are intersecting genes associated to methylation specific to introns/exons/3UTR/promoters?

#Promoter
N <- as.numeric(length(unique(DM_sites_M1$SYMBOL)))           # all DM genes
m <- as.numeric(length(prom_genes))                           # prom DM genes
k <- as.numeric(length(intersect_genes))                      # all intersecting DM genes
q <- as.numeric(length(prom_inter_genes))                     # intersecting prom on DM genes
hyp_prom_trend <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#exon
N <- as.numeric(length(unique(DM_sites_M1$SYMBOL)))           # all DM genes
m <- as.numeric(length(exon_genes))                           # exon DM genes
k <- as.numeric(length(intersect_genes))                      # all intersecting DM genes
q <- as.numeric(length(exon_inter_genes))                     # intersecting exon DM genes
hyp_ex <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#intron
N <- as.numeric(length(unique(DM_sites_M1$SYMBOL)))           # all DM genes
m <- as.numeric(length(intron_genes))                         # intron DM genes
k <- as.numeric(length(intersect_genes))                      # all intersecting DM genes
q <- as.numeric(length(intron_inter_genes))                   # intersecting intron DM genes
hyp_int <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#Intergenic
N <- as.numeric(length(unique(DM_sites_M1$SYMBOL)))           # all DM genes
m <- as.numeric(length(other_genes))                          # intergenic DM genes
k <- as.numeric(length(intersect_genes))                      # all intersecting DM genes
q <- as.numeric(length(other_inter_genes))                    # intersecting prom on DM genes
hyp_other <- phyper(q-1, m, N-m, k, lower.tail = FALSE)


perc_prom <- as.numeric(length(prom_inter_genes))*100/as.numeric(length(prom_genes))
perc_ex = as.numeric(length(exon_inter_genes))*100/as.numeric(length(exon_genes))
perc_int = as.numeric(length(intron_inter_genes))*100/as.numeric(length(intron_genes))
perc_other = as.numeric(length(other_inter_genes))*100/as.numeric(length(other_genes))


Stats_region <- data.frame(M1 = c(paste0(round(perc_prom, 1)), hyp_prom_trend, paste0(round(perc_ex, 1)), hyp_ex, 
                                  paste0(round(perc_int, 1)), hyp_int, paste0(round(perc_other, 1)), hyp_other))
row.names(Stats_region) = c("Perc ∩ prom", "Sign prom (p)", "Perc ∩ exon", "Sign exon (p)", "Perc ∩ intron", 
                            "Sign intron (p)", "Perc ∩ intergen", "Sign intergen (p)")



## Graph: Region type ##
# Building counts
region_summary <- data.frame(
  region = c("Promoter", "Exon", "Intron", "Intergenic"),
  DM_all = c(length(prom_genes), length(exon_genes), length(intron_genes), length(other_genes)),
  inters_DM = c(length(prom_inter_genes), length(exon_inter_genes), length(intron_inter_genes), length(other_inter_genes))
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
         region = factor(region, levels = c("Promoter", "Exon", "Intron", "Intergenic")))

region_colors <- c(
  "Promoter" = "#E41A1C",
  "Exon"     = "#FF7F00",
  "Intron"   = "#4393C3",
  "Intergenic"   = "#984EA3"
)

# Tot of annotation on bar
totals <- region_long %>%
  group_by(group) %>%
  summarise(total = sum(n), .groups = "drop")

p_region_1 <- ggplot(region_long, aes(x = group, y = pct, fill = region)) +
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
    title = "Method 1: Genomic region distribution",
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


#####
##### METHOD 2 (CHIPseeker): #####

#Importing DM sites
DM_sites_M2 <- read.csv("C:/Users/pierp/Desktop/THESIS PROJECT/3_Benchmark/dataset_1/DM_sites_Met2.csv")

#Intersection
intersect_genes_2 <- intersect(DE_results$SYMBOL, DM_sites_M2$SYMBOL)

## DMR on promoter/exon/intron/UTR/intergenic genes
prom_genes <- unique(DM_sites_M2$SYMBOL[grepl("^Promoter", DM_sites_M2$region_type)])
exon_genes <- unique(DM_sites_M2$SYMBOL[grepl("UTR|^Exon", DM_sites_M2$region_type)])
intron_genes <- unique(DM_sites_M2$SYMBOL[grepl("^Intron", DM_sites_M2$region_type)])
other_genes <- unique(DM_sites_M2$SYMBOL[grepl("^Distal Intergenic", DM_sites_M2$region_type)])


#Obtaining only exon/intron/UTR intersecting genes
prom_inter_genes <- intersect(prom_genes, DE_results$SYMBOL)
exon_inter_genes <- intersect(exon_genes, DE_results$SYMBOL)
intron_inter_genes <- intersect(intron_genes, DE_results$SYMBOL)
other_inter_genes <- intersect(other_genes, DE_results$SYMBOL)

## Some of the cases have very low numbers: must be careful with the test interpretation. 

## Hypergeometric tests
#Are intersecting genes associated to methylation specific to introns/exons/3UTR/promoters?

#Promoter
N <- as.numeric(length(unique(DM_sites_M2$SYMBOL)))           # all DM genes
m <- as.numeric(length(prom_genes))                           # prom DM genes
k <- as.numeric(length(intersect_genes_2))                      # all intersecting DM genes
q <- as.numeric(length(prom_inter_genes))                     # intersecting prom DM genes
hyp_prom_trend <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#exon
N <- as.numeric(length(unique(DM_sites_M2$SYMBOL)))           # all DM genes
m <- as.numeric(length(exon_genes))                           # exon DM genes
k <- as.numeric(length(intersect_genes_2))                      # all intersecting DM genes
q <- as.numeric(length(exon_inter_genes))                     # intersecting exon DM genes
hyp_ex <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#intron
N <- as.numeric(length(unique(DM_sites_M2$SYMBOL)))           # all DM genes
m <- as.numeric(length(intron_genes))                         # intron DM genes
k <- as.numeric(length(intersect_genes_2))                      # all intersecting DM genes
q <- as.numeric(length(intron_inter_genes))                   # intersecting intron DM genes
hyp_int <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#intergenic
N <- as.numeric(length(unique(DM_sites_M2$SYMBOL)))           # all DM genes
m <- as.numeric(length(other_genes))                          # intergenic DM genes
k <- as.numeric(length(intersect_genes_2))                      # all intersecting DM genes
q <- as.numeric(length(other_inter_genes))                    # intersecting prom on DM genes
hyp_other <- phyper(q-1, m, N-m, k, lower.tail = FALSE)


perc_prom <- as.numeric(length(prom_inter_genes))*100/as.numeric(length(prom_genes))
perc_ex = as.numeric(length(exon_inter_genes))*100/as.numeric(length(exon_genes))
perc_int = as.numeric(length(intron_inter_genes))*100/as.numeric(length(intron_genes))
perc_other = as.numeric(length(other_inter_genes))*100/as.numeric(length(other_genes))


Stats_region$M2 <- c(paste0(round(perc_prom, 1)), hyp_prom_trend, paste0(round(perc_ex, 1)), hyp_ex, 
                                  paste0(round(perc_int, 1)), hyp_int, paste0(round(perc_other, 1)), hyp_other)




## Graph: Region type ##
# Building counts
region_summary <- data.frame(
  region = c("Promoter", "Exon", "Intron", "Intergenic"),
  DM_all = c(length(prom_genes), length(exon_genes), length(intron_genes), length(other_genes)),
  inters_DM = c(length(prom_inter_genes), length(exon_inter_genes), length(intron_inter_genes), length(other_inter_genes))
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
         region = factor(region, levels = c("Promoter", "Exon", "Intron", "Intergenic")))

region_colors <- c(
  "Promoter" = "#E41A1C",
  "Exon"     = "#FF7F00",
  "Intron"   = "#4393C3",
  "Intergenic"   = "#984EA3"
)

# Tot of annotation on bar
totals <- region_long %>%
  group_by(group) %>%
  summarise(total = sum(n), .groups = "drop")

p_region_2 <- ggplot(region_long, aes(x = group, y = pct, fill = region)) +
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
    title = "Method 2: Genomic region distribution",
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
 
#####
##### Final outputs: Graphs and Table #####

##Plotting both graphs:
pdf("C:/Users/pierp/Desktop/THESIS PROJECT/4_Integration_results/dataset_1/Region_analysis.pdf")
plot_grid(p_region_1, p_region_2, ncol = 2)
dev.off()


#Adding metric names and writing output:
Stats_region$metric <- rownames(Stats_region)
Stats_region <- Stats_region[ , c("metric", setdiff(names(Stats_region), "metric"))]  # Puts "metric" as first column

#Writing on file
write.csv(Stats_region, "C:/Users/pierp/Desktop/THESIS PROJECT/4_Integration_results/dataset_1/Stats_table_region.csv", row.names = FALSE)
