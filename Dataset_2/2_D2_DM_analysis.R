library(data.table)
library(methylKit)
library(matrixStats)
library(ggplot2)
library(rtracklayer)
library(karyoploteR)


PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"

pdf(file.path(PATH, "Dataset_2", "2_BS-Seq", "WGBS Graphs D2.pdf"), width = 12, height = 8)


### 1. File downloading, and MethylKit object building ###

#Names of the import files
files <- c(
  #Olig2 (17)
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877183_1524_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877231_1525_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877219_1527_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877216_1532_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877220_1536_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877217_1539_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877222_1541_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877177_3545_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877184_3586_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877229_3590_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877235_3602_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877212_4615_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877238_AN03398_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877214_AN05483_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877218_AN10090_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877239_AN15240_Olig2_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "olig2", "GSM2877240_AN16799_Olig2_CpG_WGBS.txt.gz"),
  
  #NeuN (17)
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877174_1524_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877198_1525_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877200_1527_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877201_1532_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877202_1536_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877204_1539_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877205_1541_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877176_3545_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877178_3586_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877228_3590_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877206_3602_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877186_4615_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877208_AN03398_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877213_AN05483_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877215_AN10090_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877209_AN15240_NeuN_CpG_WGBS.txt.gz"),
  file.path(PATH, "Dataset_2", "2_BS-Seq", "neun", "GSM2877210_AN16799_NeuN_CpG_WGBS.txt.gz")
)

#Naming samples:
sample_names <- c(
  "rep_1_olig2",  "rep_2_olig2",  "rep_3_olig2",  "rep_4_olig2",
  "rep_5_olig2",  "rep_6_olig2",  "rep_7_olig2",  "rep_8_olig2",
  "rep_9_olig2",  "rep_10_olig2", "rep_11_olig2", "rep_12_olig2",
  "rep_13_olig2", "rep_14_olig2", "rep_15_olig2", "rep_16_olig2",
  "rep_17_olig2",
  "rep_1_neun",   "rep_2_neun",   "rep_3_neun",   "rep_4_neun",
  "rep_5_neun",   "rep_6_neun",   "rep_7_neun",   "rep_8_neun",
  "rep_9_neun",   "rep_10_neun",  "rep_11_neun",  "rep_12_neun",
  "rep_13_neun",  "rep_14_neun",  "rep_15_neun",  "rep_16_neun",
  "rep_17_neun"
)

dir.create("bismark_cov", showWarnings = FALSE)
dir.create("methylDB", showWarnings = FALSE)


# a. Object building: saving on disk in bismarkCoverage format
#    chr  start  end  meth%  count_M  count_U
for (i in seq_along(files)) {
  df <- fread(files[i], header = TRUE, col.names = c("chromosome", "CpG_position", "M_counts", "Total_counts"))
  df <- df[grepl("^chr", chromosome) & Total_counts > 0]   # removes spike-in e coverage 0
  out <- data.table(
    chr     = df$chromosome,
    start   = df$CpG_position,
    end     = df$CpG_position,
    meth    = round(100 * df$M_counts / df$Total_counts, 4),
    count_M = df$M_counts,
    count_U = df$Total_counts - df$M_counts
  )
  fwrite(out,
         file.path("bismark_cov", paste0(sample_names[i], ".cov")),
         sep = "\t", col.names = FALSE)
  rm(df, out); gc()
}

# b. Object building: methRead with backend tabix: lighter on RAM
methyl_obj <- methRead(
  location  = as.list(file.path("bismark_cov", paste0(sample_names, ".cov"))),
  sample.id = as.list(sample_names),
  assembly  = "hg19",                
  treatment = c(rep(1, 17), rep(0, 17)),    # 1 = Olig2, 0 = Neun
  context   = "CpG",
  pipeline  = "bismarkCoverage",
  mincov    = 1,
  dbtype    = "tabix",
  dbdir     = "methylDB"
)


## Stats visualization
#Methylation
 par(mfrow = c(3, 4))
 for (i in 1:34) {
   getMethylationStats(methyl_obj[[i]], plot = TRUE, both.strands = FALSE)
 }
 
#Coverage
 par(mfrow = c(3, 4))
 for (i in 1:34) {
   getCoverageStats(methyl_obj[[i]], plot = TRUE, both.strands = FALSE)
 }

#Resetting layout to normal
par(mfrow = c(1, 1))


### 2. Filtering and uniting ###

#Analyzing the data:
# 1. Checking means of coverage
cov_files <- file.path("bismark_cov", paste0(sample_names, ".cov"))

cov_means <- sapply(cov_files, function(f) {
  dt <- fread(f, col.names = c("chr","start","end","meth","count_M","count_U"))
  mean(dt$count_M + dt$count_U)
})
names(cov_means) <- sample_names
print(cov_means)

# 2. Percentiles
set.seed(42)
cov_sample <- unlist(lapply(methyl_obj, function(x) {
  cov <- getData(x)$coverage
  sample(cov, min(1e5, length(cov)))
}))
quantile(cov_sample, c(0.99, 0.999, 0.9999, 0.99999))

#Based on these quantiles, the coverage filter must be decided. Standard is 10, but if the medians are lower it must be changed
#PCR artifacts abundance will instead set ceiling threshold.

filtered_methyl_obj=filterByCoverage(methyl_obj,lo.count=10,lo.perc=NULL,     #lowest mean is at 15. Others are between 20-30
                                     hi.count=NULL,hi.perc=99.99,
                                     save.db = TRUE,
                                     suffix = "filtered")

rm(methyl_obj)

## Merging
meth=unite(filtered_methyl_obj, destrand=FALSE,
           save.db = TRUE,
           suffix = "united")

rm(filtered_methyl_obj)

## Explorative analysis on the merged:    

#This dataset is too big. Must be filtered before. Filtering:
#Extracting methyl matrix
pm <- percMethylation(meth)
#SD per CpG
sds <- rowSds(pm, na.rm = TRUE)
#Indexing top CpG most variable
top_n <- 10000
top_idx <- order(sds, decreasing = TRUE)[seq_len(min(top_n, length(sds)))]
#Subset of meth object
meth_pca <- meth[top_idx, ]

## Plotting: Clustering and PCA
clusterSamples(meth_pca, dist = "correlation", method = "ward.D2", plot = TRUE)
PCASamples(meth_pca)



### 3. Differential Analysis ###

#Adding metadata for covariates
metadata <- data.frame(
  sample = c(paste0("rep_", 1:17, "_olig2"),
             paste0("rep_", 1:17, "_neun")),
  group = c(rep("olig2", 17), rep("neun", 17)),
  donor = c(
    #Olig (17): 
    1:17,
    #Neun (17):
    1:17
  ),
  stringsAsFactors = FALSE
)


#Covariates:
covariates <- data.frame(donor = factor(metadata$donor, levels = 1:17))

#Doing the analysis with the correction for overdispersion: "MN".
myDiff <- calculateDiffMeth(meth,
                            covariates = covariates,
                            overdispersion = "MN",
                            test = "F",
                            save.db = TRUE,
                            suffix = "diff")            



### !!! From this point on, moved to another work station. Files were inserted in path: file.path(PATH, "Dataset_2", "2_BS-Seq") ###
## Reading:
#myDiff <- readMethylDB(file.path(PATH, "Dataset_2", "2_BS-Seq", "methylDiff_united_diff.txt.bgz"))


## Selecting differentially methylated bases:
#All differentially methylated bases:
myDiff25p=getMethylDiff(myDiff, difference=25, qvalue=0.01, suffix = "25p")

#Hyper and hypo-only DM:
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper", suffix = "25p_hyper")
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo", suffix = "25p_hypo")

## Graph: volcano plot

#Too many sites for the plot. Filtering high qvalue and low meth.diff values. Subsampling:
SUBSAMPLE_N <- 1000
set.seed(42)

myDiff_df <- as.data.frame(as(myDiff, "GRanges"))
myDiff_df_plot <- myDiff_df[!is.na(myDiff_df$qvalue) & 
                              myDiff_df$qvalue < 0.4 & 
                              abs(myDiff_df$meth.diff) > 4, ]

if (nrow(myDiff_df_plot) > SUBSAMPLE_N) {
  myDiff_df_plot <- myDiff_df_plot[sample(nrow(myDiff_df_plot), SUBSAMPLE_N), ]
}

#Reordering: significant on top
is_sig <- myDiff_df_plot$qvalue < 0.01 & abs(myDiff_df_plot$meth.diff) > 25
myDiff_df_plot <- rbind(myDiff_df_plot[!is_sig, ], myDiff_df_plot[is_sig, ])
df_sig <- myDiff_df_plot[myDiff_df_plot$qvalue < 0.01 & abs(myDiff_df_plot$meth.diff) > 25, ]

#Plotting:
gg <- ggplot(myDiff_df_plot, aes(x = meth.diff, y = -log10(qvalue))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-25, 25), linetype = "dashed", color = "blue") +
  geom_point(data = df_sig, color = "red") +
  theme_minimal() +
  labs(title = sprintf("Volcano plot: D2 DM analysis",
                       nrow(myDiff_df_plot), nrow(df_sig)),
       x = "meth.diff", y = "-log10(adj pvalue)")
print(gg)



### 4. Coordinate conversion, annotation and export of DM sites ###

##Liftover
## Before doing the annotation, it is important to note that the methyl calling was obtained with the hg19 genome
## But the RNA-Seq data quantification was done using the hg38. For consistency it is good to switch to the hg38 genome.
#Doing this with liftOver()

#Downloading chain and meth files
chain <- import.chain(file.path(PATH, "references", "hg19ToHg38.over.chain"))
meth <- readMethylDB(file.path(PATH, "Dataset_2", "2_BS-Seq", "methylBase_united.txt.bgz"))

#meth_dm only has DM sites
dm_hg19_GR <- as(myDiff25p, "GRanges")
meth_dm <- selectByOverlap(meth, dm_hg19_GR)

meth_dm_GR <- as(meth_dm, "GRanges")


### Now the liftovers: hg19 -> hg38

## First, liftover and export of meth25p
seqlevelsStyle(meth_dm_GR) <- "UCSC"
meth_DM_hg38 <- unlist(liftOver(meth_dm_GR, chain))

#Exporting:
meth_DM_hg38_df <- as.data.frame(meth_DM_hg38)
meth_DM_hg38_df$coord_key <- paste(meth_DM_hg38_df$seqnames, meth_DM_hg38_df$start, sep="_")
write.csv(meth_DM_hg38_df, file.path(PATH, "Dataset_2", "2_BS-Seq", "meth25p.csv"), row.names = FALSE)


## Then, liftover and export of myDiff25p:
myDiff25p_GR <- as(myDiff25p, "GRanges")
seqlevelsStyle(myDiff25p_GR) <- "UCSC"
myDiff25p_GR_hg38 <- unlist(liftOver(myDiff25p_GR, chain))
#Adding coord_key
myDiff25p_GR_hg38$coord_key <- paste(seqnames(myDiff25p_GR_hg38), start(myDiff25p_GR_hg38), sep = "_")

#Exporting DM sites:
saveRDS(myDiff25p_GR_hg38, file = file.path(PATH, "Dataset_2", "2_BS-Seq", "myDiff25p_GR_hg38.rds"))

dev.off()



## Karyo visualization
pdf(file.path(PATH, "Dataset_2", "2_BS-Seq", "Chr DM distribution D2.pdf"), width = 12, height = 8)

## Graph: Position of methylation sites on all chromosomes
## Checking if their position is clusterized around centromeres.

kp <- plotKaryotype(genome="hg38")
kp <- kpPlotDensity(kp, myDiff25p_GR_hg38)

dev.off()