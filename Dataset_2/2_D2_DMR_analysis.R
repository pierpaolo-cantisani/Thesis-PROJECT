library(data.table)
library(methylKit)
library(ggplot2)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomeInfoDb)
library(ChIPseeker)
library(rtracklayer)
library(dplyr)
library(karyoploteR)

setwd("C:/Users/pierp/Desktop/THESIS PROJECT/Dataset_2/2_BS-Seq")

#pdf("C:/Users/pierp/Desktop/THESIS PROJECT/Dataset_2/2_BS-Seq/WGBS Graphs D2.pdf", width = 12, height = 8)



### 1. File downloading, and MethylKit object building ###

#The data files are not in a typical bismark output format, but they were processed and simplified.
#Files only have 4 columns: "chromosome", "position", "counts M", "counts M+U (coverage)".

#names of the import files
files <- c(
  # NeuN Schizo (21)
  "schizo/GSM2877179_1505_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877188_1506_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877189_1507_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877191_1508_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877192_1510_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877193_1511_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877171_1512_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877194_1513_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877195_1514_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877196_1515_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877197_1518_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877181_1523_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877175_4336_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877185_4361_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877223_4395_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877227_4448_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877244_4730_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877241_4804_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877243_AN09634_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877232_AN17799_NeuN_CpG_WGBS.txt.gz",
  "schizo/GSM2877245_AN18099_NeuN_CpG_WGBS.txt.gz",
# NeuN Control (19)
  "control/GSM2877174_1524_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877198_1525_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877200_1527_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877201_1532_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877202_1536_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877203_1538_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877204_1539_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877205_1541_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877176_3545_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877178_3586_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877228_3590_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877206_3602_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877207_3611_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877186_4615_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877208_AN03398_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877213_AN05483_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877215_AN10090_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877209_AN15240_NeuN_CpG_WGBS.txt.gz",
  "control/GSM2877210_AN16799_NeuN_CpG_WGBS.txt.gz")


#Naming samples:
sample_names <- c(
  "rep_1_schizo",  "rep_2_schizo",  "rep_3_schizo",  "rep_4_schizo",
  "rep_5_schizo",  "rep_6_schizo",  "rep_7_schizo",  "rep_8_schizo",
  "rep_9_schizo",  "rep_10_schizo", "rep_11_schizo", "rep_12_schizo",
  "rep_13_schizo", "rep_14_schizo", "rep_15_schizo", "rep_16_schizo",
  "rep_17_schizo", "rep_18_schizo", "rep_19_schizo", "rep_20_schizo",
  "rep_21_schizo",
  "rep_1_ctrl",  "rep_2_ctrl",  "rep_3_ctrl",  "rep_4_ctrl",
  "rep_5_ctrl",  "rep_6_ctrl",  "rep_7_ctrl",  "rep_8_ctrl",
  "rep_9_ctrl",  "rep_10_ctrl", "rep_11_ctrl", "rep_12_ctrl",
  "rep_13_ctrl", "rep_14_ctrl", "rep_15_ctrl", "rep_16_ctrl",
  "rep_17_ctrl", "rep_18_ctrl", "rep_19_ctrl")

dir.create("bismark_cov", showWarnings = FALSE)
dir.create("methylDB", showWarnings = FALSE)


# a. Object building: saving on disk in bismarkCoverage format
#    chr  start  end  meth%  count_M  count_U
for (i in seq_along(files)) {
  df <- fread(files[i], col.names = c("chromosome", "CpG_position", "M_counts", "Total_counts"))
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
  treatment = c(rep(1, 21), rep(0, 19)),   # 1 = schizo (21), 0 = ctrl (19)      
  context   = "CpG",
  pipeline  = "bismarkCoverage",
  mincov    = 1,
  dbtype    = "tabix",
  dbdir     = "methylDB"
)


#Stats visualization
#Methylation
par(mfrow = c(3, 4))
for (i in 1:40) {
  getMethylationStats(methyl_obj[[i]], plot = TRUE, both.strands = FALSE)
}

#Coverage
par(mfrow = c(3, 4))
for (i in 1:40) {
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
##Merging
meth=unite(filtered_methyl_obj, destrand=FALSE,
           min.per.group = 4L,
           save.db = TRUE,
           suffix = "united")

rm(filtered_methyl_obj)
##Explorative analysis on the merged:    #This dataset is too big. Must be filtered before
#clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
#PCASamples(meth)





### 3. Differential Analysis ###

#Adding metadata for covariates
metadata <- data.frame(
  sample = c(paste0("rep_", 1:21, "_schizo"),
             paste0("rep_", 1:19, "_ctrl")),     
  group = c(rep("schizo", 21), rep("ctrl", 19)), 
  sex = c(#schizo
    "F","F","M","F","M","M","F","F",
    "M","F","F","F","F","F","M","M",
    "F","M","M","F","M",
    #ctrl
    "M","F","M","F","M","F","F","F",
    "M","M","M","M","M","M",
    "F","M","M","F","M"),
  stringsAsFactors = FALSE
)


# Update sample_names di conseguenza (devi farlo anche a monte!)
sample_names <- metadata$sample


# Covariates:
covariates <- data.frame(
  sex       = factor(metadata$sex,       levels = c("F", "M")))

#Doing the analysis with the correction for overdispersion: "MN".
myDiff <- calculateDiffMeth(meth,
                            covariates = covariates,
                            overdispersion = "MN",
                            test = "Chisq",           #With correction the default is the F test: must force this for the comparison
                            save.db = TRUE,
                            suffix = "diff")            

##Reading
myDiff <- readMethylDB("C:/Users/pierp/Desktop/THESIS PROJECT/Dataset_2/2_BS-Seq/methylDiff_united_diff.txt.bgz")


##Finally: selecting differentially methylated bases:
#get all differentially methylated bases
myDiff10p=getMethylDiff(myDiff, difference=10, qvalue=0.05, suffix = "10p")
#get yper and hypo-only DM
myDiff10p.hyper=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hyper", suffix = "10p_hyper")
myDiff10p.hypo=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hypo", suffix = "10p_hypo")

#myDiff10_df <- data.frame(as(myDiff10p, "GRanges"))


##Volcano plot
myDiff_df <- as.data.frame(as(myDiff, "GRanges"))
#Too many sites for the plot. Filtering high qvalue and low meth.diff values:
myDiff_df_plot <- myDiff_df[!is.na(myDiff_df$qvalue) & 
                                  myDiff_df$qvalue < 0.6 & 
                                  abs(myDiff_df$meth.diff) > 2, ]


ggplot(myDiff_df_plot, aes(x = meth.diff, y = -log10(qvalue))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-10, 10), linetype = "dashed", color = "blue") +
  geom_point(data = myDiff_df_plot[!is.na(myDiff_df_plot$qvalue) & 
                                     myDiff_df_plot$qvalue < 0.05 & 
                                     abs(myDiff_df_plot$meth.diff) > 10, ], color = "red") +
  theme_minimal() +
  labs(title = "Volcano plot: infection", x = "meth.diff", y = "-log10(adj pvalue)")



### 4. Coordinate conversion, annotation and export of the universe and DM sites ###

##Liftover
## Before doing the annotation, it is important to note that the methyl calling was obtained with the h19 genome
## But the RNA-Seq data quantification was done using the h38. For consistency it is good to switch to the h38 genome.
#Doing this with liftOver()

#Download chain and meth files
chain <- import.chain("C:/Users/pierp/Desktop/THESIS PROJECT/references/hg19ToHg38.over.chain")
meth <- readMethylDB("C:/Users/pierp/Desktop/THESIS PROJECT/Dataset_2/2_BS-Seq/methylBase_united.txt.bgz")

dm_hg19_GR <- as(myDiff10p, "GRanges")
meth_dm <- selectByOverlap(meth, dm_hg19_GR)

# meth_dm only has DM sites
meth_dm_GR <- as(meth_dm, "GRanges")


### Now the liftovers: hg19 -> hg38

##First, liftover and exporting the meth10p
seqlevelsStyle(meth_dm_GR) <- "UCSC"
meth_DM_hg38 <- unlist(liftOver(meth_dm_GR, chain))

#Now exporting
meth_DM_hg38_df <- as.data.frame(meth_DM_hg38)
meth_DM_hg38_df$coord_key <- paste(meth_DM_hg38_df$seqnames, meth_DM_hg38_df$start, sep="_")
write.csv(meth_DM_hg38_df, "C:/Users/pierp/Desktop/THESIS PROJECT/Dataset_2/2_BS-Seq/meth10p.csv", row.names = FALSE)


## Then, liftover for myDiff10p:
myDiff10p_GR <- as(myDiff10p, "GRanges")
seqlevelsStyle(myDiff10p_GR) <- "UCSC"
myDiff10p_GR_hg38 <- unlist(liftOver(myDiff10p_GR, chain))
#df and coord_key
myDiff10p_GR_hg38$coord_key <- paste(seqnames(myDiff10p_GR_hg38), start(myDiff10p_GR_hg38), sep = "_")

#Exporting DM sites:
saveRDS(myDiff10p_GR_hg38, file = "C:/Users/pierp/Desktop/THESIS PROJECT/Dataset_2/2_BS-Seq/myDiff10p_GR_hg38.rds")

#dev.off()



#Karyo visualization
pdf("C:/Users/pierp/Desktop/THESIS PROJECT/Dataset_2/2_BS-Seq/Chr DM distribution D2.pdf", width = 12, height = 8)

## Graph: Position of methylation sites on all chromosomes
## Checking if their position is clusterized around centromeres. Then considering filtering

kp <- plotKaryotype(genome="hg38")
kp <- kpPlotDensity(kp, myDiff10p_GR_hg38)

dev.off()