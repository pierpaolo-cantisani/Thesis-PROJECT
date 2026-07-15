library(data.table)
library(methylKit)
library(ggplot2)
library(rtracklayer)
library(karyoploteR)


PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"

pdf(file.path(PATH, "Dataset_1", "2_BS-Seq", "WGBS Graphs D1.pdf"), width = 12, height = 8)

### 1. File downloading, and MethylKit object building ###

#The data files are not in a typical bismark output format, but they were processed and simplified.
#Files only have 4 columns: "chromosome", "position", "counts M", "counts M+U (coverage)".

#Setting file directory
setwd(file.path(PATH, "Dataset_1", "2_BS-Seq"))

#Names of the import files
files <- c(
  "GSM1565939_DC81_MTB_5mC.txt.gz", "GSM1565940_DC81_NI_5mC.txt.gz", "GSM1565941_DC82_MTB_5mC.txt.gz", 
  "GSM1565942_DC82_NI_5mC.txt.gz", "GSM1565943_DC83_MTB_5mC.txt.gz", "GSM1565944_DC83_NI_5mC.txt.gz",
  "GSM1565945_DC87_MTB_5mC.txt.gz", "GSM1565946_DC87_NI_5mC.txt.gz", "GSM1565947_DC89_MTB_5mC.txt.gz",
  "GSM1565948_DC89_NI_5mC.txt.gz", "GSM1565949_DC91_MTB_5mC.txt.gz", "GSM1565950_DC91_NI_5mC.txt.gz")

#Naming samples:
sample_names = c("rep_1_TB", "rep_1_CTRL", "rep_2_TB", "rep_2_CTRL", 
                 "rep_3_TB", "rep_3_CTRL", "rep_4_TB", "rep_4_CTRL", 
                 "rep_5_TB", "rep_5_CTRL", "rep_6_TB", "rep_6_CTRL")


## Methylkit expects a different format: creating the MethylKit input object
methyl_obj <- new("methylRawList",                             #methylraw class is defined in the methylKit library
                  lapply(seq_along(files), function(i) {
                    df <- fread(files[i], header = TRUE, col.names = c("chr", "pos", "M", "M+U"))
                    # removing spike-ins
                    df <- df[grepl("^chr", chr)]
                    new("methylRaw",
                        data.frame(chr = df$chr,
                                   start = df$pos,
                                   end = df$pos,               #Singl CpG --> start == end
                                   strand = "*",               #No strand specifics in the original files
                                   coverage = df$'M+U',
                                   numCs = df$M,
                                   numTs = df$'M+U' - df$M
                        ),
                        sample.id = sample_names[i],
                        assembly = "hg19",                     #Data are from 2015. The ref is hg19
                        context = "CpG",
                        resolution = "base")
                  }
                  ),
                  treatment = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0)    # 1 = "TB", 0 = "CTRL"
)


## Stats visualization

#Methylation
par(mfrow = c(3, 4))
for (i in 1:12) {
  getMethylationStats(methyl_obj[[i]], plot = TRUE, both.strands = FALSE)
}

#Coverage
par(mfrow = c(3, 4))
for (i in 1:12) {
  getCoverageStats(methyl_obj[[i]], plot = TRUE, both.strands = FALSE)
}

#Resetting layout to normal
par(mfrow = c(1, 1))



### 2. Filtering and uniting ###

## Analyzing the data:
#1. Checking means of coverage
sapply(methyl_obj, function(x) mean(getData(x)$coverage)) 
#2. Percentiles
cov_max <- do.call(pmax, lapply(methyl_obj, function(x) getData(x)$coverage))
quantile(cov_max, c(0.99, 0.999, 0.9999, 0.99999))

#Based on these quantiles, the coverage filter must be decided. Standard is 10, but if the medians are lower it must be changed
#PCR artifacts abundance will instead set ceiling threshold.
#Setting coverage filter at 5 (lower median is near 6)

filtered_methyl_obj=filterByCoverage(methyl_obj,lo.count=5,lo.perc=NULL,
                                     hi.count=NULL,hi.perc=99.99)

rm(methyl_obj)

## Merging
meth=unite(filtered_methyl_obj, destrand=FALSE)

rm(filtered_methyl_obj)

## Explorative analysis on the merged:
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
PCASamples(meth)



### 3. Differential Analysis ###

covariates <- data.frame(replicate = factor(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)))

#Doing the analysis with the correction for overdispersion: "MN".
myDiff <- calculateDiffMeth(meth,
                            covariates = covariates,
                            overdispersion = "MN",
                            test = "F")


## Diff Analysis was done on another machine. Then imported:
myDiff <- readRDS(file.path(PATH, "Dataset_1", "2_BS-Seq", "myDiff.rds"))


## Selecting differentially methylated bases:
#All differentially methylated bases:
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

#Hyper and hypo-only DM:
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")



## Graph: volcano plot
myDiff_df <- as.data.frame(as(myDiff, "GRanges"))
#Too many sites for the plot. Filtering high qvalue and low meth.diff values:
myDiff_df_plot <- myDiff_df[!is.na(myDiff_df$qvalue) & 
                                  myDiff_df$qvalue < 0.6 & 
                                  abs(myDiff_df$meth.diff) > 2, ]

ggplot(myDiff_df_plot, aes(x = meth.diff, y = -log10(qvalue))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-25, 25), linetype = "dashed", color = "blue") +
  geom_point(data = myDiff_df_plot[!is.na(myDiff_df_plot$qvalue) & 
                                     myDiff_df_plot$qvalue < 0.01 & 
                                     abs(myDiff_df_plot$meth.diff) > 25, ], color = "red") +
  theme_minimal() +
  labs(title = "Volcano plot: DM analysis", x = "meth.diff", y = "-log10(adj pvalue)")



### 4. Coordinate conversion, annotation and export of the universe and DM sites ###

## Liftover
#Before doing the annotation, it is important to note that the methyl calling was obtained with the h19 genome
#But the RNA-Seq data quantification was done using the h38. For consistency it is good to switch to the h38 genome.
#Doing this with liftOver()

#Downloading chain and meth files
chain <- import.chain(file.path(PATH, "references", "hg19ToHg38.over.chain"))
#This (meth) was obtained from another machine too. Importing:
meth <- readRDS(file.path(PATH, "Dataset_1", "2_BS-Seq", "meth.rds"))

dm_hg19_GR <- as(myDiff25p, "GRanges")
meth_dm <- selectByOverlap(meth, dm_hg19_GR)

#meth_dm only has DM sites
meth_dm_GR <- as(meth_dm, "GRanges")


## Now the liftovers: hg19 -> hg38
#First, liftover and exporting the meth25p
seqlevelsStyle(meth_dm_GR) <- "UCSC"
meth_DM_hg38 <- unlist(liftOver(meth_dm_GR, chain))

#Now exporting
meth_DM_hg38_df <- as.data.frame(meth_DM_hg38)
meth_DM_hg38_df$coord_key <- paste(meth_DM_hg38_df$seqnames, meth_DM_hg38_df$start, sep="_")

write.csv(meth_DM_hg38_df, file.path(PATH, "Dataset_1", "2_BS-Seq", "meth25p.csv"), row.names = FALSE)


## Then, liftover for myDiff25p:
myDiff25p_GR <- as(myDiff25p, "GRanges")
seqlevelsStyle(myDiff25p_GR) <- "UCSC"
myDiff25p_GR_hg38 <- unlist(liftOver(myDiff25p_GR, chain))
#df and coord_key
myDiff25p_GR_hg38$coord_key <- paste(seqnames(myDiff25p_GR_hg38), start(myDiff25p_GR_hg38), sep = "_")

#Sanity check. Are liftovers the same?
stopifnot(setequal(meth_DM_hg38_df$coord_key, myDiff25p_GR_hg38$coord_key))


#Exporting DM sites:
saveRDS(myDiff25p_GR_hg38, file = file.path(PATH, "Dataset_1", "2_BS-Seq", "myDiff25p_GR_hg38.rds"))

dev.off()


#Karyo visualization
pdf(file.path(PATH, "Dataset_1", "2_BS-Seq", "Chr DM distribution D1.pdf"), width = 12, height = 8)

## Graph: Position of methylation sites on all chromosomes
## Checking if their position is clusterized around centromeres. Then considering filtering

kp <- plotKaryotype(genome="hg38")
kp <- kpPlotDensity(kp, myDiff25p_GR_hg38)

dev.off()