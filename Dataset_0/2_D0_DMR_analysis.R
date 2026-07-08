library(data.table)
library(methylKit)
library(ggplot2)
library(rtracklayer)
library(karyoploteR)


PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"

pdf(file.path(PATH, "Dataset_0", "2_BS-Seq", "WGBS Graphs D0.pdf"), width = 12, height = 8)

### 1. File downloading, and MethylKit object building ###

#The data files are not in a typical bismark output format, but they were processed and simplified.
#Files only have 4 columns: "chromosome", "position", "counts M", "counts M+U (coverage)".

#Setting file directory
setwd(file.path(PATH, "Dataset_0", "2_BS-Seq"))

#names of the import files
files_df <- read.delim(file.path(PATH, "Dataset_0", "Ground Truth", "samples.tsv"))
files_tmp <- paste0(files_df$sample, ".txt")
files <- file.path(PATH, "Dataset_0", "2_BS-Seq", "wgbs", files_tmp)
names(files) <- paste0("S", seq_len(nrow(files_df)))

#Naming samples:
sample_names = files_df$sample


##Using MethylSeq.This tool expects a different format, so creating the methylKit object:
methyl_obj <- new("methylRawList",                          #methylraw class is defined in the methylKit library
                  lapply(seq_along(files), function(i) {
                    df <- fread(files[i], header = TRUE, col.names = c("chr", "pos", "M", "M+U"))
                    new("methylRaw",
                        data.frame(chr = df$chr,
                                   start = df$pos,
                                   end = df$pos,               #all(start(BS_obj) == end(BS_obj)) is TRUE --> start and end are the same, cause it's only a position.
                                   strand = "*",               #No strand specifics in the original files
                                   coverage = df$'M+U',
                                   numCs = df$M,
                                   numTs = df$'M+U' - df$M
                        ),
                        sample.id = sample_names[i],
                        assembly = "hg38",                    
                        context = "CpG",
                        resolution = "base")
                  }
                  ),
                  treatment = as.integer(files_df$condition == "case")    # 1 = "case", 0 = "ctrl"
)


#Stats visualization
#Methylation
par(mfrow = c(5, 5))
for (i in 1:100) {
  getMethylationStats(methyl_obj[[i]], plot = TRUE, both.strands = FALSE)
}
# 
#Coverage
par(mfrow = c(5, 5))
for (i in 1:100) {
  getCoverageStats(methyl_obj[[i]], plot = TRUE, both.strands = FALSE)
}

#Resetting layout to normal
par(mfrow = c(1, 1))


### 2. Filtering and uniting ###

#Analyzing the data:
# 1. Checking means of coverage
sapply(methyl_obj, function(x) mean(getData(x)$coverage)) 
# 2. Percentiles
cov_max <- do.call(pmax, lapply(methyl_obj, function(x) getData(x)$coverage))
quantile(cov_max, c(0.99, 0.999, 0.9999, 0.99999))

#Based on these quantiles, the coverage filter must be decided. Standard is 10, but if the medians are lower it must be changed
#PCR artifacts abundance will instead set ceiling threshold.

##Filtering: NO filtering for the positive control dataset
#filtered_methyl_obj=filterByCoverage(methyl_obj,lo.count=5,lo.perc=NULL,
#                                     hi.count=NULL,hi.perc=99.99)

##Merging
meth=unite(methyl_obj, destrand=FALSE)

rm(methyl_obj)

##Explorative analysis on the merged:
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
PCASamples(meth)

#If clustering doesn't follow the expected paired design, a batch effect may be present.




### 3. Differential Analysis ###
covariates <- data.frame(pair = factor(files_df$pair))

#Doing the analysis with the correction for overdispersion: "MN".
myDiff <- calculateDiffMeth(meth,
                            covariates = covariates,
                            overdispersion = "MN",
                            test = "F")


##Finally: selecting differentially methylated bases:
#get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
#get hyper and hypo
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")


##Volcano plot
myDiff_df <- as.data.frame(as(myDiff, "GRanges"))


ggplot(myDiff_df, aes(x = meth.diff, y = -log10(qvalue))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-25, 25), linetype = "dashed", color = "blue") +
  geom_point(data = myDiff_df[!is.na(myDiff_df$qvalue) & 
                                     myDiff_df$qvalue < 0.01 & 
                                     abs(myDiff_df$meth.diff) > 25, ], color = "red") +
  theme_minimal() +
  labs(title = "Volcano plot: D0 Differential methylation", x = "meth.diff", y = "-log10(adj pvalue)")



### 4. Annotation and export of the universe and DM sites ###

dm_GR <- as(myDiff25p, "GRanges")
meth_dm <- selectByOverlap(meth, dm_GR)

# meth_dm only has DM sites
meth_dm_GR <- as(meth_dm, "GRanges")

#Now exporting
meth_DM_df <- as.data.frame(meth_dm_GR)
meth_DM_df$coord_key <- paste(meth_DM_df$seqnames, meth_DM_df$start, sep="_")
write.csv(meth_DM_df, file.path(PATH, "Dataset_0", "2_BS-Seq", "meth25p.csv"), row.names = FALSE)


#Now myDiff25p
myDiff25p_GR <- as(myDiff25p, "GRanges")
myDiff25p_GR$coord_key <- paste(seqnames(myDiff25p_GR), start(myDiff25p_GR), sep = "_")

#Exporting DM sites:
saveRDS(myDiff25p_GR, file = file.path(PATH, "Dataset_0", "2_BS-Seq", "myDiff25p_GR.rds"))

dev.off()


#Karyo visualization
pdf(file.path(PATH, "Dataset_0", "2_BS-Seq", "Chr DM distribution D0.pdf"), width = 12, height = 8)

## Graph: Position of methylation sites on all chromosomes
## Checking if their position is clusterized around centromeres. Then considering filtering

kp <- plotKaryotype(genome="hg38")
kp <- kpPlotDensity(kp, myDiff25p_GR)

dev.off()


##Output data:
Output_df <- read.csv(file.path(PATH, "Dataset_0", "Ground Truth", "Output_data_tmp_1.csv"))

#Generating new output:
Output_df <- rbind(Output_df, data.frame(metric = c("CpG sites", "CpG DM sites"), value = c(nrow(getData(meth)), nrow(getData(myDiff25p)))))

write.csv(Output_df, file.path(PATH, "Dataset_0", "Ground Truth", "Output_data_tmp_2.csv"), row.names = FALSE)