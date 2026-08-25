library(data.table)
library(stringr)
library(limma)
library(ggplot2)
library(GenomicRanges)
library(karyoploteR)

PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"

pdf(file.path(PATH, "Dataset_3", "2_BS-Seq", "WGBS Graphs D3.pdf"), height = 10, width = 15)

### 1. Importing files and transforming into M-value ###

meth_D0 <- fread(file.path(PATH, "Dataset_3", "2_BS-Seq", "GSE263782_Day0_methyl_matrix.csv.gz"))
meth_D28 <- fread(file.path(PATH, "Dataset_3", "2_BS-Seq", "GSE263782_Day28_methyl_matrix.csv.gz"))

#Sorting samples numerically:
meth_D0  <- meth_D0[order(as.integer(sub("_D0$",  "", X)))]
meth_D28 <- meth_D28[order(as.integer(sub("_D28$", "", X)))]

#Transforming into matrix (limma format has samples as columns)
D0_matrix <- t(as.matrix(meth_D0[, -"X"]))
colnames(D0_matrix) <- meth_D0$X

D28_matrix <- t(as.matrix(meth_D28[, -"X"]))
colnames(D28_matrix) <- meth_D28$X

## Merging
full_matrix <- cbind(D0_matrix, D28_matrix)


## Filtering and clipping

#Filtering for variance: So CpG with steady b-value are eliminated (threshold: 0.02)
sd_beta <- apply(full_matrix, 1, sd, na.rm = TRUE)
keep <- sd_beta > 0.02 & !is.na(sd_beta)
table(keep)  # quanti CpG sopravvivono?
full_matrix    <- full_matrix[keep, ]

#Clipping: modifying all extreme beta values to 0.999(if 1) and 0.001(if 0). This allows the transformation to M-values to run smoothly
full_matrix_clip <- pmax(pmin(full_matrix , 0.999), 0.001)


#Now transformation into M-values:
M_CpG_matrix <- log2(full_matrix_clip/(1-full_matrix_clip))



### 2. DM analysis: limma ###

## Metadata
sample_ids <- colnames(M_CpG_matrix)
metadata <- data.frame(
  condition = factor(str_extract(sample_ids, "D0$|D28$"), levels = c("D0","D28")),
  donor   = factor(str_extract(sample_ids, "^[^_]+"), levels = as.character(sort(as.integer(unique(str_extract(sample_ids, "^[^_]+")))))),
  row.names = sample_ids
)

## Design and limma
design <- model.matrix(~ donor + condition, data = metadata)

fit  <- lmFit(M_CpG_matrix, design)
fit  <- eBayes(fit, robust = TRUE, trend = TRUE)

#Extracting fit results
res <- topTable(fit,
                coef       = "conditionD28",   #last in "colnames(design)"
                number     = Inf,
                sort.by    = "none",        
                adjust.method = "BH")

head(res)


## Obtaining delta-B to assess DM sites:
D0_cols <- metadata$condition == "D0"
D28_cols <- metadata$condition == "D28"

D0_mean <- rowMeans(full_matrix_clip[, D0_cols])
D28_mean <- rowMeans(full_matrix_clip[, D28_cols])

deltaB <- D28_mean - D0_mean
res$deltaB <- deltaB


## Extracting significant results:
sign_DM <- res[res$adj.P.Val < 0.05 & abs(res$deltaB) > 0.05, ]
sign_DM$seqnames <- sub("\\..*", "", rownames(sign_DM))
sign_DM$start <- as.integer(sub(".*\\.", "", rownames(sign_DM)))

res$meth.diff <- deltaB*100
## Graph: volcano plot
ggplot(res, aes(x = meth.diff, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color ="blue") +
  geom_vline(xintercept = c(-5, +5), linetype = "dashed", color ="blue") +
  geom_point(data = res[res$adj.P.Val < 0.05 & abs(res$meth.diff) > 5, ], col = "red") +
  theme_minimal() +
  labs(title = "Volcano plot: DM analysis", x= "Delta beta value", y = "-log10(adj pvalue)")


## Exporting
sign_DM$coord_key <- paste(sign_DM$seqnames, sign_DM$start, sep = "_")
write.csv(sign_DM, file.path(PATH, "Dataset_3", "2_BS-Seq", "sign_DM.csv"), row.names = FALSE)

#Changing rownames and exporting M-values
seqnames <- sub("\\..*", "", rownames(M_CpG_matrix))
start <- as.integer(sub(".*\\.", "", rownames(M_CpG_matrix)))
rownames(M_CpG_matrix) <- paste(seqnames, start, sep="_")
write.csv(M_CpG_matrix, file.path(PATH, "Dataset_3", "2_BS-Seq", "M_CpG_matrix.csv"), row.names = TRUE)

dev.off()


## Karyo visualization
#Obtaining GRanges
sign_DM_GR <- makeGRangesFromDataFrame(data.frame(seqnames = sign_DM$seqnames,
                                                  start = sign_DM$start,
                                                  end = sign_DM$start))
seqlevelsStyle(sign_DM_GR) <- "UCSC"

pdf(file.path(PATH, "Dataset_3", "2_BS-Seq", "Chr DM distribution D3.pdf"), width = 12, height = 8)

## Graph: Position of methylation sites on all chromosomes
## Checking if their position is clusterized around centromeres. Then considering filtering

kp <- plotKaryotype(genome="hg38")
kp <- kpPlotDensity(kp, sign_DM_GR)

dev.off()