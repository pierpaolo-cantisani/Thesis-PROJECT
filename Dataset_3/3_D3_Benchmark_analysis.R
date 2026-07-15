library(genomation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(rGREAT)

PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"

### Importing data ###
#Importing myDiff
sign_DM <- read.csv(file.path(PATH, "Dataset_3", "2_BS-Seq", "sign_DM.csv"))
#General df
general_df <- sign_DM[, c("seqnames", "start", "deltaB")]
general_df$coord_key <- paste(general_df$seqnames, general_df$start, sep="_")

#For all methods the GRanges object for the annotation will be:
GR_data <- GRanges(
  seqnames = sign_DM$seqnames,
  ranges   = IRanges(start = sign_DM$start, end = sign_DM$start),   # CpG = 1 bp, start==end
  strand   = "*"                                   
)

#Ref for M1 and M5
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#Used in M4 and M5
annoData <- genes(txdb)



### Association CpG-Gene : Methods ###


### FUNCTIONS ###

write_output_file <- function(general_df, specific_df, Method_n) { 
  
  #Merging the files
  output_df <- merge(general_df, specific_df,
                     by = "coord_key",
                     sort = TRUE)
  
  file_name = sprintf(file.path(PATH, "Dataset_3", "3_Benchmark", "DM_sites_Met%d.csv"), Method_n)
  write.csv(output_df, file = file_name, row.names = FALSE)
  return(output_df)
}


##### Method 1: Nearest TSS (Genomation) #####
## annotateWithGeneParts() associate each site to the nearest TSS, with no hierarchy nor range limit

## Importing reference (RefSeq hg38)
gene.obj=readTranscriptFeatures(file.path(PATH, "references", "refseq.hg38.bed"))

#Annotation
diffCpGann <- annotateWithGeneParts(GR_data, gene.obj)

#Associated genes:
tss_df <- getAssociationWithTSS(diffCpGann)

#Cleaning names and converting them to gene_symbol:
refseq_ids <- tss_df$feature.name
refseq_ids <- sub("\\.[0-9]+$", "", refseq_ids)
#Converting
symbols <- mapIds(org.Hs.eg.db,
                  keys      = refseq_ids,
                  column    = "SYMBOL",
                  keytype   = "REFSEQ",
                  multiVals = "first")
tss_df$gene.symbol <- symbols

#Extracting the associated region_type (!! For this method this is not necessarily a region in the associated gene. It just says in what kind of region the CpG is)
members <- getMembers(diffCpGann)
region_type <- apply(members, 1, function(m) {
  if (m["prom"]   == 1) return("promoter")
  if (m["exon"]   == 1) return("exon")
  if (m["intron"] == 1) return("intron")
  return("other/intergenic")
})
tss_df$region_type <- region_type
tss_df <- tss_df %>% dplyr::rename(distanceToTSS = dist.to.feature, 
                                   SYMBOL = gene.symbol)
tss_df$distanceToTSS <- abs(tss_df$distanceToTSS)

DM_sites_M1 <- data.frame(tss_df[, c("region_type", "SYMBOL", "distanceToTSS")])
gr_for_M1 <- GR_data[tss_df$target.row]
DM_sites_M1$coord_key <- paste(seqnames(gr_for_M1), start(gr_for_M1), sep="_")

## Writing output file:
M1 <- write_output_file(general_df, DM_sites_M1, 1)




##### Methods 2 & 3: Nearest TSS with hierarchy(2), and proximal promoter(3) (CHIPseeker) #####

## With CHIPseeker annotation the final list will contain the results of method 2. And the 
## genes annotated as "Promoter" will be the results of method 3, with promoter defined as (-2000, 200).
#Now using CHIPseeker (hg38)
#Using txdb

#Annotation:
peakAnno <- annotatePeak(GR_data,
                         tssRegion = c(-2000, 200),   #This is a standard definition for Promoter. Can be arbitrarly changed
                         TxDb      = txdb,
                         annoDb    = "org.Hs.eg.db")

## Writing output file for Method 2:
peakAnno_df <- as.data.frame(peakAnno)
peakAnno_df <- peakAnno_df %>% dplyr::rename("region_type" = annotation)
peakAnno_df$coord_key <- paste(peakAnno_df$seqnames, peakAnno_df$start, sep="_")
peakAnno_df$distanceToTSS <- abs(peakAnno_df$distanceToTSS)

DM_sites_M2 <- data.frame(peakAnno_df[, c("region_type", "SYMBOL", "coord_key", "distanceToTSS")])
M2 <- write_output_file(general_df, DM_sites_M2, 2)

## Writing output file for Method 3:
DM_sites_M3 <- peakAnno_df[(peakAnno_df$region_type == 'Promoter (<=1kb)') | (peakAnno_df$region_type == 'Promoter (1-2kb)'), ]
DM_sites_M3 <- data.frame(DM_sites_M3[, c("region_type", "SYMBOL", "coord_key", "distanceToTSS")])
M3 <- write_output_file(general_df, DM_sites_M3, 3)



##### Method 4: CpG is in a range [-10 kb, + 10 kb] from TSS #####

#Defining the TSS
tss_points <- promoters(annoData, upstream = 0, downstream = 1)

#Association CpG-gene
hits <- findOverlaps(GR_data, tss_points, maxgap = 10000, select = "all")

anno_range_10kb_df <- data.frame(
  seqnames       = as.character(seqnames(GR_data)[queryHits(hits)]),
  start          = start(GR_data)[queryHits(hits)],
  feature        = names(tss_points)[subjectHits(hits)],
  tss_start       = start(tss_points)[subjectHits(hits)],
  tss_strand      = as.character(strand(tss_points)[subjectHits(hits)]),
  distanceToTSS   = distance(GR_data[queryHits(hits)], tss_points[subjectHits(hits)])
)

#Assigning SYMBOL and coord_keys
symbols_4 <- mapIds(org.Hs.eg.db,
                    keys      = anno_range_10kb_df$feature,
                    column    = "SYMBOL",
                    keytype   = "ENTREZID",
                    multiVals = "first")
anno_range_10kb_df$SYMBOL <- as.character(symbols_4)
anno_range_10kb_df$coord_key <- paste(anno_range_10kb_df$seqnames, anno_range_10kb_df$start, sep="_")
anno_range_10kb_final <- anno_range_10kb_df %>% filter(!is.na(SYMBOL))

## Writing output file:
DM_sites_M4 <- anno_range_10kb_final[, c("SYMBOL", "coord_key", "distanceToTSS")]
M4 <- write_output_file(general_df, DM_sites_M4, 4)



##### METHOD 5: rGREAT (-5 kb, 1 kb) plus extension until nearest gene up to 1 MB in both directions #####

#Defining the TSS
TSS_map <- extendTSS(annoData, gene_id_type = "ENTREZ", 
                     mode = "basalPlusExt", 
                     extend_from = "TSS", 
                     basal_upstream = 5000, 
                     basal_downstream = 1000, 
                     extension = 1000000)

#Now finding overlaps between CpG and the map
hits_M5 <- findOverlaps(GR_data, TSS_map)

#Matching TSS_map and tss_points, in order to use tss_maps to calculate distanceToTSS:
matched_idx <- match(names(TSS_map)[subjectHits(hits_M5)], names(tss_points))

great_df <- data.frame(seqnames      = as.character(seqnames(GR_data))[queryHits(hits_M5)],
                       start         = start(GR_data)[queryHits(hits_M5)],
                       feature       = names(TSS_map)[subjectHits(hits_M5)],  # gene_id ENTREZ
                       tss_start       = start(tss_points)[matched_idx],
                       tss_strand      = as.character(strand(tss_points))[matched_idx],
                       distanceToTSS   = distance(GR_data[queryHits(hits_M5)], tss_points[matched_idx])
)

#Now passing from ENTREZID to SYMBOLs
symbols_5 <- mapIds(org.Hs.eg.db,
                    keys = great_df$feature,
                    column = "SYMBOL",
                    keytype = "ENTREZID",
                    multiVals = "first")
                    
great_df$SYMBOL <- as.character(symbols_5)
great_df$coord_key <- paste(great_df$seqnames, great_df$start, sep = "_")
great_final <- great_df %>% filter(!is.na(SYMBOL))

## Final output for Method 5:
DM_sites_M5 <- great_final[, c("SYMBOL", "coord_key", "distanceToTSS")]
M5 <- write_output_file(general_df, DM_sites_M5, 5)

#####