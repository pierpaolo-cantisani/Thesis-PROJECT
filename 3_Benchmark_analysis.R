library(genomation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(rGREAT)


### Importing data ###
##Importing myDiff
myDiff25p_GR_hg38 <- readRDS("C:/Users/pierp/Desktop/THESIS PROJECT/BS-Seq/dataset/myDiff25p_GR_hg38.rds")
#General df
GR_df <- as.data.frame(myDiff25p_GR_hg38)
general_df <- GR_df[, c("seqnames", "start", "meth.diff")]
general_df$coord_key <- paste(general_df$seqnames, general_df$start, sep="_")

# In this pipeline the GRanges object for the annotation will be:
GR_data <- myDiff25p_GR_hg38



### --- Association CpG-Gene : Methods --- ###

### FUNCTIONS ###

write_output_file <- function(general_df, specific_df, Method_n) { 
  
  #Merging the files
  output_df <- merge(general_df, specific_df,
                     by = "coord_key",
                     sort = TRUE)
  
  file_name = sprintf("C:/Users/pierp/Desktop/THESIS PROJECT/benchmark/DM_sites_Met%d.csv", Method_n)
  write.csv(output_df, file = file_name, row.names = FALSE)
  return(output_df)
}

##### Method 1: Nearest TSS (Genomation) #####
## annotateWithGeneParts() associate each site to the nearest TSS, with no hierarchy nor range limit

## Importing reference (RefSeq hg38)
gene.obj=readTranscriptFeatures("C:/Users/pierp/Desktop/THESIS PROJECT/references/refseq.hg38.bed")
## Annotation
diffCpGann <- annotateWithGeneParts(GR_data, gene.obj)

## Associated genes:
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
#Now extracting the associated region_type (!! For this method this is not necessarily a region in the associated gene. It just says in what kind of region the CpG is)
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
DM_sites_M1 <- data.frame(tss_df[, c("region_type", "distanceToTSS", "SYMBOL")])
gr_for_M1 <- GR_data[tss_df$target.row]
DM_sites_M1$coord_key <- paste(seqnames(gr_for_M1), start(gr_for_M1), sep="_")
#Writing output file:
M1 <- write_output_file(general_df, DM_sites_M1, 1)




##### Methods 2 & 3: Nearest TSS with hierarchy(2), and proximal promoter(3) (CHIPseeker) #####

#Now using CHIPseeker (hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

##With CHIPseeker annotation the final list will contain the results of method 2. And the 
##genes annotated as "Promoter" will be the results of method 3, with promoter defined as (-2000, 200).
#Annotation:
peakAnno <- annotatePeak(GR_data,
                         tssRegion = c(-2000, 200),   #This is the standard definition for Promoter. Can be arbitrarly changed
                         TxDb      = txdb,
                         annoDb    = "org.Hs.eg.db")

#Writing output file for Method 2:
peakAnno_df <- as.data.frame(peakAnno)
peakAnno_df <- peakAnno_df %>% dplyr::rename("region_type" = annotation)
peakAnno_df$coord_key <- paste(peakAnno_df$seqnames, peakAnno_df$start, sep="_")
DM_sites_M2 <- data.frame(peakAnno_df[, c("region_type", "distanceToTSS", "SYMBOL", "coord_key")])
M2 <- write_output_file(general_df, DM_sites_M2, 2)

#Writing output file for Method 3:
DM_sites_M3 <- peakAnno_df[(peakAnno_df$region_type == 'Promoter (<=1kb)') | (peakAnno_df$region_type == 'Promoter (1-2kb)'), ]
DM_sites_M3 <- data.frame(DM_sites_M3[, c("region_type", "distanceToTSS", "SYMBOL", "coord_key")])
M3 <- write_output_file(general_df, DM_sites_M3, 3)



##### Method 4: CpG is in a range [-10 kb, + 10 kb] from TSS #####

annoData <- genes(txdb)
tss_points <- promoters(annoData, upstream = 0, downstream = 1)

#Association CpG-gene
hits <- findOverlaps(GR_data, tss_points, maxgap = 10000, select = "all")

anno_range_10kb_df <- data.frame(
  seqnames       = as.character(seqnames(GR_data)[queryHits(hits)]),
  start          = start(GR_data)[queryHits(hits)],
  distanceToTSS  = start(GR_data)[queryHits(hits)] - start(tss_points)[subjectHits(hits)],
  feature        = names(tss_points)[subjectHits(hits)]
)

#Assigning SYMBOL and coord_keys
symbols_4 <- mapIds(org.Hs.eg.db,
                    keys      = anno_range_10kb_df$feature,
                    column    = "SYMBOL",
                    keytype   = "ENTREZID",
                    multiVals = "first")
anno_range_10kb_df$SYMBOL <- as.character(symbols_4)
anno_range_10kb_df$coord_key <- paste(anno_range_10kb_df$seqnames, anno_range_10kb_df$start, sep="_")
anno_range_10kb_final <- anno_range_10kb_df %>% filter(!is.na(feature))

#Writing output file:
DM_sites_M4 <- anno_range_10kb_final[, c("distanceToTSS", "SYMBOL", "coord_key")]
M4 <- write_output_file(general_df, DM_sites_M4, 4)



#####
##### METHOD 5: rGREAT (-5 kb, 1 kb) plus extension until nearest gene up to 1 MB in both directions #####
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_GR <- genes(txdb)

TSS_map <- extendTSS(genes_GR, gene_id_type = "ENTREZ", 
                     mode = "basalPlusExt", 
                     extend_from = c("TSS", "gene"), 
                     basal_upstream = 5000, 
                     basal_downstream = 1000, 
                     extension = 1000000)

#Now finding overlaps between CpG and the map
hits <- findOverlaps(GR_data, TSS_map)

great_df <- data.frame(seqnames      = as.character(seqnames(GR_data))[queryHits(hits)],
                       start         = start(GR_data)[queryHits(hits)],
                       feature       = names(TSS_map)[subjectHits(hits)]  # gene_id ENTREZ
)

#Now passing from ENTREZID to SYMBOLs
symbols_5 <- mapIds(org.Hs.eg.db,
                    keys = great_df$feature,
                    column = "SYMBOL",
                    keytype = "ENTREZID",
                    multiVals = "first")
                    
great_df$SYMBOL <- as.character(symbols_5)
great_df$coord_key <- paste(great_df$seqnames, great_df$start, sep = "_")
great_final <- great_df %>% filter(!is.na(feature))

#Final output for Method 5:
DM_sites_M5 <- great_final[, c("SYMBOL", "coord_key")]
M5 <- write_output_file(general_df, DM_sites_M5, 5)


#### PROBLEMI
#- Per quanto riguarda il discorso della regione genomica, il metodo 1 ha regione genomica associata
#  non affidabile (magari cade in introne del gene A, ma il TSS più vicino è quello del gene B). Anche per
#  il metodo 4 potrebe essere così.
#  Una possibile soluzione è considerare concettualmente la regione non come regione del gene associato, bensì
#  come regione in generale. Quindi è una considerazione su come il tipo di regione influenza il gene associato, non su come
#  il tipo di regione del gene ASSOCIATO influenza il gene associato


#####