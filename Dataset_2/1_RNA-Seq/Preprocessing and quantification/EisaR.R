library(eisaR)
library(Biostrings)
library(GenomicFeatures)
library(Rsamtools)

# Path to files downloaded
gtf_path    <- "references/gencode.v49.annotation.gtf"
genome_path <- "references/GRCh38.primary_assembly.genome.fa"

# 1. Generating annotation spliced + intron from GTF
grl <- getFeatureRanges(
  gtf = gtf_path,
  featureType = c("spliced", "intron"),
  intronType = "separate",
  flankLength = 90L,
  joinOverlappingIntrons = FALSE,
  verbose = TRUE
)

# 2. Exportinga GTF + t2g
exportToGtf(grl, filepath = "references/gencode.v49.spliced_intron.gtf")
df <- getTx2Gene(grl, filepath = "references/t2g_spliced_intron.tsv")

# 3. Indexing
indexFa(genome_path)

# 4. Extracting sequeces from FASTA
genome <- FaFile(genome_path)
seqs <- extractTranscriptSeqs(x = genome, transcripts = grl)

# 5. writing FASTA output (spliced + intron)
writeXStringSet(seqs, filepath = "references/merged_spliced_intron.fa")