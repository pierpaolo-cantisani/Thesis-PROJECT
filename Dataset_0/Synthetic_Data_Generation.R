#!/usr/bin/env Rscript
# ============================================================================
# Round-trip synthetic data generator (OPTION 1, reference method = M2/ChIPseeker)
# RNA-seq (.sf, keyed by SYMBOL) + WGBS (chr pos M coverage, single-strand)
#
# M2 (ChIPseeker annotatePeak, GENCODE TxDb + org.Hs.eg.db) is run INSIDE this
# script on the generated CpG coordinates, and the DE signal is anchored to the
# gene M2 associates each CpG to. Therefore M2's own output IS the ground truth:
# running M2 in the pipeline returns 100% of the prepared associations, and any
# other method's associations are deviations measured against M2.
#
# Interface note: there are NO planted category fractions. In option 1 the
# association is DISCOVERED by M2 (which assigns every site to a nearest gene,
# far ones labelled "Distal Intergenic"), so you set:
#   p_DM      : fraction of CpGs that are DM
#   p_gene_DE : fraction of M2-associated genes that are DE
# The both/dm_only/de_only/neither split and the region-type distribution are
# EMERGENT, measured from ground truth and written to design_summary.txt --
# that measured value is the round-trip target.
#
# ChIPseeker region hierarchy: Promoter > 5'UTR > 3'UTR > Exon > Intron >
#   Downstream > Distal Intergenic. Placement below mirrors these definitions
#   (promoter within tssRegion, exon = CDS, intron, UTR, distal = far), but the
#   authoritative region_type is whatever M2 assigns.
#
# Sign convention: meth.diff = case - ctrl (>0 hyper); log2FC = log2(case/ctrl)
#   Q1(+,+) up&hyper   Q2(+,-) up&hypo   Q3(-,-) down&hypo   Q4(-,+) down&hyper
#
# Outputs (directories must already exist):
#   outdir_rna/rnaseq/<sample>.sf   quant.sf (Name = SYMBOL)
#   outdir_wgbs/wgbs/<sample>.txt   chr pos M coverage  (single-strand, M+U = coverage)
#   outdir/samples.tsv / ground_truth.tsv / design_summary.txt
# ============================================================================

library(rtracklayer)     # import()
library(GenomicRanges)
library(GenomicFeatures) # makeTxDbFromGFF  (newer Bioconductor: also needs txdbmaker)
library(ChIPseeker)      # annotatePeak
library(org.Hs.eg.db)    # annoDb for SYMBOL

## =============================== PARAMETERS ================================
set.seed(1)

PATH        <- "C:/Users/pierp/Desktop/Thesis PROJECT"   # <-- project root
gtf_path    <- file.path(PATH, "references", "gencode.v49.annotation.gtf.gz")
outdir_rna  <- file.path(PATH, "Dataset_0", "1_RNA-Seq")
outdir_wgbs <- file.path(PATH, "Dataset_0", "2_BS-Seq")
outdir      <- file.path(PATH, "Dataset_0", "Ground Truth")

n_genes     <- 12000           # GENCODE genes used only as placement scaffolding
n_cpg_total <- 15000           # total CpG positions to emit

p_DM       <- 0.50             # fraction of CpGs that are DM
p_gene_DE  <- 0.50             # fraction of M2-associated genes that are DE

tss_region <- c(-2000, 200)    # ChIPseeker promoter window (matches your pipeline)

# --- quadrant split of the "both" set (drives direction balance; MUST sum 1) -
quadrant_props <- c(Q1 = 0.10, Q2 = 0.40, Q3 = 0.10, Q4 = 0.40)

# --- experimental design (PAIRED case-control) ------------------------------
n_pairs      <- 50
subj_sd_rna  <- 0.4
subj_sd_meth <- 0.7

# --- effect sizes (LARGE => power ~1 => recovery ~ identity) ----------------
lfc_mag    <- 2.0
deltaM_mag <- 0.40

# --- RNA-seq (negative binomial) -------------------------------------------
base_mean_log   <- 6
nb_dispersion   <- 0.1
tx_length_range <- c(800, 6000)

# --- WGBS -------------------------------------------------------------------
mean_coverage  <- 20
cov_dispersion <- 0.2
base_beta      <- 0.5
bb_conc        <- 30

# --- CpG placement: BY DISTANCE FROM TSS (mimics the real process) ----------
# No region type is chosen. Each non-distal CpG is dropped at TSS +/- an offset
# drawn log-uniformly from tss_dist_range (strand-aware); M2 then labels whatever
# region it happened to fall in -- exactly as with real CpGs. Tune the three
# knobs until design_summary's promoter fraction matches your data (~20-30%):
# smaller tss_dist_range / lower p_downstream => more promoters.
frac_distal    <- 0.10           # fraction dropped far from any gene (-> Distal Intergenic)
tss_dist_range <- c(100, 1e5)    # |offset| from TSS in bp, sampled log-uniformly
p_downstream   <- 0.80           # P(offset is downstream of TSS, i.e. into the gene body)
## ===========================================================================

stopifnot(abs(sum(quadrant_props) - 1) < 1e-9)

## ------------------------- 1. genes + features from GTF (placement) --------
gtf   <- import(gtf_path)
genes <- gtf[gtf$type == "gene"]
genes <- keepStandardChromosomes(genes, pruning.mode = "coarse")
genes <- genes[!duplicated(genes$gene_id)]
if (length(genes) < n_genes) stop("GTF has fewer genes than n_genes")
genes <- genes[sample(length(genes), n_genes)]

tss <- ifelse(strand(genes) == "-", end(genes), start(genes))
gene_tab <- data.frame(
  gene_id = genes$gene_id, chr = as.character(seqnames(genes)),
  gstart = start(genes), gend = end(genes),
  strand = as.character(strand(genes)), tss = tss, stringsAsFactors = FALSE)

## ------------------------- 2. place CpG coordinates (by TSS distance) -------
n <- n_cpg_total
cpg <- data.frame(cpg_id = sprintf("cpg%05d", seq_len(n)), stringsAsFactors = FALSE)
is_distal <- runif(n) < frac_distal

chr <- character(n); pos <- integer(n)

# non-distal: random host gene, position = TSS +/- log-uniform offset (strand-aware).
# We do NOT pick a region type; M2 will label wherever the CpG landed.
k  <- sum(!is_distal)
gs <- gene_tab[sample(nrow(gene_tab), k, replace = TRUE), ]
d  <- round(10^runif(k, log10(tss_dist_range[1]), log10(tss_dist_range[2])))
downstream <- runif(k) < p_downstream
# genomic sign: downstream = +d on '+' strand, -d on '-' strand (upstream reversed)
sgn <- ifelse(gs$strand == "-", -1L, 1L) * ifelse(downstream, 1L, -1L)
chr[!is_distal] <- gs$chr
pos[!is_distal] <- pmax(1L, as.integer(gs$tss + sgn * d))

# distal: random far positions -> M2 labels them "Distal Intergenic"
chr_pool <- unique(gene_tab$chr)
chr[is_distal] <- sample(chr_pool, sum(is_distal), replace = TRUE)
pos[is_distal] <- as.integer(runif(sum(is_distal), 1e4, 2e8))

cpg$chr <- chr; cpg$pos <- as.integer(pos)
cpg$strand <- "+"    # single-strand projection (destranded CpGs)

## ------------------------- 3. M2: ChIPseeker annotatePeak ------------------
txdb <- makeTxDbFromGFF(gtf_path, format = "gtf")
GR_data <- GRanges(cpg$chr, IRanges(cpg$pos, cpg$pos), strand = cpg$strand)
peakAnno <- annotatePeak(GR_data, tssRegion = tss_region, TxDb = txdb, annoDb = "org.Hs.eg.db")
pa <- as.data.frame(peakAnno)
pa$coord_key <- paste(pa$seqnames, pa$start, sep = "_")
pa$distanceToTSS <- abs(pa$distanceToTSS)

cpg$coord_key <- paste(cpg$chr, cpg$pos, sep = "_")
idx <- match(cpg$coord_key, pa$coord_key)          # annotatePeak may reorder -> map by coord
cpg$m2_symbol <- pa$SYMBOL[idx]
cpg$m2_region <- pa$annotation[idx]                # raw ChIPseeker label = ground-truth region_type
cpg$m2_dist   <- pa$distanceToTSS[idx]


## ------------------------- 4. DE anchored to M2 gene + directions ----------
cpg$is_DM <- runif(n) < p_DM
u_sym <- unique(cpg$m2_symbol[!is.na(cpg$m2_symbol)])
p_up  <- as.numeric(quadrant_props["Q1"] + quadrant_props["Q2"])
gene_de <- data.frame(SYMBOL = u_sym, stringsAsFactors = FALSE)
gene_de$is_DE <- runif(nrow(gene_de)) < p_gene_DE
gene_de$dir   <- 0L
de_rows <- which(gene_de$is_DE)
gene_de$dir[de_rows] <- ifelse(runif(length(de_rows)) < p_up, 1L, -1L)
gene_de$log2FC <- gene_de$dir * lfc_mag

mi <- match(cpg$m2_symbol, gene_de$SYMBOL)
cpg$host_isDE <- ifelse(is.na(mi), FALSE, gene_de$is_DE[mi])
cpg$gene_dir  <- ifelse(is.na(mi), 0L,    gene_de$dir[mi])

cpg$deltaM <- 0; cpg$quadrant <- NA_character_
p_hyper_up   <- as.numeric(quadrant_props["Q1"] / (quadrant_props["Q1"] + quadrant_props["Q2"]))
p_hyper_down <- as.numeric(quadrant_props["Q4"] / (quadrant_props["Q3"] + quadrant_props["Q4"]))
both_i <- which(cpg$is_DM & cpg$host_isDE)
up_i   <- both_i[cpg$gene_dir[both_i] ==  1L]
down_i <- both_i[cpg$gene_dir[both_i] == -1L]
hu <- runif(length(up_i))   < p_hyper_up
hd <- runif(length(down_i)) < p_hyper_down
cpg$deltaM[up_i]   <- ifelse(hu, deltaM_mag, -deltaM_mag); cpg$quadrant[up_i]   <- ifelse(hu, "Q1", "Q2")
cpg$deltaM[down_i] <- ifelse(hd, deltaM_mag, -deltaM_mag); cpg$quadrant[down_i] <- ifelse(hd, "Q4", "Q3")
dmonly_i <- which(cpg$is_DM & !cpg$host_isDE)
cpg$deltaM[dmonly_i] <- sample(c(-1, 1), length(dmonly_i), TRUE) * deltaM_mag

cpg$category <- with(cpg,
                     ifelse(is_DM & host_isDE, "both",
                            ifelse(is_DM & !host_isDE, "dm_only",
                                   ifelse(!is_DM & host_isDE, "de_only", "neither"))))

## ------------------------- 5. per-SYMBOL expression ------------------------
gene_de$length  <- sample(tx_length_range[1]:tx_length_range[2], nrow(gene_de), TRUE)
gene_de$mu_ctrl <- exp(rnorm(nrow(gene_de), base_mean_log, 1))
gene_de$mu_case <- gene_de$mu_ctrl * 2^gene_de$log2FC

## ------------------------- 6. simulate + write RNA (.sf, keyed by SYMBOL) --
write_sf <- function(counts, len, ids, path) {
  eff <- pmax(len - 200, 1); rate <- counts / eff
  tpm <- if (sum(rate) > 0) rate / sum(rate) * 1e6 else rep(0, length(rate))
  write.table(data.frame(Name = ids, Length = len, EffectiveLength = eff,
                         TPM = round(tpm, 4), NumReads = round(counts, 3)),
              path, sep = "\t", quote = FALSE, row.names = FALSE)
}
for (p in seq_len(n_pairs)) {
  re <- rnorm(nrow(gene_de), 0, subj_sd_rna)                 # shared subject baseline (paired)
  write_sf(rnbinom(nrow(gene_de), mu = gene_de$mu_ctrl * exp(re), size = 1 / nb_dispersion),
           gene_de$length, gene_de$SYMBOL, file.path(outdir_rna, "rnaseq", sprintf("ctrl_%d.sf", p)))
  write_sf(rnbinom(nrow(gene_de), mu = gene_de$mu_case * exp(re), size = 1 / nb_dispersion),
           gene_de$length, gene_de$SYMBOL, file.path(outdir_rna, "rnaseq", sprintf("case_%d.sf", p)))
}

## ------------------------- 7. simulate + write WGBS (chr pos M coverage) ---
cpg$beta_ctrl <- base_beta
cpg$beta_case <- pmin(pmax(base_beta + cpg$deltaM, 0.01), 0.99)
rbetabinom <- function(k, size, mu, conc) rbinom(k, size, rbeta(k, mu * conc, (1 - mu) * conc))
logit <- function(x) log(x / (1 - x)); invlogit <- function(x) 1 / (1 + exp(-x))
clampb <- function(x) pmin(pmax(x, 0.01), 0.99)
ord <- order(cpg$chr, cpg$pos)
write_wgbs <- function(beta, cond, p) {
  cov <- rnbinom(n, mu = mean_coverage, size = 1 / cov_dispersion) + 1   # coverage = M + U
  M   <- rbetabinom(n, cov, beta, bb_conc)                              # methylated count
  write.table(data.frame(chr = cpg$chr, pos = cpg$pos, M = M, coverage = cov)[ord, ],
              file.path(outdir_wgbs, "wgbs", sprintf("%s_%d.txt", cond, p)),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
for (p in seq_len(n_pairs)) {
  re <- rnorm(n, 0, subj_sd_meth)
  bc <- clampb(invlogit(logit(clampb(cpg$beta_ctrl)) + re))
  write_wgbs(bc, "ctrl", p); write_wgbs(clampb(bc + cpg$deltaM), "case", p)
}

## ------------------------- 8. ground truth + summary -----------------------
write.table(cpg, file.path(outdir, "ground_truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
samples <- data.frame(
  sample    = c(sprintf("ctrl_%d", seq_len(n_pairs)), sprintf("case_%d", seq_len(n_pairs))),
  condition = rep(c("ctrl", "case"), each = n_pairs),
  pair      = factor(rep(seq_len(n_pairs), times = 2)), stringsAsFactors = FALSE)
write.table(samples, file.path(outdir, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

coarse <- function(x) ifelse(grepl("^Promoter", x), "promoter",
                             ifelse(grepl("5' UTR", x),   "5UTR",
                                    ifelse(grepl("3' UTR", x),   "3UTR",
                                           ifelse(grepl("^Exon", x),    "exon",
                                                  ifelse(grepl("^Intron", x),  "intron",
                                                         ifelse(grepl("^Downstream", x), "downstream",
                                                                ifelse(grepl("Intergenic", x),  "intergenic", "other")))))))
sm <- c(
  sprintf("CpGs: %d | M2 genes: %d | DE genes: %d", n, nrow(gene_de), sum(gene_de$is_DE)),
  "--- REALISED marginals (round-trip target; emergent under M2) ---",
  sprintf("DM (any)   : %.3f", mean(cpg$is_DM)),
  sprintf("both DE&DM : %.3f", mean(cpg$category == "both")),
  sprintf("dm_only    : %.3f", mean(cpg$category == "dm_only")),
  sprintf("de_only    : %.3f", mean(cpg$category == "de_only")),
  "--- quadrant split of the BOTH set (realised) ---",
  capture.output(print(round(prop.table(table(cpg$quadrant)), 3))),
  "--- M2 region_type of DM CpGs (coarse) ---",
  capture.output(print(round(prop.table(table(coarse(cpg$m2_region[cpg$is_DM]))), 3)))
)
writeLines(sm, file.path(outdir, "design_summary.txt"))
cat(paste(sm, collapse = "\n"), "\n")