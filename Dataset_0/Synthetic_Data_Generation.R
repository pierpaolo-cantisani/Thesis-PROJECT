#!/usr/bin/env Rscript
# ============================================================================
# Round-trip synthetic data generator: RNA-seq (Salmon .sf) + WGBS (chr/pos/M/coverage)
#
# Goal: emit data whose MARGINAL properties are known by construction, so you can
# check that they survive the pipeline (a technical / round-trip control, NOT a
# precision-recall test of association methods):
#   - fraction of CpGs associated to a gene
#   - fraction that are BOTH DM (CpG) and DE (host)
#   - the quadrant split of that "both" set in the log2FC (x) vs meth.diff (y) plane
#   - the region-type distribution of the CpGs
#
# Design: ONE CpG per associated gene (each categorised associated CpG gets a
# UNIQUE host gene), so a gene's DE status/direction is unambiguous and every
# marginal is exact. Consequence: the CpG-per-gene distribution is degenerate
# (all hosting genes have exactly 1 CpG) -- for a pipeline step that stratifies
# by #CpG-per-gene, the expected result is "only the 1-CpG bin is populated,
# all higher bins empty". That degenerate expectation is itself the control.
#
# Sign convention: meth.diff = case - ctrl (>0 hyper); log2FC = log2(case/ctrl)
#   Q1(+,+) up&hyper   Q2(+,-) up&hypo   Q3(-,-) down&hypo   Q4(-,+) down&hyper
#
# Outputs (directories must already exist -- creation is intentionally NOT done here):
#   outdir_rna/rnaseq/<sample>.sf   Salmon-style quant.sf (Name Length EffectiveLength TPM NumReads)
#   outdir_wgbs/wgbs/<sample>.txt   tab-separated: chr pos strand M coverage  (M = methylated, coverage = M+U)
#   outdir/samples.tsv              design table: sample condition pair  (for the PAIRED model)
#   outdir/ground_truth.tsv         per-CpG truth table
#   outdir/design_summary.txt       planted marginals = your round-trip target
# ============================================================================

library(rtracklayer)     # import()
library(GenomicRanges)

## =============================== PARAMETERS ================================
set.seed(1)

PATH        <- "C:/Users/pierp/Desktop/Thesis PROJECT"   # <-- project root
gtf_path    <- file.path(PATH, "references", "gencode.v49.annotation.gtf.gz")
outdir_rna  <- file.path(PATH, "Dataset_0", "1_RNA-Seq")
outdir_wgbs <- file.path(PATH, "Dataset_0", "2_BS-Seq")
outdir      <- file.path(PATH, "Dataset_0", "Ground Truth")

n_genes     <- 12000           # MUST be >= n_cpg_total * (1 - frac["background"]); see check below
n_cpg_total <- 15000           # total CpG positions to emit

# --- category fractions over ALL CpGs (MUST sum to 1) -----------------------
#   background : not associated to any gene (intergenic noise)
#   both       : CpG is DM AND its host gene is DE            <-- the "DE DM intersection"
#   dm_only    : CpG is DM, host gene NOT DE
#   de_only    : host gene DE, this CpG NOT DM
#   neither    : associated, no perturbation
frac <- c(background = 0.30,
          both       = 0.40,
          dm_only    = 0.10,
          de_only    = 0.10,
          neither    = 0.10)

# --- quadrant split of the "both" set (MUST sum to 1) -----------------------
quadrant_props <- c(Q1 = 0.10, Q2 = 0.40, Q3 = 0.10, Q4 = 0.40)

# --- experimental design (PAIRED case-control) ------------------------------
# n_pairs individuals, each measured in BOTH conditions => 2 * n_pairs samples.
# ctrl_p and case_p share a per-subject baseline (random intercept).
n_pairs      <- 50
subj_sd_rna  <- 0.4   # between-subject SD, log scale (RNA). 0 => unpaired-equivalent
subj_sd_meth <- 0.7   # between-subject SD, logit scale (methylation)

# --- effect sizes (LARGE => power ~1 => recovery ~ identity) ----------------
lfc_mag    <- 2.0              # |log2FC| for DE genes
deltaM_mag <- 0.40             # |delta beta| for DM CpGs (beta scale)

# --- RNA-seq (negative binomial) -------------------------------------------
base_mean_log   <- 6          # log baseline expression mean
nb_dispersion   <- 0.1
tx_length_range <- c(800, 6000)

# --- WGBS -------------------------------------------------------------------
mean_coverage  <- 20
cov_dispersion <- 0.2         # NB dispersion for coverage
base_beta      <- 0.5         # baseline methylation
bb_conc        <- 30          # beta-binomial concentration (higher = tighter)

# --- CpG placement relative to host gene (sets the region_type stat) --------
# promoter/downstream are derived from TSS/gene-end; exon/intron/utr are read
# from the GTF sub-gene features (see section 1b). If a gene lacks a chosen
# feature (e.g. a non-coding gene has no UTR), placement falls back to promoter
# and the ACTUAL region_type is what gets recorded in ground_truth.
region_props   <- c(promoter = 0.30, utr = 0.10, exon = 0.25,
                    intron = 0.25, downstream = 0.10)   # MUST sum to 1
promoter_win   <- c(-2000, 200)   # rel. to TSS, strand-aware
downstream_win <- c(1, 5000)      # past gene end
## ===========================================================================

stopifnot(abs(sum(frac) - 1)           < 1e-9)
stopifnot(abs(sum(quadrant_props) - 1) < 1e-9)
stopifnot(abs(sum(region_props) - 1)   < 1e-9)

## ------------------------- 1. genes + TSS from GTF -------------------------
gtf   <- import(gtf_path)
genes <- gtf[gtf$type == "gene"]
genes <- keepStandardChromosomes(genes, pruning.mode = "coarse")
genes <- genes[!duplicated(genes$gene_id)]
if (length(genes) < n_genes) stop("GTF has fewer genes than n_genes")
genes <- genes[sample(length(genes), n_genes)]

tss <- ifelse(strand(genes) == "-", end(genes), start(genes))
gene_tab <- data.frame(
  gene_id = genes$gene_id,
  chr     = as.character(seqnames(genes)),
  gstart  = start(genes), gend = end(genes),
  strand  = as.character(strand(genes)),
  tss     = tss,
  stringsAsFactors = FALSE
)

## ------------------------- 1b. sub-gene features ---------------------------
# Read exons and UTRs straight from the GTF (GENCODE labels them "exon"/"UTR",
# and every feature line carries gene_id). Introns are not a GTF feature type,
# so we compute them as (gene span) minus (exonic footprint), per gene.
# NOTE: GENCODE uses a generic "UTR" type. On Ensembl GTFs it is instead
#       "five_prime_utr"/"three_prime_utr" -- adjust the type filter if needed.
gids   <- gene_tab$gene_id
ex_gr  <- gtf[gtf$type == "exon" & gtf$gene_id %in% gids]
utr_gr <- gtf[gtf$type == "UTR"  & gtf$gene_id %in% gids]

exons_by <- reduce(split(ex_gr,  factor(ex_gr$gene_id,  levels = gids)))   # merged exonic footprint / gene
utr_by   <- reduce(split(utr_gr, factor(utr_gr$gene_id, levels = gids)))   # merged UTR / gene

genes_gr <- GRanges(gene_tab$chr, IRanges(gene_tab$gstart, gene_tab$gend),
                    strand = gene_tab$strand)
names(genes_gr) <- gids
introns_by <- psetdiff(genes_gr, exons_by[gids])   # gene span minus exons = introns / gene

# helper: draw a position uniformly over the total length of a gene's feature set
sample_pos_in <- function(gr) {
  if (length(gr) == 0) return(NA_integer_)
  i <- sample(seq_along(gr), 1, prob = width(gr))          # width-weighted interval
  sample(start(gr)[i]:end(gr)[i], 1)
}

## ------------------------- 2. assign CpG categories ------------------------
cats     <- names(frac)
n_by_cat <- round(frac * n_cpg_total)
n_by_cat[length(n_by_cat)] <- n_cpg_total - sum(n_by_cat[-length(n_by_cat)])  # fix rounding
category <- rep(cats, times = n_by_cat)
category <- sample(category)                                   # shuffle
cpg <- data.frame(cpg_id = sprintf("cpg%05d", seq_len(n_cpg_total)),
                  category = category, stringsAsFactors = FALSE)

is_assoc <- cpg$category != "background"
n_assoc  <- sum(is_assoc)
if (n_assoc > n_genes)
  stop(sprintf("Not enough genes: need >= %d associated genes but n_genes = %d. Raise n_genes or the background fraction.",
               n_assoc, n_genes))

## ------------------------- 3. place associated CpGs ------------------------
host_idx <- integer(n_cpg_total)
host_idx[is_assoc] <- sample(nrow(gene_tab), n_assoc)          # UNIQUE host per associated CpG (1 CpG/gene)

region_type <- rep(NA_character_, n_cpg_total)
chr <- rep(NA_character_, n_cpg_total); pos <- rep(NA_integer_, n_cpg_total)
strand_v <- rep("*", n_cpg_total)

place_one <- function(g) {
  rt <- sample(names(region_props), 1, prob = region_props)
  p <- switch(rt,
              promoter = { off <- sample(promoter_win[1]:promoter_win[2], 1)
              if (g$strand == "-") g$tss - off else g$tss + off },
              exon     = sample_pos_in(exons_by[[g$gene_id]]),
              intron   = sample_pos_in(introns_by[[g$gene_id]]),
              utr      = sample_pos_in(utr_by[[g$gene_id]]),
              downstream = { off <- sample(downstream_win[1]:downstream_win[2], 1)
              if (g$strand == "-") g$gstart - off else g$gend + off })
  if (is.na(p)) {   # fallback if the chosen feature is absent for this gene
    rt  <- "promoter"
    off <- sample(promoter_win[1]:promoter_win[2], 1)
    p   <- if (g$strand == "-") g$tss - off else g$tss + off
  }
  list(rt = rt, chr = g$chr, pos = max(1L, as.integer(p)), strand = g$strand)
}
for (i in which(is_assoc)) {
  pl <- place_one(gene_tab[host_idx[i], ])
  region_type[i] <- pl$rt; chr[i] <- pl$chr; pos[i] <- pl$pos; strand_v[i] <- pl$strand
}
# background CpGs: random intergenic-ish positions (see caveat in the notes)
chr_pool <- unique(gene_tab$chr)
for (i in which(!is_assoc)) {
  chr[i] <- sample(chr_pool, 1)
  pos[i] <- sample(1e4:2e8, 1)
  strand_v[i] <- sample(c("+", "-"), 1)
}
cpg$region_type <- region_type; cpg$chr <- chr; cpg$pos <- pos; cpg$strand <- strand_v
cpg$host_gene <- ifelse(is_assoc, gene_tab$gene_id[host_idx], NA)

## ------------------------- 4. assign directions ----------------------------
# quadrant -> (sign log2FC, sign meth.diff)
q_signs <- list(Q1 = c(lfc =  1, dm =  1), Q2 = c(lfc =  1, dm = -1),
                Q3 = c(lfc = -1, dm = -1), Q4 = c(lfc = -1, dm =  1))

cpg$is_DM     <- cpg$category %in% c("both", "dm_only")
cpg$host_isDE <- cpg$category %in% c("both", "de_only")
cpg$quadrant  <- NA_character_
cpg$log2FC    <- 0
cpg$deltaM    <- 0

# BOTH: draw quadrant, set coupled signs
both_i <- which(cpg$category == "both")
qn     <- round(quadrant_props * length(both_i))
qn[length(qn)] <- length(both_i) - sum(qn[-length(qn)])
qlab   <- sample(rep(names(quadrant_props), times = qn))
cpg$quadrant[both_i] <- qlab
cpg$log2FC[both_i]   <- vapply(qlab, function(q) q_signs[[q]]["lfc"], 0) * lfc_mag
cpg$deltaM[both_i]   <- vapply(qlab, function(q) q_signs[[q]]["dm"],  0) * deltaM_mag

# dm_only / de_only: independent random signs
dm_i <- which(cpg$category == "dm_only")
cpg$deltaM[dm_i] <- sample(c(-1, 1), length(dm_i), TRUE) * deltaM_mag
de_i <- which(cpg$category == "de_only")
cpg$log2FC[de_i] <- sample(c(-1, 1), length(de_i), TRUE) * lfc_mag

## ------------------------- 5. per-gene expression --------------------------
# collapse to gene level: a gene's log2FC comes from its (unique) primary CpG
gene_tab$log2FC <- 0
prim <- !is.na(cpg$host_gene) & cpg$host_isDE
m <- match(cpg$host_gene[prim], gene_tab$gene_id)
gene_tab$log2FC[m] <- cpg$log2FC[prim]

gene_tab$length   <- sample(tx_length_range[1]:tx_length_range[2], nrow(gene_tab), TRUE)
gene_tab$mu_ctrl  <- exp(rnorm(nrow(gene_tab), base_mean_log, 1))
gene_tab$mu_case  <- gene_tab$mu_ctrl * 2^gene_tab$log2FC

## ------------------------- 6. simulate + write RNA (.sf) -------------------
write_sf <- function(counts, len, ids, path) {
  eff <- pmax(len - 200, 1)
  rate <- counts / eff
  tpm  <- if (sum(rate) > 0) rate / sum(rate) * 1e6 else rep(0, length(rate))
  df <- data.frame(Name = ids, Length = len, EffectiveLength = eff,
                   TPM = round(tpm, 4), NumReads = round(counts, 3))
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
for (p in seq_len(n_pairs)) {
  # per-gene subject baseline, SHARED by ctrl_p and case_p (this is what makes it paired)
  re      <- rnorm(nrow(gene_tab), 0, subj_sd_rna)
  mu_ctrl <- gene_tab$mu_ctrl * exp(re)
  mu_case <- gene_tab$mu_case * exp(re)   # treatment (2^log2FC) sits on top of the shared baseline
  write_sf(rnbinom(nrow(gene_tab), mu = mu_ctrl, size = 1 / nb_dispersion),
           gene_tab$length, gene_tab$gene_id,
           file.path(outdir_rna, "rnaseq", sprintf("ctrl_%d.sf", p)))
  write_sf(rnbinom(nrow(gene_tab), mu = mu_case, size = 1 / nb_dispersion),
           gene_tab$length, gene_tab$gene_id,
           file.path(outdir_rna, "rnaseq", sprintf("case_%d.sf", p)))
}

## ------------------------- 7. simulate + write WGBS ------------------------
cpg$beta_ctrl <- base_beta
cpg$beta_case <- pmin(pmax(base_beta + cpg$deltaM, 0.01), 0.99)

rbetabinom <- function(n, size, mu, conc) {
  p <- rbeta(n, mu * conc, (1 - mu) * conc)
  rbinom(n, size, p)
}
logit    <- function(x) log(x / (1 - x))
invlogit <- function(x) 1 / (1 + exp(-x))
clampb   <- function(x) pmin(pmax(x, 0.01), 0.99)
ord <- order(cpg$chr, cpg$pos)

write_wgbs <- function(beta, cond, p) {
  cov <- rnbinom(n_cpg_total, mu = mean_coverage, size = 1 / cov_dispersion) + 1  # coverage = M + U
  M   <- rbetabinom(n_cpg_total, cov, beta, bb_conc)                              # methylated count ("M")
  # last two columns match your files: M (methylated) and coverage (M+U).
  # If a downstream branch needs M-values: Mvalue = log2((M + 1) / (cov - M + 1))
  df <- data.frame(chr = cpg$chr, pos = cpg$pos, strand = cpg$strand,
                   M = M, coverage = cov)[ord, ]
  write.table(df, file.path(outdir_wgbs, "wgbs", sprintf("%s_%d.txt", cond, p)),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

for (p in seq_len(n_pairs)) {
  # per-CpG subject baseline on the logit scale, SHARED by ctrl_p and case_p
  re          <- rnorm(n_cpg_total, 0, subj_sd_meth)
  beta_ctrl_p <- clampb(invlogit(logit(clampb(cpg$beta_ctrl)) + re))
  beta_case_p <- clampb(beta_ctrl_p + cpg$deltaM)   # delta beta treatment on top of the shared baseline
  write_wgbs(beta_ctrl_p, "ctrl", p)
  write_wgbs(beta_case_p, "case", p)
}

## ------------------------- 8. ground truth + summary -----------------------
write.table(cpg, file.path(outdir, "ground_truth.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# design table for the PAIRED model (e.g. ~ pair + condition)
samples <- data.frame(
  sample    = c(sprintf("ctrl_%d", seq_len(n_pairs)), sprintf("case_%d", seq_len(n_pairs))),
  condition = rep(c("ctrl", "case"), each = n_pairs),
  pair      = factor(rep(seq_len(n_pairs), times = 2)),
  stringsAsFactors = FALSE
)
write.table(samples, file.path(outdir, "samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sm <- c(
  sprintf("Total CpGs: %d   |   Genes: %d   |   1 CpG per hosting gene", n_cpg_total, n_genes),
  "--- planted marginals (round-trip target) ---",
  sprintf("associated : %.3f", mean(cpg$category != "background")),
  sprintf("both DE&DM : %.3f", mean(cpg$category == "both")),
  sprintf("DM (any)   : %.3f", mean(cpg$is_DM)),
  sprintf("host DE(any): %.3f", mean(cpg$host_isDE)),
  "--- quadrant split of the BOTH set ---",
  capture.output(print(round(prop.table(table(cpg$quadrant)), 3))),
  "--- region-type of DM CpGs (actual, after fallback) ---",
  capture.output(print(round(prop.table(table(cpg$region_type[cpg$is_DM])), 3)))
)
writeLines(sm, file.path(outdir, "design_summary.txt"))
cat(paste(sm, collapse = "\n"), "\n")