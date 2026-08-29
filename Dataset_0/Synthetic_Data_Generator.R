#!/usr/bin/env Rscript
# ============================================================================
# Round-trip synthetic data generator (reference method = M2/ChIPseeker)
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
# --- delta-level (paired case-ctrl) coupling, INDIPENDENTE dall'eQTM ---------
# Il canale Delta (2.a) e l'eQTM (2.b) leggono gli STESSI valori grezzi di
# metilazione/espressione simulati. Un CpG (o un gene) usato da ENTRAMBI i
# meccanismi crea interferenza: 2.a è statisticamente molto più esigente di 2.b
# (~2400 test con solo ~60 veri positivi, contro il segnale "pulito", non
# differenziato, della 2.b su 100 campioni), quindi il canale Delta richiede
# un'iniezione molto più forte -- e quel rumore extra, se applicato a un CpG/gene
# GIA' usato dall'eQTM, degrada anche il segnale eQTM (osservato empiricamente:
# con 40 CpG condivisi, l'eQTM è crollato da 0.99 a 0.57). Per eliminare questa
# interferenza, le coppie Delta sono scelte SENZA alcuna sovrapposizione di CpG
# NE' di gene host con le 100 coppie eQTM (v. sezione 5c) -- quindi 0 in comune.
rho_delta_eqtm <- 0.97     # forza del canale Delta: quasi tutta la varianza spiegata da Z (poco rumore idiosincratico)
# v4 (dopo l'isolamento 0-in-comune, eQTM tornato a 0.97 ma Delta osservato a 0.35):
# delta_meth_sd=1.0 era troppo forte per lo spazio PROBABILITA' -- a rho=0.97 il termine
# e' quasi deterministico (w_idio piccolo), quindi d_meth ~ N(0, delta_meth_sd) di fatto.
# Con baseline beta_ctrl=0.5 e clamp [0.01,0.99], una SD di 1.0 satura il clamp per una
# frazione enorme dei 50 soggetti per coppia (stimato numericamente: ~35% dei ctrl e
# ~52% dei case), appiattendo molti valori e distruggendo l'ordine dei ranghi di cui
# Spearman ha bisogno -- stesso tipo di guasto di v1, solo parziale invece che totale.
# Abbassato a 0.5: la saturazione evitabile sui ctrl scende a ~10% (quella residua sui
# case, ~35-39%, e' un pavimento preesistente dovuto a deltaM_mag=0.40 applicato a TUTTI
# i CpG "both", indipendente dal canale Delta). Il segnale deterministico via Z_delta resta
# comunque forte (coeff. ~0.49 invece di ~0.98). Va ancora validato empiricamente in R.
delta_meth_sd  <- 0.5      # SD (spazio PROBABILITA', non logit) del termine metilazione del canale Delta
n_delta_pairs  <- 60       # coppie CpG-gene per il canale Delta, tutte indipendenti dall'eQTM (0 in comune)
# --- DM-site multiplicity ("dose") effect on the host gene ------------------
# Genes hosting more DM CpGs (regardless of category, counted from cpg$is_DM)
# get: (1) higher odds of being DE, (2) a larger |log2FC| once DE. Direction
# (dir, up vs down) is still drawn at p_up=0.5 AFTER this, so the up/down
# balance from part 1 is untouched -- only DE-probability and effect magnitude
# scale with site count. This is what pipeline-2 sections 4 and 5 test for.
# CAUTION: dose_logor_beta raises P(DE) for EVERY gene with >=1 DM site, not just
# multi-site genes -- with p_DM=0.5 over 15000 CpGs, a large fraction of all genes
# have >=1 DM site, so a high value inflates the overall DE-gene count substantially.
# That larger DE universe feeds directly into the analysis script's per-method
# Spearman families (pipeline 1, sections 2.a/2.b): more tests in the same BH family
# generally means a LOWER effective power for every true signal in that family,
# including the pre-existing eQTM/Delta recall. If recall drops after raising this,
# lower it first before touching anything else -- 0.8 was too aggressive and measurably
# hurt eQTM/Delta recall; 0.3 below is a more conservative starting point.
dose_logor_beta <- 0.3     # log-odds of being DE added per DM site on the gene (0 = no effect, current behaviour)
dose_fc_beta    <- 0.5     # |log2FC| added per DM site beyond the first, for DE genes (0 = no effect)
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
is_distal <- runif(n) < frac_distal        # fixed once: preserves frac_distal exactly
chr_pool  <- unique(gene_tab$chr)
# (re)draw chr/pos for a set of row indices, respecting each row's distal flag.
# non-distal: random host gene, pos = TSS +/- log-uniform offset (strand-aware).
# distal:     random far positions -> M2 labels them "Distal Intergenic".
draw_positions <- function(idx) {
  ch <- character(length(idx)); ps <- integer(length(idx))
  nd <- !is_distal[idx]; ds <- is_distal[idx]
  if (any(nd)) {
    gs <- gene_tab[sample(nrow(gene_tab), sum(nd), replace = TRUE), ]
    d  <- round(10^runif(sum(nd), log10(tss_dist_range[1]), log10(tss_dist_range[2])))
    downstream <- runif(sum(nd)) < p_downstream
    sgn <- ifelse(gs$strand == "-", -1L, 1L) * ifelse(downstream, 1L, -1L)
    ch[nd] <- gs$chr
    ps[nd] <- pmax(1L, as.integer(gs$tss + sgn * d))
  }
  if (any(ds)) {
    ch[ds] <- sample(chr_pool, sum(ds), replace = TRUE)
    ps[ds] <- as.integer(runif(sum(ds), 1e4, 2e8))
  }
  list(chr = ch, pos = ps)
}
# initial draw, then resample ONLY the surplus copies of any colliding (chr,pos)
d0 <- draw_positions(seq_len(n)); chr <- d0$chr; pos <- d0$pos
guard <- 0L
repeat {
  key  <- paste(chr, pos)
  redo <- which(duplicated(key))           # 2nd+ occurrence; first of each key is kept
  if (!length(redo)) break
  dd <- draw_positions(redo); chr[redo] <- dd$chr; pos[redo] <- dd$pos
  guard <- guard + 1L
  if (guard > 100L) stop("dedup loop not converging: check tss_dist_range / n_cpg_total")
}
cpg$chr <- chr; cpg$pos <- as.integer(pos)
cpg$strand <- "+"    # single-strand projection (destranded CpGs)
stopifnot(!any(duplicated(paste(cpg$chr, cpg$pos))))   # closes the diagnosis: must hold
## ------------------------- 3. M2: ChIPseeker annotatePeak ------------------
txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
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
## --- DM-site multiplicity ("dose") per gene, computed from is_DM alone -----
## (available now: is_DM is per-CpG and independent of gene DE status, so this
## can be computed before is_DE/log2FC are drawn, with no circularity)
dm_count_tab <- table(cpg$m2_symbol[cpg$is_DM])
gene_de$n_DM_sites <- as.integer(dm_count_tab[gene_de$SYMBOL])
gene_de$n_DM_sites[is.na(gene_de$n_DM_sites)] <- 0L
## Odds of being DE increase with n_DM_sites (drives pipeline-2 section 4,
## "hyp_multi"): p_gene_DE is the baseline (n_DM_sites = 0), each extra DM
## site adds dose_logor_beta to the log-odds.
p_DE_i <- plogis(qlogis(p_gene_DE) + dose_logor_beta * gene_de$n_DM_sites)
gene_de$is_DE <- runif(nrow(gene_de)) < p_DE_i
gene_de$dir   <- 0L
de_rows <- which(gene_de$is_DE)
gene_de$dir[de_rows] <- ifelse(runif(length(de_rows)) < p_up, 1L, -1L)
## |log2FC| increases with n_DM_sites beyond the first (drives pipeline-2
## section 5, "sites vs |log2FC|" / "sites vs log2FC"); direction (dir) is
## still assigned above at p_up = 0.5, so the up/down balance is unaffected.
gene_de$log2FC <- gene_de$dir * (lfc_mag + dose_fc_beta * pmax(gene_de$n_DM_sites - 1, 0))
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
## ------------------------- 5b. eQTM injection setup (livelli, 100 coppie) --
stopifnot("category" %in% names(cpg), "m2_symbol" %in% names(cpg))   # dipendenze
rho_eqtm <- 0.8         # forza accoppiamento; 0 = nessuno (torna al comportamento attuale)
n_eqtm   <- 100
both_idx <- which(cpg$category == "both")
stopifnot(length(both_idx) >= n_eqtm)
eqtm_cpg <- sort(sample(both_idx, n_eqtm))          # 100 CpG target
cpg$is_eqtm <- FALSE
cpg$is_eqtm[eqtm_cpg] <- TRUE
eqtm_gene_sym <- unique(cpg$m2_symbol[eqtm_cpg])    # geni-host target (per l'espressione)
Z <- rnorm(n_pairs)                                 # latente per-soggetto, estratto UNA volta
k_expr <-  sqrt(rho_eqtm) * subj_sd_rna
k_meth <- -sqrt(rho_eqtm) * subj_sd_meth            # segno - : più metilato -> meno espresso
w_idio <- sqrt(1 - rho_eqtm)                        # peso della parte idiosincratica
## ------------------------- 5c. delta-level coupling, INDIPENDENTE dall'eQTM
## Il blocco 5b applica lo stesso termine soggetto-specifico a ctrl e case,
## quindi si cancella nel delta case-ctrl (pipeline 2.a). Qui aggiungiamo un
## SECONDO termine, con segno OPPOSTO tra ctrl (-1/2) e case (+1/2), che quindi
## sopravvive per intero nel delta. Media zero su ctrl/case -> non altera le
## medie di popolazione (log2FC, deltaM) né le proporzioni in design_summary.
## Le n_delta_pairs coppie sono scelte da "both" ESCLUDENDO sia i 100 CpG
## eQTM sia i geni host dei 100 CpG eQTM (eqtm_gene_sym) -- quindi 0 CpG e 0
## geni in comune con l'eQTM: nessuna interferenza possibile tra i due canali.
candidate_delta <- setdiff(both_idx, eqtm_cpg)
candidate_delta <- candidate_delta[!(cpg$m2_symbol[candidate_delta] %in% eqtm_gene_sym)]
stopifnot(length(candidate_delta) >= n_delta_pairs)
pick_distinct_genes <- function(idx_pool, k) {
  # campiona k indici da idx_pool assicurando geni host tutti diversi tra loro
  pool <- idx_pool[sample(length(idx_pool))]
  picked <- integer(0); used_genes <- character(0)
  for (i in pool) {
    g <- cpg$m2_symbol[i]
    if (!(g %in% used_genes)) { picked <- c(picked, i); used_genes <- c(used_genes, g) }
    if (length(picked) == k) break
  }
  if (length(picked) < k) stop("delta: candidati con gene host distinto insufficienti")
  sort(picked)
}
delta_eqtm_cpg   <- pick_distinct_genes(candidate_delta, n_delta_pairs)   # n_delta_pairs CpG totali
delta_host_genes <- unique(cpg$m2_symbol[delta_eqtm_cpg])
cpg$is_delta_eqtm <- FALSE
cpg$is_delta_eqtm[delta_eqtm_cpg] <- TRUE
gene_de$is_delta_gene <- gene_de$SYMBOL %in% delta_host_genes

Z_delta      <- rnorm(n_pairs)                          # latente del canale Delta, indipendente da Z (eQTM)
k_expr_delta <-  sqrt(rho_delta_eqtm) * subj_sd_rna
k_meth_delta <- -sqrt(rho_delta_eqtm) * delta_meth_sd   # delta_meth_sd, NON subj_sd_meth: scala diversa (v. sopra)
w_idio_delta <- sqrt(1 - rho_delta_eqtm)
is_delta_g <- gene_de$is_delta_gene   # geni host del canale Delta (per la sezione 6, RNA)
is_delta_c <- cpg$is_delta_eqtm       # CpG del canale Delta (per la sezione 7, WGBS)
## ------------------------- 6. simulate + write RNA (.sf, keyed by SYMBOL) --
write_sf <- function(counts, len, ids, path) {
  eff <- pmax(len - 200, 1); rate <- counts / eff
  tpm <- if (sum(rate) > 0) rate / sum(rate) * 1e6 else rep(0, length(rate))
  write.table(data.frame(Name = ids, Length = len, EffectiveLength = eff,
                         TPM = round(tpm, 4), NumReads = round(counts, 3)),
              path, sep = "\t", quote = FALSE, row.names = FALSE)
}
hit_rna <- gene_de$SYMBOL %in% eqtm_gene_sym                       ## <-- eQTM: geni target (una volta)
for (p in seq_len(n_pairs)) {
  re <- rnorm(nrow(gene_de), 0, subj_sd_rna)                       # shared subject baseline (paired)
  re[hit_rna] <- w_idio * rnorm(sum(hit_rna), 0, subj_sd_rna) +    ## <-- eQTM: parte condivisa Z[p]
    k_expr * Z[p]                                     ## <-- eQTM
  
  d_expr <- rep(0, nrow(gene_de))                                  ## <-- Delta: termine ASIMMETRICO
  if (any(is_delta_g)) d_expr[is_delta_g] <- w_idio_delta * rnorm(sum(is_delta_g), 0, subj_sd_rna) +
    k_expr_delta * Z_delta[p]                                      ## <-- 60 geni Delta: latente Z_delta, indipendente da Z (eQTM)
  re_ctrl <- re - 0.5 * d_expr                                     # split simmetrico: sopravvive nel delta
  re_case <- re + 0.5 * d_expr
  
  write_sf(rnbinom(nrow(gene_de), mu = gene_de$mu_ctrl * exp(re_ctrl), size = 1 / nb_dispersion),
           gene_de$length, gene_de$SYMBOL, file.path(outdir_rna, "rnaseq", sprintf("ctrl_%d.sf", p)))
  write_sf(rnbinom(nrow(gene_de), mu = gene_de$mu_case * exp(re_case), size = 1 / nb_dispersion),
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
  re[cpg$is_eqtm] <- w_idio * rnorm(sum(cpg$is_eqtm), 0, subj_sd_meth) + ## <-- eQTM: parte condivisa Z[p]
    k_meth * Z[p]                                       ## <-- eQTM (segno - per anti-corr)
  
  ## v3 (ricalibrata numericamente, v. commento sopra rho_delta_eqtm/delta_meth_sd):
  ## d_meth torna in SPAZIO PROBABILITA' (come v1) -- lo spazio logit (v2) evitava la
  ## saturazione ma smorzava troppo il segnale tramite la derivata della sigmoide, che
  ## a n=50 e su una famiglia di ~2400 test con solo 60 veri positivi non sopravviveva a
  ## BH (recall osservato 0.03). La differenza rispetto a v1 e' che ora d_meth usa
  ## delta_meth_sd (1.0, dedicato, NON subj_sd_meth) con rho_delta_eqtm=0.97: piu' forte
  ## e piu' deterministico, per compensare la doppia rumorosita' di una differenza a 50
  ## coppie contro il segnale "pulito" a 100 campioni dell'eQTM 2.b.
  d_meth <- rep(0, n)                                                   ## <-- Delta: termine ASIMMETRICO
  if (any(is_delta_c)) d_meth[is_delta_c] <- w_idio_delta * rnorm(sum(is_delta_c), 0, delta_meth_sd) +
    k_meth_delta * Z_delta[p]                                           ## <-- 60 CpG Delta: latente Z_delta, indipendente da Z (eQTM)
  
  bc <- clampb(invlogit(logit(clampb(cpg$beta_ctrl)) + re) - 0.5 * d_meth)   # ctrl: -1/2 del termine
  write_wgbs(bc, "ctrl", p)
  write_wgbs(clampb(bc + cpg$deltaM + d_meth), "case", p)                   # case: bc + deltaM + termine
  # delta (case-ctrl, pre-clamp) = deltaM + d_meth -> d_meth sopravvive intero nel delta,
  # mentre la media (ctrl+case)/2 resta invariata rispetto allo script originale.
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
  sprintf("M2 genes: %d", nrow(gene_de)),
  sprintf("DE genes: %d", sum(gene_de$is_DE)),
  sprintf("CpGs: %d",     n),
  sprintf("DM CpGs: %d",  sum(cpg$is_DM)),
  print(""),
  print("--- REALISED marginals (round-trip target; emergent under M2) ---"),
  sprintf("DM (any)   : %.3f", mean(cpg$is_DM)),
  sprintf("both DE&DM : %.3f", mean(cpg$category == "both")),
  sprintf("dm_only    : %.3f", mean(cpg$category == "dm_only")),
  "--- quadrant split of the BOTH set (realised) ---",
  capture.output(print(round(prop.table(table(cpg$quadrant)), 3))),
  print(""),
  "--- M2 region_type of DM CpGs (coarse) ---",
  capture.output(print(round(prop.table(table(coarse(cpg$m2_region[cpg$is_DM]))), 3))),
  print(""),
  "--- eQTM-like couplings injected (full list: eqtm_pairs.tsv / delta_eqtm_pairs.tsv) ---",
  sprintf("condition-pooled eQTM (pipeline 2.b, 'Spearman eQTM'): %d CpG-gene pairs expected, rho = %.2f",
          length(eqtm_cpg), rho_eqtm),
  sprintf("paired-delta eQTM   (pipeline 2.a, 'Spearman Delta') : %d CpG-gene pairs expected, rho = %.2f",
          length(delta_eqtm_cpg), rho_delta_eqtm),
  sprintf("  fully INDEPENDENT of the %d condition-pooled eQTM above: 0 shared CpGs and 0 shared",
          length(eqtm_cpg)),
  print(""),
  "--- DM-site multiplicity ('dose') effect on the host gene ---",
  sprintf("log-odds of being DE, per extra DM site on the gene : +%.2f (0 = no effect; baseline p_gene_DE = %.2f)",
          dose_logor_beta, p_gene_DE),
  sprintf("|log2FC| added per DM site beyond the first, if DE  : +%.2f (0 = no effect; baseline lfc_mag = %.2f)",
          dose_fc_beta, lfc_mag),
  "  => genes with more DM sites are more likely to be DE AND show a larger |log2FC|",
  "     (magnitude only). This does NOT set the sign of log2FC",
  print(""),
  "--- directional (up vs down) balance: NO up/down preference by design ---",
  sprintf("DE genes up   : %d (%.3f of DE genes)",
          sum(gene_de$dir ==  1L), mean(gene_de$dir[gene_de$is_DE] ==  1L)),
  sprintf("DE genes down : %d (%.3f of DE genes)",
          sum(gene_de$dir == -1L), mean(gene_de$dir[gene_de$is_DE] == -1L))
)
writeLines(sm, file.path(outdir, "Ground_truth_for_comparison.txt"))
cat(paste(sm, collapse = "\n"), "\n")
## <-- eQTM: elenco esplicito delle 100 coppie target CpG-gene
eqtm_pairs <- data.frame(
  coord_key = paste(cpg$chr[eqtm_cpg], cpg$pos[eqtm_cpg], sep = "_"),
  cpg_id    = cpg$cpg_id[eqtm_cpg],
  SYMBOL    = cpg$m2_symbol[eqtm_cpg]
)
write.table(eqtm_pairs, file.path(outdir, "eqtm_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
## <-- eQTM-delta: elenco esplicito delle 60 coppie target CpG-gene (0 in comune con l'eQTM)
delta_eqtm_pairs <- data.frame(
  coord_key = paste(cpg$chr[delta_eqtm_cpg], cpg$pos[delta_eqtm_cpg], sep = "_"),
  cpg_id    = cpg$cpg_id[delta_eqtm_cpg],
  SYMBOL    = cpg$m2_symbol[delta_eqtm_cpg]
)
write.table(delta_eqtm_pairs, file.path(outdir, "delta_eqtm_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)