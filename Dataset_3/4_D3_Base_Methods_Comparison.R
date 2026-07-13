library(dplyr)
library(DESeq2)
library(UpSetR)
library(ggplot2)
library(cowplot)
library(writexl)


### Razionale Analisi: 

##Pipeline 0: Considera quanto ogni metodo identifica il sottoinsieme di geni intersezione DE DM

##Pipeline 1: Considera quante delle associazioni nell'intersezione sono "significative" per ogni metodo. 
# Se un metodo ha una % di significativi più alta si può dire che "fa associazioni più affidabili" 

##Pipeline 2 mi restituisce quanto ogni metodo è in grado di catturare il trend generale dei siti DM. 
# Dà risultati diversi dalla pipeline 1 ma vanno interpretati correttamente:
# Non c'è intersezione DE, quindi questa mi dice: quale metodo trova più "quadranti", e "precision in generale
# !! da capire se fare intersezione con DE


#QUINDI:  PIP0 è esplorativa sulla lista di geni
#         PIP1 chiede se metilazione ed espressione co-variano a livello individuale. 
#         PIP2 chiede se la direzione dell'effetto disease è concorde tra le due omiche.


PATH <- "C:/Users/pierp/Desktop/Thesis PROJECT"

pdf(file.path(PATH, "Dataset_3", "4_Integration_results", "Base_comparison D3.pdf"))

#Importing:
#Importing RNA-Seq DE genes
sign_DE <- read.csv(file.path(PATH, "Dataset_3", "1_RNA-Seq", "DE_results.csv"))
sign_DE <- sign_DE %>% dplyr::rename(SYMBOL = hugo_symbol)  #renaming for coherence

#This will be different for the different lists:
M1_df <- read.csv(file.path(PATH, "Dataset_3", "3_Benchmark", "DM_sites_Met1.csv"))
M2_df <- read.csv(file.path(PATH, "Dataset_3", "3_Benchmark", "DM_sites_Met2.csv"))
M3_df <- read.csv(file.path(PATH, "Dataset_3", "3_Benchmark", "DM_sites_Met3.csv"))
M4_df <- read.csv(file.path(PATH, "Dataset_3", "3_Benchmark", "DM_sites_Met4.csv"))
M5_df <- read.csv(file.path(PATH, "Dataset_3", "3_Benchmark", "DM_sites_Met5.csv"))


##### --- PIPELINE 0: Simple Intersection --- #####

#Imserting BS-Seq DM sites into a list
M_list <- list()
M_list[[1]] <- M1_df
M_list[[2]] <- M2_df
M_list[[3]] <- M3_df
M_list[[4]] <- M4_df
M_list[[5]] <- M5_df


inters_table <- sapply(seq_along(M_list), function(i) {
  M_inters <- intersect(sign_DE$SYMBOL, M_list[[i]]$SYMBOL)
  c("tot genes" = as.numeric(length(unique(M_list[[i]]$SYMBOL))),
    "inters genes" = as.numeric(length(M_inters)),
    "precision" = as.numeric(length(M_inters)/length(unique(M_list[[i]]$SYMBOL))),
    "unique genes" =  length(setdiff(M_list[[i]]$SYMBOL, 
                                     unlist(lapply(M_list[-i], function(x) x$SYMBOL)))),
    "Recall genes" = length(unique(M_inters)) / nrow(sign_DE)
  )
})
inters_table <- as.data.frame(inters_table)
colnames(inters_table) <- paste0("M", seq_along(M_list))
names(M_list) <- paste0("M", seq_along(M_list))


##Upset plot pre-Statistic
M_symbol_list <- lapply(M_list, function(x) {
  s <- x$SYMBOL
  s <- sub("\\.\\d+$", "", s)    # ".1", ".2", ... from repeated genes
  unique(s)
})

#Plot intersection
upset(fromList(M_symbol_list),
      mainbar.y.label = "Intersecting genes - Pre-statistic",
      sets.x.label = "Tot genes per method")
## Last column from this Upset plot shows CORE genes



## !! Important table. Table of the upset plot with gene names
### Obtaining and exporting the names of the genes for categories in the upset plot:

#Using M_symbol_list, used for the Upset plot
all_genes  <- unique(unlist(M_symbol_list))
membership <- sapply(M_symbol_list, function(g) all_genes %in% g)
rownames(membership) <- all_genes

#Category labels (es. "M1", "M2&M3", "M1&M2&M3&M4&M5")
category_label <- apply(membership, 1, function(present) {
  paste(names(M_symbol_list)[present], collapse = "&")
})

#Grouping genes for category
genes_by_region <- split(all_genes, category_label)

# Sorting categories: first for # of Method (degree), then alphabetical
degree <- sapply(strsplit(names(genes_by_region), "&"), length)
genes_by_region <- genes_by_region[order(degree, names(genes_by_region))]

# Naming columns with category:
#   "M1"       -> "unique M1"
#   "M1&M2"    -> "M1&M2"
cate_names <- sapply(names(genes_by_region), function(r) {
  parts <- strsplit(r, "&")[[1]]
  if (length(parts) == 1)                       return(paste("unique", parts))
  if (length(parts) == length(M_symbol_list))   return("CORE (all methods)")
  paste(parts, collapse = "&")
})

#Elonging all columns to the longest one, adding empty values
max_len <- max(sapply(genes_by_region, length))
Upset_genes <- lapply(genes_by_region, function(g) {
  g <- sort(g)
  c(g, rep("", max_len - length(g)))   # empty cells
})

Upset_genes_df <- as.data.frame(Upset_genes, stringsAsFactors = FALSE, check.names = FALSE)
colnames(Upset_genes_df) <- cate_names

write_xlsx(Upset_genes_df, file.path(PATH, "Dataset_3", "4_Integration_results", "Upset_genes.xlsx"))



#####
##### --- PIPELINE 1: Vector Comparison --- #####

#This pipeline only considers the sites that were found to be differentially methylated in the condition (TB-CTRL)
#These (for each different benchmark method) are intersected with DE genes (from DESeq2).
#On this gene set a statistic is obtained, to assess how many of these are "significant" associations. Then methods will be compared on this.

### FUNCTIONS ###

create_DE_matrix <- function(symbols, expr_matrix, DE_symbols) {
  #Selecting only genes that are also DE: intersection
  symbols <- symbols[symbols %in% DE_symbols]                    #Comment this to see the case without intersection (Without this line the intersection DM-DE is not done). Substitute with: symbols <- symbols[symbols %in% rownames(expr_matrix)]
  #Filtering
  Matrix_expr <- expr_matrix[symbols, , drop = FALSE]
  #Changing col names (rimuovo prefisso "S")
  colnames(Matrix_expr) <- sub("^S", "", colnames(Matrix_expr))
  # Making duplicated rownames unique
  rownames(Matrix_expr) <- make.unique(symbols)
  
  return(Matrix_expr)
}

create_DM_matrix <- function(Method_df, M_matrix, DE_matrix) {
  #Selecting couples (coord_key, SYMBOL) in M_matrix
  keep <- Method_df$coord_key %in% rownames(M_matrix)
  Method_df <- Method_df[keep, , drop = FALSE]
  
  #Filtering: only keeping genes that were in the universe (i.e. in common for the WGBS and RNA-Seq experiment)
  keep_com <- Method_df$SYMBOL %in% rownames(DE_matrix)
  Method_df <- Method_df[keep_com, , drop = FALSE]
  
  #If the same coord_key is present more than once, all cases are kept
  Matrix_M <- M_matrix[Method_df$coord_key, , drop = FALSE]
  #Changing row.names to SYMBOL
  rownames(Matrix_M) <- make.unique(Method_df$SYMBOL)  #To reconnect to the CpG, order must be followed: the first of the duplicated genes is associated to the first CpG in the "DM_sites_M" df.
  
  return(Matrix_M)
}




### 1. Importing files and Matrix creation ###

## Importing files

#Importing RNA-Seq data
dds <- readRDS(file.path(PATH, "Dataset_3", "1_RNA-Seq", "dds.rds"))

##Methylation: Creating the DM matrix:

# Obtaining the matrix for M-values (log2((numCs + 1) / (numTs + 1)))
Matrix_Mv <- read.csv(file.path(PATH, "Dataset_3", "2_BS-Seq", "M_CpG_matrix.csv"), check.names = FALSE)
rownames(Matrix_Mv) <- Matrix_Mv[, 1]
Matrix_Mv <- Matrix_Mv[, -1]
Matrix_Mv <- as.matrix(Matrix_Mv)

##RNA expression: Creating the DE matrix
#Normalization: vst()
vst_dds <- vst(dds, blind = FALSE)
vst_expr_matrix <- assay(vst_dds)

# Reorder vst_expr_matrix from interleaved
d0_cols  <- grep("_D0$",  colnames(vst_expr_matrix), value = TRUE)
d28_cols <- grep("_D28$", colnames(vst_expr_matrix), value = TRUE)
d0_cols  <- d0_cols[order(as.integer(sub("_D0$",  "", d0_cols)))]
d28_cols <- d28_cols[order(as.integer(sub("_D28$", "", d28_cols)))]
vst_expr_matrix <- vst_expr_matrix[, c(d0_cols, d28_cols)]

### 2.a Spearman Correlation Delta ###

#Subtracting the ctrl mean to ctrl samples and case mean to case samples

#methylation
# Delta paired: subject k -> olig2_k - neun_k
Matrix_Mv_new <- Matrix_Mv[, 43:84] - Matrix_Mv[, 1:42]

#expression
vst_expr_matrix_new <- vst_expr_matrix[, 43:84] - vst_expr_matrix[, 1:42]


Matrix_exp <- list()
Matrix_exp[[1]] <- create_DE_matrix(M1_df$SYMBOL, vst_expr_matrix_new, sign_DE$SYMBOL)
Matrix_exp[[2]] <- create_DE_matrix(M2_df$SYMBOL, vst_expr_matrix_new, sign_DE$SYMBOL)
Matrix_exp[[3]] <- create_DE_matrix(M3_df$SYMBOL, vst_expr_matrix_new, sign_DE$SYMBOL)
Matrix_exp[[4]] <- create_DE_matrix(M4_df$SYMBOL, vst_expr_matrix_new, sign_DE$SYMBOL)
Matrix_exp[[5]] <- create_DE_matrix(M5_df$SYMBOL, vst_expr_matrix_new, sign_DE$SYMBOL)
## Not all initial genes are mantained here. Only the ones that were present in the "universe"


## Obtaining the Matrices Mval for each different method:
Matrix_Mval <- list()
Matrix_Mval[[1]] <- create_DM_matrix(M1_df, Matrix_Mv_new, Matrix_exp[[1]])
Matrix_Mval[[2]] <- create_DM_matrix(M2_df, Matrix_Mv_new, Matrix_exp[[2]])
Matrix_Mval[[3]] <- create_DM_matrix(M3_df, Matrix_Mv_new, Matrix_exp[[3]])
Matrix_Mval[[4]] <- create_DM_matrix(M4_df, Matrix_Mv_new, Matrix_exp[[4]])
Matrix_Mval[[5]] <- create_DM_matrix(M5_df, Matrix_Mv_new, Matrix_exp[[5]])



## Final check. Before statistics, we must check that the matrices are coherent
for (m in seq_along(Matrix_exp)) {
  stopifnot(identical(rownames(Matrix_exp[[m]]), rownames(Matrix_Mval[[m]])))
}
message("CHECK PASSED")


Spear_res_list <- list()
Spear_padj_list <- list()
for(m in seq_along(Matrix_exp)) {
  #If there is no DE DM intersection, skip
  if (nrow(Matrix_exp[[m]]) == 0) {
    res <- data.frame(rho = numeric(0), pvalue = numeric(0),
                      padj = numeric(0), SYMBOL = character(0))
    Spear_res_list[[m]]  <- res
    Spear_padj_list[[m]] <- res
    next
  }
  
  res_t <- sapply(seq_len(nrow(Matrix_exp[[m]])), function (i) {
    corr <- cor.test(Matrix_Mval[[m]][i, ], Matrix_exp[[m]][i, ], method = 'spearman')
    c(rho = unname(corr$estimate), pvalue = corr$p.value)
  })
  
  res <- t(res_t)
  res <- as.data.frame(res)
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res$SYMBOL <- row.names(Matrix_exp[[m]])
  
  Spear_res_list[[m]] <- res
  Spear_padj_list[[m]] <- res %>% filter(padj < 0.05)
}
names(Spear_res_list) <- paste0("M", seq_along(Matrix_exp))
names(Spear_padj_list) <- paste0("M", seq_along(Matrix_exp))


### 2.b Spearman Correlation eQTM ###

#Subtracting the ctrl mean to ctrl samples and case mean to case samples
#methylation
D0_mean <- rowMeans(Matrix_Mv[, 1:42, drop = FALSE])
D28_mean <- rowMeans(Matrix_Mv[, 43:84, drop = FALSE])

Matrix_Mv_Q <- Matrix_Mv
Matrix_Mv_Q[, 1:42] <- Matrix_Mv_Q[, 1:42] - D0_mean
Matrix_Mv_Q[, 43:84] <- Matrix_Mv_Q[, 43:84] - D28_mean

#expression
D0_vst_mean <- rowMeans(vst_expr_matrix[, 1:42, drop = FALSE])
D28_vst_mean <- rowMeans(vst_expr_matrix[, 43:84, drop = FALSE])

vst_expr_matrix_Q <- vst_expr_matrix
vst_expr_matrix_Q[, 1:42] <- vst_expr_matrix_Q[, 1:42] - D0_vst_mean
vst_expr_matrix_Q[, 43:84] <- vst_expr_matrix_Q[, 43:84] - D28_vst_mean

Matrix_exp_Q <- list()
Matrix_exp_Q[[1]] <- create_DE_matrix(M1_df$SYMBOL, vst_expr_matrix_Q, sign_DE$SYMBOL)
Matrix_exp_Q[[2]] <- create_DE_matrix(M2_df$SYMBOL, vst_expr_matrix_Q, sign_DE$SYMBOL)
Matrix_exp_Q[[3]] <- create_DE_matrix(M3_df$SYMBOL, vst_expr_matrix_Q, sign_DE$SYMBOL)
Matrix_exp_Q[[4]] <- create_DE_matrix(M4_df$SYMBOL, vst_expr_matrix_Q, sign_DE$SYMBOL)
Matrix_exp_Q[[5]] <- create_DE_matrix(M5_df$SYMBOL, vst_expr_matrix_Q, sign_DE$SYMBOL)
## Not all initial genes are mantained here. Only the ones that were present in the "universe"


## Obtaining the Matrices Mval for each different method:
Matrix_Mval_Q <- list()
Matrix_Mval_Q[[1]] <- create_DM_matrix(M1_df, Matrix_Mv_Q, Matrix_exp_Q[[1]])
Matrix_Mval_Q[[2]] <- create_DM_matrix(M2_df, Matrix_Mv_Q, Matrix_exp_Q[[2]])
Matrix_Mval_Q[[3]] <- create_DM_matrix(M3_df, Matrix_Mv_Q, Matrix_exp_Q[[3]])
Matrix_Mval_Q[[4]] <- create_DM_matrix(M4_df, Matrix_Mv_Q, Matrix_exp_Q[[4]])
Matrix_Mval_Q[[5]] <- create_DM_matrix(M5_df, Matrix_Mv_Q, Matrix_exp_Q[[5]])

##Spearman
Spear_res_Q_list <- list()
Spear_padj_Q_list <- list()
for(m in seq_along(Matrix_exp_Q)) {
  #If there is no DE DM intersection, skip
  if (nrow(Matrix_exp_Q[[m]]) == 0) {
    res <- data.frame(rho = numeric(0), pvalue = numeric(0),
                      padj = numeric(0), SYMBOL = character(0))
    Spear_res_Q_list[[m]]  <- res
    Spear_padj_Q_list[[m]] <- res
    next
  }
  res_t <- sapply(seq_len(nrow(Matrix_exp_Q[[m]])), function (i) {
    corr <- cor.test(Matrix_Mval_Q[[m]][i, ], Matrix_exp_Q[[m]][i, ], method = 'spearman')
    c(rho = unname(corr$estimate), pvalue = corr$p.value)
  })
  
  res <- t(res_t)
  res <- as.data.frame(res)
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res$SYMBOL <- row.names(Matrix_exp_Q[[m]])
  
  Spear_res_Q_list[[m]] <- res
  Spear_padj_Q_list[[m]] <- res %>% filter(padj < 0.05)
}
names(Spear_res_Q_list) <- paste0("M", seq_along(Matrix_exp_Q))
names(Spear_padj_Q_list) <- paste0("M", seq_along(Matrix_exp_Q))


### 3. Visualizations for METHOD COMPARISONS ###

## Identifying unique/intersecting/common genes to all methods ##

## 3.1 Upset plot Spearman Delta ##
Spear_symbol_list <- lapply(Spear_padj_list, function(x) {
  s <- x$SYMBOL
  s <- sub("\\.\\d+$", "", s)    # ".1", ".2", ... from repeated genes
  unique(s)
})
names(Spear_symbol_list) <- names(Spear_res_list)

#Plot significant Spearman Delta
if (sum(lengths(Spear_symbol_list) > 0) >= 2) {
  print(upset(fromList(Spear_symbol_list),
              mainbar.y.label = "Intersecting sign genes - Spearman Delta",
              sets.x.label = "Tot genes per method"))
} else {
  message("SKIP upset Delta: less than 2 non-empty sets")
}

#Upset plot Spearman eQTM
Spear_symbol_Q_list <- lapply(Spear_padj_Q_list, function(x) {
  s <- x$SYMBOL
  s <- sub("\\.\\d+$", "", s)    # ".1", ".2", ... from repeated genes
  unique(s)
})
names(Spear_symbol_Q_list) <- names(Spear_res_Q_list)

#Plot significant Spearman eQTM
if (sum(lengths(Spear_symbol_Q_list) > 0) >= 2) {
  print(upset(fromList(Spear_symbol_Q_list),
              mainbar.y.label = "Intersecting sign genes - Spearman eQTM",
              sets.x.label = "Tot genes per method"))
} else {
  message("SKIP upset eQTM: less than 2 non-empty sets")
}



## 3.2 Spearman correlation visualization ##

# Rho distribution in methods

#Delta
df_Spear <- bind_rows(Spear_padj_list, .id = "method")
#Density
ggplot(df_Spear, aes(x = rho, fill = method)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Spearman rho Delta", y = "Density") +
  theme_minimal()

#eQTM
df_Spear_Q <- bind_rows(Spear_padj_Q_list, .id = "method")
#Density
ggplot(df_Spear_Q, aes(x = rho, fill = method)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Spearman rho eQTM", y = "Density") +
  theme_minimal()


## 3.3 Barplot: significant vs tested pairs (stacked) ##

# Builds a long df: for each method -> Significant + Non significant counts
make_sig_barplot_df <- function(res_list, padj_list) {
  do.call(rbind, lapply(seq_along(res_list), function(i) {
    tot <- nrow(res_list[[i]])
    sig <- nrow(padj_list[[i]])
    data.frame(
      method   = names(res_list)[i],
      category = c("Significant", "Non significant"),
      count    = c(sig, tot - sig),
      stringsAsFactors = FALSE
    )
  }))
}

plot_sig_bar <- function(df, res_list, title) {
  # keep method order M1..M5 (not alphabetical surprises)
  df$method   <- factor(df$method, levels = names(res_list))
  # stacking order: Significant on top of Non significant
  df$category <- factor(df$category, levels = c("Non significant", "Significant"))
  
  # totals + significant fraction, for labels
  totals <- do.call(rbind, lapply(names(res_list), function(m) {
    sub <- df[df$method == m, ]
    tot <- sum(sub$count)
    sig <- sub$count[sub$category == "Significant"]
    data.frame(method = m, tot = tot, sig = sig,
               perc = ifelse(tot > 0, 100 * sig / tot, 0))
  }))
  totals$method <- factor(totals$method, levels = names(res_list))
  
  ggplot(df, aes(x = method, y = count, fill = category)) +
    geom_col(width = 0.7) +
    # total tested pairs on top of each bar
    geom_text(data = totals, aes(x = method, y = tot, label = tot),
              inherit.aes = FALSE, vjust = -0.4, size = 3) +
    # % significant inside the colored segment
    geom_text(data = totals,
              aes(x = method, y = tot, label = sprintf("%d (%.1f%%)", sig, perc)),
              inherit.aes = FALSE, vjust = 1.4, size = 2.8, color = "grey20") +
    scale_fill_manual(values = c("Non significant" = "grey80",
                                 "Significant"     = "#D55E00")) +
    labs(title = title, x = "Method",
         y = "Tested CpG-gene pairs", fill = NULL) +
    theme_minimal()
}

bar_delta <- make_sig_barplot_df(Spear_res_list,   Spear_padj_list)
bar_qtm   <- make_sig_barplot_df(Spear_res_Q_list, Spear_padj_Q_list)

print(plot_sig_bar(bar_delta, Spear_res_list,
                   "Significant vs tested pairs - Spearman Delta"))
print(plot_sig_bar(bar_qtm,   Spear_res_Q_list,
                   "Significant vs tested pairs - Spearman eQTM"))


## 3.4 Percentage of significant and unique/different genes ##
strip_suffix <- function(x) unique(sub("\\.\\d+$", "", x))

compare_table <- data.frame(
  "Method" = names(Spear_res_list),
  "Tot associations - intersection" = sapply(seq_along(Matrix_Mval), function(j) {
    nrow(Matrix_Mval[[j]])
  }),
  "Perc sig CpG-gene Delta" = sapply(seq_along(Spear_padj_list), function(i) {
    length(Spear_padj_list[[i]]$padj)/length(Spear_res_list[[i]]$padj)
  }),
  "Perc sig genes Delta" = sapply(seq_along(Spear_padj_list), function(i) {
    sig_genes <- length(unique(sub("\\.\\d+$", "", Spear_padj_list[[i]]$SYMBOL)))
    tot_genes <- length(unique(sub("\\.\\d+$", "", Spear_res_list[[i]]$SYMBOL)))
    sig_genes / tot_genes
  }),
  "Perc unique sig Delta" = sapply(seq_along(Spear_padj_list), function(i) {
    unique_sig <- setdiff(strip_suffix(Spear_padj_list[[i]]$SYMBOL),
                          strip_suffix(unlist(lapply(Spear_padj_list[-i], function(x) x$SYMBOL))))
    unique_all <- setdiff(strip_suffix(Spear_res_list[[i]]$SYMBOL),
                          strip_suffix(unlist(lapply(Spear_res_list[-i], function(x) x$SYMBOL))))
    length(unique_sig) / length(unique_all)
  }),
  "Perc sig CpG-gene eQTM" = sapply(seq_along(Spear_padj_Q_list), function(i) {
    length(Spear_padj_Q_list[[i]]$padj)/length(Spear_res_Q_list[[i]]$padj)
  }),
  "Perc sig genes eQTM" = sapply(seq_along(Spear_padj_Q_list), function(i) {
    sig_genes <- length(unique(sub("\\.\\d+$", "", Spear_padj_Q_list[[i]]$SYMBOL)))
    tot_genes <- length(unique(sub("\\.\\d+$", "", Spear_res_Q_list[[i]]$SYMBOL)))
    sig_genes / tot_genes
  }),
  "Perc unique sig eQTM" = sapply(seq_along(Spear_padj_Q_list), function(i) {
    unique_sig <- setdiff(strip_suffix(Spear_padj_Q_list[[i]]$SYMBOL),
                          strip_suffix(unlist(lapply(Spear_padj_Q_list[-i], function(x) x$SYMBOL))))
    unique_all <- setdiff(strip_suffix(Spear_res_Q_list[[i]]$SYMBOL),
                          strip_suffix(unlist(lapply(Spear_res_Q_list[-i], function(x) x$SYMBOL))))
    length(unique_sig) / length(unique_all)
  }),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

Pip1_table <- as.data.frame(t(compare_table))
colnames(Pip1_table) <- compare_table$Method
Pip1_table <- Pip1_table[-1, ]



#####
##### --- PIPELINE 2: Scalar comparison --- #####

#In this pipeline every DM-DE association is compared with 1 scalar value: log2FC (expr) and mean Mvalue[TB-Ctrl] (methyl) 
#What is extracted is the directionality of the whole dataset's associations: quadrants


### FUNCTIONS ###

create_integr_df <- function(Method_dataframe, Mean_mv_dataframe, rnaseqFC_dataframe) {
  df <- data.frame(
    SYMBOL = Method_dataframe$SYMBOL,
    coord_key = Method_dataframe$coord_key,
    Mv = Mean_mv_dataframe$Mv_mean[match(Method_dataframe$coord_key, Mean_mv_dataframe$coord_key)],
    log2FC = rnaseqFC_dataframe$log2FoldChange[match(Method_dataframe$SYMBOL, rnaseqFC_dataframe$SYMBOL)]
  )
  return(df)
}



### 1. Importing files and matrix/dataframes creation ###

##Methylation: Matrix_Mv already created

#Obtaining the delta matrix
#D0_mean and D28 were already called in PIP 1

Mean_Mv_df <- data.frame(coord_key = rownames(Matrix_Mv),
                         Mv_mean = D28_mean - D0_mean)


##Obtaining the integration matrices: For each CpG --> its mean Mvalue(subtracted conditions) + associated gene/s + its log2FC in the RNA-Seq
Method_final_df <- list()
Method_final_df[[1]] <- create_integr_df(M1_df, Mean_Mv_df, sign_DE)
Method_final_df[[2]] <- create_integr_df(M2_df, Mean_Mv_df, sign_DE)
Method_final_df[[3]] <- create_integr_df(M3_df, Mean_Mv_df, sign_DE)
Method_final_df[[4]] <- create_integr_df(M4_df, Mean_Mv_df, sign_DE)
Method_final_df[[5]] <- create_integr_df(M5_df, Mean_Mv_df, sign_DE)
## Not all intial genes are mantained here. Only the ones that were present in the "universe"



#Filtering rows with NA (SYMBOLS not in common between Method e rnaseq are NA after "create_integr_df()")
Method_final_df <- lapply(Method_final_df, function(df) {
  n_in  <- nrow(df)
  df    <- df[complete.cases(df), ]
  message(sprintf("Kept %d / %d pair (%.1f%%)",
                  nrow(df), n_in, 100 * nrow(df) / n_in))
  df
})
names(Method_final_df) <- paste0("M", seq_along(Method_final_df))



### 2. Comparison

# How many expected association does each method find? 
quadrant_enrichment <- sapply(Method_final_df, function(df) {
  gene_df <- df %>%
    group_by(coord_key) %>%
    summarise(Mv_med = median(Mv),
              log2FC = unique(log2FC)[1],
              .groups = "drop") %>%
    filter(Mv_med != 0, log2FC != 0)
  
  #If there is no DE DM intersection, return NA
  if (nrow(gene_df) == 0) {
    return(c(expected_perc    = NA_real_,
             odds_ratio       = NA_real_,
             pvalue_quadrants = NA_real_))
  }
  
  q1 <- sum(gene_df$Mv_med > 0 & gene_df$log2FC > 0)
  q2 <- sum(gene_df$Mv_med < 0 & gene_df$log2FC > 0)  # hypo + up
  q3 <- sum(gene_df$Mv_med < 0 & gene_df$log2FC < 0)
  q4 <- sum(gene_df$Mv_med > 0 & gene_df$log2FC < 0)  # hyper + down
  
  expected_total <- q2 + q4
  unexpected_total <- q1 + q3
  
  #And statistics: Binomial test: H0 = 50/50 split, H1 = expected > unexpected
  p_bin_test <- binom.test(expected_total, expected_total + unexpected_total, p = 0.5, alternative = "greater")$p.value
  
  #This next one (expected_perc) is the fundamental metric: the % of points in the expected quadrants.
  c(expected_perc = 100 * expected_total / (expected_total + unexpected_total),
    odds_ratio = expected_total / unexpected_total, 
    pvalue_quadrants = p_bin_test)
})

quadrant_table <- as.data.frame(quadrant_enrichment)

#Multiple test correction: BH
quadrant_table["pvalue_quadrants_BH", ] <- p.adjust(
  as.numeric(quadrant_table["pvalue_quadrants", ]), 
  method = "BH"
)

## Graph: quadrants
all_Mv <- unlist(lapply(Method_final_df, '[[', "Mv"))
all_log2FC <- unlist(lapply(Method_final_df, '[[', "log2FC"))

Mmax <- max(abs(all_Mv))
FCmax <- max(abs(all_log2FC))

xlims <- c(-Mmax, Mmax)
ylims <- c(-FCmax, FCmax)


#Adding intersection with DE genes:
gg_list <- list()
for(m in seq_along(Method_final_df)) {
  gg_list[[m]] <- ggplot(Method_final_df[[m]], aes(x = Mv, y = log2FC)) +
    geom_point(alpha = 0.5) +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    geom_point(data = Method_final_df[[m]][Method_final_df[[m]]$SYMBOL %in% sign_DE$SYMBOL, ], color = "red") +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    theme_minimal() +
    labs(title = sprintf("Method %d: quadrants", m), x = "M value", y = "log2FC")
}

grid <- plot_grid(plotlist = gg_list, nrow = 2, ncol = 3)
print(grid)

#Exporting the Method_df list
write_xlsx(Method_final_df, file.path(PATH, "Dataset_3", "4_Integration_results", "Method_final_df.xlsx"))

## Writing outputs
Stats_table <- rbind(inters_table, Pip1_table, quadrant_table)
Stats_table$metric <- rownames(Stats_table)
Stats_table <- Stats_table[ , c("metric", setdiff(names(Stats_table), "metric"))]  # Puts "metric" as first column
write.csv(Stats_table, file.path(PATH, "Dataset_3", "4_Integration_results", "Stats_table.csv"), row.names = FALSE)

dev.off()