# 08_gsva.R
suppressPackageStartupMessages({
  library(GSVA)
  library(msigdbr)
  library(limma)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(matrixStats)
  library(readr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(EnhancedVolcano)
})

# -------------------- CONFIG --------------------
base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
expr_file <- file.path(base_dir, "All_ComBat_matrix.csv")
cov_file <- file.path(base_dir, "All_Covariates_augmented.csv")
stopifnot(file.exists(expr_file), file.exists(cov_file))
out_dir <- file.path(base_dir, "Results_GSVA")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- LOAD ----------------------
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE) %>% as.matrix()
mode(expr) <- "numeric"
cov <- read.csv(cov_file, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), cov$Sample))

# -------------------- HELPERS -------------------
mk_factor <- function(x) {
  f <- factor(ifelse(is.na(x) | x == "", "Unknown", x))
  if (nlevels(f) >= 2) f else NULL
}
mk_numeric <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  if (all(is.na(v))) {
    return(NULL)
  }
  v[is.na(v)] <- median(v, na.rm = TRUE)
  if (is.na(sd(v)) || sd(v) == 0) {
    return(NULL)
  }
  as.numeric(scale(v))
}

detect_ids <- function(ids) {
  if (all(grepl("^ENSG", ids))) {
    return("ENSEMBL")
  }
  if (all(grepl("^[0-9]+$", ids))) {
    return("ENTREZID")
  }
  "SYMBOL"
}
map_to_symbol <- function(ids, kind) {
  if (kind == "SYMBOL") {
    return(ids)
  }
  keys <- ids
  if (kind == "ENSEMBL") keys <- sub("\\.\\d+$", "", keys)
  sy <- suppressMessages(mapIds(org.Hs.eg.db, keys = keys, keytype = kind, column = "SYMBOL", multiVals = "first"))
  unname(sy)
}
collapse_by_symbol <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (any(duplicated(symbols))) {
    rs <- rowsum(mat, group = symbols)
    cnt <- table(symbols)
    cnt <- cnt[rownames(rs)]
    rs <- sweep(rs, 1, as.numeric(cnt), "/") # mean duplicates
    return(rs)
  } else {
    rownames(mat) <- symbols
    return(mat)
  }
}

mk_TvN_contrast <- function(design) {
  cn_tum <- grep("^Tissue2.*Tumor$", colnames(design), value = TRUE)
  cn_nor <- grep("^Tissue2.*Normal$", colnames(design), value = TRUE)
  if (length(cn_tum) != 1 || length(cn_nor) != 1) {
    return(NULL)
  }
  C <- matrix(0, nrow = ncol(design), ncol = 1, dimnames = list(colnames(design), "TvN"))
  C[cn_tum, 1] <- 1
  C[cn_nor, 1] <- -1
  C
}

# -------------------- PREP EXPRESSION FOR GSVA -------------------
expr <- expr[rowSds(expr) > 0, , drop = FALSE]
id_kind <- detect_ids(rownames(expr))
sym <- map_to_symbol(rownames(expr), ifelse(id_kind == "SYMBOL", "SYMBOL", id_kind))
expr_sym <- collapse_by_symbol(expr, sym)

# ---------- Robust gene-set fetcher ----------
get_sets <- function(coll, sub = NULL, keyword = NULL) {
  safe_msig <- function(...) suppressWarnings(try(msigdbr(species = "Homo sapiens", ...), silent = TRUE))
  candidates <- list(
    safe_msig(collection = coll),
    if (!is.null(sub)) safe_msig(collection = paste(coll, sub, sep = ":")) else NULL,
    if (!is.null(sub)) safe_msig(collection = paste(coll, gsub("^CP:", "", sub), sep = ":")) else NULL,
    safe_msig(category = coll, subcategory = sub),
    if (!is.null(sub)) safe_msig(category = coll, subcategory = gsub("^CP:", "", sub)) else NULL
  )
  df <- NULL
  for (d in candidates) {
    if (!is.null(d) && !inherits(d, "try-error") && is.data.frame(d) && nrow(d) > 0) {
      df <- d
      break
    }
  }
  if (is.null(df)) stop(sprintf("msigdbr could not return collection '%s' (sub '%s').", coll, sub))
  if (!is.null(sub) || !is.null(keyword)) {
    key <- if (is.null(keyword)) sub else keyword
    cols <- intersect(c("gs_subcat", "gs_subcategory", "subcategory", "collection", "gs_name"), names(df))
    keep <- Reduce(`|`, lapply(cols, function(cc) grepl(key, as.character(df[[cc]]), ignore.case = TRUE)))
    if (any(keep)) df <- df[keep, , drop = FALSE]
  }
  split(df$gene_symbol, df$gs_name)
}

sets_H <- get_sets("H")
sets_C2_KEGG <- get_sets("C2", sub = "CP:KEGG", keyword = "KEGG")
sets_C2_REAC <- get_sets("C2", sub = "CP:REACTOME", keyword = "REACTOME")
sets_C5_GOBP <- get_sets("C5", sub = "GO:BP", keyword = "BP|GOBP|GO:BP")

collections <- list(
  H            = sets_H,
  C2_KEGG      = sets_C2_KEGG,
  C2_REACTOME  = sets_C2_REAC,
  C5_GOBP      = sets_C5_GOBP
)

# -------------------- GSVA (with safe fallbacks) -----------------
run_one_gsva <- function(gsets, label) {
  gsets2 <- lapply(gsets, intersect, y = rownames(expr_sym))
  lens <- vapply(gsets2, length, integer(1))
  gsets2 <- gsets2[lens >= 5]
  message(sprintf("[%s] %d gene sets", label, length(gsets2)))
  if (length(gsets2) == 0) {
    return(NULL)
  }

  gs <- try(
    GSVA::gsva(
      expr = expr_sym,
      gset.idx.list = gsets2,
      method = "gsva",
      kcdf = "Gaussian",
      mx.diff = TRUE,
      parallel.sz = 1,
      verbose = FALSE
    ),
    silent = TRUE
  )

  if (inherits(gs, "try-error")) {
    param_fun <- try(getFromNamespace("GSVAParam", "GSVA"), silent = TRUE)
    if (!inherits(param_fun, "try-error")) {
      param <- param_fun(expr_sym, gsets2, kcdf = "Gaussian", mx.diff = TRUE)
      gs <- GSVA::gsva(param)
    } else {
      message(sprintf("[%s] GSVA unavailable; using z-mean pathway scores as fallback.", label))
      Xz <- t(scale(t(expr_sym)))
      Xz[is.na(Xz)] <- 0
      gs_mat <- vapply(gsets2, function(genes) {
        if (!length(genes)) {
          return(rep(NA_real_, ncol(Xz)))
        }
        colMeans(Xz[genes, , drop = FALSE])
      }, numeric(ncol(Xz)))
      gs <- t(gs_mat)
      rownames(gs) <- names(gsets2)
      colnames(gs) <- colnames(expr_sym)
    }
  }

  write.csv(gs, file.path(out_dir, paste0("GSVA_scores_", label, ".csv")), quote = FALSE)
  gs
}

gsva_scores <- lapply(names(collections), function(nm) run_one_gsva(collections[[nm]], nm))
names(gsva_scores) <- names(collections)
gsva_scores <- Filter(Negate(is.null), gsva_scores)
stopifnot(length(gsva_scores) > 0)

# -------------------- LIMMA ON GSVA SCORES (GLOBAL) --------------
cov$Tissue2 <- factor(ifelse(is.na(cov$Tissue), "Unknown", cov$Tissue),
  levels = c("Normal", "Tumor", "Unknown")
)

limma_on_scores <- function(scores, label) {
  dfrm <- data.frame(Tissue2 = cov$Tissue2)
  opt <- list(
    Study    = mk_factor(cov$Study),
    Platform = mk_factor(cov$Platform),
    Sex      = mk_factor(cov$Sex),
    Stage    = mk_factor(cov$Stage)
  )
  for (nm in names(opt)) if (!is.null(opt[[nm]])) dfrm[[nm]] <- opt[[nm]]
  AgeZ <- mk_numeric(cov$Age)
  if (!is.null(AgeZ)) dfrm$Age <- AgeZ

  design <- model.matrix(~ 0 + Tissue2 + ., data = dfrm)
  colnames(design) <- make.names(colnames(design))
  stopifnot(nrow(design) == ncol(scores))

  fit <- lmFit(scores, design)
  C <- mk_TvN_contrast(design)
  if (is.null(C)) {
    return(NULL)
  }
  fit2 <- eBayes(contrasts.fit(fit, C))
  tt <- topTable(fit2, coef = 1, n = Inf, sort.by = "P")
  tt$Collection <- label
  write.csv(tt, file.path(out_dir, paste0("GSVA_DGE_", label, "_Tumor_vs_Normal.csv")), row.names = TRUE)

  png(file.path(out_dir, paste0("Volcano_GSVA_", label, "_Tumor_vs_Normal.png")), width = 1400, height = 1000, res = 160)
  EnhancedVolcano(tt,
    lab = rownames(tt), x = "logFC", y = "P.Value",
    pCutoff = 0.05, FCcutoff = 0.15,
    title = paste("GSVA:", label, "Tumor vs Normal"),
    subtitle = "Adjusted covariates (NA→Unknown)"
  )
  dev.off()

  topn <- head(rownames(tt[order(abs(tt$t), decreasing = TRUE), , drop = FALSE]), 30)
  if (length(topn) >= 2) {
    ha <- ComplexHeatmap::HeatmapAnnotation(
      Study = cov$Study, Tissue = cov$Tissue2, Sex = cov$Sex,
      col = list(
        Study  = structure(scales::hue_pal()(length(unique(cov$Study))), names = unique(cov$Study)),
        Tissue = c(Normal = "#7f7f7f", Tumor = "#d62728", Unknown = "#1f77b4")
      )
    )
    png(file.path(out_dir, paste0("Heatmap_GSVA_", label, "_Top30.png")), width = 1400, height = 1000, res = 160)
    ComplexHeatmap::Heatmap(scores[topn, ],
      name = "GSVA", show_row_names = TRUE, show_column_names = FALSE,
      top_annotation = ha, cluster_columns = TRUE, cluster_rows = TRUE
    )
    dev.off()
  }
  tt
}

res_global <- lapply(names(gsva_scores), function(nm) limma_on_scores(gsva_scores[[nm]], nm))
names(res_global) <- names(gsva_scores)

# --------------- PER-DEMOGRAPHIC GSVA LIMMA (TvN) ----------------
per_demo_dir <- file.path(out_dir, "per_demo")
dir.create(per_demo_dir, showWarnings = FALSE)

do_per_group <- function(scores, label, varname) {
  lv <- unique(stats::na.omit(cov[[varname]]))
  for (v in lv) {
    idx <- !is.na(cov[[varname]]) & cov[[varname]] == v
    cov_g <- cov[idx, , drop = FALSE]
    scr_g <- scores[, idx, drop = FALSE]

    T2 <- factor(ifelse(is.na(cov_g$Tissue), "Unknown", cov_g$Tissue))
    if (!all(c("Normal", "Tumor") %in% levels(droplevels(T2)))) next
    if (ncol(scr_g) < 6) next

    dfrm <- data.frame(Tissue2 = factor(ifelse(is.na(cov_g$Tissue), "Unknown", cov_g$Tissue),
      levels = c("Normal", "Tumor", "Unknown")
    ))
    opt <- list(
      Study    = mk_factor(cov_g$Study),
      Platform = mk_factor(cov_g$Platform),
      Sex      = mk_factor(cov_g$Sex),
      Stage    = mk_factor(cov_g$Stage)
    )
    for (nm in names(opt)) if (!is.null(opt[[nm]])) dfrm[[nm]] <- opt[[nm]]
    AgeZ <- mk_numeric(cov_g$Age)
    if (!is.null(AgeZ)) dfrm$Age <- AgeZ

    design <- model.matrix(~ 0 + Tissue2 + ., data = dfrm)
    colnames(design) <- make.names(colnames(design))
    if (nrow(design) != ncol(scr_g)) next

    C <- mk_TvN_contrast(design)
    if (is.null(C)) next
    fit_g <- eBayes(contrasts.fit(lmFit(scr_g, design), C))
    ttg <- topTable(fit_g, coef = 1, n = Inf, sort.by = "P")
    fn <- paste0("GSVA_", label, "_", varname, "_", gsub("[^A-Za-z0-9]+", "_", v), "_Tumor_vs_Normal.csv")
    write.csv(ttg, file.path(per_demo_dir, fn), row.names = TRUE)
  }
}
for (nm in names(gsva_scores)) {
  do_per_group(gsva_scores[[nm]], nm, "Sex")
  do_per_group(gsva_scores[[nm]], nm, "Stage")
  do_per_group(gsva_scores[[nm]], nm, "Study")
  do_per_group(gsva_scores[[nm]], nm, "Platform")
}
message("GSVA completed. Outputs in: ", out_dir)

# ==================  NETWORK PREP (append)  ==================
message("==== Network prep: starting ====")
net_dir <- file.path(base_dir, "Results_network_prep")
dir.create(net_dir, showWarnings = FALSE, recursive = TRUE)

writeLines(
  capture.output(sessionInfo()),
  file.path(net_dir, "sessionInfo_06_gsva.txt")
)

stopifnot(ncol(expr) == nrow(cov))
if (anyDuplicated(cov$Sample)) {
  new_ids <- make.unique(cov$Sample)
  message("Made unique sample IDs in cov (", sum(duplicated(cov$Sample)), " duplicates).")
  cov$Sample <- new_ids
  colnames(expr) <- new_ids
  if (exists("expr_sym")) colnames(expr_sym) <- new_ids
}
stopifnot(identical(colnames(expr), cov$Sample))

cov_network <- cov
cov_network$Tissue <- ifelse(is.na(cov_network$Tissue), "Unknown", cov_network$Tissue)
cov_network$Sex <- ifelse(is.na(cov_network$Sex), "Unknown", cov_network$Sex)
cov_network$Stage <- ifelse(is.na(cov_network$Stage), "Unknown", cov_network$Stage)
cov_network$Study <- ifelse(is.na(cov_network$Study), "Unknown", cov_network$Study)
cov_network$Platform <- ifelse(is.na(cov_network$Platform), "Unknown", cov_network$Platform)
write.csv(cov_network, file.path(net_dir, "Covariates_network_clean.csv"), row.names = FALSE)

topN <- min(10000L, nrow(expr))
rv <- matrixStats::rowVars(expr)
idx <- order(rv, decreasing = TRUE)[seq_len(topN)]
expr_var <- expr[idx, , drop = FALSE]

expr_sym_var <- NULL
if (exists("expr_sym")) {
  topN_sym <- min(10000L, nrow(expr_sym))
  rv_sym <- matrixStats::rowVars(expr_sym)
  idx_sym <- order(rv_sym, decreasing = TRUE)[seq_len(topN_sym)]
  expr_sym_var <- expr_sym[idx_sym, , drop = FALSE]
}

tumor_idx <- cov_network$Tissue == "Tumor"
normal_idx <- cov_network$Tissue == "Normal"

expr_var_tumor <- expr_var[, tumor_idx, drop = FALSE]
expr_var_normal <- expr_var[, normal_idx, drop = FALSE]

if (!is.null(expr_sym_var)) {
  expr_sym_var_tumor <- expr_sym_var[, tumor_idx, drop = FALSE]
  expr_sym_var_normal <- expr_sym_var[, normal_idx, drop = FALSE]
}

safe_log1p <- function(m, eps = 1e-6) {
  m <- as.matrix(m)
  mn <- suppressWarnings(min(m, na.rm = TRUE))
  if (is.finite(mn) && mn <= -1) m <- m - mn + eps
  log1p(m)
}

expr_var_log1p <- safe_log1p(expr_var)
expr_var_tumor_log1p <- safe_log1p(expr_var_tumor)
expr_var_normal_log1p <- safe_log1p(expr_var_normal)

if (!is.null(expr_sym_var)) {
  expr_sym_var_log1p <- safe_log1p(expr_sym_var)
  expr_sym_var_tumor_log1p <- safe_log1p(expr_sym_var_tumor)
  expr_sym_var_normal_log1p <- safe_log1p(expr_sym_var_normal)
}

write.csv(expr_var, file.path(net_dir, "Expr_varTopN_all.csv"), quote = FALSE)
write.csv(expr_var_tumor, file.path(net_dir, "Expr_varTopN_Tumor.csv"), quote = FALSE)
write.csv(expr_var_normal, file.path(net_dir, "Expr_varTopN_Normal.csv"), quote = FALSE)
write.csv(expr_var_log1p, file.path(net_dir, "Expr_varTopN_all_log1p.csv"), quote = FALSE)
write.csv(expr_var_tumor_log1p, file.path(net_dir, "Expr_varTopN_Tumor_log1p.csv"), quote = FALSE)
write.csv(expr_var_normal_log1p, file.path(net_dir, "Expr_varTopN_Normal_log1p.csv"), quote = FALSE)

if (!is.null(expr_sym_var)) {
  write.csv(expr_sym_var, file.path(net_dir, "ExprSYMBOL_varTopN_all.csv"), quote = FALSE)
  write.csv(expr_sym_var_tumor, file.path(net_dir, "ExprSYMBOL_varTopN_Tumor.csv"), quote = FALSE)
  write.csv(expr_sym_var_normal, file.path(net_dir, "ExprSYMBOL_varTopN_Normal.csv"), quote = FALSE)
  write.csv(expr_sym_var_log1p, file.path(net_dir, "ExprSYMBOL_varTopN_all_log1p.csv"), quote = FALSE)
  write.csv(expr_sym_var_tumor_log1p, file.path(net_dir, "ExprSYMBOL_varTopN_Tumor_log1p.csv"), quote = FALSE)
  write.csv(expr_sym_var_normal_log1p, file.path(net_dir, "ExprSYMBOL_varTopN_Normal_log1p.csv"), quote = FALSE)
}

if (exists("gsva_scores") && length(gsva_scores) > 0) {
  saveRDS(gsva_scores, file.path(out_dir, "GSVA_scores_list.rds"))
}

if (exists("res_global") && length(res_global) > 0) {
  comb <- do.call(rbind, lapply(names(res_global), function(nm) {
    df <- res_global[[nm]]
    if (!is.null(df)) {
      df$Collection <- nm
      df$Pathway <- rownames(df)
      df
    }
  }))
  if (!is.null(comb)) {
    comb <- comb[!is.na(comb$P.Value), , drop = FALSE]
    rownames(comb) <- NULL
    write.csv(comb, file.path(out_dir, "GSVA_TvN_allCollections_results.csv"), row.names = FALSE)
  }
}

cat(paste0(
  "Network prep outputs:\n",
  "  - ", file.path(net_dir, "Expr_varTopN_all.csv"), "\n",
  "  - ", file.path(net_dir, "Expr_varTopN_Tumor.csv"), "\n",
  "  - ", file.path(net_dir, "Expr_varTopN_Normal.csv"), "\n",
  "  - ", file.path(net_dir, "Expr_varTopN_all_log1p.csv"), "\n",
  "  - ", file.path(net_dir, "Expr_varTopN_Tumor_log1p.csv"), "\n",
  "  - ", file.path(net_dir, "Expr_varTopN_Normal_log1p.csv"), "\n",
  if (!is.null(expr_sym_var)) {
    paste0(
      "  - ", file.path(net_dir, "ExprSYMBOL_varTopN_all.csv"), "\n",
      "  - ", file.path(net_dir, "ExprSYMBOL_varTopN_Tumor.csv"), "\n",
      "  - ", file.path(net_dir, "ExprSYMBOL_varTopN_Normal.csv"), "\n",
      "  - ", file.path(net_dir, "ExprSYMBOL_varTopN_all_log1p.csv"), "\n",
      "  - ", file.path(net_dir, "ExprSYMBOL_varTopN_Tumor_log1p.csv"), "\n",
      "  - ", file.path(net_dir, "ExprSYMBOL_varTopN_Normal_log1p.csv"), "\n"
    )
  } else {
    ""
  },
  "  - ", file.path(net_dir, "Covariates_network_clean.csv"), "\n",
  if (file.exists(file.path(out_dir, "GSVA_TvN_allCollections_results.csv"))) {
    paste0("  - ", file.path(out_dir, "GSVA_TvN_allCollections_results.csv"), "\n")
  } else {
    ""
  },
  if (file.exists(file.path(out_dir, "GSVA_scores_list.rds"))) {
    paste0("  - ", file.path(out_dir, "GSVA_scores_list.rds"), "\n")
  } else {
    ""
  }
))
message("==== Network prep: done ====")
