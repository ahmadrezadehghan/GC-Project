# 07_edge.R
suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
  library(dplyr)
  library(readr)
  library(matrixStats)
  library(circlize)
})

base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
expr_file <- file.path(base_dir, "All_ComBat_matrix.csv")
cov_file <- file.path(base_dir, "All_Covariates_augmented.csv")
stopifnot(file.exists(expr_file), file.exists(cov_file))

out_dir <- file.path(base_dir, "Results_DGE")
dir.create(out_dir, showWarnings = FALSE)

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE) %>% as.matrix()
mode(expr) <- "numeric"
cov <- read.csv(cov_file, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), cov$Sample))

# ---------- helpers
nz <- function(x) x[!is.na(x)]
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

# ---------- GLOBAL: Tumor vs Normal (robust design — only informative covariates)
cov$Tissue2 <- factor(ifelse(is.na(cov$Tissue), "Unknown", cov$Tissue),
  levels = c("Normal", "Tumor", "Unknown")
)
if (!all(c("Normal", "Tumor") %in% levels(droplevels(cov$Tissue2)))) {
  stop("Need both Tumor and Normal for global DGE.")
}

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

# sanity
stopifnot(!anyNA(design))
stopifnot(nrow(design) == ncol(expr))

fit <- lmFit(expr, design)

# robust contrast (بدون makeContrasts رشته‌ای)
cn_tum <- grep("^Tissue2.*Tumor$", colnames(design), value = TRUE)
cn_nor <- grep("^Tissue2.*Normal$", colnames(design), value = TRUE)
stopifnot(length(cn_tum) == 1, length(cn_nor) == 1)

C <- matrix(0, nrow = ncol(design), ncol = 1, dimnames = list(colnames(design), "Tumor_vs_Normal"))
C[cn_tum, 1] <- 1
C[cn_nor, 1] <- -1

fit2 <- eBayes(contrasts.fit(fit, C))
tt <- topTable(fit2, coef = 1, n = Inf, sort.by = "P")
write.csv(tt, file.path(out_dir, "DGE_global_Tumor_vs_Normal.csv"), row.names = TRUE)

# Volcano
png(file.path(out_dir, "Volcano_global_Tumor_vs_Normal.png"), width = 1400, height = 1000, res = 160)
EnhancedVolcano(tt,
  lab = rownames(tt), x = "logFC", y = "P.Value", pCutoff = 0.05, FCcutoff = log2(1.5),
  title = "Global DGE: Tumor vs Normal", subtitle = "Adjusted covariates (NA→Unknown)"
)
dev.off()

# Heatmap (top 50 by |t|)
top50 <- rownames(head(tt[order(abs(tt$t), decreasing = TRUE), , drop = FALSE], 50))
if (length(top50) >= 2) {
  ha <- ComplexHeatmap::HeatmapAnnotation(
    Study = cov$Study, Tissue = cov$Tissue2, Sex = cov$Sex,
    col = list(
      Study  = structure(scales::hue_pal()(length(unique(cov$Study))), names = unique(cov$Study)),
      Tissue = c(Normal = "#7f7f7f", Tumor = "#d62728", Unknown = "#1f77b4")
    )
  )
  png(file.path(out_dir, "Heatmap_global_Top50.png"), width = 1400, height = 1000, res = 160)
  ComplexHeatmap::Heatmap(expr[top50, ],
    name = "Expr", show_row_names = TRUE, show_column_names = FALSE,
    top_annotation = ha, cluster_columns = TRUE, cluster_rows = TRUE
  )
  dev.off()
}

# ---------- PER-DEMOGRAPHIC DGE (Tumor vs Normal within each group)
per_out <- file.path(out_dir, "per_demo")
dir.create(per_out, showWarnings = FALSE)

build_design <- function(cov_sub) {
  dfrm <- data.frame(
    Tissue2 = factor(ifelse(is.na(cov_sub$Tissue), "Unknown", cov_sub$Tissue),
      levels = c("Normal", "Tumor", "Unknown")
    )
  )
  # include covariates only if informative in this subgroup
  opt <- list(
    Study    = mk_factor(cov_sub$Study),
    Platform = mk_factor(cov_sub$Platform),
    Sex      = mk_factor(cov_sub$Sex),
    Stage    = mk_factor(cov_sub$Stage)
  )
  for (nm in names(opt)) if (!is.null(opt[[nm]])) dfrm[[nm]] <- opt[[nm]]
  AgeZ <- mk_numeric(cov_sub$Age)
  if (!is.null(AgeZ)) dfrm$Age <- AgeZ

  design <- model.matrix(~ 0 + Tissue2 + ., data = dfrm)
  colnames(design) <- make.names(colnames(design))

  # need both Tumor and Normal columns to exist
  cn_tum <- grep("^Tissue2.*Tumor$", colnames(design), value = TRUE)
  cn_nor <- grep("^Tissue2.*Normal$", colnames(design), value = TRUE)
  if (length(cn_tum) != 1 || length(cn_nor) != 1) {
    return(NULL)
  }

  list(
    design = design,
    contrast = {
      C <- matrix(0, nrow = ncol(design), ncol = 1, dimnames = list(colnames(design), "Tumor_vs_Normal"))
      C[cn_tum, 1] <- 1
      C[cn_nor, 1] <- -1
      C
    }
  )
}

do_per_group <- function(varname) {
  lv <- unique(stats::na.omit(cov[[varname]]))
  for (v in lv) {
    idx <- !is.na(cov[[varname]]) & cov[[varname]] == v
    cov_g <- cov[idx, , drop = FALSE]
    expr_g <- expr[, idx, drop = FALSE]

    # need both classes present and enough samples
    T2 <- factor(ifelse(is.na(cov_g$Tissue), "Unknown", cov_g$Tissue))
    if (!all(c("Normal", "Tumor") %in% levels(droplevels(T2)))) next
    if (ncol(expr_g) < 6) next

    bd <- build_design(cov_g)
    if (is.null(bd)) next
    fit_g <- eBayes(contrasts.fit(lmFit(expr_g, bd$design), bd$contrast))
    ttg <- topTable(fit_g, coef = 1, n = Inf, sort.by = "P")
    fn <- paste0("DGE_", varname, "_", gsub("[^A-Za-z0-9]+", "_", v), "_Tumor_vs_Normal.csv")
    write.csv(ttg, file.path(per_out, fn), row.names = TRUE)
  }
}
for (vn in c("Sex", "Stage", "Study", "Platform")) do_per_group(vn)

message("DGE completed.")
