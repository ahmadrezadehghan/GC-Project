# 05_all_combat.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(sva)
  library(matrixStats)
  library(umap)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(readr)
})

base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
expr_file <- file.path(base_dir, "Full_matrix.csv")
cov_file <- file.path(base_dir, "CombinedCovariates.csv")
demo_file <- file.path(base_dir, "Demo.xlsx")
stopifnot(file.exists(expr_file), file.exists(cov_file), file.exists(demo_file))

out_dir <- file.path(base_dir, "Results_all")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- load
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cov <- read.csv(cov_file, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), cov$Sample))

# ---------- read Demo.xlsx to harvest optional covariates (Sex/Stage/Age/Site)
norm_key <- function(x) tolower(gsub("\\s+", " ", trimws(as.character(x))))
demo <- readxl::read_excel(demo_file) %>% as.data.frame()
features_raw <- demo[[1]]
features_key <- norm_key(features_raw)
demo_vals <- function(idx) {
  r <- demo[idx, -1, drop = TRUE]
  names(r) <- names(demo)[-1]
  r
}
find_row <- function(cands) {
  hit <- which(features_key %in% norm_key(cands))
  ifelse(length(hit) >= 1, hit[1], NA_integer_)
}

pull_vec <- function(keys) {
  idx <- find_row(keys)
  if (!is.na(idx)) demo_vals(idx) else setNames(rep(NA_character_, ncol(demo) - 1), names(demo)[-1])
}

sex_vec <- pull_vec(c("Sex", "Gender", "Biological Sex"))
stage_vec <- pull_vec(c("Stage", "Tumor Stage", "Pathologic Stage", "Clinical Stage"))
age_vec <- pull_vec(c("Age", "Age at Diagnosis", "Age_years", "Age (years)"))
site_vec <- pull_vec(c("Primary Site", "Cancer Type", "Disease", "Site", "Organ", "Tumor Type"))

# align to samples (add if missing in CombinedCovariates)
mapv <- function(v) {
  x <- v[cov$Sample]
  if (is.null(x)) rep(NA, length(cov$Sample)) else as.character(x)
}
cov$Sex <- if (!"Sex" %in% names(cov) || all(is.na(cov$Sex))) trimws(ifelse(is.na(mapv(sex_vec)), NA, mapv(sex_vec))) else cov$Sex
cov$Stage <- if (!"Stage" %in% names(cov) || all(is.na(cov$Stage))) trimws(ifelse(is.na(mapv(stage_vec)), NA, mapv(stage_vec))) else cov$Stage
if (!"Age" %in% names(cov)) {
  cov$Age <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", mapv(age_vec))))
}
if (!"Site" %in% names(cov)) {
  cov$Site <- trimws(ifelse(is.na(mapv(site_vec)), NA, mapv(site_vec)))
}

# tidy Stage
cov$Stage <- gsub("stage\\s*", "", tolower(cov$Stage))
cov$Stage <- toupper(gsub("(^[ivx]+).*", "\\1", cov$Stage)) # I, II, III, IV from roman
cov$Stage[!cov$Stage %in% c("I", "II", "III", "IV")] <- NA

# ---------- ensure Study column (derived from Platform if absent)
if (!"Study" %in% names(cov)) {
  cov$Study <- ifelse(cov$Platform == "RNAseq", "RNAseq", "Microarray")
}

# ---------- build a global batch across all samples
# For microarray-like samples, prefer Series or Org×Country if available; else Study
org_vec <- pull_vec(c("Organization", "Organisation", "Institute", "Laboratory", "Submitter Organization", "Center", "Centre"))
country_vec <- pull_vec(c("Country", "Origin Country", "Submitter Country", "Location", "Nation"))
series_vec <- pull_vec(c("Series", "GSE", "GEO Series", "Series ID", "Study Accession", "Project ID", "Study"))

org <- trimws(ifelse(is.na(org_vec[cov$Sample]), "NA", as.character(org_vec[cov$Sample])))
cty <- trimws(ifelse(is.na(country_vec[cov$Sample]), "NA", as.character(country_vec[cov$Sample])))
ser <- trimws(ifelse(is.na(series_vec[cov$Sample]), "NA", as.character(series_vec[cov$Sample])))

batch_all <- cov$Study
is_micro <- cov$Platform == "Microarray"
use_series <- is_micro & ser != "NA"
batch_all[use_series] <- paste0("SERIES:", ser[use_series])
use_orgcty <- is_micro & !use_series & !(org == "NA" & cty == "NA")
batch_all[use_orgcty] <- paste0("ORGCTY:", org[use_orgcty], "_", cty[use_orgcty])
batch_all <- factor(batch_all)

# ---------- split multi-sample vs singletons
tb <- table(batch_all)
multi_lvls <- names(tb)[tb >= 2]
singleton_lvls <- names(tb)[tb == 1]
idx_multi <- batch_all %in% multi_lvls
idx_single <- batch_all %in% singleton_lvls
expr_multi <- as.matrix(expr[, idx_multi, drop = FALSE])
batch_multi <- droplevels(batch_all[idx_multi])

# ---------- CLEAN Sex to canonical values (Male/Female/NA)
sx <- tolower(cov$Sex)
sx[sx %in% c("m", "male", "man", "boy")] <- "male"
sx[sx %in% c("f", "female", "woman", "girl")] <- "female"
sx[!sx %in% c("male", "female")] <- NA
cov$Sex <- ifelse(is.na(sx), NA, ifelse(sx == "male", "Male", "Female"))

# re-align multi-sample covariates after cleaning
cov_multi <- cov[idx_multi, , drop = FALSE]

# ---------- helpers for safe ComBat
mk_factor <- function(x) {
  f <- factor(ifelse(is.na(x) | x == "", "Unknown", x))
  if (nlevels(f) < 2) NULL else f
}
mk_numeric <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  if (all(is.na(v))) {
    return(NULL)
  }
  v[is.na(v)] <- median(v, na.rm = TRUE)
  if (sd(v) == 0 || is.na(sd(v))) {
    return(NULL)
  }
  v
}
build_mod_df <- function(cov_sub) {
  cols <- list(Intercept = rep(1, nrow(cov_sub)))
  sx <- mk_factor(cov_sub$Sex)
  if (!is.null(sx)) cols$Sex <- sx
  st <- mk_factor(cov_sub$Stage)
  if (!is.null(st)) cols$Stage <- st
  ag <- mk_numeric(cov_sub$Age)
  if (!is.null(ag)) cols$Age <- as.numeric(scale(ag))
  as.data.frame(cols, check.names = FALSE)
}
run_combat_subset <- function(dat, batch, cov_sub) {
  if (ncol(dat) == 0 || nlevels(batch) <= 1) {
    return(as.matrix(dat))
  }
  df_cov <- build_mod_df(cov_sub) # NOTE: Tissue is NOT in mod (we stratify by Tissue)
  drop_order <- c("Stage", "Sex", "Age") # progressively drop if confounded
  for (k in 0:length(drop_order)) {
    df_try <- df_cov
    if (k >= 1) {
      drop_cols <- drop_order[1:k]
      df_try <- df_try[, setdiff(colnames(df_try), drop_cols), drop = FALSE]
      if (ncol(df_try) == 0) df_try <- data.frame(Intercept = rep(1, ncol(dat)))
    }
    mod <- if (ncol(df_try) > 1) model.matrix(~ . - 1, data = df_try) else model.matrix(~1, data = df_try)
    if (nrow(mod) != ncol(dat)) stop("Design/sample mismatch in subset.")
    res <- try(ComBat(
      dat = as.matrix(dat), batch = batch, mod = mod,
      par.prior = TRUE, prior.plots = FALSE, mean.only = FALSE
    ), silent = TRUE)
    if (!inherits(res, "try-error")) {
      return(res)
    }
  }
  # final fallback: intercept-only
  mod0 <- model.matrix(~1, data = data.frame(Intercept = rep(1, ncol(dat))))
  ComBat(
    dat = as.matrix(dat), batch = batch, mod = mod0,
    par.prior = TRUE, prior.plots = FALSE, mean.only = FALSE
  )
}

# ---------- STRATIFY BY TISSUE and run ComBat per stratum (on multi-sample batches)
tissue_ms <- factor(ifelse(is.na(cov_multi$Tissue), "Unknown", cov_multi$Tissue),
  levels = c("Normal", "Tumor", "Unknown")
)
expr_multi_cb_list <- list()
for (lev in levels(tissue_ms)) {
  idx_lev <- which(tissue_ms == lev)
  dat <- expr_multi[, idx_lev, drop = FALSE]
  bat <- droplevels(batch_multi[idx_lev])
  cov_sub <- cov_multi[idx_lev, , drop = FALSE]
  expr_multi_cb_list[[lev]] <- run_combat_subset(dat, bat, cov_sub)
}
expr_multi_cb <- do.call(cbind, expr_multi_cb_list)
expr_multi_cb <- expr_multi_cb[, colnames(expr_multi), drop = FALSE] # restore original multi-sample order

# ---------- combine back (singletons unchanged)
expr_all_cb <- cbind(expr_multi_cb, as.matrix(expr[, idx_single, drop = FALSE]))
expr_all_cb <- expr_all_cb[, colnames(expr), drop = FALSE]

# ---------- UMAP overview (by Study / Tissue)
set.seed(42)
um_dir <- file.path(out_dir, "UMAP")
dir.create(um_dir, showWarnings = FALSE)
um_cfg <- umap.defaults
um_cfg$n_neighbors <- min(15, ncol(expr_all_cb) - 1)
um_cfg$min_dist <- 0.2

um_bef <- umap::umap(t(expr), config = um_cfg)$layout
sel <- rowVars(expr_all_cb)
topn <- min(3000, max(2, sum(sel > 0)))
um_aft <- umap::umap(t(expr_all_cb[order(sel, decreasing = TRUE)[1:topn], , drop = FALSE]), config = um_cfg)$layout

df_bef_s <- data.frame(UMAP1 = um_bef[, 1], UMAP2 = um_bef[, 2], Color = cov$Study, Title = "Before (Study)")
df_aft_s <- data.frame(UMAP1 = um_aft[, 1], UMAP2 = um_aft[, 2], Color = cov$Study, Title = "After (Study)")
df_bef_t <- data.frame(UMAP1 = um_bef[, 1], UMAP2 = um_bef[, 2], Color = ifelse(is.na(cov$Tissue), "Unknown", cov$Tissue), Title = "Before (Tissue)")
df_aft_t <- data.frame(UMAP1 = um_aft[, 1], UMAP2 = um_aft[, 2], Color = ifelse(is.na(cov$Tissue), "Unknown", cov$Tissue), Title = "After (Tissue)")

plot_panel <- function(df) {
  ggplot(df, aes(UMAP1, UMAP2, color = Color)) +
    geom_point(size = 1.8, alpha = 0.85) +
    scale_color_viridis_d() +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom", legend.box = "horizontal",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = unique(df$Title), x = NULL, y = NULL, color = NULL)
}
p <- (plot_panel(df_bef_s) | plot_panel(df_aft_s)) / (plot_panel(df_bef_t) | plot_panel(df_aft_t)) +
  patchwork::plot_annotation(title = "ALL samples — UMAP Before/After global ComBat")

ggsave(file.path(um_dir, "ALL_UMAP_4panel.png"), p, width = 12, height = 10, dpi = 300)

# ---------- save
write.csv(expr_all_cb, file.path(base_dir, "All_ComBat_matrix.csv"), quote = FALSE)
write.csv(cov, file.path(base_dir, "All_Covariates_augmented.csv"), row.names = FALSE)
message("Wrote: All_ComBat_matrix.csv, All_Covariates_augmented.csv, and UMAP panel.")
