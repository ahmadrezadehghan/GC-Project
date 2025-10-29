# 04_normalize_and_merge.R
# Batch correct Microarray (GEO/GSM-like) using ComBat (multi-sample batches only) and plot UMAP/PCA if possible

suppressPackageStartupMessages({
  library(dplyr)
  library(sva)
  library(matrixStats)
  library(umap)
  library(ggplot2)
  library(viridis)
  library(readxl)
  library(patchwork)
  library(readr)
})

# -------------------- CONFIG --------------------
base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
expr_file <- file.path(base_dir, "Full_matrix.csv")
cov_file <- file.path(base_dir, "CombinedCovariates.csv")
demo_file <- file.path(base_dir, "Demo.xlsx")
stopifnot(file.exists(expr_file), file.exists(cov_file), file.exists(demo_file))

# -------------------- LOAD ----------------------
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cov <- read.csv(cov_file, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), cov$Sample))

# ----------------- SELECT MICROARRAY -----------------
# (به‌جای cov$Study از Platform استفاده می‌کنیم)
keep <- cov$Platform == "Microarray"
expr_geo <- expr[, keep, drop = FALSE]
cov_geo <- cov[keep, , drop = FALSE]

ns <- ncol(expr_geo)
if (ns == 0) stop("No Microarray samples found (Platform == 'Microarray').")
message(sprintf("Microarray samples found: %d", ns))

# ----------------- DEMO HELPERS ------------------
norm_key <- function(x) tolower(gsub("\\s+", " ", trimws(as.character(x))))

# read demo and prepare a robust row finder (case-insensitive, trims spaces)
demo <- readxl::read_excel(demo_file) %>% as.data.frame()
features_raw <- demo[[1]]
features_key <- norm_key(features_raw)
demo_vals <- function(idx) {
  r <- demo[idx, -1, drop = TRUE]
  names(r) <- names(demo)[-1]
  r
}

find_row <- function(candidates) {
  cand <- norm_key(candidates)
  hit <- which(features_key %in% cand)
  if (length(hit) >= 1) hit[1] else NA_integer_
}

# try multiple options for organization / country / series
idx_org <- find_row(c(
  "Organization", "Organisation", "Organization Name",
  "Institute", "Institution", "Center", "Centre",
  "Laboratory", "Lab", "Submitter Organization"
))
idx_ctry <- find_row(c(
  "Country", "Country ", "Origin Country", "Submitter Country",
  "Location", "Nation"
))
idx_series <- find_row(c(
  "Series", "GSE", "Series ID", "GEO Series", "Study Accession",
  "Study", "Project", "Project ID"
))

# extract vectors (may be NA)
mk_na_vec <- function() setNames(rep(NA_character_, ncol(demo) - 1), names(demo)[-1])
org_vec <- if (!is.na(idx_org)) demo_vals(idx_org) else mk_na_vec()
country_vec <- if (!is.na(idx_ctry)) demo_vals(idx_ctry) else mk_na_vec()
series_vec <- if (!is.na(idx_series)) demo_vals(idx_series) else mk_na_vec()

# align to Microarray samples
organization <- trimws(ifelse(is.na(org_vec[cov_geo$Sample]), "NA", as.character(org_vec[cov_geo$Sample])))
country <- trimws(ifelse(is.na(country_vec[cov_geo$Sample]), "NA", as.character(country_vec[cov_geo$Sample])))
series <- trimws(ifelse(is.na(series_vec[cov_geo$Sample]), "NA", as.character(series_vec[cov_geo$Sample])))

# ----------------- BUILD BATCH -------------------
# Priority: Series/GSE if available; else Organization×Country; else Platform; else single batch
use_series <- any(series != "NA") && length(unique(series[series != "NA"])) >= 2
if (use_series) {
  batch_geo <- factor(series)
} else {
  has_org_or_ctry <- !(all(organization == "NA") & all(country == "NA"))
  if (has_org_or_ctry) {
    batch_geo <- interaction(organization, country, drop = TRUE)
  } else if (!is.null(cov_geo$Platform) && length(unique(cov_geo$Platform)) > 1) {
    batch_geo <- factor(cov_geo$Platform)
  } else {
    batch_geo <- factor(rep("Microarray", nrow(cov_geo))) # single batch
  }
}

# ----------------- BASIC FILTERS -----------------
# zero-variance genes
expr_geo <- expr_geo[rowVars(as.matrix(expr_geo)) > 0, , drop = FALSE]

# Tissue with Unknown kept (این اسکریپت Unknown را حذف نمی‌کند)
tissue3 <- factor(ifelse(is.na(cov_geo$Tissue), "Unknown", cov_geo$Tissue),
  levels = c("Normal", "Tumor", "Unknown")
)

# ----------------- COMBAT while KEEPING SINGLETON BATCHES -----------------
tb <- table(batch_geo)
multi_lvls <- names(tb)[tb >= 2]
singleton_lvls <- names(tb)[tb == 1]

idx_multi <- batch_geo %in% multi_lvls
idx_singleton <- batch_geo %in% singleton_lvls

expr_multi <- expr_geo[, idx_multi, drop = FALSE]
batch_multi <- droplevels(batch_geo[idx_multi])
cov_multi <- cov_geo[idx_multi, , drop = FALSE]
tissue3_multi <- droplevels(tissue3[idx_multi])

run_combat <- (ncol(expr_multi) > 0) && (nlevels(batch_multi) > 1)

if (run_combat) {
  # Use tissue covariate only if both classes present among multi-sample batches
  present <- levels(droplevels(tissue3_multi[tissue3_multi != "Unknown"]))
  if (length(present) >= 2) {
    mod <- model.matrix(~tissue3_multi)
  } else {
    mod <- model.matrix(~1, data = data.frame(dummy = rep(1, ncol(expr_multi))))
  }
  message("Running ComBat on multi-sample batches (mean & variance)…")
  expr_multi_cb <- ComBat(
    dat = as.matrix(expr_multi),
    batch = batch_multi, mod = mod,
    par.prior = TRUE, prior.plots = FALSE,
    mean.only = FALSE
  )
  message("ComBat done.")
} else {
  expr_multi_cb <- as.matrix(expr_multi) # nothing to correct
  if (ncol(expr_multi) > 0) message("Only one batch (or none) among multi-sample groups; skipping ComBat for them.")
}

# Combine back: corrected multi-sample + untouched singletons, in original order
expr_geo_cb <- cbind(expr_multi_cb, as.matrix(expr_geo[, idx_singleton, drop = FALSE]))
expr_geo_cb <- expr_geo_cb[, colnames(expr_geo), drop = FALSE] # restore column order

# ----------------- UMAP (guarded) ----------------
set.seed(42)
um_ready <- ncol(expr_geo_cb) >= 3
if (um_ready) {
  um_cfg <- umap.defaults
  um_cfg$n_neighbors <- min(15, ncol(expr_geo_cb) - 1)
  um_cfg$min_dist <- 0.2
  um_cfg$random_state <- 42

  um_bef <- umap::umap(t(expr_geo), config = um_cfg)$layout

  sel <- rowVars(expr_geo_cb)
  topn <- min(3000, max(2, sum(sel > 0)))
  expr_sel <- expr_geo_cb[order(sel, decreasing = TRUE)[1:topn], , drop = FALSE]
  um_aft <- umap::umap(t(expr_sel), config = um_cfg)$layout

  df_bef_b <- data.frame(UMAP1 = um_bef[, 1], UMAP2 = um_bef[, 2], Color = batch_geo, Title = "Before (Batch)")
  df_aft_b <- data.frame(UMAP1 = um_aft[, 1], UMAP2 = um_aft[, 2], Color = batch_geo, Title = "After (Batch)")
  df_bef_t <- data.frame(UMAP1 = um_bef[, 1], UMAP2 = um_bef[, 2], Color = tissue3, Title = "Before (Tissue)")
  df_aft_t <- data.frame(UMAP1 = um_aft[, 1], UMAP2 = um_aft[, 2], Color = tissue3, Title = "After (Tissue)")

  batch_levels <- levels(batch_geo)
  if (is.null(batch_levels)) batch_levels <- "Batch"
  batch_cols <- viridis(max(1, length(batch_levels)))
  names(batch_cols) <- batch_levels
  tissue_cols <- c(Normal = "#7f7f7f", Tumor = "#d62728", Unknown = "#1f77b4")

  plot_panel <- function(df, pal) {
    ggplot(df, aes(UMAP1, UMAP2, color = Color)) +
      geom_point(size = 2.2, alpha = 0.85) +
      {
        if (is.factor(df$Color) || is.character(df$Color)) {
          scale_color_manual(values = pal, drop = FALSE)
        } else {
          scale_color_viridis_c()
        }
      } +
      theme_bw(base_size = 13) +
      theme(
        legend.position = "bottom", legend.box = "horizontal",
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)
      ) +
      labs(title = unique(df$Title), x = NULL, y = NULL, color = NULL)
  }

  p1 <- plot_panel(df_bef_b, batch_cols)
  p2 <- plot_panel(df_aft_b, batch_cols)
  p3 <- plot_panel(df_bef_t, tissue_cols)
  p4 <- plot_panel(df_aft_t, tissue_cols)

  final_plot <- (p1 | p2) / (p3 | p4) +
    patchwork::plot_annotation(
      title = "Microarray (GEO-like) — UMAP Before/After",
      subtitle = if (run_combat) {
        "Top: Batch | Bottom: Tissue (Normal/Tumor/Unknown) | After = ComBat (multi-batches only)"
      } else {
        "Top: Batch | Bottom: Tissue (Normal/Tumor/Unknown) | After = no ComBat (single/one batch)"
      }
    )

  print(final_plot)
  ggsave(file.path(base_dir, "GEO_UMAP_4panel.png"), final_plot, width = 12, height = 10, dpi = 300)
} else {
  message(sprintf("Only %d Microarray samples — skipping UMAP.", ncol(expr_geo_cb)))
}

# ----------------- SAVE OUTPUTS ------------------
# Batch × Tissue table (includes Unknown)
tab_bt <- table(batch_geo, tissue3)
write.csv(as.data.frame(tab_bt), file.path(base_dir, "GEO_batchXtissue_table.csv"), row.names = FALSE)

write.csv(expr_geo_cb, file.path(base_dir, "GEO_ComBat_matrix.csv"), quote = FALSE)
message("Saved GEO_ComBat_matrix.csv and GEO_batchXtissue_table.csv")
