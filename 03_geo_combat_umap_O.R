# 03_geo_combat_umap_O.R
# Batch correct GEO (microarray) using ComBat and plot UMAP
# Removes NA/Unknown/Other tissue samples from all outputs

suppressPackageStartupMessages({
  library(dplyr)
  library(sva)
  library(matrixStats)
  library(umap)
  library(ggplot2)
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

# ------- SELECT MICROARRAY SUBSET (instead of cov$Study) -------
# Previous bug: keep_geo <- cov$Study == "GEO"
keep_geo <- cov$Platform == "Microarray"
if (!any(keep_geo)) stop("No Microarray samples found (Platform == 'Microarray'). Check CombinedCovariates.csv")

expr_geo <- expr[, keep_geo, drop = FALSE]
cov_geo <- cov[keep_geo, , drop = FALSE]
stopifnot(ncol(expr_geo) >= 3)

# ----------------- DEMO HELPERS ------------------
norm_key <- function(x) tolower(gsub("\\s+", " ", trimws(as.character(x))))

demo <- as.data.frame(readxl::read_excel(demo_file))
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
  "Organization", "Organisation", "Organization Name", "Institute", "Institution",
  "Center", "Centre", "Laboratory", "Lab", "Submitter Organization"
))
idx_ctry <- find_row(c("Country", "Country ", "Origin Country", "Submitter Country", "Location", "Nation"))
idx_series <- find_row(c("Series", "GSE", "Series ID", "GEO Series", "Study Accession", "Study", "Project", "Project ID"))

mk_na_vec <- function() setNames(rep(NA_character_, ncol(demo) - 1), names(demo)[-1])
org_vec <- if (!is.na(idx_org)) demo_vals(idx_org) else mk_na_vec()
country_vec <- if (!is.na(idx_ctry)) demo_vals(idx_ctry) else mk_na_vec()
series_vec <- if (!is.na(idx_series)) demo_vals(idx_series) else mk_na_vec()

# align to GEO samples
organization <- trimws(ifelse(is.na(org_vec[cov_geo$Sample]), "NA", as.character(org_vec[cov_geo$Sample])))
country <- trimws(ifelse(is.na(country_vec[cov_geo$Sample]), "NA", as.character(country_vec[cov_geo$Sample])))
series <- trimws(ifelse(is.na(series_vec[cov_geo$Sample]), "NA", as.character(series_vec[cov_geo$Sample])))

# ----------------- BUILD BATCH -------------------
# Priority: Series/GSE; else Organization×Country; else Platform; else single batch
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
    batch_geo <- factor(rep("GEO", nrow(cov_geo))) # single batch
  }
}

# ----------------- BASIC FILTERS -----------------
# zero-variance genes
expr_geo <- expr_geo[rowVars(as.matrix(expr_geo)) > 0, , drop = FALSE]

# Tissue; drop NA/Unknown/Other samples entirely from all downstream outputs
tissue_raw <- cov_geo$Tissue
tissue_norm <- tolower(trimws(as.character(tissue_raw)))
is_bad_tissue <- is.na(tissue_raw) |
  tissue_norm %in% c("unknown", "other", "na", "") |
  nchar(tissue_norm) == 0
# Keep only Normal / Tumor
keep_tissue <- (!is_bad_tissue) & tissue_norm %in% c("normal", "tumor")
expr_geo <- expr_geo[, keep_tissue, drop = FALSE]
cov_geo <- cov_geo[keep_tissue, , drop = FALSE]
batch_geo <- droplevels(batch_geo[keep_tissue])

# final tissue factor (no Unknown/Other)
tissue2 <- factor(ifelse(tolower(cov_geo$Tissue) == "normal", "Normal", "Tumor"),
  levels = c("Normal", "Tumor")
)

if (ncol(expr_geo) < 3) stop("Fewer than 3 samples remain after tissue filtering; cannot run UMAP.")

# Drop tiny batches (<2) AFTER tissue filtering (stabilizes ComBat)
tb <- table(batch_geo)
ok_batch <- batch_geo %in% names(tb)[tb >= 2]
expr_geo <- expr_geo[, ok_batch, drop = FALSE]
cov_geo <- cov_geo[ok_batch, , drop = FALSE]
batch_geo <- droplevels(batch_geo[ok_batch])
tissue2 <- droplevels(tissue2[ok_batch])

if (ncol(expr_geo) < 3) stop("Fewer than 3 samples remain after tiny-batch removal; cannot run UMAP.")

# ----------------- REPORT ------------------------
cat("Contingency (Batch × Tissue) AFTER filtering and tiny-batch removal:\n")
tab_bt <- table(batch_geo, tissue2)
print(tab_bt)
write.csv(as.data.frame(tab_bt),
  file.path(base_dir, "GEO_batchXtissue_table.csv"),
  row.names = FALSE
)

# ----------------- COMBAT (adaptive design) ----------
run_combat <- (nlevels(batch_geo) > 1)
if (run_combat) {
  # With only Normal/Tumor left, use tissue if both present
  present <- levels(droplevels(tissue2))
  if (length(present) >= 2) {
    mod <- model.matrix(~tissue2)
  } else {
    mod <- model.matrix(~1, data = data.frame(dummy = rep(1, ncol(expr_geo))))
  }
  cat("Running ComBat (mean & variance)…\n")
  expr_geo_cb <- ComBat(
    dat = as.matrix(expr_geo),
    batch = batch_geo, mod = mod,
    par.prior = TRUE, prior.plots = FALSE,
    mean.only = FALSE
  )
  cat("ComBat done.\n")
} else {
  cat("Only one batch detected; skipping ComBat and using original expression.\n")
  expr_geo_cb <- as.matrix(expr_geo)
}

# ----------------- UMAP --------------------------
set.seed(42)
um_cfg <- umap.defaults
um_cfg$n_neighbors <- min(15, ncol(expr_geo) - 1)
um_cfg$min_dist <- 0.2
um_cfg$random_state <- 42

um_bef <- umap::umap(t(expr_geo), config = um_cfg)$layout

sel <- rowVars(expr_geo_cb)
topn <- min(3000, max(2, sum(sel > 0)))
expr_sel <- expr_geo_cb[order(sel, decreasing = TRUE)[1:topn], , drop = FALSE]
um_aft <- umap::umap(t(expr_sel), config = um_cfg)$layout

df_bef_b <- data.frame(UMAP1 = um_bef[, 1], UMAP2 = um_bef[, 2], Color = batch_geo, Title = "Before (Batch)")
df_aft_b <- data.frame(UMAP1 = um_aft[, 1], UMAP2 = um_aft[, 2], Color = batch_geo, Title = "After (Batch)")
df_bef_t <- data.frame(UMAP1 = um_bef[, 1], UMAP2 = um_bef[, 2], Color = tissue2, Title = "Before (Tissue)")
df_aft_t <- data.frame(UMAP1 = um_aft[, 1], UMAP2 = um_aft[, 2], Color = tissue2, Title = "After (Tissue)")

# ---------- COLOR PALETTES ----------
nice_cols <- c(
  "#2ca02c", "#1f77b4", "#d62728", "#ff7f0e",
  "#9467bd", "#8c564b", "#e377c2", "#17becf",
  "#bcbd22", "#7f7f7f", "#aec7e8", "#ff9896",
  "#98df8a", "#c5b0d5", "#c49c94", "#f7b6d2", "#9edae5", "#dbdb8d"
)
batch_levels <- levels(batch_geo)
if (is.null(batch_levels)) batch_levels <- "Batch"
batch_cols <- setNames(rep(nice_cols, length.out = length(batch_levels)), batch_levels)

# Bottom (Tissue): Normal = green, Tumor = red
tissue_cols <- c(Normal = "#2ca02c", Tumor = "#d62728")

plot_panel <- function(df, pal) {
  ggplot(df, aes(UMAP1, UMAP2, color = Color)) +
    geom_point(size = 2.2, alpha = 0.85) +
    {
      if (is.factor(df$Color) || is.character(df$Color)) {
        scale_color_manual(values = pal, drop = TRUE)
      } else {
        scale_color_gradient()
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
    title = "Microarray subset — UMAP Before/After",
    subtitle = if (run_combat) {
      "Top: Batch | Bottom: Tissue (Normal/Tumor) | After = ComBat-corrected"
    } else {
      "Top: Batch | Bottom: Tissue (Normal/Tumor) | After = no ComBat (single batch)"
    }
  )

print(final_plot)
ggsave(file.path(base_dir, "GEO_UMAP_4panel.png"), final_plot, width = 12, height = 10, dpi = 300)

# ----------------- SAVE OUTPUTS ------------------
# Table includes only Normal/Tumor; matrix includes only kept samples
write.csv(expr_geo_cb, file.path(base_dir, "GEO_ComBat_matrix.csv"), quote = FALSE)
cat("Saved GEO_ComBat_matrix.csv, GEO_UMAP_4panel.png, and GEO_batchXtissue_table.csv\n")
