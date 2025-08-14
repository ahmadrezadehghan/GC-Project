# 01_geo_combat_umap.R
# Batch correct GEO (microarray) using ComBat and plot UMAP/PCA

suppressPackageStartupMessages({
  library(dplyr); library(sva); library(matrixStats)
  library(umap); library(ggplot2); library(viridis)
  library(readxl); library(patchwork); library(readr)
})

# -------------------- CONFIG --------------------
base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
expr_file <- file.path(base_dir, "Full_matrix.csv")
cov_file  <- file.path(base_dir, "CombinedCovariates.csv")
demo_file <- file.path(base_dir, "Demo.xlsx")
stopifnot(file.exists(expr_file), file.exists(cov_file), file.exists(demo_file))

# -------------------- LOAD ----------------------
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cov  <- read.csv(cov_file, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), cov$Sample))

# GEO subset with defined tissue
keep <- cov$Study == "GEO" & !is.na(cov$Tissue)
expr_geo <- expr[, keep, drop=FALSE]
cov_geo  <- cov[keep, , drop=FALSE]
stopifnot(ncol(expr_geo) >= 10)

# ----------------- DEMO HELPERS ------------------
norm_key <- function(x) tolower(gsub("\\s+"," ", trimws(as.character(x))))

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
  hit  <- which(features_key %in% cand)
  if (length(hit) >= 1) hit[1] else NA_integer_
}

# try multiple options for organization / country / series
idx_org   <- find_row(c("Organization", "Organisation", "Organization Name",
                        "Institute", "Institution", "Center", "Centre",
                        "Laboratory", "Lab", "Submitter Organization"))
idx_ctry  <- find_row(c("Country", "Country ", "Origin Country", "Submitter Country",
                        "Location", "Nation"))
idx_series<- find_row(c("Series", "GSE", "Series ID", "GEO Series", "Study Accession",
                        "Study", "Project", "Project ID"))

# extract vectors (may be NA)
org_vec    <- if (!is.na(idx_org))   demo_vals(idx_org)   else setNames(rep(NA_character_, ncol(demo)-1), names(demo)[-1])
country_vec<- if (!is.na(idx_ctry))  demo_vals(idx_ctry)  else setNames(rep(NA_character_, ncol(demo)-1), names(demo)[-1])
series_vec <- if (!is.na(idx_series))demo_vals(idx_series)else setNames(rep(NA_character_, ncol(demo)-1), names(demo)[-1])

# align to GEO samples
organization <- org_vec[cov_geo$Sample]
country      <- country_vec[cov_geo$Sample]
series       <- series_vec[cov_geo$Sample]

organization <- trimws(ifelse(is.na(organization), "NA", as.character(organization)))
country      <- trimws(ifelse(is.na(country),      "NA", as.character(country)))
series       <- trimws(ifelse(is.na(series),       "NA", as.character(series)))

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
    batch_geo <- factor(rep("GEO", nrow(cov_geo)))  # single batch; ComBat will be skipped
  }
}

# Drop tiny batches (<2) to stabilize ComBat estimates
tb <- table(batch_geo)
ok_batch <- batch_geo %in% names(tb)[tb >= 2]
expr_geo <- expr_geo[, ok_batch, drop = FALSE]
cov_geo  <- cov_geo[ok_batch, , drop = FALSE]
batch_geo <- droplevels(batch_geo[ok_batch])

# ----------------- BASIC FILTERS -----------------
# zero-variance genes
expr_geo <- expr_geo[rowVars(as.matrix(expr_geo)) > 0, , drop = FALSE]

tissue <- factor(cov_geo$Tissue, levels = c("Normal","Tumor"))
stopifnot(length(unique(na.omit(tissue))) >= 2)

# ----------------- REPORT ------------------------
cat("Contingency (Batch × Tissue) AFTER tiny-batch removal:\n")
tab_bt <- table(batch_geo, tissue)
print(tab_bt)
write.csv(as.data.frame(tab_bt), file.path(base_dir, "GEO_batchXtissue_table.csv"), row.names = FALSE)

# ----------------- COMBAT (conditional) ----------
run_combat <- (nlevels(batch_geo) > 1)
if (run_combat) {
  mod <- model.matrix(~ tissue)
  cat("Running ComBat (mean & variance)…\n")
  expr_geo_cb <- ComBat(dat = as.matrix(expr_geo),
                        batch = batch_geo, mod = mod,
                        par.prior = TRUE, prior.plots = FALSE,
                        mean.only = FALSE)
  cat("ComBat done.\n")
} else {
  cat("Only one batch detected; skipping ComBat and using original expression.\n")
  expr_geo_cb <- as.matrix(expr_geo)
}

# ----------------- UMAP & PCA --------------------
set.seed(42)
# Before-correction UMAP on all genes (or a cap)
um_cfg <- umap.defaults; um_cfg$n_neighbors <- 15; um_cfg$min_dist <- 0.2; um_cfg$random_state <- 42
um_bef <- umap::umap(t(expr_geo), config = um_cfg)$layout

# After-correction UMAP on top variable genes
sel <- rowVars(expr_geo_cb); topn <- min(3000, sum(sel > 0))
expr_sel <- expr_geo_cb[order(sel, decreasing = TRUE)[1:topn], , drop = FALSE]
um_aft <- umap::umap(t(expr_sel), config = um_cfg)$layout

df_bef_b <- data.frame(UMAP1=um_bef[,1], UMAP2=um_bef[,2], Color=batch_geo, Title="Before (Batch)")
df_aft_b <- data.frame(UMAP1=um_aft[,1], UMAP2=um_aft[,2], Color=batch_geo, Title="After (Batch)")
df_bef_t <- data.frame(UMAP1=um_bef[,1], UMAP2=um_bef[,2], Color=tissue,    Title="Before (Tissue)")
df_aft_t <- data.frame(UMAP1=um_aft[,1], UMAP2=um_aft[,2], Color=tissue,    Title="After (Tissue)")

# palettes (handle 1-level case safely)
batch_levels <- levels(batch_geo)
if (is.null(batch_levels)) batch_levels <- "Batch"
batch_cols <- viridis(max(1, length(batch_levels))); names(batch_cols) <- batch_levels
tissue_cols <- c(Normal="#7f7f7f", Tumor="#d62728")

plot_panel <- function(df, pal){
  ggplot(df, aes(UMAP1, UMAP2, color = Color)) +
    geom_point(size=2.4, alpha=0.85) +
    {
      if (is.factor(df$Color) || is.character(df$Color)) {
        scale_color_manual(values = pal)
      } else {
        scale_color_viridis_c()
      }
    } +
    theme_bw(base_size=13) +
    theme(legend.position="bottom", legend.box="horizontal",
          panel.grid.major = element_line(colour="grey90"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face="bold", hjust=0.5)) +
    labs(title=unique(df$Title), x=NULL, y=NULL, color=NULL)
}

p1 <- plot_panel(df_bef_b, batch_cols)
p2 <- plot_panel(df_aft_b, batch_cols)
p3 <- plot_panel(df_bef_t, tissue_cols)
p4 <- plot_panel(df_aft_t, tissue_cols)

final_plot <- (p1 | p2) / (p3 | p4) +
  patchwork::plot_annotation(
    title   = "GEO (Microarray) — UMAP Before/After",
    subtitle= if (run_combat) "Top: Batch | Bottom: Tissue | After = ComBat-corrected"
    else "Top: Batch | Bottom: Tissue | After = no ComBat (single batch)"
  )

print(final_plot)
ggsave(file.path(base_dir, "GEO_UMAP_4panel.png"), final_plot, width=12, height=10, dpi=300)

# Save corrected (or original) matrix used for 'After'
write.csv(expr_geo_cb, file.path(base_dir, "GEO_ComBat_matrix.csv"), quote = FALSE)
cat("Saved GEO_ComBat_matrix.csv, GEO_UMAP_4panel.png, and GEO_batchXtissue_table.csv\n")
