# 03_rnaseq_tcga_only.R
# TCGA-only ComBat (requires TCGA normals; stops if absent)
suppressPackageStartupMessages({
  library(dplyr); library(sva); library(matrixStats)
  library(umap); library(ggplot2); library(viridis); library(patchwork)
  library(readr)
})

base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
expr_file <- file.path(base_dir, "Full_matrix.csv")
cov_file  <- file.path(base_dir, "CombinedCovariates.csv")
stopifnot(file.exists(expr_file), file.exists(cov_file))

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cov  <- read.csv(cov_file, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), cov$Sample))

keep <- cov$Study == "TCGA" & !is.na(cov$Tissue)
expr_tcga <- expr[, keep, drop=FALSE]
cov_tcga  <- cov[keep, , drop=FALSE]
tissue <- factor(cov_tcga$Tissue, levels=c("Normal","Tumor"))

# Need both tissues
if (length(unique(na.omit(tissue))) < 2) {
  stop("TCGA lacks both tissue classes here. Add TCGA normals or skip this path.")
}

# Batch within TCGA: use whatever is available; fallback to one batch (no correction)
batch <- factor(rep("TCGA", ncol(expr_tcga)))
# If you have a 'Center' or 'Plate' in Demo.xlsx, replace batch with that.

expr_tcga <- expr_tcga[rowVars(as.matrix(expr_tcga)) > 0, , drop=FALSE]
mod <- model.matrix(~ tissue)

expr_tcga_cb <- if (nlevels(batch) > 1) {
  ComBat(dat=as.matrix(expr_tcga), batch=batch, mod=mod,
         par.prior=TRUE, mean.only=FALSE)
} else {
  as.matrix(expr_tcga) # nothing to correct
}

# UMAP
set.seed(42)
sel <- rowVars(expr_tcga_cb); topn <- min(3000, sum(sel>0))
expr_sel <- expr_tcga_cb[order(sel, decreasing=TRUE)[1:topn], ]
um_cfg <- umap.defaults; um_cfg$n_neighbors <- 15; um_cfg$min_dist <- 0.2; um_cfg$random_state <- 42
um_bef <- umap::umap(t(expr_tcga),     config=um_cfg)$layout
um_aft <- umap::umap(t(expr_sel),      config=um_cfg)$layout

df_bef_t <- data.frame(UMAP1=um_bef[,1], UMAP2=um_bef[,2], Color=tissue, Title="Before (Tissue)")
df_aft_t <- data.frame(UMAP1=um_aft[,1], UMAP2=um_aft[,2], Color=tissue, Title="After (Tissue)")
tissue_cols <- c(Normal="#7f7f7f", Tumor="#d62728")

plot_panel <- function(df){
  ggplot(df, aes(UMAP1, UMAP2, color = Color)) +
    geom_point(size=2.4, alpha=0.85) +
    scale_color_manual(values = tissue_cols) +
    theme_bw(base_size=13) +
    theme(legend.position="bottom", legend.box="horizontal",
          panel.grid.major = element_line(colour="grey90"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face="bold", hjust=0.5)) +
    labs(title=unique(df$Title), x=NULL, y=NULL, color=NULL)
}
final_plot <- (plot_panel(df_bef_t) | plot_panel(df_aft_t)) +
  patchwork::plot_annotation(title="TCGA — UMAP Before/After (Tissue)")
print(final_plot)
ggsave(file.path(base_dir, "TCGA_UMAP_Tissue.png"), final_plot, width=10, height=6, dpi=300)

write.csv(expr_tcga_cb, file.path(base_dir, "TCGA_ComBat_matrix.csv"), quote=FALSE)
cat("Saved TCGA_ComBat_matrix.csv and TCGA_UMAP_Tissue.png\n")
