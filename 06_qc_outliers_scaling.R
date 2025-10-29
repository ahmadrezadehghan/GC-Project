# 06_qc_outliers_scaling.R
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(matrixStats)
    library(ggplot2)
    library(umap)
    library(rrcov)
    library(patchwork)
})

# -------------------- CONFIG --------------------
base_dir <- Sys.getenv("GC_BASE_DIR", unset = "E:/1----------GC_29_4_mainDATAFRAME---------1")
full_file <- file.path(base_dir, "Full_matrix.csv")
post_file <- file.path(base_dir, "All_ComBat_matrix.csv")
cov_file <- file.path(base_dir, "All_Covariates_augmented.csv")
out_dir <- file.path(base_dir, "Results_QC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(post_file), file.exists(cov_file))
if (!file.exists(full_file)) message("⚠️ Full_matrix.csv یافت نشد؛ فقط مرحله‌ی پس از ComBat را QC می‌کنم.")

# -------------------- HELPERS -------------------
safe_log1p <- function(m, eps = 1e-6) {
    m <- as.matrix(m)
    mode(m) <- "numeric"
    mn <- suppressWarnings(min(m, na.rm = TRUE))
    if (is.finite(mn) && mn <= -1) m <- m - mn + eps
    log1p(m)
}
pca_plot <- function(X, lab, color = NULL) {
    pr <- prcomp(t(X), scale. = TRUE)
    pc <- as.data.frame(pr$x[, 1:2])
    pc$lab <- lab
    if (!is.null(color)) pc$color <- color
    gg <- ggplot(pc, aes(PC1, PC2, color = if (!is.null(color)) color else lab)) +
        geom_point(alpha = .85, size = 2) +
        theme_bw(base_size = 12) +
        labs(title = paste0("PCA: ", lab), color = NULL)
    list(plot = gg, obj = pr)
}
umap_plot <- function(X, lab, color = NULL) {
    cfg <- umap.defaults
    cfg$n_neighbors <- min(15, ncol(X) - 1)
    cfg$min_dist <- 0.2
    um <- umap::umap(t(X), config = cfg)$layout
    df <- data.frame(UMAP1 = um[, 1], UMAP2 = um[, 2], color = if (is.null(color)) lab else color)
    ggplot(df, aes(UMAP1, UMAP2, color = color)) +
        geom_point(alpha = .85, size = 2) +
        theme_bw(base_size = 12) +
        labs(title = paste0("UMAP: ", lab), color = NULL)
}
dens_plot <- function(X, title) {
    Xl <- safe_log1p(X)
    df <- reshape2::melt(as.data.frame(Xl))
    colnames(df) <- c("Gene", "Value")
    ggplot(df, aes(Value, group = Gene)) +
        geom_density(alpha = 0.05, color = "grey30") +
        theme_bw(base_size = 12) +
        labs(title = title, x = "log1p(expr)", y = "Density")
}
box_plot <- function(X, title) {
    Xl <- safe_log1p(X)
    sm <- apply(Xl, 2, function(v) quantile(v, probs = c(.25, .5, .75), na.rm = TRUE))
    df <- data.frame(Sample = colnames(Xl), Q1 = sm[1, ], Med = sm[2, ], Q3 = sm[3, ])
    ggplot(df, aes(x = reorder(Sample, Med), ymin = Q1, lower = Q1, middle = Med, upper = Q3, ymax = Q3)) +
        geom_boxplot(stat = "identity") +
        coord_flip() +
        theme_bw(base_size = 11) +
        labs(title = title, x = NULL, y = "log1p(expr) quantiles")
}

detect_outliers <- function(X, k_pcs = 5, iqr_mult = 3) {
    Xl <- t(scale(t(safe_log1p(X))))
    Xl[is.na(Xl)] <- 0
    # robust PCA
    rp <- try(rrcov::PcaHubert(t(Xl), k = min(k_pcs, ncol(Xl) - 1)), silent = TRUE)
    rd <- rep(NA_real_, ncol(Xl))
    if (!inherits(rp, "try-error")) {
        rd <- rp@distances
        cut <- sqrt(qchisq(0.975, df = rp@k))
        r_out <- rd > cut
    } else {
        r_out <- rep(FALSE, ncol(Xl))
    }
    # IQR on library size (sum) و میانه‌ی سیگنال
    lib <- colSums(safe_log1p(X))
    med <- colMedians(safe_log1p(X))
    flags <- function(v) {
        q <- quantile(v, probs = c(.25, .75), na.rm = TRUE)
        iqr <- q[2] - q[1]
        v < (q[1] - iqr_mult * iqr) | v > (q[2] + iqr_mult * iqr)
    }
    iqr_out <- flags(lib) | flags(med)
    data.frame(
        Sample = colnames(X),
        RPCA_Outlier = r_out,
        IQR_Outlier = iqr_out,
        LibSize = lib, Median = med,
        RobustDistance = rd,
        Any_Outlier = r_out | iqr_out,
        stringsAsFactors = FALSE
    )
}

# -------------------- LOAD ----------------------
post <- read.csv(post_file, row.names = 1, check.names = FALSE) %>% as.matrix()
mode(post) <- "numeric"
cov <- read.csv(cov_file, stringsAsFactors = FALSE)
stopifnot(identical(colnames(post), cov$Sample))
pre <- NULL
if (file.exists(full_file)) {
    pre <- read.csv(full_file, row.names = 1, check.names = FALSE) %>% as.matrix()
    mode(pre) <- "numeric"
    pre <- pre[, colnames(post), drop = FALSE] # align to same samples
}

# -------------------- QC PLOTS ------------------
if (!is.null(pre)) {
    g1 <- dens_plot(pre, "Density (pre-ComBat)")
    g2 <- dens_plot(post, "Density (post-ComBat)")
    g <- g1 + g2 + plot_layout(ncol = 2)
    ggsave(file.path(out_dir, "QC_density_pre_post.png"), g, width = 12, height = 5, dpi = 150)
    ggsave(file.path(out_dir, "QC_box_pre.png"), box_plot(pre, "Box (pre)"), width = 10, height = 10, dpi = 150)
}
ggsave(file.path(out_dir, "QC_box_post.png"), box_plot(post, "Box (post)"), width = 10, height = 10, dpi = 150)

# PCA/UMAP
if (!is.null(pre)) {
    p1 <- pca_plot(pre, "pre-ComBat", color = cov$Tissue)$plot
    p2 <- pca_plot(post, "post-ComBat", color = cov$Tissue)$plot
    u1 <- umap_plot(pre, "pre-ComBat", color = cov$Tissue)
    u2 <- umap_plot(post, "post-ComBat", color = cov$Tissue)
    ggsave(file.path(out_dir, "PCA_pre_post.png"), (p1 | p2), width = 12, height = 5, dpi = 150)
    ggsave(file.path(out_dir, "UMAP_pre_post.png"), (u1 | u2), width = 12, height = 5, dpi = 150)
} else {
    p2 <- pca_plot(post, "post-ComBat", color = cov$Tissue)$plot
    u2 <- umap_plot(post, "post-ComBat", color = cov$Tissue)
    ggsave(file.path(out_dir, "PCA_post.png"), p2, width = 6, height = 5, dpi = 150)
    ggsave(file.path(out_dir, "UMAP_post.png"), u2, width = 6, height = 5, dpi = 150)
}

# -------------------- OUTLIERS ------------------
ol <- detect_outliers(post, k_pcs = 5, iqr_mult = 3)
write.csv(ol, file.path(out_dir, "Outliers_report.csv"), row.names = FALSE)
keep <- ol$Sample[!ol$Any_Outlier]
drop <- ol$Sample[ol$Any_Outlier]
writeLines(keep, file.path(out_dir, "retained_samples.txt"))
writeLines(drop, file.path(out_dir, "removed_samples.txt"))
cat(sprintf(
    "QC done. Kept %d / %d samples (outliers flagged: %d)\n",
    length(keep), nrow(ol), length(drop)
))
