# 10_network_variants.R
suppressPackageStartupMessages({
    library(minet)
    library(infotheo)
    library(igraph)
    library(dplyr)
    library(readr)
    library(matrixStats)
    library(purrr)
    library(pheatmap)
    library(pROC)
    library(pcit) # CRAN: pcit
})

base_dir <- Sys.getenv("GC_BASE_DIR", unset = "E:/1----------GC_29_4_mainDATAFRAME---------1")
expr_file <- file.path(base_dir, "All_ComBat_matrix.csv")
cov_file <- file.path(base_dir, "All_Covariates_augmented.csv")
out_dir <- file.path(base_dir, "Results_Networks_Variants")
dir.create(out_dir, TRUE, TRUE)
stopifnot(file.exists(expr_file), file.exists(cov_file))

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE) %>% as.matrix()
mode(expr) <- "numeric"
cov <- read.csv(cov_file, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), cov$Sample))

# --- انتخاب ژن‌ها ---
topN <- min(2000L, nrow(expr))
rv <- matrixStats::rowVars(expr)
keep <- order(rv, decreasing = TRUE)[seq_len(topN)]
X <- expr[keep, , drop = FALSE]

groups <- list(
    Tumor  = which(tolower(cov$Tissue) == "tumor"),
    Normal = which(tolower(cov$Tissue) == "normal")
)

build_edges_minet <- function(Xsub, method = c("aracne", "clr", "mrnet"), bins = 3L) {
    disc <- infotheo::discretize(data.frame(t(Xsub)), nbins = bins)
    mim <- minet::build.mim(disc, estimator = "mi.empirical")
    adj <- switch(match.arg(method),
        "aracne" = minet::aracne(mim),
        "clr"    = minet::clr(mim),
        "mrnet"  = minet::mrnet(mim)
    )
    idx <- which(adj != 0, arr.ind = TRUE)
    idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
    if (nrow(idx) == 0) {
        return(NULL)
    }
    data.frame(
        g1 = rownames(adj)[idx[, 1]],
        g2 = colnames(adj)[idx[, 2]],
        w = mapply(function(i, j) mim[i, j], idx[, 1], idx[, 2]),
        stringsAsFactors = FALSE
    )
}
build_edges_pcit <- function(Xsub) {
    C <- cor(t(Xsub), method = "spearman", use = "pairwise.complete.obs")
    pc <- pcit(C)
    edges <- NULL
    for (i in 1:(nrow(C) - 1)) {
        for (j in (i + 1):ncol(C)) {
            if (PCIT::decision(pc, i, j)) {
                edges <- rbind(edges, data.frame(g1 = rownames(C)[i], g2 = colnames(C)[j], w = abs(C[i, j])))
            }
        }
    }
    edges
}

write_edges <- function(ed, method, grp) {
    if (is.null(ed) || nrow(ed) == 0) {
        return(NULL)
    }
    fn <- file.path(out_dir, paste0(method, "_", grp, "_edges.csv"))
    write_csv(ed, fn)
    fn
}

edge_key <- function(df) paste(pmin(df$g1, df$g2), pmax(df$g1, df$g2), sep = "|")

overlap_summary <- list()

for (g in names(groups)) {
    idx <- groups[[g]]
    if (length(idx) < 10) {
        message("Skip ", g, ": n<10")
        next
    }
    Xg <- X[, idx, drop = FALSE]
    E_ar <- build_edges_minet(Xg, "aracne")
    E_cl <- build_edges_minet(Xg, "clr")
    E_mr <- build_edges_minet(Xg, "mrnet")
    E_pc <- build_edges_pcit(Xg)

    files <- list(
        aracne = write_edges(E_ar, "aracne", g),
        clr    = write_edges(E_cl, "clr", g),
        mrnet  = write_edges(E_mr, "mrnet", g),
        pcit   = write_edges(E_pc, "pcit", g)
    )

    # همپوشانی
    keys <- lapply(list(aracne = E_ar, clr = E_cl, mrnet = E_mr, pcit = E_pc), function(ed) if (is.null(ed)) character(0) else unique(edge_key(ed)))
    meths <- names(keys)
    combs <- t(combn(meths, 2))
    ov <- apply(combs, 1, function(x) {
        a <- keys[[x[1]]]
        b <- keys[[x[2]]]
        inter <- length(intersect(a, b))
        uni <- length(union(a, b))
        data.frame(Group = g, A = x[1], B = x[2], Jaccard = ifelse(uni == 0, NA, inter / uni), Inter = inter, Union = uni)
    }) %>% bind_rows()
    overlap_summary[[g]] <- ov
}
if (length(overlap_summary)) {
    write_csv(bind_rows(overlap_summary), file.path(out_dir, "edge_overlap_jaccard.csv"))
}
cat("Network variants built.\n")
