# 13_diffcoexp.R
suppressPackageStartupMessages({
    library(DGCA) # CRAN: DGCA
    library(dplyr)
    library(readr)
    library(matrixStats)
    library(igraph)
})

base_dir <- Sys.getenv("GC_BASE_DIR", unset = "E:/1----------GC_29_4_mainDATAFRAME---------1")
expr_file <- file.path(base_dir, "All_ComBat_matrix.csv")
cov_file <- file.path(base_dir, "All_Covariates_augmented.csv")
out_dir <- file.path(base_dir, "Results_DiffCoExp")
dir.create(out_dir, TRUE, TRUE)
stopifnot(file.exists(expr_file), file.exists(cov_file))

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE) %>% as.matrix()
mode(expr) <- "numeric"
cov <- read.csv(cov_file, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), cov$Sample))

# --- ژن‌های متغیر برای کاهش محاسبات ---
topN <- min(1200L, nrow(expr))
rv <- matrixStats::rowVars(expr)
keep <- order(rv, decreasing = TRUE)[seq_len(topN)]
X <- expr[keep, , drop = FALSE]

classes <- factor(ifelse(tolower(cov$Tissue) == "tumor", "Tumor",
    ifelse(tolower(cov$Tissue) == "normal", "Normal", NA)
))
ok <- !is.na(classes)
X <- X[, ok, drop = FALSE]
classes <- droplevels(classes[ok])

# DGCA
design <- model.matrix(~ 0 + classes) # گروه‌ها
colnames(design) <- levels(classes)
res <- ddcorAll(
    inputMat = X, design = design, compare = c("Tumor", "Normal"),
    corrType = "spearman", nPairs = 1e7, nPerm = 0, adjust = "BH", splitSet = NULL
)
# significant differential correlations
sig <- res %>% filter(q.adjusted < 0.05)
write_csv(res, file.path(out_dir, "DGCA_all_pairs.csv"))
write_csv(sig, file.path(out_dir, "DGCA_sig_pairs_q05.csv"))

# ماژول‌های DiffCoExp: گراف از زوج‌های معنادار
if (nrow(sig) >= 2) {
    g <- graph_from_data_frame(sig[, c("Gene1", "Gene2")], directed = FALSE)
    cmp <- components(g)$membership
    mod_df <- data.frame(Gene = names(cmp), DiffCoExModule = paste0("DCM", cmp))
    write_csv(mod_df, file.path(out_dir, "DiffCoEx_modules.csv"))
}
cat("DGCA DiffCoExp done.\n")
