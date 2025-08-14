# 02_all_batches_combat_umap.R
# Attempt global ComBat; STOP if any batch lacks both tissues (confounding)
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

# Build batch = Study × Platform (safe default)
batch <- interaction(cov$Study, cov$Platform, drop=TRUE)
tissue <- factor(cov$Tissue, levels=c("Normal", "Tumor"))

# Keep only samples with known tissue
keep <- !is.na(tissue)

# Filter expr, batch, tissue, and cov together to ensure all have the same length
expr <- expr[, keep, drop=FALSE]
batch <- batch[keep]
tissue <- tissue[keep]
cov <- cov[keep, , drop=FALSE]

# Drop tiny batches (<2)
tb <- table(batch)
keep2 <- batch %in% names(tb)[tb >= 2]
expr <- expr[, keep2, drop=FALSE]
batch <- batch[keep2]
tissue <- tissue[keep2]
cov <- cov[keep2, , drop=FALSE]

# Check and exclude TCGA samples if necessary
tcga_tissue <- cov$Tissue[cov$Study == "TCGA"]
if (length(unique(tcga_tissue)) == 1) {
  cat("Skipping TCGA as it lacks both Tumor and Normal tissue.\n")
  
  # Exclude TCGA from the analysis
  cov <- cov[cov$Study != "TCGA", ]  # Remove TCGA from cov
  expr <- expr[, colnames(expr) %in% cov$Sample, drop=FALSE]  # Align expr with updated cov
  batch <- batch[cov$Study != "TCGA"]
  tissue <- tissue[cov$Study != "TCGA"]
}

# Exclude batches where only Normal samples are present (e.g., GTEX.RNAseq)
batch_table <- table(batch, tissue)
confounded_batches <- rowSums(batch_table > 0) == 1  # Check if any batch has only 1 tissue level
if (any(confounded_batches)) {
  bad_batches <- names(confounded_batches)[confounded_batches]
  cat("Excluding confounded batches:", paste(bad_batches, collapse=", "), "\n")
  
  # Exclude confounded batches (those with only Normal or Tumor)
  cov <- cov[!batch %in% bad_batches, ]  # Remove confounded batches from cov
  expr <- expr[, colnames(expr) %in% cov$Sample, drop=FALSE]  # Align expr with updated cov
  batch <- batch[!batch %in% bad_batches]  # Update batch
  tissue <- tissue[!batch %in% bad_batches]  # Update tissue
}

# After exclusion, ensure that batch, tissue, and expr are properly aligned
expr <- expr[, colnames(expr) %in% cov$Sample]  # Ensure expr aligns with the updated cov
batch <- batch[colnames(expr) %in% cov$Sample]  # Align batch with the updated expr
tissue <- tissue[colnames(expr) %in% cov$Sample]  # Align tissue with the updated expr

# Ensure there are no NAs in batch, tissue, and expr after filtering
stopifnot(all(!is.na(batch)))
stopifnot(all(!is.na(tissue)))
stopifnot(all(!is.na(expr)))

# Confounding check: each batch must contain both levels
tab <- table(batch, tissue)
cat("Confounding check results:\n")
print(tab)

# This should now work without errors:
# You can now proceed with ComBat or other analyses as needed
