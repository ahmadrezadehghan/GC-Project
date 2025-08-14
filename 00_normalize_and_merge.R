# 00_normalize_and_merge.R
# Build clean, aligned matrices + covariates from Demo.xlsx

suppressPackageStartupMessages({
  library(readxl)
  library(preprocessCore)
  library(matrixStats)
  library(dplyr)
  library(readr)
  library(stringr)
})

# -------------------- CONFIG --------------------
base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
cov_file        <- file.path(base_dir, "Demo.xlsx")
rna_expr_file   <- file.path(base_dir, "RNA_seq.xlsx")
micro_expr_file <- file.path(base_dir, "Micro.xlsx")
stopifnot(file.exists(cov_file), file.exists(rna_expr_file), file.exists(micro_expr_file))

# ----------------- UTILITIES --------------------
clean_ids <- function(df, id_col = 1){
  df <- df[!is.na(df[[id_col]]), , drop=FALSE]
  df <- df[!duplicated(df[[id_col]]), , drop=FALSE]
  rn <- df[[id_col]]
  df <- df[, -id_col, drop=FALSE]
  rownames(df) <- rn
  df
}

# Force expression columns numeric
force_numeric <- function(df, label = "expr") {
  coerce_col <- function(v) {
    if (is.numeric(v)) return(v)
    x <- as.character(v)
    x <- gsub(",", "", x)             # remove commas
    x <- gsub("[[:space:]]+", "", x)  # remove spaces
    x[x %in% c("", "NA", "NaN", "Inf", "-Inf", "—", "–", "-", ".")] <- NA
    suppressWarnings(as.numeric(x))
  }
  num_df <- as.data.frame(lapply(df, coerce_col), check.names = FALSE)
  
  # Drop columns with <2 non-NA
  keep_col <- colSums(!is.na(num_df)) >= 2
  if (any(!keep_col)) {
    message(sprintf("[%s] Dropping %d non-numeric/empty columns: %s",
                    label, sum(!keep_col),
                    paste(names(num_df)[!keep_col], collapse = ", ")))
  }
  num_df <- num_df[, keep_col, drop = FALSE]
  
  # Drop rows with all NA
  keep_row <- rowSums(!is.na(num_df)) > 0
  if (any(!keep_row)) {
    message(sprintf("[%s] Dropping %d genes with all NA after coercion.", label, sum(!keep_row)))
  }
  num_df <- num_df[keep_row, , drop = FALSE]
  
  colnames(num_df) <- make.unique(colnames(num_df))
  num_df
}

impute_row_median <- function(mat){
  t(apply(mat,1,function(r){
    if(anyNA(r)) r[is.na(r)] <- median(r, na.rm=TRUE)
    r
  }))
}

compute_ref <- function(mat){
  sm <- apply(mat, 2, sort, na.last = NA)
  rowMeans(sm, na.rm = TRUE)
}

zscore_rows <- function(df){
  m <- as.matrix(df)
  mu <- rowMeans(m)
  sdv <- matrixStats::rowSds(m)
  sdv[sdv==0] <- 1
  as.data.frame((m - mu)/sdv, check.names = FALSE)
}

infer_study <- function(s) {
  ifelse(grepl("^GTEX", s, ignore.case=TRUE), "GTEX",
         ifelse(grepl("^TCGA", s, ignore.case=TRUE), "TCGA",
                ifelse(grepl("^GSM",  s, ignore.case=TRUE), "GEO", "OTHER")))
}

infer_platform <- function(s) {
  ifelse(grepl("RNA|TCGA|GTEX", s, ignore.case=TRUE), "RNAseq", "Microarray")
}

# ----------------- LOAD & NORMALIZE -----------------
message("Reading RNA-seq and Microarray…")
rna_raw   <- readxl::read_excel(rna_expr_file)   %>% as.data.frame()
micro_raw <- readxl::read_excel(micro_expr_file) %>% as.data.frame()

# Clean IDs and force numeric
rna_df <- clean_ids(rna_raw,   1) %>% force_numeric("RNA-seq")
mic_df <- clean_ids(micro_raw, 1) %>% force_numeric("Microarray")

# Convert to matrices
rna_mat <- as.matrix(rna_df)
mic_mat <- as.matrix(mic_df)

# log2 + impute + quantile normalize
rna_log <- log2(rna_mat + 1)
mic_log <- log2(mic_mat + 1)

rna_imp <- impute_row_median(rna_log)
mic_imp <- impute_row_median(mic_log)

rna_qn <- normalize.quantiles(rna_imp, keep.names = TRUE); dimnames(rna_qn) <- dimnames(rna_imp)
mic_qn <- normalize.quantiles(mic_imp, keep.names = TRUE); dimnames(mic_qn) <- dimnames(mic_imp)

# Rescale RNA to microarray target
micro_target <- compute_ref(mic_qn)
rna_to_micro <- normalize.quantiles.use.target(rna_qn, target = micro_target)
dimnames(rna_to_micro) <- dimnames(rna_qn)

norm_rna   <- as.data.frame(rna_to_micro, check.names = FALSE)
norm_micro <- as.data.frame(mic_qn,       check.names = FALSE)

# Z-scores for integration
norm_rna_z   <- zscore_rows(norm_rna)
norm_micro_z <- zscore_rows(norm_micro)

# Common genes only
common_genes <- intersect(rownames(norm_rna_z), rownames(norm_micro_z))
stopifnot(length(common_genes) > 1000)
norm_rna_z   <- norm_rna_z[common_genes, , drop=FALSE]
norm_micro_z <- norm_micro_z[common_genes, , drop=FALSE]
expr_merged  <- cbind(norm_rna_z, norm_micro_z)

# ----------------- DEMO.xlsx: TISSUE CLEAN -----------------
demo <- readxl::read_excel(cov_file) %>% as.data.frame()
features <- demo[[1]]

grab_row <- function(name) {
  stopifnot(name %in% features)
  r <- demo[features==name, -1, drop=TRUE]
  names(r) <- names(demo)[-1]
  r
}

trow <- grab_row("Tissue")

clean_tissue <- function(x){
  s <- tolower(trimws(gsub("\\s+"," ", as.character(x))))
  if (grepl("normal",      s)) return("Normal")
  if (grepl("adenocarc|carcinoma|cancer|tumou?r|tumor", s)) return("Tumor")
  return(NA_character_)
}
tissue_vec <- vapply(trow, clean_tissue, character(1))

# Align to merged matrix
samples <- colnames(expr_merged)
tissue_aligned <- tissue_vec[samples]

# ----------------- PLATFORM / STUDY -----------------
platform <- ifelse(samples %in% colnames(norm_rna_z), "RNAseq",
                   ifelse(samples %in% colnames(norm_micro_z), "Microarray", infer_platform(samples)))
study    <- infer_study(samples)

# ----------------- SAVE OUTPUTS -----------------
write.csv(norm_rna,   file.path(getwd(), "norm_rna.csv"),   quote=FALSE)
write.csv(norm_micro, file.path(getwd(), "norm_micro.csv"), quote=FALSE)

write.csv(norm_rna_z,   file.path(getwd(), "norm_rna_z.csv"),   quote=FALSE)
write.csv(norm_micro_z, file.path(getwd(), "norm_micro_z.csv"), quote=FALSE)

write.csv(expr_merged,  file.path(base_dir, "Full_matrix.csv"), quote=FALSE)

cov_df <- data.frame(Sample = samples,
                     Platform = platform,
                     Study = study,
                     Tissue = tissue_aligned,
                     stringsAsFactors = FALSE)
write.csv(cov_df, file.path(base_dir, "CombinedCovariates.csv"), row.names = FALSE)

message("✅ Done: Full_matrix.csv + CombinedCovariates.csv created in base_dir.")
