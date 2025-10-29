# 00_Covfilecreator.R
# Output EXACT columns:
#   Sample, Platform, Tissue, Sex, Stage, LocationCode, LaurenCode, Age_Code
# No NA in file (missing values written as empty strings)

suppressPackageStartupMessages({ library(readxl) })

# -------------------- CONFIG --------------------
base_dir        <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
demo_file       <- file.path(base_dir, "Demo.xlsx")
rna_expr_file   <- file.path(base_dir, "RNA_seq.xlsx")
micro_expr_file <- file.path(base_dir, "Micro.xlsx")
expr_file       <- file.path(base_dir, "Full_matrix.csv")  # optional for sample order

if (!file.exists(demo_file)) stop("Demo.xlsx not found at: ", demo_file)

# -------------------- HELPERS --------------------
norm_key <- function(x) tolower(gsub("\\s+"," ", trimws(as.character(x))))
find_row <- function(features_vec, candidates){
  key <- norm_key(features_vec); cand <- norm_key(candidates)
  hit <- which(key %in% cand); if (length(hit) >= 1) hit[1] else NA_integer_
}
clean_tissue <- function(x){
  if (is.na(x)) return(NA_character_)
  s <- tolower(trimws(gsub("\\s+"," ", as.character(x))))
  if (!nzchar(s)) return(NA_character_)
  if (grepl("(^|\\b)normal\\b|adjacent normal|non[- ]tumou?r", s)) return("Normal")
  if (grepl("adenocarc|carcinoma|cancer|tumou?r|tumor", s))       return("Tumor")
  NA_character_
}
clean_sex <- function(x){
  if (is.na(x)) return(NA_character_)
  s <- tolower(trimws(as.character(x))); if (!nzchar(s)) return(NA_character_)
  if (s %in% c("m","male")) "Male" else if (s %in% c("f","female")) "Female" else NA_character_
}
clean_stage_overall <- function(x){
  if (is.na(x)) return(NA_character_); s <- toupper(trimws(as.character(x)))
  s <- sub("[ABC]$", "", s)                      # drop substage letter (IA->I etc.)
  if (s %in% c("I","1"))   "I" else
    if (s %in% c("II","2"))  "II" else
      if (s %in% c("III","3")) "III" else
        if (s %in% c("IV","4"))  "IV" else NA_character_
}
clean_T <- function(x){
  if (is.na(x)) return(NA_character_)
  s <- tolower(gsub("[^a-z0-9]", "", as.character(x)))
  if (grepl("tis", s)) "Tis" else
    if (grepl("t4b", s)) "T4b" else
      if (grepl("t4a|t4", s)) "T4a" else
        if (grepl("t3",  s)) "T3"  else
          if (grepl("t2",  s)) "T2"  else
            if (grepl("t1",  s)) "T1"  else
              if (grepl("t0",  s)) "T0"  else NA_character_
}
clean_N <- function(x){
  if (is.na(x)) return(NA_character_)
  s <- tolower(gsub("[^a-z0-9]", "", as.character(x)))
  if (grepl("n3b", s)) "N3b" else
    if (grepl("n3a|n3", s)) "N3a" else
      if (grepl("n2", s)) "N2" else
        if (grepl("n1", s)) "N1" else
          if (grepl("n0", s)) "N0" else NA_character_
}
clean_M <- function(x){
  if (is.na(x)) return(NA_character_)
  s <- tolower(gsub("[^a-z0-9]", "", as.character(x)))
  if (grepl("m1", s)) "M1" else
    if (grepl("m0", s)) "M0" else NA_character_
}
clean_location_code <- function(x){
  if (is.na(x)) return(NA_integer_)
  if (is.numeric(x)) { v <- as.integer(round(x)); if (v %in% 1:3) return(v) }
  s <- tolower(trimws(as.character(x))); if (!nzchar(s)) return(NA_integer_)
  if (grepl("gej|cardia|fundus|proximal|upper|prox(imal)? ?1/?3", s)) 1L else
    if (grepl("\\b(body|middle)\\b|mid( ?1/?3)?", s))                 2L else
      if (grepl("antrum|distal|lower|dist(al)? ?1/?3", s))            3L else NA_integer_
}
clean_lauren_code <- function(x){
  if (is.na(x)) return(NA_integer_)
  if (is.numeric(x)) { v <- as.integer(round(x)); if (v %in% 1:3) return(v) }
  s <- tolower(trimws(as.character(x))); if (!nzchar(s)) return(NA_integer_)
  if (grepl("diffuse", s)) 1L else if (grepl("intestinal", s)) 2L else
    if (grepl("mixed|indeterminate", s)) 3L else NA_integer_
}
clean_age_code <- function(x){
  if (is.na(x)) return(NA_integer_)
  if (is.numeric(x)) { a <- as.numeric(x); if (is.na(a)) return(NA_integer_);
  if (a < 40) 1L else if (a <= 60) 2L else 3L }
  s <- tolower(trimws(as.character(x))); if (!nzchar(s)) return(NA_integer_)
  if (s %in% c("1","u40","<40","under40","younger")) 1L else
    if (s %in% c("2","40to60","40-60","40–60","40<=age<=60")) 2L else
      if (s %in% c("3",">60","gt60","older","60+")) 3L else {
        suppressWarnings({ v <- as.numeric(s);
        if (!is.na(v)) return(if (v < 40) 1L else if (v <= 60) 2L else 3L) })
        NA_integer_
      }
}
merge_chr <- function(canon, raw){
  out <- canon; need <- is.na(out) | out == ""
  if (any(need)) {
    r <- as.character(raw); r[!nzchar(trimws(r))] <- NA_character_
    out[need] <- r[need]
  }
  out
}
merge_code <- function(canon, raw){
  out <- canon; need <- is.na(out)
  if (any(need)) {
    r <- suppressWarnings(as.integer(round(as.numeric(as.character(raw)))))
    out[need] <- r[need]
  }
  out
}
truthy <- function(x){
  if (is.null(x)) return(NA)
  y <- tolower(trimws(as.character(x)))
  ifelse(y %in% c("1","true","t","yes","y"), TRUE,
         ifelse(y %in% c("0","false","f","no","n",""), FALSE, NA))
}

# -------------------- SAMPLES --------------------
if (file.exists(expr_file)) {
  hdr <- read.csv(expr_file, nrows = 1, check.names = FALSE)
  cn <- colnames(hdr); if (length(cn) < 2) stop("Full_matrix.csv header appears invalid.")
  samples <- cn[-1]
} else {
  rna_samples <- micro_samples <- character(0)
  if (file.exists(rna_expr_file))   rna_samples   <- colnames(read_excel(rna_expr_file,   n_max = 0))[-1]
  if (file.exists(micro_expr_file)) micro_samples <- colnames(read_excel(micro_expr_file, n_max = 0))[-1]
  samples <- unique(c(rna_samples, setdiff(micro_samples, rna_samples)))
  if (length(samples) == 0) stop("No samples inferred from expression files.")
}

# -------------------- READ DEMO --------------------
demo <- as.data.frame(read_excel(demo_file)); stopifnot(ncol(demo) >= 2)
features <- demo[[1]]
grab_row <- function(primary, alts = character(0)){
  nm <- colnames(demo)[-1]
  idx <- find_row(features, c(primary, alts))
  if (is.na(idx)) { v <- rep(NA_character_, length(nm)); names(v) <- nm; return(v) }
  v <- as.vector(t(demo[idx, -1, drop = FALSE])); names(v) <- nm; v
}

# Main rows (raw)
Tissue_raw   <- grab_row("tissue",      c("tissue type","sample tissue"))
Sex_raw      <- grab_row("sex",         c("gender"))
Stage_raw    <- grab_row("stage",       c("tumor stage","tumour stage","stage (ajcc)","ajcc stage","pstage","cstage"))
T_raw        <- grab_row("t_stage",     c("t","pt","ct","t stage","t-stage"))
N_raw        <- grab_row("n_stage",     c("n","pn","cn","n stage","n-stage"))
M_raw        <- grab_row("m_stage",     c("m","pm","cm","m stage","m-stage"))
Location_raw <- grab_row("locationcode",c("location code","location","site","tumor location","anatomic site"))
Lauren_raw   <- grab_row("laurencode",  c("lauren","lauren classification"))
AgeCode_raw  <- grab_row("age_code",    c("agecode","age code"))
Age_raw      <- grab_row("age",         c("age (years)","patient age"))

# Optional one-hot flags (only to fill gaps)
Sex_Male_raw    <- grab_row("sex_male", c())
Sex_Female_raw  <- grab_row("sex_female", c())
Stage_I_raw     <- grab_row("stage_i", c("stage i"))
Stage_II_raw    <- grab_row("stage_ii", c("stage ii"))
Stage_III_raw   <- grab_row("stage_iii", c("stage iii"))
Stage_IV_raw    <- grab_row("stage_iv", c("stage iv"))
Tissue_Tumor_raw  <- grab_row("tissue_tumor", c())
Tissue_Normal_raw <- grab_row("tissue_normal", c())

# Align to sample order
align <- function(v) unname(v[samples])
Tissue_raw   <- align(Tissue_raw);   Sex_raw  <- align(Sex_raw);   Stage_raw <- align(Stage_raw)
T_raw        <- align(T_raw);        N_raw    <- align(N_raw);      M_raw     <- align(M_raw)
Location_raw <- align(Location_raw); Lauren_raw <- align(Lauren_raw)
AgeCode_raw  <- align(AgeCode_raw);  Age_raw  <- align(Age_raw)
Sex_Male_raw <- align(Sex_Male_raw); Sex_Female_raw <- align(Sex_Female_raw)
Stage_I_raw  <- align(Stage_I_raw);  Stage_II_raw <- align(Stage_II_raw)
Stage_III_raw<- align(Stage_III_raw);Stage_IV_raw <- align(Stage_IV_raw)
Tissue_Tumor_raw <- align(Tissue_Tumor_raw); Tissue_Normal_raw <- align(Tissue_Normal_raw)

# -------------------- CLEAN / MERGE --------------------
Tissue_c <- vapply(Tissue_raw, clean_tissue, character(1))
Sex_c    <- vapply(Sex_raw,    clean_sex,    character(1))
Stage_c  <- vapply(Stage_raw,  clean_stage_overall, character(1))

# Fill Sex/Tissue from flags if missing
fill_from_flags <- function(curr, yes_vec, no_vec = NULL, yes_val = NULL, no_val = NULL){
  out <- curr
  for (i in seq_along(out)) {
    if (is.na(out[i])) {
      y <- truthy(yes_vec[i]); n <- if (!is.null(no_vec)) truthy(no_vec[i]) else NA
      if (isTRUE(y) && !isTRUE(n))       out[i] <- yes_val
      else if (isTRUE(n) && !isTRUE(y))  out[i] <- no_val
    }
  }
  out
}
Sex_c    <- fill_from_flags(Sex_c,    Sex_Male_raw,  Sex_Female_raw,  "Male",  "Female")
Sex_c    <- fill_from_flags(Sex_c,    Sex_Female_raw, Sex_Male_raw,   "Female","Male")
Tissue_c <- fill_from_flags(Tissue_c, Tissue_Tumor_raw, Tissue_Normal_raw, "Tumor", "Normal")
Tissue_c <- fill_from_flags(Tissue_c, Tissue_Normal_raw, Tissue_Tumor_raw, "Normal","Tumor")

# Fill Stage from Stage_I..IV flags if missing
which_flag <- function(i_vec, ii_vec, iii_vec, iv_vec){
  out <- rep(NA_character_, length(i_vec))
  for (k in seq_along(out)) {
    vals <- c(I=truthy(i_vec[k]), II=truthy(ii_vec[k]), III=truthy(iii_vec[k]), IV=truthy(iv_vec[k]))
    if (sum(isTRUE(vals)) == 1) out[k] <- names(vals)[which(isTRUE(vals))]
  }
  out
}
Stage_flag <- which_flag(Stage_I_raw, Stage_II_raw, Stage_III_raw, Stage_IV_raw)
Stage_c[is.na(Stage_c)] <- Stage_flag[is.na(Stage_c)]

# Safe inference from TNM: if M1 -> Stage IV
T_c <- vapply(T_raw, clean_T, character(1))
N_c <- vapply(N_raw, clean_N, character(1))
M_c <- vapply(M_raw, clean_M, character(1))
needs_stage <- is.na(Stage_c) & !is.na(M_c) & M_c == "M1"
Stage_c[needs_stage] <- "IV"

# Normals: no stage (leave blank in file)
is_normal <- !is.na(Tissue_c) & tolower(Tissue_c) == "normal"
Stage_c[is_normal] <- NA_character_

# Location/Lauren/Age_Code (cleaned; numbers kept numeric)
LocationCode_c <- vapply(Location_raw, clean_location_code, integer(1))
LaurenCode_c   <- vapply(Lauren_raw,   clean_lauren_code,   integer(1))
AgeCode_from_row <- vapply(AgeCode_raw, clean_age_code, integer(1))
AgeCode_from_age <- vapply(Age_raw,     clean_age_code, integer(1))
Age_Code_c <- ifelse(!is.na(AgeCode_from_row), AgeCode_from_row, AgeCode_from_age)

# Final merged fields (canonical preferred, fall back to raw)
Tissue       <- merge_chr(Tissue_c, Tissue_raw)
Sex          <- merge_chr(Sex_c,    Sex_raw)
Stage        <- merge_chr(Stage_c,  Stage_raw)
LocationCode <- merge_code(LocationCode_c, Location_raw)
LaurenCode   <- merge_code(LaurenCode_c,   Lauren_raw)
Age_Code     <- merge_code(Age_Code_c,     ifelse(!is.na(AgeCode_raw), AgeCode_raw, Age_raw))

# -------------------- PLATFORM --------------------
rna_hdr2   <- if (file.exists(rna_expr_file))   read_excel(rna_expr_file,   n_max = 0) else NULL
micro_hdr2 <- if (file.exists(micro_expr_file)) read_excel(micro_expr_file, n_max = 0) else NULL
rna_samples   <- if (!is.null(rna_hdr2)   && ncol(rna_hdr2)   >= 2) colnames(rna_hdr2)[-1]   else character(0)
micro_samples <- if (!is.null(micro_hdr2) && ncol(micro_hdr2) >= 2) colnames(micro_hdr2)[-1] else character(0)
Platform <- ifelse(samples %in% rna_samples, "RNAseq",
                   ifelse(samples %in% micro_samples, "Microarray", NA_character_))

# -------------------- OUTPUT (exact columns) --------------------
cov_df <- data.frame(
  Sample = samples,
  Platform = Platform,
  Tissue = Tissue,
  Sex = Sex,
  Stage = Stage,
  LocationCode = LocationCode,
  LaurenCode = LaurenCode,
  Age_Code = Age_Code,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Write with EMPTY cells for missing (no "NA" strings)
out_path <- file.path(base_dir, "CombinedCovariates.csv")
write.csv(cov_df, out_path, row.names = FALSE, na = "")
cat(sprintf("✅ Wrote %d rows × %d cols to %s (missing values written as empty)\n",
            nrow(cov_df), ncol(cov_df), out_path))
