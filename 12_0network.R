# 12_0network.R
suppressPackageStartupMessages({
  library(minet)
  library(infotheo)
  library(igraph)
  library(dplyr)
  library(readr)
  library(matrixStats)
  library(ggplot2)
  library(stringr)
})

# -------------------- CONFIG --------------------
base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
expr_file <- file.path(base_dir, "All_ComBat_matrix.csv")
cov_file <- file.path(base_dir, "All_Covariates_augmented.csv")
stopifnot(file.exists(expr_file), file.exists(cov_file))

# در 12_0network.R
out_dir <- file.path(base_dir, "Results_Networks_simple")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- LOAD ----------------------
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE) %>% as.matrix()
mode(expr) <- "numeric"
cov <- read.csv(cov_file, stringsAsFactors = FALSE)

stopifnot("Sample" %in% names(cov))
stopifnot(identical(colnames(expr), cov$Sample))

# -------------------- RECODES -------------------
norm_txt <- function(x) tolower(trimws(as.character(x)))
site_norm <- norm_txt(cov$Site)

loc_code <- rep(NA_integer_, nrow(cov))
loc_lab <- rep(NA_character_, nrow(cov))
is_code1 <- str_detect(site_norm, "gastro[_- ]?esoph|gej|fundus|cardia|proximal.*1/3")
is_code2 <- str_detect(site_norm, "body|middle.*1/3")
is_code3 <- str_detect(site_norm, "antrum|distal.*1/3")
loc_code[is_code1] <- 1L
loc_lab[is_code1] <- "Location:Code1"
loc_code[is_code2] <- 2L
loc_lab[is_code2] <- "Location:Code2"
loc_code[is_code3] <- 3L
loc_lab[is_code3] <- "Location:Code3"
cov$LocationCode <- loc_code
cov$Location <- loc_lab

lauren_col <- c("Lauren", "LaurenClassification", "lauren", "laurenclassification")
lauren_col <- lauren_col[lauren_col %in% names(cov)]
if (length(lauren_col)) {
  la <- norm_txt(cov[[lauren_col[1]]])
  la_code <- rep(NA_integer_, nrow(cov))
  la_code[str_detect(la, "diffuse")] <- 1L
  la_code[str_detect(la, "intestinal")] <- 2L
  la_code[str_detect(la, "mix|mixed")] <- 3L
  cov$LaurenCode <- la_code
} else {
  cov$LaurenCode <- NA_integer_
}

age_num <- suppressWarnings(as.numeric(cov$Age))
cov$AgeGroup <- ifelse(!is.na(age_num) & age_num < 40, 1L,
  ifelse(!is.na(age_num) & age_num <= 60, 2L,
    ifelse(!is.na(age_num) & age_num > 60, 3L, NA_integer_)
  )
)

if ("Stage" %in% names(cov)) {
  st <- toupper(gsub("[^IVX]", "", as.character(cov$Stage)))
  st[st %in% c("", NA)] <- NA
  cov$StageClean <- factor(st, levels = c("I", "II", "III", "IV"))
} else {
  cov$StageClean <- NA
}

write.csv(cov, file.path(out_dir, "Covariates_recoded_for_networks.csv"), row.names = FALSE)

# -------------------- GENE SELECTION -------------------
sel <- rowVars(expr)
topn <- min(2000L, nrow(expr))
genes <- names(sort(sel, decreasing = TRUE))[seq_len(topn)]
X <- expr[genes, , drop = FALSE]

# -------------------- NETWORK BUILDER ------------------
build_mi_net <- function(Xsub, out_prefix) {
  if (ncol(Xsub) < 10) {
    message(out_prefix, ": <10 samples, skipped.")
    return(invisible(NULL))
  }
  disc <- infotheo::discretize(t(Xsub), nbins = 3)
  mim <- minet::build.mim(data = disc, estimator = "mi.empirical")
  net <- minet::aracne(mim)

  el <- which(net != 0, arr.ind = TRUE)
  el <- el[el[, 1] < el[, 2], , drop = FALSE]
  if (nrow(el) == 0) {
    message(out_prefix, ": no edges after ARACNE.")
    return(invisible(NULL))
  }

  edges <- data.frame(
    gene1 = rownames(net)[el[, 1]],
    gene2 = colnames(net)[el[, 2]],
    weight = mim[el],
    stringsAsFactors = FALSE
  )
  readr::write_csv(edges, paste0(out_prefix, "_edges.csv"))

  ord <- order(-edges$weight)
  keep <- edges[ord[seq_len(min(200, nrow(edges)))], , drop = FALSE]
  g <- igraph::graph_from_data_frame(keep, directed = FALSE)
  png(paste0(out_prefix, "_graph.png"), width = 1400, height = 1000, res = 140)
  plot(g,
    vertex.size = 4, vertex.label = NA,
    edge.width = scales::rescale(E(g)$weight, to = c(0.4, 3))
  )
  dev.off()
}

# -------------------- GROUPS ---------------------------
groups <- list(
  Sex_Male = which(tolower(cov$Sex) == "male"),
  Sex_Female = which(tolower(cov$Sex) == "female"),
  Stage_I = which(cov$StageClean == "I"),
  Stage_II = which(cov$StageClean == "II"),
  Stage_III = which(cov$StageClean == "III"),
  Stage_IV = which(cov$StageClean == "IV"),
  Location_1 = which(cov$LocationCode == 1L),
  Location_2 = which(cov$LocationCode == 2L),
  Location_3 = which(cov$LocationCode == 3L),
  AgeGrp_1_lt40 = which(cov$AgeGroup == 1L),
  AgeGrp_2_40to60 = which(cov$AgeGroup == 2L),
  AgeGrp_3_gt60 = which(cov$AgeGroup == 3L),
  Lauren_1_diffuse = which(cov$LaurenCode == 1L),
  Lauren_2_intestinal = which(cov$LaurenCode == 2L),
  Lauren_3_mixed = which(cov$LaurenCode == 3L)
)
if ("Tissue" %in% names(cov)) {
  groups$Tissue_Tumor <- which(tolower(cov$Tissue) == "tumor")
  groups$Tissue_Normal <- which(tolower(cov$Tissue) == "normal")
}

# -------------------- BUILD NETWORKS ------------------
for (nm in names(groups)) {
  idx <- groups[[nm]]
  if (length(idx) >= 10) {
    Xsub <- X[, idx, drop = FALSE]
    out_prefix <- file.path(out_dir, paste0("Network_", nm))
    message("Building: ", nm, " (n=", length(idx), ")")
    build_mi_net(Xsub, out_prefix)
  } else {
    message("Skipping: ", nm, " (n=", length(idx), " < 10)")
  }
}
message("Networks built (MI/ARACNE). Outputs in: ", out_dir)
