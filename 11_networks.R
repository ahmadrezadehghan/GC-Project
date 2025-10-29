# 11_networks.R
suppressPackageStartupMessages({
  library(minet)
  library(infotheo)
  library(igraph)
  library(dplyr)
  library(readr)
  library(matrixStats)
  library(ggplot2)
  library(survival)
  library(limma)
})

set.seed(42)

# -------------------- CONFIG --------------------
base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
expr_file <- file.path(base_dir, "All_ComBat_matrix.csv")
cov_file <- file.path(base_dir, "All_Covariates_augmented.csv")
stopifnot(file.exists(expr_file), file.exists(cov_file))

out_dir <- file.path(base_dir, "Results_Networks")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# MI settings
topN_genes <- 2000L # top variable genes
min_samples <- 10L # minimum samples per subgroup to build a network
mi_bins <- 3L # discretization bins
preview_edges <- 200L # edges to plot in quick preview graph

# -------------------- LOAD ----------------------
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE) |> as.matrix()
mode(expr) <- "numeric"
cov <- read.csv(cov_file, stringsAsFactors = FALSE)

# align expr <-> cov and de-duplicate sample IDs if needed
stopifnot(identical(colnames(expr), cov$Sample))
if (anyDuplicated(cov$Sample)) {
  new_ids <- make.unique(cov$Sample)
  message("Made unique sample IDs in cov (", sum(duplicated(cov$Sample)), " duplicates).")
  cov$Sample <- new_ids
  colnames(expr) <- new_ids
}
stopifnot(identical(colnames(expr), cov$Sample))

# -------------------- DEMOGRAPHIC CODINGS --------------------
# Stage: keep canonical I/II/III/IV, else NA/Unknown
cov$Stage <- toupper(gsub("stage\\s*", "", cov$Stage))
cov$Stage[!cov$Stage %in% c("I", "II", "III", "IV")] <- NA

# --- Location -> LocationCode ---
pick_first_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) df[[hit[1]]] else rep(NA_character_, nrow(df))
}
loc_raw <- pick_first_col(cov, c("Location", "location", "PrimarySite", "Site", "site", "Primary_Site"))

norm_txt <- function(x) tolower(trimws(as.character(x)))

map_location_to_code <- function(x) {
  x <- norm(x)
  out <- rep(NA_integer_, length(x))
  out[grepl("gastro.?esophageal|gastroesophageal|ge junction|junction|cardia", x)] <- 1
  out[grepl("^cardia$", x)] <- 1
  out[grepl("fundus", x)] <- 1
  out[grepl("proximal\\s*1/3|proximal third|upper third", x)] <- 1

  out[grepl("body", x)] <- 2
  out[grepl("middle\\s*1/3|middle third|mid third", x)] <- 2

  out[grepl("antrum", x)] <- 3
  out[grepl("distal\\s*1/3|distal third|lower third", x)] <- 3

  out[is.na(out) & grepl("gej", x)] <- 1
  out
}
cov$LocationText <- ifelse(is.na(loc_raw) | loc_raw == "", NA, loc_raw)
cov$LocationCode <- map_location_to_code(cov$LocationText)

# --- Lauren classification -> LaurenCode ---
laur_raw <- pick_first_col(cov, c("Lauren", "LaurenClass", "Lauren_classification", "laurenclassification", "LaurenClassification"))
map_lauren_to_code <- function(x) {
  x <- norm(x)
  out <- rep(NA_integer_, length(x))
  out[grepl("diffuse", x)] <- 1
  out[grepl("intestinal", x)] <- 2
  out[grepl("mixed|mix", x)] <- 3
  out
}
cov$LaurenText <- ifelse(is.na(laur_raw) | laur_raw == "", NA, laur_raw)
cov$LaurenCode <- map_lauren_to_code(cov$LaurenText)

# --- Age -> numeric + AgeCode ---
age_num <- suppressWarnings(as.numeric(cov$Age))
cov$Age <- age_num
cov$AgeCode <- ifelse(!is.na(age_num) & age_num < 40, 1L,
  ifelse(!is.na(age_num) & age_num <= 60, 2L,
    ifelse(!is.na(age_num) & age_num > 60, 3L, NA)
  )
)

# Sex canonical
sx <- tolower(cov$Sex)
sx[sx %in% c("m", "male", "man", "boy")] <- "male"
sx[sx %in% c("f", "female", "woman", "girl")] <- "female"
cov$Sex <- ifelse(is.na(sx), NA, ifelse(sx == "male", "Male", "Female"))

# Tissue canonical
cov$Tissue <- ifelse(is.na(cov$Tissue) | cov$Tissue == "", "Unknown", cov$Tissue)

# -------------------- SURVIVAL: detect & screen demographics --------------------
find_surv_cols <- function(cov) {
  time_cands <- c(
    "OS_time", "OS.time", "OS_days", "OS.months", "Overall_Survival_Time",
    "OverallSurvivalMonths", "OverallSurvivalDays", "SurvivalTime", "time", "Time"
  )
  status_cands <- c(
    "OS_event", "OS.event", "OS_status", "OS.status", "vital_status",
    "Status", "Event", "event", "status"
  )
  pick <- function(cands) {
    hit <- intersect(cands, names(cov))
    if (length(hit)) hit[1] else NA_character_
  }
  list(time = pick(time_cands), status = pick(status_cands))
}

make_surv <- function(cov, time_col, status_col) {
  tt <- suppressWarnings(as.numeric(cov[[time_col]]))
  ss_raw <- cov[[status_col]]
  if (is.numeric(ss_raw)) {
    ss <- as.integer(ss_raw > 0)
  } else {
    sstr <- tolower(trimws(as.character(ss_raw)))
    ss <- ifelse(sstr %in% c("1", "true", "event", "dead", "deceased", "death"), 1L, 0L)
  }
  ok <- is.finite(tt) & !is.na(ss)
  list(S = Surv(tt[ok], ss[ok]), ok = ok)
}

sc <- find_surv_cols(cov)
# بعد از find_surv_cols:
if (is.na(sc$time)) sc$time <- "MySurvivalTimeColumn"
if (is.na(sc$status)) sc$status <- "MyEventColumn"

do_surv <- !(is.na(sc$time) || is.na(sc$status))

sig_covars <- character(0)
if (do_surv) {
  sv <- make_surv(cov, sc$time, sc$status)
  demo_vars <- c("Sex", "Stage", "LocationCode", "LaurenCode", "Age")
  surv_rows <- lapply(demo_vars, function(v) {
    if (!v %in% names(cov)) {
      return(NULL)
    }
    x <- cov[[v]][sv$ok]
    if (all(is.na(x))) {
      return(NULL)
    }
    d <- data.frame(x = x)
    if (is.numeric(x)) {
      fit <- try(coxph(sv$S ~ x, data = d), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(NULL)
      }
      s <- summary(fit)
      data.frame(var = v, type = "numeric", hr = exp(coef(fit)), wald_p = s$wald["pvalue"], stringsAsFactors = FALSE)
    } else {
      d$x <- factor(d$x)
      if (nlevels(d$x) < 2) {
        return(NULL)
      }
      fit <- try(coxph(sv$S ~ x, data = d), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(NULL)
      }
      s <- summary(fit)
      data.frame(var = v, type = "factor", hr = NA_real_, wald_p = s$wald["pvalue"], stringsAsFactors = FALSE)
    }
  })
  surv_tab <- do.call(rbind, surv_rows)
  if (!is.null(surv_tab)) {
    surv_tab$padj <- p.adjust(surv_tab$wald_p, method = "BH")
    readr::write_csv(surv_tab, file.path(out_dir, "Survival_univariate_covariates.csv"))
    sig_covars <- surv_tab$var[surv_tab$padj < 0.05]
  }
  message(
    "Survival covariate screen (FDR<0.05): ",
    ifelse(length(sig_covars) == 0, "none", paste(sig_covars, collapse = ", "))
  )
} else {
  warning("Could not detect survival time/status columns. Skipping survival-based adjustment.")
}

# helpers for residualization before MI
build_cov_mm <- function(cov_sub, vars) {
  if (length(vars) == 0) {
    return(NULL)
  }
  df <- cov_sub[, vars, drop = FALSE]
  for (nm in names(df)) if (is.character(df[[nm]])) df[[nm]] <- factor(df[[nm]])
  mm <- model.matrix(~ 0 + ., data = df)
  if (ncol(mm) == 0) {
    return(NULL)
  }
  mm
}

infer_group_var <- function(group_name) {
  if (grepl("^Sex_", group_name)) {
    return("Sex")
  }
  if (grepl("^Stage_", group_name)) {
    return("Stage")
  }
  if (grepl("^Location_", group_name)) {
    return("LocationCode")
  }
  if (grepl("^Lauren_", group_name)) {
    return("LaurenCode")
  }
  if (grepl("^Age_", group_name)) {
    return("Age")
  }
  if (grepl("^Tissue_", group_name)) {
    return("Tissue")
  }
  NA_character_
}

# -------------------- GENE SELECTION --------------------
sel <- rowVars(expr)
sel[is.na(sel)] <- 0
topN <- min(topN_genes, sum(sel > 0))
genes <- names(sort(sel, decreasing = TRUE))[seq_len(topN)]
X <- expr[genes, , drop = FALSE]

# -------------------- GROUP DEFINITIONS --------------------
groups <- list(
  Sex_Male = which(cov$Sex == "Male"),
  Sex_Female = which(cov$Sex == "Female"),
  Stage_I = which(cov$Stage == "I"),
  Stage_II = which(cov$Stage == "II"),
  Stage_III = which(cov$Stage == "III"),
  Stage_IV = which(cov$Stage == "IV"),
  Location_1 = which(cov$LocationCode == 1),
  Location_2 = which(cov$LocationCode == 2),
  Location_3 = which(cov$LocationCode == 3),
  Lauren_Diffuse = which(cov$LaurenCode == 1),
  Lauren_Intestinal = which(cov$LaurenCode == 2),
  Lauren_Mixed = which(cov$LaurenCode == 3),
  Age_Code1_u40 = which(cov$AgeCode == 1),
  Age_Code2_40to60 = which(cov$AgeCode == 2),
  Age_Code3_gt60 = which(cov$AgeCode == 3),
  Tissue_Tumor = which(cov$Tissue == "Tumor"),
  Tissue_Normal = which(cov$Tissue == "Normal")
)

# save group sizes
sizes <- vapply(groups, length, integer(1))
write.csv(data.frame(Group = names(groups), N = sizes),
  file.path(out_dir, "Group_sizes.csv"),
  row.names = FALSE
)

# -------------------- NETWORK BUILDER --------------------
build_mi_net <- function(Xsub, out_prefix, bins = mi_bins, preview_n = preview_edges) {
  if (ncol(Xsub) < min_samples) {
    return(invisible(NULL))
  }
  disc <- infotheo::discretize(data.frame(t(Xsub)), nbins = bins)
  mim <- minet::build.mim(data = disc, estimator = "mi.empirical")
  adj <- minet::aracne(mim)

  idx <- which(adj != 0, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  if (nrow(idx) == 0) {
    message("No edges for: ", out_prefix)
    return(invisible(NULL))
  }

  edges <- data.frame(
    gene1 = rownames(adj)[idx[, 1]],
    gene2 = colnames(adj)[idx[, 2]],
    weight = mapply(function(i, j) mim[i, j], idx[, 1], idx[, 2]),
    stringsAsFactors = FALSE
  )

  deg <- table(c(edges$gene1, edges$gene2))
  deg_df <- data.frame(gene = names(deg), degree = as.integer(deg))
  deg_df <- deg_df[order(-deg_df$degree), ]

  write.csv(edges, paste0(out_prefix, "_edges.csv"), row.names = FALSE)
  write.csv(deg_df, paste0(out_prefix, "_node_degrees.csv"), row.names = FALSE)

  topE <- edges[order(-edges$weight), , drop = FALSE]
  topE <- topE[seq_len(min(preview_n, nrow(topE))), , drop = FALSE]
  g <- igraph::graph_from_data_frame(topE, directed = FALSE)

  png(paste0(out_prefix, "_graph.png"), width = 1400, height = 1000, res = 140)
  plot(g,
    vertex.size = 4, vertex.label = NA,
    edge.width = scales::rescale(E(g)$weight, to = c(0.2, 3))
  )
  dev.off()

  igraph::write_graph(g, paste0(out_prefix, "_preview.graphml"), format = "graphml")
  invisible(list(edges = edges, degrees = deg_df))
}

message("Building networks (MI/ARACNE) with survival-informed adjustment when available…")

for (nm in names(groups)) {
  idx <- groups[[nm]]
  if (length(idx) < min_samples) next
  Xsub <- X[, idx, drop = FALSE]

  grp_var <- infer_group_var(nm)
  adj_vars <- setdiff(sig_covars, grp_var)

  if (length(adj_vars) > 0) {
    mm <- build_cov_mm(cov[idx, , drop = FALSE], adj_vars)
    if (!is.null(mm)) {
      Xsub <- limma::removeBatchEffect(Xsub, covariates = mm)
    }
  }

  out_prefix <- file.path(out_dir, paste0("Network_", nm))
  build_mi_net(Xsub, out_prefix)
}

message("Networks built (MI/ARACNE). Outputs in: ", out_dir)
