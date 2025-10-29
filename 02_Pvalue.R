# 02_Pvalue.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# -------------------- CONFIG --------------------
base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
cov_file <- file.path(base_dir, "CombinedCovariates.csv")
if (!file.exists(cov_file)) stop("CombinedCovariates.csv not found at: ", cov_file)

out_dir <- file.path(base_dir, "Results_Covariate_Pvals")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Which covariates to consider (will use only those present)
covariate_order <- c("Tissue", "Sex", "Stage", "LocationCode", "LaurenCode", "Age_Code")

# -------------------- LOAD ----------------------
cov <- read.csv(cov_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))

if (!"Sample" %in% names(cov)) stop("CombinedCovariates.csv must contain a 'Sample' column.")

# normalize Stage to I/II/III/IV (everything else -> NA)
stage_fix <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s <- sub("[ABC]$", "", s) # drop substage letters
  s[s %in% c("1", "I")] <- "I"
  s[s %in% c("2", "II")] <- "II"
  s[s %in% c("3", "III")] <- "III"
  s[s %in% c("4", "IV")] <- "IV"
  s[!s %in% c("I", "II", "III", "IV")] <- NA
  s
}
if ("Stage" %in% names(cov)) cov$Stage <- stage_fix(cov$Stage)
if ("LocationCode" %in% names(cov)) cov$LocationCode <- suppressWarnings(as.integer(cov$LocationCode))
if ("LaurenCode" %in% names(cov)) cov$LaurenCode <- suppressWarnings(as.integer(cov$LaurenCode))
if ("Age_Code" %in% names(cov)) cov$Age_Code <- suppressWarnings(as.integer(cov$Age_Code))

use_cols <- intersect(covariate_order, names(cov))

# -------------------- PRINT INPUT SUMMARY --------------------
cat("==== INPUT SUMMARY ====\n")
cat(sprintf("File: %s\n", cov_file))
cat(sprintf("Rows (samples): %d\n", nrow(cov)))
cat(sprintf("Columns available for testing: %s\n\n", paste(use_cols, collapse = ", ")))

tab_print <- function(x, title, maxk = 20) {
  f <- addNA(factor(x, exclude = NULL), ifany = TRUE)
  tb <- sort(table(f, useNA = "ifany"), decreasing = TRUE)
  cat(sprintf("  %s:\n", title))
  if (length(tb) == 0) {
    cat("    (no data)\n")
    return(invisible(NULL))
  }
  nm <- names(tb)
  nm[nm == "<NA>"] <- "Missing"
  show_n <- min(length(tb), maxk)
  for (i in seq_len(show_n)) cat(sprintf("    - %s: %d\n", nm[i], tb[i]))
  if (length(tb) > show_n) cat(sprintf("    (+%d more)\n", length(tb) - show_n))
  invisible(tb)
}

cat("Level counts (including Missing):\n")
level_counts <- list()
for (v in use_cols) {
  tb <- tab_print(cov[[v]], v)
  if (!is.null(tb)) {
    level_counts[[v]] <- data.frame(Variable = v, Level = names(tb), Count = as.integer(tb), row.names = NULL)
  }
}
cat("==== END SUMMARY ====\n\n")

if (length(level_counts)) {
  write_csv(bind_rows(level_counts), file.path(out_dir, "level_counts.csv"))
}

# -------------------- TEST ENGINE --------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Robust chooser: Chi-square (with simulation if needed) or Fisher (2x2)
chisq_or_fisher <- function(x, y, simulate_p = TRUE) {
  # Keep Missing as an explicit level so nothing is dropped
  x <- addNA(factor(x, exclude = NULL), ifany = TRUE)
  y <- addNA(factor(y, exclude = NULL), ifany = TRUE)

  tab <- table(x, y, useNA = "ifany")
  if (nrow(tab) == 0L || ncol(tab) == 0L) {
    return(list(method = "NA", statistic = NA_real_, df = NA_integer_, p = NA_real_, n = length(na.omit(x))))
  }
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  n <- sum(tab)
  if (n == 0) {
    return(list(method = "NA", statistic = NA_real_, df = NA_integer_, p = NA_real_, n = 0L))
  }

  # 2x2 → Fisher exact
  if (nrow(tab) == 2L && ncol(tab) == 2L) {
    ft <- fisher.test(tab)
    return(list(
      method = "FisherExact(2x2)",
      statistic = unname(ft$estimate %||% NA_real_),
      df = NA_integer_, p = ft$p.value, n = n
    ))
  }

  # Chi-square; use simulation if requested or expected counts are low
  exp <- suppressWarnings((rowSums(tab) %*% t(colSums(tab))) / n)
  do_sim <- isTRUE(simulate_p) || any(exp < 5)

  ct <- try(chisq.test(tab, simulate.p.value = do_sim, B = if (do_sim) 9999 else 0), silent = TRUE)

  if (inherits(ct, "try-error")) {
    return(list(method = "ChiSqFailed", statistic = NA_real_, df = NA_integer_, p = NA_real_, n = n))
  } else {
    df_val <- if (!is.null(ct$parameter)) unname(ct$parameter) else as.integer((nrow(tab) - 1) * (ncol(tab) - 1))
    meth <- if (do_sim) sprintf("ChiSqSim(%dx%d)", nrow(tab), ncol(tab)) else "ChiSq"
    return(list(method = meth, statistic = unname(ct$statistic), df = df_val, p = ct$p.value, n = n))
  }
}

# Optional monotone trend p-value for ordered pairs (Spearman)
is_ordered_code <- function(vname) vname %in% c("Stage", "LocationCode", "LaurenCode", "Age_Code")
stage_to_num <- function(s) {
  s <- toupper(trimws(as.character(s)))
  s[s %in% c("I", "1")] <- "1"
  s[s %in% c("II", "2")] <- "2"
  s[s %in% c("III", "3")] <- "3"
  s[s %in% c("IV", "4")] <- "4"
  suppressWarnings(as.numeric(s))
}
to_numeric_order <- function(vname, x) {
  if (vname == "Stage") {
    return(stage_to_num(x))
  }
  suppressWarnings(as.numeric(x))
}
spearman_trend <- function(v1_name, v1, v2_name, v2) {
  if (!is_ordered_code(v1_name) || !is_ordered_code(v2_name)) {
    return(NA_real_)
  }
  x <- to_numeric_order(v1_name, v1)
  y <- to_numeric_order(v2_name, v2)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) {
    return(NA_real_)
  }
  suppressWarnings(cor.test(x[ok], y[ok], method = "spearman")$p.value)
}

# Core runner for a given data.frame of samples and a set of variables
run_pair_tests <- function(df, varnames, label, drop_constant = TRUE) {
  vars <- intersect(varnames, names(df))
  if (drop_constant) {
    vars <- vars[vapply(vars, function(v) {
      u <- unique(addNA(df[[v]], ifany = TRUE))
      length(u[!is.na(u)]) + if (any(is.na(df[[v]]))) 1L else 0L
    }, integer(1)) >= 2]
  }
  if (length(vars) < 2) {
    message(sprintf("[%s] Not enough variables to test (|vars|=%d).", label, length(vars)))
    out <- data.frame(
      Var1 = character(), Var2 = character(), Method = character(),
      Statistic = numeric(), DF = integer(), P.Value = numeric(),
      N = integer(), Trend_P = numeric(), stringsAsFactors = FALSE
    )
    return(out)
  }

  pairs <- t(combn(vars, 2))
  cat(sprintf("[%-12s] Testing %d pairs across %d samples…\n", label, nrow(pairs), nrow(df)))

  res <- apply(pairs, 1, function(vv) {
    a <- vv[1]
    b <- vv[2]
    test <- chisq_or_fisher(df[[a]], df[[b]])
    data.frame(
      Var1 = a,
      Var2 = b,
      Method = test$method,
      Statistic = test$statistic,
      DF = test$df,
      P.Value = test$p,
      N = test$n,
      Trend_P = spearman_trend(a, df[[a]], b, df[[b]]),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  res <- res %>%
    mutate(adj.P.Val = p.adjust(P.Value, method = "BH")) %>%
    arrange(adj.P.Val, P.Value)

  out_path <- file.path(out_dir, paste0("pairwise_pvalues_", label, ".csv"))
  write_csv(res, out_path)
  message(sprintf("[%-12s] WROTE: %s (pairs=%d)", label, out_path, nrow(res)))
  res
}

# -------------------- RUN TESTS --------------------
if (length(use_cols) < 2) stop("Not enough covariate columns found to test.")

# 1) ALL SAMPLES
res_all <- run_pair_tests(cov, use_cols, label = "All")

# 2) TUMOR / NORMAL strata (if Tissue exists)
if ("Tissue" %in% use_cols) {
  idx_tumor <- which(cov$Tissue == "Tumor")
  idx_normal <- which(cov$Tissue == "Normal")

  if (length(idx_tumor) >= 3) {
    res_tumor <- run_pair_tests(cov[idx_tumor, , drop = FALSE],
      setdiff(use_cols, "Tissue"),
      label = "TumorOnly"
    )
  } else {
    message("[TumorOnly  ] Too few Tumor samples to run (n<3); skipping.")
  }

  if (length(idx_normal) >= 3) {
    res_normal <- run_pair_tests(cov[idx_normal, , drop = FALSE],
      setdiff(use_cols, "Tissue"),
      label = "NormalOnly"
    )
  } else {
    message("[NormalOnly ] Too few Normal samples to run (n<3); skipping.")
  }
}

# 3) Write a small manifest
manifest <- data.frame(
  File = c(
    "level_counts.csv",
    "pairwise_pvalues_All.csv",
    if ("Tissue" %in% use_cols) c("pairwise_pvalues_TumorOnly.csv", "pairwise_pvalues_NormalOnly.csv") else character(0)
  ),
  Description = c(
    "Level counts per variable (including Missing)",
    "Pairwise p-values across all samples",
    if ("Tissue" %in% use_cols) {
      c(
        "Pairwise p-values within Tumor samples (Tissue excluded)",
        "Pairwise p-values within Normal samples (Tissue excluded)"
      )
    } else {
      character(0)
    }
  ),
  stringsAsFactors = FALSE
)
write_csv(manifest, file.path(out_dir, "MANIFEST.csv"))

cat("\nDone. CSVs written in:\n  ", out_dir, "\n")
