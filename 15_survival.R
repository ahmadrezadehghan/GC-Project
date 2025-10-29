# 15_survival.R
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
  library(readxl)
  library(readr)
  library(ggplot2)
})

base_dir <- "E:/1----------GC_29_4_mainDATAFRAME---------1"
cov_file <- file.path(base_dir, "All_Covariates_augmented.csv")
demo_file <- file.path(base_dir, "Demo.xlsx")
stopifnot(file.exists(cov_file), file.exists(demo_file))

out_dir <- file.path(base_dir, "Results_Survival")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cov <- read.csv(cov_file, stringsAsFactors = FALSE)

# ---------- optional GSVA scores (smart loader) ----------
# Prefer the RDS bundle from 06_gsva.R; else try per-collection CSVs; else NULL
gsva_dir <- file.path(base_dir, "Results_GSVA")
gsva_list_rds <- file.path(gsva_dir, "GSVA_scores_list.rds")
gsva_scores <- NULL
if (file.exists(gsva_list_rds)) {
  gsva_scores <- readRDS(gsva_list_rds)
} else {
  csvs <- list.files(gsva_dir, pattern = "^GSVA_scores_.*\\.csv$", full.names = TRUE)
  if (length(csvs)) {
    # load all and keep as a named list
    nm <- sub("^GSVA_scores_(.*)\\.csv$", "\\1", basename(csvs))
    lst <- setNames(vector("list", length(csvs)), nm)
    for (i in seq_along(csvs)) {
      m <- suppressWarnings(read.csv(csvs[i], row.names = 1, check.names = FALSE))
      lst[[i]] <- as.matrix(m)
    }
    names(lst) <- nm
    gsva_scores <- lst
  }
}

# ---------- pull survival from Demo.xlsx ----------
norm_key <- function(x) tolower(gsub("\\s+", " ", trimws(as.character(x))))
demo <- readxl::read_excel(demo_file) %>% as.data.frame()
features_raw <- demo[[1]]
features_key <- norm_key(features_raw)
demo_vals <- function(idx) {
  r <- demo[idx, -1, drop = TRUE]
  names(r) <- names(demo)[-1]
  r
}
find_row <- function(cands) {
  hit <- which(features_key %in% norm_key(cands))
  ifelse(length(hit) >= 1, hit[1], NA_integer_)
}

time_vec <- {
  idx <- find_row(c("OS time", "Overall Survival Time", "Survival Time", "Days to event", "OS_days", "Time_to_event", "OS_days", "OS.months"))
  if (!is.na(idx)) demo_vals(idx) else setNames(rep(NA, ncol(demo) - 1), names(demo)[-1])
}
event_vec <- {
  idx <- find_row(c("OS event", "Vital Status", "Event", "Status", "OS_status", "Death_event", "OS.event"))
  if (!is.na(idx)) demo_vals(idx) else setNames(rep(NA, ncol(demo) - 1), names(demo)[-1])
}

cov$OS_time <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", time_vec[cov$Sample])))
ev <- tolower(as.character(event_vec[cov$Sample]))
cov$OS_event <- ifelse(ev %in% c("1", "dead", "deceased", "death", "event", "true", "yes"), 1,
  ifelse(ev %in% c("0", "alive", "censored", "no", "false"), 0, NA)
)

# ---------- KM by demographics ----------
km_dir <- file.path(out_dir, "KM")
dir.create(km_dir, showWarnings = FALSE)

do_km <- function(group, name) {
  keep <- !is.na(cov$OS_time) & !is.na(cov$OS_event) & !is.na(group)
  d <- data.frame(OS_time = cov$OS_time[keep], OS_event = cov$OS_event[keep])
  d$grp <- factor(group[keep])
  if (nlevels(d$grp) < 2) {
    return(invisible(NULL))
  }
  fit <- survfit(Surv(OS_time, OS_event) ~ grp, data = d)
  p <- ggsurvplot(fit,
    risk.table = TRUE, pval = TRUE, conf.int = FALSE, ggtheme = theme_bw(),
    title = paste("KM:", name)
  )
  ggsave(file.path(km_dir, paste0("KM_", gsub("[^A-Za-z0-9]+", "_", name), ".png")),
    width = 9, height = 7, dpi = 150, plot = p$plot
  )
  invisible(TRUE)
}
do_km(cov$Sex, "Sex (Male vs Female)")
do_km(cov$Stage, "Stage (I/II/III/IV)")

# ---------- Adjusted Cox across all demographics ----------
cox_df <- cov %>%
  dplyr::select(OS_time, OS_event, Sex, Stage, Age, Study, Platform, Tissue) %>%
  na.omit()
if (nrow(cox_df) >= 50) {
  cox_df$Sex <- factor(cox_df$Sex)
  cox_df$Stage <- factor(cox_df$Stage, levels = c("I", "II", "III", "IV"))
  cox_df$Study <- factor(cox_df$Study)
  cox_df$Platform <- factor(cox_df$Platform)
  cox_df$Tissue <- factor(ifelse(is.na(cox_df$Tissue), "Unknown", cox_df$Tissue))
  fit <- coxph(Surv(OS_time, OS_event) ~ Sex + Stage + Age + Study + Platform + Tissue, data = cox_df)
  summ <- summary(fit)
  write.csv(as.data.frame(summ$coefficients), file.path(out_dir, "Cox_adjusted_coefficients.csv"))
  ctab <- as.data.frame(summ$conf.int)
  ctab$var <- rownames(ctab)
  p <- ggplot(ctab, aes(x = var, y = `exp(coef)`, ymin = `lower .95`, ymax = `upper .95`)) +
    geom_pointrange() +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(title = "Adjusted Cox model", y = "Hazard Ratio", x = NULL)
  ggsave(file.path(out_dir, "Cox_adjusted_forest.png"), width = 9, height = 7, dpi = 150, plot = p)
}

# ---------- Optional: pathway-based survival (median split) ----------
if (!is.null(gsva_scores) && length(gsva_scores) > 0) {
  # نمونه: از HALLMARK ها اگر موجود بود
  if ("H" %in% names(gsva_scores)) {
    sc <- gsva_scores[["H"]]
  } else {
    sc <- gsva_scores[[1]]
  }
  common <- intersect(colnames(sc), cov$Sample)
  if (length(common) >= 30) {
    s <- cov[match(common, cov$Sample), ]
    pick <- intersect(rownames(sc), c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ANGIOGENESIS"))
    if (length(pick) > 0) {
      score <- colMeans(sc[pick, common, drop = FALSE])
      d <- data.frame(
        OS_time = s$OS_time, OS_event = s$OS_event,
        grp = factor(ifelse(score >= median(score, na.rm = TRUE), "High", "Low"))
      )
      d <- d[complete.cases(d), ]
      if (nrow(d) >= 30) {
        fit <- survfit(Surv(OS_time, OS_event) ~ grp, data = d)
        p <- ggsurvplot(fit,
          risk.table = TRUE, pval = TRUE, ggtheme = theme_bw(),
          title = paste0("KM by pathway score: ", paste(pick, collapse = "+"))
        )
        ggsave(file.path(out_dir, "KM_pathway_score.png"), width = 9, height = 7, dpi = 150, plot = p$plot)
      }
    }
  }
}
message("Survival analyses completed. Outputs in: ", out_dir)
