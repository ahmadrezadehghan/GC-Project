# 16_survival_advanced.R
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(survival)
    library(survminer)
    library(glmnet)
    library(timeROC)
    library(pROC)
})

base_dir <- Sys.getenv("GC_BASE_DIR", unset = "E:/1----------GC_29_4_mainDATAFRAME---------1")
expr_file <- file.path(base_dir, "All_ComBat_matrix.csv")
cov_file <- file.path(base_dir, "All_Covariates_augmented.csv")
wgcna_dir <- file.path(base_dir, "Results_WGCNA")
cent_file <- file.path(base_dir, "Results_Integration", "Tumor_ARACNE_centrality.csv") # اختیاری
out_dir <- file.path(base_dir, "Results_Survival_Adv")
dir.create(out_dir, TRUE, TRUE)

stopifnot(file.exists(expr_file), file.exists(cov_file))

cov <- read.csv(cov_file, stringsAsFactors = FALSE)
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE) %>% as.matrix()
stopifnot(identical(colnames(expr), cov$Sample))

# ---- آماده‌سازی بقا ----
time_cols <- c("OS_time", "OS.days", "OS_days", "SurvivalTime", "Time")
event_cols <- c("OS_event", "Status", "Event", "status")
tcol <- intersect(time_cols, names(cov))
ecol <- intersect(event_cols, names(cov))
stopifnot(length(tcol) >= 1, length(ecol) >= 1)
tt <- suppressWarnings(as.numeric(cov[[tcol[1]]]))
ev <- cov[[ecol[1]]]
if (!is.numeric(ev)) ev <- ifelse(tolower(ev) %in% c("1", "dead", "deceased", "event", "true", "yes"), 1, 0)
ok <- is.finite(tt) & !is.na(ev)
cov <- cov[ok, ]
expr <- expr[, ok, drop = FALSE]
tt <- tt[ok]
ev <- ev[ok]

# ---- ویژگی‌ها: Module Eigengenes + (اختیاری) هاب‌ژن‌ها ----
ME_file <- file.path(wgcna_dir, "WGCNA_ModuleEigengenes.csv")
Xf <- NULL
if (file.exists(ME_file)) {
    MEs <- read.csv(ME_file, row.names = 1, check.names = FALSE)
    MEs <- MEs[rownames(MEs) %in% cov$Sample, , drop = FALSE]
    MEs <- MEs[cov$Sample, , drop = FALSE]
    Xf <- as.matrix(MEs)
}
if (file.exists(cent_file)) {
    cen <- read.csv(cent_file)
    hubs <- cen %>%
        filter(degree >= quantile(degree, 0.95, na.rm = TRUE)) %>%
        pull(Gene)
    hubs <- intersect(hubs, rownames(expr))
    if (length(hubs) >= 5) {
        Xh <- t(scale(t(expr[hubs, cov$Sample, drop = FALSE])))
        colnames(Xh) <- paste0("HUB_", colnames(Xh)) # ← اشتباه؛ باید ردیف‌ها ژن هستند
        # اصلاح: ویژگی‌های هاب = میانگین یا PC1
        Xh <- t(expr[hubs, cov$Sample, drop = FALSE])
        colnames(Xh) <- paste0("HUB_", hubs)
        Xf <- if (is.null(Xf)) as.matrix(Xh) else cbind(Xf, as.matrix(Xh))
    }
}
stopifnot(!is.null(Xf))

# ---- Cox LASSO ----
y <- Surv(tt, ev)
x <- as.matrix(Xf)
x[!is.finite(x)] <- 0
cv <- cv.glmnet(x, y, family = "cox", nfolds = 5, alpha = 1)
fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = cv$lambda.min)
beta <- as.matrix(coef(fit))
beta <- beta[beta[, 1] != 0, , drop = FALSE]
write.csv(data.frame(Feature = rownames(beta), Coef = beta[, 1]),
    file.path(out_dir, "CoxLASSO_selected_features.csv"),
    row.names = FALSE
)

risk <- as.numeric(x %*% beta[, 1])
df <- data.frame(Sample = cov$Sample, Risk = risk, Time = tt, Event = ev)
write.csv(df, file.path(out_dir, "Risk_scores_train.csv"), row.names = FALSE)

# ---- ارزیابی ROC زمان‌مند (۱۲/۳۶/۶۰ ماه یا معادل روز) ----
# تخمین مقیاس: اگر زمان بزرگتر از 1000، واحد احتمالاً روز است.
scal <- ifelse(max(tt, na.rm = TRUE) > 1000, 30.4, 1) # روز→ماه ~30.4
Tvec <- c(12, 36, 60) * scal
Tvec <- Tvec[Tvec < max(tt, na.rm = TRUE)]
if (length(Tvec)) {
    roc <- timeROC(T = tt, delta = ev, marker = risk, cause = 1, times = Tvec, iid = TRUE)
    write.csv(data.frame(Time = Tvec, AUC = roc$AUC), file.path(out_dir, "timeROC_AUC.csv"), row.names = FALSE)
}

# KM بر اساس tertile ریسک
grp <- cut(risk, breaks = quantile(risk, probs = c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE), include.lowest = TRUE, labels = c("Low", "Mid", "High"))
fitkm <- survfit(Surv(tt, ev) ~ grp)
p <- ggsurvplot(fitkm, risk.table = TRUE, pval = TRUE, ggtheme = theme_bw())
ggsave(file.path(out_dir, "KM_by_risk_tertiles.png"), p$plot, width = 8, height = 6, dpi = 150)
cat("Advanced survival done.\n")
