# 17_external_validation.R
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(survival)
    library(glmnet)
    library(timeROC)
    library(WGCNA)
})

train_dir <- Sys.getenv("GC_BASE_DIR", unset = "E:/1----------GC_29_4_mainDATAFRAME---------1")
valid_dir <- Sys.getenv("GC_VALIDATION_DIR", unset = NA)
stopifnot(!is.na(valid_dir))

# --- فایل‌های آموزش ---
ME_train_file <- file.path(train_dir, "Results_WGCNA", "WGCNA_ModuleEigengenes.csv")
coef_file <- file.path(train_dir, "Results_Survival_Adv", "CoxLASSO_selected_features.csv")
stopifnot(file.exists(ME_train_file), file.exists(coef_file))

# --- داده اعتبارسنجی ---
expr_v <- read.csv(file.path(valid_dir, "All_ComBat_matrix.csv"), row.names = 1, check.names = FALSE) %>% as.matrix()
cov_v <- read.csv(file.path(valid_dir, "All_Covariates_augmented.csv"))
stopifnot(identical(colnames(expr_v), cov_v$Sample))

# زمان/وضعیت بقا
time_cols <- c("OS_time", "OS.days", "OS_days", "SurvivalTime", "Time")
event_cols <- c("OS_event", "Status", "Event", "status")
tcol <- intersect(time_cols, names(cov_v))
ecol <- intersect(event_cols, names(cov_v))
stopifnot(length(tcol) >= 1, length(ecol) >= 1)
tt <- suppressWarnings(as.numeric(cov_v[[tcol[1]]]))
ev <- cov_v[[ecol[1]]]
if (!is.numeric(ev)) ev <- ifelse(tolower(ev) %in% c("1", "dead", "deceased", "event", "true", "yes"), 1, 0)
ok <- is.finite(tt) & !is.na(ev)
expr_v <- expr_v[, ok, drop = FALSE]
cov_v <- cov_v[ok, ]
tt <- tt[ok]
ev <- ev[ok]

# --- تولید MEs در اعتبارسنجی بر اساس رنگ‌های آموزش ---
modcol <- read.csv(file.path(train_dir, "Results_WGCNA", "WGCNA_ModuleColors.csv"))
genes_train <- modcol$Gene
common <- intersect(genes_train, rownames(expr_v))
datExpr <- t(expr_v[common, , drop = FALSE])
mod_in_valid <- modcol$Module[match(common, modcol$Gene)]
ME_v <- moduleEigengenes(datExpr, colors = mod_in_valid)$eigengenes
ME_v <- orderMEs(ME_v)
write.csv(ME_v, file.path(valid_dir, "Results_WGCNA", "WGCNA_ModuleEigengenes_VALID.csv"), row.names = TRUE)

# --- ریسک‌مدل آموزش را بخوان ---
coef_tab <- read.csv(coef_file)
feat <- coef_tab$Feature
beta <- coef_tab$Coef
names(beta) <- feat

# اگر فیچرها MEs باشند:
X <- ME_v
# اگر فیچرهای HUB_* هم داشتیم و در اعتبارسنجی می‌خواهیم، باید سازگار تولید شوند (برای سادگی اینجا فقط MEs)
X <- X[, intersect(colnames(X), names(beta)), drop = FALSE]
beta <- beta[colnames(X)]
risk <- as.numeric(as.matrix(X) %*% beta)

# ROC زمان‌مند
scal <- ifelse(max(tt, na.rm = TRUE) > 1000, 30.4, 1)
Tvec <- c(12, 36, 60) * scal
Tvec <- Tvec[Tvec < max(tt, na.rm = TRUE)]
if (length(Tvec)) {
    roc <- timeROC(T = tt, delta = ev, marker = risk, cause = 1, times = Tvec, iid = TRUE)
    write.csv(data.frame(Time = Tvec, AUC = roc$AUC), file.path(valid_dir, "Results_Survival_Adv", "timeROC_AUC_VALID.csv"), row.names = FALSE)
}
write.csv(data.frame(Sample = cov_v$Sample, Risk = risk, Time = tt, Event = ev),
    file.path(valid_dir, "Results_Survival_Adv", "Risk_scores_VALID.csv"),
    row.names = FALSE
)

cat("External validation complete.\n")
