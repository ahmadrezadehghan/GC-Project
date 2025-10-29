# 09_wgcna_modules.R
suppressPackageStartupMessages({
    library(WGCNA)
    allowWGCNAThreads()
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(patchwork)
    library(survival)
})

options(stringsAsFactors = FALSE)
base_dir <- Sys.getenv("GC_BASE_DIR", unset = "E:/1----------GC_29_4_mainDATAFRAME---------1")
expr_file <- file.path(base_dir, "All_ComBat_matrix.csv")
cov_file <- file.path(base_dir, "All_Covariates_augmented.csv")
out_dir <- file.path(base_dir, "Results_WGCNA")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(expr_file), file.exists(cov_file))

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE) %>% as.matrix()
mode(expr) <- "numeric"
cov <- read.csv(cov_file)
stopifnot(identical(colnames(expr), cov$Sample))

# --- فیلتر ژن‌های متغیر ---
topN <- min(10000L, nrow(expr))
rv <- matrixStats::rowVars(expr)
keep <- order(rv, decreasing = TRUE)[seq_len(topN)]
datExpr <- t(expr[keep, , drop = FALSE]) # نمونه‌ها×ژن‌ها
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) datExpr <- datExpr[, gsg$goodGenes]

# --- صفات کمی/کیفی ---
numize <- function(x) {
    v <- suppressWarnings(as.numeric(x))
    if (all(is.na(v))) NA else v
}
stage_to_num <- function(s) {
    s <- toupper(trimws(as.character(s)))
    s[s %in% c("I", "1")] <- 1
    s[s %in% c("II", "2")] <- 2
    s[s %in% c("III", "3")] <- 3
    s[s %in% c("IV", "4")] <- 4
    suppressWarnings(as.numeric(s))
}
traits <- data.frame(
    Tissue = ifelse(cov$Tissue %in% c("Tumor", "Normal"), ifelse(cov$Tissue == "Tumor", 1, 0), NA),
    Sex = ifelse(tolower(cov$Sex) %in% c("male", "female"), ifelse(tolower(cov$Sex) == "male", 1, 0), NA),
    Stage = stage_to_num(cov$Stage),
    LocationCode = suppressWarnings(as.numeric(cov$LocationCode)),
    LaurenCode = suppressWarnings(as.numeric(cov$LaurenCode)),
    Age = suppressWarnings(as.numeric(cov$Age)),
    row.names = cov$Sample,
    check.names = FALSE
)
traits <- traits[rownames(datExpr), , drop = FALSE]

# --- انتخاب توان ---
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0, networkType = "signed")
softPower <- if (is.finite(sft$powerEstimate)) sft$powerEstimate else 6
png(file.path(out_dir, "pickSoftThreshold.png"), 1600, 600, res = 150)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Power", ylab = "Scale Free Topology", type = "n", main = "Scale independence"
)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "blue")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Power", ylab = "Mean connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# --- ماژول‌ها ---
net <- blockwiseModules(datExpr,
    power = softPower, networkType = "signed",
    TOMType = "signed", minModuleSize = 30, reassignThreshold = 0,
    mergeCutHeight = 0.25, numericLabels = FALSE, pamRespectsDendro = TRUE,
    saveTOMs = FALSE, verbose = 2
)
moduleColors <- net$colors
MEs <- orderMEs(net$MEs)
geneInfo <- data.frame(Gene = colnames(datExpr), Module = moduleColors)
write.csv(geneInfo, file.path(out_dir, "WGCNA_ModuleColors.csv"), row.names = FALSE)
write.csv(MEs, file.path(out_dir, "WGCNA_ModuleEigengenes.csv"))

png(file.path(out_dir, "WGCNA_Dendro.png"), 1600, 800, res = 150)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
)
dev.off()

# --- ارتباط مدول‌ها با صفات ---
ME_mat <- MEs[rownames(traits), , drop = FALSE]
cors <- cor(ME_mat, traits, use = "pairwise.complete.obs")
ps <- corPvalueStudent(cors, nrow(ME_mat))
df_ct <- data.frame(
    Module = rownames(cors),
    Trait = rep(colnames(traits), each = nrow(cors)),
    Cor = as.vector(cors), P = as.vector(ps)
)
write.csv(df_ct, file.path(out_dir, "Module_Trait_correlations.csv"), row.names = FALSE)

# heatmap
library(pheatmap)
lbl <- paste0(rownames(cors), "\n(n=", table(moduleColors)[gsub("^ME", "", rownames(cors))], ")")
rownames(cors) <- lbl
png(file.path(out_dir, "Module_Trait_heatmap.png"), 1200, 900, res = 150)
pheatmap(cors, cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE, number_format = "%.2f")
dev.off()

# --- هاب‌های درون‌مدولی ---
ADJ <- adjacency(datExpr, power = softPower, type = "signed")
IMC <- intramodularConnectivity(ADJ, moduleColors)
imc_df <- data.frame(Gene = rownames(IMC), IMC, Module = moduleColors)
write.csv(imc_df, file.path(out_dir, "IntramodularConnectivity.csv"), row.names = FALSE)

# --- ارتباط با بقا (در صورت موجود) ---
cox_tab <- NULL
time_cols <- c("OS_time", "OS.days", "OS_days", "SurvivalTime", "Time")
event_cols <- c("OS_event", "Status", "Event", "status")
tcol <- intersect(time_cols, names(cov))
ecol <- intersect(event_cols, names(cov))
if (length(tcol) >= 1 && length(ecol) >= 1) {
    tt <- suppressWarnings(as.numeric(cov[[tcol[1]]]))
    ev <- cov[[ecol[1]]]
    if (!is.numeric(ev)) ev <- ifelse(tolower(ev) %in% c("1", "dead", "deceased", "event", "true", "yes"), 1, 0)
    ok <- is.finite(tt) & !is.na(ev)
    if (sum(ok) >= 50) {
        for (mn in colnames(ME_mat)) {
            d <- data.frame(t = tt[ok], e = ev[ok], x = ME_mat[ok, mn])
            fit <- try(coxph(Surv(t, e) ~ x, data = d), silent = TRUE)
            if (!inherits(fit, "try-error")) {
                s <- summary(fit)
                cox_tab <- rbind(cox_tab, data.frame(Module = mn, HR = exp(coef(fit)), P = s$wald["pvalue"]))
            }
        }
        if (!is.null(cox_tab)) write.csv(cox_tab, file.path(out_dir, "Module_Cox_univariate.csv"), row.names = FALSE)
    }
}
cat("WGCNA done.\n")
