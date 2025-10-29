# 14_integration_network_dge_gsva.R
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(igraph)
    library(purrr)
    library(pheatmap)
})

base_dir <- Sys.getenv("GC_BASE_DIR", unset = "E:/1----------GC_29_4_mainDATAFRAME---------1")
expr_file <- file.path(base_dir, "All_ComBat_matrix.csv")
dge_file <- file.path(base_dir, "Results_DGE", "DGE_global_Tumor_vs_Normal.csv")
wgcna_dir <- file.path(base_dir, "Results_WGCNA")
gsva_dir <- file.path(base_dir, "Results_GSVA")
net_dir <- file.path(base_dir, "Results_Networks_Variants")
out_dir <- file.path(base_dir, "Results_Integration")
dir.create(out_dir, TRUE, TRUE)

stopifnot(file.exists(expr_file), file.exists(dge_file))

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE) %>% as.matrix()
mode(expr) <- "numeric"
dge <- read.csv(dge_file, row.names = 1, check.names = FALSE)
dge$Gene <- rownames(dge)

# ---- 1) Hub-Driver: DE ∩ هاب‌های شبکه ----
pick_net <- function(method = "aracne", group = "Tumor") {
    f <- file.path(net_dir, sprintf("%s_%s_edges.csv", method, group))
    if (file.exists(f)) read.csv(f) else NULL
}
edges <- pick_net("aracne", "Tumor")
if (!is.null(edges) && nrow(edges)) {
    g <- graph_from_data_frame(edges[, 1:2], directed = FALSE)
    cent <- data.frame(
        Gene = names(V(g)),
        degree = degree(g),
        betweenness = betweenness(g, directed = FALSE, normalized = TRUE),
        closeness = closeness(g, normalized = TRUE),
        stringsAsFactors = FALSE
    )
    # هاب‌ها = بالاترین 5%
    thr <- quantile(cent$degree, 0.95, na.rm = TRUE)
    hubs <- cent %>% filter(degree >= thr)
    write_csv(cent, file.path(out_dir, "Tumor_ARACNE_centrality.csv"))
    write_csv(hubs, file.path(out_dir, "Tumor_ARACNE_hubs.csv"))

    # همپوشانی با DE (adj.P.Val یا P.Value)
    sigDE <- dge %>%
        mutate(adj = ifelse(!is.na(adj.P.Val), adj.P.Val, p.adjust(P.Value, "BH"))) %>%
        filter(adj < 0.05)
    hub_drivers <- hubs %>% inner_join(sigDE, by = c("Gene" = "Gene"))
    write_csv(hub_drivers, file.path(out_dir, "HubDrivers_Tumor_ARACNE.csv"))
}

# ---- 2) WGCNA × GSVA: همبستگی ME با pathway scores ----
ME_file <- file.path(wgcna_dir, "WGCNA_ModuleEigengenes.csv")
if (file.exists(ME_file)) {
    MEs <- read.csv(ME_file, row.names = 1, check.names = FALSE)
    # بارگذاری GSVA (لیست RDS یا CSVهای جدا)
    gsva_rds <- file.path(gsva_dir, "GSVA_scores_list.rds")
    gs_list <- list()
    if (file.exists(gsva_rds)) {
        gs_list <- readRDS(gsva_rds)
    } else {
        csvs <- list.files(gsva_dir, "^GSVA_scores_.*\\.csv$", full.names = TRUE)
        for (cf in csvs) {
            nm <- sub("^GSVA_scores_(.*)\\.csv$", "\\1", basename(cf))
            gs_list[[nm]] <- as.matrix(read.csv(cf, row.names = 1, check.names = FALSE))
        }
    }
    if (length(gs_list)) {
        cors_all <- list()
        for (nm in names(gs_list)) {
            sc <- gs_list[[nm]]
            common <- intersect(colnames(sc), rownames(MEs))
            if (length(common) < 10) next
            me <- MEs[common, , drop = FALSE]
            sc <- sc[, common, drop = FALSE]
            C <- cor(me, t(sc), use = "pairwise.complete.obs")
            df <- reshape2::melt(C)
            colnames(df) <- c("Module", "Pathway", "Correlation")
            df$Collection <- nm
            cors_all[[nm]] <- df
        }
        if (length(cors_all)) {
            cors_all <- bind_rows(cors_all)
            write_csv(cors_all, file.path(out_dir, "ME_GSVA_correlations.csv"))
            # heatmap جمع‌وجور از 30 مسیر برتر هر کالکشن
            top_pw <- cors_all %>%
                group_by(Collection, Pathway) %>%
                summarize(MaxAbs = max(abs(Correlation)), .groups = "drop") %>%
                arrange(desc(MaxAbs)) %>%
                group_by(Collection) %>%
                slice_head(n = 30) %>%
                pull(Pathway) %>%
                unique()
            last_coll <- tail(names(gs_list), 1)
            sc <- gs_list[[last_coll]]
            common <- intersect(colnames(sc), rownames(MEs))
            H <- cor(MEs[common, , drop = FALSE], t(sc[top_pw, common, drop = FALSE]))
            png(file.path(out_dir, "ME_GSVA_heatmap_top30.png"), 1400, 900, res = 150)
            pheatmap(H, cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = FALSE)
            dev.off()
        }
    }
}
cat("Integration done.\n")
