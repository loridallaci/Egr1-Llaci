## =====================================================================
## Supplementary Fig. 1 panels f-g: the small female EGR1-motif-HIGH subcluster
## (c4) and its reactivation of the female EGR1 regulon.
##
##  (f) EGR1 motif activity (chromVAR z) across female clusters -> c4 is the
##      small high population; c0/c1 are the big low clusters.
##  (g) Female EGR1 regulon module scores (activated targets / repressed targets)
##      in c4 vs c0/c1 -> activated up, repressed down = motif activity drives
##      the independently-defined (calling-card + KD) regulon.
##
## Regulon: 04_callingcard_analysis/output_full_summary/Egr1CC_DE_full_summary_FEMALES.csv
##   activated = down-after-KD (EGR1 maintains); repressed = up-after-KD (EGR1 represses).
##
## Runs locally: only reads the stored object (no RunChromVAR).
## Outputs (05_integration/output/):
##   SupplFig1f_EGR1motif_by_female_cluster.pdf
##   SupplFig1g_female_regulon_moduleScores.pdf
##   SupplFig1fg_c4_regulon.pdf   (combined f|g montage)
## =====================================================================
suppressMessages({ library(Seurat); library(Signac); library(ggplot2); library(patchwork) })
rds     <- "C:/Users/loril/Downloads/lot6_merged_RPCAintegrated_011524_motifsAdded_030124_chromVARadded_030124_V4_050424_2x.rds"
reg_csv <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output_full_summary/Egr1CC_DE_full_summary_FEMALES.csv"
out     <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/05_integration/output"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

big_theme <- theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), legend.title = element_text(size = 14),
        legend.text = element_text(size = 14), strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank())

obj <- readRDS(rds)
obj$sex <- factor(ifelse(obj$sex %in% c("1", 1), "female", "male"), levels = c("female","male"))
mn <- Motifs(obj[["peaks"]])@motif.names
egr1 <- names(mn)[toupper(unlist(mn)) == "EGR1"][1]
obj$egr1_act <- as.numeric(GetAssayData(obj, assay = "chromvar", slot = "data")[egr1, ])
DefaultAssay(obj) <- "RNA"
fem <- subset(obj, subset = sex == "female")

## ---- panel f: EGR1 motif activity across female clusters -----------
keep <- names(which(table(fem$seurat_clusters) >= 20))          # drop near-empty clusters
fdf <- data.frame(cluster = factor(fem$seurat_clusters, levels = keep),
                  egr1 = fem$egr1_act)[fem$seurat_clusters %in% keep, ]
role <- ifelse(levels(fdf$cluster) == "4", "EGR1-high (c4)",
        ifelse(levels(fdf$cluster) %in% c("0","1"), "EGR1-low big (c0/c1)", "other"))
fdf$role <- factor(role[match(fdf$cluster, levels(fdf$cluster))],
                   levels = c("EGR1-high (c4)","EGR1-low big (c0/c1)","other"))
pf <- ggplot(fdf, aes(cluster, egr1, fill = role)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_violin(scale = "width", colour = NA, alpha = 0.9) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("EGR1-high (c4)" = "#D6604D", "EGR1-low big (c0/c1)" = "#4A6FE3",
                               "other" = "grey75"), name = NULL) +
  labs(title = "EGR1 motif activity by female cluster",
       x = "female cluster", y = "EGR1 motif activity (chromVAR z)") +
  big_theme + theme(legend.position = "top")
ggsave(file.path(out, "SupplFig1f_EGR1motif_by_female_cluster.pdf"), pf, width = 7, height = 5.5)

## ---- panel g: regulon module scores, c4 vs big --------------------
reg <- read.csv(reg_csv, stringsAsFactors = FALSE)
rna_genes <- rownames(fem[["RNA"]])
activated <- intersect(unique(reg$SYMBOL[reg$kd3_female_sig == "Down after KD"]), rna_genes)
repressed <- intersect(unique(reg$SYMBOL[reg$kd3_female_sig == "Up after KD"]),   rna_genes)
fem <- AddModuleScore(fem, features = list(activated, repressed), name = c("act","rep"), seed = 42)
fem$grp <- factor(ifelse(fem$seurat_clusters == "4", "EGR1-high (c4)",
                  ifelse(fem$seurat_clusters %in% c("0","1"), "EGR1-low (c0/c1)", NA)),
                  levels = c("EGR1-high (c4)","EGR1-low (c0/c1)"))
sub <- fem[, !is.na(fem$grp)]; md <- sub@meta.data
gdf <- rbind(data.frame(grp = md$grp, score = md$act1, set = sprintf("Activated targets\n(n=%d)", length(activated))),
             data.frame(grp = md$grp, score = md$rep2, set = sprintf("Repressed targets\n(n=%d)", length(repressed))))
pact <- signif(wilcox.test(md$act1 ~ md$grp)$p.value, 2)
prep <- signif(wilcox.test(md$rep2 ~ md$grp)$p.value, 2)
pg <- ggplot(gdf, aes(grp, score, fill = grp)) +
  geom_violin(scale = "width", colour = NA, alpha = 0.9) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white") +
  facet_wrap(~ set) +
  scale_fill_manual(values = c("EGR1-high (c4)" = "#D6604D", "EGR1-low (c0/c1)" = "#4A6FE3"), guide = "none") +
  labs(title = "Female EGR1 regulon expression", x = NULL, y = "module score") +
  big_theme + theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(out, "SupplFig1g_female_regulon_moduleScores.pdf"), pg, width = 8, height = 5.5)

## ---- combined f | g montage ---------------------------------------
mont <- pf + pg + plot_layout(widths = c(1, 1.15)) + plot_annotation(tag_levels = list(c("f","g")))
ggsave(file.path(out, "SupplFig1fg_c4_regulon.pdf"), mont, width = 15, height = 5.5)

cat(sprintf("panel f/g done. activated=%d repressed=%d | Wilcoxon act p=%.2g rep p=%.2g\n",
            length(activated), length(repressed), pact, prep))
cat("->", out, "\n")
