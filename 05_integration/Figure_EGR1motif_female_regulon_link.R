## =====================================================================
## EGR1 motif activity  <->  female EGR1 regulon  (multiome, integration panel)
##
## Question: do the EGR1-motif-HIGH female cells actually express the female
## EGR1 regulon that was defined ORTHOGONALLY (calling-card binding + Egr1-KD
## response)? A "yes, directionally" here is a cross-modality validation:
##   ATAC motif activity (chromVAR)  ->  cell state (cluster)  ->  RNA regulon.
##
## Within females, chromVAR EGR1 motif activity is high in a small cluster (c4,
## ~249 cells, mean z +1.6) and low in the two big clusters (c0/c1, mean z <0).
## The regulon (04_callingcard_analysis/output_full_summary/
## Egr1CC_DE_full_summary_FEMALES.csv) splits by KD direction into
##   ACTIVATED targets  (down after KD; EGR1 maintains them)  and
##   REPRESSED targets  (up after KD;   EGR1 represses them).
## Prediction: activated up in c4, repressed down in c4.
##
## Outputs (05_integration/output/):
##   Fig_EGR1motif_female_regulon_moduleScores.pdf   (violins, hi vs lo)
##   Fig_EGR1motif_female_regulon_DEoverlap_OR.pdf   (directional odds ratios)
##   female_EGR1motif_regulon_link_stats.csv         (all numbers)
##   female_c4_vs_big_DE.csv                          (the DE table)
##
## Runs locally: only reads the object (Seurat attaches fine); no RunChromVAR.
## =====================================================================
suppressMessages({ library(Seurat); library(Signac); library(ggplot2); library(patchwork) })

rds     <- "C:/Users/loril/Downloads/lot6_merged_RPCAintegrated_011524_motifsAdded_030124_chromVARadded_030124_V4_050424_2x.rds"
reg_csv <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output_full_summary/Egr1CC_DE_full_summary_FEMALES.csv"
out     <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/05_integration/output"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

grp_cols <- c(`EGR1-high (c4)` = "#D6604D", `EGR1-low (c0/c1)` = "#4A6FE3")
## presentation fonts (>=14 ticks, >=16 axis titles, >=18 title, >=14 legend)
big_theme <- theme_classic(base_size = 16) +
  theme(plot.title   = element_text(size = 18, face = "bold"),
        axis.title   = element_text(size = 16), axis.text = element_text(size = 14),
        legend.title = element_text(size = 14), legend.text = element_text(size = 14),
        strip.text   = element_text(size = 16, face = "bold"), strip.background = element_blank())

## ---- object + EGR1 motif activity + regulon sets -------------------
obj <- readRDS(rds)
obj$sex <- factor(ifelse(obj$sex %in% c("1", 1), "female", "male"), levels = c("female","male"))
mn   <- Motifs(obj[["peaks"]])@motif.names
egr1 <- names(mn)[toupper(unlist(mn)) == "EGR1"][1]
obj$egr1_act <- as.numeric(GetAssayData(obj, assay = "chromvar", slot = "data")[egr1, ])

reg <- read.csv(reg_csv, stringsAsFactors = FALSE)
DefaultAssay(obj) <- "RNA"
rna_genes <- rownames(obj[["RNA"]])
activated <- intersect(unique(reg$SYMBOL[reg$kd3_female_sig == "Down after KD"]), rna_genes)
repressed <- intersect(unique(reg$SYMBOL[reg$kd3_female_sig == "Up after KD"]),   rna_genes)

## ---- females, hi/lo groups, module scores --------------------------
fem <- subset(obj, subset = sex == "female")
fem$grp <- factor(ifelse(fem$seurat_clusters == "4", "EGR1-high (c4)",
                  ifelse(fem$seurat_clusters %in% c("0","1"), "EGR1-low (c0/c1)", NA)),
                  levels = c("EGR1-high (c4)","EGR1-low (c0/c1)"))
fem <- AddModuleScore(fem, features = list(activated, repressed),
                      name = c("act_score","rep_score"), seed = 42)   # act_score1, rep_score2
sub <- fem[, !is.na(fem$grp)]
md  <- sub@meta.data

## ---- module-score violins (hi vs lo) ------------------------------
mk_violin <- function(col, ttl, sub_n){
  p_w <- wilcox.test(md[[col]] ~ md$grp)$p.value
  ggplot(md, aes(grp, .data[[col]], fill = grp)) +
    geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.85) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = grp_cols, guide = "none") +
    labs(title = ttl,
         subtitle = sprintf("%s (n=%d genes) | Wilcoxon p = %.1e", sub_n, length(get(sub_n)), p_w),
         x = NULL, y = "module score") +
    big_theme + theme(plot.subtitle = element_text(size = 13))
}
pv <- mk_violin("act_score1", "EGR1-ACTIVATED targets", "activated") +
      mk_violin("rep_score2", "EGR1-REPRESSED targets", "repressed") +
      plot_annotation(title = "Female EGR1 regulon vs EGR1 motif activity",
                      theme = theme(plot.title = element_text(size = 18, face = "bold")))
ggsave(file.path(out, "Fig_EGR1motif_female_regulon_moduleScores.pdf"), pv, width = 11, height = 5.5)

## ---- DE c4 vs big, directional Fisher overlaps --------------------
Idents(sub) <- "grp"
de <- FindMarkers(sub, ident.1 = "EGR1-high (c4)", ident.2 = "EGR1-low (c0/c1)",
                  logfc.threshold = 0.25, min.pct = 0.1)
write.csv(de, file.path(out, "female_c4_vs_big_DE.csv"))
up   <- rownames(de)[de$avg_log2FC >  0.25 & de$p_val_adj < 0.05]
down <- rownames(de)[de$avg_log2FC < -0.25 & de$p_val_adj < 0.05]
fish <- function(hits, set){
  a <- length(intersect(hits, set)); b <- length(setdiff(hits, set))
  c_ <- length(setdiff(set, hits));  d <- length(setdiff(rna_genes, union(hits, set)))
  ft <- fisher.test(matrix(c(a, b, c_, d), 2), alternative = "greater")
  data.frame(overlap = a, OR = as.numeric(ft$estimate), p = ft$p.value,
             lo = ft$conf.int[1], hi = ft$conf.int[2]) }
orr <- rbind(
  cbind(test = "Up in c4\n-> ACTIVATED",   fish(up, activated)),
  cbind(test = "Down in c4\n-> REPRESSED", fish(down, repressed)))
orr$test <- factor(orr$test, levels = orr$test)

pbar <- ggplot(orr, aes(test, OR, fill = test)) +
  geom_col(width = 0.6, colour = "black") +
  geom_errorbar(aes(ymin = pmax(lo,0.1), ymax = hi), width = 0.2) +
  geom_text(aes(label = sprintf("OR=%.1f\np=%.1e\n(%d genes)", OR, p, overlap)),
            vjust = -0.3, size = 4.6) +
  scale_fill_manual(values = c("#D6604D", "#4A6FE3"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(title = "Directional regulon enrichment in EGR1-high cluster (c4 vs c0/c1)",
       x = NULL, y = "odds ratio (Fisher)") +
  big_theme
ggsave(file.path(out, "Fig_EGR1motif_female_regulon_DEoverlap_OR.pdf"), pbar, width = 8, height = 5.5)

## ---- stats table ---------------------------------------------------
grpstat <- function(v) c(hi_mean = mean(v[md$grp=="EGR1-high (c4)"]),
                         lo_mean = mean(v[md$grp=="EGR1-low (c0/c1)"]))
stats <- data.frame(
  set        = c("activated (down-after-KD)", "repressed (up-after-KD)"),
  n_genes    = c(length(activated), length(repressed)),
  hi_c4_mean = c(grpstat(md$act_score1)["hi_mean"], grpstat(md$rep_score2)["hi_mean"]),
  lo_big_mean= c(grpstat(md$act_score1)["lo_mean"], grpstat(md$rep_score2)["lo_mean"]),
  wilcox_p   = c(wilcox.test(md$act_score1 ~ md$grp)$p.value,
                 wilcox.test(md$rep_score2 ~ md$grp)$p.value),
  DE_overlap = orr$overlap, DE_OR = orr$OR, DE_p = orr$p)
write.csv(stats, file.path(out, "female_EGR1motif_regulon_link_stats.csv"), row.names = FALSE)
print(stats)
cat("\nUp-in-c4 activated hits:", paste(intersect(up, activated), collapse = ", "), "\n")
cat("Done ->", out, "\n")
