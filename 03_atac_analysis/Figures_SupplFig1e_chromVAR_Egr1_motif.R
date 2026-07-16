## =====================================================================
## Supplementary Fig 1e - chromVAR EGR1-motif activity on the RPCA UMAP,
## faceted by sex.
##
## Plots the STORED umap.RPCA embedding + the STORED chromvar assay carried in
## the processed lot6 object (chromVAR was computed once on HTCF and saved into
## this object: ..._chromVARadded_030124...). Panel 1b (RNA) reads the SAME
## stored umap.RPCA, so 1b and 1e sit on an identical layout cell-for-cell.
## No coordinate files, no recomputation.
##
## Both sexes share ONE colour scale (a data-frame + scale_colour_gradient2),
## unlike Seurat FeaturePlot(split.by=), which rescales each facet and hides the
## male-vs-female difference.
##
## Runs locally: Seurat attaches fine for read-only; the Matrix/chromVAR
## segfault only affects RunChromVAR, which we do NOT call (the assay is stored).
## =====================================================================
suppressMessages({ library(Seurat); library(Signac); library(ggplot2) })

rds <- "C:/Users/loril/Downloads/lot6_merged_RPCAintegrated_011524_motifsAdded_030124_chromVARadded_030124_V4_050424_2x.rds"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/03_atac_analysis/output/figures/suppl_fig1"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(rds)
## sex coded 1 = female, 2 = male (Cell Ranger ARC -1/-2 aggregation order)
obj$sex <- factor(ifelse(obj$sex %in% c("1", 1), "female", "male"), levels = c("female","male"))

## EGR1 motif id from the peaks-assay motif names (MA0162.4 in JASPAR2020)
mn   <- Motifs(obj[["peaks"]])@motif.names
egr1 <- names(mn)[toupper(unlist(mn)) == "EGR1"][1]
stopifnot(!is.na(egr1))
cat("EGR1 motif id:", egr1, "\n")

emb <- Embeddings(obj, "umap.RPCA")
act <- as.numeric(GetAssayData(obj, assay = "chromvar", slot = "data")[egr1, ])
df  <- data.frame(umapRPCA_1 = emb[,1], umapRPCA_2 = emb[,2], sex = obj$sex, activity = act)
df  <- df[order(abs(df$activity)), ]                 # extreme cells drawn on top
lim <- as.numeric(ceiling(max(abs(df$activity))))    # full range -> graded red, no saturation

p <- ggplot(df, aes(umapRPCA_1, umapRPCA_2, colour = activity)) +
  geom_point(size = 0.7) + facet_wrap(~ sex) +
  scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#D6604D", midpoint = 0,
                         limits = c(-lim, lim), name = "EGR1 motif\nactivity (z)") +
  coord_fixed() +
  labs(title = paste0("EGR1 motif activity (chromVAR, ", egr1, ")"),
       x = "umapRPCA_1", y = "umapRPCA_2") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size = 14), legend.text = element_text(size = 14))
ggsave(file.path(out, "SupplFig1e_chromVAR_Egr1motif_RPCA.pdf"), p, width = 10, height = 5)

## small git-friendly per-cell table (an OUTPUT/result table, not a figure input)
write.csv(data.frame(cell = colnames(obj), sex = obj$sex,
                     umapRPCA_1 = emb[,1], umapRPCA_2 = emb[,2], EGR1_motif_activity = act),
          file.path(out, "SupplFig1e_chromVAR_Egr1motif_activity.csv"), row.names = FALSE)

cat("Done. plot + table ->", out, " | color limit +/-", lim, "\n")
