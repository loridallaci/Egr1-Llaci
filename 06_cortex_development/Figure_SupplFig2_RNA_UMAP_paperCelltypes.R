## =====================================================================
## Supplementary Fig. 2 — cortex development: RNA UMAP annotated with the
## PAPER'S OWN cell-type labels (option 3).
##
## Source = the authors' deposited CELLxGENE object (45,549 QC'd nuclei), which
## carries their `author_cell_type` calls (Fig. 1F scheme) + raw counts + sex +
## stage. We do NOT re-derive QC or re-call cell types here: cells and labels are
## the paper's. Only the UMAP embedding is re-derived, per stage, with the same
## RPCA-across-sex pipeline as the rest of Supp Fig. 2, so the geometry matches.
##
##   h5ad  : datasets.cellxgene.cziscience.com/06ee790a-...-.h5ad
##   built : build_from_h5ad.R -> Roussos_cortex_paperAnnot_seurat.rds
##
## Two outputs (NEW names; existing marker-scored plots are left untouched):
##   SupplFig2_cortex_RNA_UMAPs_paperCelltypes.pdf/.png   (re-derived, per stage)
##   SupplFig2_cortex_paperNativeUMAP_celltypes.pdf/.png  (the paper's own UMAP)
## =====================================================================
suppressMessages({ library(Seurat); library(ggplot2); library(patchwork); library(ggrepel) })
options(future.globals.maxSize = 100 * 1024^3); set.seed(1)

obj_fp <- "C:/Users/loril/Downloads/Roussos_cortex_paperAnnot_seurat.rds"
out    <- Sys.getenv("OUT_DIR", unset = "C:/Users/loril/Documents/GitHub/Egr1-Llaci/06_cortex_development/output")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

## stage display names (drop EaFet: female-only, no M-vs-F contrast — matches Supp Fig 2b)
stage_map    <- c("LaFet"="LateFetal", "Inf"="Infant", "Child"="Child", "Adol"="Adol", "Adult"="Adult")
stage_levels <- c("LateFetal","Infant","Child","Adol","Adult")

## fixed palette + legend order for the paper's 15 author_cell_type labels
ct_order <- c("RG","IPC","EN-fetal-early","EN-fetal-late","EN",
              "IN-fetal","IN-MGE","IN-CGE","OPC","Oligodendrocytes",
              "Astrocytes","Microglia","Endothelial","Pericytes","VSMC")
ct_cols <- c(RG="#8C564B", IPC="#E377C2", "EN-fetal-early"="#AEC7E8",
             "EN-fetal-late"="#3B6DB3", EN="#1F77B4", "IN-fetal"="#FF9896",
             "IN-MGE"="#D62728", "IN-CGE"="#E45756", OPC="#9467BD",
             Oligodendrocytes="#2CA02C", Astrocytes="#FF7F0E", Microglia="#17BECF",
             Endothelial="#7F7F7F", Pericytes="#4D4D4D", VSMC="#C49C94")

big_theme <- theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 16), axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 14), legend.text = element_text(size = 13))

full <- readRDS(obj_fp)
full$stage <- unname(stage_map[full$stage_raw])   # unname: else Seurat matches names to barcodes
full$author_cell_type <- factor(full$author_cell_type, levels = ct_order)

sex_cols <- c(male = "#4A6FE3", female = "#F39AC9")   # same convention as Supp Fig 2b1

## ---------- (1) re-derived per-stage RPCA-across-sex UMAPs ----------
umaps <- list(); sex_umaps <- list(); present_list <- list()
for (st in stage_levels) {
  o <- full[, which(full$stage == st)]
  cat("stage", st, ":", ncol(o), "cells,", length(unique(o$sex)), "sexes\n")
  o[["RNA"]] <- as(o[["RNA"]], "Assay5")
  o[["RNA"]] <- split(o[["RNA"]], f = o$sex)
  o <- NormalizeData(o, verbose = FALSE); o <- FindVariableFeatures(o, verbose = FALSE)
  o <- ScaleData(o, verbose = FALSE);     o <- RunPCA(o, verbose = FALSE)
  o <- IntegrateLayers(o, method = RPCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.RPCA", verbose = FALSE)
  o <- RunUMAP(o, reduction = "integrated.RPCA", dims = 1:30,
               reduction.name = "umap.rna", verbose = FALSE)
  e  <- Embeddings(o, "umap.rna")
  df <- data.frame(x = e[,1], y = e[,2],
                   celltype = factor(o$author_cell_type[rownames(e)], levels = ct_order),
                   sex = factor(o$sex[rownames(e)], levels = c("male","female")))
  cent <- aggregate(cbind(x, y) ~ celltype, df, median)
  present_list[[st]] <- as.character(unique(df$celltype))
  umaps[[st]] <- ggplot(df[sample(nrow(df)), ], aes(x, y, colour = celltype)) +
    geom_point(size = 0.35, alpha = 0.75) +
    geom_text_repel(data = cent, aes(x, y, label = celltype), inherit.aes = FALSE,
                    size = 4.2, fontface = "bold", seed = 1, box.padding = 0.4,
                    min.segment.length = 0, colour = "black", bg.color = "white", bg.r = 0.15,
                    max.overlaps = Inf) +
    coord_fixed() + labs(title = st, x = NULL, y = NULL) + big_theme
  ## same embedding, coloured by sex (companion to the existing panel b1)
  sex_umaps[[st]] <- ggplot(df[sample(nrow(df)), ], aes(x, y, colour = sex)) +
    geom_point(size = 0.35, alpha = 0.7) +
    scale_colour_manual(values = sex_cols, name = "Sex",
                        guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_fixed() + labs(title = st, x = NULL, y = NULL) + big_theme +
    theme(legend.position = if (st == stage_levels[length(stage_levels)]) "right" else "none")
  rm(o); gc()
}
## by-sex montage (paper's cells, same embedding as the cell-type montage)
sex_mont <- wrap_plots(sex_umaps[stage_levels], nrow = 1) +
  plot_annotation(title = "Cortex development — RNA UMAP (RPCA, integrated across sex), coloured by sex",
                  subtitle = "Paper's cells (CELLxGENE 45,549 nuclei); same embedding as the cell-type panel",
                  theme = theme(plot.title = element_text(size = 18, face = "bold"),
                                plot.subtitle = element_text(size = 13)))
ggsave(file.path(out, "SupplFig2_cortex_RNA_UMAPs_paperCells_bySex.pdf"), sex_mont, width = 5.2*5+2.5, height = 6.0, limitsize = FALSE)
ggsave(file.path(out, "SupplFig2_cortex_RNA_UMAPs_paperCells_bySex.png"), sex_mont, width = 5.2*5+2.5, height = 6.0, dpi = 150, limitsize = FALSE)
present_all <- ct_order[ct_order %in% unlist(present_list)]
ct_scale <- scale_colour_manual(values = ct_cols, limits = present_all, drop = FALSE,
                                name = "Cell type (paper)",
                                guide = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
umaps_s <- lapply(seq_along(stage_levels), function(i) {
  p <- umaps[[stage_levels[i]]] + ct_scale
  p + theme(legend.position = if (i == length(stage_levels)) "right" else "none")
})
mont <- wrap_plots(umaps_s, nrow = 1) +
  plot_annotation(title = "Cortex development — RNA UMAP (RPCA, integrated across sex), the paper's cell-type labels",
                  subtitle = "Cells & labels: authors' deposited annotation (Zhu et al., Sci Adv 2023; CELLxGENE 45,549 nuclei)",
                  theme = theme(plot.title = element_text(size = 18, face = "bold"),
                                plot.subtitle = element_text(size = 13)))
ggsave(file.path(out, "SupplFig2_cortex_RNA_UMAPs_paperCelltypes.pdf"), mont, width = 5.2*5+2.5, height = 6.0, limitsize = FALSE)
ggsave(file.path(out, "SupplFig2_cortex_RNA_UMAPs_paperCelltypes.png"), mont, width = 5.2*5+2.5, height = 6.0, dpi = 150, limitsize = FALSE)

## ---------- (2) the paper's OWN global UMAP, coloured by author_cell_type ----------
pe <- Embeddings(full, "paperUMAP")
pdf_df <- data.frame(x = pe[,1], y = pe[,2],
                     celltype = factor(full$author_cell_type[rownames(pe)], levels = ct_order))
cent2 <- aggregate(cbind(x, y) ~ celltype, pdf_df, median)
pnat <- ggplot(pdf_df[sample(nrow(pdf_df)), ], aes(x, y, colour = celltype)) +
  geom_point(size = 0.25, alpha = 0.6) +
  geom_text_repel(data = cent2, aes(x, y, label = celltype), inherit.aes = FALSE,
                  size = 4.5, fontface = "bold", seed = 1, box.padding = 0.5,
                  min.segment.length = 0, colour = "black", bg.color = "white", bg.r = 0.15, max.overlaps = Inf) +
  scale_colour_manual(values = ct_cols, name = "Cell type (paper)",
                      guide = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  coord_fixed() + labs(title = "Paper's own UMAP (all stages), authors' cell-type labels",
                       x = "UMAP 1", y = "UMAP 2") + big_theme +
  theme(axis.title = element_text(size = 16))
ggsave(file.path(out, "SupplFig2_cortex_paperNativeUMAP_celltypes.pdf"), pnat, width = 11, height = 9)
ggsave(file.path(out, "SupplFig2_cortex_paperNativeUMAP_celltypes.png"), pnat, width = 11, height = 9, dpi = 150)

## cell-count table for source data
tab <- as.data.frame.matrix(table(full$author_cell_type, full$stage))
write.csv(tab, file.path(out, "paperCelltype_counts_byStage.csv"))
cat("\nDone ->", out, "\n")
