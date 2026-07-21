## =====================================================================
## Supplementary Fig. 2 — cortex development: RNA UMAP (RPCA-integrated
## across sex, same geometry as Figure_SupplFig2_cellnumbers_umaps.R) with
## clusters ANNOTATED by cell type.
##
## Cell types + marker genes are taken from the source dataset paper:
##   Zhu K, Bendl J, Rahman S, ... Roussos P. "Multi-omic profiling of the
##   developing human cerebral cortex at the single-cell level." Sci Adv
##   2023;9(41):eadg3754. Cell-type marker panel = their Fig. 1F. The paper's
##   fine subtypes (EN-fetal-early/late, IN-MGE/CGE, Endo/Peri/VSMC) are
##   collapsed into 9 major classes here for a legible UMAP overlay.
##
## Pipeline per stage: split RNA by sex batch -> Normalize/HVG/Scale/PCA ->
## IntegrateLayers(RPCA) -> UMAP + FindClusters on integrated.RPCA (dims 1:30)
## -> AddModuleScore for each cell-type marker set -> assign each cluster to
## the cell type with the highest mean module score.
##
## Outputs (06_cortex_development/output/):
##   SupplFig2_cortex_RNA_UMAPs_celltypes.pdf     (UMAP montage, coloured by cell type)
##   SupplFig2_cortex_celltype_marker_dotplots.pdf (marker dotplot per stage)
##   celltype_cluster_assignment_<stage>.csv       (cluster -> module scores -> label)
##
## NOTE: only LateFetal/Infant/Child run locally (pre-motif objects in Downloads).
##   Adol/Adult objects are on the .181 workstation; set CORTEX_FULL=1 + DATA_DIR
##   there (withMotifs_V4 files) to render all 5.
## =====================================================================
suppressMessages({ library(Seurat); library(Signac); library(ggplot2); library(patchwork) })
options(future.globals.maxSize = 400 * 1024^3)
set.seed(1)

out <- Sys.getenv("OUT_DIR", unset = "C:/Users/loril/Documents/GitHub/Egr1-Llaci/06_cortex_development/output")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

## ---- paper cell-type marker panel (Zhu/Roussos 2023, Fig. 1F) --------
## 9 major classes; markers restricted to the paper's panel (+ a few canonical
## synonyms in the same class). Order = fetal/progenitor -> neuronal -> glia -> vascular.
markers <- list(
  RG    = c("HES5","VIM","PAX6","SOX2"),                # radial glia
  IPC   = c("EOMES","PPP1R17"),                         # intermediate progenitors
  ExN   = c("SATB2","SLC17A7","NEUROD2","NEUROD6"),     # excitatory neurons
  InN   = c("GAD1","GAD2","DLX2","LHX6"),               # inhibitory neurons (MGE/CGE merged)
  OPC   = c("OLIG1","OLIG2","SOX10","PDGFRA"),          # oligo precursors
  Oligo = c("MOBP","OPALIN","MBP","PLP1"),              # oligodendrocytes
  Astro = c("AQP4","GFAP","SLC1A3"),                    # astrocytes
  Micro = c("PTPRC","CX3CR1","CSF1R","P2RY12"),         # microglia
  Vasc  = c("CLDN5","FLT1","PDGFRB","COL1A2")           # endothelial/pericyte/VSMC
)
ct_levels <- names(markers)
ct_cols <- c(RG="#8C564B", IPC="#E377C2", ExN="#1F77B4", InN="#D62728",
             OPC="#9467BD", Oligo="#2CA02C", Astro="#FF7F0E", Micro="#17BECF",
             Vasc="#7F7F7F")

## ---- stage inputs ----------------------------------------------------
if (Sys.getenv("CORTEX_FULL") == "1") {
  data_dir <- Sys.getenv("DATA_DIR", unset = "/home/lllaci/data/cortex_development")
  stage_files <- c(
    LateFetal = "LaFet_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
    Infant    = "Inf_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
    Child     = "Child_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
    Adol      = "Adol_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
    Adult     = "Adult_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds")
} else {
  data_dir <- "C:/Users/loril/Downloads"
  stage_files <- c(
    LateFetal = "LaFet_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds",
    Infant    = "Inf_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds",
    Child     = "Child_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds")
}
stage_levels <- names(stage_files)

big_theme <- theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.title = element_text(size = 14), legend.text = element_text(size = 14))

umaps <- list(); dots <- list(); present_list <- list()
for (st in names(stage_files)) {
  fp <- file.path(data_dir, stage_files[[st]])
  if (!file.exists(fp)) { cat("SKIP (missing):", st, "\n"); next }
  cat("\n==== stage:", st, "====\n"); o <- readRDS(fp)
  for (a in intersect(c("peaks","ATAC"), Assays(o)))
    suppressWarnings(try(Fragments(o[[a]]) <- NULL, silent = TRUE))

  ## split by the -1/-2 barcode batch (== sex; a bijection per stage). Sex is
  ## NOT displayed here, so the per-stage suffix->sex flip is irrelevant: we only
  ## need the batch to reproduce the RPCA-across-sex integration geometry.
  o$batch <- factor(sub(".*-", "", colnames(o)))

  DefaultAssay(o) <- "RNA"
  if (!inherits(o[["RNA"]], "Assay5")) o[["RNA"]] <- as(o[["RNA"]], "Assay5")
  o[["RNA"]] <- split(o[["RNA"]], f = o$batch)
  o <- NormalizeData(o, verbose = FALSE); o <- FindVariableFeatures(o, verbose = FALSE)
  o <- ScaleData(o, verbose = FALSE);     o <- RunPCA(o, verbose = FALSE)
  o <- IntegrateLayers(o, method = RPCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.RPCA", verbose = FALSE)
  o <- RunUMAP(o, reduction = "integrated.RPCA", dims = 1:30,
               reduction.name = "umap.rna", verbose = FALSE)
  o <- FindNeighbors(o, reduction = "integrated.RPCA", dims = 1:30, verbose = FALSE)
  o <- FindClusters(o, resolution = 0.6, verbose = FALSE)
  o[["RNA"]] <- JoinLayers(o[["RNA"]])

  ## ---- module score per cell-type marker set, then argmax per cluster ----
  o <- AddModuleScore(o, features = markers, name = "ctscore", assay = "RNA",
                      nbin = 20, ctrl = 50, seed = 1)
  score_cols <- paste0("ctscore", seq_along(markers))
  clus <- as.character(o$seurat_clusters)                 # length ncells, cell order
  meanscore <- sapply(score_cols, function(sc) tapply(o[[sc]][,1], clus, mean))
  colnames(meanscore) <- ct_levels
  assign_ct <- ct_levels[max.col(meanscore, ties.method = "first")]  # per cluster
  names(assign_ct) <- rownames(meanscore)
  ## assign per cell; strip names so Seurat maps by position, not by barcode
  o$celltype <- factor(unname(assign_ct[clus]), levels = ct_levels)

  ## record the assignment for source data
  ntab <- table(clus)
  asg <- data.frame(cluster = rownames(meanscore),
                    n_cells = as.integer(ntab[rownames(meanscore)]),
                    celltype = unname(assign_ct[rownames(meanscore)]),
                    round(meanscore, 3), check.names = FALSE)
  write.csv(asg, file.path(out, paste0("celltype_cluster_assignment_", st, ".csv")), row.names = FALSE)
  cat("clusters ->", paste(sprintf("%s:%s", rownames(meanscore), as.character(assign_ct)), collapse="  "), "\n")

  ## ---- UMAP coloured by cell type, labelled at cell-type centroids ----
  e <- Embeddings(o, "umap.rna")
  df <- data.frame(x = e[,1], y = e[,2], celltype = o$celltype[rownames(e)])
  cent <- aggregate(cbind(x, y) ~ celltype, df, median)
  present_ct <- levels(droplevels(df$celltype))
  present_list[[st]] <- present_ct
  umaps[[st]] <- ggplot(df[sample(nrow(df)), ], aes(x, y, colour = celltype)) +
    geom_point(size = 0.35, alpha = 0.75) +
    ggrepel::geom_text_repel(data = cent, aes(x, y, label = celltype), inherit.aes = FALSE,
                             size = 5, fontface = "bold", seed = 1, box.padding = 0.4,
                             min.segment.length = 0, colour = "black", bg.color = "white", bg.r = 0.15) +
    ## colour scale is applied once, at montage time (below), so patchwork
    ## collects a single unified legend across all panels.
    coord_fixed() + labs(title = st, x = NULL, y = NULL) + big_theme +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  ## ---- marker dotplot: paper markers x assigned cell type -------------
  Idents(o) <- "celltype"
  feats <- unique(unlist(markers))
  dp <- DotPlot(o, features = feats, assay = "RNA") +
    coord_flip() +
    scale_colour_gradient(low = "grey85", high = "#B2182B") +
    labs(title = st, x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 13), legend.text = element_text(size = 12))
  dots[[st]] <- dp
  rm(o); gc()
}

ord <- intersect(stage_levels, names(umaps))
if (length(ord)) {
  present_all <- ct_levels[ct_levels %in% unlist(present_list)]  # union, in canonical order
  ## Apply ONE identical colour scale to every panel, then collect. Doing this in
  ## an explicit lapply (rather than via `&`) guarantees byte-identical guide grobs,
  ## which is what patchwork requires to merge them into a single legend.
  ct_scale <- scale_colour_manual(values = ct_cols, limits = present_all, drop = FALSE,
                                  name = "Cell type",
                                  guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
  ## Same colour scale on every panel; show the legend only on the last panel
  ## (which contains all cell types) — avoids patchwork's flaky legend collection.
  umaps_s <- lapply(seq_along(ord), function(i) {
    p <- umaps[[ord[i]]] + ct_scale
    p + theme(legend.position = if (i == length(ord)) "right" else "none")
  })
  mont <- wrap_plots(umaps_s, nrow = 1) +
    plot_annotation(title = "Cortex development — RNA UMAP (RPCA, integrated across sex), clusters annotated by cell type",
                    subtitle = "Cell types & markers per Zhu et al., Sci Adv 2023 (Fig. 1F)",
                    theme = theme(plot.title = element_text(size = 18, face = "bold"),
                                  plot.subtitle = element_text(size = 14)))
  ggsave(file.path(out, "SupplFig2_cortex_RNA_UMAPs_celltypes.pdf"), mont,
         width = 5.2*length(ord)+2, height = 6.0, limitsize = FALSE)
  ggsave(file.path(out, "SupplFig2_cortex_RNA_UMAPs_celltypes.png"), mont,
         width = 5.2*length(ord)+2, height = 6.0, dpi = 150, limitsize = FALSE)

  dmont <- wrap_plots(dots[ord], nrow = 1)
  ggsave(file.path(out, "SupplFig2_cortex_celltype_marker_dotplots.pdf"), dmont,
         width = 5.5*length(ord), height = 8, limitsize = FALSE)
  ggsave(file.path(out, "SupplFig2_cortex_celltype_marker_dotplots.png"), dmont,
         width = 5.5*length(ord), height = 8, dpi = 150, limitsize = FALSE)
}
cat("\nStages built:", paste(ord, collapse=", "), "\nDone ->", out, "\n")
