## =====================================================================
## Supplementary Fig 1c - ATAC-seq UMAP (MACS2 peaks), rLSI-INTEGRATED.
##   Matched analog of the RPCA-integrated RNA panel (Figures_SupplFig1b_RPCA_UMAP.R):
##   RNA integrates across sex with reciprocal PCA; ATAC integrates across sex
##   with reciprocal LSI (rLSI) - same anchor framework, LSI in place of PCA.
##
##   peaks assay -> per-sex TFIDF/FindTopFeatures/RunSVD ->
##   FindIntegrationAnchors(reduction = "rlsi") -> IntegrateEmbeddings ->
##   RunUMAP(reduction = "integrated_lsi").  LSI component 1 dropped (depth).
##
## Note: integrating across sex removes the sex signal by design - these are
## QC / cell-type UMAPs (both sexes present, cells co-cluster by type), not a
## readout of sex biology (that comes from DARs / RENIN / chromVAR).
## Runs on HTCF (Signac/Seurat v5 + data).
## =====================================================================
suppressMessages({library(Seurat); library(Signac); library(ggplot2)})

rds <- "C:/Users/loril/Documents/data/multiome_081721_aggregated_data_analysis/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks.rds"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/03_atac_analysis/output/figures/suppl_fig1"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
sex_cols <- c(female = "#F39AC9", male = "#4A6FE3")

lot6 <- readRDS(rds)
lot6$sex <- ifelse(grepl("-1$", colnames(lot6)), "female",
             ifelse(grepl("-2$", colnames(lot6)), "male", NA))

## ---- global LSI (provides the reduction to integrate) ---------------
DefaultAssay(lot6) <- "peaks"
lot6 <- RunTFIDF(lot6)
lot6 <- FindTopFeatures(lot6, min.cutoff = 5)
lot6 <- RunSVD(lot6)

## ---- reciprocal LSI integration across sex --------------------------
obj.list <- SplitObject(lot6, split.by = "sex")
obj.list <- lapply(obj.list, function(x) {
  x <- RunTFIDF(x)
  x <- FindTopFeatures(x, min.cutoff = 5)
  x <- RunSVD(x)
  x
})

anchors <- FindIntegrationAnchors(
  object.list     = obj.list,
  anchor.features = rownames(lot6),   # peaks as features
  reduction       = "rlsi",           # ATAC analog of "rpca"
  dims            = 2:30              # drop LSI-1 (sequencing-depth component)
)

integrated <- IntegrateEmbeddings(
  anchorset          = anchors,
  reductions         = lot6[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate  = 1:30
)

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30,
                      reduction.name = "umap.atac", reduction.key = "umapATAC_")

## ---- panel c : rLSI-integrated ATAC UMAP, colored by sex ------------
p <- DimPlot(integrated, reduction = "umap.atac", group.by = "sex", cols = sex_cols) +
  ggtitle("ATAC-seq UMAP (rLSI-integrated)")
ggsave(file.path(out, "SupplFig1c_ATAC_UMAP_rLSI.pdf"), p, width = 6, height = 5)

p_split <- DimPlot(integrated, reduction = "umap.atac", group.by = "sex",
                   split.by = "sex", cols = sex_cols) + ggtitle("ATAC rLSI UMAP, by sex")
ggsave(file.path(out, "SupplFig1c_ATAC_UMAP_rLSI_bySex.pdf"), p_split, width = 10, height = 5)

cat("Done. rLSI-integrated panel 1c in:\n  ", out, "\n")
