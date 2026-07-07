## =====================================================================
## Supplementary Fig 1b - RNA UMAP from the RPCA-INTEGRATED sample.
## Reproduces 02_rna_analysis/multiome_rna/rpca_integration.R:
##   split RNA by sex -> Normalize/HVG/Scale/PCA -> IntegrateLayers(RPCA)
##   -> RunUMAP(reduction = "integrated.RPCA", reduction.name = "umap.RPCA")
## Saves the RPCA UMAP (panel b, by sex) + the per-sex split (panel e style)
## as PDFs into git. Runs locally (Seurat v5).
## =====================================================================
suppressMessages({library(Seurat); library(Signac); library(ggplot2); library(dplyr)})

rds <- "C:/Users/loril/Documents/data/multiome_081721_aggregated_data_analysis/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks.rds"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/multiome_rna/output/figures/suppl_fig1"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
sex_cols <- c(female = "#F39AC9", male = "#4A6FE3")

lot6 <- readRDS(rds)
lot6$sex <- ifelse(grepl("-1$", colnames(lot6)), "female",
             ifelse(grepl("-2$", colnames(lot6)), "male", NA))
DefaultAssay(lot6) <- "RNA"

## ---- Seurat v5 RPCA integration (split layers by sex) --------------
lot6[["RNA"]] <- split(lot6[["RNA"]], f = lot6$sex)
lot6 <- NormalizeData(lot6)
lot6 <- FindVariableFeatures(lot6)
lot6 <- ScaleData(lot6)
lot6 <- RunPCA(lot6, verbose = FALSE)
lot6 <- IntegrateLayers(lot6, method = RPCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.RPCA",
                        verbose = FALSE)
lot6 <- RunUMAP(lot6, reduction = "integrated.RPCA", dims = 1:30,
                reduction.name = "umap.RPCA", reduction.key = "umapRPCA_")
lot6[["RNA"]] <- JoinLayers(lot6[["RNA"]])   # rejoin for plotting/FeaturePlot

## ---- panel b : RPCA UMAP colored by sex ---------------------------
p_b <- DimPlot(lot6, reduction = "umap.RPCA", group.by = "sex", cols = sex_cols) +
  ggtitle("RNA UMAP (RPCA-integrated)")
ggsave(file.path(out, "SupplFig1b_RNA_UMAP_RPCA.pdf"), p_b, width = 6, height = 5)

## ---- panel b/e : per-sex split + Egr1 FeaturePlot on RPCA UMAP -----
p_split <- DimPlot(lot6, reduction = "umap.RPCA", group.by = "sex",
                   split.by = "sex", cols = sex_cols) + ggtitle("RPCA UMAP, by sex")
ggsave(file.path(out, "SupplFig1b_RNA_UMAP_RPCA_bySex.pdf"), p_split, width = 10, height = 5)

p_e <- FeaturePlot(lot6, features = "Egr1", reduction = "umap.RPCA",
                   split.by = "sex", order = TRUE) & theme(legend.position = "right")
ggsave(file.path(out, "SupplFig1e_FeatureUMAP_Egr1_RPCA.pdf"), p_e, width = 10, height = 5)

cat("Done. RPCA UMAP panels in:\n  ", out, "\n")
