# Need seurat version 5 

suppressMessages(library(Seurat)) 
suppressMessages(library(Signac))
suppressMessages(library(SeuratWrappers)) 
suppressMessages(library(RENIN)) 
suppressMessages(library(harmony))
#plan(sequential)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)

lot6 = readRDS("data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds")  # original: "/home/lllaci/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds"
obj <- lot6
DefaultAssay(obj) <- 'RNA'

# Normalize 
obj <- NormalizeData(obj)
# Find variable features
obj <- FindVariableFeatures(obj)
# Scale data (important!)
obj <- ScaleData(obj)

# PCA -> neighbors -> clusters -> UMAP
obj <- RunPCA(obj, features = VariableFeatures(obj))
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.RPCA",
  verbose = FALSE
)

obj <- FindNeighbors(obj, reduction = "integrated.RPCA", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "RPCA_clusters")
obj <- RunUMAP(obj, reduction = "integrated.RPCA", dims = 1:30, reduction.name = "umap.RPCA")

pdf(paste("lot6_RPCA_integrated_umap_011524.pdf", sep = ""))
DimPlot(
  obj,
  reduction = "umap.RPCA",
  group.by = c("sex", "predicted.celltype.l2", "RPCA_clusters"),
  combine = FALSE, label.size = 2
)
dev.off()


