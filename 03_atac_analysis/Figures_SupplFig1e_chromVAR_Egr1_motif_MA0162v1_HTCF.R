## =====================================================================
## Supplementary Fig 1e - chromVAR EGR1-motif activity on the RPCA UMAP,
## split by sex.  [HTCF version, MA0162.1 - OLD EGR1 matrix]
##
## Same as Figures_SupplFig1e_chromVAR_Egr1_motif_HTCF.R but uses the OLD
## EGR1 matrix MA0162.1 instead of the latest MA0162.4. The default
## getMatrixSet(..., all_versions = FALSE) only returns the newest version,
## so MA0162.1 is fetched explicitly by full ID and added as the only motif.
##
## Purpose: check the EGR1 motif-activity readout is robust to the JASPAR
## matrix version (MA0162 is EGR1 in both; versions differ only in the
## training data, and the EGR1 GC-rich site is conserved).
##
## Paths from env vars (set by the sbatch wrapper):
##   MULTIOME_RDS - lot6 multiome Seurat object (RNA + MACS2 `peaks` assays)
##   OUTDIR       - where to write the PDF + CSV
##
## Outputs (kept separate from the MA0162.4 run):
##   SupplFig1e_chromVAR_Egr1motif_MA0162.1_RPCA.pdf
##   SupplFig1e_chromVAR_Egr1motif_MA0162.1_activity.csv
## =====================================================================
suppressMessages({
  library(Seurat); library(Signac); library(ggplot2); library(dplyr)
  library(JASPAR2020); library(TFBSTools)
  library(BSgenome.Mmusculus.UCSC.mm10); library(chromVAR)
})

rds <- Sys.getenv("MULTIOME_RDS")
out <- Sys.getenv("OUTDIR", unset = ".")
stopifnot("Set MULTIOME_RDS to the multiome .rds path" = nzchar(rds), file.exists(rds))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
sex_cols <- c(female = "#F39AC9", male = "#4A6FE3")

obj <- readRDS(rds)
obj$sex <- ifelse(grepl("-1$", colnames(obj)), "female",
            ifelse(grepl("-2$", colnames(obj)), "male", NA))

## ---- 1. RPCA-integrated RNA UMAP (umap.RPCA) -----------------------
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sex)
obj <- NormalizeData(obj); obj <- FindVariableFeatures(obj); obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- IntegrateLayers(obj, method = RPCAIntegration,
                       orig.reduction = "pca", new.reduction = "integrated.RPCA", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "integrated.RPCA", dims = 1:30,
               reduction.name = "umap.RPCA", reduction.key = "umapRPCA_")
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

## ---- 2. Add ONLY the old EGR1 matrix MA0162.1, then chromVAR -------
egr1_v1 <- getMatrixByID(JASPAR2020, ID = "MA0162.1")   # full versioned ID pins the version
egr1_id <- ID(egr1_v1)                                  # "MA0162.1"
cat("Using motif:", egr1_id, "-", name(egr1_v1), "\n")
pfm <- do.call(TFBSTools::PFMatrixList, list(egr1_v1))

DefaultAssay(obj) <- "peaks"
obj <- AddMotifs(obj, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
stopifnot(egr1_id %in% colnames(GetMotifData(obj, assay = "peaks", slot = "data")))
obj <- RunChromVAR(obj, genome = BSgenome.Mmusculus.UCSC.mm10)   # assay "chromvar"

## ---- 3. FeaturePlot of EGR1(MA0162.1) motif activity on RPCA UMAP --
DefaultAssay(obj) <- "chromvar"
feat <- gsub("_", "-", egr1_id)                          # Seurat sanitises "_" in feature names
if (!feat %in% rownames(obj)) feat <- rownames(obj)[1]
p <- FeaturePlot(obj, features = feat, reduction = "umap.RPCA",
                 split.by = "sex", order = TRUE, pt.size = 0.4) &
     scale_colour_gradientn(colours = c("#3B4CC0", "grey90", "#B40426")) &
     theme(legend.position = "right")
p <- p & patchwork::plot_annotation(title = "EGR1 motif activity (chromVAR, MA0162.1)")
ggsave(file.path(out, "SupplFig1e_chromVAR_Egr1motif_MA0162.1_RPCA.pdf"), p, width = 10, height = 5)

## ---- 4. small git-friendly table ----------------------------------
act <- GetAssayData(obj, assay = "chromvar", slot = "data")[feat, ]
emb <- Embeddings(obj, "umap.RPCA")
write.csv(data.frame(cell = colnames(obj), sex = obj$sex,
                     umapRPCA_1 = emb[, 1], umapRPCA_2 = emb[, 2],
                     EGR1_motif_activity = as.numeric(act)),
          file.path(out, "SupplFig1e_chromVAR_Egr1motif_MA0162.1_activity.csv"), row.names = FALSE)

cat("Done (MA0162.1). plot+csv ->", out, "\n")
