## =====================================================================
## Supplementary Fig 1e - chromVAR EGR1-motif activity on the RPCA UMAP,
## split by sex.
##
## Builds the chromVAR object the same way as 01_preprocessing/add_motifs_chromvar.R
## (JASPAR2020 CORE vertebrate motifs -> AddMotifs -> RunChromVAR), computes the
## RPCA-integrated RNA UMAP, finds the EGR1 motif, and FeaturePlots its per-cell
## motif-activity z-score split by sex.
##
## Outputs:
##   git : SupplFig1e_chromVAR_Egr1motif_RPCA.pdf  + Egr1 activity table (csv)
##   local: full chromVAR object (too large for git) in the data folder
## =====================================================================
## Matrix >= 1.6.4 + MatrixGenerics/chromVAR break sparse rowSums/colSums
## (`object 'CRsparse_rowSums' not found`). SeuratObject needs Matrix >= 1.6.4,
## so we can't just downgrade; instead use Matrix 1.6.4 (temp lib) and patch the
## broken methods with matrix-multiply versions just before RunChromVAR (below).
.libPaths(c("C:/Users/loril/AppData/Local/Temp/rlib_matrix_old", .libPaths()))

## NOTE: load only Seurat/Signac before the v5 layer split() below.
## Loading chromVAR/Bioconductor matrix packages first breaks the v3->v5
## assay conversion (object 'CRsparse_colSums' not found); load them AFTER.
suppressMessages({ library(Seurat); library(Signac); library(ggplot2); library(dplyr) })

data_dir <- "C:/Users/loril/Documents/data/multiome_081721_aggregated_data_analysis"
rds   <- file.path(data_dir, "female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks.rds")
gitout<- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/03_atac_analysis/output/figures/suppl_fig1"
dir.create(gitout, showWarnings = FALSE, recursive = TRUE)
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

## ---- 2. Add motifs + chromVAR (same as add_motifs_chromvar.R) ------
## load the chromVAR/Bioconductor stack only now (after the v5 split above)
suppressMessages({ library(JASPAR2020); library(TFBSTools); library(BSgenome.Mmusculus.UCSC.mm10); library(chromVAR) })
## patch the broken sparse rowSums/colSums (see header note) with matrix-multiply
## versions; coerce logical sparse (motif matches) -> numeric so crossprod works.
setMethod("rowSums", "CsparseMatrix", function(x, na.rm = FALSE, dims = 1, ...) {
  if (!methods::is(x, "dMatrix")) x <- x * 1; as.numeric(x %*% rep(1, ncol(x))) })
setMethod("colSums", "CsparseMatrix", function(x, na.rm = FALSE, dims = 1, ...) {
  if (!methods::is(x, "dMatrix")) x <- x * 1; as.numeric(crossprod(rep(1, nrow(x)), x)) })
## force SERIAL chromVAR: parallel workers OOM AND don't inherit the patch above
suppressMessages(library(BiocParallel)); register(SerialParam(), default = TRUE)
cat("BiocParallel default backend:", class(bpparam())[1], "\n")
pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
DefaultAssay(obj) <- "peaks"
obj <- AddMotifs(obj, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

## ---- chromVAR by hand, forced SERIAL, in the MAIN process ----------
## Signac::RunChromVAR spawns parallel workers that re-load the broken main-lib
## Matrix and ignore our patches -> crossprod/rowSums failures. Calling chromVAR
## directly with BPPARAM = SerialParam() keeps it all in this process, where the
## temp-lib Matrix 1.6.4 + the rowSums/colSums patches apply.
counts <- GetAssayData(obj, assay = "peaks", slot = "counts")            # peaks x cells
ranges <- granges(obj[["peaks"]])
mm     <- as(GetMotifData(obj, assay = "peaks", slot = "data") * 1, "CsparseMatrix")  # peaks x motifs (numeric)
rse <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts), rowRanges = ranges)
rse <- chromVAR::addGCBias(rse, genome = BSgenome.Mmusculus.UCSC.mm10)
bg  <- chromVAR::getBackgroundPeaks(rse)   # not parallel; no BPPARAM arg
dev <- chromVAR::computeDeviations(object = rse, annotations = mm,
                                   background_peaks = bg)   # uses registered SerialParam
zscores <- chromVAR::deviationScores(dev)                                # motifs x cells
obj[["chromvar"]] <- CreateAssayObject(data = zscores)

## ---- 3. find the EGR1 motif id ------------------------------------
mn <- Motifs(obj[["peaks"]])@motif.names           # named list: id -> TF name
egr1_id <- names(mn)[toupper(unlist(mn)) == "EGR1"][1]
cat("EGR1 motif id:", egr1_id, "\n")

## ---- 4. FeaturePlot of EGR1 motif activity on RPCA UMAP, by sex ----
DefaultAssay(obj) <- "chromvar"
p <- FeaturePlot(obj, features = egr1_id, reduction = "umap.RPCA",
                 split.by = "sex", order = TRUE, pt.size = 0.4) &
     scale_colour_gradientn(colours = c("#3B4CC0", "grey90", "#B40426")) &
     theme(legend.position = "right")
p <- p & patchwork::plot_annotation(title = "EGR1 motif activity (chromVAR)")
ggsave(file.path(gitout, "SupplFig1e_chromVAR_Egr1motif_RPCA.pdf"), p, width = 10, height = 5)

## ---- 5. save outputs ----------------------------------------------
## small, git-friendly: per-cell EGR1 activity + UMAP coords + sex
act <- GetAssayData(obj, assay = "chromvar", slot = "data")[egr1_id, ]
emb <- Embeddings(obj, "umap.RPCA")
tab <- data.frame(cell = colnames(obj), sex = obj$sex,
                  umapRPCA_1 = emb[, 1], umapRPCA_2 = emb[, 2],
                  EGR1_motif_activity = as.numeric(act))
write.csv(tab, file.path(gitout, "SupplFig1e_chromVAR_Egr1motif_activity.csv"), row.names = FALSE)

## full object (too big for git) -> local data folder, for reuse
saveRDS(obj, file.path(data_dir, "female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_RPCA_local.rds"))

cat("Done.\n  plot+csv -> ", gitout, "\n  full object -> ", data_dir, "\n")
