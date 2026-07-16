## =====================================================================
## Supplementary Fig 1e - chromVAR EGR1-motif activity on the RPCA UMAP,
## split by sex.  [MA0162.1 variant]
##
## Identical to Figures_SupplFig1e_chromVAR_Egr1_motif_fix.R EXCEPT it uses
## the OLD EGR1 matrix MA0162.1 instead of the latest MA0162.4. The default
## pipeline (getMatrixSet(..., all_versions = FALSE)) only ever returns the
## newest version, so MA0162.1 is not in the object; here we fetch it
## explicitly by full ID and add only that single motif.
##
## Purpose: sanity-check that the EGR1 motif-activity readout is robust to
## the JASPAR matrix version (MA0162 is EGR1 in both; the two versions differ
## only in the underlying training data, and EGR1's GC-rich site is conserved).
##
## Outputs (kept separate from the MA0162.4 run so the two can be compared):
##   SupplFig1e_chromVAR_Egr1motif_MA0162.1_RPCA.pdf  + activity csv
## =====================================================================
.libPaths(c("C:/Users/loril/AppData/Local/Temp/rlib_matrix_old", .libPaths()))
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

## ---- 2. Add ONLY the old EGR1 matrix MA0162.1 ----------------------
suppressMessages({ library(JASPAR2020); library(TFBSTools); library(BSgenome.Mmusculus.UCSC.mm10); library(chromVAR) })
setMethod("rowSums", "CsparseMatrix", function(x, na.rm = FALSE, dims = 1, ...) {
  if (!methods::is(x, "dMatrix")) x <- x * 1; as.numeric(x %*% rep(1, ncol(x))) })
setMethod("colSums", "CsparseMatrix", function(x, na.rm = FALSE, dims = 1, ...) {
  if (!methods::is(x, "dMatrix")) x <- x * 1; as.numeric(crossprod(rep(1, nrow(x)), x)) })
suppressMessages(library(BiocParallel)); register(SerialParam(), default = TRUE)
cat("BiocParallel default backend:", class(bpparam())[1], "\n")

## fetch EGR1 matrix version 1 explicitly (the full versioned ID pins the version)
egr1_v1 <- getMatrixByID(JASPAR2020, ID = "MA0162.1")
egr1_id <- ID(egr1_v1)                                  # "MA0162.1"
cat("Using motif:", egr1_id, "-", name(egr1_v1), "\n")
pfm <- do.call(TFBSTools::PFMatrixList, list(egr1_v1))

DefaultAssay(obj) <- "peaks"
obj <- AddMotifs(obj, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
stopifnot(egr1_id %in% colnames(GetMotifData(obj, assay = "peaks", slot = "data")))

## ---- 3. chromVAR for the SINGLE MA0162.1 motif, forced SERIAL ------
## coerce counts + annotation to explicit numeric dgCMatrix (no logical -> crossprod)
counts <- GetAssayData(obj, assay = "peaks", slot = "counts")            # peaks x cells
counts <- as(counts, "CsparseMatrix"); counts <- as(counts, "dMatrix")
if (!is.double(counts@x)) counts@x <- as.double(counts@x)

mm_full <- GetMotifData(obj, assay = "peaks", slot = "data")            # peaks x motifs (logical)
mm_egr1 <- mm_full[, egr1_id, drop = FALSE] * 1                         # peaks x 1, numeric
mm_egr1 <- as(as(mm_egr1, "CsparseMatrix"), "dMatrix")
if (!is.double(mm_egr1@x)) mm_egr1@x <- as.double(mm_egr1@x)

cat(sprintf("counts:   class=%s typeof(@x)=%s dim=%dx%d\n",
            class(counts)[1], typeof(counts@x), nrow(counts), ncol(counts)))
cat(sprintf("mm_egr1:  class=%s typeof(@x)=%s dim=%dx%d nnz=%d\n",
            class(mm_egr1)[1], typeof(mm_egr1@x), nrow(mm_egr1), ncol(mm_egr1), length(mm_egr1@x)))

ranges <- granges(obj[["peaks"]])
rse <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts), rowRanges = ranges)
rse <- chromVAR::addGCBias(rse, genome = BSgenome.Mmusculus.UCSC.mm10)
bg  <- chromVAR::getBackgroundPeaks(rse)
dev <- tryCatch(
  chromVAR::computeDeviations(object = rse, annotations = mm_egr1, background_peaks = bg),
  error = function(e) { cat("computeDeviations ERROR:", conditionMessage(e), "\n"); stop(e) })
zscores <- chromVAR::deviationScores(dev)                                # 1 x cells
if (is.null(rownames(zscores))) rownames(zscores) <- egr1_id
obj[["chromvar"]] <- CreateAssayObject(data = zscores)

## ---- 4. FeaturePlot of EGR1(MA0162.1) motif activity on RPCA UMAP --
DefaultAssay(obj) <- "chromvar"
feat <- rownames(zscores)[1]
p <- FeaturePlot(obj, features = feat, reduction = "umap.RPCA",
                 split.by = "sex", order = TRUE, pt.size = 0.4) &
     scale_colour_gradientn(colours = c("#3B4CC0", "grey90", "#B40426")) &
     theme(legend.position = "right")
p <- p & patchwork::plot_annotation(title = "EGR1 motif activity (chromVAR, MA0162.1)")
ggsave(file.path(gitout, "SupplFig1e_chromVAR_Egr1motif_MA0162.1_RPCA.pdf"), p, width = 10, height = 5)

## ---- 5. save small git-friendly table -----------------------------
act <- GetAssayData(obj, assay = "chromvar", slot = "data")[feat, ]
emb <- Embeddings(obj, "umap.RPCA")
tab <- data.frame(cell = colnames(obj), sex = obj$sex,
                  umapRPCA_1 = emb[, 1], umapRPCA_2 = emb[, 2],
                  EGR1_motif_activity = as.numeric(act))
write.csv(tab, file.path(gitout, "SupplFig1e_chromVAR_Egr1motif_MA0162.1_activity.csv"), row.names = FALSE)

cat("Done (MA0162.1). plot+csv ->", gitout, "\n")
