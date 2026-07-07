## =====================================================================
## Supplementary Fig 1e - chromVAR EGR1-motif activity on the RPCA UMAP,
## split by sex.  [PATCHED single-motif variant]
##
## Same as Figures_SupplFig1e_chromVAR_Egr1_motif.R, but computes chromVAR
## deviations for ONLY the EGR1 motif annotation (background peaks are
## motif-independent, so the EGR1 z-score is identical to the full run). This
## collapses the 745-motif BiocParallel loop - where the Matrix>=1.6.4 sparse
## crossprod bug fires - to a single serial iteration, and coerces the count /
## annotation matrices to explicit numeric dgCMatrix so crossprod never sees a
## logical. Type diagnostics are printed before the deviation call.
##
## Outputs:
##   git : SupplFig1e_chromVAR_Egr1motif_RPCA.pdf + Egr1 activity table (csv)
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

## ---- 2. Add motifs (same as add_motifs_chromvar.R) -----------------
suppressMessages({ library(JASPAR2020); library(TFBSTools); library(BSgenome.Mmusculus.UCSC.mm10); library(chromVAR) })
setMethod("rowSums", "CsparseMatrix", function(x, na.rm = FALSE, dims = 1, ...) {
  if (!methods::is(x, "dMatrix")) x <- x * 1; as.numeric(x %*% rep(1, ncol(x))) })
setMethod("colSums", "CsparseMatrix", function(x, na.rm = FALSE, dims = 1, ...) {
  if (!methods::is(x, "dMatrix")) x <- x * 1; as.numeric(crossprod(rep(1, nrow(x)), x)) })
suppressMessages(library(BiocParallel)); register(SerialParam(), default = TRUE)
cat("BiocParallel default backend:", class(bpparam())[1], "\n")
pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
DefaultAssay(obj) <- "peaks"
obj <- AddMotifs(obj, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

## ---- 3. find the EGR1 motif id FIRST ------------------------------
mn <- Motifs(obj[["peaks"]])@motif.names           # named list: id -> TF name
egr1_id <- names(mn)[toupper(unlist(mn)) == "EGR1"][1]
stopifnot(!is.na(egr1_id))
cat("EGR1 motif id:", egr1_id, "\n")

## ---- 4. chromVAR for the SINGLE EGR1 motif, forced SERIAL ----------
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

## ---- 5. EGR1 motif activity on the RPCA UMAP, by sex --------------
## Plot from a data frame so BOTH sexes share ONE colour scale and the panel
## matches Fig 1b's RPCA-UMAP-by-sex layout. (Seurat FeaturePlot(split.by=)
## rescales each facet independently, which hid the male-vs-female difference.)
DefaultAssay(obj) <- "chromvar"
feat <- rownames(zscores)[1]
emb  <- Embeddings(obj, "umap.RPCA")
dfp  <- data.frame(umapRPCA_1 = emb[,1], umapRPCA_2 = emb[,2],
                   sex = factor(obj$sex, levels = c("female","male")),
                   activity = as.numeric(GetAssayData(obj, assay="chromvar", slot="data")[feat, ]))
dfp  <- dfp[order(abs(dfp$activity)), ]                          # extreme cells on top
lim  <- as.numeric(ceiling(quantile(abs(dfp$activity), 0.99)))   # symmetric robust cap
p <- ggplot(dfp, aes(umapRPCA_1, umapRPCA_2, colour = activity)) +
  geom_point(size = 0.7) + facet_wrap(~ sex) +
  scale_colour_gradient2(low = "#3B4CC0", mid = "grey92", high = "#B40426", midpoint = 0,
                         limits = c(-lim, lim), oob = scales::squish, name = "EGR1 motif\nactivity (z)") +
  coord_fixed() +
  labs(title = paste0("EGR1 motif activity (chromVAR, ", feat, ")"),
       x = "umapRPCA_1", y = "umapRPCA_2") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size = 14), legend.text = element_text(size = 14))
ggsave(file.path(gitout, "SupplFig1e_chromVAR_Egr1motif_RPCA.pdf"), p, width = 10, height = 5)

## ---- 6. save small git-friendly table -----------------------------
act <- GetAssayData(obj, assay = "chromvar", slot = "data")[feat, ]
emb <- Embeddings(obj, "umap.RPCA")
tab <- data.frame(cell = colnames(obj), sex = obj$sex,
                  umapRPCA_1 = emb[, 1], umapRPCA_2 = emb[, 2],
                  EGR1_motif_activity = as.numeric(act))
write.csv(tab, file.path(gitout, "SupplFig1e_chromVAR_Egr1motif_activity.csv"), row.names = FALSE)

cat("Done. plot+csv ->", gitout, "\n")
