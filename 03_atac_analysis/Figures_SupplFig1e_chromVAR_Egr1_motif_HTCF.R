## =====================================================================
## Supplementary Fig 1e - chromVAR EGR1-motif activity on the RPCA UMAP,
## split by sex.  [HTCF version - canonical Signac::RunChromVAR]
##
## Run on HTCF, where the Matrix/chromVAR stack is intact (the local Windows
## install hits a Matrix >= 1.6 crossprod bug inside chromVAR). No matrix hacks
## here - just RunChromVAR.
##
## Paths are read from env vars (set by the sbatch wrapper) so nothing is
## hard-coded:
##   MULTIOME_RDS - the lot6 multiome Seurat object (RNA + MACS2 `peaks` assays)
##   OUTDIR       - where to write the PDF + CSV
##
## Outputs:
##   SupplFig1e_chromVAR_Egr1motif_RPCA.pdf
##   SupplFig1e_chromVAR_Egr1motif_activity.csv
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

## ---- 2. Motifs + chromVAR (canonical) -----------------------------
pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE",
                    tax_group = "vertebrates", all_versions = FALSE))
DefaultAssay(obj) <- "peaks"
obj <- AddMotifs(obj, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
obj <- RunChromVAR(obj, genome = BSgenome.Mmusculus.UCSC.mm10)   # assay "chromvar"

## ---- 3. EGR1 motif id ---------------------------------------------
mn <- Motifs(obj[["peaks"]])@motif.names
egr1_id <- names(mn)[toupper(unlist(mn)) == "EGR1"][1]
stopifnot(!is.na(egr1_id))
cat("EGR1 motif id:", egr1_id, "\n")

## ---- 4. EGR1 motif activity on the RPCA UMAP, by sex --------------
## Plot from a data frame so BOTH sexes share ONE colour scale and the panel
## matches Fig 1b's RPCA-UMAP-by-sex layout. (Seurat FeaturePlot(split.by=)
## rescales each facet independently, which hid the male-vs-female difference.)
DefaultAssay(obj) <- "chromvar"
emb <- Embeddings(obj, "umap.RPCA")
dfp <- data.frame(umapRPCA_1 = emb[,1], umapRPCA_2 = emb[,2],
                  sex = factor(obj$sex, levels = c("female","male")),
                  activity = as.numeric(GetAssayData(obj, assay="chromvar", slot="data")[egr1_id, ]))
dfp <- dfp[order(abs(dfp$activity)), ]                          # extreme cells on top
lim <- as.numeric(ceiling(max(abs(dfp$activity))))   # full range -> graded red, no saturation
p <- ggplot(dfp, aes(umapRPCA_1, umapRPCA_2, colour = activity)) +
  geom_point(size = 0.7) + facet_wrap(~ sex) +
  scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#D6604D", midpoint = 0,
                         limits = c(-lim, lim), name = "EGR1 motif\nactivity (z)") +
  coord_fixed() +
  labs(title = paste0("EGR1 motif activity (chromVAR, ", egr1_id, ")"),
       x = "umapRPCA_1", y = "umapRPCA_2") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size = 14), legend.text = element_text(size = 14))
ggsave(file.path(out, "SupplFig1e_chromVAR_Egr1motif_RPCA.pdf"), p, width = 10, height = 5)

## ---- 5. small git-friendly table ----------------------------------
act <- GetAssayData(obj, assay = "chromvar", slot = "data")[egr1_id, ]
emb <- Embeddings(obj, "umap.RPCA")
write.csv(data.frame(cell = colnames(obj), sex = obj$sex,
                     umapRPCA_1 = emb[, 1], umapRPCA_2 = emb[, 2],
                     EGR1_motif_activity = as.numeric(act)),
          file.path(out, "SupplFig1e_chromVAR_Egr1motif_activity.csv"), row.names = FALSE)

cat("Done. plot+csv ->", out, "\n")
