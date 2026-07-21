## =====================================================================
## Supplementary Fig. 2 — cortex development: (1) cell numbers by stage x sex,
## (2) RNA + ATAC UMAPs per stage, coloured by sex.
##
## Sex is assigned from the barcode suffix (-1 = female, -2 = male), the same
## convention used by cortex_dev_RENIN_pipeline.R. Motifs are irrelevant to
## these panels, so any of the aggregated stage objects works.
##
## Canonical run = HTCF, over the 5 withMotifs_V4 objects (set DATA_DIR).
## Local preview = the 3 pre-motif objects in Downloads (Child/Infant/LateFetal):
##   set CORTEX_LOCAL_PREVIEW=1 (Adol/Adult are HTCF-only, so preview is 3/5).
##
## Outputs (06_cortex_development/output/):
##   CellCounts_AllCells_byStageSex.csv
##   SupplFig2_cortex_cellnumbers.pdf
##   SupplFig2_cortex_UMAPs.pdf
## =====================================================================
suppressMessages({ library(Seurat); library(Signac); library(ggplot2); library(patchwork) })
options(future.globals.maxSize = 400 * 1024^3)  # rLSI anchor step ships the object list to a worker;
                                                # the withMotifs_V4 objects (JASPAR motif matrix in the
                                                # peaks assay) exceed 12 GiB, so this must be generous.
                                                # Run on a high-memory host (e.g. .181, 503 GiB).

stage_levels <- c("LateFetal","Infant","Child","Adol","Adult")   # EaFet dropped: both donors female (XIST+/Y-), no M-vs-F contrast
sex_cols <- c(male = "#4A6FE3", female = "#F39AC9")   # male plotted first
out <- Sys.getenv("OUT_DIR", unset = "C:/Users/loril/Documents/GitHub/Egr1-Llaci/06_cortex_development/output")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

if (Sys.getenv("CORTEX_LOCAL_PREVIEW") == "1") {
  data_dir <- "C:/Users/loril/Downloads"
  stage_files <- c(
    LateFetal = "LaFet_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds",
    Infant    = "Inf_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds",
    Child     = "Child_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds")
} else {
  data_dir <- Sys.getenv("DATA_DIR", unset = "/home/lllaci/data/cortex_development")
  stage_files <- c(
    LateFetal = "LaFet_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
    Infant    = "Inf_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
    Child     = "Child_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
    Adol      = "Adol_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
    Adult     = "Adult_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds")
}

big_theme <- theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), legend.title = element_text(size = 14),
        legend.text = element_text(size = 14), strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank())

counts <- list(); rna_umaps <- list(); atac_umaps <- list()
for (st in names(stage_files)) {
  fp <- file.path(data_dir, stage_files[[st]])
  if (!file.exists(fp)) { cat("SKIP (missing):", st, "\n"); next }
  cat("stage:", st, "\n")
  o <- readRDS(fp)
  ## drop Fragment objects: they point to ATAC fragment files on HTCF that aren't
  ## present locally, and any subset/split would try to validate them (GetIndexFile
  ## error). TFIDF/SVD/UMAP use the peak counts matrix, so fragments aren't needed.
  for (a in intersect(c("peaks","ATAC"), Assays(o)))
    suppressWarnings(try(Fragments(o[[a]]) <- NULL, silent = TRUE))
  ## Sex: use the per-stage assignment STORED in the object by cortex_dev_add_motifs.R.
  ## Do NOT re-derive from the barcode suffix — the suffix->sex map FLIPS per stage
  ## (LaFet/Child/Adol have -1=male, -2=female; Inf/Adult/EaFet have -1=female), so a
  ## uniform -1=female map silently swaps male/female for LaFet, Child and Adol.
  if (!"sex" %in% colnames(o[[]])) stop("stored $sex missing for stage ", st)
  o$sex <- as.character(o$sex)
  o <- o[, o$sex %in% c("female","male")]
  o$sex <- factor(o$sex, levels = c("male","female"))
  counts[[st]] <- as.data.frame(table(sex = o$sex)); counts[[st]]$stage <- st

  ## ---- RNA UMAP: RPCA-integrated across sex (matches Supp Fig 1b) ----
  DefaultAssay(o) <- "RNA"
  if (!inherits(o[["RNA"]], "Assay5")) o[["RNA"]] <- as(o[["RNA"]], "Assay5")
  o[["RNA"]] <- split(o[["RNA"]], f = o$sex)
  o <- NormalizeData(o, verbose = FALSE); o <- FindVariableFeatures(o, verbose = FALSE)
  o <- ScaleData(o, verbose = FALSE); o <- RunPCA(o, verbose = FALSE)
  o <- IntegrateLayers(o, method = RPCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.RPCA", verbose = FALSE)
  o <- RunUMAP(o, reduction = "integrated.RPCA", dims = 1:30, reduction.name = "umap.rna", verbose = FALSE)
  o[["RNA"]] <- JoinLayers(o[["RNA"]])
  e <- Embeddings(o, "umap.rna")
  rna_df <- data.frame(x = e[,1], y = e[,2],
                       sex = factor(o$sex[rownames(e)], levels = c("male","female")))

  ## ---- ATAC UMAP: rLSI-integrated across sex (matches Supp Fig 1c) ---
  DefaultAssay(o) <- "peaks"
  o <- RunTFIDF(o, verbose = FALSE); o <- FindTopFeatures(o, min.cutoff = 5); o <- RunSVD(o, verbose = FALSE)
  ol <- SplitObject(o, split.by = "sex")
  ol <- lapply(ol, function(x){ x <- RunTFIDF(x, verbose=FALSE); x <- FindTopFeatures(x, min.cutoff=5); x <- RunSVD(x, verbose=FALSE); x })
  anch  <- FindIntegrationAnchors(object.list = ol, anchor.features = rownames(o), reduction = "rlsi", dims = 2:30, verbose = FALSE)
  integ <- IntegrateEmbeddings(anchorset = anch, reductions = o[["lsi"]],
                               new.reduction.name = "integrated_lsi", dims.to.integrate = 1:30, verbose = FALSE)
  integ <- RunUMAP(integ, reduction = "integrated_lsi", dims = 2:30, reduction.name = "umap.atac", verbose = FALSE)
  ea <- Embeddings(integ, "umap.atac")
  atac_df <- data.frame(x = ea[,1], y = ea[,2],
                        sex = factor(o$sex[rownames(ea)], levels = c("male","female")))

  mk <- function(df, ttl) ggplot(df[sample(nrow(df)), ], aes(x, y, colour = sex)) +
      geom_point(size = 0.3, alpha = 0.7) +
      scale_colour_manual(values = sex_cols, guide = "none") + coord_fixed() +
      labs(title = ttl, x = NULL, y = NULL) + big_theme +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
  rna_umaps[[st]]  <- mk(rna_df,  st)
  atac_umaps[[st]] <- mk(atac_df, st)
  rm(o); gc()
}

## ---- cell-number bar (stage x sex) --------------------------------
cdf <- do.call(rbind, counts)
cdf$stage <- factor(cdf$stage, levels = stage_levels)
write.csv(cdf[, c("stage","sex","Freq")], file.path(out, "CellCounts_AllCells_byStageSex.csv"), row.names = FALSE)
pc <- ggplot(cdf, aes(stage, Freq, fill = sex)) +
  geom_col(position = position_dodge(0.75), width = 0.7, colour = "black") +
  scale_fill_manual(values = sex_cols, name = "Sex") +
  labs(title = "Cortex development: cell numbers by stage and sex",
       x = "developmental stage", y = "number of nuclei") +
  big_theme + theme(legend.position = "top", axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(out, "SupplFig2_cortex_cellnumbers.pdf"), pc, width = 8, height = 5.5)

## ---- UMAP montage (RNA row / ATAC row, columns = stage) -----------
ord <- intersect(stage_levels, names(rna_umaps))
row_rna  <- wrap_plots(rna_umaps[ord],  nrow = 1)
row_atac <- wrap_plots(atac_umaps[ord], nrow = 1)
mont <- (row_rna / row_atac) +
  plot_annotation(title = "Cortex development — RNA (RPCA, top) and ATAC (rLSI, bottom) UMAPs, integrated across sex",
                  theme = theme(plot.title = element_text(size = 18, face = "bold")))
ggsave(file.path(out, "SupplFig2_cortex_UMAPs.pdf"), mont, width = 3.2*length(ord), height = 6.6)

cat("\nStages built:", paste(ord, collapse=", "), "\n")
print(cdf[, c("stage","sex","Freq")]); cat("Done ->", out, "\n")
