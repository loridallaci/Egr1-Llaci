# ============================================================================
# Add JASPAR motifs to the PAPER-CELL objects -> *_withMotifs_V4 (100% cells)
# ============================================================================
# Variant of cortex_dev_add_motifs.R for the paper-cell rebuild
# (cortex_dev_build_paperCells.R output). KEY DIFFERENCE: sex is already set
# authoritatively from paper_cell_annotation.csv (0 NA), so we do NOT overwrite
# it from the barcode suffix — that suffix reassignment is what caused the
# earlier LaFet/Child/Adol sex-swap. Processes the 4 rebuilt stages; Adol is
# added later once its cellranger-arc run finishes.
#
# Run on .181 (RENIN env):
#   LD_LIBRARY_PATH=$HOME/local/lib /home/lllaci/R-4.2.2/bin/Rscript cortex_dev_add_motifs_paperCells.R
# ============================================================================
suppressMessages({library(Seurat); library(Signac); library(JASPAR2020); library(TFBSTools); library(BSgenome.Hsapiens.UCSC.hg38)})
set.seed(1234)

data_dir <- "/home/lllaci/data/cortex_development_paperCells"
prefix   <- c("LaFet", "Inf", "Child", "Adult")          # Adol appended later
in_name  <- function(p) paste0(p, "_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds")
out_name <- function(p) paste0(p, "_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds")

cat("Fetching JASPAR2020 CORE vertebrate PFMs...\n")
pfm <- getMatrixSet(x = JASPAR2020,
                    opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))

for (p in prefix) {
  fp <- file.path(data_dir, in_name(p))
  if (!file.exists(fp)) { cat("SKIP (missing):", p, "\n"); next }
  cat("\n=====", p, "=====\n")
  obj <- readRDS(fp)
  cat("sex (from paper annotation, kept):\n"); print(table(obj$sex, useNA = "ifany"))

  DefaultAssay(obj) <- "peaks"
  obj <- AddMotifs(object = obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
  obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay")   # v4 Assay for RENIN

  out <- file.path(data_dir, out_name(p))
  saveRDS(obj, out); cat("saved ->", out, "\n")
}
cat("\n=== Done. paper-cell *_withMotifs_V4 objects in:", data_dir, "===\n")
