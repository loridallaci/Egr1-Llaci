# =============================================================================
# 01_preprocessing: Add Motifs and Run ChromVAR
# =============================================================================
# Description:
#   Adds JASPAR2020 motif information to the peaks assay and computes
#   per-cell transcription factor motif activity scores using ChromVAR.
#   Sex labels are assigned based on Cell Ranger ARC barcode suffixes
#   (-1 = female, -2 = male from aggregation).
#
# Run order:
#   1. qc_filtering.R
#   2. macs2_peak_calling.R
#   3. add_motifs_chromvar.R   <-- this script
#
# Input:
#   female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks.rds
#
# Output:
#   female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds
# =============================================================================

# --- Libraries ----------------------------------------------------------------

library(Seurat)
library(Signac)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
library(patchwork)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# --- Load object --------------------------------------------------------------

object <- readRDS("output/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks.rds")  # original: "/home/lllaci/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks.rds"

# --- Assign sex labels --------------------------------------------------------
# Cell Ranger ARC appends barcode suffixes during aggregation:
# -1 = female sample, -2 = male sample

DefaultAssay(object) <- "RNA"
object$sex <- ifelse(grepl("-1$", colnames(object)), "female",
                     ifelse(grepl("-2$", colnames(object)), "male", NA))
object <- SetIdent(object, value = object@meta.data$sex)

# --- Add motifs ---------------------------------------------------------------
# Retrieve vertebrate CORE motif matrices from JASPAR2020

pfm <- getMatrixSet(
  x    = JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)

DefaultAssay(object) <- "peaks"
object <- AddMotifs(
  object = object,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm    = pfm
)

# --- Run ChromVAR -------------------------------------------------------------
# Computes per-cell motif activity scores stored in the 'chromvar' assay

object <- RunChromVAR(
  object = object,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# --- Save ---------------------------------------------------------------------

saveRDS(object, "output/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds")  # original: "/home/lllaci/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds"

message("Done. Object saved with motifs and ChromVAR assay.")
