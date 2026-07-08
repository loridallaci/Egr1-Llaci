# ============================================================================
# Cortex development (Roussos) — add JASPAR motifs -> *_withMotifs_V4 objects
# ============================================================================
# STEP 1 (between cortex_dev_build_objects.R and cortex_dev_RENIN_pipeline.R).
# Takes each per-stage `<Stage>_aggregated_object_..._RNAandPeakAssaysAvailable.rds`,
# assigns sex from the barcode suffix, adds motif information to the `peaks`
# assay with `AddMotifs` (JASPAR2020 CORE vertebrates), converts the RNA assay
# to Seurat v4 (`Assay`) for RENIN, and saves `<Stage>_..._withMotifs_V4_122825.rds`
# — the objects the RENIN pipeline and Supp Fig 2 a/b scripts consume.
#
# Consolidated from Roussos_cortex_development_RENIN.Rmd (one block per stage)
# and the "STAGE 0" pseudocode in cortex_dev_RENIN_pipeline.R.
#
# Compute-heavy over full peak sets — run on the HPC via srun/sbatch.
# ============================================================================

library(Seurat)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

# ============================================================================
# SET PATHS - UPDATE THESE FOR YOUR SYSTEM
# ============================================================================

# Directory holding the STEP 0 objects (output of cortex_dev_build_objects.R).
# On HTCF this is /home/lllaci/data/cortex_development.
data_dir <- "C:/Users/loril/Documents/data/roussos_cortex_development/aggregated_objects"

# Stages to process. NOTE the sex mapping flips per stage in the original .Rmd:
#   Adult / Infant / EaFet: -1 = female, -2 = male
#   Adol / Child / LaFet:   -1 = male,   -2 = female
stage_prefix <- c(EaFet = "EaFet", LaFet = "LaFet", Inf = "Inf",
                  Child = "Child", Adol = "Adol", Adult = "Adult")
female_is_suffix1 <- c(EaFet = TRUE, LaFet = FALSE, Inf = TRUE,
                       Child = FALSE, Adol = FALSE, Adult = TRUE)

in_name  <- function(p) paste0(p, "_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds")
out_name <- function(p) paste0(p, "_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds")

# ============================================================================
# JASPAR2020 CORE vertebrate motifs (shared across stages)
# ============================================================================

cat("Fetching JASPAR2020 CORE vertebrate PFMs...\n")
pfm <- getMatrixSet(x = JASPAR2020,
                    opts = list(collection = "CORE", tax_group = "vertebrates",
                                all_versions = FALSE))

# ============================================================================
# PER-STAGE: add motifs -> convert RNA to v4 -> save
# ============================================================================

for (st in names(stage_prefix)) {
  p  <- stage_prefix[[st]]
  fp <- file.path(data_dir, in_name(p))
  if (!file.exists(fp)) { cat("SKIP (missing):", st, "\n"); next }
  cat("\n=====", st, "=====\n")
  obj <- readRDS(fp)

  # sex from barcode suffix (stage-specific mapping)
  s1 <- if (female_is_suffix1[[st]]) "female" else "male"
  s2 <- if (female_is_suffix1[[st]]) "male"   else "female"
  obj$sex <- ifelse(grepl("-1$", colnames(obj)), s1,
             ifelse(grepl("-2$", colnames(obj)), s2, NA))
  print(table(obj$sex, useNA = "ifany"))

  # add motif information to the peaks assay
  DefaultAssay(obj) <- "peaks"
  obj <- AddMotifs(object = obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

  # convert RNA assay to Seurat v4 (Assay) as RENIN expects
  obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay")

  out <- file.path(data_dir, out_name(p))
  saveRDS(obj, out)
  cat("saved ->", out, "\n")
}

cat("\n=== Done. *_withMotifs_V4 objects in:", data_dir, "===\n")
cat("Next: cortex_dev_RENIN_pipeline.R (motif CSVs) and the Supp Fig 2 scripts.\n")

sessionInfo()
