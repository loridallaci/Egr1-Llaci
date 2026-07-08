# ============================================================================
# Cortex development (Roussos) — build the per-stage aggregated Seurat objects
# ============================================================================
# STEP 0 of the cortex-development arm (Supplementary Fig. 2). Rebuilds the
# `<Stage>_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds`
# objects that every downstream cortex script consumes.
#
# External data: Zhu et al., "Multi-omic profiling of the developing human
# cerebral cortex at the single-cell level", Sci Adv 2023;9(41):eadg3754
# (doi:10.1126/sciadv.adg3754). NOT processed from FASTQ — we start from the
# deposited 10X Multiome Cell Ranger ARC outputs (GEO GSE204684; also Broad
# Single Cell Portal SCP1859 / the dual-assay S3 hub). One sample per stage,
# each with a filtered_feature_bc_matrix.h5 and an atac_fragments.tsv.gz.
#
# Per stage this script: builds an RNA + ATAC Seurat/Signac object (hg38),
# QC-filters cells, re-calls peaks with MACS2, quantifies them as a `peaks`
# assay, and saves. Motifs are added later (see cortex_dev_RENIN_pipeline.R,
# "STAGE 0"), producing the *_withMotifs_V4 objects.
#
# Consolidated from Roussos_aggregated_files_112223.Rmd (the original had one
# copy-pasted block per stage). Compute-heavy (MACS2 + FeatureMatrix over full
# fragment files) — run on the HPC via srun/sbatch, not the login node.
# ============================================================================

library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(dplyr)
library(patchwork)

set.seed(1234)

# ============================================================================
# SET PATHS - UPDATE THESE FOR YOUR SYSTEM
# ============================================================================

# Root holding the Roussos Cell Ranger ARC outputs, one sub-directory per stage.
# Original (author's Mac, external drive):
#   /Volumes/DUAL_DRIVE/data_mining/roussos_human_cerebral_cortex_science23
data_root <- "C:/Users/loril/Documents/data/roussos_cortex_development"

# One 10X Multiome sample per developmental stage. Each stage directory must
# contain `filtered_feature_bc_matrix.h5` and `atac_fragments.tsv.gz`
# (+ its .tbi index). Order = Supp Fig. 2 columns (Early Fetal -> Adult).
stages <- c("EaFet", "LaFet", "Inf", "Child", "Adol", "Adult")
stage_dir <- setNames(file.path(data_root, stages), stages)

# Path to MACS2 executable (no `module load` on this HPC — point at your conda/venv).
# original (author's machine): "/path/to/macs2"
macs2_path <- "macs2"  # e.g. "~/miniconda3/envs/r_multiome/bin/macs2"

# ENCODE hg38 blacklist (download once).
# original (author's machine): "/path/to/hg38-blacklist.v2.bed"
blacklist_file <- file.path(data_root, "hg38-blacklist.v2.bed")

# Output directory for the aggregated objects.
output_dir <- file.path(data_root, "aggregated_objects")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# QC thresholds (identical to Roussos_aggregated_files_112223.Rmd)
qc <- list(nCount_RNA_max = 40000, nCount_ATAC_max = 40000,
           nucleosome_signal_min = 0.1, nucleosome_signal_max = 2,
           TSS_enrichment_min = 0.5, TSS_enrichment_max = 9)

# ============================================================================
# SHARED ANNOTATION (hg38)
# ============================================================================

cat("Loading hg38 gene annotation (EnsDb.Hsapiens.v86)...\n")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

blacklist <- rtracklayer::import(blacklist_file)

# ============================================================================
# PER-STAGE BUILD
# ============================================================================

build_stage <- function(stage) {
  cat("\n========================= ", stage, " =========================\n")
  h5   <- file.path(stage_dir[[stage]], "filtered_feature_bc_matrix.h5")
  frag <- file.path(stage_dir[[stage]], "atac_fragments.tsv.gz")
  stopifnot(file.exists(h5), file.exists(frag))

  # --- 1. RNA + ATAC object ---
  counts <- Read10X_h5(h5)
  obj <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA")
  obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks, sep = c(":", "-"),
    fragments = frag, annotation = annotation)

  # --- 2. QC metrics + sex from barcode suffix ---
  DefaultAssay(obj) <- "ATAC"
  obj <- NucleosomeSignal(obj)
  obj <- TSSEnrichment(obj)
  # NOTE: the -1/-2 -> sex mapping is stage-specific in the original .Rmd; verify
  # per stage before use (Adult/Infant/EaFet: -1 = female; Adol/Child/LaFet: -1 = male).
  obj$sex <- ifelse(grepl("-1$", colnames(obj)), "female",
             ifelse(grepl("-2$", colnames(obj)), "male", NA))
  cat("cells before filter:", ncol(obj), "\n")

  # --- 3. filter ---
  obj <- subset(obj, subset =
    nCount_ATAC < qc$nCount_ATAC_max & nCount_RNA < qc$nCount_RNA_max &
    nucleosome_signal > qc$nucleosome_signal_min & nucleosome_signal < qc$nucleosome_signal_max &
    TSS.enrichment > qc$TSS_enrichment_min & TSS.enrichment < qc$TSS_enrichment_max)
  cat("cells after filter: ", ncol(obj), "\n")

  # --- 4. MACS2 peaks -> quantify -> add `peaks` assay ---
  peaks <- CallPeaks(obj, macs2.path = macs2_path)
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(peaks, blacklist, invert = TRUE)
  saveRDS(peaks, file.path(output_dir, paste0(stage, "_filtered_macs2Peaks_112223.rds")))

  macs2_counts <- FeatureMatrix(
    fragments = Fragments(obj), features = peaks, cells = colnames(obj))
  obj[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts, fragments = frag, annotation = annotation)

  # --- 5. save ---
  out <- file.path(output_dir,
    paste0(stage, "_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds"))
  saveRDS(obj, out)
  cat("saved ->", out, "\n")
  invisible(out)
}

for (st in stages) build_stage(st)

cat("\n=== Done. Aggregated per-stage objects in:", output_dir, "===\n")
cat("Next: add JASPAR motifs to make the *_withMotifs_V4 objects",
    "(see cortex_dev_RENIN_pipeline.R, STAGE 0).\n")

sessionInfo()
