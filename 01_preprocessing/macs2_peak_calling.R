# ============================================================================
# MACS2 Peak Calling for Multiome ATAC Data
# ============================================================================
# This script calls peaks using MACS2 and adds them to the Seurat object
# Run this on a system with MACS2 installed
# ============================================================================

library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomicRanges)

# ============================================================================
# SET PATHS
# ============================================================================

# Path to filtered Seurat object
# original (author's machine): "/path/to/seurat_object_filtered.rds"
seurat_file <- "output/seurat_object_filtered.rds"

# Path to MACS2 executable
# original (author's machine): "/path/to/macs2"
macs2_path <- "macs2"  # e.g., "/usr/local/bin/macs2" or "~/anaconda3/bin/macs2"

# Output directory
# original (author's machine): "/path/to/output"
output_dir <- "output"

# Blacklist regions (download from ENCODE)
# original (author's machine): "/path/to/mm10-blacklist.v2.bed"
blacklist_file <- "data/mm10-blacklist.v2.bed"  # For mouse
# blacklist_file <- "data/hg38-blacklist.v2.bed"  # For human  # original: "/path/to/hg38-blacklist.v2.bed"

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading filtered Seurat object...\n")
seurat_obj <- readRDS(seurat_file)

# Get gene annotations
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"

# ============================================================================
# CALL PEAKS WITH MACS2
# ============================================================================

cat("Calling peaks with MACS2...\n")

DefaultAssay(seurat_obj) <- "ATAC"

peaks <- CallPeaks(
  seurat_obj,
  macs2.path = macs2_path
)

cat("Total peaks called:", length(peaks), "\n")

# ============================================================================
# FILTER PEAKS
# ============================================================================

cat("Filtering peaks...\n")

# Load blacklist regions
blacklist <- rtracklayer::import(blacklist_file)

# Remove peaks on nonstandard chromosomes
peaks_filtered <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

# Remove peaks in blacklist regions
peaks_filtered <- subsetByOverlaps(
  x = peaks_filtered,
  ranges = blacklist,
  invert = TRUE
)

cat("Peaks after filtering:", length(peaks_filtered), "\n")

# Save peaks
saveRDS(peaks_filtered, file.path(output_dir, "macs2_peaks_filtered.rds"))

# ============================================================================
# QUANTIFY PEAKS
# ============================================================================

cat("Quantifying peaks...\n")

# Quantify counts in each peak
peak_counts <- FeatureMatrix(
  fragments = Fragments(seurat_obj),
  features = peaks_filtered,
  cells = colnames(seurat_obj)
)

# ============================================================================
# ADD PEAKS TO SEURAT OBJECT
# ============================================================================

cat("Adding peaks to Seurat object...\n")

# Get fragment file path
fragpath <- seurat_obj@assays$ATAC@fragments[[1]]@path

# Create a new assay using the MACS2 peak set
seurat_obj[["peaks"]] <- CreateChromatinAssay(
  counts = peak_counts,
  fragments = fragpath,
  annotation = annotation
)

cat("Peaks added successfully\n")

# ============================================================================
# SAVE UPDATED OBJECT
# ============================================================================

cat("Saving Seurat object with peaks...\n")

saveRDS(
  seurat_obj,
  file.path(output_dir, "seurat_object_filtered_withPeaks.rds")
)

cat("\n=== Peak Calling Complete ===\n")
cat("Output saved to:", output_dir, "\n\n")

sessionInfo()

