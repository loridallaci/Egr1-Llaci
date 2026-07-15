# ============================================================================
# Re-call MACS2 peaks with the CORRECT mouse effective genome size (mm10),
# then diff the new peak set against the original (human-default) peak set.
#
# Original run used Signac's default effective.genome.size = 2.7e9 (human 'hs').
# Mouse mm10 should be 1.87e9 (MACS2 'mm' preset).
#
# This script writes the new peaks to a NEW file (does NOT overwrite the
# original) and prints a comparison so you can judge how much actually changes.
#
# Run on the cluster with: sbatch macs2_aggr_lot6_mm10gsize_check.sh
# ============================================================================

library(Seurat)
library(Signac)
library(GenomicRanges)

# --- Paths (identical to macs2_aggr_lot6_012925.R) ---------------------------
rds_in   <- "/scratch/rmlab/rmlab_shared3/l.llaci/output/multiome_new_macs2_012925/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925.rds"
fragpath <- "/scratch/rmlab/rmlab_shared3/l.llaci/output/multiome_new_macs2_012925/data/atac_fragments.tsv.gz"
macs2    <- "/ref/rmlab/software/spack-0.21.1/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/py-macs2-2.2.8-k7i7lcl5lqoabuaqsie22uzuzov7hzv2/bin/macs2"

peaks_old_rds <- "/scratch/rmlab/rmlab_shared3/l.llaci/output/multiome_new_macs2_012925/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_macs2Peaks.rds"
peaks_new_rds <- "/scratch/rmlab/rmlab_shared3/l.llaci/output/multiome_new_macs2_012925/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_macs2Peaks_mm10gsize.rds"

# --- Load object & set fragment path -----------------------------------------
aggr_filtered <- readRDS(rds_in)
aggr_filtered@assays$ATAC@fragments[[1]]@path <- fragpath

# --- Re-call peaks with mouse genome size ------------------------------------
# ONLY change vs the original: effective.genome.size = 1.87e9
peaks_new <- CallPeaks(
  aggr_filtered,
  macs2.path            = macs2,
  effective.genome.size = 1.87e9   # mm10; was 2.7e9 (human default) originally
)
saveRDS(peaks_new, peaks_new_rds)

# ============================================================================
# DIFF: new (mm10, 1.87e9) vs original (human, 2.7e9)
# ============================================================================
peaks_old <- readRDS(peaks_old_rds)

# Standardize chromosomes so the comparison is apples-to-apples
peaks_old <- keepStandardChromosomes(peaks_old, pruning.mode = "coarse")
peaks_new <- keepStandardChromosomes(peaks_new, pruning.mode = "coarse")

n_old <- length(peaks_old)
n_new <- length(peaks_new)

bp_old <- sum(width(reduce(peaks_old)))
bp_new <- sum(width(reduce(peaks_new)))

# Reciprocal overlap counts
hits_new_in_old <- sum(countOverlaps(peaks_new, peaks_old) > 0)
hits_old_in_new <- sum(countOverlaps(peaks_old, peaks_new) > 0)

new_only <- n_new - hits_new_in_old   # peaks gained
old_only <- n_old - hits_old_in_new   # peaks lost

cat("\n================ MACS2 genome-size diff ================\n")
cat(sprintf("Original (2.7e9, human) : %d peaks, %s bp covered\n", n_old, format(bp_old, big.mark=",")))
cat(sprintf("New      (1.87e9, mm10) : %d peaks, %s bp covered\n", n_new, format(bp_new, big.mark=",")))
cat(sprintf("Net change in peak count: %+d (%.2f%%)\n", n_new - n_old, 100*(n_new - n_old)/n_old))
cat("--------------------------------------------------------\n")
cat(sprintf("New peaks overlapping an original peak: %d (%.2f%% of new)\n", hits_new_in_old, 100*hits_new_in_old/n_new))
cat(sprintf("Original peaks overlapping a new peak : %d (%.2f%% of original)\n", hits_old_in_new, 100*hits_old_in_new/n_old))
cat(sprintf("Peaks present ONLY in new (gained): %d\n", new_only))
cat(sprintf("Peaks present ONLY in old (lost)  : %d\n", old_only))
cat("========================================================\n\n")

# Save a small summary table next to the peaks
summary_df <- data.frame(
  metric = c("n_peaks_old_2.7e9","n_peaks_new_1.87e9","bp_old","bp_new",
             "new_overlapping_old","old_overlapping_new","gained_only","lost_only"),
  value  = c(n_old, n_new, bp_old, bp_new,
             hits_new_in_old, hits_old_in_new, new_only, old_only)
)
write.csv(summary_df,
          "/scratch/rmlab/rmlab_shared3/l.llaci/output/multiome_new_macs2_012925/data/macs2_gsize_diff_summary.csv",
          row.names = FALSE)

sessionInfo()
