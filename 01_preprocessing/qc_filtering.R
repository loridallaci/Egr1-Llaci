# ============================================================================
# Quality Control and Filtering for Multiome Data
# ============================================================================
# This script performs QC on aggregated Cell Ranger ARC output and filters
# low-quality cells based on RNA and ATAC metrics
# ============================================================================

library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79)  # Use EnsDb.Hsapiens.v86 for human
library(BSgenome.Mmusculus.UCSC.mm10)  # Use BSgenome.Hsapiens.UCSC.hg38 for human
library(dplyr)
library(patchwork)
library(scales)

set.seed(1234)

# ============================================================================
# SET PATHS - UPDATE THESE FOR YOUR SYSTEM
# ============================================================================

# Input files from Cell Ranger ARC aggregation
data_dir <- "/path/to/cellranger_aggr/outs"
counts_file <- file.path(data_dir, "filtered_feature_bc_matrix.h5")
fragments_file <- file.path(data_dir, "atac_fragments.tsv.gz")

# Output directory
output_dir <- "/path/to/output"
figures_dir <- file.path(output_dir, "figures/qc")

# Create output directories if they don't exist
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading Cell Ranger ARC output...\n")

# Load the RNA and ATAC data
counts <- Read10X_h5(counts_file)

# Get gene annotations (mm10 for mouse, hg38 for human)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"

# ============================================================================
# CREATE SEURAT OBJECT
# ============================================================================

cat("Creating Seurat object...\n")

# Create a Seurat object containing the RNA data
seurat_obj <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# Create ATAC assay and add it to the object
seurat_obj[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragments_file,
  annotation = annotation
)

cat("Cells before filtering:", ncol(seurat_obj), "\n")

# ============================================================================
# CALCULATE QC METRICS
# ============================================================================

cat("Calculating QC metrics...\n")

# Switch to ATAC assay for ATAC-specific metrics
DefaultAssay(seurat_obj) <- "ATAC"

# Calculate nucleosome signal
seurat_obj <- NucleosomeSignal(seurat_obj)

# Calculate TSS enrichment
seurat_obj <- TSSEnrichment(seurat_obj)

# Switch to RNA assay for RNA-specific metrics
DefaultAssay(seurat_obj) <- "RNA"

# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")  # Use "^MT-" for human

# ============================================================================
# ADD SAMPLE METADATA
# ============================================================================

# Add sex information based on barcode suffix
# Adjust this based on your aggregation scheme
# -1 suffix = first sample, -2 suffix = second sample, etc.
seurat_obj$sex <- ifelse(grepl("-1$", colnames(seurat_obj)), "female", 
                         ifelse(grepl("-2$", colnames(seurat_obj)), "male", NA))

# ============================================================================
# VISUALIZE QC METRICS (PRE-FILTERING)
# ============================================================================

cat("Generating pre-filtering QC plots...\n")

# RNA QC violin plots
p1 <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
) & scale_y_continuous(labels = comma)

ggsave(
  file.path(figures_dir, "RNA_QC_prefilter.pdf"),
  plot = p1,
  width = 12,
  height = 4
)

# ATAC QC violin plots
p2 <- VlnPlot(
  seurat_obj,
  features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
  pt.size = 0.1,
  ncol = 4
) & scale_y_continuous(labels = comma)

ggsave(
  file.path(figures_dir, "ATAC_QC_prefilter.pdf"),
  plot = p2,
  width = 14,
  height = 4
)

# Combined QC plot
p3 <- VlnPlot(
  object = seurat_obj,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
) & scale_y_continuous(labels = comma)

ggsave(
  file.path(figures_dir, "Combined_QC_prefilter.pdf"),
  plot = p3,
  width = 14,
  height = 4
)

# QC by sex (if applicable)
if ("sex" %in% colnames(seurat_obj@meta.data)) {
  Idents(seurat_obj) <- 'sex'
  
  p4 <- VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    pt.size = 0.1,
    ncol = 3,
    cols = c("female" = "#F39AC9", "male" = "#4A6FE3")
  ) & scale_y_continuous(labels = comma)
  
  ggsave(
    file.path(figures_dir, "RNA_QC_bySex_prefilter.pdf"),
    plot = p4,
    width = 12,
    height = 4
  )
}

# ============================================================================
# CELL FILTERING
# ============================================================================

cat("Applying filtering thresholds...\n")

# Define filtering thresholds
# Adjust these based on your data distribution
filtering_params <- list(
  nCount_RNA_max = 35000,
  nCount_ATAC_max = 500000,
  percent_mt_max = 15,
  nucleosome_signal_min = 0.1,
  nucleosome_signal_max = 1.2,
  TSS_enrichment_min = 2.5,
  TSS_enrichment_max = 8
)

# Apply filters
seurat_filtered <- subset(
  x = seurat_obj,
  subset = nCount_ATAC < filtering_params$nCount_ATAC_max &
    nCount_RNA < filtering_params$nCount_RNA_max &
    percent.mt < filtering_params$percent_mt_max &
    nucleosome_signal > filtering_params$nucleosome_signal_min &
    nucleosome_signal < filtering_params$nucleosome_signal_max &
    TSS.enrichment > filtering_params$TSS_enrichment_min &
    TSS.enrichment < filtering_params$TSS_enrichment_max
)

cat("\nFiltering Results:\n")
cat("Cells before filtering:", ncol(seurat_obj), "\n")
cat("Cells after filtering:", ncol(seurat_filtered), "\n")
cat("Cells removed:", ncol(seurat_obj) - ncol(seurat_filtered), "\n")
cat("Percent retained:", round(ncol(seurat_filtered) / ncol(seurat_obj) * 100, 2), "%\n\n")

# ============================================================================
# VISUALIZE QC METRICS (POST-FILTERING)
# ============================================================================

cat("Generating post-filtering QC plots...\n")

# RNA QC violin plots (post-filter)
p5 <- VlnPlot(
  seurat_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
) & scale_y_continuous(labels = comma)

ggsave(
  file.path(figures_dir, "RNA_QC_postfilter.pdf"),
  plot = p5,
  width = 12,
  height = 4
)

# ATAC QC violin plots (post-filter)
p6 <- VlnPlot(
  seurat_filtered,
  features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
  pt.size = 0.1,
  ncol = 4
) & scale_y_continuous(labels = comma)

ggsave(
  file.path(figures_dir, "ATAC_QC_postfilter.pdf"),
  plot = p6,
  width = 14,
  height = 4
)

# ============================================================================
# CALL PEAKS WITH MACS2
# ============================================================================

cat("\nNote: Peak calling with MACS2 should be done separately.\n")
cat("See 'macs2_peak_calling.R' script for peak calling.\n")

# Uncomment and adjust if MACS2 is available
# DefaultAssay(seurat_filtered) <- "ATAC"
# peaks <- CallPeaks(seurat_filtered, macs2.path = '/path/to/macs2')
# saveRDS(peaks, file.path(output_dir, "macs2_peaks.rds"))

# ============================================================================
# SAVE FILTERED OBJECT
# ============================================================================

cat("Saving filtered Seurat object...\n")

# Save unfiltered object (for reference)
saveRDS(
  seurat_obj,
  file.path(output_dir, "seurat_object_unfiltered.rds")
)

# Save filtered object
saveRDS(
  seurat_filtered,
  file.path(output_dir, "seurat_object_filtered.rds")
)

# ============================================================================
# GENERATE QC SUMMARY TABLE
# ============================================================================

cat("Generating QC summary statistics...\n")

qc_summary <- data.frame(
  Metric = c(
    "Total cells (pre-filter)",
    "Total cells (post-filter)",
    "Cells removed",
    "Percent retained",
    "Median RNA features",
    "Median RNA counts",
    "Median ATAC features",
    "Median ATAC counts",
    "Median TSS enrichment",
    "Median nucleosome signal",
    "Median percent MT"
  ),
  Value = c(
    ncol(seurat_obj),
    ncol(seurat_filtered),
    ncol(seurat_obj) - ncol(seurat_filtered),
    round(ncol(seurat_filtered) / ncol(seurat_obj) * 100, 2),
    median(seurat_filtered$nFeature_RNA),
    median(seurat_filtered$nCount_RNA),
    median(seurat_filtered$nFeature_ATAC),
    median(seurat_filtered$nCount_ATAC),
    round(median(seurat_filtered$TSS.enrichment), 2),
    round(median(seurat_filtered$nucleosome_signal), 2),
    round(median(seurat_filtered$percent.mt), 2)
  )
)

write.csv(
  qc_summary,
  file.path(output_dir, "qc_summary_statistics.csv"),
  row.names = FALSE
)

print(qc_summary)

# ============================================================================
# SESSION INFO
# ============================================================================

cat("\n=== QC and Filtering Complete ===\n")
cat("\nOutput files saved to:", output_dir, "\n")
cat("Figures saved to:", figures_dir, "\n\n")

sessionInfo()