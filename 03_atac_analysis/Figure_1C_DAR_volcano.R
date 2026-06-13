# =============================================================================
# 03_atac_analysis: Differentially Accessible Regions (DARs) Volcano Plot
# =============================================================================
# Description:
#   Annotates differentially accessible peaks with nearest genes using
#   EnsDb.Mmusculus.v79, classifies peaks as male- or female-enriched,
#   and generates a publication-ready volcano plot.
#
# Input:
#   lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv
#
# Output:
#   - lot6_DAR_peaks_annotated.csv
#   - figure_DAR_volcano_male_vs_female.pdf
# =============================================================================

# --- Libraries ----------------------------------------------------------------

library(Seurat)
library(Signac)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(ensembldb)
library(stringr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(GenomicFeatures)
library(dplyr)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# --- Load data ----------------------------------------------------------------

da_peaks <- read.csv("output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks_logfc0.csv")  # original: "/home/lllaci/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks_logfc0.csv"

# --- Annotate peaks with nearest gene -----------------------------------------

# Convert peak coordinates (format: chr-start-end) to GRanges object
parts       <- str_split_fixed(da_peaks$X, "-", 3)
da_peaks_gr <- GRanges(
  seqnames = parts[, 1],
  ranges   = IRanges(
    start = as.numeric(parts[, 2]),
    end   = as.numeric(parts[, 3])
  )
)

# Remove "chr" prefix to match Ensembl-style chromosome names (e.g. "1" not "chr1")
seqlevels(da_peaks_gr) <- sub("^chr", "", seqlevels(da_peaks_gr))

# Get gene coordinates from Ensembl mouse annotation
genes_gr <- genes(EnsDb.Mmusculus.v79)

# Find nearest gene for each peak
nearest_idx      <- nearest(da_peaks_gr, genes_gr)
nearest_genes_gr <- genes_gr[nearest_idx]

# Add gene annotations to data frame
da_peaks$nearest_gene              <- mcols(nearest_genes_gr)$gene_name
da_peaks$nearest_gene_id           <- mcols(nearest_genes_gr)$gene_id
da_peaks$distance_to_nearest_gene  <- mcols(
  distanceToNearest(da_peaks_gr, genes_gr)
)$distance

# Save annotated peaks
write.csv(da_peaks, "output/lot6_DAR_peaks_annotated_logfc0.csv", row.names = FALSE)  # original: "/home/lllaci/lot6_DAR_peaks_annotated_logfc0.csv"

# --- Classify peaks -----------------------------------------------------------

da_peaks$significance <- "Not Significant"
da_peaks$significance[
  da_peaks$avg_log2FC >=  0.5 & da_peaks$p_val_adj <= 0.05
] <- "Male-enriched"
da_peaks$significance[
  da_peaks$avg_log2FC <= -0.5 & da_peaks$p_val_adj <= 0.05
] <- "Female-enriched"
da_peaks$significance <- factor(
  da_peaks$significance,
  levels = c("Male-enriched", "Female-enriched", "Not Significant")
)

# --- Summary ------------------------------------------------------------------

n_male   <- sum(da_peaks$significance == "Male-enriched")
n_female <- sum(da_peaks$significance == "Female-enriched")
cat("Male-enriched peaks:  ", n_male,           "\n")
cat("Female-enriched peaks:", n_female,          "\n")
cat("Total significant:    ", n_male + n_female, "\n")

# --- Volcano plot -------------------------------------------------------------

# Label only significant peaks
sig_to_label <- da_peaks %>%
  dplyr::filter(significance != "Not Significant")

# Cap -log10(p_val_adj) to handle p = 0 (avoids Inf on y axis)
da_peaks$log10p <- pmin(-log10(da_peaks$p_val_adj + 1e-300), 300)

p <- ggplot(da_peaks, aes(
  x     = avg_log2FC,
  y     = log10p,
  color = significance
)) +
  geom_point(
    alpha = 0.7,
    size  = 1.8,
    shape = 16
  ) +
  geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_color_manual(
    values = c(
      "Male-enriched"   = "dodgerblue3",
      "Female-enriched" = "pink",
      "Not Significant" = "grey70"
    ),
    name = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    axis.text       = element_text(color = "black"),
    axis.title      = element_text(face = "bold"),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(
    title = "Differentially Accessible Peaks: Male vs Female",
    x     = expression(log[2]~"Fold Change"),
    y     = expression(-log[10]~"(adjusted p-value)")
  )

cairo_pdf("output/DAR_volcano_male_vs_female_logfc0.pdf", width = 7, height = 6)  # original: "/home/lllaci/DAR_volcano_male_vs_female_logfc0.pdf"
print(p)
dev.off()

message("Done. Annotated peaks and volcano plot saved.")













