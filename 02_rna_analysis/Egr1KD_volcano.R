# =============================================================================
# RNA-seq Volcano Plots: Egr1 KD vs WT (Male and Female separately)
# =============================================================================
# Colors:
#   Male   Egr1 KD  -> #00BFFF (deep sky blue)
#   Female Egr1 KD  -> #E75480 (dark pink)
#   Male   Control  -> darkgreen
#   Female Control  -> lightgreen
# All points are circles. log2FC capped at +/- 3 for display.
# =============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)

# --- Directories --------------------------------------------------------------

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1_bulk_RNAseq/Egr1_KD_030525/diffexp_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers"
base_dir              <- "output"
# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output"
figureOutputDirectory <- "output"
dir.create(figureOutputDirectory, showWarnings = FALSE, recursive = TRUE)

LOG2FC_CAP <- 4

# --- Genes to label -----------------------------------------------------------

genes_to_label <- c("Egr1")

# --- Comparisons --------------------------------------------------------------

comparisons <- list(
  list(
    file       = "Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt",
    up_label   = "Male Egr1KD gRNA3",
    down_label = "Male Control gRNA1",
    up_color   = "#00BFFF",
    down_color = "darkgreen",
    title      = "Male Egr1KD gRNA3 vs Male NoTreatment gRNA1",
    out_pdf    = "Male_Egr1KD_gRNA3_vs_Male_NoTreatg1_volcano_plot_vst_filtered_Figure3C.pdf"
  ),
  list(
    file       = "Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt",
    up_label   = "Female Egr1 gRNA3",
    down_label = "Female Control gRNA1",
    up_color   = "#FF66CC",
    down_color = "#74AC64",
    title      = "Female Egr1KD gRNA3 vs Female NoTreatment gRNA1",
    out_pdf    = "Female_Egr1KD_gRNA3_vs_Female_NoTreatg1_volcano_plot_vst_filtered_Figure3C.pdf"
  )
)

# --- Loop ---------------------------------------------------------------------

for (comp in comparisons) {
  
  cat("\n============================================================\n")
  cat("Processing:", comp$title, "\n")
  
  # Load data
  de <- read.csv(
    file.path(base_dir, comp$file),
    sep = "\t"
  )
  
  # Replace p-value of 0 with very small number
  de$pvalue[de$pvalue == 0] <- 1e-300
  
  # Remove NA p-values
  de <- de[!is.na(de$pvalue), ]
  
  # --- Cap log2FC for display only (classification uses real values) ----------
  
  de$log2fc_plot <- pmax(pmin(de$log2FoldChange, LOG2FC_CAP), -LOG2FC_CAP)
  
  # --- Classify ---------------------------------------------------------------
  
  de$diffexpressed <- "NO"
  de$diffexpressed[de$log2FoldChange >=  0.5 & de$pvalue <= 0.05] <- comp$up_label
  de$diffexpressed[de$log2FoldChange <= -0.5 & de$pvalue <= 0.05] <- comp$down_label
  
  # Summary
  n_up   <- sum(de$diffexpressed == comp$up_label)
  n_down <- sum(de$diffexpressed == comp$down_label)
  cat(comp$up_label,   ":", n_up,   "\n")
  cat(comp$down_label, ":", n_down, "\n")
  
  # --- Labels -----------------------------------------------------------------
  
  de$delabel <- ifelse(de$SYMBOL %in% genes_to_label, de$SYMBOL, NA)
  
  # --- Colors -----------------------------------------------------------------
  
  mycolors <- setNames(
    c(comp$down_color, comp$up_color, "grey60"),
    c(comp$down_label, comp$up_label, "NO")
  )
  
  # --- Plot -------------------------------------------------------------------
  
  p <- ggplot(de, aes(
    x     = log2fc_plot,
    y     = -log10(pvalue),
    col   = diffexpressed,
    label = delabel
  )) +
    geom_point(alpha = 0.8, shape = 16) +
    geom_text_repel(
      box.padding  = 2.5,
      max.overlaps = Inf,
      nudge_x      = 0.1,
      nudge_y      = 0.1,
      size         = 4,
      show.legend  = FALSE,
      fontface     = "italic"
    ) +
    scale_color_manual(values = mycolors, name = "Expression Change") +
    scale_x_continuous(
      limits = c(-LOG2FC_CAP, LOG2FC_CAP),
      breaks = seq(-LOG2FC_CAP, LOG2FC_CAP, by = 1)
    ) +
    theme_minimal() +
    labs(
      title    = comp$title,
      subtitle = paste0(comp$up_label, ": ", n_up,
                        "  |  ", comp$down_label, ": ", n_down,
                        "  |  log2FC capped at \u00b1", LOG2FC_CAP),
      x        = "Log2 Fold Change",
      y        = "-Log10 P-value"
    ) +
    theme(
      plot.title    = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40")
    )
  
  # Save PDF
  out_path <- file.path(figureOutputDirectory, comp$out_pdf)
  cairo_pdf(out_path, width = 8, height = 6)
  print(p)
  dev.off()
  cat("Saved:", out_path, "\n")
  
  # Count significant genes
  de_down <- filter(de, log2FoldChange <= -0.5 & pvalue <= 0.05)
  de_up   <- filter(de, log2FoldChange >=  0.5 & pvalue <= 0.05)
  cat("Down:", nrow(de_down), "  Up:", nrow(de_up), "\n")
}

message("\nDone. All volcano plots saved.")

#CONTROLS
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- Directories --------------------------------------------------------------

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1_bulk_RNAseq/Egr1_KD_030525/diffexp_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers"
base_dir              <- "output"
# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output"
figureOutputDirectory <- "output"
dir.create(figureOutputDirectory, showWarnings = FALSE, recursive = TRUE)

LOG2FC_CAP <- 10

# --- Genes to label -----------------------------------------------------------

genes_to_label <- c("Egr1")

# --- Comparisons --------------------------------------------------------------

comparisons <- list(
  list(
    file       = "Male_NoTreatg1_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt",
    up_label   = "Male Egr1 WT",
    down_label = "Female Egr1 WT",
    up_color   = "darkgreen",
    down_color = "#74AC64",#74AC64
    title      = "Male NoTreatment gRNA1 vs Female NoTreatment gRNA1",
    out_pdf    = "Male_NoTreatg1_vs_Female_NoTreatg1_volcano_plot_vst_filtered_Figure3C.pdf"
  ),
  list(
    file       = "Male_Egr1KDg3_vs_Female_Egr1KDg3_DE_vst_filtered_091625.txt",
    up_label   = "Male Egr1 gRNA3",
    down_label = "Female Egr1 gRNA3",
    up_color   = "#00BFFF",
    down_color = "#FF66CC",
    title      = "Male Egr1KD gRNA3 vs Female Egr1KD gRNA3",
    out_pdf    = "Male_Egr1KD_gRNA3_vs_Female_Egr1KD_gRNA3_volcano_plot_vst_filtered_Figure3C.pdf"
  )
)

# --- Loop ---------------------------------------------------------------------

for (comp in comparisons) {
  
  cat("\n============================================================\n")
  cat("Processing:", comp$title, "\n")
  
  # Load data
  de <- read.csv(
    file.path(base_dir, comp$file),
    sep = "\t"
  )
  
  # Replace p-value of 0 with very small number
  de$pvalue[de$pvalue == 0] <- 1e-300
  
  # Remove NA p-values
  de <- de[!is.na(de$pvalue), ]
  
  # --- Cap log2FC for display only (classification uses real values) ----------
  
  de$log2fc_plot <- pmax(pmin(de$log2FoldChange, LOG2FC_CAP), -LOG2FC_CAP)
  
  # --- Classify ---------------------------------------------------------------
  
  de$diffexpressed <- "NO"
  de$diffexpressed[de$log2FoldChange >=  0.5 & de$pvalue <= 0.05] <- comp$up_label
  de$diffexpressed[de$log2FoldChange <= -0.5 & de$pvalue <= 0.05] <- comp$down_label
  
  # Summary
  n_up   <- sum(de$diffexpressed == comp$up_label)
  n_down <- sum(de$diffexpressed == comp$down_label)
  cat(comp$up_label,   ":", n_up,   "\n")
  cat(comp$down_label, ":", n_down, "\n")
  
  # --- Labels -----------------------------------------------------------------
  
  de$delabel <- ifelse(de$SYMBOL %in% genes_to_label, de$SYMBOL, NA)
  
  # --- Colors -----------------------------------------------------------------
  
  mycolors <- setNames(
    c(comp$down_color, comp$up_color, "grey60"),
    c(comp$down_label, comp$up_label, "NO")
  )
  
  # --- Plot -------------------------------------------------------------------
  
  p <- ggplot(de, aes(
    x     = log2fc_plot,
    y     = -log10(pvalue),
    col   = diffexpressed,
    label = delabel
  )) +
    geom_point(alpha = 0.8, shape = 16) +
    geom_text_repel(
      box.padding  = 2.5,
      max.overlaps = Inf,
      nudge_x      = 0.1,
      nudge_y      = 0.1,
      size         = 4,
      show.legend  = FALSE,
      fontface     = "italic"
    ) +
    scale_color_manual(values = mycolors, name = "Expression Change") +
    scale_x_continuous(
      limits = c(-LOG2FC_CAP, LOG2FC_CAP),
      breaks = seq(-LOG2FC_CAP, LOG2FC_CAP, by = 1)
    ) +
    theme_minimal() +
    labs(
      title    = comp$title,
      subtitle = paste0(comp$up_label, ": ", n_up,
                        "  |  ", comp$down_label, ": ", n_down,
                        "  |  log2FC capped at \u00b1", LOG2FC_CAP),
      x        = "Log2 Fold Change",
      y        = "-Log10 P-value"
    ) +
    theme(
      plot.title    = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40")
    )
  
  # Save PDF
  out_path <- file.path(figureOutputDirectory, comp$out_pdf)
  cairo_pdf(out_path, width = 8, height = 6)
  print(p)
  dev.off()
  cat("Saved:", out_path, "\n")
  
  # Count significant genes
  de_down <- filter(de, log2FoldChange <= -0.5 & pvalue <= 0.05)
  de_up   <- filter(de, log2FoldChange >=  0.5 & pvalue <= 0.05)
  cat("Down:", nrow(de_down), "  Up:", nrow(de_up), "\n")
}

message("\nDone. All volcano plots saved.")
