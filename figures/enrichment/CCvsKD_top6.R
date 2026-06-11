# =============================================================================
# Enrichr Analysis: CC 20kb NearestGene overlap with KD DE genes
# =============================================================================
# For Male and Female:
#   1. Read CC peaks file -> get Gene Name1 column
#   2. Read KD vs WT DE file -> split into KD-enriched and WT-enriched
#   3. Overlap CC genes with KD-enriched and WT-enriched separately
#   4. Run Enrichr on each overlap list
#   5. Save TXT + PDF to output/enrichment/Egr1CCvsKD
# =============================================================================

library(enrichR)
library(ggplot2)
library(dplyr)

# --- Directories --------------------------------------------------------------

peaks_dir  <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output"
de_dir     <- "C:/Users/loril/Documents/Egr1/Egr1_bulk_RNAseq/Egr1_KD_030525/diffexp_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers"
output_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/enrichment/Egr1CCvsKD"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Enrichr databases --------------------------------------------------------

dbs <- c("GO_Biological_Process_2023")

# --- Comparisons --------------------------------------------------------------

comparisons <- list(
  list(
    peaks_file = "Male_Egr1CC_peaks_20kbThreshhold_091125_111225.txt",
    de_file    = "Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt",
    sex        = "Male"
  ),
  list(
    peaks_file = "Female_Egr1CC_peaks_20kbThreshhold_091125_111225.txt",
    de_file    = "Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt",
    sex        = "Female"
  )
)

# =============================================================================
# Helper: run Enrichr and save TXT + PDF
# =============================================================================

run_enrichr_and_save <- function(genes, label, output_dir, dbs) {
  
  if (length(genes) == 0) {
    cat("  No genes for", label, "- skipping\n")
    return(invisible(NULL))
  }
  
  cat("  Running Enrichr for:", label, "(", length(genes), "genes )\n")
  enriched <- enrichr(genes, dbs)
  
  for (enrichr_db in dbs) {
    result <- enriched[[enrichr_db]]
    
    if (is.null(result) || nrow(result) == 0) {
      cat("    Skipping", enrichr_db, "- no results\n"); next
    }
    
    # Save TXT
    txt_path <- file.path(output_dir, paste0(label, "_", enrichr_db, ".txt"))
    write.table(result, txt_path, quote = FALSE, row.names = TRUE, sep = '\t')
    cat("    Saved:", basename(txt_path), "\n")
    
    # Save PDF
    pdf_path <- file.path(output_dir, paste0(label, "_", enrichr_db, ".pdf"))
    p <- plotEnrich(result, showTerms = 6, numChar = 60, y = "Count", orderBy = "P.value") +
      ggtitle(paste0(label, "\n", enrichr_db)) +
      theme(plot.title = element_text(hjust = 0.5, size = 8))
    pdf(file = pdf_path, width = 9, height = 4, onefile = TRUE, useDingbats = FALSE)
    print(p)
    invisible(dev.off())
    cat("    Saved:", basename(pdf_path), "\n")
  }
}

# =============================================================================
# Loop over Male and Female
# =============================================================================

for (comp in comparisons) {
  
  cat("\n============================================================\n")
  cat("Processing:", comp$sex, "\n")
  
  # --- Load CC peaks file and get Gene Name1 ----------------------------------
  
  peaks <- read.csv(
    file.path(peaks_dir, comp$peaks_file),
    sep = '\t', header = TRUE, check.names = FALSE
  )
  
  cc_genes <- peaks[["Gene Name1"]]
  cc_genes <- as.character(cc_genes[!is.na(cc_genes) & cc_genes != ""])
  cc_genes <- unique(cc_genes)
  cat("CC peak nearest genes:", length(cc_genes), "\n")
  
  # --- Load DE file -----------------------------------------------------------
  
  de <- read.csv(
    file.path(de_dir, comp$de_file),
    sep = '\t'
  )
  de$pvalue[de$pvalue == 0] <- 1e-300
  de <- de[!is.na(de$pvalue), ]
  
  de_KD <- filter(de, log2FoldChange >= 0.5  & pvalue <= 0.05)
  de_WT <- filter(de, log2FoldChange <= -0.5 & pvalue <= 0.05)
  
  cat("KD-enriched DE genes:", nrow(de_KD), "\n")
  cat("WT-enriched DE genes:", nrow(de_WT), "\n")
  
  # --- Overlap ----------------------------------------------------------------
  
  overlap_KD <- intersect(cc_genes, de_KD$SYMBOL)
  overlap_WT <- intersect(cc_genes, de_WT$SYMBOL)
  
  cat("Overlap with KD-enriched:", length(overlap_KD), "\n")
  cat("Overlap with WT-enriched:", length(overlap_WT), "\n")
  
  # Save overlap gene lists
  write.table(overlap_KD,
              file.path(output_dir, paste0(comp$sex, "_CC20kb_overlap_KDenriched_top6.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  write.table(overlap_WT,
              file.path(output_dir, paste0(comp$sex, "_CC20kb_overlap_WTenriched_top6.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # --- Run Enrichr on overlaps ------------------------------------------------
  
  run_enrichr_and_save(
    genes      = overlap_KD,
    label      = paste0(comp$sex, "_CC20kb_overlap_KDenriched_top6"),
    output_dir = output_dir,
    dbs        = dbs
  )
  
  run_enrichr_and_save(
    genes      = overlap_WT,
    label      = paste0(comp$sex, "_CC20kb_overlap_WTenriched_top6"),
    output_dir = output_dir,
    dbs        = dbs
  )
}

message("\nDone. All results saved to: ", output_dir)