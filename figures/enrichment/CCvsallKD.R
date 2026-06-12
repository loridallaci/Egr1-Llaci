# =============================================================================
# Enrichr Analysis: Egr1 KD vs WT DE genes (Male and Female separately)
# =============================================================================
# For each comparison, runs Enrichr on:
#   - KD-enriched genes  (log2FC >= 0.5, pvalue <= 0.05)
#   - WT-enriched genes  (log2FC <= -0.5, pvalue <= 0.05)
# Exports .txt results and high-quality PDF plots
# =============================================================================

library(enrichR)
library(ggplot2)
library(dplyr)

# --- Directories --------------------------------------------------------------

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1_bulk_RNAseq/Egr1_KD_030525/diffexp_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers"
base_dir    <- "data"
# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output"
output_dir  <- "output"

# --- Enrichr databases --------------------------------------------------------

dbs <- c("WikiPathways_2019_Mouse", "GO_Molecular_Function_2017b", "GO_Biological_Process_2023")

# --- Comparisons --------------------------------------------------------------

comparisons <- list(
  list(
    file        = "Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt",
    up_label    = "Male_Egr1KD",       # KD-enriched (positive log2FC)
    down_label  = "Male_WT",           # WT-enriched (negative log2FC)
    sample_name = "Male_Egr1KDg3_vs_Male_NoTreatg1"
  ),
  list(
    file        = "Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt",
    up_label    = "Female_Egr1KD",
    down_label  = "Female_WT",
    sample_name = "Female_Egr1KDg3_vs_Female_NoTreatg1"
  )
)

# =============================================================================
# Loop over comparisons
# =============================================================================

for (comp in comparisons) {
  
  cat("\n============================================================\n")
  cat("Processing:", comp$sample_name, "\n")
  
  # Load data
  de <- read.csv(
    file.path(base_dir, comp$file),
    sep         = "\t",
    row.names   = 1
  )
  de$SYMBOL <- rownames(de)
  
  # Replace p-value of 0
  de$pvalue[de$pvalue == 0] <- 1e-300
  de <- de[!is.na(de$pvalue), ]
  
  # --- Split into KD-enriched and WT-enriched ---------------------------------
  
  de_kd <- de %>% filter(log2FoldChange >=  0.5 & pvalue <= 0.05)
  de_wt <- de %>% filter(log2FoldChange <= -0.5 & pvalue <= 0.05)
  
  cat("KD-enriched genes:", nrow(de_kd), "\n")
  cat("WT-enriched genes:", nrow(de_wt), "\n")
  
  # --- Save gene lists --------------------------------------------------------
  
  write.table(de_kd$SYMBOL,
              file.path(output_dir, paste0(comp$sample_name, "_KDenriched_genes.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  write.table(de_wt$SYMBOL,
              file.path(output_dir, paste0(comp$sample_name, "_WTenriched_genes.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # ===========================================================================
  # Run Enrichr for each gene set x each database
  # ===========================================================================
  
  gene_sets <- list(
    list(genes = de_kd$SYMBOL, label = paste0(comp$sample_name, "_KDenriched")),
    list(genes = de_wt$SYMBOL, label = paste0(comp$sample_name, "_WTenriched"))
  )
  
  for (gs in gene_sets) {
    
    if (length(gs$genes) == 0) {
      cat("Skipping", gs$label, "- no genes\n")
      next
    }
    
    cat("\n  Running Enrichr for:", gs$label, "\n")
    enriched <- enrichr(gs$genes, dbs)
    
    for (db in dbs) {
      
      cat("    Database:", db, "\n")
      result <- enriched[[db]]
      
      if (is.null(result) || nrow(result) == 0) {
        cat("    No results for", db, "\n")
        next
      }
      
      # Save txt
      txt_path <- file.path(output_dir, paste0(gs$label, "_", db, ".txt"))
      write.table(result, txt_path, quote = FALSE, row.names = TRUE, sep = "\t")
      cat("    Saved txt:", txt_path, "\n")
      
      # Save PDF
      pdf_path <- file.path(output_dir, paste0(gs$label, "_", db, ".pdf"))
      p <- plotEnrich(result, showTerms = 20, numChar = 40,
                      y = "Count", orderBy = "P.value") +
        ggtitle(paste0(gs$label, "\n", db)) +
        theme(plot.title = element_text(hjust = 0.5, size = 9))
      
      pdf(file = pdf_path, width = 9, height = 7, onefile = TRUE, useDingbats = FALSE)
      print(p)
      invisible(dev.off())
      cat("    Saved pdf:", pdf_path, "\n")
    }
  }
}

message("\nDone. All Enrichr results saved.")