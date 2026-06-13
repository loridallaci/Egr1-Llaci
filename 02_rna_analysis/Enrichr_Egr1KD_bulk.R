# =============================================================================
# Enrichr Analysis: Egr1 KD vs WT DE genes (Male and Female separately)
# =============================================================================

library(enrichR)
library(ggplot2)
library(dplyr)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}

if (websiteLive) dbs <- listEnrichrDbs()
if (websiteLive) head(dbs)
dbs

# --- Directories --------------------------------------------------------------

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1_bulk_RNAseq/Egr1_KD_030525/diffexp_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers"
base_dir         <- "output"
# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/enrichment/Egr1KD"
OutputDirectory  <- "output/enrichment/Egr1KD"
dir.create(OutputDirectory, showWarnings = FALSE, recursive = TRUE)

dbs <- c("WikiPathways_2019_Mouse", "GO_Biological_Process_2023", "GO_Molecular_Function_2023")

# =============================================================================
# MALE
# =============================================================================

de <- read.csv(
  file.path(base_dir, 'Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt'),
  sep = '\t'
)

de_KD <- filter(de, log2FoldChange >= 0.5  & pvalue <= 0.05)   # KD-enriched
de_WT <- filter(de, log2FoldChange <= -0.5 & pvalue <= 0.05)   # WT-enriched

cat("Male KD-enriched genes:", nrow(de_KD), "\n")
cat("Male WT-enriched genes:", nrow(de_WT), "\n")

# --- Male KD-enriched ---------------------------------------------------------

enriched <- enrichr(de_KD$SYMBOL, dbs)

for (enrichr_db in dbs) {
  result <- enriched[[enrichr_db]]
  if (is.null(result) || nrow(result) == 0) {
    cat("  Skipping", enrichr_db, "- no results\n"); next
  }
  enrichr_output_df <- paste0(OutputDirectory, "/Male_Egr1KD_enriched_", enrichr_db, ".txt")
  write.table(result, enrichr_output_df, quote = FALSE, row.names = TRUE, sep = '\t')
  
  enrichr_output_figure <- paste0(OutputDirectory, "/Male_Egr1KD_enriched_", enrichr_db, ".pdf")
  p <- plotEnrich(result, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value") +
    ggtitle(paste0("Male Egr1 KD-enriched\n", enrichr_db)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  pdf(file = enrichr_output_figure, width = 7, height = 7, onefile = TRUE, useDingbats = FALSE)
  print(p)
  invisible(dev.off())
}

# --- Male WT-enriched ---------------------------------------------------------

enriched <- enrichr(de_WT$SYMBOL, dbs)

for (enrichr_db in dbs) {
  result <- enriched[[enrichr_db]]
  if (is.null(result) || nrow(result) == 0) {
    cat("  Skipping", enrichr_db, "- no results\n"); next
  }
  enrichr_output_df <- paste0(OutputDirectory, "/Male_WT_enriched_", enrichr_db, ".txt")
  write.table(result, enrichr_output_df, quote = FALSE, row.names = TRUE, sep = '\t')
  
  enrichr_output_figure <- paste0(OutputDirectory, "/Male_WT_enriched_", enrichr_db, ".pdf")
  p <- plotEnrich(result, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value") +
    ggtitle(paste0("Male WT-enriched\n", enrichr_db)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  pdf(file = enrichr_output_figure, width = 7, height = 7, onefile = TRUE, useDingbats = FALSE)
  print(p)
  invisible(dev.off())
}

# =============================================================================
# FEMALE
# =============================================================================

de <- read.csv(
  file.path(base_dir, 'Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt'),
  sep = '\t'
)

de_KD <- filter(de, log2FoldChange >= 0.5  & pvalue <= 0.05)
de_WT <- filter(de, log2FoldChange <= -0.5 & pvalue <= 0.05)

cat("Female KD-enriched genes:", nrow(de_KD), "\n")
cat("Female WT-enriched genes:", nrow(de_WT), "\n")

# --- Female KD-enriched -------------------------------------------------------

enriched <- enrichr(de_KD$SYMBOL, dbs)

for (enrichr_db in dbs) {
  result <- enriched[[enrichr_db]]
  if (is.null(result) || nrow(result) == 0) {
    cat("  Skipping", enrichr_db, "- no results\n"); next
  }
  enrichr_output_df <- paste0(OutputDirectory, "/Female_Egr1KD_enriched_", enrichr_db, ".txt")
  write.table(result, enrichr_output_df, quote = FALSE, row.names = TRUE, sep = '\t')
  
  enrichr_output_figure <- paste0(OutputDirectory, "/Female_Egr1KD_enriched_", enrichr_db, ".pdf")
  p <- plotEnrich(result, showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value") +
    ggtitle(paste0("Female Egr1 KD-enriched\n", enrichr_db)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  pdf(file = enrichr_output_figure, width = 7, height = 7, onefile = TRUE, useDingbats = FALSE)
  print(p)
  invisible(dev.off())
}

# --- Female WT-enriched -------------------------------------------------------

enriched <- enrichr(de_WT$SYMBOL, dbs)

for (enrichr_db in dbs) {
  result <- enriched[[enrichr_db]]
  if (is.null(result) || nrow(result) == 0) {
    cat("  Skipping", enrichr_db, "- no results\n"); next
  }
  enrichr_output_df <- paste0(OutputDirectory, "/Female_WT_enriched_", enrichr_db, ".txt")
  write.table(result, enrichr_output_df, quote = FALSE, row.names = TRUE, sep = '\t')
  
  enrichr_output_figure <- paste0(OutputDirectory, "/Female_WT_enriched_", enrichr_db, ".pdf")
  p <- plotEnrich(result, showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value") +
    ggtitle(paste0("Female WT-enriched\n", enrichr_db)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  pdf(file = enrichr_output_figure, width = 7, height = 7, onefile = TRUE, useDingbats = FALSE)
  print(p)
  invisible(dev.off())
}

message("\nDone. All Enrichr results saved to: ", OutputDirectory)