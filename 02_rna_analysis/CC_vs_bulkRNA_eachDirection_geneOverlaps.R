# =============================================================================
# CC Peaks vs DE genes overlap analysis
# =============================================================================

library(dplyr)

base_dir   <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"
output_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/enrichment/Egr1CCvsKD/"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

de_dir <- "C:/Users/loril/Documents/Egr1/Egr1_bulk_RNAseq/Egr1_KD_030525/diffexp_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers"

# =============================================================================
# STEP 1: Load CC peaks and extract Gene Name1
# =============================================================================

cc_peaks <- read.table(
  file.path(base_dir, "Male_Egr1CC_peaks_20kbThreshhold_091125_111225.txt"),
  sep = "\t", header = TRUE
)

cc_genes <- unique(toupper(cc_peaks$Gene.Name1))
cc_genes <- cc_genes[!is.na(cc_genes) & cc_genes != ""]
cat("Unique CC peak genes:", length(cc_genes), "\n")

# =============================================================================
# STEP 2: Load DE genes and split by direction
# =============================================================================

de <- read.csv(
  file.path(de_dir, "Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt"),
  sep = "\t"
)

de$SYMBOL_UPPER <- toupper(de$SYMBOL)

de_KD <- filter(de, log2FoldChange >= 0.5  & pvalue <= 0.05)  # upregulated in KD
de_WT <- filter(de, log2FoldChange <= -0.5 & pvalue <= 0.05)  # upregulated in WT

cat("Male KD-enriched genes:", nrow(de_KD), "\n")
cat("Male WT-enriched genes:", nrow(de_WT), "\n")

# =============================================================================
# STEP 3: Find overlaps
# =============================================================================

overlap_KD <- intersect(cc_genes, toupper(de_KD$SYMBOL))
overlap_WT <- intersect(cc_genes, toupper(de_WT$SYMBOL))

cat("\nCC peaks near KD-enriched DE genes:", length(overlap_KD), "\n")
print(overlap_KD)

cat("\nCC peaks near WT-enriched DE genes:", length(overlap_WT), "\n")
print(overlap_WT)

# =============================================================================
# STEP 4: Save overlap tables with full DE stats
# =============================================================================

overlap_KD_df <- de_KD %>%
  filter(toupper(SYMBOL) %in% overlap_KD) %>%
  left_join(
    cc_peaks %>%
      filter(toupper(Gene.Name1) %in% overlap_KD) %>%
      select(Chr, Start, End, Gene.Name1, Experiment.Insertions,
             TPH.Experiment, pvalue_adj.Reference),
    by = c("SYMBOL" = "Gene.Name1")
  )

overlap_WT_df <- de_WT %>%
  filter(toupper(SYMBOL) %in% overlap_WT) %>%
  left_join(
    cc_peaks %>%
      filter(toupper(Gene.Name1) %in% overlap_WT) %>%
      select(Chr, Start, End, Gene.Name1, Experiment.Insertions,
             TPH.Experiment, pvalue_adj.Reference),
    by = c("SYMBOL" = "Gene.Name1")
  )

write.csv(overlap_KD_df,
          file.path(output_dir, "overlap_MaleCC_vs_MaleKD_enriched.csv"),
          row.names = FALSE)

write.csv(overlap_WT_df,
          file.path(output_dir, "overlap_MaleCC_vs_MaleWT_enriched.csv"),
          row.names = FALSE)

cat("\nSaved overlap tables to:", output_dir, "\n")

# =============================================================================
# STEP 5: Gene sets for survival z-score analysis
# =============================================================================

geneset_KD <- overlap_KD  # CC peaks + upregulated when Egr1 KD → Egr1-repressed
geneset_WT <- overlap_WT  # CC peaks + downregulated when Egr1 KD → direct Egr1 targets

cat("\nGene set (CC + KD-enriched):", paste(geneset_KD, collapse = ", "), "\n")
cat("Gene set (CC + WT-enriched):", paste(geneset_WT, collapse = ", "), "\n")


# Do the same for females
# =============================================================================
# CC Peaks vs DE genes overlap analysis - FEMALE
# =============================================================================

library(dplyr)

base_dir   <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"
output_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/enrichment/Egr1CCvsKD/"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

de_dir <- "C:/Users/loril/Documents/Egr1/Egr1_bulk_RNAseq/Egr1_KD_030525/diffexp_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers"

# =============================================================================
# STEP 1: Load CC peaks and extract Gene Name1
# =============================================================================

cc_peaks <- read.table(
  file.path(base_dir, "Female_Egr1CC_peaks_20kbThreshhold_091125_111225.txt"),
  sep = "\t", header = TRUE
)

cc_genes <- unique(toupper(cc_peaks$Gene.Name1))
cc_genes <- cc_genes[!is.na(cc_genes) & cc_genes != ""]
cat("Unique CC peak genes:", length(cc_genes), "\n")

# =============================================================================
# STEP 2: Load DE genes and split by direction
# =============================================================================

de <- read.csv(
  file.path(de_dir, "Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt"),
  sep = "\t"
)

de$SYMBOL_UPPER <- toupper(de$SYMBOL)

de_KD <- filter(de, log2FoldChange >= 0.5  & pvalue <= 0.05)  # upregulated in KD
de_WT <- filter(de, log2FoldChange <= -0.5 & pvalue <= 0.05)  # upregulated in WT

cat("Female KD-enriched genes:", nrow(de_KD), "\n")
cat("Female WT-enriched genes:", nrow(de_WT), "\n")

# =============================================================================
# STEP 3: Find overlaps
# =============================================================================

overlap_KD <- intersect(cc_genes, toupper(de_KD$SYMBOL))
overlap_WT <- intersect(cc_genes, toupper(de_WT$SYMBOL))

cat("\nCC peaks near KD-enriched DE genes:", length(overlap_KD), "\n")
print(overlap_KD)

cat("\nCC peaks near WT-enriched DE genes:", length(overlap_WT), "\n")
print(overlap_WT)

# =============================================================================
# STEP 4: Save overlap tables with full DE stats
# =============================================================================

overlap_KD_df <- de_KD %>%
  filter(toupper(SYMBOL) %in% overlap_KD) %>%
  left_join(
    cc_peaks %>%
      filter(toupper(Gene.Name1) %in% overlap_KD) %>%
      select(Chr, Start, End, Gene.Name1, Experiment.Insertions,
             TPH.Experiment, pvalue_adj.Reference),
    by = c("SYMBOL" = "Gene.Name1")
  )

overlap_WT_df <- de_WT %>%
  filter(toupper(SYMBOL) %in% overlap_WT) %>%
  left_join(
    cc_peaks %>%
      filter(toupper(Gene.Name1) %in% overlap_WT) %>%
      select(Chr, Start, End, Gene.Name1, Experiment.Insertions,
             TPH.Experiment, pvalue_adj.Reference),
    by = c("SYMBOL" = "Gene.Name1")
  )

write.csv(overlap_KD_df,
          file.path(output_dir, "overlap_FemaleCC_vs_FemaleKD_enriched.csv"),
          row.names = FALSE)

write.csv(overlap_WT_df,
          file.path(output_dir, "overlap_FemaleCC_vs_FemaleWT_enriched.csv"),
          row.names = FALSE)

cat("\nSaved overlap tables to:", output_dir, "\n")

# =============================================================================
# STEP 5: Gene sets for survival z-score analysis
# =============================================================================

geneset_KD_F <- overlap_KD  # CC peaks + upregulated when Egr1 KD → Egr1-repressed
geneset_WT_F <- overlap_WT  # CC peaks + downregulated when Egr1 KD → direct Egr1 targets

cat("\nGene set (CC + KD-enriched):", paste(geneset_KD_F, collapse = ", "), "\n")
cat("Gene set (CC + WT-enriched):", paste(geneset_WT_F, collapse = ", "), "\n")

