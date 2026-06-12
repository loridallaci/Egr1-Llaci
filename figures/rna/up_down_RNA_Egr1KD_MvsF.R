# =============================================================================
# Overlap of DE genes between Male and Female Egr1 KD
# Genes that go UP in both, and genes that go DOWN in both
# Males LEFT | Females RIGHT
# =============================================================================

library(dplyr)
library(ggVennDiagram)
library(ggplot2)

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1_bulk_RNAseq/Egr1_KD_030525/diffexp_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers"
de_dir <- "data"
# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/enrichment/Egr1KD"
output_dir <- "output/enrichment/Egr1KD"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD DE DATA
# =============================================================================

de_male <- read.csv(
  file.path(de_dir, "Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt"),
  sep = "\t"
)

de_female <- read.csv(
  file.path(de_dir, "Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt"),
  sep = "\t"
)

# =============================================================================
# SPLIT BY DIRECTION
# =============================================================================

male_up     <- filter(de_male,   log2FoldChange >= 0.5  & pvalue <= 0.05)
female_up   <- filter(de_female, log2FoldChange >= 0.5  & pvalue <= 0.05)
male_down   <- filter(de_male,   log2FoldChange <= -0.5 & pvalue <= 0.05)
female_down <- filter(de_female, log2FoldChange <= -0.5 & pvalue <= 0.05)

cat("Male UP:    ", nrow(male_up),   "\n")
cat("Female UP:  ", nrow(female_up), "\n")
cat("Male DOWN:  ", nrow(male_down), "\n")
cat("Female DOWN:", nrow(female_down), "\n")

# =============================================================================
# OVERLAPS
# =============================================================================

overlap_up   <- intersect(toupper(male_up$SYMBOL),   toupper(female_up$SYMBOL))
overlap_down <- intersect(toupper(male_down$SYMBOL), toupper(female_down$SYMBOL))

cat("\nGenes UP in both Male and Female KD:  ", length(overlap_up),   "\n")
print(overlap_up)

cat("\nGenes DOWN in both Male and Female KD:", length(overlap_down), "\n")
print(overlap_down)

# =============================================================================
# SAVE OVERLAP TABLES
# =============================================================================

overlap_up_df <- de_male %>%
  filter(toupper(SYMBOL) %in% overlap_up) %>%
  select(SYMBOL, log2FoldChange, pvalue, padj) %>%
  rename(log2FC_Male = log2FoldChange,
         pvalue_Male = pvalue,
         padj_Male = padj) %>%
  left_join(
    de_female %>%
      filter(toupper(SYMBOL) %in% overlap_up) %>%
      select(SYMBOL, log2FoldChange, pvalue, padj) %>%
      rename(log2FC_Female = log2FoldChange,
             pvalue_Female = pvalue,
             padj_Female = padj),
    by = "SYMBOL"
  ) %>%
  arrange(pvalue_Male)

overlap_down_df <- de_male %>%
  filter(toupper(SYMBOL) %in% overlap_down) %>%
  select(SYMBOL, log2FoldChange, pvalue, padj) %>%
  rename(log2FC_Male = log2FoldChange,
         pvalue_Male = pvalue,
         padj_Male = padj) %>%
  left_join(
    de_female %>%
      filter(toupper(SYMBOL) %in% overlap_down) %>%
      select(SYMBOL, log2FoldChange, pvalue, padj) %>%
      rename(log2FC_Female = log2FoldChange,
             pvalue_Female = pvalue,
             padj_Female = padj),
    by = "SYMBOL"
  ) %>%
  arrange(pvalue_Male)

write.csv(overlap_up_df,
          file.path(output_dir, "overlap_MaleFemale_KD_UP_inBoth.csv"),
          row.names = FALSE)

write.csv(overlap_down_df,
          file.path(output_dir, "overlap_MaleFemale_KD_DOWN_inBoth.csv"),
          row.names = FALSE)

cat("\nSaved overlap tables to:", output_dir, "\n")

# =============================================================================
# VENN PLOT FUNCTION (FIXED — deterministic Male LEFT / Female RIGHT)
# =============================================================================

library(ggplot2)
library(ggforce)

plot_venn_fixed <- function(set_male, set_female,
                            title_text) {
  
  male   <- unique(toupper(na.omit(set_male)))
  female <- unique(toupper(na.omit(set_female)))
  
  overlap <- intersect(male, female)
  only_m  <- setdiff(male, female)
  only_f  <- setdiff(female, male)
  
  df <- data.frame(
    x = c(-1, 1, 0),
    y = c(0, 0, 0),
    group = c("Male", "Female", "Overlap"),
    size = c(length(only_m), length(only_f), length(overlap))
  )
  
  ggplot() +
    
    # circles
    geom_circle(aes(x0 = -1, y0 = 0, r = 1.2),
                color = "#4DA6FF", fill = "#4DA6FF", alpha = 0.3) +
    geom_circle(aes(x0 =  1, y0 = 0, r = 1.2),
                color = "#FF69B4", fill = "#FF69B4", alpha = 0.3) +
    
    # labels
    annotate("text", x = -1, y = 1.4, label = "Male", fontface = "bold") +
    annotate("text", x =  1, y = 1.4, label = "Female", fontface = "bold") +
    
    annotate("text", x = -1.5, y = 0, label = length(only_m), size = 6) +
    annotate("text", x =  1.5, y = 0, label = length(only_f), size = 6) +
    annotate("text", x =  0, y = 0, label = length(overlap), size = 6, fontface = "bold") +
    
    coord_fixed() +
    xlim(-2.5, 2.5) +
    ylim(-2, 2) +
    
    ggtitle(title_text) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

# =============================================================================
# PLOT 1: UP in KD
# =============================================================================

pdf(file.path(output_dir, "venn_MaleFemale_Egr1KDBulk_KD_UP.pdf"), width = 6, height = 6)

print(
  plot_venn_fixed(male_up$SYMBOL,
                  female_up$SYMBOL,
                  "Upregulated after Egr1 KD")
)

dev.off()

cat("Saved: venn_MaleFemale_Egr1KDBulk_KD_UP.pdf\n")

# =============================================================================
# PLOT 2: DOWN in KD
# =============================================================================

pdf(file.path(output_dir, "venn_MaleFemale_Egr1KDBulk_KD_DOWN.pdf"),
    width = 6, height = 6)

print(
  plot_venn_fixed(male_down$SYMBOL,
                  female_down$SYMBOL,
                  "Downregulated after Egr1 KD")
)

dev.off()

cat("Saved: venn_MaleFemale_Egr1KDBulk_KD_DOWN.pdf\n")

# =============================================================================
# FISHER TESTS — overlap significance
# =============================================================================

universe_conservative <- intersect(
  toupper(filter(de_male,   pvalue <= 0.1)$SYMBOL),
  toupper(filter(de_female, pvalue <= 0.1)$SYMBOL)
)

cat("\nConservative universe:", length(universe_conservative), "genes\n")

run_fisher <- function(setA, setB, universe, labelA, labelB) {
  
  setA <- intersect(setA, universe)
  setB <- intersect(setB, universe)
  
  overlap <- length(intersect(setA, setB))
  a_only  <- length(setdiff(setA, setB))
  b_only  <- length(setdiff(setB, setA))
  neither <- length(universe) - overlap - a_only - b_only
  
  mat <- matrix(c(overlap, a_only, b_only, neither), nrow = 2)
  test <- fisher.test(mat, alternative = "greater")
  
  cat("\n---", labelA, "vs", labelB, "---\n")
  cat("Overlap:", overlap,
      "| OR:", round(test$estimate, 3),
      "| p =", format(test$p.value, scientific = TRUE, digits = 3), "\n")
  
  invisible(list(overlap = overlap,
                 odds_ratio = test$estimate,
                 p_value = test$p.value))
}

res_up   <- run_fisher(toupper(male_up$SYMBOL),
                       toupper(female_up$SYMBOL),
                       universe_conservative,
                       "Male UP", "Female UP")

res_down <- run_fisher(toupper(male_down$SYMBOL),
                       toupper(female_down$SYMBOL),
                       universe_conservative,
                       "Male DOWN", "Female DOWN")

cat("\nAll done. Files saved to:", output_dir, "\n")

library(enrichR)
library(dplyr)

dbs <- c("GO_Biological_Process_2023")

run_enrichr <- function(genes, label, output_dir) {
  
  genes <- unique(na.omit(toupper(genes)))
  
  if (length(genes) < 5) {
    cat("Skipping", label, "- too few genes\n")
    return(NULL)
  }
  
  cat("\nRunning Enrichr:", label, "(", length(genes), "genes)\n")
  
  enriched <- enrichr(genes, dbs)
  
  result <- enriched[[dbs[1]]]   # ✔ FIXED
  
  if (is.null(result) || nrow(result) == 0) {
    cat("No enrichment for", label, "\n")
    return(NULL)
  }
  
  # Save table
  write.table(
    result,
    file = file.path(output_dir, paste0(label, "_GO_BP_2023.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Plot
  p <- plotEnrich(result,
                  showTerms = 20,
                  numChar = 50,
                  y = "Count",
                  orderBy = "P.value") +
    ggtitle(label) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10)
    )
  
  pdf(file.path(output_dir, paste0(label, "_GO_BP_2023.pdf")),
      width = 9, height = 7)
  print(p)
  dev.off()
  
  return(result)
}

# =========================
# UP sets (from Venn)
# =========================

res_up_shared <- run_enrichr(
  up_shared,
  "UP_Shared_MaleFemale",
  output_dir
)

res_up_male <- run_enrichr(
  up_male_only,
  "UP_MaleOnly",
  output_dir
)

res_up_female <- run_enrichr(
  up_female_only,
  "UP_FemaleOnly",
  output_dir
)

# =========================
# DOWN sets (from Venn)
# =========================

res_down_shared <- run_enrichr(
  down_shared,
  "DOWN_Shared_MaleFemale",
  output_dir
)

res_down_male <- run_enrichr(
  down_male_only,
  "DOWN_MaleOnly",
  output_dir
)

res_down_female <- run_enrichr(
  down_female_only,
  "DOWN_FemaleOnly",
  output_dir
)

cat("\n✔ Enrichr complete using Venn-derived gene sets\n")


