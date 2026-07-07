## =====================================================================
## STEP 2 of the chain:  DE  +  CallingCard(20kb) intersection  ->  regulon summary
##
## Produces the two summary CSVs that Regulon_TCGA_permutation_052726.R reads:
##   Egr1CC_allGenes_DE_summary_with_GO_MALES_gRNA3_20kb_rerun_052726.csv
##   Egr1CC_allGenes_DE_summary_with_GO_FEMALES_gRNA3_20kb_rerun_052726.csv
##
## Minimal version: keeps only the columns the permutation needs
##   (SYMBOL + kd3_male_sig / kd3_female_sig).
## The full Rmd version additionally annotates Brd4/JQ1/GO/biotype — not needed
## for the survival permutation, so omitted here.
##
## FULL CHAIN:
##   1) Egr1KD_gRNA3_vs_Neg1_DE_CONSOLIDATED.R   -> the two DE tables
##   2) THIS SCRIPT                              -> the two regulon summaries
##   3) Regulon_TCGA_permutation_052726.R        -> gene sets + Cox + permutation (panel f)
## =====================================================================

library(dplyr)
library(tibble)
library(readr)

## ---- Paths ----------------------------------------------------------
out_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
de_dir  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_Egr1KD_gRNA3_vs_Neg1"

## CallingCard 20kb nearest-gene  x  gRNA3-vs-Neg1 DE  (the CC∩DE intersection, from step 1b)
cc_male_file   <- file.path(out_dir, "overlap_MaleCC_20kb_NearestGene_vs_Egr1g3vsNeg1.csv")
cc_female_file <- file.path(out_dir, "overlap_FemaleCC_20kb_NearestGene_vs_Egr1g3vsNeg1.csv")

## DE tables (output of step 1)
de_male_file   <- file.path(de_dir, "Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt")
de_female_file <- file.path(de_dir, "Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt")

## DE thresholds for the direction call (same as the Rmd)
lfc_cut <- 0.5
p_cut   <- 0.05

## ---- Helper: build one sex's regulon summary ------------------------
build_summary <- function(cc_file, de_file, sig_colname) {
  ## regulon genes = nearest gene of each CC peak within 20kb (the 'Gene' column)
  cc <- read.delim(cc_file, sep = "\t", stringsAsFactors = FALSE)
  regulon_genes <- unique(na.omit(cc$Gene))

  ## DE table (has SYMBOL, log2FoldChange, pvalue, padj)
  de <- read.delim(de_file, sep = "\t", stringsAsFactors = FALSE) %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::select(SYMBOL, log2FoldChange, pvalue, padj)

  tibble(SYMBOL = regulon_genes) %>%
    dplyr::left_join(de, by = "SYMBOL") %>%
    dplyr::mutate(
      !!sig_colname := dplyr::case_when(
        !is.na(log2FoldChange) & log2FoldChange >=  lfc_cut & pvalue <= p_cut ~ "Up after KD",
        !is.na(log2FoldChange) & log2FoldChange <= -lfc_cut & pvalue <= p_cut ~ "Down after KD",
        TRUE ~ "Not DE"
      )
    )
}

## ---- Build + write both summaries -----------------------------------
male_summary   <- build_summary(cc_male_file,   de_male_file,   "kd3_male_sig")
female_summary <- build_summary(cc_female_file, de_female_file, "kd3_female_sig")

cat("Male regulon genes:",   nrow(male_summary),
    "| DE-sig:",  sum(male_summary$kd3_male_sig   != "Not DE"), "\n")
cat("Female regulon genes:", nrow(female_summary),
    "| DE-sig:",  sum(female_summary$kd3_female_sig != "Not DE"), "\n")

write_csv(male_summary,
  file.path(out_dir, "Egr1CC_allGenes_DE_summary_with_GO_MALES_gRNA3_20kb_rerun_052726.csv"))
write_csv(female_summary,
  file.path(out_dir, "Egr1CC_allGenes_DE_summary_with_GO_FEMALES_gRNA3_20kb_rerun_052726.csv"))

message("Done. Now run Regulon_TCGA_permutation_052726.R for the survival + permutation (panel f).")
