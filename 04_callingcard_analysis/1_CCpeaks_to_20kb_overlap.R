## =====================================================================
## STEP "1b" of the chain:  raw CallingCard peaks  ->  20kb nearest-gene
##                          ->  intersect with DE  ->  overlap gene list
##
## Produces the CC∩DE overlap files that Egr1_CCxDE_to_regulonSummary_CONSOLIDATED.R reads:
##   overlap_MaleCC_20kb_NearestGene_vs_Egr1g3vsNeg1.csv
##   overlap_FemaleCC_20kb_NearestGene_vs_Egr1g3vsNeg1.csv
## ...and the 20kb-filtered peak tables:
##   {Male,Female}_Egr1CC_peaks_20kbThreshhold_090825.txt
##
## Note: the CallingCard peak caller already annotates each peak with its
## nearest RefSeq gene + distance ("Gene Name1" / "Distance1"), so the
## "20kb nearest gene" step is simply |Distance1| <= 20000. The overlap is
## the intersection of those genes with the DE-significant gene symbols.
##
## FULL CHAIN:
##   1)  Egr1KD_gRNA3_vs_Neg1_DE_CONSOLIDATED.R       -> DE tables
##   1b) THIS SCRIPT                                  -> 20kb peaks + CC∩DE overlap
##   2)  Egr1_CCxDE_to_regulonSummary_CONSOLIDATED.R  -> regulon summaries
##   3)  Regulon_TCGA_permutation_052726.R            -> survival + permutation (panel f)
## =====================================================================

library(dplyr)
library(readr)

## ---- Paths ----------------------------------------------------------
cc_dir <- "C:/Users/loril/Documents/Egr1/Egr1CC_vs_Egr1KDBulkRNA_FINAL_July2025"   # raw CC peaks (input)
de_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_Egr1KD_gRNA3_vs_Neg1"

## Outputs -> calling-card stage folder (git)
out_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cc_raw_male   <- file.path(cc_dir, "Male_Egr1CC_peaks_072825.txt")   # raw peaks = input (read only)
cc_raw_female <- file.path(cc_dir, "Female_Egr1CC_peaks_072825.txt")

de_male   <- file.path(de_dir, "Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt")
de_female <- file.path(de_dir, "Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt")

## ---- Parameters -----------------------------------------------------
tss_window <- 20000     # 20 kb nearest-gene window
lfc_cut    <- 0.5       # DE significance (matches the regulon summary step)
p_cut      <- 0.05

## ---- Helper ---------------------------------------------------------
make_overlap <- function(cc_file, de_file, peaks20_out, overlap_out) {
  ## raw CC peaks (already have 'Gene Name1' + 'Distance1' from the caller)
  cc <- read.delim(cc_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

  ## 20kb nearest-gene filter
  cc_20kb <- cc[!is.na(cc$Distance1) & abs(cc$Distance1) <= tss_window, ]
  write.table(cc_20kb, peaks20_out, sep = "\t", quote = FALSE, row.names = FALSE)

  cc_genes <- unique(na.omit(cc_20kb$`Gene Name1`))
  cc_genes <- cc_genes[cc_genes != ""]

  ## DE-significant gene symbols (gRNA3 vs Neg1)
  de <- read.delim(de_file, sep = "\t", stringsAsFactors = FALSE)
  de_sig <- de$SYMBOL[!is.na(de$log2FoldChange) & !is.na(de$pvalue) &
                        abs(de$log2FoldChange) >= lfc_cut & de$pvalue <= p_cut]
  de_sig <- unique(na.omit(de_sig))

  ## intersection = CC∩DE regulon gene list
  overlap <- intersect(cc_genes, de_sig)
  write_csv(data.frame(Gene = overlap), overlap_out)

  cat(basename(overlap_out), ":", length(cc_genes), "CC(20kb) genes  ∩ ",
      length(de_sig), "DE-sig  =>", length(overlap), "overlap\n")
  invisible(overlap)
}

## ---- Run for both sexes ---------------------------------------------
make_overlap(
  cc_raw_male, de_male,
  file.path(out_dir, "Male_Egr1CC_peaks_20kbThreshhold_090825.txt"),
  file.path(out_dir, "overlap_MaleCC_20kb_NearestGene_vs_Egr1g3vsNeg1.csv"))

make_overlap(
  cc_raw_female, de_female,
  file.path(out_dir, "Female_Egr1CC_peaks_20kbThreshhold_090825.txt"),
  file.path(out_dir, "overlap_FemaleCC_20kb_NearestGene_vs_Egr1g3vsNeg1.csv"))

message("Done. Next run Egr1_CCxDE_to_regulonSummary_CONSOLIDATED.R")
