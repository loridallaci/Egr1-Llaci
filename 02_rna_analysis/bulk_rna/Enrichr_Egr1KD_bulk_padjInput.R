# =============================================================================
# Enrichr Analysis: Egr1 KD vs WT DE genes (Male and Female separately)
# INPUT DEG CUTOFF: adjusted P (padj) <= 0.05 & |log2FC| >= 0.5
#   (padj variant of Enrichr_Egr1KD_bulk.R, which used raw pvalue <= 0.05)
# Outputs mirror the raw-p set but land in .../enrichr_padjInput/ so both are kept.
# =============================================================================
suppressMessages({library(enrichR); library(ggplot2); library(dplyr)})
setEnrichrSite("Enrichr")

indir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_Egr1KD_gRNA3_vs_Neg1"
OutputDirectory <- file.path(indir, "enrichr_padjInput")
dir.create(OutputDirectory, showWarnings = FALSE, recursive = TRUE)
dbs <- c("WikiPathways_2019_Mouse", "GO_Biological_Process_2023", "GO_Molecular_Function_2023")

files <- c(Male   = "Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt",
           Female = "Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt")

run_one <- function(sex, dir_lab, tag, title_lab) {
  de <- read.csv(file.path(indir, files[sex]), sep = "\t")
  keep <- !is.na(de$padj) & de$padj <= 0.05                     # <-- padj cutoff
  if (dir_lab == "up")   genes <- de$SYMBOL[keep & de$log2FoldChange >=  0.5]
  if (dir_lab == "down") genes <- de$SYMBOL[keep & de$log2FoldChange <= -0.5]
  genes <- unique(na.omit(genes))
  cat(sprintf("%-6s %-14s %d genes\n", sex, tag, length(genes)))
  if (length(genes) < 2) { cat("  <2 genes, skipping\n"); return(invisible()) }
  enriched <- enrichr(genes, dbs)
  for (db in dbs) {
    result <- enriched[[db]]
    if (is.null(result) || nrow(result) == 0) { cat("  skip", db, "\n"); next }
    write.table(result, file.path(OutputDirectory, paste0(sex, "_", tag, "_", db, ".txt")),
                quote = FALSE, row.names = TRUE, sep = "\t")
    p <- plotEnrich(result, showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value") +
      ggtitle(paste0(sex, " ", title_lab, "\n", db, "  (padj<=0.05 & |log2FC|>=0.5)")) +
      theme(plot.title = element_text(hjust = 0.5, size = 10))
    tryCatch({ pdf(file.path(OutputDirectory, paste0(sex, "_", tag, "_", db, ".pdf")),
                   width = 7, height = 7, onefile = TRUE, useDingbats = FALSE); print(p); dev.off() },
             error = function(e) message("  PDF fail ", db, ": ", conditionMessage(e)))
  }
}

# KD-enriched = Egr1-repressed (up in KD);  WT-enriched = Egr1-activated (down in KD)
run_one("Male",   "up",   "Egr1KD_enriched", "Egr1 KD-enriched (Egr1-repressed)")
run_one("Male",   "down", "WT_enriched",     "WT-enriched (Egr1-activated)")
run_one("Female", "up",   "Egr1KD_enriched", "Egr1 KD-enriched (Egr1-repressed)")
run_one("Female", "down", "WT_enriched",     "WT-enriched (Egr1-activated)")

message("\nDone. padj-cutoff Enrichr results in: ", OutputDirectory)
