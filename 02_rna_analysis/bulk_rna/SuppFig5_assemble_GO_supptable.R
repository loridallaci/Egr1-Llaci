## Assemble Supplementary Data workbook: GO Biological Process enrichment for the
## four Supp Fig 5 sex-DE gene sets (control/KD x female-/male-higher), one tab each.
## DE tables are kept separate as Source Data (not in this workbook).
suppressMessages(library(openxlsx))

dir_in <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_sexComparison/enrichr"
out    <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_sexComparison/SupplementaryData_SuppFig5_GO_enrichment.xlsx"

sets <- c(
  Control_Female = "Control_FemaleHigher",
  Control_Male   = "Control_MaleHigher",
  KD_Female      = "Egr1KD_FemaleHigher",
  KD_Male        = "Egr1KD_MaleHigher")

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", halign = "center", border = "bottom", fgFill = "#D9E1F2")
for (tab in names(sets)) {
  f <- file.path(dir_in, sprintf("%s_enrichrResults_GO_Biological_Process_2023.txt", sets[[tab]]))
  d <- read.delim(f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  d$N_genes <- lengths(strsplit(d$Genes, ";"))
  keep <- c("Term","Overlap","N_genes","P.value","Adjusted.P.value","Odds.Ratio","Combined.Score","Genes")
  d <- d[order(d$Adjusted.P.value), keep]
  addWorksheet(wb, tab)
  writeData(wb, tab, d, headerStyle = hdr)
  freezePane(wb, tab, firstRow = TRUE)
  setColWidths(wb, tab, cols = 1:ncol(d), widths = c(46,10,9,12,16,11,14,60))
  cat(sprintf("%-15s %d terms (%d at padj<=0.05)\n", tab, nrow(d), sum(d$Adjusted.P.value <= 0.05)))
}
saveWorkbook(wb, out, overwrite = TRUE)
cat("\nWrote:", out, "\n")
