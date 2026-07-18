## Supplementary Table 7 : GO Biological Process enrichment of the sex-biased
## genes, one combined sheet (paper-standard). Significant terms only (padj<=0.05),
## comparison/direction as COLUMNS. DE tables kept separate as Source Data.
suppressMessages(library(openxlsx))

dir_in <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_sexComparison/enrichr"
out    <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_sexComparison/SupplementaryTable7_GO_enrichment_sexDE.xlsx"

sets <- list(
  c("Control_FemaleHigher","Control","Female-higher"),
  c("Control_MaleHigher",  "Control","Male-higher"),
  c("Egr1KD_FemaleHigher", "Egr1-KD","Female-higher"),
  c("Egr1KD_MaleHigher",   "Egr1-KD","Male-higher"))

rows <- list()
for (s in sets) {
  f <- file.path(dir_in, sprintf("%s_enrichrResults_GO_Biological_Process_2023.txt", s[1]))
  d <- read.delim(f, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
  d <- d[d$Adjusted.P.value <= 0.05, ]
  d <- d[order(d$Adjusted.P.value), ]
  d$N_genes <- lengths(strsplit(d$Genes, ";"))
  d <- data.frame(Comparison = s[2], Direction = s[3],
                  d[, c("Term","Overlap","N_genes","P.value","Adjusted.P.value",
                        "Odds.Ratio","Combined.Score","Genes")],
                  check.names = FALSE)
  cat(sprintf("%-8s %-14s %d significant terms\n", s[2], s[3], nrow(d)))
  rows[[length(rows)+1]] <- d
}
tab <- do.call(rbind, rows)

wb <- createWorkbook()
addWorksheet(wb, "Supplementary Table 7")
title <- "Supplementary Table 7. GO Biological Process enrichment of genes differentially expressed between male and female astrocytes (control and Egr1-KD cells; enrichR, adjusted P <= 0.05)."
writeData(wb, 1, title, startRow = 1)
mergeCells(wb, 1, cols = 1:10, rows = 1)
addStyle(wb, 1, createStyle(textDecoration="bold", wrapText=TRUE), rows=1, cols=1)
writeData(wb, 1, tab, startRow = 3,
          headerStyle = createStyle(textDecoration="bold", halign="center", border="bottom", fgFill="#D9E1F2"))
freezePane(wb, 1, firstActiveRow = 4)
setColWidths(wb, 1, cols=1:10, widths=c(11,14,46,10,9,12,16,11,14,60))
saveWorkbook(wb, out, overwrite = TRUE)
cat(sprintf("\nTotal rows: %d\nWrote: %s\n", nrow(tab), out))
