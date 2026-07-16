## Per-regulon-direction EnrichR figures in the SAME style as the rest of the paper
## (enrichR::plotEnrich, y="Count", coloured by P.value with side colourbar).
suppressMessages({library(enrichR); library(ggplot2)})
options(timeout = 120); setEnrichrSite("Enrichr")
DB  <- "GO_Biological_Process_2023"
IN  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/07_tcga_survival/output_regulon_survival"
OUT <- file.path(IN, "enrichr_regulon_directions")

titles <- c(MaleUnique_Down="Male regulon (Down after KD)",
            MaleUnique_Up  ="Male regulon (Up after KD)",
            FemaleUnique_Down="Female regulon (Down after KD)",
            FemaleUnique_Up  ="Female regulon (Up after KD)")

for (s in names(titles)) {
  genes <- toupper(readLines(file.path(IN, paste0("regulon_genes_", s, "_052726.txt"))))
  genes <- genes[nzchar(genes)]
  result <- enrichr(genes, DB)[[1]]
  p <- plotEnrich(result, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value") +
    ggtitle(paste0(titles[[s]], "\n", DB)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  outf <- file.path(OUT, paste0("enrichr_GO_BP2023_", s, ".pdf"))
  pdf(file = outf, width = 7, height = 7, onefile = TRUE, useDingbats = FALSE)
  print(p); invisible(dev.off())
  cat("wrote:", basename(outf), "\n")
}
