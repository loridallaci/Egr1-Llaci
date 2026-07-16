## EnrichR (GO Biological Process 2023) on the 4 directional UNIQUE regulons.
suppressMessages(library(enrichR))
options(timeout = 120)
setEnrichrSite("Enrichr")
DB  <- "GO_Biological_Process_2023"
IN  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/07_tcga_survival/output_regulon_survival"
OUT <- file.path(IN, "enrichr_regulon_directions")
dir.create(OUT, showWarnings = FALSE)

sets <- c("MaleUnique_Down","MaleUnique_Up","FemaleUnique_Down","FemaleUnique_Up")
for (s in sets) {
  genes <- toupper(readLines(file.path(IN, paste0("regulon_genes_", s, "_052726.txt"))))
  genes <- genes[nzchar(genes)]
  res <- tryCatch(enrichr(genes, DB)[[1]], error = function(e){message("FAIL ",s,": ",conditionMessage(e)); NULL})
  if (is.null(res) || !nrow(res)) { cat(sprintf("\n### %s (%d genes): no result\n", s, length(genes))); next }
  res <- res[order(res$Adjusted.P.value), ]
  outf <- file.path(OUT, paste0("enrichr_GO_BP2023_", s, ".csv"))
  write.csv(res, outf, row.names = FALSE)
  cat(sprintf("\n### %s  (%d genes) -> %d terms; top by adj.P:\n", s, length(genes), nrow(res)))
  top <- head(res[, c("Term","Overlap","P.value","Adjusted.P.value")], 8)
  for (i in seq_len(nrow(top))) cat(sprintf("  %-55s  %s  adjP=%.2e\n",
      substr(top$Term[i],1,55), top$Overlap[i], top$Adjusted.P.value[i]))
}
cat("\nWrote tables to:", OUT, "\n")
