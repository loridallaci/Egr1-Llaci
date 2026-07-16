# =============================================================================
# Enrichr: DE between the female EGR1-motif-HIGH small cluster (c4) and the
# big EGR1-low clusters (c0/c1), in the multiome.
#
# Same workflow + databases as Enrichr_Egr1KD_bulk.R (enrichR + plotEnrich).
# DE input is the list used by the regulon panel, so the two are consistent:
#   05_integration/output/female_c4_vs_big_DE.csv   (rownames = gene;
#     cols: p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
#
# Two gene sets:
#   EGR1high_c4_up   = up in the EGR1-high small cluster (avg_log2FC >= 0.5)
#   EGR1low_big_up   = up in the big EGR1-low clusters   (avg_log2FC <= -0.5)
#
# Needs the Enrichr website (internet) - runs on the author's machine.
# =============================================================================
library(enrichR)
library(ggplot2)
library(dplyr)

websiteLive <- getOption("enrichR.live")
if (isTRUE(websiteLive)) setEnrichrSite("Enrichr")   # mouse symbols resolve here too

base <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/05_integration"
OutputDirectory <- file.path(base, "output/enrichment")
dir.create(OutputDirectory, showWarnings = FALSE, recursive = TRUE)

dbs <- c("GO_Molecular_Function_2023")

# presentation-legible fonts (>=14 ticks, >=16 axis titles, >=18 title, >=14 legend)
big_fonts <- theme(plot.title  = element_text(hjust = 0.5, size = 18, face = "bold"),
                   axis.title  = element_text(size = 16),
                   axis.text   = element_text(size = 14),
                   legend.title= element_text(size = 14),
                   legend.text = element_text(size = 14))

# --- DE lists ----------------------------------------------------------------
de <- read.csv(file.path(base, "output/female_c4_vs_big_DE.csv"), row.names = 1)
de$SYMBOL <- rownames(de)
up_c4  <- filter(de, avg_log2FC >=  0.5 & p_val_adj <= 0.05)   # EGR1-high small cluster
up_big <- filter(de, avg_log2FC <= -0.5 & p_val_adj <= 0.05)   # EGR1-low big clusters
cat("EGR1high_c4_up genes:", nrow(up_c4), " | EGR1low_big_up genes:", nrow(up_big), "\n")

run_enrichr <- function(genes, label, title_lab) {
  enriched <- enrichr(genes, dbs)
  for (db in dbs) {
    result <- enriched[[db]]
    if (is.null(result) || nrow(result) == 0) { cat("  skip", db, "(no results)\n"); next }
    write.table(result, file.path(OutputDirectory, paste0(label, "_", db, ".txt")),
                quote = FALSE, row.names = TRUE, sep = "\t")
    p <- plotEnrich(result, showTerms = 10, numChar = 45, y = "Count", orderBy = "P.value") +
      ggtitle(paste0(title_lab, "\n", db)) + big_fonts
    pdf(file.path(OutputDirectory, paste0(label, "_", db, ".pdf")),
        width = 12, height = 6.5, onefile = TRUE, useDingbats = FALSE)
    print(p); invisible(dev.off())
  }
}

run_enrichr(up_c4$SYMBOL,  "EGR1high_c4_up", "EGR1-high small cluster (c4) up")
run_enrichr(up_big$SYMBOL, "EGR1low_big_up", "EGR1-low big clusters (c0/c1) up")

message("\nDone. Enrichr results -> ", OutputDirectory)
