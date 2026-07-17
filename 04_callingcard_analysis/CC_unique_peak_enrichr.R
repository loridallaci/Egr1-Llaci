## =====================================================================
## Enrichr GO on the nearest genes of the male-unique / female-unique / shared
## CC peaks (Supp Fig 4). Mirrors enrichr_overlaps.R:
##   dbs = GO MF/CC/BP 2023 (+ WikiPathways mouse); plot GO_Biological_Process_2023.
##
## Gene sets = nearest gene (Gene Name1) of the PEAK-CENTRIC sets:
##   male-unique   8,265 peaks
##   female-unique 8,053 peaks
##   shared        3,712 peaks (male-side coords)
## Requires internet (Enrichr API).
## =====================================================================

suppressMessages({
  library(enrichR); library(GenomicRanges); library(ggplot2); library(dplyr); library(stringr)
})
setEnrichrSite("Enrichr")

cc_dir  <- "C:/Users/loril/Documents/Egr1/Egr1CC_vs_Egr1KDBulkRNA_FINAL_July2025"
out_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output/peak_overlap/enrichr"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dbs        <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
                "GO_Biological_Process_2023", "WikiPathways_2019_Mouse")
plot_db    <- "GO_Biological_Process_2023"

read_cc <- function(f) {
  d <- read.delim(file.path(cc_dir, f), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  d$gr <- GRanges(d$Chr, IRanges(d$Start, d$End)); d
}
M <- read_cc("Male_Egr1CC_peaks_072825.txt")
F <- read_cc("Female_Egr1CC_peaks_072825.txt")

genes_of <- function(df) { g <- unique(df$`Gene Name1`); sort(g[!is.na(g) & g != ""]) }
sets <- list(
  Male_unique   = genes_of(M[!overlapsAny(M$gr, F$gr), ]),
  Female_unique = genes_of(F[!overlapsAny(F$gr, M$gr), ]),
  Shared        = genes_of(M[ overlapsAny(M$gr, F$gr), ])
)
for (nm in names(sets)) cat(sprintf("%-14s %d nearest genes\n", nm, length(sets[[nm]])))

run_one <- function(genes, tag) {
  res <- tryCatch(enrichr(genes, dbs), error = function(e) { message("enrichr failed: ", e$message); NULL })
  if (is.null(res)) return(invisible())
  for (db in names(res)) {
    write.table(res[[db]], file.path(out_dir, sprintf("%s_enrichrResults_%s.txt", tag, db)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  ## ---- house-style enrichR plot (as in enrichr_overlaps.R / drug-screen) ----
  ## x = gene count, y = term, fill = Adjusted.P.value on a blue(low)->red(high)
  ## log10 gradient; top 10 terms by padj. Fonts bumped to the project presentation
  ## minimums (ticks >=14, axis titles >=16, title >=18, legend >=14).
  plt <- res[[plot_db]] %>%
    mutate(GeneCount = str_count(Genes, ";") + 1) %>%
    arrange(Adjusted.P.value) %>% slice_head(n = 10) %>%
    mutate(Term = str_wrap(Term, width = 40), Term = factor(Term, levels = rev(Term)))
  ph <- min(max(7, 0.9 * nrow(plt)), 16)
  p <- ggplot(plt, aes(x = GeneCount, y = Term, fill = Adjusted.P.value)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "blue", high = "red", name = "Adj. P-Value", trans = "log10") +
    theme_minimal(base_size = 14) +
    labs(title = sprintf("GO Enrichment - %s", sub("_", "-", tag)),
         x = "Number of Genes", y = "GO Term") +
    theme(axis.text.y  = element_text(size = 14, hjust = 1),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          legend.text  = element_text(size = 14),
          legend.title = element_text(size = 14),
          plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.margin  = margin(10, 10, 30, 10))
  ggsave(file.path(out_dir, sprintf("%s_enrichrResults_%s.pdf", tag, plot_db)), p, width = 10, height = ph, bg = "white")
  ggsave(file.path(out_dir, sprintf("%s_enrichrResults_%s.png", tag, plot_db)), p, width = 10, height = ph, dpi = 300, bg = "white")
  cat(sprintf("\n== %s : top GO-BP ==\n", tag))
  print(head(data.frame(Term = as.character(plt$Term), padj = signif(plt$Adjusted.P.value, 3),
                        genes = plt$Overlap), 8), row.names = FALSE)
}

for (nm in names(sets)) run_one(sets[[nm]], nm)
cat(sprintf("\nDone. Results in %s\n", out_dir))
