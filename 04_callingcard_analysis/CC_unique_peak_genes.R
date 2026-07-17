## =====================================================================
## Nearest genes of the male-unique and female-unique CC peaks (Supp Fig 4b)
##
## "Unique" = peak-centric: whole peaks from one sex that do NOT overlap any peak
## in the other sex (== bedtools intersect -v). These retain the caller's nearest-
## gene annotation (Gene Name1 / Distance1), unlike the merged region-centric sets.
##   male-unique  = 8,265 peaks
##   female-unique = 8,053 peaks
##
## Outputs (to output/peak_overlap/):
##   CC_male_unique_peaks_withGenes.txt    - full peak table + nearest gene/dist
##   CC_female_unique_peaks_withGenes.txt
##   CC_male_unique_geneList.txt           - unique gene symbols (one per line)
##   CC_female_unique_geneList.txt
## =====================================================================

suppressMessages(library(GenomicRanges))

cc_dir  <- "C:/Users/loril/Documents/Egr1/Egr1CC_vs_Egr1KDBulkRNA_FINAL_July2025"
out_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output/peak_overlap"

read_cc <- function(f) {
  d <- read.delim(file.path(cc_dir, f), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  d$gr <- GRanges(d$Chr, IRanges(d$Start, d$End))
  d
}
M <- read_cc("Male_Egr1CC_peaks_072825.txt")
F <- read_cc("Female_Egr1CC_peaks_072825.txt")

## peak-centric unique = does NOT overlap the other sex
M_uni <- M[!overlapsAny(M$gr, F$gr), ]
F_uni <- F[!overlapsAny(F$gr, M$gr), ]
cat(sprintf("male-unique %d peaks   female-unique %d peaks\n", nrow(M_uni), nrow(F_uni)))

## full peak tables (drop the helper GRanges column)
write.table(M_uni[, setdiff(names(M_uni), "gr")],
            file.path(out_dir, "CC_male_unique_peaks_withGenes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(F_uni[, setdiff(names(F_uni), "gr")],
            file.path(out_dir, "CC_female_unique_peaks_withGenes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## unique gene symbol lists (nearest gene = "Gene Name1"), blanks/NA dropped
gene_list <- function(df) {
  g <- unique(df$`Gene Name1`)
  g <- g[!is.na(g) & g != ""]
  sort(g)
}
gm <- gene_list(M_uni)
gf <- gene_list(F_uni)
writeLines(gm, file.path(out_dir, "CC_male_unique_geneList.txt"))
writeLines(gf, file.path(out_dir, "CC_female_unique_geneList.txt"))

cat(sprintf("\nunique nearest genes:  male %d   female %d\n", length(gm), length(gf)))
cat(sprintf("shared between the two unique lists: %d\n", length(intersect(gm, gf))))
cat("\n(genes appear in both lists when a male-unique and a female-unique PEAK",
    "\n sit near the same gene - the peaks differ, the nearest gene coincides.)\n")

cat("\n--- first 25 male-unique nearest genes ---\n");  print(head(gm, 25))
cat("\n--- first 25 female-unique nearest genes ---\n"); print(head(gf, 25))

## ---- targeted lookup: Ptk2b / Frzb / Nab1 across all three sets -------------
## shared set = peaks present in BOTH sexes (peak-centric, male-side coords)
M_shr <- M[overlapsAny(M$gr, F$gr), ]
gene_hits <- function(df, gene) {
  h <- df[!is.na(df$`Gene Name1`) & df$`Gene Name1` == gene, ]
  if (nrow(h) == 0) return(NULL)
  data.frame(peak = sprintf("%s:%d-%d", h$Chr, h$Start, h$End),
             dist = h$Distance1, stringsAsFactors = FALSE)
}
cat("\n\n================ Ptk2b / Frzb / Nab1 lookup ================\n")
for (g in c("Ptk2b", "Frzb", "Nab1")) {
  cat(sprintf("\n%s:\n", g))
  for (setnm in c("male-unique", "female-unique", "shared")) {
    df <- switch(setnm, `male-unique` = M_uni, `female-unique` = F_uni, shared = M_shr)
    h <- gene_hits(df, g)
    if (is.null(h)) cat(sprintf("   %-14s -- absent\n", setnm))
    else for (i in seq_len(nrow(h)))
      cat(sprintf("   %-14s %s  (nearest-gene dist %s bp)\n", setnm, h$peak[i], format(h$dist[i], big.mark=",")))
  }
}
