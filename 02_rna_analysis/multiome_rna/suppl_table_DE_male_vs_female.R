## Supplementary table: multiome RNA differential expression, Male vs Female.
## Reads the full DE result written by gene_de.R and writes a tidy table of
## ALL tested genes, with a significance flag (so reviewers can see everything
## that was tested; the paper text cites the significant count).
##
## Test: Seurat FindMarkers, Wilcoxon; ident.1 = Male, so avg_log2FC > 0 = higher
##   in male. A gene is flagged Significant when p_val_adj <= 0.05 AND
##   |avg_log2FC| >= 0.5 -- the same cutoffs used to define the male-higher /
##   female-higher DE sets (adjusted p, consistent with the rest of the paper).
base   <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci"
in_csv <- file.path(base, "data output", "DE_male_vs_female_allcells_allgenes.csv")
out_csv<- file.path(base, "02_rna_analysis", "multiome_rna", "SupplTable_multiome_DE_male_vs_female.csv")

de <- read.csv(in_csv, stringsAsFactors = FALSE, check.names = FALSE)
is_sig <- de$p_val_adj <= 0.05 & abs(de$avg_log2FC) >= 0.5
de$Significant <- ifelse(is_sig, "Yes", "No")
de$Direction   <- ifelse(!is_sig, "NS",
                    ifelse(de$avg_log2FC > 0, "Male-higher", "Female-higher"))

tab <- data.frame(
  Gene        = de$gene,
  avg_log2FC  = round(de$avg_log2FC, 4),
  p_val       = de$p_val,
  p_val_adj   = de$p_val_adj,
  pct_Male    = de$pct.1,   # fraction of male cells expressing
  pct_Female  = de$pct.2,   # fraction of female cells expressing
  Significant = de$Significant,
  Direction   = de$Direction,
  stringsAsFactors = FALSE)
## significant genes first (male-higher desc, then female-higher), then the rest
tab <- tab[order(tab$Significant != "Yes", -tab$avg_log2FC), ]

write.csv(tab, out_csv, row.names = FALSE)
cat("Wrote:", out_csv, "\n")
cat("Total tested:", nrow(tab),
    "| Significant:", sum(tab$Significant=="Yes"),
    "(Male-higher:", sum(tab$Direction=="Male-higher"),
    ", Female-higher:", sum(tab$Direction=="Female-higher"), ")\n")
