## Supplementary table: ATAC differentially accessible regions (DARs), Male vs Female.
## Reads the annotated DAR result (from DARs.R / Figure_1C_DAR_volcano.R) and writes
## a tidy table of ALL tested peaks with a significance flag, so reviewers can see
## everything tested; the paper text cites the significant count.
##
## Test: Seurat FindMarkers on the peak assay, Wilcoxon; ident.1 = male, ident.2 =
##   female, so avg_log2FC > 0 = MORE accessible in male. A peak is flagged
##   Significant when p_val_adj <= 0.05 AND |avg_log2FC| >= 0.5 -- the same cutoffs
##   used for the Fig 1C volcano (Male-enriched / Female-enriched), and matching
##   the RNA DE table (adjusted p, consistent across the paper).
base   <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci"
in_csv <- file.path(base, "data output", "lot6_DAR_peaks_annotated.csv")
out_csv<- file.path(base, "03_atac_analysis", "SupplTable_ATAC_DAR_male_vs_female.csv")

da <- read.csv(in_csv, stringsAsFactors = FALSE, check.names = FALSE)
is_sig <- da$p_val_adj <= 0.05 & abs(da$avg_log2FC) >= 0.5
da$Significant <- ifelse(is_sig, "Yes", "No")
da$Direction   <- ifelse(!is_sig, "NS",
                    ifelse(da$avg_log2FC > 0, "Male-enriched", "Female-enriched"))

tab <- data.frame(
  Peak                     = da$X,          # chr-start-end
  avg_log2FC               = round(da$avg_log2FC, 4),
  p_val                    = da$p_val,
  p_val_adj                = da$p_val_adj,
  pct_Male                 = da$pct.1,       # fraction of male cells with the peak
  pct_Female               = da$pct.2,       # fraction of female cells with the peak
  nearest_gene             = da$nearest_gene,
  distance_to_nearest_gene = da$distance_to_nearest_gene,
  Significant              = da$Significant,
  Direction                = da$Direction,
  stringsAsFactors = FALSE)
## significant peaks first (male-enriched desc, then female-enriched), then the rest
tab <- tab[order(tab$Significant != "Yes", -tab$avg_log2FC), ]

write.csv(tab, out_csv, row.names = FALSE)
cat("Wrote:", out_csv, "\n")
cat("Total peaks:", nrow(tab),
    "| Significant:", sum(tab$Significant=="Yes"),
    "(Male-enriched:", sum(tab$Direction=="Male-enriched"),
    ", Female-enriched:", sum(tab$Direction=="Female-enriched"), ")\n")
