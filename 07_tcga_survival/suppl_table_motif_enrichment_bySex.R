## Supplementary table: TF motif enrichment in sex-biased cis-regulatory elements.
## Combines the male and female Signac FindMotifs outputs into ONE labelled table
## of ALL motifs tested in each sex, with an enrichment flag (so reviewers can see
## everything that was tested; the paper text cites the enriched counts).
##
## Test: Signac FindMotifs, hypergeometric, against a GC-content-matched background
##   of accessible peaks. Motifs are called Enriched at Benjamini-Hochberg
##   adjusted P <= 0.05 -- reproducing the 567 (male) / 350 (female) in the text.
##   Rank_within_sex is on RAW p (BH clamps ties), consistent with the rest of the
##   paper; the Enriched call still uses the adjusted p.
##
## NOTE: the same 746 motifs are tested in both sexes, so the two blocks are
##   directly comparable. A motif can be enriched in BOTH sexes (e.g. EGR1:
##   male rank 12, 2.18-fold; female rank 66, 1.67-fold) -- the "male/female"
##   TF lists in the paper come from the top 30 per sex, not from exclusivity.
base    <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci"
in_dir  <- file.path(base, "07_tcga_survival", "data_motifs")
out_dir <- file.path(base, "07_tcga_survival", "output_motif_enrichment")
out_csv <- file.path(out_dir, "SupplTable_RENIN_motif_enrichment_bySex.csv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_one <- function(file, sex) {
  d <- read.csv(file.path(in_dir, file), stringsAsFactors = FALSE, row.names = 1)
  data.frame(
    Sex                = sex,
    Motif_ID           = d$motif,
    TF                 = d$motif.name,
    Observed           = d$observed,
    Background         = d$background,
    Percent_observed   = round(d$percent.observed, 3),
    Percent_background = round(d$percent.background, 3),
    Fold_enrichment    = round(d$fold.enrichment, 4),
    P_value            = d$pvalue,
    P_adj_BH           = d$p.adjust,
    stringsAsFactors   = FALSE)
}

rank_and_flag <- function(d) {
  d <- d[order(d$P_value, -d$Fold_enrichment), ]
  d$Rank_within_sex <- seq_len(nrow(d))
  d$Enriched <- ifelse(d$P_adj_BH <= 0.05, "Yes", "No")
  d
}

m <- rank_and_flag(read_one("M_all_motifs_updated.csv", "Male"))
f <- rank_and_flag(read_one("F_all_motifs_updated.csv", "Female"))

tab <- rbind(m, f)[, c("Sex", "Rank_within_sex", "Motif_ID", "TF", "Observed",
                       "Background", "Percent_observed", "Percent_background",
                       "Fold_enrichment", "P_value", "P_adj_BH", "Enriched")]
write.csv(tab, out_csv, row.names = FALSE)

cat("Wrote:", out_csv, "\n")
cat("Motifs tested per sex:", nrow(m), "male,", nrow(f), "female",
    "| identical motif sets:", identical(sort(m$Motif_ID), sort(f$Motif_ID)), "\n")
cat("Enriched (BH <= 0.05):", sum(m$Enriched == "Yes"), "male,",
    sum(f$Enriched == "Yes"), "female\n")
