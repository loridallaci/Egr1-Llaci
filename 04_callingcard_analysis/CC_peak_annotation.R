## =====================================================================
## Genomic-feature distribution (promoter / UTR / exon / intron /
##   downstream / distal-intergenic) of the Egr1 Calling-Card peaks.
##
## Same recipe as the multiome DAR annotation (Suppl_Fig1d_peak_annotation.R):
##   ChIPseeker::annotatePeak -> plotAnnoPie, mm10 TxDb, tssRegion (-2000,+500).
##
## UNLIKE the ATAC multiome (one shared MACS2 peak set for all cells), the
## Calling-Card peaks ARE called separately per sex, so the full male vs full
## female peak distributions are genuinely comparable, not identical by
## construction. We annotate:
##   - Male  (all)  = male-unique + shared
##   - Female(all)  = female-unique + shared
##   - Male-unique  (sex-biased)
##   - Female-unique(sex-biased)
##   - Shared
##
## Inputs : the CORRECTED peak-centric beds in output/peak_overlap/, produced by
##   CC_peak_overlap_venn.py|R from the raw MACCs calls (male 11,977 / female 11,772).
##
## 2026-07-17 -- REPOINTED. This script previously read
##   {male,female}_unique_peaks_originalPeaks.bed / shared_peaks_originalPeaks.bed
##   (7,357 / 7,495 / 3,187). Those were produced by a notebook using pyranges
##   .subtract()/.intersect(), which do INTERVAL ARITHMETIC: they trim and split
##   peaks into fragments rather than classifying whole peaks. They also ran on the
##   20 kb-gene-filtered subset, not the raw calls. Net effect: ~2,500 peaks missing,
##   997 "unique" peaks that actually overlap the other sex, and clipped coordinates.
##   Since a feature pie measures where a peak sits relative to a TSS, trimmed
##   coordinates directly bias it -- hence this re-run.
##
## Peak-centric sets are used (whole peaks, ORIGINAL coordinates) because the unit
## being annotated should be a real peak, not a merged construct. They sum exactly:
##   male_unique 8,265 + male_shared 3,712 = 11,977 raw male
##   female_unique 8,053 + female_shared 3,719 = 11,772 raw female
## "Shared" is annotated from BOTH sides: overlap is many-to-many, so the male-side
## and female-side shared sets are different peaks (3,712 vs 3,719), not one set.
## =====================================================================

suppressMessages({
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)   # mm10
  library(org.Mm.eg.db)
  library(GenomicRanges)
})

txdb       <- TxDb.Mmusculus.UCSC.mm10.knownGene
tss_region <- c(-2000, 500)                     # same window as the DAR / multiome annotation

## ---- Paths ----------------------------------------------------------
cc_dir  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output/peak_overlap"
fig_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output/figures/peak_annotation"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

read_bed <- function(f) {
  b <- read.table(file.path(cc_dir, f), sep = "\t", header = FALSE,
                  stringsAsFactors = FALSE, colClasses = c("character","integer","integer"))[, 1:3]
  GRanges(b[[1]], IRanges(b[[2]] + 1L, b[[3]]))  # BED is 0-based half-open -> 1-based
}

male_uni    <- read_bed("CC_male_unique_peakCentric.bed")                    # 8,265
female_uni  <- read_bed("CC_female_unique_peakCentric.bed")                  # 8,053
male_shr    <- read_bed("CC_male_shared_peakCentric_maleCoords.bed")         # 3,712
female_shr  <- read_bed("CC_female_shared_peakCentric_femaleCoords.bed")     # 3,719

male_all   <- c(male_uni,   male_shr)      # == 11,977 raw male peaks
female_all <- c(female_uni, female_shr)    # == 11,772 raw female peaks
stopifnot(length(male_all) == 11977, length(female_all) == 11772)

sets <- list(
  Male_all      = list(gr = male_all,   title = "Male Egr1 CC peaks (all)"),
  Female_all    = list(gr = female_all, title = "Female Egr1 CC peaks (all)"),
  Male_unique   = list(gr = male_uni,   title = "Male-unique Egr1 CC peaks"),
  Female_unique = list(gr = female_uni, title = "Female-unique Egr1 CC peaks"),
  Shared_maleCoords   = list(gr = male_shr,   title = "Shared Egr1 CC peaks (male coords)"),
  Shared_femaleCoords = list(gr = female_shr, title = "Shared Egr1 CC peaks (female coords)")
)

## ---- annotate + pie, and collect feature % for a summary table ------
feat_rows <- list()
for (tag in names(sets)) {
  gr <- sets[[tag]]$gr
  pa <- annotatePeak(gr, TxDb = txdb, tssRegion = tss_region,
                     annoDb = "org.Mm.eg.db", verbose = FALSE)

  pdf(file.path(fig_dir, paste0("CC_peaks_annoPie_", tag, ".pdf")), width = 8, height = 8)
  plotAnnoPie(pa, main = sets[[tag]]$title, cex = 1.4)   # legible from a distance
  dev.off()

  st <- pa@annoStat                       # Feature / Frequency (%)
  st$Set <- tag
  feat_rows[[tag]] <- st
  cat(sprintf("%-14s %6d peaks -> CC_peaks_annoPie_%s.pdf\n", tag, length(gr), tag))
}

## ---- tidy summary table + grouped barplot (easier to compare) -------
feat <- do.call(rbind, feat_rows)
wide <- reshape(feat, idvar = "Feature", timevar = "Set", direction = "wide")
names(wide) <- sub("^Frequency\\.", "", names(wide))
write.csv(wide, file.path(fig_dir, "CC_peak_feature_distribution.csv"), row.names = FALSE)

## grouped bar: feature composition side-by-side across the five sets
mat <- tapply(feat$Frequency, list(feat$Feature, factor(feat$Set, levels = names(sets))),
              sum)
mat[is.na(mat)] <- 0
pdf(file.path(fig_dir, "CC_peaks_feature_barplot.pdf"), width = 11, height = 7)
op <- par(mar = c(6, 5, 4, 12), xpd = TRUE)
cols <- grDevices::hcl.colors(nrow(mat), "Spectral")
bp <- barplot(mat, beside = FALSE, col = cols, las = 2,
              ylab = "% of peaks", cex.lab = 1.4, cex.axis = 1.2, cex.names = 1.2,
              main = "Egr1 CC peak genomic-feature distribution", cex.main = 1.5)
legend("topright", inset = c(-0.28, 0), legend = rownames(mat), fill = cols,
       bty = "n", cex = 1.1, title = "Feature")
par(op)
dev.off()

cat("\nDone. Pies + barplot + CSV in:\n  ", fig_dir, "\n\n")
print(wide, row.names = FALSE)
