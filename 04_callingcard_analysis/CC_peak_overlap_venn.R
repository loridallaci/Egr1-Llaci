## =====================================================================
## Male vs Female Calling-Card peak overlap  ->  Venn + corrected peak sets
##
## WHY THIS SCRIPT EXISTS
##   The previously committed sets (male_unique_peaks_originalPeaks.bed 7,357 /
##   female_unique_peaks_originalPeaks.bed 7,495 / shared_peaks_originalPeaks.bed
##   3,187) were produced by a notebook that is not in this repo, and they do not
##   reproduce from the raw per-sex calls. Audit (2026-07-17) found:
##     - 1,433 male + 1,090 female raw peaks silently dropped
##       (7,357+3,187 = 10,544 vs 11,977 raw male; 7,495+3,187 = 10,682 vs 11,772 raw female)
##     - 997 "male_unique" peaks in fact overlap a raw female peak (and vice versa)
##     - only ~2,122/3,187 "shared" rows exactly match a raw peak despite the
##       "originalPeaks" filename -> coordinates had been merged
##     - a single symmetric shared count (3,187) is not obtainable from a
##       peak-centric intersect, which is inherently many-to-many/asymmetric
##   This script regenerates the split reproducibly from the raw calls.
##
## bedtools NOTE
##   bedtools is not installed on this Windows box. GenomicRanges::findOverlaps
##   with default settings is semantically identical to `bedtools intersect`
##   (>= 1 bp overlap on the same chromosome, strand-agnostic). The equivalent
##   bedtools calls are given per block below so this can be checked on HTCF.
##
## TWO WAYS TO COUNT (both reported; they answer different questions)
##   (1) PEAK-CENTRIC  = `bedtools intersect -u/-v`, keeps each sex's own peaks.
##       Honest but the two "shared" counts DISAGREE (3,712 male-side vs 3,719
##       female-side) because overlap is many-to-many. Cannot draw a 2-circle Venn.
##   (2) REGION-CENTRIC = merge the union, then ask which sexes hit each region.
##       Categories are disjoint and sum exactly -> this is what the Venn uses.
##
## Inputs : raw per-sex MACCs calls (pycallingcards), unaffected by the bug above.
## Outputs: 04_callingcard_analysis/output/peak_overlap/
## =====================================================================

suppressMessages({
  library(GenomicRanges)
  library(VennDiagram)
})
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")  # silence per-plot logs

## ---- Paths ----------------------------------------------------------
cc_dir  <- "C:/Users/loril/Documents/Egr1/Egr1CC_vs_Egr1KDBulkRNA_FINAL_July2025"
out_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output/peak_overlap"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

male_file   <- file.path(cc_dir, "Male_Egr1CC_peaks_072825.txt")
female_file <- file.path(cc_dir, "Female_Egr1CC_peaks_072825.txt")
stopifnot(file.exists(male_file), file.exists(female_file))

read_cc <- function(f) {
  d <- read.delim(f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(d$Chr, IRanges(d$Start, d$End))
  mcols(gr)$gene <- d$`Gene Name1`
  mcols(gr)$dist <- d$Distance1
  sort(gr)
}
write_bed <- function(gr, f) {
  write.table(data.frame(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr)),
              file.path(out_dir, f), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

M <- read_cc(male_file)
F <- read_cc(female_file)
cat(sprintf("Raw peaks:  male %d   female %d\n\n", length(M), length(F)))

## ---- (1) PEAK-CENTRIC intersect -------------------------------------
## bedtools intersect -a male.bed -b female.bed -v   (unique)
## bedtools intersect -a male.bed -b female.bed -u   (shared, male coords)
m_hit <- overlapsAny(M, F)
f_hit <- overlapsAny(F, M)

cat("(1) PEAK-CENTRIC (each sex's own peaks, >=1bp overlap):\n")
cat(sprintf("    male   %5d = %5d unique + %5d overlapping-female\n", length(M), sum(!m_hit), sum(m_hit)))
cat(sprintf("    female %5d = %5d unique + %5d overlapping-male\n", length(F), sum(!f_hit), sum(f_hit)))
cat(sprintf("    NB shared is ASYMMETRIC (%d vs %d) - expected, overlap is many-to-many.\n\n",
            sum(m_hit), sum(f_hit)))

write_bed(M[!m_hit], "CC_male_unique_peakCentric.bed")
write_bed(F[!f_hit], "CC_female_unique_peakCentric.bed")
write_bed(M[m_hit],  "CC_male_shared_peakCentric_maleCoords.bed")
write_bed(F[f_hit],  "CC_female_shared_peakCentric_femaleCoords.bed")

## ---- (2) REGION-CENTRIC merge -> disjoint categories -----------------
## bedtools cat male.bed female.bed | bedtools sort | bedtools merge
u   <- reduce(c(granges(M), granges(F)))
inM <- overlapsAny(u, M)
inF <- overlapsAny(u, F)

male_only   <- u[ inM & !inF]
shared      <- u[ inM &  inF]
female_only <- u[!inM &  inF]

stopifnot(length(male_only) + length(shared) + length(female_only) == length(u))

cat("(2) REGION-CENTRIC (merged union; disjoint, sums exactly) <- used for the Venn:\n")
cat(sprintf("    male-only %d + shared %d + female-only %d = %d union regions\n",
            length(male_only), length(shared), length(female_only), length(u)))
cat(sprintf("    male circle  = %d   female circle = %d\n\n",
            length(male_only) + length(shared), length(female_only) + length(shared)))

write_bed(male_only,   "CC_male_only_regions.bed")
write_bed(female_only, "CC_female_only_regions.bed")
write_bed(shared,      "CC_shared_regions.bed")

## ---- Summary table ---------------------------------------------------
summ <- data.frame(
  Method   = c(rep("peak-centric", 4), rep("region-centric", 4)),
  Set      = c("male total", "male unique", "female total", "female unique",
               "male-only", "shared", "female-only", "union total"),
  N        = c(length(M), sum(!m_hit), length(F), sum(!f_hit),
               length(male_only), length(shared), length(female_only), length(u))
)
write.csv(summ, file.path(out_dir, "CC_peak_overlap_summary.csv"), row.names = FALSE)

## ---- Venn (region-centric) ------------------------------------------
## Font sizes: counts ~17pt, category labels ~18pt, title ~22pt on a 12pt device.
draw_venn <- function() {
grid.newpage()
vp <- draw.pairwise.venn(
  area1      = length(male_only) + length(shared),
  area2      = length(female_only) + length(shared),
  cross.area = length(shared),
  category   = c("Male", "Female"),
  fill       = c("#4C7FB8", "#C4622D"),
  alpha      = c(0.55, 0.55),
  lwd        = 2.5,
  col        = c("#2E4F73", "#7A3B18"),
  cex        = 1.4,     # the three counts
  cat.cex    = 1.5,     # "Male" / "Female"
  cat.col    = c("#2E4F73", "#7A3B18"),
  cat.pos    = c(-30, 30),
  cat.dist   = 0.05,
  fontfamily     = "sans",
  cat.fontfamily = "sans",
  cat.fontface   = "bold",
  margin     = 0.10,
  ind        = FALSE
)
## comma-separate the three counts (grobs carry raw labels e.g. "8204")
lab_map <- setNames(format(c(length(male_only), length(shared), length(female_only)), big.mark = ",", trim = TRUE),
                    as.character(c(length(male_only), length(shared), length(female_only))))
vp <- lapply(vp, function(g) {
  if (inherits(g, "text") && !is.null(g$label) && as.character(g$label) %in% names(lab_map))
    g$label <- lab_map[[as.character(g$label)]]
  g
})
class(vp) <- "gList"

grid.draw(vp)
grid.text("Egr1 Calling-Card peak overlap",
          y = 0.955, gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
grid.text("peaks called separately per sex",
          y = 0.905, gp = gpar(fontsize = 14, fontfamily = "sans", col = "grey30"))
grid.text(sprintf("merged union regions, ≥ 1 bp overlap  |  n = %s total", format(length(u), big.mark = ",")),
          y = 0.045, gp = gpar(fontsize = 14, fontfamily = "sans", col = "grey30"))
}

venn_pdf <- file.path(out_dir, "SupplFig4b_CC_peak_overlap_Venn.pdf")
venn_png <- file.path(out_dir, "SupplFig4b_CC_peak_overlap_Venn.png")
pdf(venn_pdf, width = 7.5, height = 7);                      draw_venn(); dev.off()
png(venn_png, width = 7.5, height = 7, units = "in", res = 200, pointsize = 12); draw_venn(); dev.off()

cat("Wrote:\n  ", venn_pdf, "\n  ", venn_png, "\n  ",
    file.path(out_dir, "CC_peak_overlap_summary.csv"), "\n  + 7 bed files in ", out_dir, "\n", sep = "")
