## =====================================================================
## Supplementary Fig 1d : genomic-feature distribution (pie) of the
##   sex-biased DAR peaks - Male-biased vs Female-biased.
##
## Why DAR peaks and not "all peaks per sex": there is ONE MACS2 `peaks`
## set for the whole multiome object, shared across all cells. Subsetting
## the object to male vs female cells does NOT change the peak coordinates,
## so an "all peaks" pie is identical for M and F by construction. The
## meaningful per-sex comparison is the sex-biased DAR peaks (Male-up /
## Female-up, from DARs.R), whose genomic-feature distributions can differ.
##
## ChIPseeker plotAnnoPie -> % promoter / 5'UTR / 3'UTR / exon / intron /
##   downstream / distal-intergenic, for the Male-up and Female-up peak sets.
##
## Inputs : the DAR result csv (from DARs.R) only - no Seurat object needed,
##   so this runs anywhere (no Signac/Matrix stack required).
##   ident.1 = male, ident.2 = female  =>  avg_log2FC > 0 is Male-biased.
##   Cutoffs match the Fig 1C DAR volcano: |avg_log2FC| >= 0.5 & p_adj <= 0.05.
## =====================================================================

suppressMessages({
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)   # mm10
  library(org.Mm.eg.db)
  library(GenomicRanges)
})

txdb       <- TxDb.Mmusculus.UCSC.mm10.knownGene
tss_region <- c(-2000, 500)                      # same window as the original Rmd

## ---- Paths ----------------------------------------------------------
dar_csv <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/data output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv"  # from DARs.R
fig_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/03_atac_analysis/output/figures/suppl_fig1"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## peak id "chr-start-end" -> GRanges (regex is robust to chr names; no Signac)
ids_to_granges <- function(ids) {
  m   <- regmatches(ids, regexec("^(.*)-([0-9]+)-([0-9]+)$", ids))
  chr <- vapply(m, `[`, character(1), 2)
  s   <- as.integer(vapply(m, `[`, character(1), 3))
  e   <- as.integer(vapply(m, `[`, character(1), 4))
  GRanges(chr, IRanges(s, e))
}

annotate_and_pie <- function(ids, tag, title) {
  gr <- ids_to_granges(ids)
  pa <- annotatePeak(gr, TxDb = txdb, tssRegion = tss_region,
                     annoDb = "org.Mm.eg.db", verbose = FALSE)
  pdf(file.path(fig_dir, paste0("ATAC_peaks_annoPie_", tag, ".pdf")), width = 8, height = 8)
  plotAnnoPie(pa, main = title, cex = 1.4)   # cex bump: legible from a distance
  dev.off()
  cat(sprintf("%-20s %5d peaks -> %s\n", tag, length(gr),
              paste0("ATAC_peaks_annoPie_", tag, ".pdf")))
  invisible(pa)
}

## ---- sex-biased DAR peaks ------------------------------------------
stopifnot(file.exists(dar_csv))
da      <- read.csv(dar_csv, stringsAsFactors = FALSE)
peak_id <- da[[1]]                                  # first (unnamed) col = chr-start-end
male_up   <- peak_id[da$avg_log2FC >=  0.5 & da$p_val_adj <= 0.05]
female_up <- peak_id[da$avg_log2FC <= -0.5 & da$p_val_adj <= 0.05]

annotate_and_pie(male_up,   "DAR_Male_up_lot6",   "Male-biased DAR peaks")
annotate_and_pie(female_up, "DAR_Female_up_lot6", "Female-biased DAR peaks")

cat("\nDone. Panel-1d sex-biased DAR pies in:\n  ", fig_dir, "\n")
