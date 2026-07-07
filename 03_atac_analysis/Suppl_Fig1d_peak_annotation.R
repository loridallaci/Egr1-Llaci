## =====================================================================
## Supplementary Fig 1d : genomic-feature distribution of ATAC peaks
##   "Promoter, gene body etc distribution of data"
##
## ChIPseeker annotation of the lot6 multiome peaks:
##   plotAnnoBar  -> % promoter / 5'UTR / 3'UTR / exon / intron /
##                   downstream / distal-intergenic
## (bar plot only; the pie and distance-to-TSS panels are not used.)
##
## Done for: all peaks (per sex) AND the sex-biased DAR peaks
##   (Male-up / Female-up, from DARs.R).
##
## Extracted/cleaned from Multiome_update_01292025.Rmd.
## NOTE on peak set: the original Rmd annotated the "ATAC" (Cell Ranger ARC)
##   assay for the overall distribution. Here we default to the MACS2 "peaks"
##   assay so panel d uses the SAME peak set as the DAR volcano (DARs.R) -
##   set peak_assay <- "ATAC" to reproduce the original exactly.
## =====================================================================

suppressMessages({
  library(Seurat); library(Signac)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)   # mm10
  library(org.Mm.eg.db)
  library(GenomicRanges); library(ggplot2)
})

txdb       <- TxDb.Mmusculus.UCSC.mm10.knownGene
tss_region <- c(-2000, 500)                      # same window as the original Rmd
peak_assay <- "peaks"                            # MACS2 peaks (set "ATAC" for Cell Ranger peaks)

## ---- Paths (mirror the other 03_atac_analysis scripts) --------------
## original (author's machine): "/home/lllaci/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds"
rds_path  <- "C:/Users/loril/Documents/data/multiome_081721_aggregated_data_analysis/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks.rds"
dar_csv   <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/data output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv"  # from DARs.R
fig_dir   <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/03_atac_analysis/output/figures/suppl_fig1"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

lot6 <- readRDS(rds_path)
lot6$sex <- ifelse(grepl("-1$", colnames(lot6)), "female",
             ifelse(grepl("-2$", colnames(lot6)), "male", NA))

## ---- helper: annotate a GRanges + write the bar plot ----------------
annotate_and_plot <- function(gr, tag, title = tag) {
  pa <- annotatePeak(gr, TxDb = txdb, tssRegion = tss_region,
                     annoDb = "org.Mm.eg.db", verbose = FALSE)
  pdf(file.path(fig_dir, paste0("ATAC_peaks_annoBar_",   tag, ".pdf")), width = 10, height = 6)
  print(plotAnnoBar(pa, title = title));      dev.off()
  invisible(pa)
}

## ---- 1. all peaks, per sex -----------------------------------------
for (sx in c("male", "female")) {
  sub <- subset(lot6, subset = sex == sx)
  gr  <- granges(sub[[peak_assay]])
  annotate_and_plot(gr, paste0(peak_assay, "_all_", sx, "_lot6"),
                    title = paste0(sx, " - all ", peak_assay, " peaks"))
  cat("annotated", length(gr), peak_assay, "peaks for", sx, "\n")
}

## ---- 2. sex-biased DAR peaks (Male-up / Female-up) -----------------
if (file.exists(dar_csv)) {
  da <- read.csv(dar_csv)
  if (!"significance" %in% colnames(da)) {
    da$significance <- "Not Sig"
    da$significance[da$avg_log2FC >=  0.5 & da$p_val_adj <= 0.05] <- "Male-up"
    da$significance[da$avg_log2FC <= -0.5 & da$p_val_adj <= 0.05] <- "Female-up"
  }
  peak_id <- if ("X" %in% colnames(da)) da$X else da[[1]]
  male_gr   <- StringToGRanges(peak_id[da$significance == "Male-up"])
  female_gr <- StringToGRanges(peak_id[da$significance == "Female-up"])
  if (length(male_gr))   annotate_and_plot(male_gr,   "DAR_Male_up_lot6",   "Male-biased DAR peaks")
  if (length(female_gr)) annotate_and_plot(female_gr, "DAR_Female_up_lot6", "Female-biased DAR peaks")
  cat("annotated DAR peaks:", length(male_gr), "Male-up |", length(female_gr), "Female-up\n")
} else {
  cat("NOTE: DAR csv not found (", dar_csv, ") - run DARs.R first for the sex-biased panels.\n")
}

cat("\nDone. Panel-1d figures in:\n  ", fig_dir, "\n")
