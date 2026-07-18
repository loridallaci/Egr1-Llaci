## =====================================================================
## Supp Fig 4a - HOMER known-motif enrichment in the per-sex Egr1 CC peaks.
## Egr1 is highlighted among the top enriched motifs, per sex + shared.
##
## Data = existing valid HOMER runs (window-1000 per-sex calls + shared 20kb).
## The window-300 full-set (11,977/11,772) run could not be produced: HOMER's
## compiled binaries need GLIBC 2.29+ (cluster login/interactive nodes are older)
## and its genome extraction falls into a broken custom-directory mode. These
## runs give the same conclusion and are used instead. State the peak call in the
## caption.
##
## Egr1 rank: Male #4 (p=1e-161), Female #3 (p=1e-192), Shared #1 (p=1e-192).
## =====================================================================

suppressMessages({ library(ggplot2); library(dplyr); library(stringr); library(patchwork) })

eg <- "C:/Users/loril/Documents/Egr1"
files <- list(
  Male   = file.path(eg, "CallingCard_Analyses/041725/Homer/Separate_peaks/GBM_Male_Peaks_Egr1_MACC2_window1000_041725/knownResults.txt"),
  Female = file.path(eg, "CallingCard_Analyses/041725/Homer/Separate_peaks/GBM_Female_Peaks_Egr1_MACC2_window1000_041725/knownResults.txt"),
  Shared = file.path(eg, "Egr1 manuscript/Final Submission/CC/GBM_Peaks_Egr1_CalledSeparately_SharedMandF_20kb_091925/knownResults.txt")
)
npeaks <- c(Male = 2842, Female = 3803, Shared = 3187)
out_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output/figures/suppfig4a_homer"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

parse_p <- function(x) {                                   # "1e-170" -> 170 (= -log10 p)
  x <- gsub("^1e", "", x); as.numeric(x) * -1
}
top_motifs <- function(f, n = 8) {
  d <- read.delim(f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  d$neglogp <- parse_p(d$`P-value`)
  d$tf   <- sub("/.*$", "", d[[1]])                        # "Fra1(bZIP)/BT549-..." -> "Fra1(bZIP)"
  d$rank <- seq_len(nrow(d))
  d$isEgr1 <- grepl("^Egr1\\b", d[[1]])
  head(d[order(-d$neglogp), c("tf","neglogp","rank","isEgr1")], n)
}

panel <- function(setnm) {
  d <- top_motifs(files[[setnm]])
  d$tf <- factor(d$tf, levels = rev(d$tf))
  egr_rank <- d$rank[d$isEgr1][1]
  ggplot(d, aes(x = neglogp, y = tf, fill = isEgr1)) +
    geom_col() +
    scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#B8C4D0"), guide = "none") +
    labs(title = sprintf("%s CC peaks (n=%s)", setnm, format(npeaks[setnm], big.mark=",")),
         subtitle = sprintf("Egr1 motif rank #%d", egr_rank),
         x = expression(-log[10]~italic(P)), y = NULL) +
    theme_minimal(base_size = 14) +
    theme(plot.title    = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 15, colour = "#C0392B", face = "bold"),
          axis.title.x  = element_text(size = 16),
          axis.text.y   = element_text(size = 14),
          axis.text.x   = element_text(size = 14),
          panel.grid.major.y = element_blank())
}

p <- panel("Male") + panel("Female") + panel("Shared") +
  plot_annotation(title = "Egr1 motif enrichment in Calling-Card peaks (HOMER known motifs)",
                  theme = theme(plot.title = element_text(size = 18, face = "bold")))

ggsave(file.path(out_dir, "SuppFig4a_HOMER_Egr1_motif_panel.pdf"), p, width = 16, height = 6)
ggsave(file.path(out_dir, "SuppFig4a_HOMER_Egr1_motif_panel.png"), p, width = 16, height = 6, dpi = 200, bg = "white")
cat("Wrote SuppFig4a_HOMER_Egr1_motif_panel.{pdf,png} to\n  ", out_dir, "\n")
