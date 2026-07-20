## =====================================================================
## Supp Fig 4a - Egr1 motif enrichment in the per-sex Egr1 CC peaks.
##
## Egr1-ONLY panel (not the full ranked motif list): the point is that the
## Egr1 motif is strongly enriched in each sex's Calling Card peaks over a
## matched genomic background, confirming insertions are Egr1-directed.
##
## Data = final v5 HOMER run on the FULL per-sex peaks (Egr1 vs matched WT):
##   suppfig4a_homer_perSex.sh -> output/suppfig4a_homer_v5/*_knownResults.txt
##   Male   n = 11,978 peaks ; Egr1(Zf) P = 1e-291 ; 41.64% peaks vs 26.22% bg
##   Female n = 11,772 peaks ; Egr1(Zf) P = 1e-300 ; 43.78% peaks vs 27.78% bg
## The Egr1 consensus (TGCGTGGGYG) is identical in both sexes.
## PWM below is HOMER's known Egr1 motif (known11.motif from the v5 run).
##
## NOTE: the sex-vs-sex / depth-matched (downsampling) differential analysis
## is a SEPARATE panel/figure - do not merge it into 4a.
## Figure fonts follow the lab presentation minimums (ticks>=14, titles>=16/18).
## =====================================================================

suppressMessages({ library(ggplot2) })
have_logo <- requireNamespace("ggseqlogo", quietly = TRUE) &&
             requireNamespace("patchwork", quietly = TRUE)

out_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output/figures/suppfig4a_homer"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- enrichment numbers (from output/suppfig4a_homer_v5/*_knownResults.txt) ----
enr <- data.frame(
  Sex   = factor(rep(c("Male", "Female"), each = 2), levels = c("Male", "Female")),
  Set   = factor(rep(c("Egr1 CC peaks", "Genomic background"), 2),
                 levels = c("Egr1 CC peaks", "Genomic background")),
  Pct   = c(41.64, 26.22, 43.78, 27.78)
)
pann <- data.frame(                         # p-value labels, one per sex
  Sex = factor(c("Male", "Female"), levels = c("Male", "Female")),
  lab = c("italic(P)==1%*%10^-291", "italic(P)==1%*%10^-300"),
  y   = c(52, 52)
)

## ---- HOMER known Egr1 PWM (rows = position, cols = A C G T) ----
homer <- matrix(c(
  0.128, 0.072, 0.142, 0.658,
  0.078, 0.036, 0.882, 0.004,
  0.154, 0.523, 0.023, 0.300,
  0.001, 0.001, 0.997, 0.001,
  0.001, 0.001, 0.282, 0.716,
  0.027, 0.001, 0.971, 0.001,
  0.001, 0.002, 0.973, 0.024,
  0.001, 0.001, 0.997, 0.001,
  0.153, 0.415, 0.010, 0.422,
  0.034, 0.002, 0.940, 0.024), ncol = 4, byrow = TRUE)
colnames(homer) <- c("A", "C", "G", "T")

## ---- bar: % of regions with the Egr1 motif, CC peaks vs background, per sex ----
bar <- ggplot(enr, aes(Sex, Pct, fill = Set)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = sprintf("%.1f%%", Pct)),
            position = position_dodge(width = 0.72), vjust = -0.35, size = 5) +
  geom_text(data = pann, aes(Sex, y, label = lab), parse = TRUE,
            inherit.aes = FALSE, size = 5.2) +
  scale_fill_manual(values = c("Egr1 CC peaks" = "#C0392B",
                               "Genomic background" = "#B8C4D0")) +
  scale_y_continuous(limits = c(0, 56), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "% of regions with Egr1 motif", fill = NULL) +
  theme_classic(base_size = 16) +
  theme(plot.title   = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 16),
        axis.text    = element_text(size = 14, colour = "black"),
        legend.text  = element_text(size = 14),
        legend.position = "top")

if (have_logo) {
  pwm <- t(homer); rownames(pwm) <- c("A", "C", "G", "T")
  logo <- ggseqlogo::ggseqlogo(pwm, method = "bits") +
    ggtitle("Egr1 motif (TGCGTGGGYG)") +
    theme(plot.title = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 16),
          axis.text  = element_text(size = 14, colour = "black"))
  p <- patchwork::wrap_plots(logo, bar, ncol = 1, heights = c(1, 1.6))
} else {
  message("ggseqlogo/patchwork not installed - writing the bar panel only. ",
          "install.packages(c('ggseqlogo','patchwork')) to add the motif logo.")
  p <- bar + ggtitle("Egr1 motif enrichment in Calling Card peaks")
}

ggsave(file.path(out_dir, "SuppFig4a_Egr1_motif.pdf"), p, width = 7, height = 7)
ggsave(file.path(out_dir, "SuppFig4a_Egr1_motif.png"), p, width = 7, height = 7,
       dpi = 300, bg = "white")
cat("Wrote SuppFig4a_Egr1_motif.{pdf,png} to\n  ", out_dir, "\n")
