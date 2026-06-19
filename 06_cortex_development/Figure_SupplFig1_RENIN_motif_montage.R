# =============================================================================
# Cortex-development RENIN motif enrichment — Males x Females across 5 stages
# Builds the 2 x 5 volcano montage (Fold enrichment vs -log10 p-value).
#
# This is the PLOT-ONLY step: it reads the per-stage / per-sex motif-enrichment
# tables (*_AllCells_all_motifs.csv, produced by cortex_dev_RENIN_pipeline.R on
# the HTCF cluster) and regenerates every panel + the assembled figure. It runs
# locally in seconds — no Seurat/Signac/HTCF data required.
#
# Inputs : data_motifs/<Stage>_<Sex>_AllCells_all_motifs.csv  (10 files)
# Outputs: output/SupplFig1_RENIN_motif_montage.{png,pdf}
#          output/panels/RENIN_<Stage>_<Sex>.pdf              (per-panel)
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# Resolve paths relative to this script so it runs from anywhere
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()

in_dir  <- file.path(script_dir, "data_motifs")
out_dir <- file.path(script_dir, "output")
panel_dir <- file.path(out_dir, "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

# Figure layout: rows = sex, columns = developmental stage (left -> right)
stages      <- c("LateFetal", "Infant", "Child", "Adol", "Adult")
stage_label <- c(LateFetal = "Late Fetal", Infant = "Infant", Child = "Child",
                 Adol = "Adolescence", Adult = "Adult")
sexes       <- c("Male", "Female")
num_top_motifs <- 20

# ---- single-panel volcano (matches analyze_and_plot_motifs in the pipeline) ----
make_panel <- function(stage, sex) {
  csv <- file.path(in_dir, paste0(stage, "_", sex, "_AllCells_all_motifs.csv"))
  if (!file.exists(csv)) {
    warning("Missing: ", csv); return(NULL)
  }
  m <- read.csv(csv, row.names = 1, check.names = FALSE)
  if (is.null(m$logp)) m$logp <- -log10(m$pvalue)
  m <- m[order(m$logp, decreasing = TRUE), ]

  m$significant <- m$p.adjust <= 0.05
  m$color  <- ifelse(m$significant, "green3", "grey80")
  m$border <- ifelse(m$significant, 0.4, 0.1)
  m$label  <- ifelse(rank(-m$logp) <= num_top_motifs, m$motif.name, "")

  max_y <- max(m$logp, na.rm = TRUE) * 1.05
  max_x <- max(m$fold.enrichment, na.rm = TRUE) * 1.05

  # Top row carries the stage title; first column carries the sex row-label
  is_top  <- sex == sexes[1]
  is_left <- stage == stages[1]
  ylab <- if (is_left) paste0(ifelse(sex == "Male", "Males", "Females"),
                              "\n-log10(p-value)") else "-log10(p-value)"

  ggplot(m, aes(x = fold.enrichment, y = logp, label = label)) +
    geom_point(aes(fill = color), color = "black",
               stroke = m$border, pch = 21, size = 1.6) +
    scale_fill_identity() +
    geom_text_repel(max.overlaps = 500, size = 2.6, point.padding = 0.5,
                    force = 2, box.padding = 0.4, min.segment.length = 0,
                    segment.color = "grey40", segment.size = 0.3,
                    segment.alpha = 0.8) +
    theme_classic() +
    labs(x = "Fold enrichment", y = ylab,
         title = if (is_top) stage_label[[stage]] else NULL) +
    xlim(c(0, max_x)) + ylim(c(0, max_y)) +
    # Presentation-legible fonts (>=14 ticks, >=16 axis titles, >=18 title)
    theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 16),
          axis.text  = element_text(size = 14))
}

# ---- build all 10 panels, save each, then assemble ----
panels <- list()
for (sex in sexes) for (stage in stages) {
  p <- make_panel(stage, sex)
  panels[[paste(sex, stage)]] <- p
  if (!is.null(p)) {
    pdf(file.path(panel_dir, paste0("RENIN_", stage, "_", sex, ".pdf")),
        width = 5, height = 5)
    print(p); dev.off()
  }
}

# Row 1 = Males, Row 2 = Females; columns follow `stages`
ordered <- c(lapply(stages, function(s) panels[[paste("Male", s)]]),
             lapply(stages, function(s) panels[[paste("Female", s)]]))
montage <- wrap_plots(ordered, nrow = 2, ncol = 2 * 0 + 5) +
  plot_annotation(
    title = "RENIN motif enrichment across cortex development (Males vs Females)",
    theme = theme(plot.title = element_text(size = 18, face = "bold")))

ggsave(file.path(out_dir, "SupplFig1_RENIN_motif_montage.png"),
       montage, width = 25, height = 10, dpi = 300, limitsize = FALSE)
ggsave(file.path(out_dir, "SupplFig1_RENIN_motif_montage.pdf"),
       montage, width = 25, height = 10, limitsize = FALSE)

cat("Done. Montage + per-panel PDFs written to:\n  ", out_dir, "\n")
