## One PDF per regulon direction (top GO-BP-2023 terms), full wrapped names.
suppressMessages({library(ggplot2); library(dplyr); library(stringr)})
IN   <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/07_tcga_survival/output_regulon_survival/enrichr_regulon_directions"
NTOP <- 10

sets <- list(
  MaleUnique_Down   = list(title="Male regulon (Down after KD)",   col="#1E6FBF"),
  MaleUnique_Up     = list(title="Male regulon (Up after KD)",     col="#8FC1EC"),
  FemaleUnique_Down = list(title="Female regulon (Down after KD)", col="#C2185B"),
  FemaleUnique_Up   = list(title="Female regulon (Up after KD)",   col="#F48FB1"))

for (s in names(sets)) {
  d <- read.csv(file.path(IN, paste0("enrichr_GO_BP2023_", s, ".csv")))
  d <- d[order(d$Adjusted.P.value), ][seq_len(min(NTOP, nrow(d))), ]
  d$Term <- sub("\\s*\\(GO:\\d+\\)$", "", d$Term)           # drop GO id
  d$Term <- str_wrap(d$Term, width = 40)                    # wrap long names
  d$neglog10 <- -log10(d$Adjusted.P.value)
  d$row <- seq_len(nrow(d))
  d$uid <- factor(sprintf("%02d", d$row), levels = sprintf("%02d", rev(d$row)))
  labmap <- setNames(d$Term, as.character(d$uid))

  p <- ggplot(d, aes(neglog10, uid)) +
    geom_col(fill = sets[[s]]$col, width = 0.72) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.6, color = "grey40") +
    scale_y_discrete(labels = labmap) +
    labs(x = expression(-log[10]~"(adjusted "*italic(P)*"-value)"), y = NULL,
         title = sets[[s]]$title, subtitle = "GO Biological Process 2023") +
    theme_bw(base_size = 15) +
    theme(axis.text.y  = element_text(size = 14, color = "black"),
          axis.text.x  = element_text(size = 14, color = "black"),
          axis.title.x = element_text(size = 16),
          plot.title   = element_text(size = 18, face = "bold"),
          plot.subtitle= element_text(size = 14, color = "grey35"),
          panel.grid.minor = element_blank())

  outf <- file.path(IN, paste0("enrichr_GO_BP2023_", s, ".pdf"))
  ggsave(outf, p, width = 10, height = 6, device = "pdf")
  cat("wrote:", basename(outf), "\n")
}
