## =====================================================================
## Supp Fig 2 (extra panel): EGR1 chromVAR motif activity by cell type x sex,
## pooled across the 5 developmental stages (paper's cells + author labels).
## Hypothesis-generating: shows the compartment where EGR1 activity is
## sex-biased differs by lineage (myeloid vs glial-precursor).
##
## Input : output/paperCells_panels/Egr1motif_activity_<stage>.csv
##         (per-cell: barcode, sex, celltype, atac coords, egr1_z)
## Fonts follow the lab presentation minimums.
## =====================================================================
suppressMessages({ library(ggplot2) })

d_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/06_cortex_development/output/paperCells_panels"
files <- list.files(d_dir, "Egr1motif_activity_.*csv", full.names = TRUE)
rd <- function(f) { x <- read.csv(f); x$stage <- gsub("^.*activity_|[.]csv$", "", basename(f)); x }
all <- do.call(rbind, lapply(files, rd))

cts <- unique(all$celltype)
summ <- do.call(rbind, lapply(cts, function(ct) {
  s <- all[all$celltype == ct, ]
  m <- s$egr1_z[s$sex == "male"]; f <- s$egr1_z[s$sex == "female"]
  if (length(m) < 3 || length(f) < 3) return(NULL)
  data.frame(celltype = ct, n_M = length(m), n_F = length(f),
             mean_M = mean(m), mean_F = mean(f),
             sem_M = sd(m)/sqrt(length(m)), sem_F = sd(f)/sqrt(length(f)),
             wilcox_p = suppressWarnings(wilcox.test(m, f)$p.value))
}))
summ$p_adj  <- p.adjust(summ$wilcox_p, "BH")
summ$star   <- cut(summ$p_adj, c(-Inf,1e-3,1e-2,5e-2,Inf), labels = c("***","**","*","ns"))
summ$overall <- (summ$mean_M * summ$n_M + summ$mean_F * summ$n_F) / (summ$n_M + summ$n_F)
summ <- summ[order(summ$overall), ]
summ$celltype <- factor(summ$celltype, levels = summ$celltype)
write.csv(summ[order(-summ$overall), ], file.path(d_dir, "EGR1_activity_byCelltypeSex_summary.csv"), row.names = FALSE)

long <- rbind(
  data.frame(celltype = summ$celltype, sex = "male",   z = summ$mean_M, sem = summ$sem_M),
  data.frame(celltype = summ$celltype, sex = "female", z = summ$mean_F, sem = summ$sem_F))
long$sex <- factor(long$sex, levels = c("male","female"))
sex_cols <- c(male = "#4A6FE3", female = "#F39AC9")
star_x <- pmax(summ$mean_M, summ$mean_F) + 0.06

p <- ggplot(long, aes(z, celltype, colour = sex)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_segment(data = summ, inherit.aes = FALSE,
               aes(x = mean_M, xend = mean_F, y = celltype, yend = celltype),
               colour = "grey70", linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(data = data.frame(celltype = summ$celltype, star = summ$star, x = star_x),
            inherit.aes = FALSE, aes(x = x, y = celltype, label = star),
            hjust = 0, size = 5, colour = "grey25") +
  scale_colour_manual(values = sex_cols, name = "Sex",
                      guide = guide_legend(override.aes = list(size = 4))) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  labs(x = "EGR1 motif activity (mean chromVAR z)", y = NULL,
       title = "EGR1 motif activity by cell type and sex",
       subtitle = "Developing human cortex (5 stages pooled); dashed line = mean accessibility") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 13),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        legend.position = "top")

ggsave(file.path(d_dir, "SuppFig2_EGR1_activity_byCelltypeSex.pdf"), p, width = 8, height = 7)
ggsave(file.path(d_dir, "SuppFig2_EGR1_activity_byCelltypeSex.png"), p, width = 8, height = 7,
       dpi = 300, bg = "white")
cat("Wrote SuppFig2_EGR1_activity_byCelltypeSex.{pdf,png} + summary CSV\n")
