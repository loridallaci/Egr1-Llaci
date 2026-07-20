## Figure 5 panel b — sex-biased drug sensitivity (Male vs Female combined Z),
## WT | KD facets, hits colored by direction, with the 3 KD-only drugs labeled.
## Reads the confirmation-screen combined-Z tab directly (thioguanine kept as 2 runs).
suppressMessages({library(readxl); library(ggplot2); library(ggrepel)})

wb  <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/Cell Titer Glo/Maxene Drug Screen/85compoundsCheck_analyses_July2026.xlsx"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/09_drug_screen"

d <- suppressMessages(read_excel(wb, sheet = "Zc combined", col_names = FALSE))[-1, ]
df <- data.frame(name = as.character(d[[9]]),
                 M_WT = as.numeric(d[[3]]), M_KD = as.numeric(d[[4]]),
                 F_WT = as.numeric(d[[5]]), F_KD = as.numeric(d[[6]]))

q2  <- function(a, b) (a - b) / sqrt(2)          # <0 => male more sensitive
q2WT <- q2(df$M_WT, df$F_WT); qWT <- p.adjust(2 * pnorm(-abs(q2WT)), "BH")
q2KD <- q2(df$M_KD, df$F_KD); qKD <- p.adjust(2 * pnorm(-abs(q2KD)), "BH")
sigWT <- qWT <= 0.05; sigKD <- qKD <= 0.05

wt <- data.frame(xM = df$M_WT, yF = df$F_WT, q2 = q2WT, q = qWT, name = df$name,
                 geno = "WT cells", kdonly = FALSE)
kd <- data.frame(xM = df$M_KD, yF = df$F_KD, q2 = q2KD, q = qKD, name = df$name,
                 geno = "KD cells", kdonly = sigKD & !sigWT)
all <- rbind(wt, kd)
all$geno <- factor(all$geno, levels = c("WT cells", "KD cells"))

all$cat <- "n.s. (q > 0.05)"
all$cat[all$q <= 0.05 & all$q2 < 0] <- "Male-biased hit (q ≤ 0.05)"
all$cat[all$q <= 0.05 & all$q2 > 0] <- "Female-biased hit (q ≤ 0.05)"
all$cat <- factor(all$cat, levels = c("Male-biased hit (q ≤ 0.05)",
                                      "Female-biased hit (q ≤ 0.05)", "n.s. (q > 0.05)"))

## labels: the 3 KD-only drugs, on the KD facet only
clean <- function(n) sub(" (Citrate|Acetate|Nitrate|Hydrochloride|Sulfate|Sodium|Mesylate|Pamoate|Pivalate|Calcium|Hemisulfate)$",
                         "", tools::toTitleCase(tolower(n)))
all$label <- NA_character_
sel <- all$geno == "KD cells" & all$kdonly
all$label[sel] <- clean(all$name[sel])
all$label[!is.na(all$label) & grepl("Prasterone", all$label)] <- "Prasterone (DHEA)"
cat("KD-only labeled:", paste(na.omit(all$label), collapse = ", "), "\n")

cols <- c("Male-biased hit (q ≤ 0.05)" = "#1874CD",
          "Female-biased hit (q ≤ 0.05)" = "#C71585",
          "n.s. (q > 0.05)" = "grey72")
lo <- min(c(all$xM, all$yF), na.rm = TRUE) * 1.03; lim <- 1

p <- ggplot(all, aes(xM, yF, fill = cat)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(aes(size = cat), shape = 21, color = "grey25", stroke = 0.3, alpha = 0.9) +
  geom_text_repel(aes(label = label), data = subset(all, !is.na(label)),
                  size = 4.4, fontface = "italic", color = "black", na.rm = TRUE,
                  box.padding = 0.7, point.padding = 0.4, min.segment.length = 0,
                  segment.color = "grey40", segment.size = 0.3,
                  nudge_x = 2.5, nudge_y = -3, seed = 1) +
  facet_wrap(~geno) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_size_manual(values = c(3.2, 3.2, 2.2), guide = "none") +
  coord_equal(xlim = c(lo, lim), ylim = c(lo, lim)) +
  labs(x = "Male Egr1 combined Z score", y = "Female Egr1 combined Z score",
       title = "Sex-biased drug sensitivity (Male vs Female Z score)") +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        panel.grid.minor = element_blank(), aspect.ratio = 1)

ggsave(file.path(out, "Figure5_panelB_sexdiff_scatter.png"), p, width = 11, height = 6.4, dpi = 200, bg = "white")
ggsave(file.path(out, "Figure5_panelB_sexdiff_scatter.pdf"), p, width = 11, height = 6.4)
cat(sprintf("WT hits: male %d female %d | KD hits: male %d female %d\nsaved to %s\n",
            sum(sigWT & q2WT < 0), sum(sigWT & q2WT > 0),
            sum(sigKD & q2KD < 0), sum(sigKD & q2KD > 0), out))
