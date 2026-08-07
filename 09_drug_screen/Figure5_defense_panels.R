## =============================================================================
## Figure 5 / defense drug-screen panels - regenerated from the canonical data,
## emitting BOTH .pdf (vector) and .png for every panel.
##
## Source of truth:
##   confirmation_screen_sexdiff_emmeans/output/emmeans_sexdiff_results.xlsx
##     sheet "emmeans_results"  - per-compound emmeans contrasts (85 compounds)
##     sheet "PerRep_Zscores"   - per-replicate z per genotype
##
## SIGNIFICANCE CONVENTION: raw emmeans p <= 0.05 (NOT BH q).
##   This is what the defense deck and Fig5_revision panels use, and it gives
##   WT 63 / KD 53 hits and 37 Egr1-dependent. Filtering on BH q instead would
##   give 61 / 52 / 28. The two must not be mixed.
##
## Panels written to 09_drug_screen/defense_panels/:
##   panelB_scatter_MvF_z        male vs female z, WT and KD facets
##   panelB_barchart_hitcounts   sex-biased hits per genotype
##   panelC_pie_Egr1dependence   % of WT hits that are Egr1-dependent
##   panelD_pctEgr1dep_byClass   same, broken down by drug class
##
## Fonts follow the lab presentation standard (ticks >=14, axis titles >=16,
## plot titles >=18, legend >=14).
## =============================================================================
suppressMessages({ library(readxl); library(ggplot2) })

ROOT <- if (basename(getwd()) == "09_drug_screen") "." else "09_drug_screen"
SRC  <- file.path(ROOT, "confirmation_screen_sexdiff_emmeans", "output",
                  "emmeans_sexdiff_results.xlsx")
OUT  <- file.path(ROOT, "defense_panels")
dir.create(OUT, showWarnings = FALSE)

stopifnot(file.exists(SRC))
res <- as.data.frame(read_excel(SRC, sheet = "emmeans_results"))
zz  <- as.data.frame(read_excel(SRC, sheet = "PerRep_Zscores"))
cat(sprintf("loaded %d compounds\n", nrow(res)))

ALPHA <- 0.05
save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat("  wrote", name, ".pdf/.png\n")
}
base_thm <- theme_bw(base_size = 16) +
  theme(plot.title    = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 13),
        axis.title    = element_text(size = 16),
        axis.text     = element_text(size = 14, colour = "black"),
        legend.text   = element_text(size = 14),
        legend.title  = element_blank(),
        strip.text    = element_text(size = 16, face = "bold"),
        panel.grid.minor = element_blank())

## ---- hit definitions (raw p) ------------------------------------------------
wt_hit <- res$Q2WT_p <= ALPHA
kd_hit <- res$Q2KD_p <= ALPHA
dep    <- res$Q4_p  <= ALPHA          # sex difference significantly changed by KD
cat(sprintf("WT hits %d | KD hits %d | Q4 (Egr1-dependent) %d\n",
            sum(wt_hit), sum(kd_hit), sum(dep)))

## ---- Panel B: male vs female z, per genotype --------------------------------
mz <- function(pre) rowMeans(zz[, paste0(pre, "_z", 1:3)], na.rm = TRUE)
sc <- rbind(
  data.frame(geno = "WT cells", male = mz("M6_WT"), female = mz("F6_WT"),
             hit = ifelse(!wt_hit, "n.s.", ifelse(res$Q2WT_est < 0, "Male-biased hit", "Female-biased hit"))),
  data.frame(geno = "KD cells", male = mz("M6_KD"), female = mz("F6_KD"),
             hit = ifelse(!kd_hit, "n.s.", ifelse(res$Q2KD_est < 0, "Male-biased hit", "Female-biased hit"))))
sc$geno <- factor(sc$geno, c("WT cells", "KD cells"))
sc$hit  <- factor(sc$hit, c("Male-biased hit", "Female-biased hit", "n.s."))
COL <- c("Male-biased hit" = "#3B4CC0", "Female-biased hit" = "#C0392B", "n.s." = "grey65")

p <- ggplot(sc, aes(male, female, colour = hit)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.6, alpha = 0.9) +
  facet_wrap(~geno) + scale_colour_manual(values = COL) +
  labs(title = "Sex-biased drug sensitivity (male vs female z score)",
       x = "Male z score", y = "Female z score") +
  base_thm + theme(legend.position = "bottom")
save2(p, "panelB_scatter_MvF_z", 9.5, 5.8)

## ---- Panel B companion: hit counts ------------------------------------------
mk <- function(hit, est, lab) data.frame(geno = lab,
  cat = c("Male-biased hit", "Female-biased hit", "n.s."),
  n = c(sum(hit & est < 0), sum(hit & est > 0), sum(!hit)))
bc <- rbind(mk(wt_hit, res$Q2WT_est, "WT cells"), mk(kd_hit, res$Q2KD_est, "KD cells"))
bc$geno <- factor(bc$geno, c("WT cells", "KD cells"))
bc$cat  <- factor(bc$cat, c("Male-biased hit", "Female-biased hit", "n.s."))

p <- ggplot(bc, aes(geno, n, fill = cat)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.4, size = 5) +
  scale_fill_manual(values = COL) + expand_limits(y = max(bc$n) * 1.12) +
  labs(title = "Sex-biased hits per genotype", x = NULL, y = "number of drugs") +
  base_thm + theme(legend.position = "bottom")
save2(p, "panelB_barchart_hitcounts", 7.0, 5.4)

## ---- Panel C: Egr1-dependence pie (of WT hits) ------------------------------
nd <- sum(wt_hit & dep); ni <- sum(wt_hit & !dep); tot <- nd + ni
pie <- data.frame(lab = c("Egr1-dependent\n(sig. interaction)", "Egr1-independent"),
                  n = c(nd, ni))
pie$pct <- round(100 * pie$n / tot)
p <- ggplot(pie, aes("", n, fill = lab)) +
  geom_col(width = 1, colour = "white") + coord_polar("y", start = 0) +
  geom_text(aes(label = sprintf("%d\n(%d%%)", n, pct)),
            position = position_stack(vjust = 0.5), size = 6, fontface = "bold") +
  scale_fill_manual(values = c("#3B4CC0", "grey75")) +
  labs(title = sprintf("Egr1 dependence (n = %d WT sex-biased hits)", tot)) +
  theme_void(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 14), legend.title = element_blank())
save2(p, "panelC_pie_Egr1dependence", 7.2, 5.2)

## ---- Panel D: % Egr1-dependent by drug class --------------------------------
w <- res[wt_hit, ]
cl <- aggregate(cbind(dep = dep[wt_hit], tot = rep(1, sum(wt_hit))) ~ Class, data = w, FUN = sum)
cl$pct <- 100 * cl$dep / cl$tot
cl <- cl[order(cl$pct, cl$tot), ]
cl$Class <- factor(cl$Class, levels = cl$Class)

p <- ggplot(cl, aes(pct, Class)) +
  geom_col(fill = "#3B5488", width = 0.75) +
  geom_text(aes(label = sprintf("%.0f%% (%d/%d)", pct, dep, tot)),
            hjust = -0.08, size = 5, fontface = "bold") +
  scale_x_continuous(limits = c(0, 128), breaks = seq(0, 100, 25),
                     labels = paste0(seq(0, 100, 25), "%")) +
  labs(title = "Sex-biased drug vulnerability is Egr1-dependent",
       subtitle = sprintf("%d of %d WT sex-biased hits are Egr1-dependent", nd, tot),
       x = "% of WT sex-biased hits that are Egr1-dependent", y = NULL) +
  base_thm
save2(p, "panelD_pctEgr1dep_byClass", 10.5, 7.0)

cat("\nAll panels written to", normalizePath(OUT), "\n")
