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
##   panelB_scatter_MvF_z         male vs female z, WT and KD facets
##   panelB_barchart_hitcounts    sex-biased hits per genotype (WT 61+2, KD 51+2)
##   panelC_pie_Egr1dependence    Egr1-dependence of the WT-or-KD UNION (n=64: 37/27)
##   panelD_pctEgr1dep_byClass    Egr1-dependence of WT hits ONLY (n=63: 37/26), by class
##   panelE_dumbbell_WTvsKD_n64   per-compound WT->KD shift for the union (n=64)
##
## NOTE the two denominators are deliberately different and must not be swapped:
##   panel C = WT-or-KD union      = 64 rows -> 37 dependent / 27 independent
##   panel D = WT sex-biased hits  = 63 rows -> 37 dependent / 26 independent
## They differ because 1 compound (Vanillin) is a KD-only hit, and THIOGUANINE
## appears twice in the 85-compound plate, so union-by-row (64) exceeds
## union-by-name (63).
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
  ok <- c(pdf = FALSE, png = FALSE)
  ok["pdf"] <- tryCatch({ ggsave(file.path(OUT, paste0(name, ".pdf")), p,
                                 width = w, height = h, device = cairo_pdf); TRUE },
                        error = function(e) FALSE)
  ok["png"] <- tryCatch({ ggsave(file.path(OUT, paste0(name, ".png")), p,
                                 width = w, height = h, dpi = 300); TRUE },
                        error = function(e) FALSE)
  cat(sprintf("  %-32s pdf:%s png:%s%s\n", name,
              ifelse(ok["pdf"], "ok", "LOCKED"), ifelse(ok["png"], "ok", "LOCKED"),
              if (all(ok)) "" else "   <-- close the open file and re-run"))
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
wt_hit  <- res$Q2WT_p <= ALPHA
kd_hit  <- res$Q2KD_p <= ALPHA
dep     <- res$Q4_p  <= ALPHA         # sex difference significantly changed by KD
uni_hit <- wt_hit | kd_hit            # sex-biased in EITHER genotype

## Direction across the union: a compound's sex is taken from whichever
## genotype it is a hit in. No compound flips direction between genotypes.
uni_dir <- ifelse(wt_hit, ifelse(res$Q2WT_est < 0, "M", "F"),
                          ifelse(res$Q2KD_est < 0, "M", "F"))
uni_M <- sum(uni_hit & uni_dir == "M")
uni_F <- sum(uni_hit & uni_dir == "F")

cat(sprintf("WT hits %d (%dM/%dF) | KD hits %d (%dM/%dF) | union %d (%dM/%dF) | Q4 %d\n",
            sum(wt_hit), sum(wt_hit & res$Q2WT_est < 0), sum(wt_hit & res$Q2WT_est > 0),
            sum(kd_hit), sum(kd_hit & res$Q2KD_est < 0), sum(kd_hit & res$Q2KD_est > 0),
            sum(uni_hit), uni_M, uni_F, sum(dep)))

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
## Label each genotype with its TOTAL number of sex-biased hits, so the
## male-biased bar (61) is not misread as the WT hit total (63). The union
## across genotypes is 64, of which 62 are male-biased - that is the
## "62 of 64" figure, a different denominator from either bar here.
bc$geno <- factor(bc$geno, c("WT cells", "KD cells"))
bc$cat  <- factor(bc$cat, c("Male-biased hit", "Female-biased hit", "n.s."))
tot_lab <- data.frame(
  geno = factor(c("WT cells", "KD cells"), c("WT cells", "KD cells")),
  n    = c(sum(wt_hit), sum(kd_hit)))
tot_lab$txt <- sprintf("%d sex-biased hits", tot_lab$n)

p <- ggplot(bc, aes(geno, n, fill = cat)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.4, size = 5) +
  geom_text(data = tot_lab, aes(geno, y = max(bc$n) * 1.16, label = txt),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = COL) + expand_limits(y = max(bc$n) * 1.24) +
  labs(title = "Sex-biased hits per genotype",
       subtitle = sprintf("%d drugs sex-biased in WT or KD; %d of those are male-biased",
                          sum(uni_hit), uni_M),
       x = NULL, y = "number of drugs") +
  base_thm + theme(legend.position = "bottom")
save2(p, "panelB_barchart_hitcounts", 7.0, 5.4)

## ---- Panel C: Egr1-dependence pie -------------------------------------------
## Denominator = the UNION of WT and KD sex-biased hits (64), not WT alone (63).
## Both thioguanine runs are retained as separate rows, so 64 = 63 unique
## compounds + the second thioguanine run. All 37 Q4-significant compounds fall
## inside this union, so the split is 37 dependent / 27 independent.
nd <- sum(uni_hit & dep); ni <- sum(uni_hit & !dep); tot <- nd + ni
cat(sprintf("  panel C denominator = WT-or-KD union: %d (%d dependent / %d independent)\n",
            tot, nd, ni))
pie <- data.frame(lab = c("Egr1-dependent\n(sig. interaction)", "Egr1-independent"),
                  n = c(nd, ni))
pie$pct <- round(100 * pie$n / tot)
p <- ggplot(pie, aes("", n, fill = lab)) +
  geom_col(width = 1, colour = "white") + coord_polar("y", start = 0) +
  geom_text(aes(label = sprintf("%d\n(%d%%)", n, pct)),
            position = position_stack(vjust = 0.5), size = 6, fontface = "bold") +
  scale_fill_manual(values = c("#3B4CC0", "grey75")) +
  labs(title = "Egr1 dependence of sex-biased hits",
       subtitle = sprintf("n = %d hits sex-biased in WT or KD", tot)) +
  theme_void(base_size = 16) +
  theme(plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 6)),
        legend.text   = element_text(size = 14), legend.title = element_blank())
save2(p, "panelC_pie_Egr1dependence", 8.0, 5.6)

## ---- Panel D: % Egr1-dependent by drug class --------------------------------
## Panel D's denominator is WT hits only (63) - NOT panel C's union (64).
## Compute it locally so the subtitle can never drift from the bars.
nd_wt  <- sum(wt_hit & dep)
tot_wt <- sum(wt_hit)
cat(sprintf("  panel D denominator = WT hits only: %d (%d dependent)\n", tot_wt, nd_wt))

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
       subtitle = sprintf("%d of %d WT sex-biased hits are Egr1-dependent", nd_wt, tot_wt),
       x = "% of WT sex-biased hits that are Egr1-dependent", y = NULL) +
  base_thm
save2(p, "panelD_pctEgr1dep_byClass", 10.5, 7.0)

## ---- Panel D2: class composition of ALL Egr1-dependent drugs ----------------
## Panel D above asks "what fraction of hits in each class are Egr1-dependent"
## and is therefore restricted to WT hits. This panel instead takes the full set
## of Egr1-dependent drugs (Q4 raw p <= 0.05) and shows how they distribute
## across drug classes - i.e. how many pathways the Egr1-dependent vulnerability
## spans. All Q4-significant compounds are sex-biased hits, so this set is the
## same 37 either way; what changes is that the denominator is the 37, not the
## hits per class.
dd <- res[dep, ]
cd <- as.data.frame(table(dd$Class), stringsAsFactors = FALSE)
names(cd) <- c("Class", "n")
cd <- cd[cd$n > 0, ]
cd$pct <- 100 * cd$n / sum(cd$n)
cd <- cd[order(cd$n, cd$Class, decreasing = c(FALSE, TRUE), method = "radix"), ]
cd$Class <- factor(cd$Class, levels = cd$Class)
cat(sprintf("  panel D2: %d Egr1-dependent drugs spanning %d classes\n",
            sum(cd$n), nrow(cd)))

p <- ggplot(cd, aes(n, Class)) +
  geom_col(fill = "#3B4CC0", width = 0.75) +
  geom_text(aes(label = sprintf("%d  (%.0f%%)", n, pct)),
            hjust = -0.12, size = 5, fontface = "bold") +
  scale_x_continuous(limits = c(0, max(cd$n) * 1.42),
                     breaks = scales::breaks_pretty(5)) +
  labs(title = "Egr1-dependent drug vulnerabilities span multiple classes",
       subtitle = sprintf("all %d Egr1-dependent drugs (Q4 raw p <= 0.05), across %d drug classes",
                          sum(cd$n), nrow(cd)),
       x = "number of Egr1-dependent drugs", y = NULL) +
  base_thm
save2(p, "panelD2_Egr1dependent_byClass", 10.5, 7.0)

## ---- Panel E: dumbbell, WT vs KD sex difference for all union hits ----------
## Grey = Egr1 WT, blue = Egr1 KD; arrow shows the WT->KD shift. Bold label =
## significant Egr1 x sex interaction (Q4 raw p <= 0.05). Duplicate compound
## names (THIOGUANINE appears twice) are disambiguated as "(run 1)/(run 2)".
lab <- res$Compound
dupn <- names(which(table(lab) > 1))
for (d in dupn) {
  ix <- which(lab == d)
  lab[ix] <- sprintf("%s (run %d)", d, seq_along(ix))
}
lab <- tools::toTitleCase(tolower(lab))

db <- data.frame(drug = lab, wt = res$Q2WT_est, kd = res$Q2KD_est,
                 sig = dep, stringsAsFactors = FALSE)[uni_hit, ]
db <- db[order(db$wt), ]
db$drug <- factor(db$drug, levels = db$drug)
long <- rbind(data.frame(drug = db$drug, z = db$wt, geno = "Egr1 WT"),
              data.frame(drug = db$drug, z = db$kd, geno = "Egr1 KD"))
long$geno <- factor(long$geno, c("Egr1 WT", "Egr1 KD"))

p <- ggplot(db) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_segment(aes(x = wt, xend = kd, y = drug, yend = drug, linewidth = sig),
               colour = "grey25", arrow = arrow(length = unit(0.10, "cm"), type = "closed")) +
  geom_point(data = long, aes(z, drug, colour = geno), size = 2.4) +
  scale_linewidth_manual(values = c("FALSE" = 0.25, "TRUE" = 0.75), guide = "none") +
  scale_colour_manual(values = c("Egr1 WT" = "grey55", "Egr1 KD" = "#1E90FF")) +
  labs(title = sprintf("Sex-biased drug-screen hits (n = %d): WT vs KD", nrow(db)),
       subtitle = sprintf("ordered by WT male-vs-female z\nbold = significant Egr1 x sex interaction (raw p <= 0.05, n = %d)", sum(db$sig)),
       x = "Male-vs-female sensitivity in WT (Q2), z units", y = NULL) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 8,
                       face = ifelse(db$sig, "bold", "plain"), colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        legend.position = "top", legend.title = element_blank(),
        panel.grid.minor = element_blank())
save2(p, "panelE_dumbbell_WTvsKD_n64", 7.5, 12.5)

cat("\nAll panels written to", normalizePath(OUT), "\n")
