# ============================================================
#  Sex difference (M/F AUC fold change) per genotype.
#  SEPARATE model per figure. Between-genotype comparison vs
#  Egr1 WT shown as the ACTUAL p-value (per Rosy), on a bracket.
# ============================================================
library(dplyr); library(ggplot2); library(emmeans); library(ggsignif); library(conflicted)
conflict_prefer("select","dplyr"); conflict_prefer("filter","dplyr")

# ---- Data -------------------------------------------------
df <- data.frame(
  condition = c(rep("EgrWT",4), rep("EgrKD",4), rep("Nrp1",8),
                rep("Hmox1",8), rep("Mast4",8), rep("Ptk2b",8)),
  sex = c(rep(c("M","F"),each=2), rep(c("M","F"),each=2),
          rep(c("M","F"),each=4), rep(c("M","F"),each=4),
          rep(c("M","F"),each=4), rep(c("M","F"),each=4)),
  value = c(1843,1806, 691.1,726.6,
            1612,1477, 1014,1853,
            1468,1749,1485,1549, 1065,1977,1054,1167,
            1575,1719,1749,1663, 1376,1386,1104,683.3,
            1583,1446,1758,1827, 1414,1443,939.9,1327,
            1815,1519,1430,1839, 1069,751,769.1,1031)
)
df$log_value <- log2(df$value)
df$condition <- factor(df$condition, levels=c("EgrWT","EgrKD","Nrp1","Hmox1","Mast4","Ptk2b"))
df$sex       <- factor(df$sex, levels=c("M","F"))
cond_levels  <- levels(df$condition)

outdir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# helper: format p the way Rosy asked (p = 0.054, and a sensible floor for tiny values)
pfmt <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "f", digits = 3)))

run_panel <- function(conds, plot_title, file_name, w, ref = "EgrWT") {
  dd <- df |> filter(condition %in% conds) |>
    mutate(condition = droplevels(factor(condition, levels = cond_levels)),
           sex       = factor(sex, levels = c("M","F")))
  
  fit <- lm(log_value ~ condition * sex, data = dd)
  cat("\n", plot_title, "  (residual df = ", df.residual(fit), ")\n", sep = "")
  
  fc <- contrast(emmeans(fit, ~ sex | condition), "pairwise", by = "condition") |>
    summary(infer = TRUE, adjust = "none") |> as.data.frame() |>
    mutate(fold_change = 2^estimate, fold_lower = 2^lower.CL, fold_upper = 2^upper.CL)
  
  it <- contrast(emmeans(fit, ~ sex * condition),
                 interaction = c("pairwise","trt.vs.ctrl")) |>
    summary(infer = TRUE, adjust = "none") |> as.data.frame()
  cc <- names(it)[vapply(it, function(x) any(grepl(ref, as.character(x))), logical(1))][1]
  it$comparison <- as.character(it[[cc]])
  it <- it |>
    mutate(target = trimws(sub("-.*","",comparison)),
           RoR = 2^estimate, RoR_lo = 2^lower.CL, RoR_hi = 2^upper.CL,
           p_adj = p.adjust(p.value, "fdr"),
           plabel = pfmt(p.value))          # <-- actual p-value text (raw, two-sided)
  print(it[, c("comparison","RoR","p.value","p_adj")], digits = 3)
  
  comps <- it |> filter(target %in% conds)
  ymax  <- max(fc$fold_upper)
  yp    <- ymax * 1.04 * 1.10 ^ (seq_len(nrow(comps)) - 1)
  
  p <- ggplot(fc, aes(condition, fold_change)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_errorbar(aes(ymin = fold_lower, ymax = fold_upper), width = 0.15, linewidth = 0.5) +
    geom_point(size = 3.2, colour = "#5E3C99") +
    scale_y_continuous(trans = "log2", breaks = c(0.5,0.75,1,1.5,2,3,4),
                       expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = NULL, y = "Male / Female AUC  (fold change)", title = plot_title) +
    theme_classic(base_size = 13)
  
  if (nrow(comps) > 0)
    p <- p + geom_signif(annotations = comps$plabel,        # <-- p = X.XXX instead of */ns
                         xmin = rep(ref, nrow(comps)),
                         xmax = comps$target, y_position = yp,
                         tip_length = 0.005, textsize = 3.4)
  
  ggsave(file.path(outdir, file_name), plot = p, width = w, height = 4.5,
         units = "in", device = cairo_pdf)
  p
}

p_egr <- run_panel(c("EgrWT","EgrKD"),
                   "AUC: Egr1 WT vs Egr1 KD",
                   "sex_foldchange_EgrWT_vs_EgrKD_5.pdf", w = 4)

p_kd  <- run_panel(c("EgrWT","Nrp1","Hmox1","Mast4","Ptk2b"),
                   "AUC: Egr1 WT vs other knockdowns",
                   "sex_foldchange_EgrWT_vs_otherKDs_5.pdf", w = 7)

print(p_egr); print(p_kd)