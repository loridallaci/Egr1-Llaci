## =====================================================================
## Forest plot per RENIN motif list (male-biased / female-biased RENIN
## motif TFs): each TF's gene-expression HR (95% CI) in MALE vs FEMALE
## TCGA patients, on one plot. TFs ordered by the multivariate overall
## p-value (most significant on top). One PDF per RENIN motif set.
## Reads the COMBINED_allTFs tables made by motif_multivariate_TCGA.R.
## =====================================================================
suppressMessages({library(dplyr); library(ggplot2)})

base   <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/05_tcga_survival/motif_multivariate_TCGA"

make_forest <- function(set) {
  d <- read.csv(file.path(base, set, paste0("cox_wide_", set, "_COMBINED_allTFs.csv")), check.names = FALSE)
  long <- bind_rows(
    data.frame(TF = d$TF, Patient = "Male TCGA",
               HR = d[["Gene expression HR male"]],
               Lo = d[["Gene expression HR lower CI male"]],
               Hi = d[["Gene expression HR upper CI male"]],
               p  = d[["Gene expression p value male"]]),
    data.frame(TF = d$TF, Patient = "Female TCGA",
               HR = d[["Gene expression HR female"]],
               Lo = d[["Gene expression HR lower CI female"]],
               Hi = d[["Gene expression HR upper CI female"]],
               p  = d[["Gene expression p value female"]]))
  long <- long[is.finite(long$HR) & is.finite(long$Lo) & is.finite(long$Hi), ]
  long$Patient <- factor(long$Patient, levels = c("Male TCGA", "Female TCGA"))
  ## order TFs by the own-sex multivariate (overall-model) p-value:
  ## most significant (smallest p) on TOP. Male motifs -> male p, female motifs -> female p.
  own_p_col <- if (set == "MaleMotifs") "Multivariate p value male" else "Multivariate p value female"
  ord <- d$TF[order(d[[own_p_col]], decreasing = TRUE)]   # largest p first (bottom) -> smallest last (top)
  long$TF  <- factor(long$TF, levels = ord)
  long$sig <- ifelse(!is.na(long$p) & long$p < 0.05, "p < 0.05", "ns")
  set_label <- if (set == "MaleMotifs") "Male-biased RENIN motifs" else "Female-biased RENIN motifs"
  p <- ggplot(long, aes(HR, TF, color = Patient)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = Lo, xmax = Hi), height = 0.3, position = position_dodge(width = 0.6)) +
    geom_point(size = 2.8, position = position_dodge(width = 0.6)) +   # all filled
    geom_text(aes(x = Hi, label = ifelse(sig == "p < 0.05", "*", "")), # * just past each bar's upper CI, same row
              position = position_dodge(width = 0.6), hjust = -0.35, vjust = 0.78, size = 7, show.legend = FALSE) +
    scale_color_manual(values = c("Male TCGA" = "#1E90FF", "Female TCGA" = "#FF69B4"), name = NULL) +
    scale_x_log10() +
    scale_y_discrete(expand = expansion(add = 0.7)) +   # padding so top/bottom dodged points aren't clipped
    labs(x = "RENIN motif TF gene-expression HR (95% CI, log scale)",
         y = "RENIN motif TF",
         title = paste0(set_label, " - multivariate Cox in TCGA GBM (M vs F patients)"),
         subtitle = "ordered by multivariate overall p-value (most significant on top)",
         caption = "* = gene-expression HR p < 0.05 (Wald test on the Expression term)") +
    theme_bw() + theme(legend.position = "top")
  out <- file.path(base, set, paste0(set, "_RENINmotifs_forest_M_and_F.pdf"))
  ggsave(out, p, width = 8, height = 8)
  cat(set, ":", nlevels(long$TF), "TFs ->", out, "\n")
}

for (s in c("MaleMotifs", "FemaleMotifs")) make_forest(s)
cat("Done.\n")
