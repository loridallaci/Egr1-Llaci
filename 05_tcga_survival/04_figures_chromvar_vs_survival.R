# =============================================================================
# 05_tcga_survival: 04 - ChromVAR Activity vs Survival Figures
# =============================================================================
# Description:
#   Integrates RENIN ChromVAR motif activity scores (male - female) with
#   multivariate Cox survival results. Generates:
#     1. Scatter: ChromVAR difference vs male HR
#     2. Scatter: ChromVAR difference vs female HR
#     3. Quadrant plot: ChromVAR diff vs HR sex difference
#     4. Heatmap: top 30 sex-differential motifs (ChromVAR + HR)
#     5. Bar plot: top 20 motifs with HR significance annotation
#     6. Forest plots: per gene for male, female, and all patients
#
# Run order:
#   1. 01_load_and_prepare_tcga_data.R
#   2. 02_multivariate_cox_regression.R
#   3. 03_kaplan_meier_plots.R
#   4. 04_figures_chromvar_vs_survival.R  <-- this script
#
# Input:
#   - RENINmotif_multivariate_male_and_female.csv
#   - cox_results_MALES/FEMALES/ALL.csv
#
# Output:
#   - chromvar_diff_vs_male_HR.pdf
#   - chromvar_diff_vs_female_HR.pdf
#   - chromvar_vs_HR_quadrant.pdf
#   - chromvar_HR_combined_heatmap.pdf
#   - chromvar_top20_with_HR_annotation.pdf
#   - TCGA_[sex]_multivariatePlot_[gene].pdf  (forest plots)
#   - gene_HR_summary_by_sex.csv
# =============================================================================

# --- Libraries ----------------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(patchwork)
library(survival)

source("05_tcga_survival/utils.R")

# original (author's machine): "/home/lllaci/data/tcga_survival_results"
output_dir <- "output/tcga_survival_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load integrated data -----------------------------------------------------

chromvar_diff      <- read.csv(
  file.path(output_dir, "RENINmotif_multivariate_male_and_female.csv"),
  row.names = 1
)
survival_table_male   <- read.csv(file.path(output_dir, "cox_results_MALES.csv"))
survival_table_female <- read.csv(file.path(output_dir, "cox_results_FEMALES.csv"))
survival_table_all    <- read.csv(file.path(output_dir, "cox_results_ALL.csv"))

# Keep only first instance of each TF (handles duplicates from "::" splitting)
chromvar_diff <- chromvar_diff %>% distinct(TF_name, .keep_all = TRUE)

# --- Classify ChromVAR and HR significance ------------------------------------

chromvar_diff <- chromvar_diff %>%
  mutate(
    abs_diff          = abs(diff),
    chromvar_category = case_when(
      diff >  1 ~ "Male-enriched",
      diff < -1 ~ "Female-enriched",
      TRUE      ~ "No difference"
    ),
    # HR significance: HR > 1.1 and p <= 0.05
    male_HR_sig = ifelse(
      !is.na(Expression.pvalue_M) & !is.na(Expression.coef_M) &
        Expression.pvalue_M <= 0.05 & Expression.coef_M > 1.1,
      "Significant", "Not significant"
    ),
    female_HR_sig = ifelse(
      !is.na(Expression.pvalue_F) & !is.na(Expression.coef_F) &
        Expression.pvalue_F <= 0.05 & Expression.coef_F > 1.1,
      "Significant", "Not significant"
    ),
    combined_category = case_when(
      chromvar_category == "Male-enriched"   & male_HR_sig   == "Significant" ~
        "Male-enriched + Male HR sig",
      chromvar_category == "Female-enriched" & female_HR_sig == "Significant" ~
        "Female-enriched + Female HR sig",
      chromvar_category == "Male-enriched"   ~ "Male-enriched only",
      chromvar_category == "Female-enriched" ~ "Female-enriched only",
      TRUE ~ "No ChromVAR difference"
    ),
    HR_diff = Expression.coef_M - Expression.coef_F
  )

# =============================================================================
# FIGURE 1: ChromVAR diff vs Male HR
# =============================================================================

sig_male <- chromvar_diff %>%
  filter(male_HR_sig == "Significant" & abs(diff) > 1)

p1 <- ggplot(chromvar_diff, aes(x = diff, y = Expression.coef_M)) +
  geom_hline(yintercept = 1,      linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "grey50") +
  geom_point(aes(color = male_HR_sig, size = abs_diff), alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "#D32F2F",
                                "Not significant" = "grey60")) +
  scale_size_continuous(range = c(2, 6)) +
  { if (nrow(sig_male) > 0)
    geom_text_repel(data = sig_male, aes(label = TF_name),
                    size = 3, max.overlaps = 20, color = "black") } +
  labs(title    = "ChromVAR Sex Difference vs Male Survival HR",
       subtitle = "Each point is a transcription factor motif",
       x        = "ChromVAR difference (male - female)",
       y        = "Hazard Ratio in Males",
       color    = "Male HR significance",
       size     = "Absolute ChromVAR diff") +
  theme_classic(base_size = 14)

ggsave(file.path(output_dir, "chromvar_diff_vs_male_HR.pdf"),
       p1, width = 12, height = 8)

# =============================================================================
# FIGURE 2: ChromVAR diff vs Female HR
# =============================================================================

sig_female <- chromvar_diff %>%
  filter(female_HR_sig == "Significant" & abs(diff) > 1)

p2 <- ggplot(chromvar_diff, aes(x = diff, y = Expression.coef_F)) +
  geom_hline(yintercept = 1,       linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_point(aes(color = female_HR_sig, size = abs_diff), alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "#D32F2F",
                                "Not significant" = "grey60")) +
  scale_size_continuous(range = c(2, 6)) +
  { if (nrow(sig_female) > 0)
    geom_text_repel(data = sig_female, aes(label = TF_name),
                    size = 3, max.overlaps = 20, color = "black") } +
  labs(title    = "ChromVAR Sex Difference vs Female Survival HR",
       subtitle = "Each point is a transcription factor motif",
       x        = "ChromVAR difference (male - female)",
       y        = "Hazard Ratio in Females",
       color    = "Female HR significance",
       size     = "Absolute ChromVAR diff") +
  theme_classic(base_size = 14)

ggsave(file.path(output_dir, "chromvar_diff_vs_female_HR.pdf"),
       p2, width = 12, height = 8)

# =============================================================================
# FIGURE 3: Quadrant plot - ChromVAR diff vs HR sex difference
# =============================================================================

interesting <- chromvar_diff %>%
  filter((male_HR_sig == "Significant" | female_HR_sig == "Significant") &
           abs(diff) > 1)

p3 <- ggplot(chromvar_diff, aes(x = diff, y = HR_diff)) +
  geom_hline(yintercept = 0, color = "grey30") +
  geom_vline(xintercept = 0, color = "grey30") +
  geom_point(aes(color = combined_category, size = abs_diff), alpha = 0.7) +
  scale_color_manual(values = c(
    "Male-enriched + Male HR sig"     = "#D32F2F",
    "Female-enriched + Female HR sig" = "#E24A90",
    "Male-enriched only"              = "#4A90E2",
    "Female-enriched only"            = "#90CAF9",
    "No ChromVAR difference"          = "grey70"
  )) +
  scale_size_continuous(range = c(2, 6)) +
  { if (nrow(interesting) > 0)
    geom_text_repel(data = interesting, aes(label = TF_name),
                    size = 3, max.overlaps = 15, color = "black") } +
  labs(title    = "ChromVAR Sex Difference vs HR Sex Difference",
       x        = "ChromVAR difference (male - female)",
       y        = "HR difference (male HR - female HR)",
       color    = "Category",
       size     = "Absolute ChromVAR diff") +
  theme_classic(base_size = 14)

ggsave(file.path(output_dir, "chromvar_vs_HR_quadrant.pdf"),
       p3, width = 13, height = 9)

# =============================================================================
# FIGURE 4: Heatmap - top 30 sex-differential motifs
# =============================================================================

top30 <- chromvar_diff %>%
  arrange(desc(abs_diff)) %>%
  head(30)

heatmap_data <- data.frame(
  ChromVAR_Female = top30$female_mean,
  ChromVAR_Male   = top30$male_mean,
  ChromVAR_Diff   = top30$diff,
  log2HR_Male     = log2(top30$Expression.coef_M),
  log2HR_Female   = log2(top30$Expression.coef_F),
  row.names       = paste0(top30$TF_name, " (", top30$motif, ")")
)

annotation_col <- data.frame(
  Male_HR_sig   = ifelse(top30$male_HR_sig   == "Significant", "Yes", "No"),
  Female_HR_sig = ifelse(top30$female_HR_sig == "Significant", "Yes", "No"),
  row.names     = rownames(heatmap_data)
)

ann_colors <- list(
  Male_HR_sig   = c("Yes" = "#D32F2F", "No" = "grey90"),
  Female_HR_sig = c("Yes" = "#E24A90", "No" = "grey90")
)

pheatmap(
  t(as.matrix(heatmap_data)),
  cluster_rows    = FALSE,
  cluster_cols    = TRUE,
  scale           = "row",
  color           = colorRampPalette(c("#4A90E2", "white", "#D32F2F"))(100),
  main            = "Top 30 Sex-Differential Motifs: ChromVAR + Survival HR",
  fontsize_row    = 10,
  fontsize_col    = 7,
  angle_col       = 45,
  cellheight      = 25,
  cellwidth       = 15,
  annotation_col  = annotation_col,
  annotation_colors = ann_colors,
  filename        = file.path(output_dir, "chromvar_HR_combined_heatmap.pdf"),
  width           = 16,
  height          = 8
)

# =============================================================================
# FIGURE 5: Bar plot - top 20 motifs with HR annotation
# =============================================================================

top20 <- chromvar_diff %>%
  arrange(desc(abs_diff)) %>%
  head(20) %>%
  mutate(
    TF_label     = paste0(TF_name, " (", motif, ")"),
    HR_indicator = case_when(
      male_HR_sig == "Significant" & female_HR_sig == "Significant" ~ "Both sig",
      male_HR_sig == "Significant"   ~ "Male sig",
      female_HR_sig == "Significant" ~ "Female sig",
      TRUE ~ "Neither sig"
    )
  )

p4 <- ggplot(top20, aes(x = reorder(TF_label, diff), y = diff)) +
  geom_bar(aes(fill = chromvar_category), stat = "identity") +
  geom_point(aes(shape = HR_indicator, color = HR_indicator),
             size = 4, position = position_nudge(y = 0.1)) +
  coord_flip() +
  scale_fill_manual(values = c("Male-enriched"   = "#4A90E2",
                               "Female-enriched" = "#E24A90",
                               "No difference"   = "grey70")) +
  scale_shape_manual(values = c("Both sig"    = 16, "Male sig" = 17,
                                "Female sig"  = 15, "Neither sig" = 1)) +
  scale_color_manual(values = c("Both sig"    = "#D32F2F", "Male sig"    = "#F44336",
                                "Female sig"  = "#E91E63", "Neither sig" = "grey50")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title    = "Top 20 Sex-Differential Motifs",
       subtitle = "Bar = ChromVAR difference; Points = Survival HR significance",
       x        = "",
       y        = "ChromVAR difference (male - female)",
       fill     = "ChromVAR category",
       shape    = "HR significance",
       color    = "HR significance") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.y     = element_text(size = 9))

ggsave(file.path(output_dir, "chromvar_top20_with_HR_annotation.pdf"),
       p4, width = 12, height = 10)

# =============================================================================
# FIGURE 6: Forest plots per gene (male, female, all)
# =============================================================================

message("Generating forest plots...")
for (gene_to_plot in geneID_included) {

  for (group_info in list(
    list(table = survival_table_male,   label = "MALES"),
    list(table = survival_table_female, label = "FEMALES"),
    list(table = survival_table_all,    label = "ALL")
  )) {

    if (!gene_to_plot %in% group_info$table$Gene) next

    p <- create_forest_plot(
      gene_name     = gene_to_plot,
      survival_data = group_info$table,
      title_suffix  = paste0(" - ", group_info$label)
    )

    ggsave(
      filename = file.path(output_dir,
                           paste0("TCGA_", group_info$label,
                                  "_multivariatePlot_", gene_to_plot, ".pdf")),
      plot   = p,
      width  = 7,
      height = 6
    )
  }
}

message("All figures saved to: ", output_dir)
