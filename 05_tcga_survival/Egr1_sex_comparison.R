# =============================================================================
# EGR1 Expression Comparison: Male vs Female (IDH wild-type GBM)
# =============================================================================
# Run after: 01_load_and_prepare_tcga_data.R
# Uses: tcga_pheno (contains EGR1 column and Gender column)
# =============================================================================

library(dplyr)
library(ggplot2)
library(ggpubr)

if (!exists("tcga_pheno")) stop("Run 01_load_and_prepare_tcga_data.R first.")

output_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"

gene      <- "EGR1"   # column name in tcga_pheno (uppercased by script 01)
plot_file <- file.path(output_dir, "EGR1_expression_TCGA_IDHwt_male_vs_female.pdf")

# --- Sanity check -------------------------------------------------------------

if (!gene %in% colnames(tcga_pheno)) {
  stop(sprintf("'%s' not found in tcga_pheno. Check geneID_included.", gene))
}

# --- Subset and label ---------------------------------------------------------

dat <- tcga_pheno %>%
  filter(!is.na(Gender), !is.na(.data[[gene]])) %>%
  mutate(Sex = factor(Gender, levels = c("Male", "Female")))

# --- Descriptive statistics ---------------------------------------------------

cat("=== EGR1 Expression by Sex ===\n")
desc <- dat %>%
  group_by(Sex) %>%
  summarise(
    n      = n(),
    mean   = round(mean(.data[[gene]]), 3),
    median = round(median(.data[[gene]]), 3),
    sd     = round(sd(.data[[gene]]), 3),
    IQR    = round(IQR(.data[[gene]]), 3),
    .groups = "drop"
  )
print(desc)

# --- Statistical tests --------------------------------------------------------

male_vals   <- dat %>% filter(Sex == "Male")   %>% pull(gene)
female_vals <- dat %>% filter(Sex == "Female") %>% pull(gene)

# Wilcoxon rank-sum (Mann-Whitney) — non-parametric, robust for expression data
wt  <- wilcox.test(male_vals, female_vals, exact = FALSE)

# t-test — for reference / normality assumption
tt  <- t.test(male_vals, female_vals)

# Normality check (Shapiro-Wilk, max 5000 samples)
sw_male   <- if (length(male_vals)   <= 5000) shapiro.test(male_vals)   else NULL
sw_female <- if (length(female_vals) <= 5000) shapiro.test(female_vals) else NULL

cat("\n=== Statistical Tests ===\n")
cat(sprintf("Wilcoxon rank-sum test:  W = %.1f,  p = %.4g\n", wt$statistic, wt$p.value))
cat(sprintf("Welch t-test:            t = %.3f,  p = %.4g\n", tt$statistic, tt$p.value))
if (!is.null(sw_male))   cat(sprintf("Shapiro-Wilk (Male):     W = %.4f, p = %.4g\n", sw_male$statistic,   sw_male$p.value))
if (!is.null(sw_female)) cat(sprintf("Shapiro-Wilk (Female):   W = %.4f, p = %.4g\n", sw_female$statistic, sw_female$p.value))

# Effect size: rank-biserial correlation (r = 1 - 2W/(n1*n2))
r_effect <- 1 - (2 * wt$statistic) / (length(male_vals) * length(female_vals))
cat(sprintf("Effect size (rank-biserial r): %.3f\n", r_effect))

# --- Plot: Violin + boxplot + jitter ------------------------------------------

p <- ggplot(dat, aes(x = Sex, y = .data[[gene]], fill = Sex)) +
  geom_violin(trim = FALSE, alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Sex), width = 0.08, size = 1.2, alpha = 0.5) +
  stat_compare_means(
    method       = "wilcox.test",
    label        = "p.format",
    label.x      = 1.5,
    label.y      = max(dat[[gene]], na.rm = TRUE) * 1.05,
    size         = 4.5,
    fontface     = "bold"
  ) +
  scale_fill_manual(values  = c("Male" = "#4E79A7", "Female" = "#F28E2B")) +
  scale_color_manual(values = c("Male" = "#4E79A7", "Female" = "#F28E2B")) +
  labs(
    title    = sprintf("%s Expression in IDH Wild-Type GBM", gene),
    subtitle = sprintf("Male n=%d  |  Female n=%d  |  Wilcoxon p=%.3g",
                       length(male_vals), length(female_vals), wt$p.value),
    x        = NULL,
    y        = sprintf("%s Expression (log2 CPM)", gene)
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40")
  )

print(p)
ggsave(plot_file, plot = p, width = 5, height = 6, dpi = 300)
cat(sprintf("\nPlot saved to: %s\n", plot_file))