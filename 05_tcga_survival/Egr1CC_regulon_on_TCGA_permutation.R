# =============================================================================
# NEGATIVE CONTROL + EMPIRICAL P-VALUES: scrambled gene sets vs TCGA survival
# -----------------------------------------------------------------------------
# Mirrors Egr1CC_regulon_on_TCGA.R, replacing each real regulon with a
# SIZE-MATCHED random gene set. Same z-score -> Cox pipeline.
#
# PART A : one scrambled draw per group  -> forest plot (same style as real)
# PART B : N random draws per group      -> permutation null + EMPIRICAL P
#          one-sided p  = fraction of random sets as harmful as the real HR
#          two-sided p  = fraction of random sets as far from HR=1 (log scale)
#
# Run 01_load_and_prepare_tcga_data.R first (tcga_pheno objects must exist).
# =============================================================================

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

if (!exists("tcga_pheno")) stop("Run 01_load_and_prepare_tcga_data.R first.")

set.seed(42)   # reproducible

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/enrichment/Egr1CCvsKD/"
output_dir <- "output/enrichment/Egr1CCvsKD/"

# "random_genome" = random genes from the whole TCGA matrix (standard control)
# "reshuffle_within_regulons" = random genes from the pooled regulon universe
scramble_mode <- "random_genome"

n_perm <- 1000   # permutations per group. Raise to 10000 to resolve p < 0.001
# (10000 takes a few minutes; min reportable p = 1/(n_perm+1)).

# =============================================================================
# STEP 1: Read real gene sets (for their SIZES only)
# =============================================================================

overlap_Male_KD   <- read.csv(file.path(output_dir, "overlap_MaleCC_vs_MaleKD_enriched.csv"))
overlap_Male_WT   <- read.csv(file.path(output_dir, "overlap_MaleCC_vs_MaleWT_enriched.csv"))
overlap_Female_KD <- read.csv(file.path(output_dir, "overlap_FemaleCC_vs_FemaleKD_enriched.csv"))
overlap_Female_WT <- read.csv(file.path(output_dir, "overlap_FemaleCC_vs_FemaleWT_enriched.csv"))

geneset_Male_KD   <- unique(toupper(overlap_Male_KD$SYMBOL))
geneset_Male_WT   <- unique(toupper(overlap_Male_WT$SYMBOL))
geneset_Female_KD <- unique(toupper(overlap_Female_KD$SYMBOL))
geneset_Female_WT <- unique(toupper(overlap_Female_WT$SYMBOL))

# =============================================================================
# STEP 2: Gene pool to scramble FROM (numeric expression cols, no metadata)
# =============================================================================

non_gene_cols <- c("survival", "status", "Recurrence", "Age", "Subtype",
                   "MGMT_status", "Sex", "sample", "SampleID", "barcode",
                   "patient", "OS", "OS_days", "OS_event", "vital_status")

gene_pool <- setdiff(colnames(tcga_pheno), non_gene_cols)
gene_pool <- gene_pool[vapply(tcga_pheno[gene_pool], is.numeric, logical(1))]
gene_pool <- gene_pool[!grepl("_binary$|_score$", gene_pool)]
cat("Expression gene pool:", length(gene_pool), "genes\n")

if (scramble_mode == "random_genome") {
  scramble_pool <- gene_pool
} else if (scramble_mode == "reshuffle_within_regulons") {
  scramble_pool <- intersect(
    unique(c(geneset_Male_WT, geneset_Male_KD,
             geneset_Female_WT, geneset_Female_KD)), gene_pool)
} else {
  stop("scramble_mode must be 'random_genome' or 'reshuffle_within_regulons'")
}
cat("Scramble pool (", scramble_mode, "):", length(scramble_pool), "genes\n")

set_sizes <- c(
  MaleCC_WT_enriched   = sum(geneset_Male_WT   %in% gene_pool),
  MaleCC_KD_enriched   = sum(geneset_Male_KD   %in% gene_pool),
  FemaleCC_WT_enriched = sum(geneset_Female_WT %in% gene_pool),
  FemaleCC_KD_enriched = sum(geneset_Female_KD %in% gene_pool)
)
print(set_sizes)

runs <- list(
  list(geneset = "MaleCC_WT_enriched",   dataset = "Males"),
  list(geneset = "MaleCC_KD_enriched",   dataset = "Males"),
  list(geneset = "FemaleCC_WT_enriched", dataset = "Females"),
  list(geneset = "FemaleCC_KD_enriched", dataset = "Females")
)
datasets <- list(All = tcga_pheno, Males = tcga_pheno_male, Females = tcga_pheno_female)

# =============================================================================
# STEP 3: z-score gene-set Cox (identical math to the real script)
# =============================================================================

run_zscore_survival_scrambled <- function(pheno, geneset, dataset_name, geneset_name) {
  genes_available <- geneset[geneset %in% colnames(pheno)]
  if (length(genes_available) < 2) {
    cat("  Skipping", geneset_name, "x", dataset_name, "- <2 genes\n"); return(NULL)
  }
  gene_mat_z <- scale(as.matrix(pheno[, genes_available, drop = FALSE]))
  pheno$geneset_score <- rowMeans(gene_mat_z, na.rm = TRUE)
  
  cox_res <- coxph(
    Surv(survival, status) ~ geneset_score + Recurrence + Age + Subtype + MGMT_status,
    data = pheno)
  cox_sum <- summary(cox_res)
  
  hr    <- round(cox_sum$coefficients["geneset_score", "exp(coef)"], 3)
  lower <- round(cox_sum$conf.int["geneset_score",    "lower .95"],  3)
  upper <- round(cox_sum$conf.int["geneset_score",    "upper .95"],  3)
  pval  <- signif(cox_sum$coefficients["geneset_score", "Pr(>|z|)"], 3)
  
  cat(sprintf("  SCRAMBLED %-22s x %-8s HR=%.3f (%.3f-%.3f) p=%s\n",
              geneset_name, dataset_name, hr, lower, upper, pval))
  data.frame(Geneset = geneset_name, Dataset = dataset_name,
             N_genes = length(genes_available),
             HR = hr, Lower95 = lower, Upper95 = upper, Pvalue = pval,
             stringsAsFactors = FALSE)
}

# =============================================================================
# PART A: ONE scrambled draw per group -> forest plot
# =============================================================================

scrambled_sets <- lapply(set_sizes, function(n) sample(scramble_pool, n))

scrambled_results <- list()
for (r in runs) {
  res <- tryCatch(
    run_zscore_survival_scrambled(datasets[[r$dataset]], scrambled_sets[[r$geneset]],
                                  r$dataset, r$geneset),
    error = function(e) { warning(e$message); NULL })
  if (!is.null(res)) scrambled_results[[paste0(r$geneset, "_", r$dataset)]] <- res
}

scrambled_df <- bind_rows(scrambled_results)
print(scrambled_df)
write.csv(scrambled_df, file.path(output_dir, "TCGA_zscore_survival_SCRAMBLED.csv"),
          row.names = FALSE)

scrambled_df$Label <- paste0("SCRAMBLED ", scrambled_df$Geneset,
                             "\n(", scrambled_df$Dataset, ")")
scrambled_df$Label <- factor(scrambled_df$Label, levels = rev(scrambled_df$Label))
scrambled_df$Significant <- ifelse(scrambled_df$Pvalue < 0.05, "p < 0.05", "p >= 0.05")

p_forest <- ggplot(scrambled_df, aes(x = HR, y = Label, color = Significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = Lower95, xmax = Upper95), height = 0.2, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(x = max(Upper95) + 1, label = paste0("p = ", signif(Pvalue, 3))),
            hjust = 0, size = 3.5, color = "black") +
  scale_color_manual(values = c("p < 0.05" = "#E41A1C", "p >= 0.05" = "grey40")) +
  labs(x = "Hazard Ratio (95% CI)", y = NULL, color = NULL,
       title = "Scrambled Gene Set Survival Analysis - TCGA GBM (negative control)") +
  scale_x_continuous(limits = c(0, max(scrambled_df$Upper95) + 4)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 10), legend.position = "top",
        panel.grid.minor = element_blank())
print(p_forest)
ggsave(file.path(output_dir, "forest_plot_genesets_SCRAMBLED.pdf"),
       plot = p_forest, width = 8, height = 5)

# =============================================================================
# PART B: permutation null  ->  EMPIRICAL P-VALUES
# =============================================================================

# real HRs (and Cox p) from the original analysis
real_csv <- file.path(output_dir, "TCGA_zscore_survival_summary.csv")
if (file.exists(real_csv)) {
  real_df <- read.csv(real_csv)
} else {
  real_df <- data.frame(
    Geneset = c("MaleCC_WT_enriched","MaleCC_KD_enriched",
                "FemaleCC_WT_enriched","FemaleCC_KD_enriched"),
    Dataset = c("Males","Males","Females","Females"),
    HR      = c(2.734, 1.370, 2.358, 2.744),
    Pvalue  = c(0.00727, 0.614, 0.146, 0.135))
}

# one random-set HR
perm_one <- function(pheno, n_genes) {
  genes <- sample(scramble_pool[scramble_pool %in% colnames(pheno)], n_genes)
  pheno$geneset_score <- rowMeans(
    scale(as.matrix(pheno[, genes, drop = FALSE])), na.rm = TRUE)
  fit <- tryCatch(
    coxph(Surv(survival, status) ~ geneset_score + Recurrence + Age + Subtype + MGMT_status,
          data = pheno),
    error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  unname(summary(fit)$coefficients["geneset_score", "exp(coef)"])
}

null_summary <- list()
null_long    <- list()

for (r in runs) {
  key      <- paste0(r$geneset, "_", r$dataset)
  n_genes  <- set_sizes[[r$geneset]]
  pheno    <- datasets[[r$dataset]]
  row_idx  <- which(real_df$Geneset == r$geneset & real_df$Dataset == r$dataset)[1]
  real_hr  <- real_df$HR[row_idx]
  real_cox_p <- if ("Pvalue" %in% names(real_df)) real_df$Pvalue[row_idx] else NA
  
  null_hr <- replicate(n_perm, perm_one(pheno, n_genes))
  null_hr <- null_hr[is.finite(null_hr)]
  n_valid <- length(null_hr)
  
  # one-sided: as extreme as the real HR IN ITS OBSERVED DIRECTION
  # (+1 numerator/denominator = standard correction, p never exactly 0)
  if (real_hr >= 1) {
    emp_p_one <- (sum(null_hr >= real_hr) + 1) / (n_valid + 1)
  } else {
    emp_p_one <- (sum(null_hr <= real_hr) + 1) / (n_valid + 1)
  }
  # two-sided: as far from HR = 1 on the log scale
  emp_p_two <- (sum(abs(log(null_hr)) >= abs(log(real_hr))) + 1) / (n_valid + 1)
  
  cat(sprintf("%-22s x %-8s | real HR=%.3f | null median=%.3f | 1-sided p=%.4f | 2-sided p=%.4f\n",
              r$geneset, r$dataset, real_hr, median(null_hr), emp_p_one, emp_p_two))
  
  null_summary[[key]] <- data.frame(
    Geneset = r$geneset, Dataset = r$dataset,
    Real_HR = real_hr, Real_Cox_P = real_cox_p,
    Null_median_HR = round(median(null_hr), 3),
    Empirical_P_onesided = round(emp_p_one, 4),
    Empirical_P_twosided = round(emp_p_two, 4),
    N_perm = n_valid)
  null_long[[key]] <- data.frame(
    Pairing = paste0(r$geneset, "\n(", r$dataset, ")"),
    HR = null_hr, Real_HR = real_hr)
}

null_summary_df <- bind_rows(null_summary)
print(null_summary_df)
write.csv(null_summary_df, file.path(output_dir, "TCGA_zscore_PERMUTATION_null.csv"),
          row.names = FALSE)

# ---- Null distribution plot, annotated with the empirical p-values ----------
null_long_df <- bind_rows(null_long)

pval_labels <- null_summary_df %>%
  mutate(Pairing = paste0(Geneset, "\n(", Dataset, ")"),
         label = paste0("1-sided p = ", signif(Empirical_P_onesided, 3),
                        "\n2-sided p = ", signif(Empirical_P_twosided, 3)))

p_null <- ggplot(null_long_df, aes(x = HR)) +
  geom_histogram(bins = 40, fill = "grey70", color = "white") +
  geom_vline(aes(xintercept = Real_HR), color = "#E41A1C", linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(data = pval_labels, aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.3, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ Pairing, scales = "free", ncol = 2) +
  labs(x = "Hazard Ratio of random gene sets", y = "Count",
       title = paste0("Permutation null vs real regulon HR  (",
                      n_perm, " random sets, red = real)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        strip.text = element_text(size = 9))
print(p_null)
ggsave(file.path(output_dir, "permutation_null_distributions.pdf"),
       plot = p_null, width = 9, height = 6)

cat("\nDone. Empirical p-values in TCGA_zscore_PERMUTATION_null.csv\n")





# =============================================================================
# PART C: Forest-plot representation of the 1000-permutation test
# -----------------------------------------------------------------------------
# Run AFTER Egr1CC_regulon_SCRAMBLED_control_pvalues.R, in the SAME R session
# (uses null_long, null_summary_df, and output_dir from its PART B).
#
# Each row summarizes ALL 1000 random gene sets - NOT a single draw:
#   grey bar    = 2.5th-97.5th percentile of the 1000 random-set HRs
#   grey point  = median random-set HR (centre of the null)
#   red diamond = real Egr1 regulon HR
#   label       = real HR + one-sided empirical p
# If the red diamond sits outside the grey bar, the regulon beats random.
# =============================================================================

library(dplyr)
library(ggplot2)

if (!exists("null_long_df") || !exists("null_summary_df"))
  stop("Run Egr1CC_regulon_SCRAMBLED_control_pvalues.R first (same session).")
if (!exists("output_dir")) output_dir <- getwd()

use_log_axis <- FALSE   # TRUE = log10 x-axis (better if the female bars dominate)

# ---- summarise the 1000-permutation null per group -------------------------
perm_forest <- null_long_df %>%
  group_by(Pairing) %>%
  summarise(null_median = median(HR),
            null_lo     = quantile(HR, 0.025),
            null_hi     = quantile(HR, 0.975),
            real_hr     = first(Real_HR),
            .groups = "drop") %>%
  left_join(
    null_summary_df %>%
      transmute(Pairing   = paste0(Geneset, "\n(", Dataset, ")"),
                emp_p_one = Empirical_P_onesided),
    by = "Pairing")

# row order: MaleCC_WT top -> FemaleCC_KD bottom
row_order <- c("MaleCC_WT_enriched\n(Males)",   "MaleCC_KD_enriched\n(Males)",
               "FemaleCC_WT_enriched\n(Females)", "FemaleCC_KD_enriched\n(Females)")
perm_forest$Pairing <- factor(perm_forest$Pairing, levels = rev(row_order))

perm_forest$plabel <- sprintf("real HR = %.2f\np = %.3g",
                              perm_forest$real_hr, perm_forest$emp_p_one)

print(perm_forest)

x_text <- max(perm_forest$null_hi, perm_forest$real_hr) * 1.05

# ---- forest plot -----------------------------------------------------------
p_perm_forest <- ggplot(perm_forest, aes(y = Pairing)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = null_lo, xmax = null_hi),
                 height = 0.18, color = "grey55", linewidth = 1.1) +
  geom_point(aes(x = null_median, color = "Random gene sets (n = 1000)"), size = 3.5) +
  geom_point(aes(x = real_hr,     color = "Real Egr1 regulon"), shape = 18, size = 6) +
  geom_text(aes(x = x_text, label = plabel), hjust = 0, size = 3.2,
            color = "black", lineheight = 0.9) +
  scale_color_manual(name = NULL,
                     values = c("Random gene sets (n = 1000)" = "grey40",
                                "Real Egr1 regulon"           = "#E41A1C")) +
  labs(x = "Hazard Ratio   (bar = 2.5-97.5th percentile of 1000 random sets)",
       y = NULL,
       title = "Egr1 regulon vs permutation null - TCGA GBM",
       subtitle = "Grey = size-matched random gene sets | Red = real regulon") +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        axis.text.y   = element_text(size = 10),
        legend.position = "top",
        panel.grid.minor = element_blank())

if (use_log_axis) {
  p_perm_forest <- p_perm_forest +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.40)))
} else {
  p_perm_forest <- p_perm_forest +
    scale_x_continuous(limits = c(0, x_text * 1.45))
}

print(p_perm_forest)
ggsave(file.path(output_dir, "forest_plot_permutation_null.pdf"),
       plot = p_perm_forest, width = 9, height = 5)
cat("Saved forest_plot_permutation_null.pdf\n")
