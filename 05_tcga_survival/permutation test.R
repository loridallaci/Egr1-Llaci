# =============================================================================
# 05_tcga_survival: 02b - Permutation Test for Multivariate Cox Regression
# =============================================================================
# Description:
#   For each TF/gene, tests whether its prognostic effect in the multivariate
#   Cox model is INDEPENDENT of the covariates (Recurrence, Age, Subtype,
#   MGMT_status). This is done by permuting the gene expression column many
#   times and comparing the observed log-HR to the resulting null distribution.
#
#   A significant permutation p-value means the gene's prognostic signal is
#   NOT explained by the covariates alone — i.e., it adds independent value.
#
# Run order:
#   1. 01_load_and_prepare_tcga_data.R
#   2. 02_multivariate_cox_regression.R   (provides survival_table_male/female/all)
#   3. 02b_permutation_test.R             <-- this script
#
# Output:
#   - permutation_results_MALES.csv
#   - permutation_results_FEMALES.csv
#   - permutation_results_ALL.csv
#   - permutation_results_COMBINED.csv
# =============================================================================

library(survival)
library(dplyr)

source("05_tcga_survival/utils.R")

# Ensure upstream results are available
if (!exists("tcga_pheno_male") || !exists("survival_table_male")) {
  stop("Run 01_load_and_prepare_tcga_data.R and 02_multivariate_cox_regression.R first.")
}

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"
output_dir <- "output/"
# --- Settings -----------------------------------------------------------------

N_PERM    <- 1000   # Number of permutations (increase to 5000 for publication)
set.seed(42)        # Reproducibility

OS_var    <- "survival"
OS_event  <- "status"
all_covarID <- c("Recurrence", "Age", "Subtype", "MGMT_status")

# --- Run permutation tests ----------------------------------------------------

message(sprintf("Running permutation tests (%d permutations per gene)...", N_PERM))
message("This may take several minutes depending on the number of genes.\n")

perm_male   <- run_permutation_test(
  pheno_data  = tcga_pheno_male,
  genes       = geneID_included,
  obs_results = survival_table_male,
  OS_var      = OS_var,
  OS_event    = OS_event,
  covariates  = all_covarID,
  n_perm      = N_PERM
)

perm_female <- run_permutation_test(
  pheno_data  = tcga_pheno_female,
  genes       = geneID_included,
  obs_results = survival_table_female,
  OS_var      = OS_var,
  OS_event    = OS_event,
  covariates  = all_covarID,
  n_perm      = N_PERM
)

perm_all    <- run_permutation_test(
  pheno_data  = tcga_pheno,
  genes       = geneID_included,
  obs_results = survival_table_all,
  OS_var      = OS_var,
  OS_event    = OS_event,
  covariates  = all_covarID,
  n_perm      = N_PERM
)

# --- Save individual group results --------------------------------------------

write.csv(perm_male,   file.path(output_dir, "permutation_results_MALES.csv"),   row.names = FALSE)
write.csv(perm_female, file.path(output_dir, "permutation_results_FEMALES.csv"), row.names = FALSE)
write.csv(perm_all,    file.path(output_dir, "permutation_results_ALL.csv"),     row.names = FALSE)

# --- Combine and save ---------------------------------------------------------

perm_combined <- bind_rows(
  mutate(perm_male,   Group = "Males"),
  mutate(perm_female, Group = "Females"),
  mutate(perm_all,    Group = "All")
)

write.csv(perm_combined,
          file.path(output_dir, "permutation_results_COMBINED.csv"),
          row.names = FALSE)

# --- Merge permutation p-values into Cox summary ------------------------------
# Attaches perm_pvalue and perm_sig columns to the existing combined_summary

combined_summary_with_perm <- combined_summary %>%
  left_join(
    perm_combined %>%
      dplyr::select(Gene, Group, perm_pvalue, perm_padj, perm_sig),
    by = c("Gene", "Group")
  )

write.csv(combined_summary_with_perm,
          file.path(output_dir, "gene_HR_summary_with_permutation.csv"),
          row.names = FALSE)

# --- Summary ------------------------------------------------------------------

cat("\n=== Permutation Test Complete ===\n")
cat(sprintf("Permutations per gene: %d\n", N_PERM))
cat(sprintf("Genes tested: %d\n", length(geneID_included)))

for (grp in c("Males", "Females", "All")) {
  sub <- perm_combined %>% filter(Group == grp)
  cat(sprintf(
    "%s — perm p<0.05: %d | BH-adjusted p<0.05: %d\n",
    grp,
    sum(sub$perm_pvalue < 0.05, na.rm = TRUE),
    sum(sub$perm_padj   < 0.05, na.rm = TRUE)
  ))
}

cat("Results saved to:", output_dir, "\n")