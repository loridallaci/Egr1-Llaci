# =============================================================================
# 05_tcga_survival: 03 - Kaplan-Meier Plots
# =============================================================================
# Description:
#   Generates Kaplan-Meier survival plots for each TF in geneID_included,
#   split by median expression (High vs Low). Plots are made separately for
#   male, female, and all IDH wild-type GBM patients. Each plot is annotated
#   with multivariate Cox HR and p-value.
#
# Run order:
#   1. 01_load_and_prepare_tcga_data.R
#   2. 02_multivariate_cox_regression.R
#   3. 03_kaplan_meier_plots.R           <-- this script
#
# Input (from environment):
#   - tcga_pheno_male2, tcga_pheno_female2, tcga_pheno
#   - geneID_included
#   - survival_table_male, survival_table_female, survival_table_all
#     (from 02_multivariate_cox_regression.R)
#
# Output:
#   - TCGA_KM_[sex]_[gene].pdf  (one per gene per group)
# =============================================================================

# --- Libraries ----------------------------------------------------------------

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

source("05_tcga_survival/utils.R")

if (!exists("tcga_pheno_male2")) {
  stop("Run 01_load_and_prepare_tcga_data.R first.")
}
if (!exists("survival_table_male")) {
  stop("Run 02_multivariate_cox_regression.R first.")
}

# original (author's machine): "/home/lllaci/data/tcga_survival_results"
output_dir <- "output/tcga_survival_results"

# --- Helper: KM plot with Cox annotation -------------------------------------

make_km_with_cox <- function(gene_name, pheno_data, cox_table,
                             group_label, output_dir) {

  binary_col <- paste0(gene_name, "_binary")

  if (!binary_col %in% colnames(pheno_data)) {
    warning(paste("Binary column not found for", gene_name, "- skipping"))
    return(invisible(NULL))
  }

  # Extract Cox stats for this gene
  cox_row  <- cox_table %>%
    filter(Gene == gene_name, Variable == gene_name)

  if (nrow(cox_row) == 0) {
    warning(paste("No Cox results for", gene_name, "- skipping"))
    return(invisible(NULL))
  }

  hr_gene     <- round(cox_row$HR[1], 2)
  p_gene      <- signif(cox_row$Pvalue[1], 3)
  x_pos       <- max(pheno_data$survival, na.rm = TRUE) * 0.5

  # KM fit
  fit <- survfit(
    as.formula(paste0("Surv(survival, status) ~ ", binary_col)),
    data = pheno_data
  )

  p1 <- ggsurvplot(
    fit,
    data             = pheno_data,
    risk.table       = TRUE,
    pval             = FALSE,
    conf.int         = FALSE,
    xlab             = "OS (months)",
    legend.title     = gene_name,
    legend.labs      = c("Low", "High"),
    legend           = "top",
    surv.median.line = "hv",
    palette          = "npg",
    ggtheme          = theme_bw(base_size = 12),
    title            = paste0(gene_name, " - ", group_label, " GBM")
  )

  # Annotate with multivariate Cox results
  p1$plot <- p1$plot +
    annotate("text", x = x_pos, y = 0.90, hjust = 0,
             label = paste0("HR = ", hr_gene)) +
    annotate("text", x = x_pos, y = 0.83, hjust = 0,
             label = paste0("Cox p = ", p_gene))

  # Save
  pdf_file <- file.path(output_dir,
                        paste0("TCGA_KM_", group_label, "_", gene_name, ".pdf"))
  ggsave(filename = pdf_file, plot = p1$plot,
         device = "pdf", width = 8, height = 6)

  invisible(p1$plot)
}

# --- Generate KM plots for all genes ------------------------------------------

message("Generating KM plots for MALES...")
for (gene in geneID_included) {
  make_km_with_cox(gene, tcga_pheno_male2,   survival_table_male,   "Male",   output_dir)
}

message("Generating KM plots for FEMALES...")
for (gene in geneID_included) {
  make_km_with_cox(gene, tcga_pheno_female2, survival_table_female, "Female", output_dir)
}

message("Generating KM plots for ALL samples...")
for (gene in geneID_included) {
  make_km_with_cox(gene, tcga_pheno,         survival_table_all,    "All",    output_dir)
}

message("Done. KM plots saved to: ", output_dir)
