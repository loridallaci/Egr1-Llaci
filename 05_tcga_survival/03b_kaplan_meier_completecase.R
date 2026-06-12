# =============================================================================
# 05_tcga_survival: 03b - Kaplan-Meier Plots on the COMPLETE-CASE (adjusted) sample
# =============================================================================
# Description:
#   Standard KM plots (03_kaplan_meier_plots.R) are drawn on every patient with
#   non-missing survival + expression, so the curve uses MORE patients than the
#   multivariable Cox model, which drops any patient missing a covariate
#   (Recurrence, Age, Subtype, MGMT_status). That makes the curve (e.g. n = 368)
#   and the annotated Cox HR (e.g. n = 271) describe different samples.
#
#   This script restricts each KM plot to the SAME complete-case rows the Cox
#   model uses, recomputes the median High/Low split within that subset, and
#   annotates each plot with the adjusted Cox HR / 95% CI / p-value AND the n.
#   The curve and the HR therefore describe one identical population.
#
# Run order:
#   1. 01_load_and_prepare_tcga_data_updated.R
#   2. 02_multivariate_cox_regression_updated_fixedForest.R
#   3. 03b_kaplan_meier_completecase.R        <-- this script
#
# Input (from environment, created by script 01):
#   - tcga_pheno, tcga_pheno_male, tcga_pheno_female
#   - geneID_included  (or geneID_included_M / geneID_included_F)
#
# Output:
#   - TCGA_KM_completecase_[group]_[gene].pdf  (one per gene per group)
# =============================================================================

# --- Libraries ----------------------------------------------------------------

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

if (!exists("tcga_pheno_male")) {
  stop("Run 01_load_and_prepare_tcga_data_updated.R first.")
}

# --- Settings -----------------------------------------------------------------

# Covariates MUST match those used in the Cox model (utils.R::run_cox_for_group)
covariates <- c("Recurrence", "Age", "Subtype", "MGMT_status")
OS_var     <- "survival"
OS_event   <- "status"

# --- Paths (edit for your system) ---------------------------------------------
# GitHub default: relative path, so the repo runs anywhere.
output_dir <- "output/km_completecase"
# Lori's laptop run:
# output_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/km_completecase"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Genes to plot: prefer a single combined list, else the union of the
# sex-specific motif gene lists.
if (exists("geneID_included")) {
  genes_to_plot <- geneID_included
} else {
  genes_to_plot <- unique(c(
    if (exists("geneID_included_M")) geneID_included_M else character(0),
    if (exists("geneID_included_F")) geneID_included_F else character(0)
  ))
}

# Make sure covariates are factors (script 01 already does this, but be safe)
factorize_cols <- function(df) {
  for (v in c("Recurrence", "Subtype", "MGMT_status")) {
    if (v %in% colnames(df)) df[[v]] <- factor(df[[v]])
  }
  df
}

# --- Helper: complete-case KM plot with adjusted Cox annotation ---------------

make_km_completecase <- function(gene_name, pheno_data, group_label, output_dir) {

  if (!gene_name %in% colnames(pheno_data)) {
    warning(paste("Gene", gene_name, "not in", group_label, "- skipping"))
    return(invisible(NULL))
  }

  pheno_data <- factorize_cols(pheno_data)

  # --- Restrict to the SAME rows coxph() would use (complete cases) ----------
  model_vars <- c(OS_var, OS_event, gene_name, covariates)
  cc         <- complete.cases(pheno_data[, model_vars])
  dat        <- pheno_data[cc, , drop = FALSE]

  n_cc <- nrow(dat)
  if (n_cc < 10) {
    warning(paste("Too few complete cases for", gene_name, "in", group_label,
                  "(n =", n_cc, ") - skipping"))
    return(invisible(NULL))
  }

  # --- Recompute median High/Low split WITHIN the complete-case subset -------
  # (so the split reflects the analyzed population; switch to the precomputed
  #  <gene>_binary column here if you prefer the whole-cohort median instead)
  med <- median(dat[[gene_name]], na.rm = TRUE)
  dat$.expr_group <- factor(ifelse(dat[[gene_name]] >= med, "High", "Low"),
                            levels = c("Low", "High"))

  if (length(unique(dat$.expr_group)) < 2) {
    warning(paste("Only one expression group for", gene_name, "in", group_label,
                  "- skipping"))
    return(invisible(NULL))
  }

  # --- Adjusted Cox model on the SAME subset (for HR / CI / p) ----------------
  cox_formula <- as.formula(
    paste("Surv(", OS_var, ",", OS_event, ") ~", gene_name, "+",
          paste(covariates, collapse = "+"))
  )

  hr_lab <- NULL
  p_lab  <- NULL
  tryCatch({
    cox  <- coxph(cox_formula, data = dat)
    csum <- summary(cox)
    ci   <- csum$conf.int
    cf   <- csum$coefficients
    if (gene_name %in% rownames(ci)) {
      hr   <- ci[gene_name, "exp(coef)"]
      lo   <- ci[gene_name, "lower .95"]
      hi   <- ci[gene_name, "upper .95"]
      pval <- cf[gene_name, "Pr(>|z|)"]
      hr_lab <- paste0("Adjusted HR = ", round(hr, 2),
                       " (", round(lo, 2), "–", round(hi, 2), ")")
      p_lab  <- paste0("Cox p = ", signif(pval, 3))
    }
  }, error = function(e) {
    warning(paste("Cox failed for", gene_name, "in", group_label, ":", e$message))
  })

  # --- KM fit on the complete-case subset ------------------------------------
  fit <- survfit(Surv(survival, status) ~ .expr_group, data = dat)

  x_pos <- max(dat[[OS_var]], na.rm = TRUE) * 0.45

  p1 <- ggsurvplot(
    fit,
    data             = dat,
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
    title            = paste0(gene_name, " - ", group_label,
                              " GBM (complete-case, n = ", n_cc, ")")
  )

  # Annotate with adjusted Cox results + group sizes
  y0 <- 0.92
  if (!is.null(hr_lab)) {
    p1$plot <- p1$plot +
      annotate("text", x = x_pos, y = y0,        hjust = 0, size = 3.5, label = hr_lab) +
      annotate("text", x = x_pos, y = y0 - 0.07, hjust = 0, size = 3.5, label = p_lab)
  }
  p1$plot <- p1$plot +
    annotate("text", x = x_pos, y = y0 - 0.14, hjust = 0, size = 3.2,
             label = paste0("Low: ", sum(dat$.expr_group == "Low"),
                            "  High: ", sum(dat$.expr_group == "High")))

  pdf_file <- file.path(output_dir,
                        paste0("TCGA_KM_completecase_", group_label, "_", gene_name, ".pdf"))
  ggsave(filename = pdf_file, plot = p1$plot, device = "pdf", width = 8, height = 6)

  invisible(p1$plot)
}

# --- Generate complete-case KM plots for all groups ---------------------------

datasets <- list(
  Male   = tcga_pheno_male,
  Female = tcga_pheno_female,
  All    = tcga_pheno
)

for (group_label in names(datasets)) {
  message("Generating complete-case KM plots for ", group_label, "...")
  for (gene in genes_to_plot) {
    make_km_completecase(gene, datasets[[group_label]], group_label, output_dir)
  }
}

message("Done. Complete-case KM plots saved to: ", output_dir)
