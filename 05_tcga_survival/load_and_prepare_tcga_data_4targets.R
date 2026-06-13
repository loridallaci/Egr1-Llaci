# =============================================================================
# 05_tcga_survival: 01 - Load and Prepare TCGA Data
# =============================================================================
# Description:
#   Loads TCGA GBM phenotype and expression data, filters to IDH wild-type
#   samples, matches samples between datasets, adds gene expression to
#   phenotype data frame, and splits into sex-stratified subsets.
#
# Run order (within 05_tcga_survival):
#   1. 01_load_and_prepare_tcga_data.R   <-- this script
#   2. 02_multivariate_cox_regression.R
#   3. 03_kaplan_meier_plots.R
#   4. 04_figures_chromvar_vs_survival.R
#
# Input:
#   - 2024-06-04_TCGA_GBM_pheno.txt
#   - 2024-06-04_TCGA_GBM_expression.txt
#   - RENINmotif_Male_minus_Female_ActivityScore_all_motifs.csv
#     (or Macs2_Peaks version)
#
# Output (saved to environment for downstream scripts):
#   - tcga_pheno        : all IDH wild-type samples with gene expression
#   - tcga_pheno_male   : male subset
#   - tcga_pheno_female : female subset
#   - tcga_exp          : expression matrix
#   - geneID_included   : genes found in TCGA expression data
#   - df_split          : RENIN motif data with TF names split by "::"
# =============================================================================

# --- Libraries ----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(kableExtra)
library(dplyr)
library(rlang)
library(knitr)

#source("05_tcga_survival/utils.R")

# --- Paths --------------------------------------------------------------------

# original (author's machine): "C:/Users/loril/Documents/multivariate_analysis/glioVis"
dat_dir ="data"
# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"
output_dir <- "output/"
# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/TCGA_Custom4Genes/"
motif_directory <- "output/TCGA_Custom4Genes/"

# --- Paths --------------------------------------------------------------------
library(dplyr)

# original (author's machine): "C:/Users/loril/Documents/multivariate_analysis/glioVis"
dat_dir <- "data"
#output_dir <- "output/"   # original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"

# ---------------- LOAD PHENO ----------------

tcga_pheno <- read.table(
  file.path(dat_dir, "2024-06-04_TCGA_GBM_pheno.txt"),
  sep = "\t", header = TRUE, check.names = FALSE
)

tcga_pheno$Sample <- as.character(tcga_pheno$Sample)
rownames(tcga_pheno) <- tcga_pheno$Sample

# ---------------- LOAD EXPRESSION ----------------

tcga_exp <- read.table(
  file.path(dat_dir, "2024-06-04_TCGA_GBM_expression.txt"),
  sep = "\t", header = TRUE, check.names = FALSE
) %>% as.data.frame()

rownames(tcga_exp) <- tcga_exp$Sample
tcga_exp <- tcga_exp %>% dplyr::select(-Sample) %>% as.matrix()

# ---------------- ALIGN ----------------

common_samples <- intersect(tcga_pheno$Sample, rownames(tcga_exp))

tcga_pheno <- tcga_pheno[common_samples, ]
tcga_exp   <- tcga_exp[common_samples, ]

stopifnot(all(rownames(tcga_exp) == tcga_pheno$Sample))

# ---------------- FILTER ----------------

tcga_pheno <- tcga_pheno %>%
  filter(IDH1_status == "Wild-type" & !is.na(Gender))

tcga_exp <- tcga_exp[tcga_pheno$Sample, ]

# ---------------- MERGE ----------------

tcga_pheno <- cbind(tcga_pheno, tcga_exp)

tcga_pheno_male   <- tcga_pheno %>% filter(Gender == "Male")
tcga_pheno_female <- tcga_pheno %>% filter(Gender == "Female")

# ---------------- FACTOR COVARIATES ----------------

factorize_cols <- function(df) {
  df$Recurrence  <- factor(df$Recurrence)
  df$Subtype     <- factor(df$Subtype)
  df$MGMT_status <- factor(df$MGMT_status)
  df
}

tcga_pheno        <- factorize_cols(tcga_pheno)
tcga_pheno_male   <- factorize_cols(tcga_pheno_male)
tcga_pheno_female <- factorize_cols(tcga_pheno_female)

# ---------------- FIXED 4-GENE SET ----------------

geneID_included <- c("PTK2B", "NRP1", "HMOX1", "MAST4")
geneID_included <- geneID_included[geneID_included %in% colnames(tcga_pheno)]

# ---------------- BINARY ONLY FOR 4 GENES ----------------

add_binary_cols <- function(df, genes) {
  for (gene in genes) {
    df[[paste0(gene, "_binary")]] <- factor(
      ifelse(df[[gene]] >= median(df[[gene]], na.rm = TRUE),
             "High", "Low"),
      levels = c("Low", "High")
    )
  }
  df
}

tcga_pheno        <- add_binary_cols(tcga_pheno, geneID_included)
tcga_pheno_male   <- add_binary_cols(tcga_pheno_male, geneID_included)
tcga_pheno_female <- add_binary_cols(tcga_pheno_female, geneID_included)

# ---------------- OUTPUT OBJECTS ----------------

cat("FINAL 4-GENE SET:\n")
print(geneID_included)

cat("Samples total:", nrow(tcga_pheno), "\n")

library(survival)
library(dplyr)
library(ggplot2)
library(survminer)

#source("05_tcga_survival/utils.R")

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"
output_dir <- "output/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------- INPUT CHECK ----------------

stopifnot(exists("tcga_pheno"))

# ---------------- SETTINGS ----------------

geneID_included <- c("PTK2B", "NRP1", "HMOX1", "MAST4")

all_covarID <- c("Recurrence", "Age", "Subtype", "MGMT_status")
OS_var   <- "survival"
OS_event <- "status"

datasets <- list(
  Males   = tcga_pheno_male,
  Females = tcga_pheno_female,
  All     = tcga_pheno
)

# ---------------- COX FUNCTION ----------------

run_cox <- function(pheno_data, genes) {
  
  res <- list()
  
  for (gene in genes) {
    
    if (!gene %in% colnames(pheno_data)) next
    
    f <- as.formula(
      paste0("Surv(", OS_var, ",", OS_event, ") ~ ",
             gene, " + ", paste(all_covarID, collapse = " + "))
    )
    
    fit <- tryCatch(coxph(f, data = pheno_data), error = function(e) NULL)
    if (is.null(fit)) next
    
    s <- summary(fit)
    
    res[[gene]] <- data.frame(
      Gene = gene,
      HR = s$coefficients[gene, "exp(coef)"],
      L95 = s$conf.int[gene, "lower .95"],
      U95 = s$conf.int[gene, "upper .95"],
      P   = s$coefficients[gene, "Pr(>|z|)"]
    )
  }
  
  bind_rows(res)
}

# ---------------- RUN ALL ----------------

all_results <- list()

for (d in names(datasets)) {
  
  pheno <- datasets[[d]]
  
  out <- run_cox(pheno, geneID_included)
  
  write.csv(out,
            file.path(output_dir, paste0("cox_4gene_", d, ".csv")),
            row.names = FALSE)
  
  all_results[[d]] <- out
}

saveRDS(all_results,
        file.path(output_dir, "cox_4gene_all_results.rds"))

cat("DONE — 4-gene Cox pipeline complete\n")

******************************************************************************
  # =============================================================================
# 05_tcga_survival: 02 - Multivariate Cox Regression (4-GENE FIXED VERSION)
# =============================================================================
create_forest_plot <- function(gene_name, survival_data,
                               title_suffix = "",
                               var_order    = NULL,
                               xlim_range   = NULL) {
  
  gene_data <- survival_data %>% dplyr::filter(Gene == gene_name)
  
  if (nrow(gene_data) == 0) return(NULL)
  
  unique_vars <- unique(gene_data$Variable)
  
  if (is.null(var_order)) {
    var_order <- c(gene_name,
                   "Age",
                   "MGMT (unmethylated)",
                   "Recurrence (recurrent)",
                   "Recurrence (secondary)",
                   "Mesenchymal Subtype",
                   "Proneural Subtype")
  }
  
  var_order <- var_order[var_order %in% unique_vars]
  
  plot_df <- gene_data %>%
    mutate(
      signif   = ifelse(Pvalue < 0.05, "p<0.05", "ns"),
      Variable = factor(Variable, levels = rev(var_order))
    )
  
  # 🔥 AUTO SCALE AXIS BASED ON DATA
  if (is.null(xlim_range)) {
    xmin <- min(plot_df$Lower95, na.rm = TRUE)
    xmax <- max(plot_df$Upper95, na.rm = TRUE)
    
    xlim_range <- c(
      max(0, xmin * 0.9),
      xmax * 1.1
    )
  }
  
  ggplot(plot_df, aes(x = HR, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_point(aes(color = signif), size = 4) +
    geom_errorbarh(aes(xmin = Lower95, xmax = Upper95), height = 0.3) +
    scale_color_manual(values = c("p<0.05" = "red", "ns" = "black")) +
    coord_cartesian(xlim = xlim_range) +   # ✅ FIXED
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title.y    = element_blank(),
      axis.text.y     = element_text(size = 12),
      plot.title      = element_text(hjust = 0.5, face = "bold")
    ) +
    xlab("Hazard Ratio") +
    ggtitle(paste0(gene_name, title_suffix))
}


library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survminer)

#source("05_tcga_survival/utils.R")

#if (!exists("tcga_pheno_male")) {
#  stop("Run 01_load_and_prepare_tcga_data.R first.")
#}

# --- Paths --------------------------------------------------------------------

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/TCGA_Custom4Genes/"
output_dir <- "output/TCGA_Custom4Genes/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Covariates --------------------------------------------------------------

all_covarID <- c("Recurrence", "Age", "Subtype", "MGMT_status")
OS_var   <- "survival"
OS_event <- "status"

factorize_cols <- function(df) {
  df$Recurrence  <- factor(df$Recurrence)
  df$Subtype     <- factor(df$Subtype)
  df$MGMT_status <- factor(df$MGMT_status)
  df
}

tcga_pheno        <- factorize_cols(tcga_pheno)
tcga_pheno_male   <- factorize_cols(tcga_pheno_male)
tcga_pheno_female <- factorize_cols(tcga_pheno_female)

# =============================================================================
# FIXED GENE SET (ONLY YOUR 4 GENES)
# =============================================================================

geneID_included <- c("PTK2B", "NRP1", "HMOX1", "MAST4", "EGR1")

cat("Using ONLY 4 genes:\n")
print(geneID_included)

# keep only genes that exist in dataset
geneID_included <- geneID_included[geneID_included %in% colnames(tcga_pheno)]

cat("Genes found in TCGA:\n")
print(geneID_included)

# =============================================================================
# CREATE BINARY EXPRESSION ONLY FOR 4 GENES
# =============================================================================

add_binary_cols <- function(df, genes) {
  for (gene in genes) {
    if (!gene %in% colnames(df)) next
    
    df[[paste0(gene, "_binary")]] <- factor(
      ifelse(df[[gene]] >= median(df[[gene]], na.rm = TRUE),
             "High", "Low"),
      levels = c("Low", "High")
    )
  }
  df
}

tcga_pheno        <- add_binary_cols(tcga_pheno, geneID_included)
tcga_pheno_male   <- add_binary_cols(tcga_pheno_male, geneID_included)
tcga_pheno_female <- add_binary_cols(tcga_pheno_female, geneID_included)

cat("Binary columns created for 4 genes only\n")

# =============================================================================
# DATASETS
# =============================================================================

datasets <- list(
  Males   = tcga_pheno_male,
  Females = tcga_pheno_female,
  All     = tcga_pheno
)

# =============================================================================
# WIDE OUTPUT FUNCTION (UNCHANGED BUT SAFE)
# =============================================================================

make_wide_output <- function(pheno_data, genes, OS_var, OS_event,
                             covariates, suffix) {
  
  wide_list <- list()
  
  for (gene_name in genes) {
    
    if (!gene_name %in% colnames(pheno_data)) next
    
    formula <- as.formula(
      paste("Surv(", OS_var, ",", OS_event, ") ~", gene_name, "+",
            paste(covariates, collapse = "+"))
    )
    
    tryCatch({
      cox_obj <- coxph(formula, data = pheno_data)
      cox_sum <- summary(cox_obj)
      
      overall_p <- cox_sum$logtest["pvalue"]
      df_coef   <- as.data.frame(cox_sum$coefficients)
      df_ci     <- as.data.frame(cox_sum$conf.int)
      
      row_out <- data.frame(
        SYMBOL_UPPER = gene_name,
        Overall_pvalue = as.numeric(overall_p)
      )
      
      for (var in rownames(df_coef)) {
        
        var_label <- var
        var_label <- gsub(gene_name, "Expression", var_label)
        
        row_out[[paste0(var_label, "_coef_", suffix)]] <- df_ci[var, "exp(coef)"]
        row_out[[paste0(var_label, "_Lower95_", suffix)]] <- df_ci[var, "lower .95"]
        row_out[[paste0(var_label, "_Upper95_", suffix)]] <- df_ci[var, "upper .95"]
        row_out[[paste0(var_label, "_SE_", suffix)]] <- df_coef[var, "se(coef)"]
        row_out[[paste0(var_label, "_pvalue_", suffix)]] <- df_coef[var, "Pr(>|z|)"]
      }
      
      wide_list[[gene_name]] <- row_out
      
    }, error = function(e) {
      warning(paste("Error with gene", gene_name, ":", e$message))
    })
  }
  
  bind_rows(wide_list)
}

# =============================================================================
# MAIN LOOP
# =============================================================================

all_long_results <- list()

for (dataset_name in names(datasets)) {
  
  pheno <- datasets[[dataset_name]]
  
  cat("\n=== Running dataset:", dataset_name, "===\n")
  
  out_sub  <- file.path(output_dir, "4gene_analysis", dataset_name)
  plot_sub <- file.path(out_sub, "forest_plots")
  km_sub   <- file.path(out_sub, "km_plots")
  
  dir.create(plot_sub, recursive = TRUE, showWarnings = FALSE)
  dir.create(km_sub, recursive = TRUE, showWarnings = FALSE)
  
  # --- LONG COX --------------------------------------------------------------
  
  long_res <- run_cox_for_group(
    pheno_data = pheno,
    genes      = geneID_included,
    OS_var     = OS_var,
    OS_event   = OS_event,
    covariates = all_covarID
  )
  
  write.csv(long_res,
            file.path(out_sub, "cox_long_4genes.csv"),
            row.names = FALSE)
  
  all_long_results[[dataset_name]] <- long_res
  
  # --- WIDE ------------------------------------------------------------------
  
  wide_res <- make_wide_output(
    pheno, geneID_included,
    OS_var, OS_event,
    all_covarID,
    suffix = dataset_name
  )
  
  write.csv(wide_res,
            file.path(out_sub, "cox_wide_4genes.csv"),
            row.names = FALSE)
  
  # --- FOREST PLOTS (FIXED & GUARANTEED RUN) -------------------------------
  
  for (gene in geneID_included) {
    
    if (!gene %in% colnames(pheno)) next
    
    tryCatch({
      
      p <- create_forest_plot(
        gene_name = gene,
        survival_data = long_res,
        title_suffix = paste0(" | ", dataset_name)
      )
      
      ggsave(
        filename = file.path(plot_sub,
                             paste0("forest_", gene, "_", dataset_name, ".pdf")),
        plot = p,
        width = 6,
        height = 4
      )
      
    }, error = function(e) {
      warning(paste("Forest plot failed:", gene, e$message))
    })
  }
  
  # --- KM PLOTS -------------------------------------------------------------
  
  for (gene in geneID_included) {
    
    binary_col <- paste0(gene, "_binary")
    if (!binary_col %in% colnames(pheno)) next
    
    tryCatch({
      
      fit <- survfit(as.formula(
        paste0("Surv(", OS_var, ",", OS_event, ") ~ ", binary_col)
      ), data = pheno)
      
      p1 <- ggsurvplot(
        fit,
        data = pheno,
        risk.table = TRUE,
        pval = FALSE,
        conf.int = FALSE,
        xlab = "OS (months)",
        legend.title = gene,
        legend.labs = c("Low", "High"),
        ggtheme = theme_bw()
      )
      
      ggsave(
        file.path(km_sub,
                  paste0("KM_", gene, "_", dataset_name, ".pdf")),
        p1$plot,
        width = 7,
        height = 5
      )
      
    }, error = function(e) {
      warning(paste("KM failed:", gene, e$message))
    })
  }
}

# =============================================================================
# SAVE
# =============================================================================

saveRDS(all_long_results,
        file.path(output_dir, "4gene_all_long_results.rds"))

cat("\nDONE: 4-gene Cox + forest + KM pipeline complete\n")


