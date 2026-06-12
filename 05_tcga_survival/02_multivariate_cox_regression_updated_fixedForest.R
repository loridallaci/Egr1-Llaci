# =============================================================================
# 05_tcga_survival: 02 - Multivariate Cox Regression
# (Forest plots use AUTO-SCALING x-axis from 4-gene script)
# =============================================================================
# Logic:
#   - Male motif genes   -> Cox on Males, Females, All
#   - Female motif genes -> Cox on Males, Females, All
#   = 6 combinations total
#
# Outputs per combination (motif_sex x dataset):
#   - Long CSV (one row per gene per variable)
#   - Wide CSV (one row per gene, all covariates spread)
#   - Forest plot PDFs   <-- now auto-scaled per gene
#   - Kaplan-Meier plot PDFs
# =============================================================================

library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survminer)

source("05_tcga_survival/utils.R")

if (!exists("tcga_pheno_male")) {
  stop("Run 01_load_and_prepare_tcga_data.R first.")
}

# =============================================================================
# OVERRIDE: AUTO-SCALING FOREST PLOT (taken from 4-gene script)
# This replaces whatever create_forest_plot() is defined in utils.R
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
  
  # AUTO SCALE AXIS BASED ON DATA (no fixed xlim)
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
    coord_cartesian(xlim = xlim_range) +
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

# --- Paths --------------------------------------------------------------------

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/tcga_all_motifs/"
output_dir <- "output/tcga_all_motifs/"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Covariates ---------------------------------------------------------------

all_covarID <- c("Recurrence", "Age", "Subtype", "MGMT_status")
OS_var      <- "survival"
OS_event    <- "status"

factorize_cols <- function(df) {
  df$Recurrence  <- factor(df$Recurrence)
  df$Subtype     <- factor(df$Subtype)
  df$MGMT_status <- factor(df$MGMT_status)
  df
}

tcga_pheno        <- factorize_cols(tcga_pheno)
tcga_pheno_male   <- factorize_cols(tcga_pheno_male)
tcga_pheno_female <- factorize_cols(tcga_pheno_female)

# --- Gene lists per motif sex -------------------------------------------------

if (!exists("geneID_included_M") || !exists("geneID_included_F")) {
  geneID_included_M <- df_split %>%
    filter(sex == "Male") %>%
    pull(SYMBOL_UPPER) %>%
    unique() %>%
    .[. %in% colnames(tcga_pheno)]
  
  geneID_included_F <- df_split %>%
    filter(sex == "Female") %>%
    pull(SYMBOL_UPPER) %>%
    unique() %>%
    .[. %in% colnames(tcga_pheno)]
}

cat("Male motif genes found in TCGA:  ", length(geneID_included_M), "\n")
cat("Female motif genes found in TCGA:", length(geneID_included_F), "\n")

# --- Dataset list -------------------------------------------------------------

datasets <- list(
  Males   = tcga_pheno_male,
  Females = tcga_pheno_female,
  All     = tcga_pheno
)

# --- Motif gene lists ---------------------------------------------------------

motif_genesets <- list(
  MaleMotifs   = geneID_included_M,
  FemaleMotifs = geneID_included_F
)

# --- Helper: wide output ------------------------------------------------------

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
      
      row_out <- data.frame(SYMBOL_UPPER   = gene_name,
                            Overall_pvalue = as.numeric(overall_p),
                            stringsAsFactors = FALSE)
      
      for (var in rownames(df_coef)) {
        
        var_label <- var
        var_label <- gsub(gene_name,                 "Expression",           var_label)
        var_label <- gsub("MGMT_statusUnmethylated", "MGMT_Unmethylated",    var_label)
        var_label <- gsub("RecurrenceRecurrent",     "Recurrence_Recurrent", var_label)
        var_label <- gsub("RecurrenceSecondary",     "Recurrence_Secondary", var_label)
        var_label <- gsub("SubtypeMesenchymal",      "Subtype_Mesenchymal",  var_label)
        var_label <- gsub("SubtypeProneural",        "Subtype_Proneural",    var_label)
        
        row_out[[paste0(var_label, "_coef_",    suffix)]] <- df_ci[var,   "exp(coef)"]
        row_out[[paste0(var_label, "_Lower95_", suffix)]] <- df_ci[var,   "lower .95"]
        row_out[[paste0(var_label, "_Upper95_", suffix)]] <- df_ci[var,   "upper .95"]
        row_out[[paste0(var_label, "_SE_",      suffix)]] <- df_coef[var, "se(coef)"]
        row_out[[paste0(var_label, "_pvalue_",  suffix)]] <- df_coef[var, "Pr(>|z|)"]
      }
      
      wide_list[[gene_name]] <- row_out
      
    }, error = function(e) {
      warning(paste("Error with gene", gene_name, ":", e$message))
    })
  }
  
  bind_rows(wide_list)
}

# =============================================================================
# MAIN LOOP: motif geneset x dataset
# =============================================================================

all_long_results <- list()

for (motif_sex in names(motif_genesets)) {
  
  genes <- motif_genesets[[motif_sex]]
  cat("\n===", motif_sex, "- genes:", length(genes), "===\n")
  
  out_sub  <- file.path(output_dir, motif_sex)
  plot_sub <- file.path(out_sub, "forest_plots")
  km_sub   <- file.path(out_sub, "km_plots")
  
  dir.create(out_sub,  showWarnings = FALSE, recursive = TRUE)
  dir.create(plot_sub, showWarnings = FALSE, recursive = TRUE)
  dir.create(km_sub,   showWarnings = FALSE, recursive = TRUE)
  
  wide_results <- list()
  
  for (dataset_name in names(datasets)) {
    
    pheno <- datasets[[dataset_name]]
    cat("  Running Cox:", motif_sex, "x", dataset_name, "\n")
    
    suffix <- paste0(substr(motif_sex, 1, 1), substr(dataset_name, 1, 1))
    
    # --- Long Cox -----------------------------------------------------------
    long_res <- run_cox_for_group(
      pheno_data = pheno,
      genes      = genes,
      OS_var     = OS_var,
      OS_event   = OS_event,
      covariates = all_covarID
    )
    
    write.csv(long_res,
              file.path(out_sub, paste0("cox_long_", motif_sex, "_", dataset_name, "_updated.csv")),
              row.names = FALSE)
    
    all_long_results[[paste0(motif_sex, "_", dataset_name)]] <- long_res
    
    # --- Wide ---------------------------------------------------------------
    wide_res <- make_wide_output(
      pheno_data = pheno,
      genes      = genes,
      OS_var     = OS_var,
      OS_event   = OS_event,
      covariates = all_covarID,
      suffix     = suffix
    )
    
    write.csv(wide_res,
              file.path(out_sub, paste0("cox_wide_", motif_sex, "_", dataset_name, "_updated.csv")),
              row.names = FALSE)
    
    wide_results[[dataset_name]] <- wide_res %>%
      rename(!!paste0("Overall_pvalue_", dataset_name) := Overall_pvalue)
    
    # --- Forest plots (AUTO-SCALED x-axis) ----------------------------------
    for (gene in genes) {
      
      if (!gene %in% colnames(pheno)) next
      
      tryCatch({
        p <- create_forest_plot(
          gene_name     = gene,
          survival_data = long_res,
          title_suffix  = paste0(" | ", motif_sex, " | ", dataset_name)
          # NOTE: xlim_range omitted -> auto-scales from data
        )
        
        if (!is.null(p)) {
          ggsave(
            filename = file.path(plot_sub,
                                 paste0("forest_", gene, "_", dataset_name, ".pdf")),
            plot = p, width = 6, height = 4
          )
        }
      }, error = function(e) {
        warning(paste("Forest plot failed for", gene, ":", e$message))
      })
    }
    
    # --- Kaplan-Meier plots -------------------------------------------------
    for (gene in genes) {
      
      binary_col <- paste0(gene, "_binary")
      if (!binary_col %in% colnames(pheno)) next
      
      obs_row <- long_res %>% filter(Gene == gene, Variable == gene)
      
      if (nrow(obs_row) > 0) {
        hr_label <- paste0("HR (High vs Low) = ", round(obs_row$HR[1], 2),
                           " (", round(obs_row$Lower95[1], 2),
                           "\u2013", round(obs_row$Upper95[1], 2), ")")
        p_label  <- paste0("Multivariate Cox p = ", signif(obs_row$Pvalue[1], 3))
      } else {
        hr_label <- NULL
        p_label  <- NULL
      }
      
      tryCatch({
        fit <- survfit(as.formula(
          paste0("Surv(", OS_var, ", ", OS_event, ") ~ ", binary_col)
        ), data = pheno)
        
        p1 <- ggsurvplot(
          fit,
          data             = pheno,
          risk.table       = TRUE,
          pval             = FALSE,
          conf.int         = FALSE,
          xlab             = "OS (months)",
          legend.title     = gene,
          legend.labs      = c("Low", "High"),
          legend           = "top",
          surv.median.line = "hv",
          palette          = "npg",
          ggtheme          = theme_bw(),
          title            = paste0(gene, " - ", motif_sex, " | ", dataset_name)
        )
        
        if (!is.null(hr_label)) {
          p1$plot <- p1$plot +
            ggplot2::annotate("text", x = 30, y = 0.90, hjust = 0,
                              size = 3.5, label = hr_label) +
            ggplot2::annotate("text", x = 30, y = 0.80, hjust = 0,
                              size = 3.5, label = p_label)
        }
        
        ggsave(
          filename = file.path(km_sub,
                               paste0("KM_", gene, "_", dataset_name, ".pdf")),
          plot  = p1$plot,
          width = 7, height = 5
        )
      }, error = function(e) {
        warning(paste("KM plot failed for", gene, "in", dataset_name, ":", e$message))
      })
    }
    
  } # end dataset loop
  
  # --- Merge wide tables across all 3 datasets --------------------------------
  wide_combined <- wide_results[["Males"]] %>%
    full_join(wide_results[["Females"]], by = "SYMBOL_UPPER") %>%
    full_join(wide_results[["All"]],     by = "SYMBOL_UPPER")
  
  write.csv(wide_combined,
            file.path(out_sub, paste0("cox_wide_", motif_sex, "_COMBINED_updated.csv")),
            row.names = FALSE)
  
  cat("  Saved combined wide table for", motif_sex, "\n")
  
} # end motif sex loop

saveRDS(all_long_results, file.path(output_dir, "all_long_cox_results_updated.rds"))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== DONE ===\n")
cat("Male motif genes:   ", length(geneID_included_M), "\n")
cat("Female motif genes: ", length(geneID_included_F), "\n")
cat("Male rows:          ", nrow(tcga_pheno_male), "\n")
cat("Female rows:        ", nrow(tcga_pheno_female), "\n")
cat("All rows:           ", nrow(tcga_pheno), "\n")
