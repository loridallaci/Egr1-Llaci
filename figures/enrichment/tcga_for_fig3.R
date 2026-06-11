# =============================================================================
# 05_tcga_survival: CUSTOM GENES — MULTIVARIATE COX + FOREST + KM
# =============================================================================

library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survminer)

source("05_tcga_survival/utils.R")

# -----------------------------------------------------------------------------
# USER INPUT (YOUR GENES HERE)
# -----------------------------------------------------------------------------

genes_of_interest <- c("EGR1", "PTK2B", "NRP1", "HMOX1", "MAST4")
genes_of_interest <- toupper(genes_of_interest)

# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------

dat_dir <- "C:/Users/loril/Documents/multivariate_analysis/glioVis"

output_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/custom_genes/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------

tcga_pheno <- read.table(
  file.path(dat_dir, "2024-06-04_TCGA_GBM_pheno.txt"),
  sep = "\t", header = TRUE, check.names = FALSE
)

tcga_exp <- read.table(
  file.path(dat_dir, "2024-06-04_TCGA_GBM_expression.txt"),
  sep = "\t", header = TRUE, check.names = FALSE
) %>% data.frame()

rownames(tcga_pheno) <- tcga_pheno$Sample
rownames(tcga_exp)   <- tcga_exp$Sample

# -----------------------------------------------------------------------------
# ALIGN SAMPLES
# -----------------------------------------------------------------------------

common_samples <- intersect(tcga_pheno$Sample, rownames(tcga_exp))

tcga_pheno <- tcga_pheno[common_samples, ]
tcga_exp   <- tcga_exp[common_samples, ]

tcga_exp$Sample <- NULL
tcga_exp <- as.matrix(tcga_exp)

stopifnot(all.equal(rownames(tcga_exp), tcga_pheno$Sample))

# -----------------------------------------------------------------------------
# FILTER CLINICAL DATA
# -----------------------------------------------------------------------------

tcga_pheno <- tcga_pheno %>%
  filter(IDH1_status == "Wild-type" & !is.na(Gender))

tcga_exp <- tcga_exp[tcga_pheno$Sample, ]

# -----------------------------------------------------------------------------
# MERGE EXPRESSION INTO PHENOTYPE
# -----------------------------------------------------------------------------

tcga_pheno <- cbind(tcga_pheno, tcga_exp)

# -----------------------------------------------------------------------------
# KEEP ONLY GENES THAT EXIST
# -----------------------------------------------------------------------------

genes_found <- genes_of_interest[genes_of_interest %in% colnames(tcga_pheno)]
missing     <- setdiff(genes_of_interest, genes_found)

cat("Found genes:\n")
print(genes_found)

if (length(missing) > 0) {
  cat("Missing genes:\n")
  print(missing)
}

# -----------------------------------------------------------------------------
# FACTORIZE COVARIATES
# -----------------------------------------------------------------------------

factorize_cols <- function(df) {
  df$Recurrence  <- factor(df$Recurrence)
  df$Subtype     <- factor(df$Subtype)
  df$MGMT_status <- factor(df$MGMT_status)
  df
}

tcga_pheno <- factorize_cols(tcga_pheno)

# -----------------------------------------------------------------------------
# SPLIT BY SEX
# -----------------------------------------------------------------------------

tcga_male   <- tcga_pheno %>% filter(Gender == "Male")
tcga_female <- tcga_pheno %>% filter(Gender == "Female")

# -----------------------------------------------------------------------------
# BINARY EXPRESSION (MEDIAN SPLIT PER GROUP)
# -----------------------------------------------------------------------------

add_binary <- function(df, genes) {
  for (g in genes) {
    df[[paste0(g, "_binary")]] <- factor(
      ifelse(df[[g]] >= median(df[[g]], na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")
    )
  }
  df
}

tcga_pheno <- add_binary(tcga_pheno, genes_found)
tcga_male  <- add_binary(tcga_male, genes_found)
tcga_female<- add_binary(tcga_female, genes_found)

# -----------------------------------------------------------------------------
# DATASETS
# -----------------------------------------------------------------------------

datasets <- list(
  All     = tcga_pheno,
  Males   = tcga_male,
  Females = tcga_female
)

# -----------------------------------------------------------------------------
# MODEL SETTINGS
# -----------------------------------------------------------------------------

covariates <- c("Recurrence", "Age", "Subtype", "MGMT_status")

OS_var   <- "survival"
OS_event <- "status"

# -----------------------------------------------------------------------------
# STORAGE
# -----------------------------------------------------------------------------

all_long_results <- list()

# =============================================================================
# MAIN LOOP
# =============================================================================

gene_sets <- list(
  SelectedGenes = genes_found
)

for (set_name in names(gene_sets)) {
  
  genes <- gene_sets[[set_name]]
  
  cat("\n===", set_name, "— genes:", length(genes), "===\n")
  
  out_sub  <- file.path(output_dir, set_name)
  plot_sub <- file.path(out_sub, "forest_plots")
  km_sub   <- file.path(out_sub, "km_plots")
  
  dir.create(out_sub,  showWarnings = FALSE, recursive = TRUE)
  dir.create(plot_sub, showWarnings = FALSE, recursive = TRUE)
  dir.create(km_sub,   showWarnings = FALSE, recursive = TRUE)
  
  for (dataset_name in names(datasets)) {
    
    pheno <- datasets[[dataset_name]]
    
    cat("Running Cox:", set_name, "x", dataset_name, "\n")
    
    suffix <- paste0(substr(set_name, 1, 1), substr(dataset_name, 1, 1))
    
    # -------------------------------------------------------------------------
    # LONG COX
    # -------------------------------------------------------------------------
    
    long_res <- run_cox_for_group(
      pheno_data = pheno,
      genes      = genes,
      OS_var     = OS_var,
      OS_event   = OS_event,
      covariates = covariates
    )
    
    write.csv(long_res,
              file.path(out_sub, paste0("cox_long_", set_name, "_", dataset_name, ".csv")),
              row.names = FALSE
    )
    
    all_long_results[[paste0(set_name, "_", dataset_name)]] <- long_res
    
    
    # -------------------------------------------------------------------------
    # FOREST PLOTS (UNCHANGED FUNCTION USAGE)
    # -------------------------------------------------------------------------
    
    for (gene in genes) {
      tryCatch({
        p <- create_forest_plot(
          gene_name     = gene,
          survival_data = long_res,
          title_suffix  = paste0(" — ", set_name, " | ", dataset_name)
        )
        
        ggsave(
          filename = file.path(plot_sub,
                               paste0("forest_", gene, "_", dataset_name, ".pdf")),
          plot = p,
          width = 6,
          height = 4
        )
      }, error = function(e) {
        warning(paste("Forest failed:", gene, dataset_name))
      })
    }
    
    # -------------------------------------------------------------------------
    # KM PLOTS
    # -------------------------------------------------------------------------
    
    for (gene in genes) {
      
      binary_col <- paste0(gene, "_binary")
      if (!binary_col %in% colnames(pheno)) next
      
      tryCatch({
        
        fit <- survfit(as.formula(
          paste0("Surv(", OS_var, ", ", OS_event, ") ~ ", binary_col)
        ), data = pheno)
        
        p1 <- ggsurvplot(
          fit,
          data = pheno,
          risk.table = TRUE,
          legend.labs = c("Low", "High"),
          title = paste0(gene, " — ", set_name, " | ", dataset_name),
          ggtheme = theme_bw()
        )
        
        ggsave(
          filename = file.path(km_sub,
                               paste0("KM_", gene, "_", dataset_name, ".pdf")),
          plot = p1$plot,
          width = 7,
          height = 5
        )
        
      }, error = function(e) {
        warning(paste("KM failed:", gene, dataset_name))
      })
    }
  }
}

# -----------------------------------------------------------------------------
# SAVE ALL RESULTS
# -----------------------------------------------------------------------------

saveRDS(all_long_results,
        file.path(output_dir, "all_long_cox_results.rds"))

cat("\n=== DONE ===\n")