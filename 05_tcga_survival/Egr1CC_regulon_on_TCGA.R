# =============================================================================
# TCGA Survival Analysis - Z-score gene set method
# CC peaks x DE genes overlap sets (Male and Female separately)
# =============================================================================

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# Must run 01_load_and_prepare_tcga_data.R first
if (!exists("tcga_pheno")) stop("Run 01_load_and_prepare_tcga_data.R first.")

base_dir   <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"
output_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/enrichment/Egr1CCvsKD/"

# =============================================================================
# STEP 1: Read in overlap gene sets
# =============================================================================

overlap_Male_KD <- read.csv(file.path(output_dir, "overlap_MaleCC_vs_MaleKD_enriched.csv"))
overlap_Male_WT <- read.csv(file.path(output_dir, "overlap_MaleCC_vs_MaleWT_enriched.csv"))
overlap_Female_KD <- read.csv(file.path(output_dir, "overlap_FemaleCC_vs_FemaleKD_enriched.csv"))
overlap_Female_WT <- read.csv(file.path(output_dir, "overlap_FemaleCC_vs_FemaleWT_enriched.csv"))

# Extract gene symbols and convert to uppercase to match TCGA columns
geneset_Male_KD   <- unique(toupper(overlap_Male_KD$SYMBOL))
geneset_Male_WT   <- unique(toupper(overlap_Male_WT$SYMBOL))
geneset_Female_KD <- unique(toupper(overlap_Female_KD$SYMBOL))
geneset_Female_WT <- unique(toupper(overlap_Female_WT$SYMBOL))

cat("Male KD gene set:  ", length(geneset_Male_KD),   "genes\n")
cat("Male WT gene set:  ", length(geneset_Male_WT),   "genes\n")
cat("Female KD gene set:", length(geneset_Female_KD), "genes\n")
cat("Female WT gene set:", length(geneset_Female_WT), "genes\n")

# =============================================================================
# STEP 2: Check which genes are present in TCGA expression data
# =============================================================================

check_genes <- function(geneset, label) {
  found   <- geneset[geneset %in% colnames(tcga_pheno)]
  missing <- setdiff(geneset, colnames(tcga_pheno))
  cat("\n", label, "\n")
  cat("  Found:  ", length(found),   "—", paste(found,   collapse = ", "), "\n")
  cat("  Missing:", length(missing), "—", paste(missing, collapse = ", "), "\n")
  return(found)
}

geneset_Male_KD_found   <- check_genes(geneset_Male_KD,   "Male CC x KD-enriched:")
geneset_Male_WT_found   <- check_genes(geneset_Male_WT,   "Male CC x WT-enriched:")
geneset_Female_KD_found <- check_genes(geneset_Female_KD, "Female CC x KD-enriched:")
geneset_Female_WT_found <- check_genes(geneset_Female_WT, "Female CC x WT-enriched:")

# =============================================================================
# STEP 3: Z-score function
# =============================================================================

run_zscore_survival <- function(pheno, geneset, dataset_name, geneset_name, output_dir) {
  
  # Check enough genes available
  genes_available <- geneset[geneset %in% colnames(pheno)]
  if (length(genes_available) < 2) {
    cat("Skipping", geneset_name, "x", dataset_name, "— fewer than 2 genes available\n")
    return(NULL)
  }
  
  cat("\nRunning:", geneset_name, "x", dataset_name,
      "(", length(genes_available), "genes )\n")
  
  # Z-score each gene then take mean across gene set
  gene_mat  <- as.matrix(pheno[, genes_available, drop = FALSE])
  gene_mat_z <- scale(gene_mat)
  pheno$geneset_score <- rowMeans(gene_mat_z, na.rm = TRUE)
  
  # Median split
  pheno$geneset_score_binary <- factor(
    ifelse(pheno$geneset_score >= median(pheno$geneset_score, na.rm = TRUE),
           "High", "Low"),
    levels = c("Low", "High")
  )
  
  # Cox regression
  cox_res <- coxph(
    Surv(survival, status) ~ geneset_score + Recurrence + Age + Subtype + MGMT_status,
    data = pheno
  )
  cox_sum <- summary(cox_res)
  
  hr    <- round(cox_sum$coefficients["geneset_score", "exp(coef)"], 3)
  lower <- round(cox_sum$conf.int["geneset_score",    "lower .95"],  3)
  upper <- round(cox_sum$conf.int["geneset_score",    "upper .95"],  3)
  pval  <- signif(cox_sum$coefficients["geneset_score", "Pr(>|z|)"], 3)
  
  cat("  HR:", hr, "(", lower, "-", upper, ") p =", pval, "\n")
  
  # KM plot
  fit <- survfit(Surv(survival, status) ~ geneset_score_binary, data = pheno)
  
  p <- ggsurvplot(
    fit,
    data             = pheno,
    risk.table       = TRUE,
    pval             = FALSE,
    conf.int         = FALSE,
    xlab             = "OS (months)",
    legend.title     = "Gene Set Score",
    legend.labs      = c("Low", "High"),
    legend           = "top",
    surv.median.line = "hv",
    palette          = "npg",
    ggtheme          = theme_bw(),
    title            = paste0(geneset_name, " — ", dataset_name)
  )
  
  # Annotate HR and p-value on plot
  p$plot <- p$plot +
    ggplot2::annotate("text", x = 30, y = 0.90, hjust = 0, size = 3.5,
                      label = paste0("HR = ", hr, " (", lower, "–", upper, ")")) +
    ggplot2::annotate("text", x = 30, y = 0.80, hjust = 0, size = 3.5,
                      label = paste0("Cox p = ", pval))
  
  ggsave(
    filename = file.path(output_dir, paste0("KM_", geneset_name, "_", dataset_name, ".pdf")),
    plot = p$plot, width = 7, height = 5
  )
  
  # Return summary row
  data.frame(
    Geneset      = geneset_name,
    Dataset      = dataset_name,
    N_genes      = length(genes_available),
    Genes        = paste(genes_available, collapse = ", "),
    HR           = hr,
    Lower95      = lower,
    Upper95      = upper,
    Pvalue       = pval
  )
}

# =============================================================================
# STEP 4: Run all combinations
# =============================================================================

datasets <- list(
  All     = tcga_pheno,
  Males   = tcga_pheno_male,
  Females = tcga_pheno_female
)

genesets <- list(
  MaleCC_KD_enriched   = geneset_Male_KD_found,
  MaleCC_WT_enriched   = geneset_Male_WT_found,
  FemaleCC_KD_enriched = geneset_Female_KD_found,
  FemaleCC_WT_enriched = geneset_Female_WT_found
)

all_results <- list()

for (geneset_name in names(genesets)) {
  for (dataset_name in names(datasets)) {
    
    result <- tryCatch(
      run_zscore_survival(
        pheno        = datasets[[dataset_name]],
        geneset      = genesets[[geneset_name]],
        dataset_name = dataset_name,
        geneset_name = geneset_name,
        output_dir   = output_dir
      ),
      error = function(e) {
        warning(paste("Failed:", geneset_name, "x", dataset_name, "—", e$message))
        NULL
      }
    )
    
    if (!is.null(result)) {
      all_results[[paste0(geneset_name, "_", dataset_name)]] <- result
    }
  }
}

# =============================================================================
# STEP 5: Save summary table
# =============================================================================

summary_df <- bind_rows(all_results)
print(summary_df)

write.csv(summary_df,
          file.path(output_dir, "TCGA_zscore_survival_summary.csv"),
          row.names = FALSE)

cat("\nDone. Results saved to:", output_dir, "\n")

# Make forest plot
# =============================================================================
# Forest plot — gene set level results (4 gene sets as rows)
# =============================================================================

library(ggplot2)
library(dplyr)

output_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/enrichment/Egr1CCvsKD/"

# =============================================================================
# INPUT: your results in the exact order you specified
# =============================================================================

results_df <- data.frame(
  Geneset  = c("MaleCC_WT_enriched", "MaleCC_KD_enriched", 
               "FemaleCC_WT_enriched", "FemaleCC_KD_enriched"),
  Dataset  = c("Males", "Males", "Females", "Females"),
  HR       = c(2.734, 1.370, 2.358, 2.744),
  Lower95  = c(1.312, 0.403, 0.743, 0.731),
  Upper95  = c(5.700, 4.658, 7.488, 10.293),
  Pvalue   = c(0.00727, 0.614, 0.146, 0.135),
  stringsAsFactors = FALSE
)

# Label for plot — combine geneset and dataset
results_df$Label <- paste0(results_df$Geneset, "\n(", results_df$Dataset, ")")

# Fix order for plot (bottom to top in ggplot)
results_df$Label <- factor(results_df$Label, levels = rev(results_df$Label))

# Significance marker
results_df$Significant <- ifelse(results_df$Pvalue < 0.05, "p < 0.05", "p ≥ 0.05")

# =============================================================================
# FOREST PLOT
# =============================================================================

p <- ggplot(results_df, aes(x = HR, y = Label, color = Significant)) +
  
  # Vertical line at HR = 1 (no effect)
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  
  # Confidence interval lines
  geom_errorbarh(aes(xmin = Lower95, xmax = Upper95),
                 height = 0.2, linewidth = 0.8) +
  
  # HR point
  geom_point(size = 4) +
  
  # P-value annotation on the right
  geom_text(aes(x = max(Upper95) + 1,
                label = paste0("p = ", signif(Pvalue, 3))),
            hjust = 0, size = 3.5, color = "black") +
  
  # Colors
  scale_color_manual(values = c("p < 0.05" = "#E41A1C", "p ≥ 0.05" = "grey40")) +
  
  # Axis labels
  labs(
    x     = "Hazard Ratio (95% CI)",
    y     = NULL,
    title = "Gene Set Survival Analysis — TCGA GBM",
    color = NULL
  ) +
  
  # Expand x axis to make room for p-value text
  scale_x_continuous(limits = c(0, max(results_df$Upper95) + 4)) +
  
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 12),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 10),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

print(p)

ggsave(
  file.path(output_dir, "forest_plot_genesets.pdf"),
  plot = p, width = 8, height = 5
)

cat("Saved to:", output_dir, "\n")

