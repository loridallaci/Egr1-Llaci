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

source("05_tcga_survival/utils.R")

# --- Paths --------------------------------------------------------------------

# original (author's machine): "C:/Users/loril/Documents/multivariate_analysis/glioVis"
dat_dir ="data"
# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"
output_dir <- "output/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"
motif_directory <- "output/"

# --- Paths --------------------------------------------------------------------
#dat_dir    <- "data/"           # user puts their data here
#output_dir <- "tcga_survival_results_GitHub/"        # outputs go here
#dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load phenotype ------------------------------------------------------------

tcga_pheno <- read.table(
  file.path(dat_dir, "2024-06-04_TCGA_GBM_pheno.txt"),
  sep = "\t", header = TRUE, check.names = FALSE
)

rownames(tcga_pheno) <- tcga_pheno$Sample

# --- Load expression -----------------------------------------------------------

tcga_exp <- read.table(
  file.path(dat_dir, "2024-06-04_TCGA_GBM_expression.txt"),
  sep = "\t", header = TRUE, check.names = FALSE
) %>% data.frame()

rownames(tcga_exp) <- tcga_exp$Sample

cat("Expression matrix dimensions:", dim(tcga_exp), "\n")

# --- Remove Sample column from expression -------------------------------------

tcga_exp <- tcga_exp %>%
  select(-Sample) %>%
  as.matrix()

# --- Align samples (CRITICAL) --------------------------------------------------

common_samples <- intersect(tcga_pheno$Sample, rownames(tcga_exp))

tcga_pheno <- tcga_pheno[common_samples, ]
tcga_exp   <- tcga_exp[common_samples, ]

stopifnot(all.equal(rownames(tcga_exp), tcga_pheno$Sample))

cat("Matched samples:", nrow(tcga_pheno), "\n")

# --- Filter clinical group (IDH WT + known sex) --------------------------------

tcga_pheno <- tcga_pheno %>%
  filter(IDH1_status == "Wild-type" & !is.na(Gender))

tcga_exp <- tcga_exp[tcga_pheno$Sample, ]

stopifnot(all.equal(rownames(tcga_exp), tcga_pheno$Sample))

cat("Samples after filtering:", nrow(tcga_pheno), "\n")
print(table(tcga_pheno$Gender))

# --- ADD ALL GENES INTO PHENOTYPE (THIS FIXES YOUR MAIN ERROR) -----------------

tcga_pheno <- cbind(tcga_pheno, tcga_exp)

genes <- colnames(tcga_exp)

# --- OPTIONAL: Summary of ALL genes (WARNING: large) --------------------------

# apply(tcga_exp, 2, summary)

# --- Sex subsets ---------------------------------------------------------------

tcga_pheno_female <- tcga_pheno %>% filter(Gender == "Female")
tcga_pheno_male   <- tcga_pheno %>% filter(Gender == "Male")

cat("Female samples:", nrow(tcga_pheno_female), "\n")
cat("Male samples:", nrow(tcga_pheno_male), "\n")

# --- Binary expression (ALL GENES) --------------------------------------------

# WARNING: this creates ~12k new columns → heavy but matches your request

for (gene in genes) {
  tcga_pheno[[paste0(gene, "_binary")]] <-
    factor(ifelse(tcga_pheno[[gene]] >= median(tcga_pheno[[gene]], na.rm = TRUE),
                  "High", "Low"),
           levels = c("Low", "High"))
}

# --- Factor covariates --------------------------------------------------------

factorize_cols <- function(df) {
  df$Recurrence  <- factor(df$Recurrence)
  df$Subtype     <- factor(df$Subtype)
  df$MGMT_status <- factor(df$MGMT_status)
  return(df)
}

tcga_pheno        <- factorize_cols(tcga_pheno)
tcga_pheno_male   <- factorize_cols(tcga_pheno_male)
tcga_pheno_female <- factorize_cols(tcga_pheno_female)

message("Data preparation complete. Ready for multivariate analysis.")

saveRDS(tcga_pheno, file.path(output_dir, "tcga_pheno_maleANDfemale_allgenes_updated.rds"))
saveRDS(tcga_pheno_male, file.path(output_dir, "tcga_pheno_male_allgenes_updated.rds"))
saveRDS(tcga_pheno_female, file.path(output_dir, "tcga_pheno_female_allgenes_updated.rds"))

#then filter genes I need

# --- Load motif files (Female + Male) + take only top 30 motifs from each dataset -----------------------------------------

motif_female <- read.csv(file.path(motif_directory, "F_all_motifs_updated.csv")) %>%
  arrange(p.adjust) %>%
  slice(1:30)

motif_male <- read.csv(file.path(motif_directory, "M_all_motifs_updated.csv")) %>%
  arrange(p.adjust) %>%
  slice(1:30)

# Add sex label
motif_female$sex <- "Female"
motif_male$sex   <- "Male"

# Combine
motif_all <- bind_rows(motif_female, motif_male)

# Use motif.name as TF name
motif_all$TF_name <- motif_all$motif.name

# Clean TF names (remove parentheses, trim)
motif_all$TF_name <- gsub("\\s*\\([^)]*\\)", "", motif_all$TF_name)
motif_all$TF_name <- trimws(motif_all$TF_name)

# Split composite TFs (e.g. Zic1::Zic2)
df_split <- motif_all %>%
  mutate(original_row = row_number()) %>%
  separate_rows(TF_name, sep = "::") %>%
  group_by(original_row) %>%
  mutate(duplicated = n() > 1) %>%
  ungroup() %>%
  mutate(
    dup_group    = if_else(duplicated,
                           paste0("dup_", match(original_row,
                                                unique(original_row[duplicated]))),
                           NA_character_),
    SYMBOL_UPPER = toupper(TF_name)
  )

# --- Gene list ----------------------------------------------------------------

geneID <- unique(toupper(df_split$TF_name))
cat("Total unique TFs:", length(geneID), "\n")

# Keep only TFs that exist in TCGA expression data
genes_found <- geneID[geneID %in% colnames(tcga_pheno)]
missing     <- geneID[!geneID %in% colnames(tcga_pheno)]

cat("TFs found in TCGA:", length(genes_found), "\n")
if (length(missing) > 0) cat("TFs NOT found:", paste(missing, collapse = ", "), "\n")

# Set geneID_included for all downstream scripts
geneID_included <- genes_found

# --- Binary expression columns (median split, per sex-specific dataset) -------
# Done separately per group so median reflects that group's distribution

add_binary_cols <- function(df, genes) {
  for (gene in genes) {
    if (gene %in% colnames(df)) {
      df[[paste0(gene, "_binary")]] <- factor(
        ifelse(df[[gene]] >= median(df[[gene]], na.rm = TRUE), "High", "Low"),
        levels = c("Low", "High")
      )
    }
  }
  return(df)
}

tcga_pheno        <- add_binary_cols(tcga_pheno,        geneID_included)
tcga_pheno_male   <- add_binary_cols(tcga_pheno_male,   geneID_included)
tcga_pheno_female <- add_binary_cols(tcga_pheno_female, geneID_included)

cat("Genes included in analysis:", length(geneID_included), "\n")
cat("Columns in tcga_pheno:",        ncol(tcga_pheno), "\n")
cat("Columns in tcga_pheno_male:",   ncol(tcga_pheno_male), "\n")
cat("Columns in tcga_pheno_female:", ncol(tcga_pheno_female), "\n")

# --- Sex-specific gene lists --------------------------------------------------

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

cat("Male motif genes found in TCGA:  ", length(geneID_included_M), "\n")
cat("Female motif genes found in TCGA:", length(geneID_included_F), "\n")

cat("Male motif genes:  ", paste(geneID_included_M, collapse=", "), "\n")
cat("Female motif genes:", paste(geneID_included_F, collapse=", "), "\n")
cat("Overlap:", length(intersect(geneID_included_M, geneID_included_F)), "\n")

