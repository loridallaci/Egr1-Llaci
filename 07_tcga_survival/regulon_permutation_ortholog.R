# =============================================================================
# Venn-UNIQUE regulons split by KD direction (Up / Down)  ->  TCGA GBM survival
#   multivariate gene-set z-score Cox  +  permutation test
# -----------------------------------------------------------------------------
# Gene-set construction (Venn-biased ∩ KD direction):
#   MaleUnique_Up    = Male_only   ∩  kd3_male_sig   == "Up after KD"     (151)
#   MaleUnique_Down  = Male_only   ∩  kd3_male_sig   == "Down after KD"   ( 51)
#   FemaleUnique_Up  = Female_only ∩  kd3_female_sig == "Up after KD"     (160)
#   FemaleUnique_Down= Female_only ∩  kd3_female_sig == "Down after KD"   ( 64)
#
# Pipeline same as Egr1CC_regulon_on_TCGA.R / _permutation.R, mirrored on
# four matched pairs (each set against its sex-matched TCGA subset).
#
# Self-bootstraps TCGA if not already loaded.
# Hardened numeric filter for the permutation gene pool.
# =============================================================================

library(survival)
library(survminer)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

set.seed(42)

# =============================================================================
# 0) Paths and settings
# =============================================================================
in_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output"

male_summary_file   <- file.path(in_dir,
                                 "Egr1CC_allGenes_DE_summary_with_GO_MALES_gRNA3_20kb_rerun_052726.csv")
female_summary_file <- file.path(in_dir,
                                 "Egr1CC_allGenes_DE_summary_with_GO_FEMALES_gRNA3_20kb_rerun_052726.csv")

# >>> OUTPUT DIRECTORY <<<
output_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/07_tcga_survival/output_regulon_survival"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("All outputs will be saved to:\n  ", output_dir, "\n\n")

# Paths from 01_load_and_prepare_tcga_data_updated.R
glio_dir     <- "C:/Users/loril/Documents/multivariate_analysis/glioVis"
tcga_rds_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output"

# Permutation settings
n_perm        <- 1000              # raise to 10000 for p < 0.001 resolution
scramble_mode <- "random_genome"   # "random_genome" or "reshuffle_within_regulons"

# Which _sig column defines Up/Down in each summary CSV
male_sig_col   <- "kd3_male_sig"
female_sig_col <- "kd3_female_sig"

# =============================================================================
# 1) Self-bootstrap TCGA pheno objects
# =============================================================================
prepare_tcga_data <- function() {
  rds_all <- file.path(tcga_rds_dir, "tcga_pheno_maleANDfemale_allgenes_updated.rds")
  rds_m   <- file.path(tcga_rds_dir, "tcga_pheno_male_allgenes_updated.rds")
  rds_f   <- file.path(tcga_rds_dir, "tcga_pheno_female_allgenes_updated.rds")
  if (all(file.exists(c(rds_all, rds_m, rds_f)))) {
    cat("Loading TCGA from RDS in\n  ", tcga_rds_dir, "\n")
    tp   <- readRDS(rds_all); tp_m <- readRDS(rds_m); tp_f <- readRDS(rds_f)
  } else {
    pheno_file <- file.path(glio_dir, "2024-06-04_TCGA_GBM_pheno.txt")
    exp_file   <- file.path(glio_dir, "2024-06-04_TCGA_GBM_expression.txt")
    if (!all(file.exists(c(pheno_file, exp_file))))
      stop("RDS missing in '", tcga_rds_dir, "' AND raw text files missing in '", glio_dir,
           "'.\nRun 01_load_and_prepare_tcga_data_updated.R first.")
    cat("Rebuilding TCGA from raw text files (slower)...\n")
    pheno   <- read.table(pheno_file, sep = "\t", header = TRUE, check.names = FALSE)
    rownames(pheno) <- pheno$Sample
    exp_mat <- read.table(exp_file,   sep = "\t", header = TRUE, check.names = FALSE)
    rownames(exp_mat) <- exp_mat$Sample
    exp_mat <- as.matrix(exp_mat[, setdiff(colnames(exp_mat), "Sample"), drop = FALSE])
    common  <- intersect(pheno$Sample, rownames(exp_mat))
    pheno   <- pheno[common, ]; exp_mat <- exp_mat[common, ]
    pheno   <- pheno[pheno$IDH1_status == "Wild-type" & !is.na(pheno$Gender), ]
    exp_mat <- exp_mat[pheno$Sample, ]
    tp   <- cbind(pheno, as.data.frame(exp_mat))
    tp_f <- tp[tp$Gender == "Female", ]; tp_m <- tp[tp$Gender == "Male", ]
  }
  fact <- function(df) {
    df$Recurrence  <- factor(df$Recurrence)
    df$Subtype     <- factor(df$Subtype)
    df$MGMT_status <- factor(df$MGMT_status); df
  }
  list(all = fact(tp), male = fact(tp_m), female = fact(tp_f))
}

needs_load <- !exists("tcga_pheno") || ncol(tcga_pheno) < 100 ||
  !exists("tcga_pheno_male") || !exists("tcga_pheno_female")
if (needs_load) {
  cur_cols <- if (exists("tcga_pheno")) ncol(tcga_pheno) else 0
  cat("Bootstrapping TCGA (current tcga_pheno has", cur_cols, "cols).\n")
  .tcga <- prepare_tcga_data()
  tcga_pheno        <- .tcga$all
  tcga_pheno_male   <- .tcga$male
  tcga_pheno_female <- .tcga$female
}
cat("TCGA ready: tcga_pheno", ncol(tcga_pheno), "cols |",
    nrow(tcga_pheno_male), "M |", nrow(tcga_pheno_female), "F\n\n")

# =============================================================================
# 2) Derive 4 gene sets: MaleUnique/FemaleUnique x Up/Down
# =============================================================================
male_df   <- read_csv(male_summary_file,   show_col_types = FALSE)
female_df <- read_csv(female_summary_file, show_col_types = FALSE)

clean_symbols <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & x != ""])
}
male_symbols   <- clean_symbols(male_df$SYMBOL)
female_symbols <- clean_symbols(female_df$SYMBOL)

male_only_mouse   <- setdiff(male_symbols,   female_symbols)
female_only_mouse <- setdiff(female_symbols, male_symbols)

# Build SYMBOL -> sig direction lookups
m_sig_lookup <- setNames(male_df[[male_sig_col]],   male_df$SYMBOL)
f_sig_lookup <- setNames(female_df[[female_sig_col]], female_df$SYMBOL)

geneset_MaleUnique_Up     <- male_only_mouse[  m_sig_lookup[male_only_mouse]   == "Up after KD"   ]
geneset_MaleUnique_Down   <- male_only_mouse[  m_sig_lookup[male_only_mouse]   == "Down after KD" ]
geneset_FemaleUnique_Up   <- female_only_mouse[f_sig_lookup[female_only_mouse] == "Up after KD"   ]
geneset_FemaleUnique_Down <- female_only_mouse[f_sig_lookup[female_only_mouse] == "Down after KD" ]

cat("Gene sets (mouse SYMBOLs):\n")
cat("  MaleUnique_Up    :", length(geneset_MaleUnique_Up),     "\n")
cat("  MaleUnique_Down  :", length(geneset_MaleUnique_Down),   "\n")
cat("  FemaleUnique_Up  :", length(geneset_FemaleUnique_Up),   "\n")
cat("  FemaleUnique_Down:", length(geneset_FemaleUnique_Down), "\n")

# Save gene lists to txt
.write <- function(x, nm)
  writeLines(x, file.path(output_dir, paste0("regulon_genes_", nm, "_052726.txt")))
.write(geneset_MaleUnique_Up,     "MaleUnique_Up")
.write(geneset_MaleUnique_Down,   "MaleUnique_Down")
.write(geneset_FemaleUnique_Up,   "FemaleUnique_Up")
.write(geneset_FemaleUnique_Down, "FemaleUnique_Down")

# =============================================================================
# 3) Match to TCGA columns (case-insensitive)
# =============================================================================
library(babelgene)
tcga_cols       <- colnames(tcga_pheno)
tcga_cols_upper <- toupper(tcga_cols)
upper_to_orig   <- setNames(tcga_cols, tcga_cols_upper)
upper_to_orig   <- upper_to_orig[!duplicated(names(upper_to_orig))]

## Proper mouse -> human ortholog conversion (babelgene), with an uppercase
## fallback for any mouse symbols babelgene cannot map; then map to TCGA columns.
mouse_to_human <- function(sym) {
  o <- tryCatch(orthologs(genes = sym, species = "mouse", human = FALSE),
                error = function(e) NULL)
  human  <- if (!is.null(o) && nrow(o)) toupper(o$human_symbol) else character(0)
  mapped <- if (!is.null(o) && nrow(o)) o$symbol               else character(0)
  c(human, toupper(setdiff(sym, mapped)))   # fallback: uppercase the unmapped
}
match_genes <- function(sym) {
  u <- unique(mouse_to_human(sym))
  unique(unname(upper_to_orig[u[u %in% names(upper_to_orig)]]))
}

genesets <- list(
  MaleUnique_Up     = match_genes(geneset_MaleUnique_Up),
  MaleUnique_Down   = match_genes(geneset_MaleUnique_Down),
  FemaleUnique_Up   = match_genes(geneset_FemaleUnique_Up),
  FemaleUnique_Down = match_genes(geneset_FemaleUnique_Down)
)

cat("\nGenes available in TCGA (after case-insensitive match):\n")
for (nm in names(genesets))
  cat(sprintf("  %-18s: %3d\n", nm, length(genesets[[nm]])))

if (sum(lengths(genesets)) == 0)
  stop("0 genes matched in TCGA - check that tcga_pheno has the expression matrix.")

# =============================================================================
# 4) z-score gene-set multivariate Cox
# =============================================================================
run_zscore_survival <- function(pheno, geneset, dataset_name, geneset_name,
                                make_km = TRUE) {
  genes_available <- geneset[geneset %in% colnames(pheno)]
  if (length(genes_available) < 2) {
    cat("  Skipping", geneset_name, "x", dataset_name, "- <2 genes\n"); return(NULL)
  }
  cat("  Running", geneset_name, "x", dataset_name,
      "(", length(genes_available), "genes )\n")
  
  gene_mat_z <- scale(as.matrix(pheno[, genes_available, drop = FALSE]))
  pheno$geneset_score <- rowMeans(gene_mat_z, na.rm = TRUE)
  pheno$geneset_score_binary <- factor(
    ifelse(pheno$geneset_score >= median(pheno$geneset_score, na.rm = TRUE),
           "High", "Low"), levels = c("Low", "High"))
  
  cox_res <- coxph(
    Surv(survival, status) ~ geneset_score + Recurrence + Age + Subtype + MGMT_status,
    data = pheno)
  cs <- summary(cox_res)
  hr    <- round(cs$coefficients["geneset_score", "exp(coef)"], 3)
  lower <- round(cs$conf.int["geneset_score",    "lower .95"],  3)
  upper <- round(cs$conf.int["geneset_score",    "upper .95"],  3)
  pval  <- signif(cs$coefficients["geneset_score", "Pr(>|z|)"], 3)
  cat("    HR =", hr, "(", lower, "-", upper, ") Cox p =", pval, "\n")
  
  if (make_km) {
    fit <- survfit(Surv(survival, status) ~ geneset_score_binary, data = pheno)
    p <- ggsurvplot(fit, data = pheno, risk.table = TRUE, pval = FALSE,
                    conf.int = FALSE, xlab = "OS (months)",
                    legend.title = "Gene Set Score",
                    legend.labs = c("Low", "High"), legend = "top",
                    surv.median.line = "hv", palette = "npg",
                    ggtheme = theme_bw(),
                    title = paste0(geneset_name, " - ", dataset_name))
    p$plot <- p$plot +
      ggplot2::annotate("text", x = 30, y = 0.90, hjust = 0, size = 3.5,
                        label = paste0("HR = ", hr, " (", lower, "-", upper, ")")) +
      ggplot2::annotate("text", x = 30, y = 0.80, hjust = 0, size = 3.5,
                        label = paste0("Cox p = ", pval))
    ggsave(file.path(output_dir,
                     paste0("KM_", geneset_name, "_", dataset_name, ".pdf")),
           plot = p$plot, width = 7, height = 5)
  }
  
  data.frame(Geneset = geneset_name, Dataset = dataset_name,
             N_genes = length(genes_available),
             Genes = paste(genes_available, collapse = ", "),
             HR = hr, Lower95 = lower, Upper95 = upper, Pvalue = pval,
             stringsAsFactors = FALSE)
}

datasets <- list(All = tcga_pheno, Males = tcga_pheno_male, Females = tcga_pheno_female)

all_results <- list()
for (gn in names(genesets)) {
  cat("\n===", gn, "===\n")
  for (dn in names(datasets)) {
    res <- tryCatch(run_zscore_survival(datasets[[dn]], genesets[[gn]], dn, gn),
                    error = function(e) { warning(e$message); NULL })
    if (!is.null(res)) all_results[[paste0(gn, "_", dn)]] <- res
  }
}

summary_df <- bind_rows(all_results)
print(summary_df[, c("Geneset","Dataset","N_genes","HR","Lower95","Upper95","Pvalue")])
write.csv(summary_df,
          file.path(output_dir,
                    "TCGA_regulon_UNIQUE_UpDown_zscore_survival_summary.csv"),
          row.names = FALSE)

# =============================================================================
# 5) Forest plot of REAL results (sex-matched pairs)
# =============================================================================
forest_order <- c("MaleUnique_Down", "MaleUnique_Up",
                  "FemaleUnique_Down", "FemaleUnique_Up")

real_forest <- summary_df %>%
  filter((Geneset %in% c("MaleUnique_Up","MaleUnique_Down")    & Dataset == "Males") |
           (Geneset %in% c("FemaleUnique_Up","FemaleUnique_Down") & Dataset == "Females")) %>%
  mutate(Geneset = factor(Geneset, levels = forest_order)) %>%
  arrange(Geneset) %>%
  mutate(Label = paste0(Geneset, "\n(", Dataset, ")"),
         Label = factor(Label, levels = rev(Label)),
         Significant = ifelse(Pvalue < 0.05, "p < 0.05", "p >= 0.05"))

p_real <- ggplot(real_forest, aes(x = HR, y = Label, color = Significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = Lower95, xmax = Upper95), height = 0.2, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(x = max(Upper95) + 1, label = paste0("p = ", signif(Pvalue, 3))),
            hjust = 0, size = 3.5, color = "black") +
  scale_color_manual(values = c("p < 0.05" = "#E41A1C", "p >= 0.05" = "grey40")) +
  scale_x_continuous(limits = c(0, max(real_forest$Upper95) + 4)) +
  labs(x = "Hazard Ratio (95% CI)", y = NULL, color = NULL,
       title = "Venn-UNIQUE x KD direction regulons - TCGA GBM survival") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "top", panel.grid.minor = element_blank())
print(p_real)
ggsave(file.path(output_dir, "forest_plot_regulon_UNIQUE_UpDown.pdf"),
       plot = p_real, width = 8, height = 5)

# =============================================================================
# 6) Permutation test  ->  empirical p-values  (hardened numeric filter)
# =============================================================================
is_numeric_strict <- function(x) is.numeric(x) && !is.factor(x) && !is.character(x)
is_numeric_in <- function(df, cols)
  vapply(df[, cols, drop = FALSE], is_numeric_strict, logical(1))

if (exists("tcga_exp")) {
  gene_pool0 <- intersect(colnames(tcga_exp), colnames(tcga_pheno))
} else {
  non_gene <- c("survival","status","Recurrence","Age","Subtype","MGMT_status",
                "Sample","Gender","IDH1_status","Sex",
                "Histology","CIMP_status","Subtype_Verhaak_2010","Therapy_Class")
  gene_pool0 <- setdiff(colnames(tcga_pheno), non_gene)
  gene_pool0 <- gene_pool0[!grepl("_binary$|_score$", gene_pool0)]
}
in_all <- Reduce(intersect, list(colnames(tcga_pheno),
                                 colnames(tcga_pheno_male),
                                 colnames(tcga_pheno_female)))
gene_pool0 <- intersect(gene_pool0, in_all)

keep <- is_numeric_in(tcga_pheno,        gene_pool0) &
  is_numeric_in(tcga_pheno_male,   gene_pool0) &
  is_numeric_in(tcga_pheno_female, gene_pool0)
gene_pool <- gene_pool0[keep]
cat("\nGene pool: ", length(gene_pool0), "candidates ->",
    length(gene_pool), "numeric in all three pheno objects\n")

if (scramble_mode == "random_genome") {
  scramble_pool <- gene_pool
} else if (scramble_mode == "reshuffle_within_regulons") {
  scramble_pool <- intersect(unlist(genesets, use.names = FALSE), gene_pool)
} else stop("scramble_mode must be 'random_genome' or 'reshuffle_within_regulons'")
cat("Scramble pool (", scramble_mode, "):", length(scramble_pool), "genes\n")

set_sizes <- sapply(genesets, length)
print(set_sizes)

runs <- list(
  list(geneset = "MaleUnique_Up",     dataset = "Males"),
  list(geneset = "MaleUnique_Down",   dataset = "Males"),
  list(geneset = "FemaleUnique_Up",   dataset = "Females"),
  list(geneset = "FemaleUnique_Down", dataset = "Females")
)

perm_one <- function(pheno, n_genes) {
  pool <- scramble_pool[scramble_pool %in% colnames(pheno)]
  if (n_genes < 2 || length(pool) < n_genes) return(NA_real_)
  genes <- sample(pool, n_genes)
  mat <- as.matrix(pheno[, genes, drop = FALSE])
  if (!is.numeric(mat)) return(NA_real_)   # defensive
  pheno$geneset_score <- rowMeans(scale(mat), na.rm = TRUE)
  fit <- tryCatch(
    coxph(Surv(survival, status) ~ geneset_score + Recurrence + Age + Subtype + MGMT_status,
          data = pheno),
    error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  unname(summary(fit)$coefficients["geneset_score", "exp(coef)"])
}

null_summary <- list(); null_long <- list()

for (r in runs) {
  key     <- paste0(r$geneset, "_", r$dataset)
  n_genes <- set_sizes[[r$geneset]]
  pheno   <- datasets[[r$dataset]]
  ix <- which(summary_df$Geneset == r$geneset & summary_df$Dataset == r$dataset)[1]
  if (is.na(ix)) { cat("No real result for", key, "- skipping\n"); next }
  real_hr    <- summary_df$HR[ix]
  real_cox_p <- summary_df$Pvalue[ix]
  
  cat("\nPermuting", key, "( n_perm =", n_perm, ", set size =", n_genes, ")...\n")
  null_hr <- replicate(n_perm, perm_one(pheno, n_genes))
  null_hr <- null_hr[is.finite(null_hr)]
  n_valid <- length(null_hr)
  if (n_valid == 0) { cat("  no valid permutations\n"); next }
  
  emp_p_one <- if (real_hr >= 1)
    (sum(null_hr >= real_hr) + 1) / (n_valid + 1)
  else
    (sum(null_hr <= real_hr) + 1) / (n_valid + 1)
  emp_p_two <- (sum(abs(log(null_hr)) >= abs(log(real_hr))) + 1) / (n_valid + 1)
  
  cat(sprintf("  real HR = %.3f | null median = %.3f | 1-sided p = %.4f | 2-sided p = %.4f\n",
              real_hr, median(null_hr), emp_p_one, emp_p_two))
  
  null_summary[[key]] <- data.frame(
    Geneset = r$geneset, Dataset = r$dataset,
    Real_HR = real_hr, Real_Cox_P = real_cox_p,
    Null_median_HR = round(median(null_hr), 3),
    Empirical_P_onesided = round(emp_p_one, 4),
    Empirical_P_twosided = round(emp_p_two, 4),
    N_perm = n_valid, stringsAsFactors = FALSE)
  null_long[[key]] <- data.frame(
    Pairing = paste0(r$geneset, "\n(", r$dataset, ")"),
    HR = null_hr, Real_HR = real_hr, stringsAsFactors = FALSE)
}

null_summary_df <- bind_rows(null_summary)
print(null_summary_df)
write.csv(null_summary_df,
          file.path(output_dir,
                    "TCGA_regulon_UNIQUE_UpDown_PERMUTATION_null.csv"),
          row.names = FALSE)

# =============================================================================
# 7) Null distribution + permutation forest plot
# =============================================================================
null_long_df <- bind_rows(null_long)

pval_labels <- null_summary_df %>%
  mutate(Pairing = paste0(Geneset, "\n(", Dataset, ")"),
         label   = paste0("1-sided p = ", signif(Empirical_P_onesided, 3),
                          "\n2-sided p = ", signif(Empirical_P_twosided, 3)))

p_null <- ggplot(null_long_df, aes(x = HR)) +
  geom_histogram(bins = 40, fill = "grey70", color = "white") +
  geom_vline(aes(xintercept = Real_HR), color = "#E41A1C", linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(data = pval_labels, aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.3, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ Pairing, scales = "free", ncol = 2) +
  labs(x = "Hazard Ratio of random gene sets", y = "Count",
       title = paste0("Permutation null vs real UNIQUE-UpDown regulon HR  (",
                      n_perm, " random sets, red = real)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(p_null)
ggsave(file.path(output_dir, "permutation_null_distributions.pdf"),
       plot = p_null, width = 9, height = 6)

perm_forest <- null_long_df %>%
  group_by(Pairing) %>%
  summarise(null_median = median(HR),
            null_lo     = quantile(HR, 0.025),
            null_hi     = quantile(HR, 0.975),
            real_hr     = Real_HR[1], .groups = "drop") %>%
  left_join(null_summary_df %>%
              transmute(Pairing   = paste0(Geneset, "\n(", Dataset, ")"),
                        emp_p_one = Empirical_P_onesided),
            by = "Pairing") %>%
  # Pull the real HR's 95% CI from the multivariate Cox summary so the paper
  # figure shows both the null distribution AND the precision of the real estimate
  left_join(summary_df %>%
              transmute(Pairing = paste0(Geneset, "\n(", Dataset, ")"),
                        real_lo = Lower95, real_hi = Upper95),
            by = "Pairing")

row_order <- paste0(forest_order, "\n(",
                    c("Males","Males","Females","Females"), ")")
perm_forest$Pairing <- factor(perm_forest$Pairing, levels = rev(row_order))

# Pretty row labels for the paper figure
relabel <- c("MaleUnique_Down\n(Males)"     = "Male regulon\n(Down after KD)",
             "MaleUnique_Up\n(Males)"       = "Male regulon\n(Up after KD)",
             "FemaleUnique_Up\n(Females)"   = "Female regulon\n(Up after KD)",
             "FemaleUnique_Down\n(Females)" = "Female regulon\n(Down after KD)")
perm_forest$PrettyLabel <- factor(
  relabel[as.character(perm_forest$Pairing)],
  levels = rev(relabel[row_order]))

# Significance stars from one-sided empirical p
perm_forest$sig_star <- with(perm_forest, dplyr::case_when(
  emp_p_one < 0.001 ~ "***",
  emp_p_one < 0.01  ~ "**",
  emp_p_one < 0.05  ~ "*",
  TRUE              ~ ""))

# Right-side annotation: HR + 95% CI + 1-sided p + star
perm_forest$plabel <- with(perm_forest, sprintf(
  "HR = %.2f (%.2f\u2013%.2f)\np = %.3g %s",
  real_hr, real_lo, real_hi, emp_p_one, sig_star))

x_text <- max(c(perm_forest$null_hi, perm_forest$real_hi), na.rm = TRUE) * 1.06

# ---- Paper-quality permutation forest plot ----
p_perm_forest <- ggplot(perm_forest, aes(y = PrettyLabel)) +
  # Reference at HR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55",
             linewidth = 0.5) +
  # Null distribution: 2.5-97.5 percentile bar
  geom_errorbarh(aes(xmin = null_lo, xmax = null_hi),
                 height = 0.22, color = "grey65", linewidth = 1.2) +
  # Real HR's 95% CI (prominent red error bar with caps, drawn for every row)
  geom_errorbarh(aes(xmin = real_lo, xmax = real_hi),
                 height = 0.18, color = "#C0392B", linewidth = 1.1) +
  # Null median (small grey dot)
  geom_point(aes(x = null_median,
                 fill = paste0("Permutation null (n = ", n_perm, ")")),
             shape = 21, size = 3, color = "grey35", stroke = 0.4) +
  # Real regulon HR (red diamond, larger, prominent)
  geom_point(aes(x = real_hr, fill = "Real Egr1 regulon"),
             shape = 23, size = 4.5, color = "black", stroke = 0.5) +
  # HR + CI + p label on the right
  geom_text(aes(x = x_text, label = plabel), hjust = 0,
            size = 3.2, color = "black", lineheight = 0.95,
            family = "sans") +
  # Sig star, slightly bigger and red, just past the diamond
  geom_text(aes(x = real_hr, y = PrettyLabel, label = sig_star),
            nudge_y = 0.30, size = 5.5, color = "#C0392B",
            fontface = "bold", family = "sans") +
  
  scale_fill_manual(
    name = NULL,
    values = setNames(c("grey65", "#C0392B"),
                      c(paste0("Permutation null (n = ", n_perm, ")"),
                        "Real Egr1 regulon")),
    guide = guide_legend(override.aes = list(
      shape  = c(21, 23),
      size   = c(3, 4.5),
      stroke = c(0.4, 0.5)))) +
  
  scale_x_continuous(
    name = "Hazard Ratio (95% CI)",
    limits = c(0, x_text * 1.55),
    breaks = scales::pretty_breaks(n = 6),
    expand = c(0, 0)) +
  
  labs(y = NULL,
       title = "Sex-biased Egr1 regulons vs. survival in TCGA GBM",
       subtitle = paste0("Multivariate Cox (z-score gene-set score + Recurrence + Age + Subtype + MGMT)  |  ",
                         "grey bar = 2.5\u201397.5th percentile of ", n_perm,
                         " size-matched random gene sets")) +
  
  theme_classic(base_family = "sans", base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0,
                                    margin = margin(b = 3)),
    plot.subtitle    = element_text(size = 8.5, color = "grey25",
                                    margin = margin(b = 12)),
    axis.text.y      = element_text(size = 10.5, color = "black",
                                    lineheight = 1.05),
    axis.text.x      = element_text(size = 10, color = "black"),
    axis.title.x     = element_text(size = 10.5, color = "black",
                                    margin = margin(t = 8)),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length= unit(0.18, "cm"),
    panel.grid       = element_blank(),
    legend.position  = "top",
    legend.justification = "left",
    legend.text      = element_text(size = 9.5),
    legend.key       = element_blank(),
    legend.box.margin = margin(b = 4),
    plot.margin      = margin(t = 14, r = 14, b = 14, l = 14))

print(p_perm_forest)

# Save both PDF (vector, for the paper) and high-res PNG (for slides)
ggsave(file.path(output_dir, "forest_plot_permutation_null.pdf"),
       plot = p_perm_forest, width = 9, height = 4.5)
ggsave(file.path(output_dir, "forest_plot_permutation_null.png"),
       plot = p_perm_forest, width = 9, height = 4.5, dpi = 600)

cat("\n=== DONE ===\nAll outputs saved to:\n  ", output_dir, "\n")