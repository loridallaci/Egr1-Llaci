## =====================================================================
## STEP 3:  Full per-gene CC∩DE summary table (one CSV per sex)
##
## Rebuilds the master summary like
##   Egr1CC_allGenes_DE_summary_with_GO_<SEX>_gRNA3_20kb_..._ADDED_HR_MANUALLY.xlsx
## but fully programmatically. Row set = the CC∩DE genes for that sex
##   (Male = 250, Female = 272).
##
## Columns produced (all auto-generated, no manual steps):
##   SYMBOL
##   in_male_CC, in_female_CC                         (CC∩DE membership flags)
##   kd3_male_log2FoldChange / pvalue / padj / sig    (Egr1 gRNA3 KD vs Neg1, Males)
##   kd3_female_log2FoldChange / pvalue / padj / sig  (Egr1 gRNA3 KD vs Neg1, Females)
##   description, gene_biotype, GO_Term_Description    (biomaRt mouse annotation)
##   Multivariate Cox in TCGA MALE  patients: Overall_pvalue_M + per-term HR/pvalue
##   Multivariate Cox in TCGA FEMALE patients: Overall_pvalue_F + per-term HR/pvalue
##     model: Surv(survival,status) ~ Expression + Recurrence + Age + Subtype + MGMT_status
##     genes matched mouse->human via babelgene ortholog (uppercase fallback)
##
## DROPPED vs the manual xlsx (not reproducible from this pipeline):
##   - manual literature curation (Proliferation/Migration/Invasion/Stem.Cells/
##     PMID/Score/Overall.Summary/Interesting.papers)
##   - DrugList (manual Mayan drug join)
##   - nt_m_vs_f and JQ1 DE blocks (separate experiments, not in this repo's pipeline)
## =====================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(readr)
  library(survival); library(babelgene); library(biomaRt)
})

## ---- Paths ----------------------------------------------------------
g        <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci"
cc_dir   <- file.path(g, "04_callingcard_analysis", "output")
de_dir   <- file.path(g, "02_rna_analysis", "bulk_rna", "output_DE_Egr1KD_gRNA3_vs_Neg1")
out_dir  <- file.path(g, "04_callingcard_analysis", "output_full_summary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

glio_dir     <- "C:/Users/loril/Documents/multivariate_analysis/glioVis"
tcga_rds_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output"

lfc_cut <- 0.5; p_cut <- 0.05
all_covarID <- c("Recurrence", "Age", "Subtype", "MGMT_status")
OS_var <- "survival"; OS_event <- "status"

## ---- 1. CC∩DE gene sets (nearest gene within 20kb, intersected w/ DE) ----
read_cc <- function(f) unique(na.omit(read.delim(file.path(cc_dir, f))$Gene))
male_cc   <- read_cc("overlap_MaleCC_20kb_NearestGene_vs_Egr1g3vsNeg1.csv")
female_cc <- read_cc("overlap_FemaleCC_20kb_NearestGene_vs_Egr1g3vsNeg1.csv")
all_genes <- union(male_cc, female_cc)
cat("CC∩DE genes:  Male =", length(male_cc), " | Female =", length(female_cc),
    " | union =", length(all_genes), "\n")

## ---- 2. DE tables (both sexes) --------------------------------------
read_de <- function(f, tag) read.delim(file.path(de_dir, f)) %>%
  filter(!is.na(SYMBOL)) %>%
  transmute(SYMBOL,
            !!paste0(tag, "_log2FoldChange") := log2FoldChange,
            !!paste0(tag, "_pvalue") := pvalue,
            !!paste0(tag, "_padj")   := padj) %>%
  distinct(SYMBOL, .keep_all = TRUE)
de_m <- read_de("Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt",   "kd3_male")
de_f <- read_de("Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt", "kd3_female")

sig_call <- function(lfc, p) case_when(
  !is.na(lfc) & lfc >=  lfc_cut & p <= p_cut ~ "Up after KD",
  !is.na(lfc) & lfc <= -lfc_cut & p <= p_cut ~ "Down after KD",
  TRUE ~ "Not DE")

## ---- 3. biomaRt annotation: description, gene_biotype, GO (BP) -------
cat("Querying biomaRt (mouse) for description/biotype/GO ...\n")
mart <- useEnsembl("ensembl", "mmusculus_gene_ensembl")
anno <- getBM(attributes = c("external_gene_name", "description", "gene_biotype"),
              filters = "external_gene_name", values = all_genes, mart = mart) %>%
  distinct(external_gene_name, .keep_all = TRUE)
go <- getBM(attributes = c("external_gene_name", "name_1006", "namespace_1003"),
            filters = "external_gene_name", values = all_genes, mart = mart) %>%
  filter(namespace_1003 == "biological_process", name_1006 != "") %>%
  group_by(external_gene_name) %>%
  summarise(GO_Term_Description = paste(unique(name_1006), collapse = "; "), .groups = "drop")
annotation <- tibble(SYMBOL = all_genes) %>%
  left_join(anno, by = c("SYMBOL" = "external_gene_name")) %>%
  left_join(go,   by = c("SYMBOL" = "external_gene_name"))

## ---- 4. TCGA male/female pheno (RDS fast-path, else rebuild) ---------
fact <- function(df){ for(c in c("Recurrence","Subtype","MGMT_status")) df[[c]] <- factor(df[[c]]); df }
rds_m <- file.path(tcga_rds_dir, "tcga_pheno_male_allgenes_updated.rds")
rds_f <- file.path(tcga_rds_dir, "tcga_pheno_female_allgenes_updated.rds")
if (file.exists(rds_m) && file.exists(rds_f)) {
  tcga_m <- fact(readRDS(rds_m)); tcga_f <- fact(readRDS(rds_f))
} else {
  pheno <- read.table(file.path(glio_dir,"2024-06-04_TCGA_GBM_pheno.txt"), sep="\t", header=TRUE, check.names=FALSE)
  rownames(pheno) <- pheno$Sample
  ex <- read.table(file.path(glio_dir,"2024-06-04_TCGA_GBM_expression.txt"), sep="\t", header=TRUE, check.names=FALSE)
  rownames(ex) <- ex$Sample; ex <- as.matrix(ex[, setdiff(colnames(ex),"Sample")])
  common <- intersect(pheno$Sample, rownames(ex)); pheno <- pheno[common,]; ex <- ex[common,]
  pheno <- pheno[pheno$IDH1_status=="Wild-type" & !is.na(pheno$Gender),]; ex <- ex[pheno$Sample,]
  tp <- cbind(pheno, as.data.frame(ex))
  tcga_m <- fact(tp[tp$Gender=="Male",]); tcga_f <- fact(tp[tp$Gender=="Female",])
}
cat("TCGA:", nrow(tcga_m), "M |", nrow(tcga_f), "F patients\n")

## ---- 5. ortholog mouse->human gene matcher (uppercase fallback) ------
tcga_cols     <- colnames(tcga_m)
upper_to_orig <- setNames(tcga_cols, toupper(tcga_cols)); upper_to_orig <- upper_to_orig[!duplicated(names(upper_to_orig))]
match_gene <- function(sym) {
  o   <- tryCatch(orthologs(genes = sym, species = "mouse", human = FALSE), error = function(e) NULL)
  hum <- if (!is.null(o) && nrow(o)) toupper(o$human_symbol) else character(0)
  cand <- unique(c(hum, toupper(sym)))
  hit <- cand[cand %in% names(upper_to_orig)]
  if (length(hit)) unname(upper_to_orig[hit[1]]) else NA_character_
}
tcga_gene <- vapply(all_genes, match_gene, character(1))
cat("Genes matched to TCGA (ortholog):", sum(!is.na(tcga_gene)), "/", length(all_genes), "\n")

## ---- 6. per-gene multivariate Cox -> named HR/pvalue per term --------
lbl_term <- function(v, gene) {
  v <- gsub(gene, "Expression", v, fixed = TRUE)
  v <- gsub("MGMT_statusUnmethylated", "MGMT_Unmethylated",    v)
  v <- gsub("RecurrenceRecurrent",     "Recurrence_Recurrent", v)
  v <- gsub("RecurrenceSecondary",     "Recurrence_Secondary", v)
  v <- gsub("Recurrence2nd Recurrence","Recurrence_2nd",       v)
  v <- gsub("SubtypeMesenchymal",      "Subtype_Mesenchymal",  v)
  v <- gsub("SubtypeProneural",        "Subtype_Proneural",    v)
  v
}
cox_row <- function(pheno, gene, suffix) {
  if (is.na(gene) || !gene %in% colnames(pheno)) return(NULL)
  fm <- as.formula(paste0("Surv(",OS_var,",",OS_event,") ~ `", gene, "` + ", paste(all_covarID, collapse="+")))
  fit <- tryCatch(coxph(fm, data = pheno), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit); ci <- s$conf.int; co <- s$coefficients
  row <- setNames(list(as.numeric(s$logtest["pvalue"])), paste0("Overall_pvalue_", suffix))
  for (v in rownames(co)) {
    lab <- lbl_term(v, gene)
    row[[paste0(lab, "_HR_",     suffix)]] <- ci[v, "exp(coef)"]
    row[[paste0(lab, "_pvalue_", suffix)]] <- co[v, "Pr(>|z|)"]
  }
  as_tibble(row)
}
multivar <- function(genes_tcga) {
  rows <- lapply(seq_along(genes_tcga), function(i) {
    g <- genes_tcga[i]
    m <- cox_row(tcga_m, g, "M"); f <- cox_row(tcga_f, g, "F")
    if (is.null(m) && is.null(f)) return(NULL)
    bind_cols(tibble(.gene_key = g), if(!is.null(m)) m, if(!is.null(f)) f)
  })
  bind_rows(rows)
}

## ---- 7. assemble one sex's summary ----------------------------------
build <- function(symbols) {
  base <- tibble(SYMBOL = symbols) %>%
    mutate(in_male_CC   = SYMBOL %in% male_cc,
           in_female_CC = SYMBOL %in% female_cc) %>%
    left_join(de_m, by = "SYMBOL") %>%
    left_join(de_f, by = "SYMBOL") %>%
    mutate(kd3_male_sig   = sig_call(kd3_male_log2FoldChange,   kd3_male_pvalue),
           kd3_female_sig = sig_call(kd3_female_log2FoldChange, kd3_female_pvalue)) %>%
    relocate(kd3_male_sig,   .after = kd3_male_padj) %>%
    relocate(kd3_female_sig, .after = kd3_female_padj) %>%
    left_join(annotation, by = "SYMBOL") %>%
    mutate(.gene_key = tcga_gene[match(SYMBOL, all_genes)])
  mv <- multivar(unique(na.omit(base$.gene_key)))
  base %>% left_join(mv, by = ".gene_key") %>% dplyr::select(-.gene_key)
}

male_summary   <- build(male_cc)
female_summary <- build(female_cc)

mp <- file.path(out_dir, "Egr1CC_DE_full_summary_MALES.csv")
fp <- file.path(out_dir, "Egr1CC_DE_full_summary_FEMALES.csv")
write_csv(male_summary,   mp)
write_csv(female_summary, fp)

cat("\nDone.\n  Male  :", nrow(male_summary),  "genes ,", ncol(male_summary),  "cols ->", mp, "\n")
cat("  Female:", nrow(female_summary),"genes ,", ncol(female_summary),"cols ->", fp, "\n")
