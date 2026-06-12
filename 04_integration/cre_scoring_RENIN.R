# =============================================================================
# 04_integration_analysis: CRE Scoring using RENIN
# =============================================================================
# Description:
#   Identifies cis-regulatory elements (CREs) associated with sex-biased
#   gene expression using the RENIN package. Computes pseudocell matrices
#   from harmonized RNA and ATAC embeddings, runs peak-gene association
#   via elastic net regression, and scores CREs as female-enriched
#   (female_cres) or male-enriched (male_cres).
#
# Run order (within 03_atac_analysis):
#   1. cre_scoring_RENIN.R         <-- this script
#   2. figure_motif_enrichment.R
#
# Input:
#   female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925
#   _withPeaks_chromVARadded_111425_V4.rds
#
# Output:
#   - cre_lists_female_male_RENIN.rds  (female_cres and male_cres peak lists)
#   - obj_RENIN_processed.rds          (object with RegionStats added)
# =============================================================================

# --- Libraries ----------------------------------------------------------------

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(RENIN))
suppressMessages(library(harmony))
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)

# --- Load object --------------------------------------------------------------

lot6 <- readRDS("data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds")  # original: "/home/lllaci/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds"

# -----------------------------
# Setup
# -----------------------------
obj <- lot6
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# Rename RNA to SCT for compatibility with RENIN pseudocell pipeline
DefaultAssay(obj) <- 'peaks' # change default assay as you can't delete the default assay
obj[["SCT"]] <- obj[["RNA"]]
obj[["RNA"]] <- NULL
# Add sex metadata
obj$sex <- ifelse(
  grepl("-1$", colnames(obj)), "female",
  ifelse(grepl("-2$", colnames(obj)), "male", NA)
)

obj <- RunHarmony(
  object = obj,
  group.by.vars = 'sex',
  reduction = 'pca',
  assay.use = 'SCT',
  project.dim = FALSE,
  reduction.save = 'harmony_SCT'
)


set.seed(1234)
DefaultAssay(obj) <- "peaks"

## RunTFIDF
obj <- RunTFIDF(obj, assay = "peaks")

## FindTopFeatures
obj <- FindTopFeatures(obj, min.cutoff = 'q0', assay = "peaks")

## RunSVD
obj <- RunSVD(obj, assay = "peaks")

obj <- RunHarmony(
  object = obj,
  group.by.vars = 'sex',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE,
  reduction.save = 'harmony_peaks'
)

# -----------------------------
# Build pseudocell matrices
# -----------------------------
Sys.time()

mats <- prepare_pseudocell_matrix(
  obj,
  assay = c("peaks", "SCT"),
  cells_per_partition = 100,
  reduction1 = "harmony_peaks",
  reduction2 = "harmony_SCT"
)

Sys.time()

expr_mat <- mats[["SCT"]]
peak_mat <- mats[["peaks"]]

# -----------------------------
# Define identity
# -----------------------------
Idents(obj) <- obj$sex

# -----------------------------
# Load external DE table
# -----------------------------
de_df <- read.csv("data/DE_male_vs_female_allcells_allgenes.csv")  # original: "/home/lllaci/DE_male_vs_female_allcells_allgenes.csv"
#de_df <- read.csv("data/DE_male_vs_female_allcells_allgenes_log025.csv")  # original: "/home/lllaci/DE_male_vs_female_allcells_allgenes_log025.csv"
#de_df <- read.csv("data/DE_male_vs_female_allcells_allgenes_nologfcthreshhold.csv")  # original: "/home/lllaci/DE_male_vs_female_allcells_allgenes_nologfcthreshhold.csv"
gene_list <- de_df$gene[1:100]


# -----------------------------
# Peak-gene association analysis
# -----------------------------
peak_results <- run_peak_aen(
  obj,
  expr_mat,
  peak_mat,
  gene_list,
  lambda2 = 0.5,
  max_distance = 5e5,
  num_bootstraps = 100
)

aen_lists <- make_aen_lists(peak_results)

# -----------------------------
# Define female-enriched genes
# (negative logFC = higher in female)
# -----------------------------
female_genes <- de_df$X[which(de_df$avg_log2FC < 0)]

# -----------------------------
# Compute CRE scores
# -----------------------------
cre_scores <- lapply(peak_results, function(x) {
  x[[4]][
    union(1, which(x[[4]][, "coef_if_kept"] != 0)),
    "coef_if_kept"
  ] *
    ifelse(x[[1]] %in% female_genes, -1, 1)
})

# Filter valid entries
cre_scores <- cre_scores[which(unlist(lapply(cre_scores, length)) > 1)]

# Combine into matrix
cre_total_scores <- bind_rows(cre_scores)

# Replace NA with 0
cre_total_scores[is.na(cre_total_scores)] <- 0

# Remove first column (gene IDs)
cre_total_scores <- cre_total_scores[, -1]

# Column-wise aggregation
cre_total_scores_sums <- colSums(cre_total_scores)

# Count number of genes per CRE
cre_num_genes <- apply(cre_total_scores, 2, function(x) {
  length(which(x != 0))
})

# -----------------------------
# Split CREs by enrichment
# -----------------------------
female_cres <- names(cre_total_scores_sums)[cre_total_scores_sums < 0]
male_cres   <- names(cre_total_scores_sums)[cre_total_scores_sums > 0]

cat("Female-enriched CREs:", length(female_cres), "\n")
cat("Male-enriched CREs:  ", length(male_cres), "\n")

# -----------------------------
# Genome annotation (optional)
# -----------------------------
library(BSgenome.Mmusculus.UCSC.mm10)

obj <- RegionStats(obj, genome = BSgenome.Mmusculus.UCSC.mm10)


# --- Save ---------------------------------------------------------------------

saveRDS(
  list(female_cres = female_cres, male_cres = male_cres),
  "output/cre_lists_female_male_RENIN.rds"  # original: "/home/lllaci/data/cre_lists_female_male_RENIN.rds"
)
saveRDS(obj, "output/obj_RENIN_processed.rds")  # original: "/home/lllaci/data/obj_RENIN_processed.rds"

message("Done. CRE lists and processed object saved.")

