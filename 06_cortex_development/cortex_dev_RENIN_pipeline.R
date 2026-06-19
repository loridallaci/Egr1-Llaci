# =============================================================================
# Supplementary Figure 1 — RENIN motif enrichment across cortex development
# Consolidated reproduction script (Males x Females  ×  5 developmental stages)
#
# Author: L. Llaci
# Gathered from:
#   - Roussos_cortex_development_RENIN.Rmd              (preprocessing + AddMotifs)
#   - automated_pipeline_allCells_RENINonly_010626.sh  (RENIN run + volcano plots)
#   - RENIN on cortex development Male vs Female.pptx   (manual 2x5 montage)
#
# Figure = 10 volcano panels: Fold enrichment (x) vs -log10(p-value) (y),
#   green = significant (p.adjust <= 0.05), top-20 motifs labeled (ggrepel).
#   Male panels  = male-enriched CREs (fr_cres, summed CRE score < 0)
#   Female panels = female-enriched CREs (h_cres, summed CRE score > 0)
#
# Run on the HPC via srun/sbatch (NOT the login node) — RENIN is compute-heavy.
# =============================================================================

suppressMessages({
  library(Seurat)
  library(Signac)
  library(SeuratWrappers)
  library(RENIN)
  library(harmony)
  library(GenomeInfoDb)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# STAGE 0 (one-time, done in the .Rmd) — build the *_withMotifs_V4 objects
# -----------------------------------------------------------------------------
# For each stage's aggregated object you:
#   1. assign sex from the barcode suffix (-1 = female, -2 = male)
#   2. run AddMotifs() with JASPAR2020 CORE vertebrates
#   3. saveRDS() as *_aggregated_object_..._withMotifs_V4_122825.rds
# Re-run only if the *_withMotifs_V4_122825.rds files do not already exist.
#
# library(JASPAR2020); library(TFBSTools)
# pfm <- getMatrixSet(JASPAR2020,
#          opts = list(collection = "CORE", tax_group = "vertebrates",
#                      all_versions = FALSE))
# object <- readRDS("<stage>_aggregated_object_Roussos_112223_...Available.rds")
# object$sex <- ifelse(grepl("-1$", colnames(object)), "female",
#                ifelse(grepl("-2$", colnames(object)), "male", NA))
# DefaultAssay(object) <- "peaks"
# object <- AddMotifs(object, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
# saveRDS(object, "<stage>_aggregated_object_..._withMotifs_V4_122825.rds")

# -----------------------------------------------------------------------------
# Volcano plotting function (identical to the pipeline's analyze_and_plot_motifs)
# -----------------------------------------------------------------------------
analyze_and_plot_motifs <- function(seurat, cres, output_dir, stage_name,
                                    label, num_top_motifs = 20) {

  motifs <- FindMotifs(object = seurat, features = cres)
  motifs$logp <- -log10(motifs$pvalue)
  motifs <- motifs[order(motifs$logp, decreasing = TRUE), ]

  # Save the full motif table (source data for the panel)
  write.csv(motifs,
            file.path(output_dir, paste0(label, "_all_motifs.csv")),
            row.names = TRUE)

  # --- plot aesthetics ---
  motifs$significant <- motifs$p.adjust <= 0.05
  motifs$color  <- ifelse(motifs$significant, "green3", "grey80")
  motifs$border <- ifelse(motifs$significant, 0.4, 0.1)
  motifs$label  <- ifelse(rank(-motifs$logp) <= num_top_motifs,
                          motifs$motif.name, "")

  max_y <- max(motifs$logp, na.rm = TRUE) * 1.05
  max_x <- max(motifs$fold.enrichment, na.rm = TRUE) * 1.05

  g <- ggplot(motifs, aes(x = fold.enrichment, y = logp, label = label)) +
    geom_point(aes(fill = color), color = "black",
               stroke = motifs$border, pch = 21, size = 2) +
    scale_fill_identity() +
    geom_text_repel(max.overlaps = 500, size = 3, point.padding = 0.5,
                    force = 2, box.padding = 0.4, min.segment.length = 0,
                    segment.color = "grey40", segment.size = 0.4,
                    segment.alpha = 0.8) +
    theme_classic() +
    ylab("-log10(p-value)") + xlab("Fold enrichment") +
    ylim(c(0, max_y)) + xlim(c(0, max_x)) +
    ggtitle(paste0("Motif Enrichment: ", label)) +
    # Presentation-legible fonts (>=14 ticks, >=16 axis titles, >=18 title)
    theme(plot.title  = element_text(size = 18, hjust = 0.5),
          axis.title  = element_text(size = 16),
          axis.text   = element_text(size = 14))

  pdf(file.path(output_dir,
        paste0("RENIN_", stage_name,
               "_MotifEnrichment_SignificantLabeled_", label, ".pdf")))
  print(g)
  dev.off()

  return(g)   # return the ggplot so we can assemble the montage
}

# -----------------------------------------------------------------------------
# Full per-stage RENIN run -> returns the two ggplot panels (male, female)
# -----------------------------------------------------------------------------
process_all_cells <- function(seurat, stage_name, output_dir,
                              cells_per_partition = 100,
                              num_top_degs = 100, num_top_motifs = 20) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat(paste0("\n========== ", stage_name, " ==========\n"))

  # 1. Match chromosome naming between peaks and annotation (UCSC "chr*")
  peak_gr <- granges(seurat[["peaks"]]);  anno_gr <- Annotation(seurat[["peaks"]])
  seqlevelsStyle(peak_gr) <- "UCSC";      seqlevelsStyle(anno_gr) <- "UCSC"
  seqlevels(anno_gr) <- paste0("chr", seqlevels(anno_gr))
  common <- intersect(seqlevels(peak_gr), seqlevels(anno_gr))
  peak_gr <- keepSeqlevels(peak_gr, common, pruning.mode = "coarse")
  anno_gr <- keepSeqlevels(anno_gr, common, pruning.mode = "coarse")
  seurat[["peaks"]]@ranges <- peak_gr;    Annotation(seurat[["peaks"]]) <- anno_gr

  # 2. RNA -> Harmony (batch = sex), move RNA into "SCT" slot for RENIN
  DefaultAssay(seurat) <- "RNA"
  seurat <- NormalizeData(seurat); seurat <- FindVariableFeatures(seurat)
  seurat <- ScaleData(seurat);     seurat <- RunPCA(seurat)
  DefaultAssay(seurat) <- "peaks"
  seurat[["SCT"]] <- seurat[["RNA"]]; seurat[["RNA"]] <- NULL
  seurat <- RunHarmony(seurat, group.by.vars = "sex", reduction = "pca",
                       assay.use = "SCT", project.dim = FALSE,
                       reduction.save = "harmony_SCT")

  # 3. Peaks -> LSI -> Harmony
  set.seed(1234)
  seurat <- RunTFIDF(seurat, assay = "peaks")
  seurat <- FindTopFeatures(seurat, min.cutoff = "q0", assay = "peaks")
  seurat <- RunSVD(seurat, assay = "peaks")
  seurat <- RunHarmony(seurat, group.by.vars = "sex", reduction = "lsi",
                       assay.use = "peaks", project.dim = FALSE,
                       reduction.save = "harmony_peaks")

  # 4. Differential expression: Female vs Male
  Idents(seurat) <- seurat$sex
  mpt <- subset(seurat, sex %in% c("female", "male"))
  de.genes <- prepare_degs(mpt, ident.1 = "female", ident.2 = "male")
  write.csv(de.genes,
            file.path(output_dir, paste0(stage_name, "_AllCells_harmonySCT_FvsM_DEG.csv")))

  # 5. Pseudocells -> peak AEN (RENIN core)
  mats <- prepare_pseudocell_matrix(seurat, assay = c("peaks", "SCT"),
                                    cells_per_partition = cells_per_partition,
                                    reduction1 = "harmony_peaks",
                                    reduction2 = "harmony_SCT")
  expr_mat <- mats[["SCT"]]; peak_mat <- mats[["peaks"]]
  gene_list <- rownames(de.genes)[1:min(num_top_degs, nrow(de.genes))]
  peak_results <- run_peak_aen(seurat, expr_mat, peak_mat, gene_list,
                               lambda2 = 0.5, max_distance = 5e+05,
                               num_bootstraps = 100)
  aen_lists <- make_aen_lists(peak_results)

  # 6. Sex-directional CRE scores: <0 = male-enriched, >0 = female-enriched
  fr_genes <- rownames(de.genes)[which(de.genes$avg_log2FC < 0)]
  cre_scores <- lapply(peak_results, function(x)
    x[[4]][union(1, which(x[[4]][, "coef_if_kept"] != 0)), "coef_if_kept"] *
      ifelse(x[[1]] %in% fr_genes, -1, 1))
  cre_scores <- cre_scores[which(lengths(cre_scores) > 1)]
  cre_total  <- bind_rows(cre_scores); cre_total[is.na(cre_total)] <- 0
  cre_total  <- cre_total[, -1]; sums <- colSums(cre_total)
  fr_cres <- names(sums)[sums < 0]   # male-enriched
  h_cres  <- names(sums)[sums > 0]   # female-enriched

  # 7. Motif enrichment + volcano panels
  seurat <- RegionStats(seurat, genome = BSgenome.Hsapiens.UCSC.hg38)
  g_male   <- analyze_and_plot_motifs(seurat, fr_cres, output_dir, stage_name,
                                      "Male_AllCells",   num_top_motifs)
  g_female <- analyze_and_plot_motifs(seurat, h_cres,  output_dir, stage_name,
                                      "Female_AllCells", num_top_motifs)

  list(male = g_male, female = g_female)
}

# =============================================================================
# DRIVER — run all 5 stages and assemble the 2x5 montage
# =============================================================================
# Edit this base to wherever the *_withMotifs_V4_122825.rds objects live.
data_dir    <- "/home/lllaci/data/cortex_development"
output_base <- file.path(data_dir, "AllCells_Analysis3")

# Stage order matches the figure columns (Late Fetal -> Adult)
stage_files <- list(
  LateFetal = "LaFet_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
  Infant    = "Inf_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
  Child     = "Child_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
  Adol      = "Adol_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
  Adult     = "Adult_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds"
)
col_titles <- c(LateFetal = "Late Fetal", Infant = "Infant", Child = "Child",
                Adol = "Adolescence", Adult = "Adult")

panels <- list()
for (stage in names(stage_files)) {
  obj <- readRDS(file.path(data_dir, stage_files[[stage]]))
  panels[[stage]] <- process_all_cells(
    seurat = obj, stage_name = stage,
    output_dir = file.path(output_base, stage),
    cells_per_partition = 100, num_top_degs = 100, num_top_motifs = 20)
}

# --- Build montage: row 1 = Males, row 2 = Females, columns = stages ---
male_row   <- lapply(names(stage_files), function(s)
                panels[[s]]$male   + ggtitle(col_titles[[s]]))
female_row <- lapply(names(stage_files), function(s)
                panels[[s]]$female + ggtitle(NULL))

montage <- (wrap_plots(male_row,   nrow = 1)) /
           (wrap_plots(female_row, nrow = 1)) +
           plot_annotation(title = "Supplementary Figure 1 — RENIN motif enrichment across cortex development")

ggsave(file.path(output_base, "SupplementaryFigure1_RENIN_montage.pdf"),
       montage, width = 25, height = 10, limitsize = FALSE)
ggsave(file.path(output_base, "SupplementaryFigure1_RENIN_montage.png"),
       montage, width = 25, height = 10, dpi = 300, limitsize = FALSE)

cat("\nDone. Per-panel PDFs + assembled montage written to:\n  ", output_base, "\n")
