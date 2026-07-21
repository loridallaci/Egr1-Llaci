# =============================================================================
# RENIN motif enrichment across cortex development, restricted to the PAPER'S
# cells (authors' deposited CELLxGENE QC set + sex), run in the ORIGINAL RENIN
# environment (R 4.2.2 / Seurat 4.1.1 / Signac 1.10 / RENIN 1.0).
#
# Difference from cortex_dev_RENIN_pipeline.R: before running RENIN, each
# withMotifs_V4 object is SUBSET to the paper's cells via barcode reconciliation
# and its $sex is set from the paper's authoritative donor annotation.
#
#   Local barcode  : {16bp_core}-{1|2}   (suffix = donor within stage)
#   Paper barcode  : {prefix}_{16bp_core}-1  (prefix = donor; sex in metadata)
#   Reconcile      : match 16bp core within the stage's two donors.
#
# Run in the original env:
#   LD_LIBRARY_PATH=$HOME/local/lib \
#   /home/lllaci/R-4.2.2/bin/Rscript cortex_dev_RENIN_paperCells.R
# =============================================================================
suppressMessages({
  library(Seurat); library(Signac); library(SeuratWrappers); library(RENIN)
  library(harmony); library(GenomeInfoDb); library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2); library(ggrepel); library(dplyr); library(patchwork)
})

data_dir    <- "/home/lllaci/data/cortex_development"
output_base <- "/home/lllaci/renin_paperCells_output"
annot_csv   <- "/home/lllaci/paper_cell_annotation.csv"     # scp'd here
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

stage_files <- list(
  LateFetal = "LaFet_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
  Infant    = "Inf_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
  Child     = "Child_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
  Adol      = "Adol_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds",
  Adult     = "Adult_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds")
col_titles <- c(LateFetal="Late Fetal", Infant="Infant", Child="Child",
                Adol="Adolescence", Adult="Adult")
# suffix (-1 / -2) -> paper donor, per stage (verified by core-barcode overlap)
stage_donor <- list(
  LateFetal = c("1"="LaFet1","2"="LaFet2"), Infant = c("1"="Inf1","2"="Inf2"),
  Child = c("1"="Child1","2"="Child2"),     Adol   = c("1"="Adol1","2"="Adol2"),
  Adult = c("1"="Adult1","2"="Adult2"))

paper <- read.csv(annot_csv, stringsAsFactors = FALSE)   # donor, core, sex, celltype
paper$key <- paste(paper$donor, paper$core, sep = ":")

# ---- subset a withMotifs object to the paper's cells + set authoritative sex ----
subset_to_paper <- function(obj, stage) {
  bc   <- colnames(obj)
  core <- sub("-[0-9]+$", "", bc)
  suf  <- sub(".*-", "", bc)
  donor <- stage_donor[[stage]][suf]
  key  <- paste(donor, core, sep = ":")
  m    <- match(key, paper$key)
  keep <- !is.na(m)
  cat(sprintf("  %s: local %d cells, paper-matched %d (%.0f%%)\n",
              stage, length(bc), sum(keep), 100*mean(keep)))
  obj <- obj[, keep]
  obj$sex      <- paper$sex[m[keep]]
  obj$celltype <- paper$celltype[m[keep]]
  obj
}

# ----- RENIN volcano plotter (verbatim from cortex_dev_RENIN_pipeline.R) -----
analyze_and_plot_motifs <- function(seurat, cres, output_dir, stage_name, label,
                                    num_top_motifs = 20) {
  motifs <- FindMotifs(object = seurat, features = cres)
  motifs$logp <- -log10(motifs$pvalue)
  motifs <- motifs[order(motifs$logp, decreasing = TRUE), ]
  write.csv(motifs, file.path(output_dir, paste0(label, "_all_motifs.csv")), row.names = TRUE)
  motifs$significant <- motifs$p.adjust <= 0.05
  motifs$color  <- ifelse(motifs$significant, "green3", "grey80")
  motifs$border <- ifelse(motifs$significant, 0.4, 0.1)
  motifs$label  <- ifelse(rank(-motifs$logp) <= num_top_motifs, motifs$motif.name, "")
  max_y <- max(motifs$logp, na.rm = TRUE) * 1.05
  max_x <- max(motifs$fold.enrichment, na.rm = TRUE) * 1.05
  g <- ggplot(motifs, aes(x = fold.enrichment, y = logp, label = label)) +
    geom_point(aes(fill = color), color = "black", stroke = motifs$border, pch = 21, size = 2) +
    scale_fill_identity() +
    geom_text_repel(max.overlaps = 500, size = 3, point.padding = 0.5, force = 2,
                    box.padding = 0.4, min.segment.length = 0, segment.color = "grey40",
                    segment.size = 0.4, segment.alpha = 0.8) +
    theme_classic() + ylab("-log10(p-value)") + xlab("Fold enrichment") +
    ylim(c(0, max_y)) + xlim(c(0, max_x)) + ggtitle(paste0("Motif Enrichment: ", label)) +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.title = element_text(size = 16), axis.text = element_text(size = 14))
  pdf(file.path(output_dir, paste0("RENIN_", stage_name, "_MotifEnrichment_SignificantLabeled_", label, ".pdf")))
  print(g); dev.off()
  g
}

# ----- per-stage RENIN (verbatim logic from cortex_dev_RENIN_pipeline.R) -----
process_all_cells <- function(seurat, stage_name, output_dir,
                              cells_per_partition = 100, num_top_degs = 100, num_top_motifs = 20) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat(paste0("\n========== ", stage_name, " ==========\n"))
  peak_gr <- granges(seurat[["peaks"]]);  anno_gr <- Annotation(seurat[["peaks"]])
  seqlevelsStyle(peak_gr) <- "UCSC";      seqlevelsStyle(anno_gr) <- "UCSC"
  seqlevels(anno_gr) <- paste0("chr", seqlevels(anno_gr))
  common <- intersect(seqlevels(peak_gr), seqlevels(anno_gr))
  peak_gr <- keepSeqlevels(peak_gr, common, pruning.mode = "coarse")
  anno_gr <- keepSeqlevels(anno_gr, common, pruning.mode = "coarse")
  seurat[["peaks"]]@ranges <- peak_gr;    Annotation(seurat[["peaks"]]) <- anno_gr
  DefaultAssay(seurat) <- "RNA"
  seurat <- NormalizeData(seurat); seurat <- FindVariableFeatures(seurat)
  seurat <- ScaleData(seurat);     seurat <- RunPCA(seurat)
  DefaultAssay(seurat) <- "peaks"
  seurat[["SCT"]] <- seurat[["RNA"]]; seurat[["RNA"]] <- NULL
  seurat <- RunHarmony(seurat, group.by.vars = "sex", reduction = "pca",
                       assay.use = "SCT", project.dim = FALSE, reduction.save = "harmony_SCT")
  set.seed(1234)
  seurat <- RunTFIDF(seurat, assay = "peaks")
  seurat <- FindTopFeatures(seurat, min.cutoff = "q0", assay = "peaks")
  seurat <- RunSVD(seurat, assay = "peaks")
  seurat <- RunHarmony(seurat, group.by.vars = "sex", reduction = "lsi",
                       assay.use = "peaks", project.dim = FALSE, reduction.save = "harmony_peaks")
  Idents(seurat) <- seurat$sex
  mpt <- subset(seurat, sex %in% c("female", "male"))
  de.genes <- prepare_degs(mpt, ident.1 = "female", ident.2 = "male")
  write.csv(de.genes, file.path(output_dir, paste0(stage_name, "_AllCells_harmonySCT_FvsM_DEG.csv")))
  mats <- prepare_pseudocell_matrix(seurat, assay = c("peaks", "SCT"),
                                    cells_per_partition = cells_per_partition,
                                    reduction1 = "harmony_peaks", reduction2 = "harmony_SCT")
  expr_mat <- mats[["SCT"]]; peak_mat <- mats[["peaks"]]
  gene_list <- rownames(de.genes)[1:min(num_top_degs, nrow(de.genes))]
  peak_results <- run_peak_aen(seurat, expr_mat, peak_mat, gene_list,
                               lambda2 = 0.5, max_distance = 5e+05, num_bootstraps = 100)
  aen_lists <- make_aen_lists(peak_results)
  fr_genes <- rownames(de.genes)[which(de.genes$avg_log2FC < 0)]
  cre_scores <- lapply(peak_results, function(x)
    x[[4]][union(1, which(x[[4]][, "coef_if_kept"] != 0)), "coef_if_kept"] *
      ifelse(x[[1]] %in% fr_genes, -1, 1))
  cre_scores <- cre_scores[which(lengths(cre_scores) > 1)]
  cre_total  <- bind_rows(cre_scores); cre_total[is.na(cre_total)] <- 0
  cre_total  <- cre_total[, -1]; sums <- colSums(cre_total)
  fr_cres <- names(sums)[sums < 0]; h_cres <- names(sums)[sums > 0]
  seurat <- RegionStats(seurat, genome = BSgenome.Hsapiens.UCSC.hg38)
  g_male   <- analyze_and_plot_motifs(seurat, fr_cres, output_dir, stage_name, "Male_AllCells",   num_top_motifs)
  g_female <- analyze_and_plot_motifs(seurat, h_cres,  output_dir, stage_name, "Female_AllCells", num_top_motifs)
  list(male = g_male, female = g_female)
}

# ============================= DRIVER =============================
panels <- list()
for (stage in names(stage_files)) {
  cat("\n#### loading", stage, "####\n")
  obj <- readRDS(file.path(data_dir, stage_files[[stage]]))
  obj <- subset_to_paper(obj, stage)
  panels[[stage]] <- process_all_cells(obj, stage_name = stage,
                       output_dir = file.path(output_base, stage),
                       cells_per_partition = 100, num_top_degs = 100, num_top_motifs = 20)
  rm(obj); gc()
}
male_row   <- lapply(names(stage_files), function(s) panels[[s]]$male   + ggtitle(col_titles[[s]]))
female_row <- lapply(names(stage_files), function(s) panels[[s]]$female + ggtitle(NULL))
montage <- (wrap_plots(male_row, nrow = 1)) / (wrap_plots(female_row, nrow = 1)) +
  plot_annotation(title = "RENIN motif enrichment across cortex development (paper's cells)")
ggsave(file.path(output_base, "SupplFig1_RENIN_montage_paperCells.pdf"), montage, width = 25, height = 10, limitsize = FALSE)
ggsave(file.path(output_base, "SupplFig1_RENIN_montage_paperCells.png"), montage, width = 25, height = 10, dpi = 300, limitsize = FALSE)
cat("\nDone. Output ->", output_base, "\n")
