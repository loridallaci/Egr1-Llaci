## =====================================================================
## Consolidated: make ONLY the two gRNA3-vs-Neg1 DE files
##   - Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt
##   - Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt
## Extracted verbatim (logic) from:
##   Egr1KD_again_ButCombineEgr1KDsamples_only0filteredout_removeOutliers.Rmd
## Keeps Egr1 guides SEPARATE (combine step NOT applied), drops Neg gRNA2,
## removes the same outlier samples, DESeq2 design ~0+group.
## =====================================================================

library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)

## ---- 1. Paths -------------------------------------------------------
## Raw inputs read from the copy committed to git (02_rna_analysis/data_raw/)
dataDirectory      <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/data_raw"
metadataDirectory  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/data_raw"
## Outputs -> bulk_rna stage folder (git)
deOutputDirectory  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_Egr1KD_gRNA3_vs_Neg1"
dir.create(deOutputDirectory, showWarnings = FALSE, recursive = TRUE)

figureOutputDirectory <- file.path(deOutputDirectory, "qc")
dir.create(figureOutputDirectory, showWarnings = FALSE, recursive = TRUE)

## ---- QC plot helper: sample-distance heatmap + PCA on vst ------------
## Called BEFORE and AFTER outlier/Neg2 removal so the removals are visualized.
group_palette <- c(
  Egr1KD_gRNA2_Female      = "#6A0DAD", Egr1KD_gRNA3_Female      = "#B266FF",
  Egr1KD_gRNA2_Male        = "#FF8C00", Egr1KD_gRNA3_Male        = "#FFB84D",
  NoTreatment_gRNA1_Male   = "#1E90FF", NoTreatment_gRNA1_Female = "#FF69B4",
  NoTreatment_gRNA2_Male   = "#00CED1", NoTreatment_gRNA2_Female = "#C71585")
sex_palette <- c(Female = "#F08080", Male = "darkgreen")

run_qc_plots <- function(cd, sd, tag) {
  sd <- sd[colnames(cd), , drop = FALSE]
  d  <- DESeqDataSetFromMatrix(round(as.matrix(cd)), colData = sd, design = ~1)
  d  <- d[rowSums(counts(d)) > 1, ]
  v  <- vst(d, blind = TRUE)

  ## sample-sample Euclidean-distance heatmap
  sdist <- dist(t(assay(v)), method = "euclidean")
  ann   <- as.data.frame(colData(v)[, c("group", "Sex")])
  ann$group <- droplevels(factor(ann$group)); ann$Sex <- droplevels(factor(ann$Sex))
  pdf(file.path(figureOutputDirectory, paste0("SampleSample_vst_Heatmap_", tag, ".pdf")),
      width = 10, height = 8)
  pheatmap(as.matrix(sdist),
           col = colorRampPalette(c("blue", "white", "red"))(20),
           border_color = "black", fontsize = 10, fontsize_col = 7,
           show_rownames = FALSE, show_colnames = TRUE, annotation_col = ann,
           annotation_colors = list(Sex = sex_palette, group = group_palette),
           clustering_distance_rows = sdist, clustering_distance_cols = sdist,
           main = paste0("Euclidean distance (vst) - ", tag))
  invisible(dev.off())

  ## PCA (top 2000 variable genes): all + per sex
  ## NOTE: use intgroup="group" only. plotPCA(intgroup=c("group","Sex")) returns a
  ## `group` column that is the "group : Sex" interaction, which won't match
  ## group_palette -> all points render grey. Derive Sex from the group name instead.
  pd <- plotPCA(v, intgroup = "group", ntop = 2000, returnData = TRUE)
  pd$group <- factor(pd$group, levels = names(group_palette))
  pd$Sex   <- ifelse(grepl("Female", pd$group), "Female", "Male")
  pv <- round(100 * attr(pd, "percentVar"))
  gg <- function(df, ttl)
    ggplot(df, aes(PC1, PC2, fill = group)) +
      geom_point(shape = 21, size = 4, color = "black") +
      scale_fill_manual(values = group_palette) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
      xlab(paste0("PC1: ", pv[1], "%")) + ylab(paste0("PC2: ", pv[2], "%")) +
      ggtitle(ttl) + coord_fixed() + theme_bw()
  pdf(file.path(figureOutputDirectory, paste0("PCA_vst_", tag, ".pdf")), width = 7, height = 7)
  print(gg(pd, paste0(tag, " - all samples")))
  print(gg(pd[pd$Sex == "Male", ],   paste0(tag, " - Males")))
  print(gg(pd[pd$Sex == "Female", ], paste0(tag, " - Females")))
  invisible(dev.off())
}

countDataFile      <- file.path(dataDirectory, "Dedup_Counts.txt")
geneDataFile       <- file.path(dataDirectory, "GeneInfo_Egr1KD.csv")
sampleMetadataFile <- file.path(metadataDirectory, "Meta_data_Egr1KD.csv")

## ---- 2. Load count, gene, and sample data ---------------------------
countData <- read.delim(countDataFile, header = TRUE,
                        stringsAsFactors = FALSE, row.names = 1)

geneData <- read.delim(geneDataFile, header = TRUE, stringsAsFactors = FALSE,
                       row.names = 1, sep = ',')
geneData <- geneData[!duplicated(geneData$ensembl_gene_id), ]

sampleData <- read.delim(sampleMetadataFile, header = TRUE, row.names = 1, sep = ',')
sampleData <- as.data.frame(t(sampleData))
colnames(sampleData)[3] <- 'group'          # 3rd col = Treatment_Sex -> "group"

countData <- countData[, rownames(sampleData)]

## QC BEFORE removal (all samples). Wrapped so a locked/open PDF doesn't halt the DE.
tryCatch(run_qc_plots(countData, sampleData, "BEFORE_removal"),
         error = function(e) message("QC BEFORE_removal skipped: ", conditionMessage(e)))

## ---- 3. Remove outliers + Neg gRNA2; KEEP Egr1KD_gRNA2 in the model --
## Egr1KD_gRNA2 stays in the DESeq2 model; the gRNA3-vs-Neg1 contrast is
## extracted from that model (original design -> 579 Male / 720 Female).
samples_to_remove <- c(
  "F6_Egr1_gRNA2_Rep2",
  "F6_Egr1_gRNA2_Rep3",
  "F6_Egr1_gRNA2_Rep4"
)

sampleData <- sampleData %>%
  filter(!(rownames(.) %in% samples_to_remove)) %>%   # drop outliers
  filter(Treatment != "NoTreatment_gRNA2")            # drop Neg gRNA2 (already absent from raw data)

countData <- countData[, colnames(countData) %in% rownames(sampleData)]
geneData  <- geneData[geneData$ensembl_gene_id %in% rownames(countData), ]

sampleData[, c('Treatment','Sex','group')] <-
  lapply(sampleData[, c('Treatment','Sex','group')], factor)

## ---- 4. Build DESeq2 object -----------------------------------------
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData   = sampleData,
                              design    = ~0 + group)

## attach gene annotation (for SYMBOL column)
rownames(geneData) <- geneData$ensembl_gene_id
geneData <- geneData[rownames(dds), ]
mcols(dds) <- DataFrame(mcols(dds), geneData)

## minimal prefilter (only genes with all-zero counts removed) + size factors
dds2 <- dds[rowSums(counts(dds)) > 1, ]
dds2 <- estimateSizeFactors(dds2)

## ---- 4b. QC plots AFTER removal (clean set used for DE) -------------
tryCatch(run_qc_plots(countData, sampleData, "AFTER_removal"),
         error = function(e) message("QC AFTER_removal skipped: ", conditionMessage(e)))

## ---- 5. Differential expression -------------------------------------
dds2 <- DESeq(dds2)

## Male: Egr1 gRNA3 KD vs Neg gRNA1  (positive log2FC = UP in KD)
res <- results(dds2, contrast = c("group", "Egr1KD_gRNA3_Male", "NoTreatment_gRNA1_Male"))
res$SYMBOL <- mcols(dds2)[, colnames(mcols(dds2)) %in% "external_gene_name"]
write.table(res,
  file = file.path(deOutputDirectory, "Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt"),
  quote = FALSE, sep = "\t", col.names = NA)

## Female: Egr1 gRNA3 KD vs Neg gRNA1
res <- results(dds2, contrast = c("group", "Egr1KD_gRNA3_Female", "NoTreatment_gRNA1_Female"))
res$SYMBOL <- mcols(dds2)[, colnames(mcols(dds2)) %in% "external_gene_name"]
write.table(res,
  file = file.path(deOutputDirectory, "Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt"),
  quote = FALSE, sep = "\t", col.names = NA)

message("Done. Wrote the two gRNA3-vs-Neg1 DE tables to:\n", deOutputDirectory)
