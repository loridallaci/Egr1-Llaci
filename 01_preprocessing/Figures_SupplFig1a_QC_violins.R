## =====================================================================
## Supplementary Fig 1a - multiome QC violins (by sex), as PDFs into git.
##
## Fast render: loads the already-processed lot6 filtered+peaks Seurat object
## (QC metrics + MACS2 peaks already computed) and re-draws the three violin
## panels exactly as in Multiome_update_01292025.Rmd, but saved as PDF in the
## git repo. (qc_filtering.R is the from-raw method; this just plots.)
## =====================================================================
suppressMessages({library(Seurat); library(Signac); library(ggplot2); library(scales); library(dplyr)})

rds <- "C:/Users/loril/Documents/data/multiome_081721_aggregated_data_analysis/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks.rds"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/01_preprocessing/output/figures/qc"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
sex_cols <- c(female = "#F39AC9", male = "#4A6FE3")

lot6 <- readRDS(rds)
lot6$sex <- ifelse(grepl("-1$", colnames(lot6)), "female",
             ifelse(grepl("-2$", colnames(lot6)), "male", NA))
Idents(lot6) <- "sex"
DefaultAssay(lot6) <- "RNA"
lot6[["percent.mt"]] <- PercentageFeatureSet(lot6, pattern = "^mt-")

## only plot features that actually exist (meta columns or assay features)
have <- function(fs) fs[fs %in% c(colnames(lot6@meta.data), rownames(lot6))]

vln <- function(features, file, w = 12, h = 4) {
  features <- have(features)
  p <- VlnPlot(lot6, features = features, pt.size = 0.1, ncol = length(features),
               cols = sex_cols) & scale_y_continuous(labels = comma)
  ggsave(file.path(out, file), p, width = w, height = h, units = "in")
  cat("wrote", file, "(", paste(features, collapse = ", "), ")\n")
}

vln(c("nFeature_RNA", "nCount_RNA", "percent.mt"),                      "SupplFig1a_RNA_QC_violins_bySex.pdf",  9)
vln(c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "pct_reads_in_peaks"), "SupplFig1a_ATAC_QC_violins_bySex.pdf")

## median summary table (matches the Rmd's group_by summary)
summ <- lot6@meta.data %>% group_by(sex) %>%
  summarise(median_nFeature_RNA = median(nFeature_RNA),
            median_nCount_RNA   = median(nCount_RNA),
            median_percent_mt   = round(median(percent.mt), 2),
            n_cells = dplyr::n(), .groups = "drop")
write.csv(summ, file.path(out, "SupplFig1a_QC_median_summary_bySex.csv"), row.names = FALSE)
print(summ)
cat("\nDone. Panel-1a PDFs in:\n  ", out, "\n")
