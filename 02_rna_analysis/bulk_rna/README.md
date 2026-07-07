# Bulk RNA-seq — Egr1 CRISPRi knockdown

Bulk (BRB-seq) RNA-seq of male & female GBM astrocytes with CRISPRi knockdown of **Egr1**
(guides gRNA2 / gRNA3) vs non-targeting control (**Neg gRNA1**).

## Scripts
- `DE_Egr1KD_gRNA3_vs_Neg1.R` — DESeq2 differential expression (Egr1 gRNA3 vs Neg gRNA1,
  per sex) + QC (sample-distance heatmap & PCA, **before and after** outlier removal).
- `Egr1KD_bulkRNA_forPaper.Rmd` — full paper analysis (R Markdown, knits to HTML).
- `Egr1KD_volcano.R` — volcano plots of the DE results.
- `Enrichr_Egr1KD_bulk.R` — pathway / GO enrichment (Enrichr) of the DEGs.
- `CC_vs_bulkRNA_eachDirection_geneOverlaps.R` — overlap of bulk DEGs with CallingCard targets.

## Data
- `data_raw/` — raw inputs: `Dedup_Counts.txt`, `GeneInfo_Egr1KD.csv`, `Meta_data_Egr1KD.csv`.
- `output_DE_Egr1KD_gRNA3_vs_Neg1/` — DE tables + `qc/` (PCA & heatmap PDFs, color-coded by group).

## Pipeline
`DE_Egr1KD_gRNA3_vs_Neg1.R`  →  `../../04_callingcard_analysis/`  →  `../../07_tcga_survival/`
