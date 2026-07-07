# Multiome RNA — single-nuclei (10x Multiome, lot6 081721)

Single-nuclei **RNA** analysis from the 10x Multiome dataset (male & female GBM astrocytes).
(The ATAC side of the multiome is in `../../03_atac_analysis/`.)

## Scripts
- `rpca_integration.R` — Seurat / Signac RPCA integration of the multiome (RNA assay).
- `Figures_SupplFig1b_RPCA_UMAP.R` — Suppl. Fig. 1b/e: RNA UMAP on the RPCA-integrated sample
  (colored by sex + per-sex split) and the Egr1 `FeaturePlot`. Writes PDFs to
  `output/figures/suppl_fig1/`. (Moved here from `03_atac_analysis/`; it is an RNA panel.)
- `gene_de.R` — single-cell differential expression (Seurat `FindMarkers`, male vs female)
  on the RNA assay.
- `suppl_table_DE_male_vs_female.R` — writes the multiome RNA DE supplementary table
  to `output/tables/SupplTable_multiome_DE_male_vs_female.csv`.

## Data
These scripts load the processed multiome Seurat object
(`female_male_aggregated_081722_..._chromVARadded_111425.rds`), which is **not committed**
to git due to size. Edit the `readRDS()` path to point to your local / cluster copy.
