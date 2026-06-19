# Supplementary Fig. 1 — "Egr1 activity is enriched in male GBM cells"

Panel-by-panel map of scripts, inputs, and status. Panels **a–e** use the lot6 mouse
multiome (081721); panels **f–h** use the cortex-development dataset. All multiome panels
need Signac/Seurat + the HTCF data, so they run **on the HTCF cluster**, not locally.

Raw multiome data (HTCF, aggregated lot6/081721):
`/lts/rmlab/rmlab_shared3/l.llaci/output/seurat_results/multiome_081721_analysis/multiome_081721_aggregated/outs/`
( `filtered_feature_bc_matrix.h5` 268M + `atac_fragments.tsv.gz` 6.4G )

| Panel | Content | Script | Status |
|-------|---------|--------|--------|
| **a** | QC violins — Genes, UMI, %Mito, Peaks, ATAC fragments, TSS enrichment (M vs F) | `01_preprocessing/qc_filtering.R` | ✅ script ready; needs HTCF data (paths still placeholder `data/`) |
| **b** | RNA UMAP | `03_atac_analysis/Suppl_Fig1bce_multiome_UMAPs.R` | ✅ written; needs HTCF data |
| **c** | ATAC-seq UMAP | `03_atac_analysis/Suppl_Fig1bce_multiome_UMAPs.R` (also computed in `DARs.R`) | ✅ written; needs HTCF data |
| **d** | Promoter / gene-body / intergenic distribution of peaks | `03_atac_analysis/Suppl_Fig1d_peak_annotation.R` | ✅ written (ChIPseeker on MACS2 `peaks`); needs HTCF data |
| **e** | Male / Female feature UMAPs (e.g. Egr1) | `03_atac_analysis/Suppl_Fig1bce_multiome_UMAPs.R` | ✅ written; needs HTCF data |
| **f** | QC for cortex development | `Roussos_cortex_development_RENIN.Rmd` | ⬜ to extract |
| **g** | Cell numbers by age × sex (bar) | cortex-dev cell-number CSVs (Adol/Adult/Child/Infant/LateFetal × M/F) | ⬜ to extract (per-stage `CellCounts_*` CSVs live on HTCF) |
| **h** | Cortex development RNA + ATAC UMAPs | `Roussos_cortex_development_RENIN.Rmd` | ⬜ to extract |
| **(RENIN motif montage)** | TF motif enrichment, M vs F × 5 stages (2×5 volcano grid) | `06_cortex_development/Figure_SupplFig1_RENIN_motif_montage.R` (pipeline: `06_cortex_development/cortex_dev_RENIN_pipeline.R`) | ✅ extracted; CSVs committed, figure rebuilds locally |

## Notes
- **DAR / peak analyses use the MACS2 `peaks` assay** (built by `01_preprocessing/macs2_peak_calling.R`),
  not the Cell Ranger ARC `ATAC` assay. The DAR volcano is `03_atac_analysis/DARs.R` /
  `Figure_1C_DAR_volcano.R`; panel **d** annotates that same MACS2 peak set.
- Panel **d** original (`Multiome_update_01292025.Rmd`) annotated the `ATAC` assay; the git
  script defaults to `peaks` for consistency with the DAR — set `peak_assay <- "ATAC"` to
  reproduce the original exactly.
- ChIPseeker window: `tssRegion = c(-2000, 500)`, `TxDb.Mmusculus.UCSC.mm10.knownGene`, `org.Mm.eg.db`.
