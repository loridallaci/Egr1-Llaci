# Egr1-Llaci

Analysis code for the EGR1 sex-differences GBM study (single-nucleus multiome + bulk RNA-seq + TCGA survival).

## Repository structure

| Folder | Contents |
|--------|----------|
| `01_preprocessing/` | Single-nucleus multiome QC/filtering, MACS2 peak calling, and motif + chromVAR setup (`qc_filtering.R`, `macs2_peak_calling.R`, `add_motifs_chromvar.R`). |
| `02_rna_analysis/` | Bulk RNA-seq / EGR1-knockdown differential expression, Enrichr, and multiome RNA (RPCA) analysis (`gene_de.R`, `Egr1KD_bulkRNA*.R`, `rpca_integration.R`). |
| `03_atac_analysis/` | Differentially accessible regions and ATAC figures (`DARs.R`, `Figure_1B_ATACtracks.R`, `Figure_1C_DAR_volcano.R`). |
| `04_callingcard_analysis/` | EGR1 calling-card binding: CC peak–gene overlap, CC×DE regulon summary, and multivariate integration (`1_CCpeaks_to_20kb_overlap.R`, `2_CCxDE_to_regulonSummary.R`, `3_CCxDE_full_summary_with_multivariate.R`). |
| `05_integration/` | RENIN cis-regulatory element scoring and TF motif-enrichment figures (`cre_scoring_RENIN.R`, `Figure_1D_RENIN_motifs_updated_FINAL.R`, `Figure_1e_RENIN_motifs_TCGA_VennDiagram.R`). |
| `06_cortex_development/` | RENIN TF motif enrichment across human cortex development (Roussos snMultiome), males vs females over five developmental stages — Suppl. Fig. 1f–h (`cortex_dev_RENIN_pipeline.R`, `Figure_SupplFig1_RENIN_motif_montage.R`). |
| `07_tcga_survival/` | TCGA-GBM multivariable Cox survival analysis (`01_load_and_prepare_tcga_data_updated.R`, `02_multivariate_cox_regression_updated_fixedForest.R`, `03_kaplan_meier_plots.R`, `04_figures_chromvar_vs_survival.R`, `utils.R`). |
| `08_growth_curves/` | Cell growth-curve / AUC comparisons between genotypes, with p-values (`AUC_plots_Controls_2reps_comparison_between_genotypes_withpvalues.R`). |
| `data/` | TCGA-GBM expression and phenotype tables (from GlioVis, accessed 2024-06-04). |
| `data output/` | Derived result tables and figure PDFs. |
| `environment/` | `requirements.R` — R/Bioconductor package install script. |
| `figures/` | Figure-generation scripts and outputs. |
| `archive/` | Superseded earlier versions of scripts, retained for reference. |

## Pipeline order

`01_preprocessing` → `03_atac_analysis` → `05_integration` → `07_tcga_survival`
(`02_rna_analysis` covers the bulk RNA-seq / knockdown arm.)

Within `07_tcga_survival`, run `01_load_and_prepare_tcga_data_updated.R` before
`02_multivariate_cox_regression_updated_fixedForest.R`.

## Data

The TCGA-GBM expression and clinical data were obtained from the
[GlioVis](http://gliovis.bioinfo.cnio.es/) portal and are publicly available.
See `01_preprocessing/DATA_AVAILABILITY.md` for raw sequencing data availability.
