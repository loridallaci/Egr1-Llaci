# Supplementary Fig. 2 — EGR1 / TF motif activity across human cortex development

Panel-by-panel map of scripts, inputs, and status. All panels use the Roussos
cortex-development snMultiome (RNA + ATAC), split by sex across five developmental
stages (Late Fetal → Infant → Child → Adolescence → Adult). The upstream Seurat
objects (`*_withMotifs_V4_122825.rds`) are large and live on **HTCF**; the motif
tables and figures are committed so the montage rebuilds locally.

(Moved here from Supplementary Fig. 1 — see `Suppl_Fig1_panel_map.md`.)

| Panel | Content | Script | Status |
|-------|---------|--------|--------|
| **a** | Cell numbers by stage × sex (bar) | `06_cortex_development/Figure_SupplFig2_cellnumbers_umaps.R` | ⏳ 3/5 stages (local preview: Child/Infant/LateFetal). Full 5-stage run: `sbatch 06_cortex_development/run_SupplFig2_cellnumbers_umaps_HTCF.sbatch` (adds Adol/Adult from HTCF) |
| **b** | RNA + ATAC UMAPs per stage, coloured by sex | `06_cortex_development/Figure_SupplFig2_cellnumbers_umaps.R` | ⏳ 3/5 stages (local preview). Full 5-stage run via same HTCF sbatch as panel a |
| **c** | RENIN TF motif enrichment, M vs F × 5 stages (2×5 volcano montage) | `06_cortex_development/Figure_SupplFig1_RENIN_motif_montage.R` (pipeline: `cortex_dev_RENIN_pipeline.R`) | ✅ extracted; 10 motif CSVs committed, figure rebuilds locally |

**No QC panel:** this is published data (Roussos et al.), so QC is stated in Methods, not shown.
Suggested Methods sentence: *"We reanalyzed the published Roussos et al. cortical-development
snMultiome; nuclei with nCount_RNA and nCount_ATAC < 40,000, nucleosome signal 0.1–2, and TSS
enrichment 0.5–9 were retained."* Sex is assigned from the barcode suffix (−1 = female, −2 = male).

## Notes
- The RENIN montage arm is fully committed under `06_cortex_development/` (scripts,
  10 `<Stage>_<Sex>_AllCells_all_motifs.csv` tables, and the assembled figure);
  see `06_cortex_development/README.md`.
- Panels a/b are extracted into `Figure_SupplFig2_cellnumbers_umaps.R`. The local
  preview (`CORTEX_LOCAL_PREVIEW=1`) covers 3 stages; the canonical 5-stage figure
  is produced on HTCF via `run_SupplFig2_cellnumbers_umaps_HTCF.sbatch`, which adds
  the HTCF-only Adol/Adult stages. Panel c (RENIN montage) is fully committed.
- The montage script keeps its original filename (`Figure_SupplFig1_RENIN_motif_montage.R`)
  for git history continuity even though the figure is now Supplementary Fig. 2.
