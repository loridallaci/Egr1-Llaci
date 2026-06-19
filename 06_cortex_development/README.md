# 06_cortex_development — RENIN motif enrichment across human cortex development

RENIN cis-regulatory element scoring + TF motif-enrichment on the Roussos
cortex-development snMultiome (RNA + ATAC), split by sex across five
developmental stages (Late Fetal → Infant → Child → Adolescence → Adult).
This is the cortex-development arm referenced by Suppl. Fig. 1 panels **f–h**
(see `../Suppl_Fig1_panel_map.md`).

## Scripts

| Script | Step | Runs |
|--------|------|------|
| `cortex_dev_RENIN_pipeline.R` | Full pipeline: per stage, Harmony (batch = sex) → DEGs (F vs M) → pseudocells → `run_peak_aen` → split CREs into male-/female-enriched → `FindMotifs` → per-panel volcano. Writes `<Stage>_<Sex>_AllCells_all_motifs.csv`. | **HTCF cluster** (Seurat/Signac/RENIN + large `*_withMotifs_V4` Seurat objects). Submit with `srun`/`sbatch`. |
| `Figure_SupplFig1_RENIN_motif_montage.R` | Plot-only: reads the 10 motif tables and assembles the 2×5 volcano montage. | **Locally**, in seconds. No HTCF data needed. |

## Data

`data_motifs/` — the 10 motif-enrichment tables (`<Stage>_<Sex>_AllCells_all_motifs.csv`,
746 motifs each), the committed output of the pipeline's `FindMotifs` step. These
are the inputs to the montage figure; committing them lets the figure be rebuilt
without re-running the cluster pipeline.

The upstream `*_withMotifs_V4_122825.rds` Seurat objects (one per stage) are large
and live on HTCF — they are **not** committed. Motif PFMs: JASPAR2020 CORE
vertebrates (`AddMotifs`, done in `Roussos_cortex_development_RENIN.Rmd`).

## Outputs

`output/SupplFig1_RENIN_motif_montage.{png,pdf}` — assembled 2×5 figure
(rows = Males / Females, columns = developmental stage). Per-panel PDFs in
`output/panels/`. Green = significant (`p.adjust ≤ 0.05`); top-20 motifs per
panel labeled.

## Reproduce the figure

```sh
Rscript 06_cortex_development/Figure_SupplFig1_RENIN_motif_montage.R
```
