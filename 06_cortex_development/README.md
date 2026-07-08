# 06_cortex_development — RENIN motif enrichment across human cortex development

RENIN cis-regulatory element scoring + TF motif-enrichment on the Roussos
cortex-development snMultiome (RNA + ATAC), split by sex across five
developmental stages (Late Fetal → Infant → Child → Adolescence → Adult).
This is the cortex-development arm of **Supplementary Fig. 2** (the RENIN motif
montage is panel **c**; cell-number bar = a, RNA/ATAC UMAPs = b — see
`../Suppl_Fig2_panel_map.md`).

## Scripts

| Script | Step | Runs |
|--------|------|------|
| `cortex_dev_build_objects.R` | **STEP 0.** Per stage, from the deposited Roussos Cell Ranger ARC outputs: build RNA + ATAC object (hg38) → QC filter → MACS2 peaks → add `peaks` assay → save `<Stage>_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds`. | **HTCF cluster** (MACS2 + full fragment files). Submit with `srun`/`sbatch`. |
| `cortex_dev_RENIN_pipeline.R` | Full RENIN pipeline: per stage, Harmony (batch = sex) → DEGs (F vs M) → pseudocells → `run_peak_aen` → split CREs into male-/female-enriched → `FindMotifs` → per-panel volcano. Writes `<Stage>_<Sex>_AllCells_all_motifs.csv`. (Motifs added to the STEP 0 objects first — see its "STAGE 0" header.) | **HTCF cluster** (Seurat/Signac/RENIN + large `*_withMotifs_V4` Seurat objects). Submit with `srun`/`sbatch`. |
| `Figure_SupplFig1_RENIN_motif_montage.R` | Plot-only: reads the 10 motif tables and assembles the 2×5 volcano montage (Supp Fig. 2c). | **Locally**, in seconds. No HTCF data needed. |
| `Figure_SupplFig2_cellnumbers_umaps.R` | Supp Fig. 2 panels a/b: cell numbers by stage×sex + per-stage RNA (RPCA) / ATAC (rLSI) UMAPs coloured by sex. Full 5-stage run on HTCF via `run_SupplFig2_cellnumbers_umaps_HTCF.sbatch`. | Local preview = 3 stages; full 5 stages on **HTCF**. |

## Data

`data_motifs/` — the 10 motif-enrichment tables (`<Stage>_<Sex>_AllCells_all_motifs.csv`,
746 motifs each), the committed output of the pipeline's `FindMotifs` step. These
are the inputs to the montage figure; committing them lets the figure be rebuilt
without re-running the cluster pipeline.

The upstream `*_withMotifs_V4_122825.rds` Seurat objects (one per stage) are large
and live on HTCF — they are **not** committed. Motif PFMs: JASPAR2020 CORE
vertebrates (`AddMotifs`, done in `Roussos_cortex_development_RENIN.Rmd`).

## Data provenance

Published external data: **Zhu K, Bendl J, Rahman S, et al. "Multi-omic profiling of
the developing human cerebral cortex at the single-cell level." Sci Adv 2023;9(41):eadg3754**
(doi:10.1126/sciadv.adg3754). 10X Multiome (snRNA + snATAC), hg38, six developmental
stages (Early Fetal → Late Fetal → Infant → Child → Adolescence → Adult).
Deposits: GEO **GSE204684** (raw + Cell Ranger outputs), Broad Single Cell Portal
**SCP1859**, CELLxGENE, UCSC hub `dual-assay.s3.amazonaws.com`, code Zenodo
10.5281/zenodo.7703253. We start from the deposited Cell Ranger ARC outputs (not FASTQ);
`cortex_dev_build_objects.R` turns them into the per-stage objects used here.

## Outputs

`output/SupplFig1_RENIN_motif_montage.{png,pdf}` — assembled 2×5 figure
(rows = Males / Females, columns = developmental stage). Per-panel PDFs in
`output/panels/`. Green = significant (`p.adjust ≤ 0.05`); top-20 motifs per
panel labeled.

## Reproduce the figure

```sh
Rscript 06_cortex_development/Figure_SupplFig1_RENIN_motif_montage.R
```
