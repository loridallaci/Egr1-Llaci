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
| `cortex_dev_add_motifs.R` | **STEP 1.** Per stage, add JASPAR2020 CORE vertebrate motifs to the `peaks` assay (`AddMotifs`), assign sex, convert RNA to Seurat v4 → save `<Stage>_..._withMotifs_V4_122825.rds`. | **HTCF cluster.** Submit with `srun`/`sbatch`. |
| `cortex_dev_RENIN_pipeline.R` | **STEP 2.** Full RENIN pipeline: per stage, Harmony (batch = sex) → DEGs (F vs M) → pseudocells → `run_peak_aen` → split CREs into male-/female-enriched → `FindMotifs` → per-panel volcano. Writes `<Stage>_<Sex>_AllCells_all_motifs.csv`. | **HTCF cluster** (Seurat/Signac/RENIN + large `*_withMotifs_V4` Seurat objects). Submit with `srun`/`sbatch`. |
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

## Full pipeline (download → figures), all six stages incl. Adol/Adult

Every script below loops over all stages (`EaFet, LaFet, Inf, Child, Adol, Adult`);
Adol/Adult are handled identically — they are only "special" in that their large
objects live on HTCF, not locally.

1. **Download** (manual): fetch the per-stage 10X Multiome Cell Ranger ARC outputs
   (`filtered_feature_bc_matrix.h5` + `atac_fragments.tsv.gz`) from GEO **GSE204684**
   (or Broad SCP1859 / the `dual-assay` S3). One sample per stage.
2. `cortex_dev_build_objects.R` → `<Stage>_..._RNAandPeakAssaysAvailable.rds` (HTCF).
3. `cortex_dev_add_motifs.R` → `<Stage>_..._withMotifs_V4_122825.rds` (HTCF).
4. `cortex_dev_RENIN_pipeline.R` → 10 `<Stage>_<Sex>_AllCells_all_motifs.csv` (HTCF).
5. Figures: `Figure_SupplFig1_RENIN_motif_montage.R` (panel c, local) and
   `Figure_SupplFig2_cellnumbers_umaps.R` (panels a/b; all 5 stages via
   `sbatch run_SupplFig2_cellnumbers_umaps_HTCF.sbatch`).

## Reproduce the montage figure (fast, local)

```sh
Rscript 06_cortex_development/Figure_SupplFig1_RENIN_motif_montage.R
```
