# Supplementary Fig. 2 (cortex development) — HTCF handoff / runbook

Self-contained context to continue this work from an HTCF session. Everything
below was reconstructed from the code in this repo + the original `.Rmd`s.

---

## Goal & current status

Supp Fig. 2 = EGR1/TF motif activity across human cortex development (Roussos
snMultiome, split by sex, 5 stages: LateFetal → Infant → Child → Adol → Adult).

| Panel | What | Status |
|-------|------|--------|
| **a** | cell numbers by stage×sex | ⏳ 3/5 stages done locally (Child/Infant/LateFetal). **Adol+Adult still needed.** |
| **b** | RNA (RPCA) + ATAC (rLSI) UMAPs by sex | ⏳ same — Adol+Adult needed |
| **c** | RENIN motif montage (2×5 volcanoes) | ✅ done; 10 motif CSVs committed, montage rebuilds locally |

**The only thing outstanding: run panels a/b for all 5 stages** (adds Adol+Adult).
Their objects are HTCF-only, which is why they're missing from the local preview.

---

## Data provenance (for manuscript methods/citation)

- **Paper:** Zhu K, Bendl J, Rahman S, … Roussos P. "Multi-omic profiling of the
  developing human cerebral cortex at the single-cell level." *Sci Adv* 2023;
  9(41):eadg3754. doi:10.1126/sciadv.adg3754.
- **Deposits:** GEO **GSE204684** (raw + Cell Ranger ARC outputs), Broad Single
  Cell Portal **SCP1859**, CELLxGENE, UCSC hub `dual-assay.s3.amazonaws.com`,
  code Zenodo 10.5281/zenodo.7703253. hg38, 6 stages, ~45,549 nuclei.
- **How the objects were made** (NOT from FASTQ): downloaded the deposited Cell
  Ranger outputs (`filtered_feature_bc_matrix.h5` + `atac_fragments.tsv.gz`) via
  Cyberduck to the Mac drive `/Volumes/DUAL_DRIVE/data_mining/roussos_human_cerebral_cortex_science23/`,
  then locally: build RNA+ATAC object → QC filter → MACS2 peaks → add motifs →
  upload to HTCF. Sex from barcode suffix (mapping FLIPS per stage — see gotchas).

---

## File locations on HTCF

All cortex scripts read from the **home dir**, not LTS:

```
/home/lllaci/data/cortex_development/
    <Stage>_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds
    (Stage = LaFet, Inf, Child, Adol, Adult, EaFet)
```

To confirm / locate on the cluster (incl. any LTS copy):
```sh
find /home/lllaci /lts/rmlab/rmlab_shared3/l.llaci /scratch/rmlab/rmlab_shared3/l.llaci \
  \( -iname '*Adol*Roussos*' -o -iname '*Adult*Roussos*' \) 2>/dev/null
```

---

## Pipeline scripts (in this repo, `06_cortex_development/`), full order

Every script loops over all 6 stages incl. Adol/Adult.

1. **Download** (manual): Cell Ranger outputs from GSE204684, one sample/stage.
2. `cortex_dev_build_objects.R` → `<Stage>_..._RNAandPeakAssaysAvailable.rds`
3. `cortex_dev_add_motifs.R` → `<Stage>_..._withMotifs_V4_122825.rds`
4. `cortex_dev_RENIN_pipeline.R` → 10 `<Stage>_<Sex>_AllCells_all_motifs.csv`
5. Figures: `Figure_SupplFig1_RENIN_motif_montage.R` (panel c, local) +
   `Figure_SupplFig2_cellnumbers_umaps.R` (panels a/b)

---

## TO FINISH IT — run panels a/b for all 5 stages on HTCF

The `_withMotifs_V4` objects already exist on HTCF, so **no rebuild is needed** —
just run the figure. From a copy of `06_cortex_development/` on HTCF:

```sh
# 1. make R available (NO `module load` — lab rule). Edit the sbatch's commented
#    conda line, e.g.:  source ~/miniconda3/etc/profile.d/conda.sh && conda activate r_multiome
# 2. submit (128G, 8h, all 5 stages):
sbatch run_SupplFig2_cellnumbers_umaps_HTCF.sbatch

# monitor:
squeue -u l.llaci

# outputs land in ./output/ :
#   SupplFig2_cortex_cellnumbers.pdf, SupplFig2_cortex_UMAPs.pdf, CellCounts_AllCells_byStageSex.csv
```

The sbatch sets `DATA_DIR=/home/lllaci/data/cortex_development`,
`OUT_DIR=$PWD/output`, leaves `CORTEX_LOCAL_PREVIEW` unset (→ all 5 stages), and
fails loudly if `Rscript` isn't on PATH.

Pull the results back to the laptop:
```sh
scp 'login.htcf.wustl.edu:<remote>/output/SupplFig2_*' \
    'C:/Users/loril/Documents/GitHub/Egr1-Llaci/06_cortex_development/output/'
```

### Fully-local alternative (no HTCF compute)
`scp` the two Adol/Adult objects down, then run
`Figure_SupplFig2_cellnumbers_umaps.R` locally (needs ~128G RAM for rLSI). Or
rebuild from GSE204684 Cell Ranger outputs via steps 2→3→5 above.

---

## Gotchas

- **No `module load`** on this HPC — activate R via conda/venv in the sbatch.
- **Sex mapping flips per stage** (barcode suffix): Adult/Infant/EaFet → `-1`=female;
  Adol/Child/LaFet → `-1`=male. Encoded in `cortex_dev_add_motifs.R`.
- Motifs are irrelevant to panels a/b (they use only RNA + peaks assays), so any
  aggregated object works there; the `withMotifs_V4` ones are just what's on HTCF.
- SSH from Windows: use plain PowerShell `ssh`/`scp` (Windows OpenSSH). Do NOT use
  `ssh -fN` (fails) or Git Bash's ssh with ControlMaster (mux errors on Windows).

---

## Git state (repo `Egr1-Llaci`, branch `main`)

Two local commits added this work (NOT yet pushed):
- `165f592` — cortex data-processing + Supp Fig 2 a/b code
- `d0d95d7` — STEP 1 motif script + full-pipeline README

Interim 3-stage preview PDFs left uncommitted (to be replaced by the 5-stage run).
