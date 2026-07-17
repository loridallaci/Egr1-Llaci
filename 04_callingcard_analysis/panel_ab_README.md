# Fig. 3a/b — Egr1 Calling Card binding

How the Calling Card (CC) panels of Figure 3 are produced. The peak-calling and
signal-track steps run on HTCF; the scripts here are committed for provenance.
Large raw inputs (qbeds, bigWigs) live on HTCF/LTS and are listed under
**Data availability** below.

## Pipeline

| Step | Script / notebook | Output |
|------|-------------------|--------|
| 1. Call peaks per sex | `calling_card_peak_pipeline.py` (`pipeline` branch); pycallingcards `cc.pp.call_peaks` (MACCs) | `Egr1CC_peak_{Male,Female}Egr1_VS_{Male,Female}WT_MACC2_*.bed` (male 11,977 / female 11,772) |
| 2a. Filter to 20 kb of a gene (**regulon step ONLY**) | `1_CCpeaks_to_20kb_overlap.R` | `{Male,Female}_Egr1CC_peaks_20kbThreshhold_*.txt` (9,456 / 9,603) |
| 2b. Classify shared/unique (**from the FULL calls, not the 20 kb subset**) | **`CC_peak_overlap_venn.py`** (or `.R`) | `output/peak_overlap/CC_{male_only,shared,female_only}_regions.bed` (8,204 / 3,613 / 7,986) + `*_peakCentric.bed` + Supp Fig 4a Venn |
| 3. CPM bigWigs | **`panel_a_1_make_cpm_bigwigs.sh`** | `cpm_bigwigs/{Male,Female}_Egr1_CPM.bw` |
| 4. **Panel a** heatmap — one panel per sex, over **all** that sex's peaks | **`panel_a_2_plot_heatmap.sh`** (deepTools) | `output/figures/fig3a_heatmap/Egr1CC_originalPeaks_CPM_matrix_{male,female}Only_Fig3.pdf` |
| 5. **Panel b** browser tracks | `Final5_…for_paper_April2026.ipynb` (`pipeline` branch) | `{Male,Female}_Egr1_WT_tracks_050726_*_{Ptk2b,Frzb,Nab1}.pdf` |
| 6. **Panels d/e** CC ∩ DE regulon | `2_CCxDE_to_regulonSummary.R`, `3_CCxDE_full_summary_with_multivariate.R` | `overlap_{Male,Female}CC_20kb_NearestGene_vs_Egr1g3vsNeg1_*.csv` |

### Panel a — what it is, and two traps (corrected 2026-07-17)
Panel a is **two heatmaps side by side (Male | Female)**, each over that sex's **full**
peak set (11,977 / 11,772), centered ±1 kb. It is **not** grouped into
shared / male-unique / female-unique, and it is **not** 20 kb-filtered — verified on
LTS (`wc -l centered_peaks_{males,females}.bed` → 11977 / 11772).

- **Trap 1:** the output filenames say `originalPeaks`. That is copy-paste residue from
  earlier attempts in the working-notes file — it does **not** describe the input.
- **Trap 2:** an earlier version of `panel_a_2_plot_heatmap.sh` plotted a three-group
  split from the `*_originalPeaks.bed` fragment sets. That **did not reproduce panel a**.
  It has been replaced.

The `*_originalPeaks.bed` sets (7,357 / 7,495 / 3,187) are **superseded and should not be
used**: they were built by pyranges `.subtract()`/`.intersect()` (interval arithmetic —
trims and splits peaks into fragments) applied to the 20 kb-filtered subset. Use the
`output/peak_overlap/` sets from step 2b instead.

## Peak calling — parameters
`cc.pp.call_peaks(method="MACCs", reference="mm10", window_size=300,
step_size=150, pvalue_cutoffbg=0.05, pvalue_cutoffTTAA=0.001, lam_win_size=1e6,
pseudocounts=0.1)`, run per sex (Egr1 vs matched WT), chrY-filtered inputs.

The **same** parameters were used for both sexes (`pvalue_cutoffbg = 0.05`,
`pvalue_cutoffTTAA = 0.001`); verified by reproducing the used peak beds
(male 11,977 / female 11,772 peaks) byte-for-byte from these settings. There is
no sex-specific background-cutoff asymmetry.

## Peak sets — which is which
- **Panel a heatmap and panels d/e use the SEPARATELY-CALLED unique/shared
  peaks** (`*_originalPeaks.bed`, `*_calledSeparately_*.txt`) — consistent with the
  Methods ("peaks were called on each sex separately").
- A combined-call *biased* set also exists (`Egr1_*_biased_peaks_calledtogether_*`)
  and an older `heatmap_CC.sh` used it; **panel a should use the separate-call set**
  (this is what `panel_a_2_plot_heatmap.sh` does).

## Data availability (HTCF/LTS — not in repo)
- CC insertions (qbeds): `/lts/rmlab/rmlab_shared3/l.llaci/Egr1_paper/Egr1_CC_Lori/101524_Egr1CC_analysis/output_and_analysis/insertions/`
- CPM bigWigs: `/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/cpm_bigwigs/`
- Working dir (full `pipeline`-branch clone): `/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper/`

The full working code (heatmap scripts, notebooks, peak beds) is on the
**`pipeline`** branch of this repo.
