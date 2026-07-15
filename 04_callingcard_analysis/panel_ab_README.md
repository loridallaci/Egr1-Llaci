# Fig. 3a/b — Egr1 Calling Card binding

How the Calling Card (CC) panels of Figure 3 are produced. The peak-calling and
signal-track steps run on HTCF; the scripts here are committed for provenance.
Large raw inputs (qbeds, bigWigs) live on HTCF/LTS and are listed under
**Data availability** below.

## Pipeline

| Step | Script / notebook | Output |
|------|-------------------|--------|
| 1. Call peaks per sex | `calling_card_peak_pipeline.py` (`pipeline` branch); pycallingcards `cc.pp.call_peaks` (MACCs) | `Egr1CC_peak_{Male,Female}Egr1_VS_{Male,Female}WT_MACC2_*.bed` |
| 2. Filter to 20 kb of a gene, classify shared/unique | `1_CCpeaks_to_20kb_overlap.R` (+ peak-comparison notebooks, `pipeline` branch) | `{Male,Female}_Egr1CC_peaks_20kbThreshhold_*.txt`; `{male,female}_unique_peaks_originalPeaks.bed`, `shared_peaks_originalPeaks.bed` |
| 3. CPM bigWigs | **`panel_a_1_make_cpm_bigwigs.sh`** | `cpm_bigwigs/{Male,Female}_Egr1_CPM.bw` |
| 4. **Panel a** heatmap | **`panel_a_2_plot_heatmap.sh`** (deepTools) | `Figure3a_Egr1CC_separateCall_CPM.pdf` |
| 5. **Panel b** browser tracks | `Final5_…for_paper_April2026.ipynb` (`pipeline` branch) | `{Male,Female}_Egr1_WT_tracks_050726_*_{Ptk2b,Frzb,Nab1}.pdf` |
| 6. **Panels d/e** CC ∩ DE regulon | `2_CCxDE_to_regulonSummary.R`, `3_CCxDE_full_summary_with_multivariate.R` | `overlap_{Male,Female}CC_20kb_NearestGene_vs_Egr1g3vsNeg1_*.csv` |

## Peak calling — parameters
`cc.pp.call_peaks(method="MACCs", reference="mm10", window_size=300,
step_size=150, pvalue_cutoffTTAA=0.001, lam_win_size=1e6, pseudocounts=0.1)`,
run per sex (Egr1 vs matched WT), chrY-filtered inputs.

> **Note (report in Methods):** the background p-value cutoff differs by sex —
> male `pvalue_cutoffbg = 0.005`, female `pvalue_cutoffbg = 0.05`. Document/justify.

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
