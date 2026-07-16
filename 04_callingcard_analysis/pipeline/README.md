# Calling Cards pipeline scripts (reads → hops → peaks)

Read-processing and peak-calling scripts for the Egr1 Calling Cards experiment
(Fig. 3a/b, d–f). These run on HTCF; large inputs/outputs (fastqs, BAMs, qbeds,
HOMER output) are **not** in the repo. Provided here for provenance/reproducibility.

## Pipeline overview

```
raw reads
  └─ Cutadapt v2.9         adapter/quality trim
  └─ Bowtie2 v2.3.5        align to mm10
  └─ SAMtools v1.13        sort / index
  └─ tag + annotate        sample / i7 / i5 / SRT barcodes; annotate piggyBac (PB)
                           insertion sites at TTAA; collapse to unique insertions ("hops")
  └─ qbed                  per-sample insertion files  (UMI-tools v1.0.0)
  └─ chrY filter           remove Y-chromosome hops (mismapping artifact; present in both sexes)
  └─ peak calling          pycallingcards cc.pp.call_peaks, MACCs method (see note below)
```

## Folder contents

### `call_hops/` — reads → qbed insertions
- `call_hops.sh` — main SLURM driver (Cutadapt → Bowtie2/mm10 → SAMtools → tag → annotate → qbed).
- `TagBam.py`, `TagBamWithSrtBC.py`, `AnnotateInsertionSites.py`, `BamToCallingCard.py` —
  mammalian Calling Cards read-tagging / insertion-annotation tools (Mitra lab; Moudgil et al., 2020).
- `call_hops_manifest.csv` — sample → barcode/index manifest.
- `srt_barcode_whitelist.txt` — self-reporting-transposon (SRT) barcode whitelist.
- `downsample_unique_alignments.sh` — optional alignment downsampling (template; **not used for the paper's final peaks**).

### `call_peaks/` — HTCF peak-calling scripts (exploratory)
Various `call_peaks_*` scripts (window1000, downsampled, `unf_macs`) + `CCFtools.py`.
**NOTE:** these are HTCF exploratory variants. The peaks used in the paper were
called in the analysis notebook with **pycallingcards `cc.pp.call_peaks` (MACCs;
window = 300, step = 150, `pvalue_cutoffbg` = 0.05, `pvalue_cutoffTTAA` = 0.001,
`lam_win_size` = 1e6, pseudocounts = 0.1), identically for both sexes** — see the
notebook on the `pipeline` branch and `../panel_ab_README.md`.

### `find_motifs/` — HOMER motif analysis
- `tf_homer.sh`, `motif_calling_from_qbed.sh`, `homer_list.csv`.

## Software versions
Cutadapt 2.9 · Bowtie2 2.3.5 · SAMtools 1.13 · UMI-tools 1.0.0 · pycallingcards (Guo et al., 2024) · HOMER.

## References
- Moudgil A, et al. Self-Reporting Transposons… *Cell* 2020;182(4):992–1008.e21. (CC method / read-tagging tools)
- Guo J, et al. Pycallingcards… *Bioinformatics* 2024;40(2):btae070. (peak calling / analysis)
