#!/usr/bin/env python
"""
Male vs Female Calling-Card peak overlap -> Venn + peak sets   [HTCF version]

Run on HTCF:
    srun --mem=16G --cpus-per-task=2 -p interactive --pty /bin/bash -l
    conda activate <your pycallingcards env>       # needs pandas, pyranges, matplotlib, matplotlib-venn
    python CC_peak_overlap_venn_htcf.py

or non-interactively:
    sbatch --mem=16G --cpus-per-task=2 --wrap "python CC_peak_overlap_venn_htcf.py"

Paths default to the HTCF layout below; override with --male / --female / --out.

-----------------------------------------------------------------------------
WHY THIS REPLACES THE NOTEBOOK CELL

The original used, on `peak_annotation_{Male,Female}_20000`:
    common        = male_gr.intersect(female_gr)
    unique_peaks1 = male_gr.subtract(female_gr)
    unique_peaks2 = female_gr.subtract(male_gr)

Two independent bugs:

1. WRONG INPUT. `*_20000` is the 20 kb-nearest-gene-FILTERED subset
   (male 9,456 / female 9,603), not the raw calls (male 11,977 / female 11,772).
   ~2,300 peaks never entered the Venn. The 20 kb filter belongs to the CC-vs-DE
   regulon step and has nothing to do with sharing between sexes.

2. `.intersect()` / `.subtract()` are INTERVAL ARITHMETIC, not classification.
   `.subtract()` erases the overlapping portion and keeps the leftover STUB --
   still counting it as "unique" even though the peak is shared -- and can SPLIT
   one peak into two stubs counted as two peaks. `.intersect()` keeps only the
   clipped overlap sliver, one row per overlapping PAIR. So a single male peak
   overlapping a female peak was counted TWICE (once unique, once shared), with
   neither piece retaining the real peak's coordinates.

   Verified: the old numbers (7,357 / 7,495 / 3,187) are reproduced EXACTLY by
   feeding the 20 kb subset through .subtract()/.intersect(). They are fragment
   counts, not peak counts, which is also why "shared" came out symmetric --
   impossible for a real many-to-many intersect.

THE FIX: `.overlap()` / `.overlap(invert=True)` == `bedtools intersect -u` / `-v`.
Whole peaks, nothing cut, every peak in exactly one category.
-----------------------------------------------------------------------------
"""

import argparse
import os

import pandas as pd
import pyranges as pr
import matplotlib
matplotlib.use("Agg")                      # headless: no X11 on compute nodes
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# ---- HTCF defaults ---------------------------------------------------------
# NOTE 2026-07-17: verified by ls on HTCF -- the working dir is on /lts/, NOT
# /scratch/, and Male_Egr1CC_peaks_072825.txt is NOT there. The MACCs calls that
# ARE present are the headerless *_window300_p05.bed files below (male 11,977 /
# female 11,772). Confirmed byte-identical to the *_p05_TTAAp001_091325_071526_p05.bed
# variant on the `pipeline` branch -- same file, two names, so either is fine.
BASE = "/lts/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper"
DEF_MALE = f"{BASE}/Egr1CC_peak_MaleEgr1_VS_MaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed"
DEF_FEMALE = f"{BASE}/Egr1CC_peak_FemaleEgr1_VS_FemaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed"
DEF_OUT = f"{BASE}/peak_overlap"

ap = argparse.ArgumentParser()
ap.add_argument("--male", default=DEF_MALE, help="raw male MACCs peak .txt (NOT the *_20kb file)")
ap.add_argument("--female", default=DEF_FEMALE, help="raw female MACCs peak .txt")
ap.add_argument("--out", default=DEF_OUT, help="output directory")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)


def load_peaks(path, expected=None):
    """Raw MACCs call -> PyRanges. Accepts the caller's .txt (header 'Chr/Start/End') or a 3-col BED."""
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: not found: {path}")
    head = open(path).readline()
    if head.lower().startswith(("chr\t", "chrom")):          # caller .txt with header
        df = pd.read_csv(path, sep="\t")
        df = df.rename(columns={"Chr": "Chromosome", "chrom": "Chromosome",
                                "start": "Start", "end": "End"})
    else:                                                     # headerless BED
        df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                         names=["Chromosome", "Start", "End"])
    df = df[["Chromosome", "Start", "End"]]
    gr = pr.PyRanges(df)
    if expected is not None and len(gr) != expected:
        print(f"  WARNING: {os.path.basename(path)} has {len(gr)} peaks, expected {expected} "
              f"-- is this the raw call and not the 20 kb-filtered file?")
    return gr


def write_bed(gr_or_df, name):
    df = gr_or_df.df if hasattr(gr_or_df, "df") else gr_or_df
    df[["Chromosome", "Start", "End"]].sort_values(["Chromosome", "Start"]).to_csv(
        os.path.join(args.out, name + ".bed"), sep="\t", header=False, index=False)


male = load_peaks(args.male, expected=11977)
female = load_peaks(args.female, expected=11772)
print(f"Raw peaks:  male {len(male)}   female {len(female)}\n")

# ---- (1) PEAK-CENTRIC  ==  bedtools intersect -u / -v ----------------------
male_unique = male.overlap(female, invert=True)
male_shared = male.overlap(female)
female_unique = female.overlap(male, invert=True)
female_shared = female.overlap(male)

assert len(male_unique) + len(male_shared) == len(male), "peak-centric must partition male"
assert len(female_unique) + len(female_shared) == len(female), "peak-centric must partition female"

print("(1) PEAK-CENTRIC (whole peaks, >=1 bp overlap):")
print(f"    male   {len(male):5d} = {len(male_unique):5d} unique + {len(male_shared):5d} overlapping-female")
print(f"    female {len(female):5d} = {len(female_unique):5d} unique + {len(female_shared):5d} overlapping-male")
print(f"    shared is ASYMMETRIC ({len(male_shared)} vs {len(female_shared)}) -- correct: overlap is many-to-many.\n")

write_bed(male_unique, "CC_male_unique_peakCentric")
write_bed(female_unique, "CC_female_unique_peakCentric")
write_bed(male_shared, "CC_male_shared_peakCentric_maleCoords")
write_bed(female_shared, "CC_female_shared_peakCentric_femaleCoords")

# ---- (2) REGION-CENTRIC  ==  merge union, then classify --------------------
union = pr.concat([male, female]).merge()
udf = union.df.copy()
udf["in_male"] = union.count_overlaps(male).df["NumberOverlaps"].values > 0
udf["in_female"] = union.count_overlaps(female).df["NumberOverlaps"].values > 0

male_only = udf[udf.in_male & ~udf.in_female]
shared = udf[udf.in_male & udf.in_female]
female_only = udf[~udf.in_male & udf.in_female]
assert len(male_only) + len(shared) + len(female_only) == len(udf), "regions must partition the union"

print("(2) REGION-CENTRIC (merged union; disjoint, sums exactly) <- the Venn:")
print(f"    male-only {len(male_only)} + shared {len(shared)} + female-only {len(female_only)} = {len(udf)}")
print(f"    male circle = {len(male_only) + len(shared)}   female circle = {len(female_only) + len(shared)}\n")

write_bed(male_only, "CC_male_only_regions")
write_bed(shared, "CC_shared_regions")
write_bed(female_only, "CC_female_only_regions")

pd.DataFrame({
    "Method": ["peak-centric"] * 4 + ["region-centric"] * 4,
    "Set": ["male total", "male unique", "female total", "female unique",
            "male-only", "shared", "female-only", "union total"],
    "N": [len(male), len(male_unique), len(female), len(female_unique),
          len(male_only), len(shared), len(female_only), len(udf)],
}).to_csv(os.path.join(args.out, "CC_peak_overlap_summary.csv"), index=False)

# ---- Venn ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 7))
v = venn2(subsets=(len(male_only), len(female_only), len(shared)),
          set_labels=("Male", "Female"), ax=ax)

for sid, colour in {"10": "#4C7FB8", "01": "#C4622D", "11": "#8A5344"}.items():
    p = v.get_patch_by_id(sid)
    if p:
        p.set_color(colour); p.set_edgecolor("black"); p.set_linewidth(2); p.set_alpha(0.55)

# Font sizes per project figure rules: counts/labels >= 14 pt, title >= 18 pt.
for sid in ("10", "01", "11"):
    lbl = v.get_label_by_id(sid)
    if lbl:
        lbl.set_fontsize(17)
        lbl.set_text(f"{int(lbl.get_text()):,}")
for lbl in v.set_labels:
    if lbl:
        lbl.set_fontsize(18); lbl.set_fontweight("bold")

ax.set_title("Egr1 Calling-Card peak overlap", fontsize=18, fontweight="bold", pad=28)
ax.text(0.5, 1.02, "peaks called separately per sex", transform=ax.transAxes,
        ha="center", fontsize=14, color="grey")
ax.text(0.5, -0.04, f"merged union regions, $\\geq$1 bp overlap  |  n = {len(udf):,} total",
        transform=ax.transAxes, ha="center", fontsize=14, color="grey")

plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(os.path.join(args.out, f"SupplFig4b_CC_peak_overlap_Venn.{ext}"),
                dpi=300, bbox_inches="tight")

print("Expected (verified on the raw calls):")
print("  peak-centric   male 8,265 uniq + 3,712 shared | female 8,053 uniq + 3,719 shared")
print("  region-centric male-only 8,204 + shared 3,613 + female-only 7,986 = 19,803")
print(f"\nWrote figure + 7 beds + summary to:\n  {args.out}")
