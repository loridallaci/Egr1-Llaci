"""
Male vs Female Calling-Card peak overlap -> Venn + peak sets   [CORRECTED]

WHAT WAS WRONG WITH THE ORIGINAL NOTEBOOK CELL
----------------------------------------------
The original used:
    common       = peak_data_Male_gr.intersect(peak_data_Female_gr)
    unique_peaks1 = peak_data_Male_gr.subtract(peak_data_Female_gr)
    unique_peaks2 = peak_data_Female_gr.subtract(peak_data_Male_gr)
on inputs `peak_annotation_Male_20000` / `peak_annotation_Female_20000`.

Three independent problems:

1. WRONG INPUT. `*_20000` is the 20 kb-nearest-gene-FILTERED set (male 9,456 /
   female 9,603), not the raw calls (male 11,977 / female 11,772). Peaks >20 kb
   from any gene were silently excluded from the Venn. The 20 kb filter belongs to
   the CC-vs-DE regulon step; it has nothing to do with whether a peak is shared
   between sexes.

2. `.intersect()` and `.subtract()` are INTERVAL ARITHMETIC, not peak classification.
     - `.intersect()` returns one row per overlapping PAIR, CLIPPED to the overlap
       segment. So "shared" coordinates are intersection segments, not real peaks,
       and one male peak overlapping two female peaks yields TWO rows.
     - `.subtract()` returns the LEFTOVER FRAGMENTS after erasing the overlap. A
       peak that partly overlaps the other sex is not removed — it is trimmed, and
       can even be SPLIT into two fragments, each counted as a "unique peak".
   The counts are therefore fragment counts, not peak counts, and the categories
   do not partition the input.

3. The counts cannot be summed into a Venn. `venn2(subsets=(u1, u2, common))`
   assumes the three are disjoint counts of the SAME kind of object. Here they are
   fragments produced by three different operations.

THE FIX
-------
Use whole-peak classification (`.overlap()` / `.overlap(invert=True)`, which are
`bedtools intersect -u` / `-v`), and report BOTH honest framings:

  (1) PEAK-CENTRIC  - keeps each sex's own peaks intact. The two "shared" counts
      DISAGREE (3,712 male-side vs 3,719 female-side) because overlap is
      many-to-many. This is correct and expected; it just cannot be drawn as a
      2-circle Venn.
  (2) REGION-CENTRIC - merge the union, then ask which sexes hit each region.
      Categories are disjoint and sum exactly. THIS is what the Venn shows.

Requires: pandas, pyranges, matplotlib, matplotlib-venn
"""

import os
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# ---- Paths -----------------------------------------------------------------
CC_DIR = r"C:\Users\loril\Documents\Egr1\Egr1CC_vs_Egr1KDBulkRNA_FINAL_July2025"
OUT_DIR = r"C:\Users\loril\Documents\GitHub\Egr1-Llaci\04_callingcard_analysis\output\peak_overlap"
os.makedirs(OUT_DIR, exist_ok=True)

MALE_RAW = os.path.join(CC_DIR, "Male_Egr1CC_peaks_072825.txt")
FEMALE_RAW = os.path.join(CC_DIR, "Female_Egr1CC_peaks_072825.txt")


def load_peaks(path):
    """Raw MACCs calls -> PyRanges. NOTE: the raw file, not the *_20000 filtered one."""
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"Chr": "Chromosome"})[["Chromosome", "Start", "End"]]
    return pr.PyRanges(df)


male = load_peaks(MALE_RAW)
female = load_peaks(FEMALE_RAW)
print(f"Raw peaks:  male {len(male)}   female {len(female)}")

# ---- (1) PEAK-CENTRIC : bedtools intersect -u / -v --------------------------
# .overlap() keeps WHOLE peaks that overlap (unlike .intersect(), which clips).
male_shared = male.overlap(female)
male_unique = male.overlap(female, invert=True)
female_shared = female.overlap(male)
female_unique = female.overlap(male, invert=True)

print("\n(1) PEAK-CENTRIC (each sex's own peaks, >=1 bp overlap):")
print(f"    male   {len(male):5d} = {len(male_unique):5d} unique + {len(male_shared):5d} overlapping-female")
print(f"    female {len(female):5d} = {len(female_unique):5d} unique + {len(female_shared):5d} overlapping-male")
print(f"    NB shared is ASYMMETRIC ({len(male_shared)} vs {len(female_shared)}) - expected, overlap is many-to-many.")

for gr, name in [(male_unique, "CC_male_unique_peakCentric"),
                 (female_unique, "CC_female_unique_peakCentric"),
                 (male_shared, "CC_male_shared_peakCentric_maleCoords"),
                 (female_shared, "CC_female_shared_peakCentric_femaleCoords")]:
    gr.df[["Chromosome", "Start", "End"]].to_csv(
        os.path.join(OUT_DIR, name + ".bed"), sep="\t", header=False, index=False)

# ---- (2) REGION-CENTRIC : merge union -> disjoint categories -----------------
union = pr.concat([male, female]).merge()
in_male = union.count_overlaps(male).df["NumberOverlaps"] > 0
in_female = union.count_overlaps(female).df["NumberOverlaps"] > 0

udf = union.df
male_only = udf[in_male & ~in_female]
shared = udf[in_male & in_female]
female_only = udf[~in_male & in_female]

assert len(male_only) + len(shared) + len(female_only) == len(udf), "categories must partition the union"

print("\n(2) REGION-CENTRIC (merged union; disjoint, sums exactly) <- used for the Venn:")
print(f"    male-only {len(male_only)} + shared {len(shared)} + female-only {len(female_only)} = {len(udf)} union regions")
print(f"    male circle = {len(male_only) + len(shared)}   female circle = {len(female_only) + len(shared)}")

for df_, name in [(male_only, "CC_male_only_regions"),
                  (shared, "CC_shared_regions"),
                  (female_only, "CC_female_only_regions")]:
    df_[["Chromosome", "Start", "End"]].to_csv(
        os.path.join(OUT_DIR, name + ".bed"), sep="\t", header=False, index=False)

pd.DataFrame({
    "Method": ["peak-centric"] * 4 + ["region-centric"] * 4,
    "Set": ["male total", "male unique", "female total", "female unique",
            "male-only", "shared", "female-only", "union total"],
    "N": [len(male), len(male_unique), len(female), len(female_unique),
          len(male_only), len(shared), len(female_only), len(udf)],
}).to_csv(os.path.join(OUT_DIR, "CC_peak_overlap_summary.csv"), index=False)

# ---- Venn (region-centric) --------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 7))
v = venn2(subsets=(len(male_only), len(female_only), len(shared)),
          set_labels=("Male", "Female"), ax=ax)

for sid, colour in {"10": "#4C7FB8", "01": "#C4622D", "11": "#8A5344"}.items():
    p = v.get_patch_by_id(sid)
    if p:
        p.set_color(colour)
        p.set_edgecolor("black")
        p.set_linewidth(2)
        p.set_alpha(0.55)

# Font sizes per project figure rules: counts/labels >=14 pt, title >=18 pt.
for sid in ("10", "01", "11"):
    lbl = v.get_label_by_id(sid)
    if lbl:
        lbl.set_fontsize(17)
        lbl.set_text(f"{int(lbl.get_text()):,}")   # thousands separators
for lbl in v.set_labels:
    if lbl:
        lbl.set_fontsize(18)
        lbl.set_fontweight("bold")

ax.set_title("Egr1 Calling-Card peak overlap", fontsize=18, fontweight="bold", pad=28)
ax.text(0.5, 1.02, "peaks called separately per sex", transform=ax.transAxes,
        ha="center", fontsize=14, color="grey")
ax.text(0.5, -0.04, f"merged union regions, >=1 bp overlap  |  n = {len(udf):,} total",
        transform=ax.transAxes, ha="center", fontsize=14, color="grey")

plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(os.path.join(OUT_DIR, f"SupplFig4b_CC_peak_overlap_Venn_py.{ext}"),
                dpi=300, bbox_inches="tight")
print(f"\nWrote figure + 7 beds + summary to:\n  {OUT_DIR}")
