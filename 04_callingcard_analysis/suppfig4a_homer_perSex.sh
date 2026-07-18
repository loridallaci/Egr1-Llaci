#!/bin/bash
#SBATCH --job-name=homer_cc_persex
#SBATCH --output=homer_cc_persex_%j.log
#SBATCH --error=homer_cc_persex_%j.err
#SBATCH -c 12
#SBATCH --mem=32000
#SBATCH --time=12:00:00
# =============================================================================
# Supp Fig 4a - HOMER known-motif enrichment on the FULL per-sex Egr1 CC peaks
# =============================================================================
# Runs HOMER separately on ALL male peaks (11,977) and ALL female peaks (11,772)
# -- NOT the unique/shared split, NOT the 20 kb-gene-filtered subset. Purpose:
# show that the Egr1 motif is a top-enriched known motif in each sex's binding.
#
# No prior HOMER run used these exact full sets (verified 2026-07-17): existing
# runs were on the window-1000 call (2,842 / 3,803), the unique/shared 20 kb sets
# (7,357 / 7,495 / 3,187), or the 20 kb-filtered full set (9,457). This fills that
# gap on the standardized window-300 calls.
#
# Read Egr1's rank + p-value from knownResults.txt (the "(of N)" in the column-6
# header records the peak count, so provenance is self-documenting).
# =============================================================================

set -euo pipefail

# HOMER needs to be on PATH. Do NOT `module load` (not available on this HPC);
# activate the env / prepend the bin dir that already has findMotifsGenome.pl.
export PATH=/ref/rmlab/software/lori/homer/bin:$PATH
GENOME=/ref/rmlab/software/lori/homer/data/genomes/mm10   # or just "mm10" if preinstalled

DIR=/lts/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper
OUT=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/Homer
cd "$DIR"
mkdir -p "$OUT"

MALE_RAW=Egr1CC_peak_MaleEgr1_VS_MaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed      # 11,977
FEMALE_RAW=Egr1CC_peak_FemaleEgr1_VS_FemaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed # 11,772

# ---- make HOMER-NATIVE peak files: id / chr / start / end / strand ----------
# MUST be 5 columns in HOMER-native order (id FIRST, explicit strand). A 3- or
# 4-column BED FAILS: HOMER's peakfile parser needs >=5 columns and its BED
# auto-detect does not fire on a 4-col file -> it silently drops ALL peaks and
# reports "0 total sequences read" (verified 2026-07-17). Do NOT use a .bed here.
awk 'BEGIN{OFS="\t"} {print "m"NR, $1, $2, $3, "+"}' "$MALE_RAW"   > male_homer.txt
awk 'BEGIN{OFS="\t"} {print "f"NR, $1, $2, $3, "+"}' "$FEMALE_RAW" > female_homer.txt
echo "male peaks:   $(wc -l < male_homer.txt)   (expect 11977)"
echo "female peaks: $(wc -l < female_homer.txt) (expect 11772)"

# ---- HOMER known-motif enrichment (-size 1000; matches the window-1000 CC runs)
# We only need KNOWN-motif enrichment (Egr1 rank, from knownResults.txt). The de
# novo step ("De novo motif finding") may error with "Filtered out all motifs" in
# this HOMER install -- that does NOT affect knownResults, which is what we report.
# -preparsedDir is not needed: mm10.1000.* background already exists and is readable.
findMotifsGenome.pl male_homer.txt   "$GENOME" "$OUT/SuppFig4a_Male_allPeaks_Egr1"   -size 1000 -p 12
findMotifsGenome.pl female_homer.txt "$GENOME" "$OUT/SuppFig4a_Female_allPeaks_Egr1" -size 1000 -p 12

echo
echo "=== Egr1 rank in each (from knownResults.txt) ==="
for s in Male Female; do
    f="$OUT/SuppFig4a_${s}_allPeaks_Egr1/knownResults.txt"
    n=$(head -1 "$f" | grep -oE "of [0-9]+" | head -1 | grep -oE "[0-9]+")
    echo "--- $s (N=$n)"
    awk -F'\t' 'NR>1 && $1 ~ /^Egr/ {print "   #"NR-1"  "$1"  p="$3"  tgt="$7}' "$f" | head -3
done

echo
echo "Done. Supp Fig 4a inputs:"
echo "  $OUT/SuppFig4a_Male_allPeaks_Egr1/knownResults.{txt,html}"
echo "  $OUT/SuppFig4a_Female_allPeaks_Egr1/knownResults.{txt,html}"
