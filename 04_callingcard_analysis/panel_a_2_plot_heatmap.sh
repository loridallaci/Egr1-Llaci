#!/bin/bash
#SBATCH --job-name=cc_heatmap_fig3a
#SBATCH --output=cc_heatmap_fig3a_%j.log
#SBATCH --error=cc_heatmap_fig3a_%j.err
#SBATCH -c 8
#SBATCH --mem=60000
#SBATCH --time=12:00:00
# =============================================================================
# Fig. 3a - Egr1 Calling Card CPM signal heatmap, ONE PANEL PER SEX
# =============================================================================
# Panel a of Llaci_Fig3 is TWO heatmaps side by side (Male | Female), each over
# that sex's OWN full peak set, centered +/-1 kb. It is NOT grouped into
# shared / male-unique / female-unique.
#
# ---------------------------------------------------------------------------
# 2026-07-17 -- REWRITTEN. The previous version of this script plotted a
# THREE-GROUP split (Shared / Male-unique / Female-unique) read from
# {shared,male_unique,female_unique}_peaks_originalPeaks.bed. That script did
# NOT reproduce the published panel a -- it generated a different figure, and
# those beds are the known-broken fragment sets (see CC_peak_overlap_venn.py).
#
# The real recipe was recovered from the working notes file
#   'Egr1 manuscript/Final Submission/CC/Egr1CC_CPM_Heatmap_Figure3'
# (final block, after the note "Let's keep it on each sex separately").
#
# VERIFIED on HTCF 2026-07-17:
#   wc -l /lts/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper/centered_peaks_{males,females}.bed
#     -> 11977 / 11772   == the FULL raw per-sex peak calls. Panel a is clean:
#        it uses ALL peaks, and is unaffected by the 20 kb filter / fragment bug.
#
# TRAP: the output filenames contain "originalPeaks". That string is copy-paste
# residue from earlier attempts in the notes file -- it does NOT describe the
# input. The input is the full peak set. Do not read the filename as provenance.
#
# NOT-FULLY-VERIFIED: the command that built centered_peaks_{males,females}.bed
# is absent from the notes file. STEP 1 below reconstructs it as center+/-250
# (500 bp), by analogy with the awk blocks used for the other bed sets in that
# same file, and because the counts match the raw calls exactly. Run STEP 0 to
# confirm against the beds already on LTS before trusting a regenerated set.
# Also note the notes file had a real mismatch (computeMatrix wrote
# ..._maleOnly.gz but plotHeatmap read centered_peaks_maleOnly.gz); the matrix
# names are made consistent here.
# =============================================================================

set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate deeptools_env
THREADS=${SLURM_CPUS_PER_TASK:-8}

DIR=/lts/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper
BW_DIR=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/cpm_bigwigs
cd "$DIR"

# The MACCs calls present on LTS (headerless BED; 11,977 / 11,772). Verified
# byte-identical to the *_p05_TTAAp001_091325_071526_p05.bed on the pipeline
# branch -- same content, two names. Male_Egr1CC_peaks_072825.txt is NOT on LTS.
RAW_MALE=Egr1CC_peak_MaleEgr1_VS_MaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed
RAW_FEMALE=Egr1CC_peak_FemaleEgr1_VS_FemaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed

BED_MALE=centered_peaks_males.bed
BED_FEMALE=centered_peaks_females.bed

BW_MALE=$BW_DIR/Male_Egr1_CPM.bw
BW_FEMALE=$BW_DIR/Female_Egr1_CPM.bw

# ---- STEP 0: verify the existing centered beds are the full peak sets -------
for b in "$BED_MALE" "$BED_FEMALE"; do
    [[ -s "$b" ]] || { echo "NOTE: $b absent -- STEP 1 will rebuild it." >&2; continue; }
    n=$(wc -l < "$b")
    w=$(awk 'NR==1{print $3-$2}' "$b")
    echo "$b: $n intervals, first width ${w}bp"
done
echo "Expected: 11977 / 11772 intervals, 500 bp wide."

# ---- STEP 1: (re)build centered beds -- ONLY if missing ---------------------
# center +/-250 -> 500 bp intervals, from the FULL raw peak calls (NOT the
# *_20kbThreshhold / *_originalPeaks sets).
# NB: the MACCs .bed files are HEADERLESS -- do NOT use NR>1 here, it would drop
# the first peak. (The *_072825.txt tables DO have a header; those are not on LTS.)
if [[ ! -s "$BED_MALE" ]]; then
    awk 'BEGIN{OFS="\t"} {c=int(($2+$3)/2); print $1, c-250, c+250}' "$RAW_MALE" > "$BED_MALE"
fi
if [[ ! -s "$BED_FEMALE" ]]; then
    awk 'BEGIN{OFS="\t"} {c=int(($2+$3)/2); print $1, c-250, c+250}' "$RAW_FEMALE" > "$BED_FEMALE"
fi

for f in "$BED_MALE" "$BED_FEMALE" "$BW_MALE" "$BW_FEMALE"; do
    [[ -s "$f" ]] || { echo "ERROR: missing/empty input: $f" >&2; exit 1; }
done

# ---- STEP 2: MALE panel ----------------------------------------------------
computeMatrix reference-point \
    --referencePoint center \
    -R "$BED_MALE" \
    -S "$BW_MALE" \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    --binSize 10 \
    -p "$THREADS" \
    -out Egr1CC_originalPeaks_CPM_matrix_maleOnly.gz

plotHeatmap \
    -m Egr1CC_originalPeaks_CPM_matrix_maleOnly.gz \
    -out Egr1CC_originalPeaks_CPM_matrix_maleOnly_Fig3.pdf \
    --colorMap "Blues" \
    --heatmapWidth 3 \
    --heatmapHeight 9

# ---- STEP 3: FEMALE panel --------------------------------------------------
computeMatrix reference-point \
    --referencePoint center \
    -R "$BED_FEMALE" \
    -S "$BW_FEMALE" \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    --binSize 10 \
    -p "$THREADS" \
    -out Egr1CC_originalPeaks_CPM_matrix_femaleOnly.gz

plotHeatmap \
    -m Egr1CC_originalPeaks_CPM_matrix_femaleOnly.gz \
    -out Egr1CC_originalPeaks_CPM_matrix_femaleOnly_Fig3.pdf \
    --colorMap "Blues" \
    --heatmapWidth 3 \
    --heatmapHeight 9

echo "Done -> Egr1CC_originalPeaks_CPM_matrix_{male,female}Only_Fig3.pdf"
echo "These two PDFs are assembled side-by-side as Fig 3a in Illustrator."
