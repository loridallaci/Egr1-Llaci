#!/bin/bash
#SBATCH --job-name=cc_heatmap
#SBATCH --output=cc_heatmap_%j.log
#SBATCH --error=cc_heatmap_%j.err
#SBATCH -c 8
#SBATCH --mem=60000
#SBATCH --time=12:00:00
# =============================================================================
# Fig. 3a, step 2 — Egr1 Calling Card signal heatmap over peaks (deepTools)
# =============================================================================
# Plots CPM-normalized Egr1 CC signal (male, female) over Egr1-bound peaks,
# centered +/-1 kb. Peaks are grouped by the SEPARATELY-CALLED classification
# (shared / male-unique / female-unique) to match the Methods ("peaks were
# called on each sex separately") and panels d/e. Requires the CPM bigWigs from
# panel_a_1_make_cpm_bigwigs.sh.
#
# Tools: deepTools (computeMatrix, plotHeatmap) in the 'deeptools_env' conda env.
# =============================================================================

set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate deeptools_env
THREADS=${SLURM_CPUS_PER_TASK:-8}

DIR=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper
BW_DIR=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/cpm_bigwigs
cd "$DIR"

# ---- SEPARATELY-CALLED peak sets (canonical; match Methods + panels d/e) ----
SHARED=$DIR/shared_peaks_originalPeaks.bed
MALE=$DIR/male_unique_peaks_originalPeaks.bed
FEMALE=$DIR/female_unique_peaks_originalPeaks.bed

BW_MALE=$BW_DIR/Male_Egr1_CPM.bw
BW_FEMALE=$BW_DIR/Female_Egr1_CPM.bw

MATRIX=$DIR/Egr1CC_separateCall_matrix.gz
HEATMAP=$DIR/Figure3a_Egr1CC_separateCall_CPM.pdf

for f in "$SHARED" "$MALE" "$FEMALE" "$BW_MALE" "$BW_FEMALE"; do
    [[ -s "$f" ]] || { echo "ERROR: missing/empty input: $f" >&2; exit 1; }
done

# ---- signal matrix, centered +/-1 kb ---------------------------------------
computeMatrix reference-point \
    --referencePoint center \
    -b 1000 -a 1000 --binSize 10 \
    -R "$SHARED" "$MALE" "$FEMALE" \
    -S "$BW_MALE" "$BW_FEMALE" \
    --missingDataAsZero -p "$THREADS" \
    -o "$MATRIX"

# ---- heatmap (region-label order MUST match -R order) ----------------------
plotHeatmap \
    -m "$MATRIX" \
    -o "$HEATMAP" \
    --colorMap Blues RdPu \
    --refPointLabel center \
    --regionsLabel "Shared" "Male-unique" "Female-unique" \
    --samplesLabel "Male Egr1 CC" "Female Egr1 CC" \
    --xAxisLabel "distance from peak center (bp)" \
    --zMin 0 --zMax 5 --yMax 5 \
    --heatmapHeight 12

echo "Done -> $HEATMAP"
