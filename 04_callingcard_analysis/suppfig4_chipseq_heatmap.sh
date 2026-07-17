#!/bin/bash
#SBATCH --job-name=cc_chip_heatmap
#SBATCH --output=cc_chip_heatmap_%j.log
#SBATCH --error=cc_chip_heatmap_%j.err
#SBATCH -c 8
#SBATCH --mem=60000
#SBATCH --time=8:00:00
# =============================================================================
# Supp Fig 4 - public EGR1 ChIP-seq signal over the FULL per-sex CC peak lists
#              (Fig-3a-style deepTools heatmap)
# =============================================================================
# Rows   = your two full CC peak lists (male 11,977 / female 11,772), the SAME
#          centered regions Fig 3a uses (centered_peaks_{males,females}.bed).
# Signal = public EGR1 ChIP-seq (Sun et al. 2019, Nat Commun; GSE67482;
#          male mouse prefrontal cortex; mm10), merged over its 2 replicates.
#          Your own CC CPM tracks are included as reference columns.
#
# Reads external EGR1 signal AT your peaks -> if it piles up at peak centers,
# the public ChIP corroborates the CC binding sites. (Peak-level overlap already
# shows 12-47x enrichment vs shuffled; this is the visual version.)
#
# CUT&RUN NOTE: GSE218312 (Yin 2025) deposits narrowPeak ONLY - no signal track -
# so it cannot be added here without reprocessing raw reads. See the companion
# script suppfig4_cutrun_signal_heatmap.sh for that pipeline.
# =============================================================================

set -euo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate deeptools_env
THREADS=${SLURM_CPUS_PER_TASK:-8}

DIR=/lts/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper
BW_DIR=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/cpm_bigwigs
CHIP_DIR=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/public_egr1_chipseq
OUT=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/suppfig4_chip
mkdir -p "$OUT"; cd "$DIR"

# ---- regions: the Fig 3a centered full per-sex peak lists ------------------
MALE_REG=centered_peaks_males.bed       # 11,977
FEMALE_REG=centered_peaks_females.bed   # 11,772
for f in "$MALE_REG" "$FEMALE_REG"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f (Fig3a region file)"; exit 1; }
done
echo "male regions:   $(wc -l < $MALE_REG)   (expect 11977)"
echo "female regions: $(wc -l < $FEMALE_REG) (expect 11772)"

# ---- signal tracks ---------------------------------------------------------
CHIP_BW=$CHIP_DIR/Egr1_chipseq_merged.bw    # made by chipseq_comparison_Egr1.sh STEP 5-6
CC_MALE=$BW_DIR/Male_Egr1_CPM.bw
CC_FEMALE=$BW_DIR/Female_Egr1_CPM.bw
for f in "$CHIP_BW" "$CC_MALE" "$CC_FEMALE"; do
    [[ -s "$f" ]] || { echo "ERROR: missing bigWig $f"; echo "  (ChIP bigWig comes from chipseq_comparison_Egr1.sh)"; exit 1; }
done

# ---- matrix: 2 region groups x 3 signal tracks, centered +/-1 kb -----------
computeMatrix reference-point --referencePoint center \
    -R "$MALE_REG" "$FEMALE_REG" \
    -S "$CC_MALE" "$CC_FEMALE" "$CHIP_BW" \
    -b 1000 -a 1000 --binSize 10 --missingDataAsZero -p "$THREADS" \
    -o "$OUT/CC_fullPeaks_vs_publicChIP_matrix.gz"

# ---- heatmap ---------------------------------------------------------------
plotHeatmap -m "$OUT/CC_fullPeaks_vs_publicChIP_matrix.gz" \
    -o "$OUT/SuppFig4_CC_fullPeaks_vs_publicChIP_heatmap.pdf" \
    --colorMap Blues Blues Greens \
    --refPointLabel center \
    --regionsLabel "Male CC peaks" "Female CC peaks" \
    --samplesLabel "Male Egr1 CC" "Female Egr1 CC" "Public Egr1 ChIP (Sun 2019)" \
    --xAxisLabel "distance from peak center (bp)" \
    --sortRegions descend --sortUsingSamples 1 \
    --heatmapHeight 12 --heatmapWidth 4 --dpi 300

# ---- same plot as PNG (for quick viewing) ----------------------------------
plotHeatmap -m "$OUT/CC_fullPeaks_vs_publicChIP_matrix.gz" \
    -o "$OUT/SuppFig4_CC_fullPeaks_vs_publicChIP_heatmap.png" \
    --plotFileFormat png \
    --colorMap Blues Blues Greens \
    --refPointLabel center \
    --regionsLabel "Male CC peaks" "Female CC peaks" \
    --samplesLabel "Male Egr1 CC" "Female Egr1 CC" "Public Egr1 ChIP (Sun 2019)" \
    --xAxisLabel "distance from peak center (bp)" \
    --sortRegions descend --sortUsingSamples 1 \
    --heatmapHeight 12 --heatmapWidth 4 --dpi 300

echo "Done -> $OUT/SuppFig4_CC_fullPeaks_vs_publicChIP_heatmap.pdf (+ .png)"
