#!/bin/bash
#SBATCH --job-name=heatmap_CC_biased
#SBATCH --output=heatmap_CC_biased_%j.log
#SBATCH --error=heatmap_CC_biased_%j.err
#SBATCH --ntasks=1
#SBATCH -c 8
#SBATCH --mem=60000
#SBATCH --time=24:00:00
# =============================================================================
# deepTools heatmap -- Egr1 CC signal over MALE- and FEMALE-biased peaks only
# =============================================================================
set -euo pipefail

# --- activate the conda env --------------------------------------------------
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate deeptools_env

THREADS=${SLURM_CPUS_PER_TASK:-8}

# --- working directory: outputs land here ------------------------------------
DIR=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper
cd "$DIR"

# --- inputs (biased peaks only -- no shared) ---------------------------------
MALE=$DIR/Egr1_male_biased_peaks_calledtogether_logfc1.bed
FEMALE=$DIR/Egr1_female_biased_peaks_calledtogether_logfc1.bed

BW_MALE=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/cpm_bigwigs/Male_Egr1_CPM.bw
BW_FEMALE=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/cpm_bigwigs/Female_Egr1_CPM.bw

MATRIX=$DIR/Egr1CC_biased_only_matrix.gz
HEATMAP=$DIR/Egr1CC_biased_only_heatmap.pdf

# --- fail early with a clear message if any input is missing/empty -----------
for f in "$MALE" "$FEMALE" "$BW_MALE" "$BW_FEMALE"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty input: $f" >&2; exit 1; }
done
echo "All inputs found. Threads: $THREADS"

# --- compute the matrix ------------------------------------------------------
computeMatrix reference-point \
  --referencePoint center \
  -b 1000 -a 1000 --binSize 10 \
  -R "$MALE" "$FEMALE" \
  -S "$BW_MALE" "$BW_FEMALE" \
  --missingDataAsZero -p "$THREADS" \
  -o "$MATRIX"

# --- plot --------------------------------------------------------------------
# region-label order MUST match the -R order: male, female
plotHeatmap \
  -m Egr1CC_biased_only_matrix.gz \
  -o Egr1CC_biased_only_heatmap_autoscale.pdf \
  --colorMap Blues RdPu \
  --refPointLabel center \
  --regionsLabel "Male-biased" "Female-biased" \
  --samplesLabel "Male Egr1 CC" "Female Egr1 CC" \
  --xAxisLabel "distance from peak center (bp)" \
  --heatmapHeight 12

echo "Done -> $HEATMAP"
