#!/bin/bash
#SBATCH --job-name=cc_extra_quicklook
#SBATCH --output=cc_extra_quicklook_%j.log
#SBATCH --error=cc_extra_quicklook_%j.err
#SBATCH -c 8
#SBATCH --mem=32000
#SBATCH --time=6:00:00
# =============================================================================
# QUICK LOOK: do other public mouse EGR1 datasets show signal at your CC peaks?
# Uses the GEO-deposited bigWigs directly (NO realignment) -> profiles over the
# two full CC peak lists, each on its OWN auto-scaled y-axis so weak-but-real
# central enrichment is visible (the trap that made the Sun ChIP look "flat").
#
# NORMALIZATION VARIES across these bigWigs (100K-downsampled / rpm / raw-sorted),
# so this is a SHAPE triage (center vs flank), not a magnitude comparison. Whatever
# shows a clean center bump, we later reprocess to CPM for the final figure.
#
# Datasets (all mm10, EGR1 genome-binding), ranked by GBM-relevance:
#   GSE172224 endogenous EGR1 MOWChIP, mouse PFC + cerebellum (normal brain) - BEST new comparator
#   GSE310510 FLAG-EGR1 ChIP, Neuro-2a neuroblastoma (neural TUMOR) - caveat: FLAG/overexpression
#   GSE164906 endogenous Egr1 ChIP, activated B cells (immune, very deep) - non-neural benchmark
# No mouse Egr1 dataset exists in glioma/GBM/astrocyte (verified negative).
# =============================================================================

set -euo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate deeptools_env
THREADS=${SLURM_CPUS_PER_TASK:-8}

BW=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/extra_egr1_bigwigs
OUT=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/suppfig4_extra
CC=/lts/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper
mkdir -p "$BW" "$OUT"; cd "$BW"

FTP=https://ftp.ncbi.nlm.nih.gov/geo/samples
get() { [[ -s "$2" ]] || wget -q -O "$2" "$1"; }   # download if missing

get "$FTP/GSM5244nnn/GSM5244362/suppl/GSM5244362_PFC-EGR1-100K-1.bw"        PFC_EGR1_rep1.bw
get "$FTP/GSM5244nnn/GSM5244363/suppl/GSM5244363_PFC-EGR1-100K-2.bw"        PFC_EGR1_rep2.bw
get "$FTP/GSM5244nnn/GSM5244366/suppl/GSM5244366_Cerebellum-EGR1-100K-1.bw" Cereb_EGR1_rep1.bw
get "$FTP/GSM5244nnn/GSM5244367/suppl/GSM5244367_Cerebellum-EGR1-100K-2.bw" Cereb_EGR1_rep2.bw
get "$FTP/GSM9302nnn/GSM9302417/suppl/GSM9302417_IPA.sorted.bw"             Neuro2a_FLAGEGR1.bw
get "$FTP/GSM5593nnn/GSM5593301/suppl/GSM5593301_aBEgr1.mm10.tag.rpm.bw"    Bcell_Egr1.bw
ls -la *.bw

MALE_REG=$CC/centered_peaks_males.bed       # 11,977
FEMALE_REG=$CC/centered_peaks_females.bed   # 11,772

# one matrix over both CC peak lists with all external tracks
computeMatrix reference-point --referencePoint center \
    -R "$MALE_REG" "$FEMALE_REG" \
    -S PFC_EGR1_rep1.bw PFC_EGR1_rep2.bw Cereb_EGR1_rep1.bw Cereb_EGR1_rep2.bw \
       Neuro2a_FLAGEGR1.bw Bcell_Egr1.bw \
    -b 1000 -a 1000 --binSize 10 --missingDataAsZero -p "$THREADS" \
    -o "$OUT/CC_vs_extra_matrix.gz"

# per-sample profile on its OWN axis (this is what reveals weak central signal)
for s in PFC_EGR1_rep1 PFC_EGR1_rep2 Cereb_EGR1_rep1 Cereb_EGR1_rep2 Neuro2a_FLAGEGR1 Bcell_Egr1; do
    computeMatrixOperations subset -m "$OUT/CC_vs_extra_matrix.gz" --samples "$s" -o "$OUT/${s}_only.gz"
    plotProfile -m "$OUT/${s}_only.gz" -o "$OUT/${s}_profile.png" --plotFileFormat png \
        --regionsLabel "Male CC" "Female CC" --plotTitle "$s at CC peaks"
done

echo "Done. Per-dataset profiles (own-axis) in $OUT/*_profile.png"
echo "Look for a symmetric bump at 'center' = that dataset's EGR1 binds your CC peaks."
