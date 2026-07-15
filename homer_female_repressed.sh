#!/bin/bash
#SBATCH --job-name=homer_female_repressed
#SBATCH --output=homer_female_repressed_%j.log
#SBATCH --error=homer_female_repressed_%j.err
#SBATCH --ntasks=1
#SBATCH -c 4
#SBATCH --mem=60000
#SBATCH --time=24:00:00

#!/bin/bash
# =============================================================================
# HOMER motif analysis -- Egr1 CC peaks near female-repressed genes
# =============================================================================

# make HOMER tools visible (drop this line if findMotifsGenome.pl is already on PATH)
export PATH=/ref/rmlab/software/lori/homer/bin:$PATH

PEAKS=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper/Female_repressed_Egr1CC_peaks.bed
GENOME=/ref/rmlab/software/lori/homer/data/genomes/mm10
OUTDIR=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/Homer/Female_repressed_Egr1CC/

mkdir -p "$OUTDIR"

echo "Peaks : $PEAKS  ($(wc -l < "$PEAKS") regions)"
echo "Genome: $GENOME"
echo "Output: $OUTDIR"

findMotifsGenome.pl \
    "$PEAKS" \
    "$GENOME" \
    "$OUTDIR" \
    -size 1000 \
    -p 12

echo "Done -> ${OUTDIR}homerResults.html"