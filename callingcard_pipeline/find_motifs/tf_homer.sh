#!/bin/bash
#
#SBATCH --job-name=homer
#SBATCH --output=logs/homer_%a.out
#SBATCH --error=logs/homer_%a.err
#SBATCH --array=2-5
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ihreiss@go.wustl.edu

# Initialize settings
GENOME=mm10
PROJECT_DIR=/scratch/rmlab/rmlab_shared/ihreiss/Bmal1_CCs/gladiator/FINAL_GOOD_PREP/all_combined
PROJECT_PEAKS=$PROJECT_DIR"/output_and_analysis/peaks"
PROJECT_OUT=$PROJECT_DIR"/output_and_analysis/motifs"
SAMPLE_LINE=$( sed -n ${SLURM_ARRAY_TASK_ID}p homer_list.csv )

# Initialize important global variables
FILENAME=
BACKGROUND=

# Read local variables from $SAMPLE_LINE in subshell,
# then pass them along to the global versions
while IFS=, read -r filename background
do
	{ FILENAME=$filename; }
	{ BACKGROUND=$background; }
done <<< "$SAMPLE_LINE"

module use /opt/htcf/modules
module use /opt/htcf/modules-legacy
module use /opt/apps/labs/rmlab/modules

module load homer

# Create the names for the output files and make them variables (for easier typing for me)
FILE_STEM=`basename --suffix=.bed $FILENAME`
HOMER_FILE_NO_BG=$PROJECT_OUT/$FILE_STEM"_NO_BG.homer"
HOMER_OUT_NO_BG=$PROJECT_OUT/"homer_NO_BG_"$FILE_STEM

if [ -n "$BACKGROUND" ]
then
	echo "Background file is "$BACKGROUND
	BG_BASE=`basename $BACKGROUND`
	BG_STEM="${BG_BASE%%.*}"
	BG_HOMER_FILE=$PROJECT_OUT/"BG_"$BG_STEM"_FOR_FILE_"$FILE_STEM".homer"
	echo "BG_HOMER_FILE is "$BG_HOMER_FILE
	HOMER_FILE_WHEN_BG=$PROJECT_OUT/$FILE_STEM"_USING_BG_"$BG_STEM".homer"
	HOMER_OUT_WHEN_BG=$PROJECT_OUT/$FILE_STEM"_USING_BG_"$BG_STEM 
else
	echo "Background file not supplied. Will run homer without a background file."
fi

if [[ $FILENAME == *".bed" ]] && [[ $BACKGROUND == *".bed" ]]
then
	cat $PROJECT_PEAKS/$FILENAME | grep -v "Start" | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"peak_"i++,"."}' | sort -k 1,1 -k2,2n |  awk -F"\t" -v OFS="\t" '{print $4,$1,$2,$3,$5}' > $HOMER_FILE_WHEN_BG
	cat $PROJECT_PEAKS/$BACKGROUND | grep -v "Start" | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"peak_"i++,"."}' | sort -k 1,1 -k2,2n |  awk -F"\t" -v OFS="\t" '{print $4,$1,$2,$3,$5}' > $BG_HOMER_FILE
	findMotifsGenome.pl $HOMER_FILE_WHEN_BG $GENOME $HOMER_OUT_WHEN_BG -size 1000 -bg $BG_HOMER_FILE -preparsedDir preparsed -homer2 -p 4
elif [[ ! -n "$BACKGROUND" ]] && [[ $FILENAME == *".bed" ]]
then
	cat $PROJECT_PEAKS/$FILENAME | grep -v "Start" | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"peak_"i++,"."}' | sort -k 1,1 -k2,2n |  awk -F"\t" -v OFS="\t" '{print $4,$1,$2,$3,$5}' > $HOMER_FILE_NO_BG
	findMotifsGenome.pl $HOMER_FILE_NO_BG $GENOME $HOMER_OUT_NO_BG -size 1000 -preparsedDir preparsed -homer2 -p 4
else
	echo "File name is "$FILENAME"."
	echo "File format not recognized. Please submit either a .sig or .bed output file from blockify or ccftools peak calling."
fi




