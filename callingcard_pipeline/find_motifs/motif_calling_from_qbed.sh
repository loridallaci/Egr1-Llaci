#!/bin/bash
#
#SBATCH --job-name=homer_qbed
#SBATCH --output=logs/homer_qbed.out
#SBATCH --error=logs/homer_qbed.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ihreiss@go.wustl.edu

#Assign these variables!
GENOME=mm10 #mm10 or hg38 only
#Change the name of your experiment directory (the directory that contains scripts, raw, and output_and_analysis)
PROJECT_DIR=/scratch/rmlab/rmlab_shared/ihreiss/Bmal1_CCs/gladiator/FINAL_GOOD_PREP/all_combined
PROJECT_INS=$PROJECT_DIR"/output_and_analysis/insertions"
PROJECT_OUT=$PROJECT_DIR"/output_and_analysis/motifs"
FILENAME=BMAL1_gladiator_2_10_2_6_9_combined.qbed
BACKGROUND=Unfused_gladiator_4_7_combined.qbed

# Create the names for the output files and make them variables (for easier typing for me)
FILE_STEM=`basename --suffix=.qbed $FILENAME`
HOMER_FILE=$PROJECT_OUT/$FILE_STEM"_qbed.homer"
HOMER_OUT=$PROJECT_OUT/"homer_qbed_"$FILE_STEM

BG_STEM=`basename --suffix=.qbed $BACKGROUND`
HOMER_BG=$PROJECT_OUT/$BG_STEM"_qbed.homer"

cat $PROJECT_INS/$FILENAME | grep -v "Start" | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"peak_"i++,$5}' | sort -k 1,1 -k2,2n |  awk -F"\t" -v OFS="\t" '{print $4,$1,$2,$3,$5}' > $HOMER_FILE
cat $PROJECT_INS/$BACKGROUND | grep -v "Start" | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"peak_"i++,$5}' | sort -k 1,1 -k2,2n |  awk -F"\t" -v OFS="\t" '{print $4,$1,$2,$3,$5}' > $HOMER_BG

module use /opt/htcf/modules
module use /opt/htcf/modules-legacy
module use /opt/apps/labs/rmlab/modules

module load homer

findMotifsGenome.pl $HOMER_FILE $GENOME $HOMER_OUT -bg $HOMER_BG -size 200 -preparsedDir preparsed -homer2 -p 4
