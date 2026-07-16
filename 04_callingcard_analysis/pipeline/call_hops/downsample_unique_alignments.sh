#!/bin/bash
# 
#SBATCH --job-name=downsamp_bulkRNACallingCardsBarcodes
#SBATCH --output=logs/downsamp_bulkRNACallingCardsBarcodes.out
#SBATCH --error=logs/downsamp_bulkRNACallingCardsBarcodes.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ihreiss@go.wustl.edu

#Assign these variables!
GENOME=mm10 #mm10 or hg38 only
#Change the name of your experiment directory (the directory that contains scripts, raw, and output_and_analysis)
PROJECT_DIR=/scratch/rmlab/rmlab_shared/ihreiss/Bmal1_CCs/gladiator_redo/spike_in

eval $(spack load --sh samtools@1.13 )

FILENAME=
P1=
P2=
P3=
P4=
P5=

#END OF SECTION TO UN-COMMENT FOR HTCF CLUSTER USE
#IF YOU ARE USING HTCF REMEMBER TO MAKE A LOGS FOLDER!

#Set up stem for naming output files.
BASE=`basename $FILENAME`
STEM="${BASE%%.*}"

PROJECT_OUT=$PROJECT_DIR"/output_and_analysis"
PROJECT_DOWNSAMP=$PROJECT_OUT"/downsample"
DOWNSAMP_STEM=$PROJECT_DOWNSAMP/$STEM

for ((i=1; i<6; i++))
do
    pvar="P$i"
    echo "${!pvar}"
    samtools view -b -s ${!pvar} $PROJECT_OUT/$FILENAME > $DOWNSAMP_STEM"_"${!pvar}".bam"
done

