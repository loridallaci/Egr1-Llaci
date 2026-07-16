#!/bin/bash
#
#SBATCH --job-name=unf_macs_call_peaks
#SBATCH --output=logs/unf_macs_call_peaks.out
#SBATCH --error=logs/unf_macs_call_peaks.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ihreiss@go.wustl.edu

. /ref/rmlab/software/spack/share/spack/setup-env.sh
spack load python@3.8.12
#spack load py-pip@21.1.2/yhamiia
spack load bedtools2
#cd ~/.local/lib/python3.8/site-packages
PATH_DIR=~/.local/lib/python3.8/site-packages 
PATH=$PATH_DIR:$PATH
HOMER_PATH=/home/ihreiss/homer/bin/
PATH=$HOMER_PATH:$PATH

#python3 -V
#which python3
# Assign these variables!
GENOME=mm10
PROJECT_DIR=/scratch/rmlab/rmlab_shared/ihreiss/Bmal1_CCs/032924_macs/041624_spikein
#FOR DESKTOP USE: Set filenames for tf or unfused files
UNF_BOUND=whole_brain_acsa2_bound.qbed
UNF_FT=whole_brain_acsa2_FT.qbed

PROJECT_INS=$PROJECT_DIR"/output_and_analysis/insertions"
PROJECT_OUT=$PROJECT_DIR"/output_and_analysis/peaks"

#set stem for file names
UNF_BOUND_STEM=`basename --suffix=.qbed $UNF_BOUND`
UNF_FT_STEM=`basename --suffix=.qbed $UNF_FT`

#call peaks
python3 $PROJECT_DIR/scripts/call_peaks/unf_macs_call_peaks.py --unf_bound_input $PROJECT_INS/$UNF_BOUND --unf_FT_input $PROJECT_INS/$UNF_FT --unf_bound_stem $UNF_BOUND_STEM --unf_FT_stem $UNF_FT_STEM --out_stem $PROJECT_OUT --genome $GENOME --homer_path $HOMER_PATH


