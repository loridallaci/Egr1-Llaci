#!/bin/bash
#
#SBATCH --job-name=new_call_peaks_M_102824_window1000
#SBATCH --output=logs/new_call_peaks_M_102824_window1000.out
#SBATCH --error=logs/new_call_peaks_M_102824_window1000.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=l.llaci@wustl.edu

. /ref/rmlab/software/spack/share/spack/setup-env.sh
spack load python@3.8.12
#spack load py-pip@21.1.2/yhamiia
spack load bedtools2
#cd ~/.local/lib/python3.8/site-packages
PATH_DIR=~/.local/lib/python3.8/site-packages
PATH=$PATH_DIR:$PATH
HOMER_PATH=/home/ihreiss/homer/bin/
PATH=$HOMER_PATH:$PATH


# Assign these variables!
GENOME=mm10
PROJECT_DIR=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/Egr1_CC_Lori/101524_Egr1CC_analysis
#FOR DESKTOP USE: Set filenames for tf or unfused files
UNF=TCGTG_Male_WT_22CNMGLT4_S249_L008_R1_001_map_sort_final.qbed
TF=TCGTG_Male_Egr1_22CNMGLT4_S249_L008_R1_001_map_sort_final.qbed

PROJECT_INS=$PROJECT_DIR"/output_and_analysis/insertions"
PROJECT_OUT=$PROJECT_DIR"/output_and_analysis/peaks"

#set stem for file names
UNF_STEM=`basename --suffix=.qbed $UNF`
TF_STEM=`basename --suffix=.qbed $TF`

#call peaks
python3 call_peaks_window1000.py --unf_input $PROJECT_INS/$UNF --tf_input $PROJECT_INS/$TF --tf_stem $TF_STEM --unf_stem $UNF_STEM --out_stem $PROJECT_OUT --genome $GENOME  

#annotate peaks and go analysis
PEAK_CALLERS="cctools cccaller"
for peak_caller in ${PEAK_CALLERS}; do
	HOMER_FILE=$PROJECT_OUT/$TF_STEM"_"$UNF_STEM"_"$peak_caller"_peaks.homer"
	HOMER_BG=$PROJECT_OUT/$UNF_STEM"_"$peak_caller"_peak_bg.homer"
	ANNOTATE_FILE=$PROJECT_OUT/$TF_STEM"_"$UNF_STEM"_"$peak_caller"_peaks_annotate.txt"
	ANNOTATE_BG=$PROJECT_OUT/$UNF_STEM"_"$peak_caller"_BG_peaks_annotate.txt"
	GO_FILE_DIRECTORY=$PROJECT_OUT/$TF_STEM"_"$UNF_STEM"_"$peak_caller"_peaks_go_analysis"
	GO_BG_DIRECTORY=$PROJECT_OUT/$UNF_STEM"_"$peak_caller"_BG_peaks_go_analysis"
	HOMER_DIRECTORY=$PROJECT_OUT/$TF_STEM"_"$UNF_STEM"_"$peak_caller"_motif_calling"
	cat $PROJECT_OUT/$TF_STEM"_"$UNF_STEM"_"$peak_caller"_peaks.bed" | grep -v "Start" | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"peak_"i++,"."}' | sort -k 1,1 -k2,2n |  awk -F"\t" -v OFS="\t" '{print $4,$1,$2,$3,$5}' > $HOMER_FILE
	cat $PROJECT_OUT/$UNF_STEM"_"$peak_caller"_peak_bg.bed" | grep -v "Start" | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"peak_"i++,"."}' | sort -k 1,1 -k2,2n |  awk -F"\t" -v OFS="\t" '{print $4,$1,$2,$3,$5}' > $HOMER_BG
	annotatePeaks.pl  $HOMER_FILE $GENOME -go $GO_FILE_DIRECTORY > $ANNOTATE_FILE
	annotatePeaks.pl  $HOMER_BG $GENOME -go $GO_BG_DIRECTORY > $ANNOTATE_BG
	findMotifsGenome.pl $HOMER_FILE $GENOME $HOMER_DIRECTORY -size 200 -bg $HOMER_BG 
done
