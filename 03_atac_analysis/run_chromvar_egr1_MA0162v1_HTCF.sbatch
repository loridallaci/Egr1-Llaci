#!/bin/bash
#SBATCH --job-name=chromvar_egr1_v1
#SBATCH --output=chromvar_egr1_v1_%j.out
#SBATCH --error=chromvar_egr1_v1_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
# ---------------------------------------------------------------------
# Suppl Fig 1e - chromVAR EGR1-motif activity, OLD matrix MA0162.1 (HTCF).
# Submit from 03_atac_analysis/ :   sbatch run_chromvar_egr1_MA0162v1_HTCF.sbatch
#
# Same setup as run_chromvar_egr1_HTCF.sbatch, but runs the MA0162.1 variant.
# NOTE: no `module load` (lab rule). Make R available your usual way, e.g.
#   source ~/miniconda3/etc/profile.d/conda.sh && conda activate <your_r_env>
# or prepend your R install to PATH:  export PATH=/path/to/R/bin:$PATH
# The R env must have: Seurat, Signac, chromVAR, JASPAR2020, TFBSTools,
# BSgenome.Mmusculus.UCSC.mm10, ggplot2, dplyr, patchwork.
# ---------------------------------------------------------------------

set -euo pipefail

## >>> EDIT THESE TWO PATHS <<<
export MULTIOME_RDS="/lts/rmlab/rmlab_shared3/l.llaci/PATH/TO/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks.rds"
export OUTDIR="/lts/rmlab/rmlab_shared3/l.llaci/output/suppl_fig1"

## >>> make R available here (see note above) <<<
# source ~/miniconda3/etc/profile.d/conda.sh && conda activate r_multiome

Rscript Figures_SupplFig1e_chromVAR_Egr1_motif_MA0162v1_HTCF.R
