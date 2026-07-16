#!/bin/bash
#SBATCH --job-name=chromvar_egr1
#SBATCH --output=chromvar_egr1_%j.out
#SBATCH --error=chromvar_egr1_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
# ---------------------------------------------------------------------
# Suppl Fig 1e - chromVAR EGR1-motif activity on the RPCA UMAP (HTCF).
# Submit from 03_atac_analysis/ :   sbatch run_chromvar_egr1_HTCF.sbatch
#
# NOTE: no `module load` (lab rule). Make R available your usual way, e.g.
#   source ~/miniconda3/etc/profile.d/conda.sh && conda activate <your_r_env>
# or prepend your R install to PATH:
#   export PATH=/path/to/R/bin:$PATH
# The R env must have: Seurat, Signac, chromVAR, JASPAR2020, TFBSTools,
# BSgenome.Mmusculus.UCSC.mm10, ggplot2, dplyr, patchwork.
# ---------------------------------------------------------------------

set -euo pipefail

## >>> EDIT THESE TWO PATHS <<<
# Processed lot6 multiome object (RNA + MACS2 `peaks` assays) on HTCF.
# Local copy is female_male_aggregated_081722_...filtered_012925_withPeaks.rds
export MULTIOME_RDS="/lts/rmlab/rmlab_shared3/l.llaci/PATH/TO/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks.rds"
export OUTDIR="/lts/rmlab/rmlab_shared3/l.llaci/output/suppl_fig1"
# Panel 1b (RNA UMAP) target - point at the repo's RNA multiome figure dir so
# 1b lands in 02_rna_analysis and 1e in 03_atac_analysis. Defaults to OUTDIR.
export RNA_OUTDIR="$OUTDIR"

## >>> make R available here (see note above) <<<
# source ~/miniconda3/etc/profile.d/conda.sh && conda activate r_multiome

Rscript Figures_SupplFig1e_chromVAR_Egr1_motif_HTCF.R
