#!/bin/bash
#SBATCH --job-name=macs2_mm10gsize_check
#SBATCH --output=logs/macs2_mm10gsize_check_%j.out
#SBATCH --error=logs/macs2_mm10gsize_check_%j.err
#SBATCH --mem-per-cpu=50G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=l.llaci@wustl.edu
#SBATCH --export=NONE          # do NOT inherit the (base) conda env that broke spack

# --- Stop conda base from shadowing spack's python (the six.moves error) -----
conda deactivate 2>/dev/null || true

# --- Initialize the CORRECT spack (0.22), not the broken default 0.17.2 ------
# .bashrc is not sourced in a batch shell, so set spack up explicitly here.
source /ref/rmlab/software/spack-0.22/share/spack/setup-env.sh

# --- Activate the same env that called the original peaks --------------------
spack env activate r25
spack load r-rhtslib r-rsamtools libiconv libgit2

# Rscript (not `R < ... --save`) gives clean batch logging
Rscript macs2_aggr_lot6_mm10gsize_check.R
