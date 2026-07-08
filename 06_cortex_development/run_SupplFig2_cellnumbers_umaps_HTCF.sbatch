#!/bin/bash
#SBATCH --job-name=suppfig2_cortex
#SBATCH --output=suppfig2_cortex_%j.out
#SBATCH --error=suppfig2_cortex_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
# ---------------------------------------------------------------------
# Suppl Fig 2 panels a/b - cortex development: cell numbers by stage x sex
# and per-stage RNA (RPCA) + ATAC (rLSI) UMAPs coloured by sex (HTCF).
#
# This is the canonical FULL run over all 5 stages. The local preview
# (CORTEX_LOCAL_PREVIEW=1) only covers 3 stages because the Adol/Adult
# *_withMotifs_V4 objects live on HTCF, not locally.
#
# Submit from 06_cortex_development/ :
#   sbatch run_SupplFig2_cellnumbers_umaps_HTCF.sbatch
#
# NOTE: no `module load` (lab rule). Make R available your usual way, e.g.
#   source ~/miniconda3/etc/profile.d/conda.sh && conda activate <your_r_env>
# or prepend your R install to PATH:
#   export PATH=/path/to/R/bin:$PATH
# The R env must have: Seurat (v5), Signac, ggplot2, patchwork.
# ---------------------------------------------------------------------

set -euo pipefail

# Directory holding the 5 *_withMotifs_V4_122825.rds stage objects on HTCF.
# This is already the script's built-in default; exported here to be explicit.
export DATA_DIR="/home/lllaci/data/cortex_development"

# Where the figures/CSV land on HTCF (the R script's default is a Windows path).
export OUT_DIR="${OUT_DIR:-$PWD/output}"
mkdir -p "$OUT_DIR"

# Leave CORTEX_LOCAL_PREVIEW unset so all 5 stages are built.
unset CORTEX_LOCAL_PREVIEW || true

## >>> make R available here (see note above) <<<
# source ~/miniconda3/etc/profile.d/conda.sh && conda activate r_multiome

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript not on PATH. Uncomment/edit the 'make R available' line" >&2
  echo "       above (e.g. conda activate <your_r_env>) and resubmit." >&2
  exit 127
fi

Rscript Figure_SupplFig2_cellnumbers_umaps.R
