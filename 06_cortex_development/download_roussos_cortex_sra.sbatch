#!/bin/bash
# =============================================================================
# Parallel download of the Zhu/Roussos developing-cortex multiome FASTQ
# (PRJNA841952 / GSE204684) from ENA — one SLURM array task per SRA run, so all
# 24 runs (12 donors x GEX+ATAC) download concurrently.
#
# Reads sra_download_manifest.tsv (columns: SRR  sample  assay  URL_R1  URL_R2  MD5_R1  MD5_R2).
# Files land in  $OUTDIR/<sample>/<assay>/  e.g.  roussos_cortex_fastq/LaFet1/gex/SRR19374512_1.fastq.gz
#
# Submit (from the folder holding this script + the manifest):
#   mkdir -p logs
#   sbatch download_roussos_cortex_sra.sbatch
# (OUTDIR defaults to the path below; override with --export=ALL,OUTDIR=/some/other/path)
#
# Total download ≈ 200 GB. No `module load` (not available on this HPC); uses wget only.
# =============================================================================
#SBATCH --job-name=roussos_dl
#SBATCH --array=1-24          # all 24 runs at once; append %8 (i.e. 1-24%8) to cap concurrency if the network/quota complains
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --output=logs/dl_%A_%a.out
#SBATCH --error=logs/dl_%A_%a.err

set -euo pipefail

MANIFEST="${MANIFEST:-sra_download_manifest.tsv}"
OUTDIR="${OUTDIR:-/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/roussos_cortex_fastq}"

# pull this task's row from the manifest
IFS=$'\t' read -r SRR SAMPLE ASSAY URL1 URL2 MD5_1 MD5_2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
DEST="$OUTDIR/$SAMPLE/$ASSAY"
mkdir -p "$DEST"
echo "[$(date)] task ${SLURM_ARRAY_TASK_ID}: ${SRR}  ${SAMPLE}/${ASSAY}  ->  ${DEST}"

download_one() {   # url  expected_md5
  local url="$1" md5="$2" f; f="$(basename "$url")"
  cd "$DEST"
  if [[ -f "$f" ]] && [[ "$(md5sum "$f" | awk '{print $1}')" == "$md5" ]]; then
    echo "  already complete: $f"; return 0
  fi
  echo "  downloading $f"
  wget -c -q --tries=10 --timeout=120 -O "$f" "$url"
  local got; got="$(md5sum "$f" | awk '{print $1}')"
  if [[ "$got" != "$md5" ]]; then
    echo "  MD5 MISMATCH for $f (got $got, expected $md5)" >&2; exit 1
  fi
  echo "  ok (md5 verified): $f"
}

download_one "$URL1" "$MD5_1"
download_one "$URL2" "$MD5_2"
echo "[$(date)] done ${SRR}"
