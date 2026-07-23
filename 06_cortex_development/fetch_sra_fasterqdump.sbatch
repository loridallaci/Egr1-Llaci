#!/bin/bash
# =============================================================================
# CORRECT way to re-download the Zhu/Roussos cortex multiome raw data (PRJNA841952)
# for cellranger-arc reprocessing. Use THIS, not download_roussos_cortex_sra.sbatch
# (that one uses ENA's published FASTQ, which drops 2 of the 4 ATAC reads — verified
# 2026-07-20 via vdb-dump: ATAC spots = 50bp genomic / 8bp i7 / 16bp barcode / 49bp
# genomic, but ENA exports only the 50 + 16, losing the 2nd genomic mate that
# cellranger-arc requires). fasterq-dump --include-technical recovers all reads.
#
# ONE-TIME SETUP on HTCF (login node), from the folder you'll submit in:
#   BASE=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/roussos_cortex_sra; mkdir -p "$BASE/logs"; cd "$BASE"
#   # SRA toolkit (prebuilt, no module load):
#   wget -q "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz" -O sratk.tgz
#   mkdir -p sratoolkit && tar xzf sratk.tgz -C sratoolkit --strip-components=1
#   mkdir -p ~/.ncbi && printf '/LIBS/GUID = "%s"\n/repository/user/cache-disabled = "true"\n' "$(cat /proc/sys/kernel/random/uuid)" > ~/.ncbi/user-settings.mkfg
#   # runs list (SRR  sample  assay):
#   curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA841952&result=read_run&fields=run_accession,sample_title&format=tsv" \
#    | awk -F'\t' 'NR>1{a=($2~/ATAC/)?"atac":"gex"; s=$2; sub(/^sn(RNA|ATAC)seq_/,"",s); print $1"\t"s"\t"a}' | sort -k2,2 -k3,3 > runs.tsv
#   sbatch fetch_sra_fasterqdump.sbatch
#
# Output: fastq/<sample>/{gex,atac}/<sample>_<assay>_S1_L001_{R1,R2,R3,I1}_001.fastq.gz
#         (already named for cellranger-arc). NEXT STEP: cellranger-arc count per donor.
# =============================================================================
#SBATCH --job-name=roussos_fq
#SBATCH --array=1-24
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fq_%A_%a.out
#SBATCH --error=logs/fq_%A_%a.err
set -euo pipefail
B="$SLURM_SUBMIT_DIR"; export PATH="$B/sratoolkit/bin:$PATH"
read -r SRR SAMPLE ASSAY < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$B/runs.tsv")
OUT="$B/fastq/$SAMPLE/$ASSAY"; TMP="$B/tmp/$SRR"; SRADIR="$B/sra"
mkdir -p "$OUT" "$TMP" "$SRADIR"; echo "[$(date)] $SRR $SAMPLE/$ASSAY"
prefetch --max-size 100g -O "$SRADIR" "$SRR"
if [ "$ASSAY" = atac ]; then
  fasterq-dump --split-files --include-technical -e 6 -t "$TMP" -O "$OUT" "$SRADIR/$SRR/$SRR.sra"
  mv "$OUT/${SRR}_1.fastq" "$OUT/${SAMPLE}_atac_S1_L001_R1_001.fastq"   # 50bp genomic
  mv "$OUT/${SRR}_2.fastq" "$OUT/${SAMPLE}_atac_S1_L001_I1_001.fastq"   # 8bp i7 index
  mv "$OUT/${SRR}_3.fastq" "$OUT/${SAMPLE}_atac_S1_L001_R2_001.fastq"   # 16bp cell barcode
  mv "$OUT/${SRR}_4.fastq" "$OUT/${SAMPLE}_atac_S1_L001_R3_001.fastq"   # 49bp genomic
else
  fasterq-dump --split-files -e 6 -t "$TMP" -O "$OUT" "$SRADIR/$SRR/$SRR.sra"
  mv "$OUT/${SRR}_1.fastq" "$OUT/${SAMPLE}_gex_S1_L001_R1_001.fastq"    # 28bp CB+UMI
  mv "$OUT/${SRR}_2.fastq" "$OUT/${SAMPLE}_gex_S1_L001_R2_001.fastq"    # 90bp cDNA
fi
gzip -f "$OUT"/*.fastq; rm -rf "$TMP" "$SRADIR/$SRR"; echo "[$(date)] done $SRR"
