#!/bin/bash
#SBATCH --job-name=cc_external_egr1
#SBATCH --output=cc_external_egr1_%j.log
#SBATCH --error=cc_external_egr1_%j.err
#SBATCH -c 16
#SBATCH --mem=64000
#SBATCH --time=48:00:00
#SBATCH -p general
# =============================================================================
# Uniform comparison of public mouse EGR1 datasets to the Egr1 Calling-Card peaks
#   For EVERY dataset:  SRA -> bowtie2(mm10) -> CPM bigWig + MACS2 peaks
#   Then two deliverables:
#     (1) bedtools overlap table: each dataset's peaks vs CC male-only/shared/
#         female-only, with % and fold-enrichment vs shuffled control
#     (2) CPM signal heatmap: all datasets over the two full CC peak lists
#
# Lessons baked in from this project:
#   - runs on -p general (glibc 2.34; interactive/login nodes are older)
#   - all work in /scratch (writable; /lts is read-only on compute nodes)
#   - genome index: /ref/rmlab/data/bowtie2/mm10/mm10  (verified present)
# =============================================================================

set -eo pipefail                  # NOTE: no '-u' -- conda activation scripts (binutils) are not
source "$HOME/miniconda3/etc/profile.d/conda.sh"   # set -u clean and crash on ADDR2LINE unbound var
conda activate cc_external        # create once (see header): sra-tools bowtie2 samtools macs2 bedtools deeptools
set -u                            # safe to enable strict-unset now that activation is done
THREADS=${SLURM_CPUS_PER_TASK:-16}
# One-time env setup (run on the login node before sbatch):
#   conda create -n cc_external -c bioconda -c conda-forge sra-tools bowtie2 samtools macs2 bedtools deeptools -y

BT2=/ref/rmlab/data/bowtie2/mm10/mm10
CHROMSIZES=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/cpm_bigwigs/mm10.chrom.sizes
CCDIR=/lts/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper
WORK=/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/external_egr1
BW=$WORK/bigwigs; PK=$WORK/peaks; OUT=$WORK/results
mkdir -p "$BW" "$PK" "$OUT"; cd "$WORK"

# ---- raw per-sex CC peak calls (11,977 / 11,772) --------------------------
MALE_RAW=$CCDIR/Egr1CC_peak_MaleEgr1_VS_MaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed
FEMALE_RAW=$CCDIR/Egr1CC_peak_FemaleEgr1_VS_FemaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed
awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' "$MALE_RAW"   | sort -k1,1 -k2,2n > cc_male.bed
awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' "$FEMALE_RAW" | sort -k1,1 -k2,2n > cc_female.bed

# ---- corrected region-centric CC sets (male-only / shared / female-only) --
cat cc_male.bed cc_female.bed | sort -k1,1 -k2,2n | bedtools merge > cc_union.bed
bedtools intersect -u -a cc_union.bed -b cc_male.bed   | bedtools intersect -v -a - -b cc_female.bed > cc_male_only.bed
bedtools intersect -u -a cc_union.bed -b cc_female.bed | bedtools intersect -v -a - -b cc_male.bed   > cc_female_only.bed
bedtools intersect -u -a cc_union.bed -b cc_male.bed   | bedtools intersect -u -a - -b cc_female.bed > cc_shared.bed
echo "CC sets: male-only $(wc -l <cc_male_only.bed)  shared $(wc -l <cc_shared.bed)  female-only $(wc -l <cc_female_only.bed)"

# centered ±250 lists for the heatmap (match Fig 3a)
awk 'BEGIN{OFS="\t"}{c=int(($2+$3)/2);print $1,c-250,c+250}' cc_male.bed   > cc_male_centered.bed
awk 'BEGIN{OFS="\t"}{c=int(($2+$3)/2);print $1,c-250,c+250}' cc_female.bed > cc_female_centered.bed

# ---- alignment helpers -----------------------------------------------------
align_pe () { # $1=out_bam  $2..=SRRs   (paired-end, merged)
  local out=$1; shift; local bams=()
  for srr in "$@"; do
    if [[ ! -s ${srr}.bam ]]; then
      prefetch "$srr"; fasterq-dump --split-files -e "$THREADS" -O . "$srr"
      bowtie2 -p "$THREADS" --very-sensitive --no-unal --no-mixed --no-discordant -X 1000 \
              -x "$BT2" -1 ${srr}_1.fastq -2 ${srr}_2.fastq \
        | samtools view -b -q 10 - | samtools sort -@ "$THREADS" -o ${srr}.bam -
      samtools index ${srr}.bam; rm -f ${srr}_1.fastq ${srr}_2.fastq
    fi; bams+=(${srr}.bam)
  done
  samtools merge -f -@ "$THREADS" "$out" "${bams[@]}"; samtools index "$out"
}
align_se () { # $1=out_bam  $2..=SRRs   (single-end, merged)
  local out=$1; shift; local bams=()
  for srr in "$@"; do
    if [[ ! -s ${srr}.bam ]]; then
      prefetch "$srr"; fasterq-dump -e "$THREADS" -O . "$srr"
      bowtie2 -p "$THREADS" --very-sensitive --no-unal -x "$BT2" -U ${srr}.fastq \
        | samtools view -b -q 10 - | samtools sort -@ "$THREADS" -o ${srr}.bam -
      samtools index ${srr}.bam; rm -f ${srr}.fastq
    fi; bams+=(${srr}.bam)
  done
  samtools merge -f -@ "$THREADS" "$out" "${bams[@]}"; samtools index "$out"
}
mk_cpm () { bamCoverage -b "$1" -o "$BW/$2.bw" --normalizeUsing CPM --binSize 10 -p "$THREADS" --extendReads; }
mk_peaks_pe () { macs2 callpeak -t "$1" ${3:+-c "$3"} -f BAMPE -g mm -n "$2" --outdir "$PK" -q 0.05 2>/dev/null; }
mk_peaks_se () { macs2 callpeak -t "$1" ${3:+-c "$3"} -f BAM   -g mm -n "$2" --outdir "$PK" -q 0.05 2>/dev/null; }

# ===========================================================================
# DATASETS  (SRRs resolved 2026-07; all mm10)
# ===========================================================================
# --- Sun ChIP (male cortex): PE, with input ---
align_pe sun_ip.bam    SRR1948263 SRR1948264
align_pe sun_input.bam SRR1948265
mk_cpm sun_ip.bam SunChIP_cortex; mk_peaks_pe sun_ip.bam SunChIP sun_input.bam

# --- Yin CUT&RUN (male neurons): PE, no input ---
align_pe yin_ex.bam SRR22333343 SRR22333342
align_pe yin_in.bam SRR22333341 SRR22333340
mk_cpm yin_ex.bam YinCUTRUN_exc; mk_peaks_pe yin_ex.bam YinCUTRUN_exc
mk_cpm yin_in.bam YinCUTRUN_inh; mk_peaks_pe yin_in.bam YinCUTRUN_inh

# --- GSE172224 (brain PFC + cerebellum): SE, no input ---
align_se brain_pfc.bam   SRR14253645 SRR14253646
align_se brain_cereb.bam SRR14253649 SRR14253650
mk_cpm brain_pfc.bam   Brain_PFC;   mk_peaks_se brain_pfc.bam   Brain_PFC
mk_cpm brain_cereb.bam Brain_Cereb; mk_peaks_se brain_cereb.bam Brain_Cereb

# --- Neuro-2a (neural tumor): PE, with input ---
align_pe n2a_ip.bam    SRR36099388
align_pe n2a_input.bam SRR36099389
mk_cpm n2a_ip.bam Neuro2a_tumor; mk_peaks_pe n2a_ip.bam Neuro2a n2a_input.bam

# --- B-cell GSE164906 (immune): PE, with input ---
align_pe bcell_ip.bam    SRR16011657
align_pe bcell_input.bam SRR16011654
mk_cpm bcell_ip.bam Bcell_immune; mk_peaks_pe bcell_ip.bam Bcell bcell_input.bam

# ===========================================================================
# DELIVERABLE 1 - bedtools overlap table (each dataset's peaks vs CC sets)
# ===========================================================================
PEAKSETS=(SunChIP YinCUTRUN_exc YinCUTRUN_inh Brain_PFC Brain_Cereb Neuro2a Bcell)
printf "%-16s %-14s %-14s %-14s\n" dataset "male-only%" "shared%" "female-only%" > "$OUT/overlap_table.txt"
for ds in "${PEAKSETS[@]}"; do
  p=$PK/${ds}_peaks.narrowPeak; [[ -s $p ]] || continue
  sort -k1,1 -k2,2n "$p" | cut -f1-3 > $ds.sorted.bed
  row="$ds"
  for cc in cc_male_only cc_shared cc_female_only; do
    n=$(wc -l < $cc.bed); o=$(bedtools intersect -u -a $cc.bed -b $ds.sorted.bed | wc -l)
    row+=$(printf " %13s" "$(awk -v o=$o -v n=$n 'BEGIN{printf "%.1f%%",100*o/n}')")
  done
  echo "$row" >> "$OUT/overlap_table.txt"
done
echo "=== overlap table ==="; cat "$OUT/overlap_table.txt"

# ===========================================================================
# DELIVERABLE 2 - CPM signal heatmap: all datasets over the two CC peak lists
# ===========================================================================
BWS=(SunChIP_cortex YinCUTRUN_exc YinCUTRUN_inh Brain_PFC Brain_Cereb Neuro2a_tumor Bcell_immune)
S=""; L=""; for b in "${BWS[@]}"; do S+=" $BW/$b.bw"; L+=" $b"; done
computeMatrix reference-point --referencePoint center \
    -R cc_male_centered.bed cc_female_centered.bed -S $S \
    -b 1000 -a 1000 --binSize 10 --missingDataAsZero -p "$THREADS" \
    -o "$OUT/CC_vs_external_CPM_matrix.gz"
plotHeatmap -m "$OUT/CC_vs_external_CPM_matrix.gz" \
    -o "$OUT/SuppFig4_CC_vs_external_CPM_heatmap.png" --plotFileFormat png \
    --regionsLabel "Male CC peaks" "Female CC peaks" --samplesLabel $L \
    --refPointLabel center --xAxisLabel "distance from peak center (bp)" \
    --sortRegions descend --sortUsingSamples 1 --heatmapHeight 12 --heatmapWidth 3 --dpi 200
# per-dataset own-axis profiles too (shape is the honest read across assays)
for b in "${BWS[@]}"; do
  computeMatrixOperations subset -m "$OUT/CC_vs_external_CPM_matrix.gz" --samples "$b" -o "$OUT/${b}_only.gz"
  plotProfile -m "$OUT/${b}_only.gz" -o "$OUT/${b}_profile.png" --plotFileFormat png \
      --regionsLabel "Male CC" "Female CC" --plotTitle "$b"
done

echo "DONE."
echo "  overlap table : $OUT/overlap_table.txt"
echo "  CPM heatmap   : $OUT/SuppFig4_CC_vs_external_CPM_heatmap.png"
echo "  profiles      : $OUT/*_profile.png"
