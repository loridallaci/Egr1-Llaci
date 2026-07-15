#!/usr/bin/env python3
import os
from collections import defaultdict

BASE_DIR = '/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper'
PUBLIC_CHIP = '/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/public_egr1_chipseq/Egr1_chipseq_peaks.bed'
PUBLIC_CUTRUN_EX = '/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/yin2024_egr1_cutrun/EGR1_ex_reproducible.bed'
PUBLIC_CUTRUN_IN = '/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/yin2024_egr1_cutrun/EGR1_in_reproducible.bed'
MALE_20KB = os.path.join(BASE_DIR, 'Male_Egr1CC_peaks_20kbThreshhold_091125_111225.txt')
FEMALE_20KB = os.path.join(BASE_DIR, 'Female_Egr1CC_peaks_20kbThreshhold_091125_111225.txt')

OUTPUT_SUMMARY = os.path.join(BASE_DIR, 'egr1_public_comparison_summary.txt')


def parse_annotated_20kb(path):
    intervals = []
    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')
        col_idx = {name: idx for idx, name in enumerate(header)}
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[col_idx['Chr']]
            start = int(parts[col_idx['Start']])
            end = int(parts[col_idx['End']])
            intervals.append((chrom, start, end, parts))
    return intervals


def parse_bed(path):
    intervals = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()[:3]
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            intervals.append((chrom, start, end, parts))
    return intervals


def group_by_chrom(intervals):
    grouped = defaultdict(list)
    for chrom, start, end, payload in intervals:
        grouped[chrom].append((start, end, payload))
    for chrom in grouped:
        grouped[chrom].sort(key=lambda x: (x[0], x[1]))
    return grouped


def intersect_intervals(a_intervals, b_intervals):
    """Count overlapping intervals and return overlap pairs."""
    overlaps = []
    i = 0
    j = 0
    while i < len(a_intervals) and j < len(b_intervals):
        a_start, a_end, a_payload = a_intervals[i]
        b_start, b_end, b_payload = b_intervals[j]
        if a_end <= b_start:
            i += 1
            continue
        if b_end <= a_start:
            j += 1
            continue
        # overlap exists
        overlaps.append((a_start, a_end, b_start, b_end, a_payload, b_payload))
        # advance the interval that ends first
        if a_end <= b_end:
            i += 1
        else:
            j += 1
    return overlaps


def compare_sets(set_name, set_intervals, public_intervals, public_name, output_prefix):
    stats = []
    combined = []
    for chrom in sorted(set_intervals.keys()):
        a_ints = set_intervals[chrom]
        b_ints = public_intervals.get(chrom, [])
        if not b_ints:
            continue
        overlaps = intersect_intervals(a_ints, b_ints)
        stats.append((chrom, len(a_ints), len(b_ints), len(overlaps)))
        combined.extend(overlaps)
    out_path = os.path.join(BASE_DIR, f'{output_prefix}.txt')
    with open(out_path, 'w') as out:
        out.write('#chrom	start	end	public_chrom	public_start	public_end\tset_region\tpublic_region\n')
        for a_start, a_end, b_start, b_end, a_payload, b_payload in combined:
            out.write(f'{a_payload[0]}\t{a_start}\t{a_end}\t{a_payload[0]}\t{b_start}\t{b_end}\t"{a_start}-{a_end}"\t"{b_start}-{b_end}"\n')
    return len(combined), out_path


def main():
    male20 = parse_annotated_20kb(MALE_20KB)
    female20 = parse_annotated_20kb(FEMALE_20KB)
    chip = parse_bed(PUBLIC_CHIP)
    cutrun_ex = parse_bed(PUBLIC_CUTRUN_EX)
    cutrun_in = parse_bed(PUBLIC_CUTRUN_IN)

    chip_group = group_by_chrom(chip)
    cutrun_group = group_by_chrom(cutrun_ex + cutrun_in)
    male_group = group_by_chrom(male20)
    female_group = group_by_chrom(female20)

    with open(OUTPUT_SUMMARY, 'w') as out:
        out.write('Dataset\tQuery_count\tPublic_count\tOverlap_count\tOverlap_pct_of_query\n')

        count, path = compare_sets('Male20kb', male_group, chip_group, 'Public_ChIP', 'Male20kb_vs_PublicEgr1ChIP_overlap')
        out.write(f'Male20kb_vs_PublicEgr1ChIP\t{len(male20)}\t{sum(len(v) for v in chip_group.values())}\t{count}\t{100*count/len(male20):.2f}\n')

        count, path = compare_sets('Female20kb', female_group, chip_group, 'Public_ChIP', 'Female20kb_vs_PublicEgr1ChIP_overlap')
        out.write(f'Female20kb_vs_PublicEgr1ChIP\t{len(female20)}\t{sum(len(v) for v in chip_group.values())}\t{count}\t{100*count/len(female20):.2f}\n')

        count, path = compare_sets('Male20kb', male_group, cutrun_group, 'Public_CUTRUN', 'Male20kb_vs_PublicEgr1CUTRUN_overlap')
        out.write(f'Male20kb_vs_PublicEgr1CUTRUN\t{len(male20)}\t{sum(len(v) for v in cutrun_group.values())}\t{count}\t{100*count/len(male20):.2f}\n')

        count, path = compare_sets('Female20kb', female_group, cutrun_group, 'Public_CUTRUN', 'Female20kb_vs_PublicEgr1CUTRUN_overlap')
        out.write(f'Female20kb_vs_PublicEgr1CUTRUN\t{len(female20)}\t{sum(len(v) for v in cutrun_group.values())}\t{count}\t{100*count/len(female20):.2f}\n')

    print('Saved summary to', OUTPUT_SUMMARY)
    print('Overlap files:')
    print('  Male20kb_vs_PublicEgr1ChIP_overlap.txt')
    print('  Female20kb_vs_PublicEgr1ChIP_overlap.txt')
    print('  Male20kb_vs_PublicEgr1CUTRUN_overlap.txt')
    print('  Female20kb_vs_PublicEgr1CUTRUN_overlap.txt')


if __name__ == '__main__':
    main()
