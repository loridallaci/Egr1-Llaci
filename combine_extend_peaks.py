#!/usr/bin/env python3
"""
Combine male and female WT Egr1 peaks and intersect with sex-biased peaks
"""
import pandas as pd
import os
from pybedtools import BedTool

os.chdir('/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper')

# Peak file column names
peak_columns = ['Chr', 'Start', 'End', 'Center', 'Experiment Insertions', 'Background insertions', 
                'Reference Insertions', 'pvalue Reference', 'pvalue Background', 'Fraction Experiment', 
                'TPH Experiment', 'Fraction background', 'TPH background', 'TPH background subtracted', 
                'pvalue_adj Reference']

print("=" * 80)
print("LOADING PEAK FILES")
print("=" * 80)

# Load peak files
male_WT = pd.read_csv("Egr1CC_peak_MaleEgr1_VS_MaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed", 
                       sep='\t', header=None)
male_WT.columns = peak_columns
print(f"Male Egr1 vs Male WT: {len(male_WT)} peaks")

female_WT = pd.read_csv("Egr1CC_peak_FemaleEgr1_VS_FemaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed", 
                         sep='\t', header=None)
female_WT.columns = peak_columns
print(f"Female Egr1 vs Female WT: {len(female_WT)} peaks")

male_vs_female = pd.read_csv("Egr1CC_peak_MaleEgr1_VS_FemaleEgr1_MACC2_YchromFiltered_window300_bg005_p05_TTAAp001_052125.bed", 
                              sep='\t', header=None)
male_vs_female.columns = peak_columns
print(f"Male Egr1 vs Female Egr1: {len(male_vs_female)} peaks")

female_vs_male = pd.read_csv("Egr1CC_peak_FemaleEgr1_VS_MaleEgr1_MACC2_YchromFiltered_window300_bg005_p05_TTAAp001_052125.bed", 
                              sep='\t', header=None)
female_vs_male.columns = peak_columns
print(f"Female Egr1 vs Male Egr1: {len(female_vs_male)} peaks")

print("\n" + "=" * 80)
print("COMBINING MALE AND FEMALE WT PEAKS")
print("=" * 80)

# Combine male and female WT peaks
all_peaks = pd.concat([male_WT, female_WT], ignore_index=True)
print(f"Total peaks (combined): {len(all_peaks)}")

print("\n" + "=" * 80)
print("OVERLAPPING COMBINED PEAKS WITH SEX-BIASED PEAKS")
print("=" * 80)

# Create bedtools objects
combined_bed = BedTool.from_dataframe(all_peaks[['Chr', 'Start', 'End']])
male_vs_female_bed = BedTool.from_dataframe(male_vs_female[['Chr', 'Start', 'End']])
female_vs_male_bed = BedTool.from_dataframe(female_vs_male[['Chr', 'Start', 'End']])

# 1. Combined peaks vs Male Egr1 vs Female Egr1
print("\n1. Combined WT peaks vs Male Egr1 vs Female Egr1...")
overlap_vs_male_biased = combined_bed.intersect(male_vs_female_bed, wo=True)
overlap_vs_male_biased_df = overlap_vs_male_biased.to_dataframe()
overlap_vs_male_biased_df.columns = ['combined_chrom', 'combined_start', 'combined_end', 
                                      'male_chrom', 'male_start', 'male_end', 'overlap_bp']
print(f"   Overlapping regions: {len(overlap_vs_male_biased_df)}")
print(f"   {100*len(overlap_vs_male_biased_df)/len(all_peaks):.1f}% of combined peaks")
overlap_vs_male_biased_df.to_csv("Combined_WT_peaks_vs_MaleEgr1vsFemaleEgr1_overlaps.txt", 
                                  sep='\t', index=False, header=True)
print("   Saved: Combined_WT_peaks_vs_MaleEgr1vsFemaleEgr1_overlaps.txt")

# 2. Combined peaks vs Female Egr1 vs Male Egr1
print("\n2. Combined WT peaks vs Female Egr1 vs Male Egr1...")
overlap_vs_female_biased = combined_bed.intersect(female_vs_male_bed, wo=True)
overlap_vs_female_biased_df = overlap_vs_female_biased.to_dataframe()
overlap_vs_female_biased_df.columns = ['combined_chrom', 'combined_start', 'combined_end', 
                                        'female_chrom', 'female_start', 'female_end', 'overlap_bp']
print(f"   Overlapping regions: {len(overlap_vs_female_biased_df)}")
print(f"   {100*len(overlap_vs_female_biased_df)/len(all_peaks):.1f}% of combined peaks")
overlap_vs_female_biased_df.to_csv("Combined_WT_peaks_vs_FemaleEgr1vsMaleEgr1_overlaps.txt", 
                                    sep='\t', index=False, header=True)
print("   Saved: Combined_WT_peaks_vs_FemaleEgr1vsMaleEgr1_overlaps.txt")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"Total combined WT peaks: {len(all_peaks)}")
print(f"Peaks overlapping Male Egr1 > Female: {len(overlap_vs_male_biased_df)}")
print(f"Peaks overlapping Female Egr1 > Male: {len(overlap_vs_female_biased_df)}")
print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
