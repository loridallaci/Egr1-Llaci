#!/usr/bin/env python3
"""
Generate Egr1 CC heatmaps from bigWig files and BED peak files.
No external dependencies required beyond pyBigWig, matplotlib, numpy.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Try to import pyBigWig
try:
    import pyBigWig
except ImportError:
    print("Installing pyBigWig...")
    os.system('pip install -q pyBigWig')
    import pyBigWig

os.chdir('/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper')
print(f"Working in: {os.getcwd()}\n")

# ============================================================================
# STEP 1: Load peak files
# ============================================================================
print("="*80)
print("STEP 1: Loading peak files")
print("="*80)

shared_peaks = pd.read_csv('combined_peaks.bed', sep='\t', header=None, names=['chr', 'start', 'end'])
male_biased = pd.read_csv('male_biased_overlaps.bed', sep='\t', header=None, names=['chr', 'start', 'end'])
female_biased = pd.read_csv('female_biased_overlaps.bed', sep='\t', header=None, names=['chr', 'start', 'end'])

print(f"Shared peaks: {len(shared_peaks)}")
print(f"Male-biased peaks: {len(male_biased)}")
print(f"Female-biased peaks: {len(female_biased)}")
print(f"Total: {len(shared_peaks) + len(male_biased) + len(female_biased)}\n")

# ============================================================================
# STEP 2: Define signal extraction function
# ============================================================================
def extract_signal_from_peaks(peak_df, bw_file, window_bp=1000, bin_size=10):
    """Extract signal from bigWig file for peaks."""
    bw = pyBigWig.open(bw_file)
    
    n_bins = (window_bp * 2) // bin_size
    matrix = np.zeros((len(peak_df), n_bins))
    
    print(f"Extracting from {os.path.basename(bw_file)}")
    print(f"  Matrix shape: {matrix.shape} ({len(peak_df)} peaks × {n_bins} bins)")
    
    for i, row in peak_df.iterrows():
        chrom = row['chr']
        peak_start = row['start']
        peak_end = row['end']
        peak_center = (peak_start + peak_end) // 2
        
        for j in range(n_bins):
            bin_start = peak_center - window_bp + (j * bin_size)
            bin_end = bin_start + bin_size
            
            try:
                stats = bw.stats(chrom, bin_start, bin_end, type='mean')
                matrix[i, j] = stats[0] if stats and stats[0] else 0.0
            except:
                matrix[i, j] = 0.0
        
        if (i + 1) % 500 == 0:
            print(f"    Processed {i + 1}/{len(peak_df)}")
    
    bw.close()
    return matrix

# ============================================================================
# STEP 3: Extract signal
# ============================================================================
print("="*80)
print("STEP 2: Extracting signal")
print("="*80)

print("\nMALE Egr1 CC:")
male_shared = extract_signal_from_peaks(shared_peaks, '../cpm_bigwigs/Male_Egr1_CPM.bw')
male_male_biased = extract_signal_from_peaks(male_biased, '../cpm_bigwigs/Male_Egr1_CPM.bw')
male_female_biased = extract_signal_from_peaks(female_biased, '../cpm_bigwigs/Male_Egr1_CPM.bw')
male_matrix = np.vstack([male_shared, male_male_biased, male_female_biased])

print("\nFEMALE Egr1 CC:")
female_shared = extract_signal_from_peaks(shared_peaks, '../cpm_bigwigs/Female_Egr1_CPM.bw')
female_male_biased = extract_signal_from_peaks(male_biased, '../cpm_bigwigs/Female_Egr1_CPM.bw')
female_female_biased = extract_signal_from_peaks(female_biased, '../cpm_bigwigs/Female_Egr1_CPM.bw')
female_matrix = np.vstack([female_shared, female_male_biased, female_female_biased])

print(f"\nCombined matrices:")
print(f"  Male: {male_matrix.shape}")
print(f"  Female: {female_matrix.shape}\n")

# ============================================================================
# STEP 4: Calculate statistics
# ============================================================================
print("="*80)
print("STEP 3: Data statistics")
print("="*80)

print(f"\nMale signal: min={male_matrix.min():.4f}, mean={male_matrix.mean():.4f}, "
      f"99th%ile={np.percentile(male_matrix, 99):.4f}, max={male_matrix.max():.4f}")

print(f"Female signal: min={female_matrix.min():.4f}, mean={female_matrix.mean():.4f}, "
      f"99th%ile={np.percentile(female_matrix, 99):.4f}, max={female_matrix.max():.4f}")

z_max = np.ceil(max(np.percentile(male_matrix, 99), 
                     np.percentile(female_matrix, 99)) * 2) / 2
print(f"\nUnified z-scale: 0 to {z_max}\n")

# ============================================================================
# STEP 5: Generate heatmaps
# ============================================================================
print("="*80)
print("STEP 4: Generating heatmaps")
print("="*80)

window_bp = 1000
n_bins = male_matrix.shape[1]
region_boundaries = [0, len(shared_peaks), len(shared_peaks) + len(male_biased)]

def plot_heatmap(matrix, title, colormap, filename, region_names, region_sizes):
    """Plot and save a heatmap."""
    fig, ax = plt.subplots(figsize=(16, 14), dpi=150)
    
    # Heatmap
    im = ax.imshow(matrix, aspect='auto', cmap=colormap, vmin=0, vmax=z_max, interpolation='bilinear')
    
    # Region separators
    for boundary in region_boundaries[1:-1]:
        ax.axhline(y=boundary - 0.5, color='white', linewidth=2, linestyle='--')
    
    # Labels
    ax.set_xlabel('Distance from peak center (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Peaks', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # X-axis
    x_positions = np.linspace(0, n_bins - 1, 5)
    x_labels = [f"{int(-(window_bp) + i * (2 * window_bp) / (n_bins - 1))}" for i in range(len(x_positions))]
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels)
    
    # Y-axis
    y_ticks = []
    y_labels = []
    pos = 0
    for name, size in zip(region_names, region_sizes):
        y_ticks.append(pos + size // 2)
        y_labels.append(name)
        pos += size
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=11, fontweight='bold')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, label='Signal (CPM)')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"✓ Saved: {filename}")
    plt.close()

region_names = ['Shared', 'Male-biased', 'Female-biased']
region_sizes = [len(shared_peaks), len(male_biased), len(female_biased)]

plot_heatmap(male_matrix, 'Male Egr1 CC - Heatmap', 'Blues', 
             'Male_Egr1CC_heatmap_FIXED.pdf', region_names, region_sizes)

plot_heatmap(female_matrix, 'Female Egr1 CC - Heatmap', 'RdPu', 
             'Female_Egr1CC_heatmap_FIXED.pdf', region_names, region_sizes)

# ============================================================================
# STEP 6: Generate profile plots
# ============================================================================
print("\nGenerating profile plots...")

x_axis = np.linspace(-window_bp, window_bp, n_bins)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

fig, axes = plt.subplots(1, 2, figsize=(16, 5), dpi=150)

for idx, (ax, matrix, title) in enumerate([(axes[0], male_matrix, 'Male'),
                                            (axes[1], female_matrix, 'Female')]):
    pos = 0
    for color, name, size in zip(colors, region_names, region_sizes):
        region_matrix = matrix[pos:pos+size, :]
        mean_signal = region_matrix.mean(axis=0)
        std_signal = region_matrix.std(axis=0)
        
        ax.plot(x_axis, mean_signal, linewidth=2, label=name, color=color)
        ax.fill_between(x_axis, mean_signal - std_signal, mean_signal + std_signal, 
                        alpha=0.2, color=color)
        pos += size
    
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Distance from peak center (bp)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Mean Signal (CPM)', fontsize=11, fontweight='bold')
    ax.set_title(f'{title} Egr1 CC - Profile Plot', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('Egr1CC_profile_comparison_FIXED.pdf', dpi=150, bbox_inches='tight')
print(f"✓ Saved: Egr1CC_profile_comparison_FIXED.pdf")
plt.close()

# ============================================================================
# STEP 7: Generate comprehensive figure
# ============================================================================
print("\nGenerating comprehensive figure...")

fig = plt.figure(figsize=(18, 12), dpi=150)
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

ax_male_line = fig.add_subplot(gs[0, 0])
ax_female_line = fig.add_subplot(gs[0, 1])
ax_male_heat = fig.add_subplot(gs[1, 0])
ax_female_heat = fig.add_subplot(gs[1, 1])

# Line plots
for ax, matrix, title in [(ax_male_line, male_matrix, 'Male'),
                          (ax_female_line, female_matrix, 'Female')]:
    pos = 0
    for color, name, size in zip(colors, region_names, region_sizes):
        region_matrix = matrix[pos:pos+size, :]
        mean_signal = region_matrix.mean(axis=0)
        std_signal = region_matrix.std(axis=0)
        ax.plot(x_axis, mean_signal, linewidth=2.5, label=name, color=color)
        ax.fill_between(x_axis, mean_signal - std_signal, mean_signal + std_signal, alpha=0.2, color=color)
        pos += size
    
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
    ax.set_xlabel('Distance from peak center (bp)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Signal (CPM)', fontsize=11, fontweight='bold')
    ax.set_title(f'{title} Egr1 CC', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)

# Heatmaps
for ax, matrix, cmap in [(ax_male_heat, male_matrix, 'Blues'),
                         (ax_female_heat, female_matrix, 'RdPu')]:
    im = ax.imshow(matrix, aspect='auto', cmap=cmap, vmin=0, vmax=z_max, interpolation='bilinear')
    for boundary in region_boundaries[1:-1]:
        ax.axhline(y=boundary - 0.5, color='white', linewidth=2, linestyle='--')
    
    ax.set_xlabel('Distance from peak center (bp)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Peaks', fontsize=10, fontweight='bold')
    
    x_positions = np.linspace(0, n_bins - 1, 5)
    x_labels = [f"{int(-(window_bp) + i * (2 * window_bp) / (n_bins - 1))}" for i in range(len(x_positions))]
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels)
    
    y_ticks = []
    y_labels = []
    pos = 0
    for name, size in zip(region_names, region_sizes):
        y_ticks.append(pos + size // 2)
        y_labels.append(name)
        pos += size
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=9, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Signal (CPM)')

fig.suptitle('Egr1 CC Sex-Biased Peak Analysis', fontsize=14, fontweight='bold', y=0.995)
plt.savefig('Egr1CC_comprehensive_analysis_FIXED.pdf', dpi=150, bbox_inches='tight')
print(f"✓ Saved: Egr1CC_comprehensive_analysis_FIXED.pdf")
plt.close()

# ============================================================================
# FINAL: Verify outputs
# ============================================================================
print("\n" + "="*80)
print("VERIFICATION")
print("="*80)

output_files = [
    'Male_Egr1CC_heatmap_FIXED.pdf',
    'Female_Egr1CC_heatmap_FIXED.pdf',
    'Egr1CC_profile_comparison_FIXED.pdf',
    'Egr1CC_comprehensive_analysis_FIXED.pdf'
]

all_exist = True
for f in output_files:
    if os.path.exists(f):
        size = os.path.getsize(f) / 1024 / 1024
        print(f"✓ {f} ({size:.2f} MB)")
    else:
        print(f"✗ {f} (NOT FOUND)")
        all_exist = False

if all_exist:
    print("\n✓ ALL OUTPUT FILES GENERATED SUCCESSFULLY!")
else:
    print("\n✗ Some files were not generated")
    sys.exit(1)
