#!/usr/bin/env python3
"""Simplified heatmap generation without pandas dependency."""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Ensure pyBigWig is available
try:
    import pyBigWig
    print("✓ pyBigWig available\n")
except ImportError:
    print("Installing pyBigWig...")
    os.system('pip install -q pyBigWig 2>/dev/null')
    import pyBigWig

os.chdir('/scratch/rmlab/rmlab_shared3/l.llaci/Egr1_paper/testing_CC_forPaper')

# ============================================================================
# Load BED files (simple format: chr start end)
# ============================================================================
def load_bed(filename):
    """Load BED file into list of (chr, start, end) tuples."""
    peaks = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                peaks.append((parts[0], int(parts[1]), int(parts[2])))
    return peaks

print("="*80)
print("LOADING PEAK FILES")
print("="*80)

shared_peaks = load_bed('combined_peaks.bed')
male_biased = load_bed('male_biased_overlaps.bed')
female_biased = load_bed('female_biased_overlaps.bed')

print(f"Shared peaks: {len(shared_peaks)}")
print(f"Male-biased peaks: {len(male_biased)}")
print(f"Female-biased peaks: {len(female_biased)}")
print(f"Total: {len(shared_peaks) + len(male_biased) + len(female_biased)}\n")

# ============================================================================
# Extract signal from bigWig files
# ============================================================================
def extract_signal(peaks, bw_file, window_bp=1000, bin_size=10):
    """Extract signal from bigWig for given peaks."""
    bw = pyBigWig.open(bw_file)
    n_bins = (window_bp * 2) // bin_size
    matrix = np.zeros((len(peaks), n_bins))
    
    print(f"Extracting from: {os.path.basename(bw_file)}")
    print(f"  Shape: {len(peaks)} peaks × {n_bins} bins")
    
    for i, (chrom, peak_start, peak_end) in enumerate(peaks):
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
            print(f"  Processed {i + 1}/{len(peaks)}")
    
    bw.close()
    return matrix

print("="*80)
print("EXTRACTING SIGNAL")
print("="*80)

print("\nMALE Egr1 CC:")
male_shared = extract_signal(shared_peaks, '../cpm_bigwigs/Male_Egr1_CPM.bw')
male_male_b = extract_signal(male_biased, '../cpm_bigwigs/Male_Egr1_CPM.bw')
male_female_b = extract_signal(female_biased, '../cpm_bigwigs/Male_Egr1_CPM.bw')
male_matrix = np.vstack([male_shared, male_male_b, male_female_b])

print("\nFEMALE Egr1 CC:")
female_shared = extract_signal(shared_peaks, '../cpm_bigwigs/Female_Egr1_CPM.bw')
female_male_b = extract_signal(male_biased, '../cpm_bigwigs/Female_Egr1_CPM.bw')
female_female_b = extract_signal(female_biased, '../cpm_bigwigs/Female_Egr1_CPM.bw')
female_matrix = np.vstack([female_shared, female_male_b, female_female_b])

print(f"\nMatrices created:")
print(f"  Male: {male_matrix.shape}")
print(f"  Female: {female_matrix.shape}\n")

# ============================================================================
# Calculate statistics and z-scale
# ============================================================================
print("="*80)
print("STATISTICS")
print("="*80)

male_p99 = np.percentile(male_matrix, 99)
female_p99 = np.percentile(female_matrix, 99)
z_max = np.ceil(max(male_p99, female_p99) * 2) / 2

print(f"Male 99th percentile: {male_p99:.4f}")
print(f"Female 99th percentile: {female_p99:.4f}")
print(f"Z-scale: 0 to {z_max}\n")

# ============================================================================
# Generate heatmaps
# ============================================================================
print("="*80)
print("GENERATING HEATMAPS")
print("="*80)

window_bp = 1000
n_bins = male_matrix.shape[1]
region_sizes = [len(shared_peaks), len(male_biased), len(female_biased)]
region_names = ['Shared', 'Male-biased', 'Female-biased']

def save_heatmap(matrix, filename, title, colormap):
    """Create and save a heatmap."""
    fig, ax = plt.subplots(figsize=(16, 14), dpi=150)
    
    im = ax.imshow(matrix, aspect='auto', cmap=colormap, vmin=0, vmax=z_max, interpolation='bilinear')
    
    # Region separators
    pos = region_sizes[0]
    ax.axhline(y=pos-0.5, color='white', linewidth=2, linestyle='--')
    pos += region_sizes[1]
    ax.axhline(y=pos-0.5, color='white', linewidth=2, linestyle='--')
    
    ax.set_xlabel('Distance from peak center (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Peaks', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # X-axis labels
    x_pos = np.linspace(0, n_bins-1, 5)
    x_labels = [str(int(-window_bp + i*(2*window_bp)/(n_bins-1))) for i in range(5)]
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    
    # Y-axis labels
    y_ticks = []
    y_labels = []
    pos = 0
    for name, size in zip(region_names, region_sizes):
        y_ticks.append(pos + size//2)
        y_labels.append(name)
        pos += size
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=11, fontweight='bold')
    
    plt.colorbar(im, ax=ax, label='Signal (CPM)')
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✓ {filename}")

save_heatmap(male_matrix, 'Male_Egr1CC_heatmap_FIXED.pdf', 'Male Egr1 CC - Heatmap', 'Blues')
save_heatmap(female_matrix, 'Female_Egr1CC_heatmap_FIXED.pdf', 'Female Egr1 CC - Heatmap', 'RdPu')

# ============================================================================
# Generate profile plots
# ============================================================================
print("\nGenerating profile plots...")

x_axis = np.linspace(-window_bp, window_bp, n_bins)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

fig, axes = plt.subplots(1, 2, figsize=(16, 5), dpi=150)

for ax_idx, (ax, matrix, title) in enumerate([(axes[0], male_matrix, 'Male'),
                                               (axes[1], female_matrix, 'Female')]):
    pos = 0
    for color, name, size in zip(colors, region_names, region_sizes):
        region = matrix[pos:pos+size, :]
        mean_signal = region.mean(axis=0)
        std_signal = region.std(axis=0)
        
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
plt.close()
print(f"✓ Egr1CC_profile_comparison_FIXED.pdf")

# ============================================================================
# Generate comprehensive figure
# ============================================================================
print("\nGenerating comprehensive figure...")

fig = plt.figure(figsize=(18, 12), dpi=150)
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

ax_ml = fig.add_subplot(gs[0, 0])
ax_fl = fig.add_subplot(gs[0, 1])
ax_mh = fig.add_subplot(gs[1, 0])
ax_fh = fig.add_subplot(gs[1, 1])

# Line plots
for ax, matrix, title in [(ax_ml, male_matrix, 'Male'),
                          (ax_fl, female_matrix, 'Female')]:
    pos = 0
    for color, name, size in zip(colors, region_names, region_sizes):
        region = matrix[pos:pos+size, :]
        mean_sig = region.mean(axis=0)
        std_sig = region.std(axis=0)
        ax.plot(x_axis, mean_sig, linewidth=2.5, label=name, color=color)
        ax.fill_between(x_axis, mean_sig - std_sig, mean_sig + std_sig, alpha=0.2, color=color)
        pos += size
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
    ax.set_xlabel('Distance from peak center (bp)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Signal (CPM)', fontsize=11, fontweight='bold')
    ax.set_title(f'{title} Egr1 CC', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)

# Heatmaps
for ax, matrix, cmap in [(ax_mh, male_matrix, 'Blues'),
                         (ax_fh, female_matrix, 'RdPu')]:
    im = ax.imshow(matrix, aspect='auto', cmap=cmap, vmin=0, vmax=z_max, interpolation='bilinear')
    pos = region_sizes[0]
    ax.axhline(y=pos-0.5, color='white', linewidth=2, linestyle='--')
    pos += region_sizes[1]
    ax.axhline(y=pos-0.5, color='white', linewidth=2, linestyle='--')
    
    ax.set_xlabel('Distance from peak center (bp)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Peaks', fontsize=10, fontweight='bold')
    
    x_pos = np.linspace(0, n_bins-1, 5)
    x_labels = [str(int(-window_bp + i*(2*window_bp)/(n_bins-1))) for i in range(5)]
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    
    y_ticks = []
    y_labels = []
    pos = 0
    for name, size in zip(region_names, region_sizes):
        y_ticks.append(pos + size//2)
        y_labels.append(name)
        pos += size
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=9, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Signal (CPM)')

fig.suptitle('Egr1 CC Sex-Biased Peak Analysis', fontsize=14, fontweight='bold', y=0.995)
plt.savefig('Egr1CC_comprehensive_analysis_FIXED.pdf', dpi=150, bbox_inches='tight')
plt.close()
print(f"✓ Egr1CC_comprehensive_analysis_FIXED.pdf")

# ============================================================================
# Verify outputs
# ============================================================================
print("\n" + "="*80)
print("VERIFICATION")
print("="*80 + "\n")

output_files = [
    'Male_Egr1CC_heatmap_FIXED.pdf',
    'Female_Egr1CC_heatmap_FIXED.pdf',
    'Egr1CC_profile_comparison_FIXED.pdf',
    'Egr1CC_comprehensive_analysis_FIXED.pdf'
]

success = True
for f in output_files:
    if os.path.exists(f):
        size = os.path.getsize(f) / (1024*1024)
        print(f"✓ {f} ({size:.2f} MB)")
    else:
        print(f"✗ {f} NOT FOUND")
        success = False

if success:
    print("\n✓ ALL OUTPUT FILES GENERATED SUCCESSFULLY!")
    sys.exit(0)
else:
    print("\n✗ Some files missing!")
    sys.exit(1)
