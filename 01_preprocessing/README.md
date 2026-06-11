# 01_preprocessing

Quality control and filtering of multiome data.

## Overview

This folder contains scripts for:
1. Running Cell Ranger ARC on raw FASTQ files
2. Aggregating multiple samples for joint analysis
3. Quality control of cells and features
4. Filtering and creating Seurat objects

---

## Scripts

### `cellranger_arc_count.sh`

SLURM batch script for running Cell Ranger ARC on HPC cluster to process individual samples.

**Requirements:**
- Cell Ranger ARC v2.0.0+
- Reference genome (mm10 or GRCh38)
- Libraries CSV file (see `libraries_template.csv`)

**Usage:**
```bash
# Update paths in script first
sbatch cellranger_arc_count.sh
```

**Outputs:**
- `outs/filtered_feature_bc_matrix.h5` - Filtered count matrix
- `outs/atac_fragments.tsv.gz` - ATAC fragment file
- `outs/atac_peak_annotation.tsv` - Peak annotations
- `outs/per_barcode_metrics.csv` - Cell-level QC metrics
- `outs/gex_molecule_info.h5` - RNA molecule information
- `outs/web_summary.html` - QC report

---

### `cellranger_arc_aggregate.sh`

SLURM batch script for aggregating multiple Cell Ranger ARC samples.

**Purpose:** Combine multiple samples (e.g., males and females, different replicates) for joint analysis.

**Requirements:**
- Cell Ranger ARC v2.0.0+
- All samples must be processed with the same reference genome
- Individual Cell Ranger ARC outputs must be complete

**Usage:**
```bash
# 1. Create aggregation CSV (see aggr_template.csv)
# 2. Update paths in script
# 3. Submit job
sbatch cellranger_arc_aggregate.sh
```

**Normalization:**
- We use `--normalize=none` because normalization will be performed in Seurat
- This preserves raw counts for proper downstream analysis

**Inputs (per sample):**
- `atac_fragments.tsv.gz` - ATAC fragment file
- `per_barcode_metrics.csv` - Cell-level metrics
- `gex_molecule_info.h5` - RNA molecule info

**Outputs:**
- `outs/filtered_feature_bc_matrix.h5` - Combined count matrix
- `outs/atac_fragments.tsv.gz` - Combined fragments
- `outs/web_summary.html` - Aggregation QC report

---

### `qc_filtering.R`

R script for quality control and cell filtering.

**Inputs:**
- Cell Ranger ARC output directories (aggregated or individual)

**QC Metrics:**

**RNA:**
- `nCount_RNA`: Total UMI counts per cell
- `nFeature_RNA`: Number of genes detected
- `percent.mt`: Mitochondrial gene percentage

**ATAC:**
- `nCount_ATAC`: Total ATAC fragments
- `nFeature_ATAC`: Number of accessible peaks
- `TSS.enrichment`: TSS enrichment score (>2 is good)
- `nucleosome_signal`: Nucleosome banding pattern (<2 is good)

**Filtering Thresholds:**
```r
# RNA
nFeature_RNA: 200 - 5000
nCount_RNA: 500 - 25000
percent.mt: < 15%

# ATAC
nCount_ATAC: > 1000
nFeature_ATAC: > 500
TSS.enrichment: > 2
nucleosome_signal: < 2
```

**Outputs:**
- Filtered Seurat object (`.rds` file)
- QC plots in `figures/qc/`

---

## Template Files

### `libraries_template.csv`
Template for Cell Ranger ARC count input. Format:
```csv
fastqs,sample,library_type
/path/to/RNA_fastqs,SampleName,Gene Expression
/path/to/ATAC_fastqs,SampleName,Chromatin Accessibility
```

### `aggr_template.csv`
Template for Cell Ranger ARC aggregation input. Format:
```csv
library_id,fragments,per_barcode_metrics,gex_molecule_info
Sample1,/path/to/sample1/outs/atac_fragments.tsv.gz,/path/to/sample1/outs/per_barcode_metrics.csv,/path/to/sample1/outs/gex_molecule_info.h5
Sample2,/path/to/sample2/outs/atac_fragments.tsv.gz,/path/to/sample2/outs/per_barcode_metrics.csv,/path/to/sample2/outs/gex_molecule_info.h5
```

---

## Samples in This Study

| Sample ID | Sex | Condition | Batch | Notes |
|-----------|-----|-----------|-------|-------|
| Sample_M | Male | GBM | 1 |
| Sample_F | Female | GBM | 2 |


---

## Run Order

### Complete Workflow:

1. **Process individual samples on HPC cluster:**
```bash
   # Run once per sample
   sbatch cellranger_arc_count.sh
```

2. **Aggregate all samples (optional but recommended):**
```bash
   # After all individual samples complete
   sbatch cellranger_arc_aggregate.sh
```

3. **QC and filtering locally in R:**
```r
   source("01_preprocessing/qc_filtering.R")
```

---

## Reference Genome

**Mouse (mm10):** `refdata-cellranger-arc-mm10-2020-A-2.0.0`  
**Human (GRCh38):** `refdata-cellranger-arc-GRCh38-2020-A-2.0.0`

Download from: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest

---

## Notes

- Cell Ranger ARC requires paired RNA and ATAC libraries from the same cells
- Minimum recommended sequencing depth:
  - RNA: 20,000 reads per cell
  - ATAC: 25,000 reads per cell
- All samples must use the same reference genome for aggregation
- Aggregation with `--normalize=none` is recommended for Seurat downstream analysis

---

## Files in This Folder
```
01_preprocessing/
├── README.md                          # This file
├── cellranger_arc_count.sh           # Process individual samples
├── cellranger_arc_aggregate.sh       # Aggregate multiple samples
├── qc_filtering.R                     # QC and filtering in R
├── libraries_template.csv             # Template for count input
├── aggr_template.csv                  # Template for aggregation input
└── DATA_AVAILABILITY.md               # Where to access raw data
```


