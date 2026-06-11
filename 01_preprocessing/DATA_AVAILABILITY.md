# Data Availability

Raw sequencing data for this project is available at:

- **GEO Accession**: [GSE######] (or "Available upon publication")
- **SRA Project**: [PRJNA######] (if applicable)

## For Reproducibility

To reproduce this analysis:

1. Download raw FASTQ files from GEO/SRA
2. Update paths in `cellranger_arc_commands.sh`:
   - `REFERENCE`: Path to Cell Ranger ARC reference genome
   - `LIBRARIES`: Path to your libraries.csv file
   - `OUTPUT_DIR`: Where to save outputs
3. Run Cell Ranger ARC: `sbatch cellranger_arc_commands.sh`
4. Download Cell Ranger outputs and proceed to `qc_filtering.R`

## Required Files (not in repository)

Due to size constraints, the following files are not included in this repository:

- Raw FASTQ files (download from GEO/SRA)
- Cell Ranger ARC outputs (`.h5` files, fragment files)
- Processed Seurat objects (`.rds` files)

These can be regenerated using the scripts provided or obtained by contacting the authors.

## Contact

For access to processed data files, contact: [your email]