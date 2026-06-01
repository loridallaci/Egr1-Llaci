Egr1 Calling Card Peak Plot Pipeline

This repository contains a minimal pipeline and notebook to generate Venn diagrams of nearest-gene overlaps for Egr1 calling card peaks.

Included files:
- `calling_card_peak_pipeline.py`: CLI pipeline for annotation, filtering, gene extraction, and plotting.
- `cc_peak_plots.ipynb`: Minimal notebook to run the plot workflow and save a Venn diagram.

Quick start

1. Install dependencies (preferably in a virtual environment):

   pip install -r requirements.txt

2. Run notebook (recommended):

   jupyter notebook cc_peak_plots.ipynb

3. Or run the pipeline script directly:

   python calling_card_peak_pipeline.py --male-bed <male.bed> --female-bed <female.bed> --output-dir output/cc_peak_pipeline

Notes

- The pipeline expects input BED files formatted similarly to the provided Egr1 calling card peak files in the repository root.
- `pycallingcards` and `pyranges` are used for annotation and overlap analysis; ensure `bedtools` (used by pycallingcards) is available and set `--bedtools-path` if different.

License: Your choice. This repository contains your analysis code and a minimal notebook for reproducible plotting.
