#!/usr/bin/env python3
"""
Calling card peak analysis pipeline for Egr1 male/female CC peaks.

This script:
  1. loads male/female peak BED files
  2. annotates each peak set with pycallingcards
  3. filters peaks by nearest gene distance (default 20kb)
  4. extracts unique nearest and 2nd-nearest gene lists
  5. computes shared and unique peaks between male and female
  6. writes a Venn diagram for the nearest gene sets

Usage:
  python calling_card_peak_pipeline.py \
    --male-bed Egr1CC_peak_MaleEgr1_VS_MaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed \
    --female-bed Egr1CC_peak_FemaleEgr1_VS_FemaleWT_MACC2_window1000_YchromFiltered_window300_p05.bed \
    --output-dir output/cc_peak_pipeline \
    --bedtools-path /ref/rmlab/software/pycallingcards/bin
"""

import argparse
import logging
from pathlib import Path
import re

import pandas as pd

try:
    import pycallingcards as cc
except ImportError:
    cc = None

try:
    import pyranges as pr
except ImportError:
    pr = None

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

try:
    from matplotlib_venn import venn2
except ImportError:
    venn2 = None

PEAK_COLUMNS = [
    'Chr', 'Start', 'End', 'Center', 'Experiment Insertions', 'Background insertions',
    'Reference Insertions', 'pvalue Reference', 'pvalue Background', 'Fraction Experiment',
    'TPH Experiment', 'Fraction background', 'TPH background', 'TPH background subtracted',
    'pvalue_adj Reference'
]

DEFAULT_BEDTOOLS_PATH = '/ref/rmlab/software/pycallingcards/bin'
DEFAULT_REFERENCE = 'mm10'
DEFAULT_DISTANCE = 20000


def read_peak_bed(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', header=None, names=PEAK_COLUMNS, dtype={0: str})
    logging.info('Loaded %s with %d rows', path.name, len(df))
    return df


def find_annotation_file(prefix: str, cwd: Path) -> Path | None:
    candidates = list(cwd.glob(f'{prefix}*annotat*.*'))
    return candidates[0] if candidates else None


def annotate_peaks(prefix: str, peak_df: pd.DataFrame, reference: str, bedtools_path: str, outdir: Path, force: bool = False) -> pd.DataFrame:
    outdir.mkdir(parents=True, exist_ok=True)
    existing = find_annotation_file(prefix, outdir)
    if existing and not force:
        logging.info('Loading precomputed annotation: %s', existing.name)
        return pd.read_csv(existing, sep='\t', engine='python')

    if cc is None:
        raise RuntimeError('pycallingcards is not installed. Install it or provide precomputed annotation files.')

    logging.info('Annotating %s peaks using pycallingcards', prefix)
    ann = cc.pp.annotation(peak_df, reference=reference, bedtools_path=bedtools_path)
    ann2 = cc.pp.combine_annotation(peak_df, ann)
    out_path = outdir / f'{prefix}_annotated.tsv'
    ann2.to_csv(out_path, sep='\t', index=False)
    logging.info('Wrote annotation to %s', out_path)
    return ann2


def filter_by_distance(df: pd.DataFrame, column: str = 'Distance1', threshold: int = DEFAULT_DISTANCE) -> pd.DataFrame:
    if column not in df.columns:
        raise ValueError(f'Missing required column: {column}')
    filtered = df[df[column].abs() <= threshold].reset_index(drop=True)
    logging.info('Filtered to %d peaks within %d bp of nearest gene', len(filtered), threshold)
    return filtered


def detect_gene_columns(df: pd.DataFrame) -> tuple[str | None, str | None]:
    lower_cols = {c.lower(): c for c in df.columns}
    gene1 = None
    gene2 = None
    for name in df.columns:
        low = name.lower()
        if low.endswith('1') and 'gene' in low:
            gene1 = name
        if low.endswith('2') and 'gene' in low:
            gene2 = name
    if gene1 is None:
        for candidate in ['gene name1', 'gene1', 'gene name', 'gene']:
            if candidate in lower_cols:
                gene1 = lower_cols[candidate]
                break
    if gene2 is None and 'gene name2' in lower_cols:
        gene2 = lower_cols['gene name2']
    return gene1, gene2


def collect_gene_set(series: pd.Series) -> set[str]:
    genes = set()
    if series is None:
        return genes
    for value in series.dropna().astype(str):
        value = value.strip()
        if not value or value.lower() == 'nan':
            continue
        parts = re.split(r'[;,\|/]+', value)
        for part in parts:
            part = part.strip()
            if part:
                genes.add(part)
    return genes


def save_gene_list(genes: set[str], out_path: Path, header: str = 'gene') -> None:
    sorted_genes = sorted(genes)
    pd.Series(sorted_genes, name=header).to_csv(out_path, index=False)
    logging.info('Wrote %d genes to %s', len(sorted_genes), out_path)


def extract_unique_genes(df: pd.DataFrame, prefix: str, outdir: Path) -> tuple[set[str], set[str], Path, Path, Path]:
    gene1_col, gene2_col = detect_gene_columns(df)
    if gene1_col is None and gene2_col is None:
        raise ValueError('Could not detect nearest/2nd-nearest gene columns in annotation output')

    genes1 = collect_gene_set(df[gene1_col]) if gene1_col else set()
    genes2 = collect_gene_set(df[gene2_col]) if gene2_col else set()
    union_genes = genes1 | genes2

    outdir.mkdir(parents=True, exist_ok=True)
    g1_file = outdir / f'{prefix}_nearest_unique_genes.csv'
    g2_file = outdir / f'{prefix}_2ndnearest_unique_genes.csv'
    g_all_file = outdir / f'{prefix}_nearestplus2nd_unique_genes.csv'

    save_gene_list(genes1, g1_file)
    save_gene_list(genes2, g2_file)
    save_gene_list(union_genes, g_all_file)

    logging.info('%s: nearest=%d, 2nd=%d, union=%d', prefix, len(genes1), len(genes2), len(union_genes))
    return genes1, genes2, union_genes, g1_file, g2_file, g_all_file


def make_pyranges(df: pd.DataFrame) -> pr.PyRanges:
    if pr is None:
        raise RuntimeError('pyranges is required for overlap analysis')
    if 'Chromosome' not in df.columns:
        df = df.rename(columns={'Chr': 'Chromosome'})
    return pr.PyRanges(df[['Chromosome', 'Start', 'End']])


def save_bed(df: pd.DataFrame, path: Path) -> None:
    df[['Chromosome', 'Start', 'End']].to_csv(path, sep='\t', header=False, index=False)
    logging.info('Wrote BED file %s (%d intervals)', path.name, len(df))


def compute_peak_overlap(male_df: pd.DataFrame, female_df: pd.DataFrame, outdir: Path) -> tuple[Path, Path, Path, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    male_pr = make_pyranges(male_df)
    female_pr = make_pyranges(female_df)

    common = male_pr.intersect(female_pr)
    only_male = male_pr.subtract(female_pr)
    only_female = female_pr.subtract(male_pr)

    common_df = common.df
    male_unique_df = only_male.df
    female_unique_df = only_female.df

    common_bed = outdir / 'common_peaks_unsorted.bed'
    male_unique_bed = outdir / 'unique_peaks_male.bed'
    female_unique_bed = outdir / 'unique_peaks_female.bed'

    save_bed(common_df, common_bed)
    save_bed(male_unique_df, male_unique_bed)
    save_bed(female_unique_df, female_unique_bed)

    common_sorted = common_df.sort_values(['Chromosome', 'Start', 'End']).reset_index(drop=True)
    sorted_common_bed = outdir / 'common_peaks_sorted.bed'
    save_bed(common_sorted, sorted_common_bed)

    logging.info('Overlap counts: common=%d, male_unique=%d, female_unique=%d', len(common_df), len(male_unique_df), len(female_unique_df))
    return common_bed, male_unique_bed, female_unique_bed, common_df, male_unique_df, female_unique_df


def annotate_common_peaks(common_bed: Path, reference: str, bedtools_path: str, outdir: Path) -> pd.DataFrame:
    if cc is None:
        raise RuntimeError('pycallingcards is required for common peak annotation')
    if not common_bed.exists():
        raise FileNotFoundError(common_bed)
    common_df = pd.read_csv(common_bed, sep='\t', header=None, names=['Chromosome', 'Start', 'End'])
    common_peak_path = outdir / common_bed.name
    common_df.to_csv(common_peak_path, sep='\t', header=False, index=False)
    ann = cc.pp.annotation(peaks_path=str(common_peak_path), reference=reference, bedtools_path=bedtools_path)
    ann_path = outdir / 'shared_peaks_annotated.tsv'
    ann.to_csv(ann_path, sep='\t', index=False)
    logging.info('Annotated shared peaks to %s', ann_path)
    return ann


def plot_venn(set_a: set[str], set_b: set[str], labels: tuple[str, str], out_path: Path) -> None:
    if plt is None or venn2 is None:
        logging.warning('matplotlib and matplotlib-venn are required to draw Venn diagrams')
        return
    plt.figure(figsize=(6, 6))
    venn2([set_a, set_b], set_labels=labels)
    plt.title(f'{labels[0]} vs {labels[1]} unique nearest genes')
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    logging.info('Saved venn diagram to %s', out_path)


def main(args: argparse.Namespace) -> None:
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    male_bed = Path(args.male_bed)
    female_bed = Path(args.female_bed)

    male_df = read_peak_bed(male_bed)
    female_df = read_peak_bed(female_bed)

    male_ann = annotate_peaks('Male_Egr1CC', male_df, args.reference, args.bedtools_path, outdir, args.force_annotate)
    female_ann = annotate_peaks('Female_Egr1CC', female_df, args.reference, args.bedtools_path, outdir, args.force_annotate)

    male_20kb = filter_by_distance(male_ann, 'Distance1', args.distance_threshold)
    female_20kb = filter_by_distance(female_ann, 'Distance1', args.distance_threshold)

    male_20kb_path = outdir / 'Male_Egr1CC_peaks_20kbThreshhold.tsv'
    female_20kb_path = outdir / 'Female_Egr1CC_peaks_20kbThreshhold.tsv'
    male_20kb.to_csv(male_20kb_path, sep='\t', index=False)
    female_20kb.to_csv(female_20kb_path, sep='\t', index=False)
    logging.info('Saved filtered annotations: %s, %s', male_20kb_path.name, female_20kb_path.name)

    male_genes1, male_genes2, male_genes_union, male_g1, male_g2, male_gall = extract_unique_genes(male_20kb, 'Male_Egr1CC', outdir)
    female_genes1, female_genes2, female_genes_union, female_g1, female_g2, female_gall = extract_unique_genes(female_20kb, 'Female_Egr1CC', outdir)

    common_bed, male_unique_bed, female_unique_bed, common_df, male_unique_df, female_unique_df = compute_peak_overlap(male_20kb, female_20kb, outdir)

    if args.annotate_shared:
        annotate_common_peaks(common_bed, args.reference, args.bedtools_path, outdir)

    if args.plot_venn:
        plot_venn(male_genes1, female_genes1, ('Male nearest', 'Female nearest'), outdir / 'venn_nearest_male_vs_female.png')
        plot_venn(male_genes_union, female_genes_union, ('Male nearest+2nd', 'Female nearest+2nd'), outdir / 'venn_nearestplus2nd_male_vs_female.png')

    logging.info('Pipeline complete. Outputs are in %s', outdir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calling card peak analysis pipeline for Egr1 male/female data')
    parser.add_argument('--male-bed', required=True, help='Male peak BED file')
    parser.add_argument('--female-bed', required=True, help='Female peak BED file')
    parser.add_argument('--output-dir', default='output/cc_peak_pipeline', help='Directory to write results')
    parser.add_argument('--reference', default=DEFAULT_REFERENCE, help='Genome reference for annotation')
    parser.add_argument('--bedtools-path', default=DEFAULT_BEDTOOLS_PATH, help='Path to bedtools used by pycallingcards')
    parser.add_argument('--distance-threshold', type=int, default=DEFAULT_DISTANCE, help='Distance threshold for nearest gene filtering')
    parser.add_argument('--force-annotate', action='store_true', help='Force annotation instead of loading precomputed results')
    parser.add_argument('--plot-venn', action='store_true', help='Write Venn diagrams for nearest gene overlaps')
    parser.add_argument('--annotate-shared', action='store_true', help='Annotate shared male/female peaks after overlap')
    parser.add_argument('--verbose', action='store_true', help='Enable debug logging')
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')
    main(args)
