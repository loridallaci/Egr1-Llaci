#!/usr/bin/env python
# coding: utf-8

# This program takes as input an aligned BAM file of species-mixed Drop-seq reads.
# It then tags each read with the SRT barcode that was appended to read name by umi-tools 
# Prefix for SRT barcode = XQ

import argparse
import pysam
import sys

parser = argparse.ArgumentParser()
parser.add_argument("input", type = str)
parser.add_argument("output", type = str)
parser.add_argument("-t", "--tag", type = str, required = True)
args = parser.parse_args()

if __name__ == "__main__":
    # Open files
    bamfile = pysam.AlignmentFile(args.input, "rb")
    outfile = pysam.AlignmentFile(args.output, "wb", template = bamfile)
    tag = args.tag.split(':')
    for read in bamfile:
        umi = "".join(read.query_name.split('_')[-1:])
        if umi == 'G':
            srt_bc = "".join(read.query_name.split('_')[-2:-1]) ## Last 2 elements of read name are from UMI tools (3bp UMI and 1bp CellBarcode, separated by _. I merge them here as a 4bp srt barcode
            read.set_tag(tag[0], srt_bc, tag[1])
            outfile.write(read)
            continue
        else:
            pass
	
    # Close files
    outfile.close()
    bamfile.close()
