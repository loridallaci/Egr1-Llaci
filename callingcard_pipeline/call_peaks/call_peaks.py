#!/usr/bin/env python
import pycallingcards as cc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--unf_input', type = str)
parser.add_argument('--tf_input', type = str)
parser.add_argument('--tf_stem', type = str)
parser.add_argument('--unf_stem', type = str)
parser.add_argument('--out_stem', type = str)
parser.add_argument('--genome', type = str)
args = parser.parse_args()

cctools_pvalue_cutoffbg = 0.001 #lower means needs more significance to pass
cctools_window_size = 1500
cctools_step_size = 500
cctools_pvalue_cutoffTTAA = 0.001
cctools_lam_win_size = None
cctools_pseudocounts = 0.4

cccaller_pvalue_cutoffbg = 0.001
cccaller_maxbetween = 1500
cccaller_pvalue_cutoffTTAA = 0.001
ccaller_lam_win_size = None
cccaller_pseudocounts = 0.4

TF = cc.rd.read_qbed(args.tf_input)

bg = cc.rd.read_qbed(args.unf_input)

cctools_peaks = args.out_stem + "/" + args.tf_stem + "_" + args.unf_stem + "_cctools_peaks.bed"
cctools_peak_bg_for_homer = args.out_stem + "/" + args.unf_stem + "_cctools_peak_bg.bed"
cccaller_peaks = args.out_stem + "/" + args.tf_stem + "_" + args.unf_stem + "_cccaller_peaks.bed"
cccaller_peak_bg_for_homer = args.out_stem + "/" + args.unf_stem + "_cccaller_peak_bg.bed"

peak_data_cctools = cc.pp.call_peaks(TF, bg,  method = "MACCs", reference = args.genome, pvalue_cutoffbg = cctools_pvalue_cutoffbg, window_size=cctools_window_size, step_size = cctools_step_size, pvalue_cutoffTTAA = cctools_pvalue_cutoffTTAA, lam_win_size = cctools_lam_win_size,  pseudocounts = cctools_pseudocounts, record = True, save = cctools_peaks)
bg_for_homer_ccftools =  cc.pp.call_peaks(bg,  method = "MACCs", reference = args.genome, pvalue_cutoffbg = cctools_pvalue_cutoffbg, window_size=cctools_window_size, step_size = cctools_step_size, pvalue_cutoffTTAA = cctools_pvalue_cutoffTTAA, lam_win_size = cctools_lam_win_size,  pseudocounts = cctools_pseudocounts, record = True, save = cctools_peak_bg_for_homer)

peak_data_cccaller = cc.pp.call_peaks(TF, bg,  method = "CCcaller", reference = args.genome, pvalue_cutoffbg = cccaller_pvalue_cutoffbg, maxbetween = cccaller_maxbetween, pvalue_cutoffTTAA = cccaller_pvalue_cutoffTTAA, lam_win_size = ccaller_lam_win_size, pseudocounts = cccaller_pseudocounts, record = True, save = cccaller_peaks)
bg_for_homer_cccaller = cc.pp.call_peaks(bg,  method = "CCcaller", reference = args.genome, pvalue_cutoffbg = cccaller_pvalue_cutoffbg, maxbetween = cccaller_maxbetween, pvalue_cutoffTTAA = cccaller_pvalue_cutoffTTAA, lam_win_size = ccaller_lam_win_size, pseudocounts = cccaller_pseudocounts, record = True, save = cccaller_peak_bg_for_homer)

cctools_graph = args.out_stem + "/" + args.tf_stem + "_" + args.unf_stem + "cctools_peaks.png"
cccaller_graph = args.out_stem + "/" + args.tf_stem + "_" + args.unf_stem + "cccaller_peaks.png"

cc.pl.draw_area("chr1",4856929,4863861,100000,peak_data_cctools, TF, args.genome, bg, figsize = (30,8),peak_line = 2,save = cctools_graph, bins = 500, example_length = 5000)
cc.pl.draw_area("chr1",4856929,4863861,100000,peak_data_cccaller, TF, args.genome, bg, figsize = (30,8),peak_line = 2,save = cccaller_graph, bins = 500, example_length = 5000)

qbed= {args.tf_stem:TF, args.unf_stem:bg}
bed = {'PEAK_cc_tools':peak_data_cctools,'PEAK_CCcaller':peak_data_cccaller}
cc.pl.WashU_browser_url(qbed,bed,genome = args.genome)

