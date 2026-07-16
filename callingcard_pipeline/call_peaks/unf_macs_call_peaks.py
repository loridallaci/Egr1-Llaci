#!/usr/bin/env python
import pycallingcards as cc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--unf_bound_input', type = str)
parser.add_argument('--unf_FT_input', type = str)
parser.add_argument('--unf_bound_stem', type = str)
parser.add_argument('--unf_FT_stem', type = str)
parser.add_argument('--out_stem', type = str)
parser.add_argument('--genome', type = str)
parser.add_argument('--homer_path', type = str)
args = parser.parse_args()

##SET PEAK CALLING PARAMS##
cctools_pvalue_cutoff = 0.01 #lower means needs more significance to pass
cctools_window_size = 5000
cctools_step_size = 2500
cctools_pvalue_cutoffTTAA = 0.001
cctools_lam_win_size = 1000000
cctools_pseudocounts = 0.1

cccaller_pvalue_cutoff = 0.01
cccaller_maxbetween = 1100
cccaller_pvalue_cutoffTTAA = 0.001
ccaller_lam_win_size = 1000000
cccaller_pseudocounts = 0.1

###READ QBED FILES##
unf_bound = cc.rd.read_qbed(args.unf_bound_input)
unf_bound['group'] = "unf_acsa2_bound"

unf_FT = cc.rd.read_qbed(args.unf_FT_input)
unf_FT['group'] = "unf_acsa2_FT"

##COMBINE BOUND AND FT AND CALL PEAKS##
unf = cc.rd.combine_qbed([unf_bound, unf_FT])
unf = unf[unf['Reads'] > 2]

cctools_peaks = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "_combined_cctools_peaks.bed"
cccaller_peaks = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "_combined_cccaller_peaks.bed"

peak_data_cctools = cc.pp.call_peaks(unf,  method = "MACCs", reference = args.genome, pvalue_cutoff = cctools_pvalue_cutoff, window_size=cctools_window_size, step_size = cctools_step_size, pvalue_cutoffTTAA = cctools_pvalue_cutoffTTAA, lam_win_size = cctools_lam_win_size,  pseudocounts = cctools_pseudocounts, record = True, save = cctools_peaks)

peak_data_cccaller = cc.pp.call_peaks(unf,  method = "CCcaller", reference = args.genome, pvalue_cutoff = cccaller_pvalue_cutoff, maxbetween = cccaller_maxbetween, pvalue_cutoffTTAA = cccaller_pvalue_cutoffTTAA, lam_win_size = ccaller_lam_win_size, pseudocounts = cccaller_pseudocounts, record = True, save = cccaller_peaks)

#cctools_graph = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "_cctools_peaks.png"
#cccaller_graph = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "_cccaller_peaks.png"

#cc.pl.draw_area("chr1",4856929,4863861,100000,peak_data_cctools, unf, args.genome, figsize = (30,8),peak_line = 2,save = cctools_graph, bins = 500, example_length = 5000)
#cc.pl.draw_area("chr1",4856929,4863861,100000,peak_data_cccaller, unf, args.genome, figsize = (30,8),peak_line = 2,save = cccaller_graph, bins = 500, example_length = 5000)
'''
##MAKE BROWSER VIEW OF COMBINED PEAKS##
stem = args.unf_bound_stem + "_" + args.unf_FT_stem
qbed= {stem:unf}
bed = {'PEAK_cc_tools':peak_data_cctools,'PEAK_CCcaller':peak_data_cccaller}
cc.pl.WashU_browser_url(qbed,bed,genome = args.genome)
'''
##CALL PEAKS SEPARATELY##
#import sys
# the mock-0.3.1 dir contains testcase.py, testutils.py & mock.py
#sys.path.append('/ref/rmlab/software/spack/opt/spack/linux-rocky8-x86_64/gcc-8.5.0')
import pybedtools

cctools_bound_peaks = args.out_stem + "/" + args.unf_bound_stem + "_bound_cctools_peaks.bed"
cccaller_bound_peaks = args.out_stem + "/" + args.unf_bound_stem + "_bound_cccaller_peaks.bed"

peak_data_bound_cctools = cc.pp.call_peaks(unf_bound,  method = "MACCs", reference = args.genome, pvalue_cutoff = cctools_pvalue_cutoff, window_size=cctools_window_size, step_size = cctools_step_size, pvalue_cutoffTTAA = cctools_pvalue_cutoffTTAA, lam_win_size = cctools_lam_win_size,  pseudocounts = cctools_pseudocounts, record = True, save = cctools_bound_peaks)
peak_data_bound_cccaller = cc.pp.call_peaks(unf_bound,  method = "CCcaller", reference = args.genome, pvalue_cutoff = cccaller_pvalue_cutoff, maxbetween = cccaller_maxbetween, pvalue_cutoffTTAA = cccaller_pvalue_cutoffTTAA, lam_win_size = ccaller_lam_win_size, pseudocounts = cccaller_pseudocounts, record = True, save = cccaller_bound_peaks)

cctools_FT_peaks = args.out_stem + "/" + args.unf_FT_stem + "_FT_cctools_peaks.bed"
cccaller_FT_peaks = args.out_stem + "/" + args.unf_FT_stem + "_FT_cccaller_peaks.bed"

peak_data_FT_cctools = cc.pp.call_peaks(unf_FT,  method = "MACCs", reference = args.genome, pvalue_cutoff = cctools_pvalue_cutoff, window_size=cctools_window_size, step_size = cctools_step_size, pvalue_cutoffTTAA = cctools_pvalue_cutoffTTAA, lam_win_size = cctools_lam_win_size,  pseudocounts = cctools_pseudocounts, record = True, save = cctools_FT_peaks)
peak_data_FT_cccaller = cc.pp.call_peaks(unf_FT,  method = "CCcaller", reference = args.genome, pvalue_cutoff = cccaller_pvalue_cutoff, maxbetween = cccaller_maxbetween, pvalue_cutoffTTAA = cccaller_pvalue_cutoffTTAA, lam_win_size = ccaller_lam_win_size, pseudocounts = cccaller_pseudocounts, record = True, save = cccaller_FT_peaks)

##THIS IS FROM JUANRU'S TUTORIAL##
peak_cctools = cc.rd.combine_qbed([peak_data_bound_cctools, peak_data_FT_cctools])
peak_cccaller = cc.rd.combine_qbed([peak_data_bound_cccaller, peak_data_FT_cccaller])

peak_cctools = pybedtools.BedTool.from_dataframe(peak_cctools).merge().to_dataframe()
peak_cccaller = pybedtools.BedTool.from_dataframe(peak_cccaller).merge().to_dataframe()

peak_data_cctools = peak_cctools.rename(columns={"chrom":"Chr", "start":"Start", "end":"End"})
peak_data_cccaller = peak_cccaller.rename(columns={"chrom":"Chr", "start":"Start", "end":"End"})

peak_annotation_cctools = cc.pp.annotation(peak_data_cctools, reference = "mm10")
peak_annotation_cccaller = cc.pp.annotation(peak_data_cccaller, reference = "mm10")

peak_annotation_cctools = cc.pp.combine_annotation(peak_data_cctools, peak_annotation_cctools)
peak_annotation_cccaller = cc.pp.combine_annotation(peak_data_cccaller, peak_annotation_cccaller)

adata_cctools = cc.pp.make_Anndata(unf, peak_annotation_cctools, ["unf_acsa2_bound", "unf_acsa2_FT"], key = 'group')
adata_cccaller = cc.pp.make_Anndata(unf, peak_annotation_cccaller, ["unf_acsa2_bound", "unf_acsa2_FT"], key = 'group')

adata_cctools = cc.tl.liftover(adata_cctools)
adata_cccaller = cc.tl.liftover(adata_cccaller)

print(adata_cctools.obs)
print(adata_cccaller.obs)

##RANK PEAK GROUPS##
cctools_rank = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "_cctools_rank.png"
cccaller_rank = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "_cccaller_rank.png"

cc.tl.rank_peak_groups(adata_cctools, 'Index', method = 'fisher_exact', key_added = 'fisher_exact', rankby = 'logfoldchanges')
cc.tl.rank_peak_groups(adata_cccaller, 'Index', method = 'fisher_exact', key_added = 'fisher_exact', rankby = 'logfoldchanges')

cc.pl.rank_peak_groups(adata_cctools, key = 'fisher_exact',rankby = 'logfoldchanges', save = cctools_rank)
cc.pl.rank_peak_groups(adata_cccaller, key = 'fisher_exact',rankby = 'logfoldchanges', save = cccaller_rank)
##DRAW AREA FOR BOUND AND FT##
unf_acsa2_bound_cctools_graph = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "acsa2_bound_cctools_peaks.png"
unf_acsa2_bound_cccaller_graph = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "acsa2_bound_cccaller_peaks.png"
unf_acsa2_FT_cctools_graph = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "acsa2_FT_cctools_peaks.png"
unf_acsa2_FT_cccaller_graph = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "acsa2_FT_cccaller_peaks.png"


cc.pl.draw_area("chr12", 56516453, 56538107, 50000, peak_data_cctools, unf, "mm10", adata = adata_cctools, font_size=2, plotsize = [1,1,3], name = "unf_acsa2_bound", key = "Index",insertionkey = "group", bins = 200, name_insertion1 = 'ACSA2 Bound Insertions', name_density1 = 'ACSA2 Bound Insertion Density', name_insertion2 = 'Total Insertions', name_density2 = 'Total Insertion Density', figsize = (30,8), peak_line = 2, title = "ACSA2 Bound CCtools", save = unf_acsa2_bound_cctools_graph)

cc.pl.draw_area("chr12", 56516453, 56538107, 50000, peak_data_cccaller, unf, "mm10", adata = adata_cccaller, font_size=2, plotsize = [1,1,3], name = "unf_acsa2_bound", key = "Index",insertionkey = "group", bins = 200, name_insertion1 = 'ACSA2 Bound Insertions', name_density1 = 'ACSA2 Bound Insertion Density', name_insertion2 = 'Total Insertions', name_density2 = 'Total Insertion Density', figsize = (30,8), peak_line = 2, title = "ACSA2 Bound CCCaller", save = unf_acsa2_bound_cccaller_graph)

cc.pl.draw_area("chr12", 56516453, 56538107, 50000, peak_data_cctools, unf, "mm10", adata = adata_cctools, font_size=2, plotsize = [1,1,3], name = "unf_acsa2_FT", key = "Index",insertionkey = "group", bins = 200, name_insertion1 = 'ACSA2 FT Insertions', name_density1 = 'ACSA2 FT Insertion Density', name_insertion2 = 'Total Insertions', name_density2 = 'Total Insertion Density', figsize = (30,8), peak_line = 2, title = "ACSA2 FT CCtools", save = unf_acsa2_FT_cctools_graph)

cc.pl.draw_area("chr12", 56516453, 56538107, 50000, peak_data_cccaller, unf, "mm10", adata = adata_cccaller, font_size=2, plotsize = [1,1,3], name = "unf_acsa2_FT", key = "Index",insertionkey = "group", bins = 200, name_insertion1 = 'ACSA2 FT Insertions', name_density1 = 'ACSA2 FT Insertion Density', name_insertion2 = 'Total Insertions', name_density2 = 'Total Insertion Density', figsize = (30,8), peak_line = 2, title = "ACSA2 FT CCCaller", save = unf_acsa2_FT_cccaller_graph)


##VOLCANO PLOT BOUND VS FT##
cctools_volcano = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "_cctools_volcano.png"
cccaller_volcano = args.out_stem + "/" + args.unf_bound_stem + "_" + args.unf_FT_stem + "_cccaller_volcano.png"

cc.pl.volcano_plot(adata_cctools, pvalue_name = 'pvalues_adj', pvalue_cutoff = 0.01, lfc_cutoff = 4, figsize = (5,4), labelright = (3,250), labelleft = (-9,250), save = cctools_volcano)
cc.pl.volcano_plot(adata_cccaller, pvalue_name = 'pvalues_adj', pvalue_cutoff = 0.01, lfc_cutoff = 4, figsize = (5,4), labelright = (3,250), labelleft = (-9,250), save = cccaller_volcano)

##motif calling##
homer_bound_cctools = args.out_stem + "/" + args.unf_bound_stem + "_cctools_bound_homer.bed"
homer_bound_cccaller = args.out_stem + "/" + args.unf_bound_stem + "_cccaller_bound_homer.bed"
homer_FT_cctools = args.out_stem + "/" + args.unf_FT_stem + "_cctools_FT_homer.bed"
homer_FT_cccaller = args.out_stem + "/" + args.unf_FT_stem + "_cccaller_FT_homer.bed"
cc.tl.call_motif(peaks_frame=cc.tl.rank_peak_groups_df(adata_cctools, key = 'fisher_exact', pval_cutoff = 0.01, logfc_min = 3, group = ['unf_acsa2_bound'])["names"].str.split("_", expand=True), reference ="mm10",save_homer = "Homer/unf_acsa2_bound_cctools", homer_path = args.homer_path, num_cores=12, save_name = homer_bound_cctools)
cc.tl.call_motif(peaks_frame=cc.tl.rank_peak_groups_df(adata_cccaller, key = 'fisher_exact', pval_cutoff = 0.01, logfc_min = 3, group = ['unf_acsa2_bound'])["names"].str.split("_", expand=True), reference ="mm10",save_homer = "Homer/unf_acsa2_bound_cccaller", homer_path = args.homer_path, num_cores=12, save_name = homer_bound_cccaller)
cc.tl.call_motif(peaks_frame=cc.tl.rank_peak_groups_df(adata_cctools, key = 'fisher_exact', pval_cutoff = 0.01, logfc_min = 3, group = ['unf_acsa2_FT'])["names"].str.split("_", expand=True), reference ="mm10",save_homer = "Homer/unf_acsa2_FT_cctools", homer_path = args.homer_path, num_cores=12, save_name = homer_FT_cctools)
cc.tl.call_motif(peaks_frame=cc.tl.rank_peak_groups_df(adata_cccaller, key = 'fisher_exact', pval_cutoff = 0.01, logfc_min = 3, group = ['unf_acsa2_FT'])["names"].str.split("_", expand=True), reference ="mm10",save_homer = "Homer/unf_acsa2_FT_cccaller", homer_path = args.homer_path, num_cores=12, save_name = homer_FT_cccaller)

bound_cctools_compare_motif = cc.tl.compare_motif("Homer/unf_acsa2_bound_cctools", "Homer/unf_acsa2_FT_cctools")
bound_cccaller_compare_motif = cc.tl.compare_motif("Homer/unf_acsa2_bound_cccaller", "Homer/unf_acsa2_FT_cccaller")
FT_cctools_compare_motif = cc.tl.compare_motif("Homer/unf_acsa2_FT_cctools", "Homer/unf_acsa2_bound_cctools")
FT_cccaller_compare_motif = cc.tl.compare_motif("Homer/unf_acsa2_FT_cccaller", "Homer/unf_acsa2_bound_cccaller")

bound_cctools_compare_motif_name = args.out_stem + "/in_" + args.unf_bound_stem + "_and_not_in_" + args.unf_FT_stem + "_cctools_compare_motif.txt"
bound_cccaller_compare_motif_name = args.out_stem + "/in_" + args.unf_bound_stem + "_and_not_in_" + args.unf_FT_stem + "_cccaller_compare_motif.txt"
FT_cctools_compare_motif_name = args.out_stem + "/in_" + args.unf_FT_stem + "_vs_" + args.unf_bound_stem + "_cctools_compare_motif.txt"
FT_cccaller_compare_motif_name = args.out_stem + "/in_" + args.unf_FT_stem + "_vs_" + args.unf_bound_stem + "_cccaller_compare_motif.txt"

bound_cctools_compare_motif.to_csv(bound_cctools_compare_motif_name, sep='\t')
bound_cccaller_compare_motif.to_csv(bound_cccaller_compare_motif_name, sep='\t')
FT_cctools_compare_motif.to_csv(FT_cctools_compare_motif_name, sep='\t')
FT_cccaller_compare_motif.to_csv(FT_cccaller_compare_motif_name, sep='\t')
print("Done!")

