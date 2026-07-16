#!/usr/bin/env python
"""
# Author: 
# Created Time : 
# File Name: 
# Description:
"""
import numpy as np
import pandas as pd
import scipy.stats as scistat
import tqdm
from scipy.sparse import lil_matrix
#import anndata as ad
from scipy.sparse import csr_matrix


def findinsertionslen(Chrom, start, end):
    # function to calculate the number of hops in the spcific area of chromosomes
    return len(Chrom[(Chrom >= start-3) & (Chrom <= end )])

def findinsertions(Chrom, start, end):
    # function returns of hops in the spcific area of chromosomes
    return Chrom[(Chrom >= start-3) & (Chrom <= end)]

def compute_cumulative_poisson(exp_hops_region,bg_hops_region,total_exp_hops,total_bg_hops,pseudocounts):
    
    if total_bg_hops >= total_exp_hops:
        return(1-scistat.poisson.cdf((exp_hops_region+pseudocounts),bg_hops_region * (float(total_exp_hops)/float(total_bg_hops)) + pseudocounts))
    else:
        return(1-scistat.poisson.cdf(((exp_hops_region *(float(total_bg_hops)/float(total_exp_hops)) )+pseudocounts),bg_hops_region + pseudocounts))

def testCompare(bound, curChromnp, curframe, boundnew, scaleFactor, pvalue_cutoff, chrom, record):
    
    # test whether the potiential peaks are true peaks by comparing to other data
    for i in range(len(bound)):
        TTAAnum = len(curframe[(curframe >= bound[i][0]) & (curframe <= bound[i][1])])
        boundnum = bound[i][2]
        #TTAAnumlam = len(curframe[(curframe >= bound[i][0] - lam_size/2) & (curframe <= bound[i][1] + lam_size/2)])
        #boundnumlam = len(curChromnp[(curChromnp >= bound[i][0] - lam_size/2) & (curChromnp <= bound[i][1] + lam_size/2)])
        lam = TTAAnum * scaleFactor
        pvalue = scistat.poisson.sf(boundnum - 1, lam)
        
        if pvalue < pvalue_cutoff:
            if record:
                boundnew.append([chrom, bound[i][0], bound[i][1], pvalue])
            else:
                boundnew.append([chrom, bound[i][0], bound[i][1]])
                
    return boundnew
    



def blockifyCompare(bound, curChrom, curframe, boundnew, scaleFactor, pvalue_cutoff, chrom, record):
# test whether the potiential peaks are true peaks by comparing to TTAAs
    last = -1
    Chrnumtotal = 0
    TTAAnumtotal = 0

    for i in range(len(bound)):
        TTAAnum = len(curframe[(curframe >= bound[i][0]) & (curframe <= bound[i][1])])
        boundnum = bound[i][2]

        lam = TTAAnum * scaleFactor
        pValue = scistat.poisson.sf(boundnum - 1, lam)

        if pValue <= pvalue_cutoff and last == -1 :
            
            last = i
            Chrnumtotal += boundnum
            TTAAnumtotal += TTAAnum

        elif pValue > pvalue_cutoff and last != -1 :
            
            if record:
                boundnew.append([chrom, bound[last][0], bound[i-1][1], scistat.poisson.sf(Chrnumtotal - 1, TTAAnumtotal*scaleFactor)])
            else:
                boundnew.append([chrom, bound[last][0], bound[i-1][1]])
                
            last = -1
            Chrnumtotal = 0
            TTAAnumtotal = 0


    if last != -1:
        
        if record:
            boundnew.append([chrom, bound[last][0], bound[i-1][1], scistat.poisson.sf(Chrnumtotal - 1, TTAAnumtotal*scaleFactor)])
        else:
            boundnew.append([chrom, bound[last][0], bound[i-1][1]])

    return boundnew



def callpeaksMACS2_bf(expdata, TTAAframe, 
                      window_size = 1000,lam_win_size=100000,step_size = 500,pseudocounts = 0.2,pvalue_cutoff = 0.01,
                     record = False):
    
    # function for MACS2 background free
    
    # The chromosomes we need to consider
    chrm = list(expdata[0].unique())
    
    chr_list = []
    start_list = []
    end_list = []
    
    if record:
        center_list = []
        num_exp_hops_list = []
        frac_exp_list = []
        tph_exp_list = []
        lambda_type_list =[]
        lambda_list = []
        pvalue_list = []
        l = []
        
    total_experiment_hops = len(expdata)
    list_of_l_names = lam_win_size
    
    
    for chrom in tqdm.tqdm(chrm):
    #for chrom in chrm:

        # find out the insertion data of our current chromosome
        curChrom = list(expdata[expdata[0] == chrom][1]) 
        curChromnp = np.array(curChrom)
        # sort it so that we could accelarate the searching afterwards
        curChrom.sort()
        
        max_pos = curChrom[-1]
        sig_start = 0
        sig_end = 0
        sig_flag = 0
        
        curTTAAframe = np.array(list(TTAAframe[TTAAframe[0]==chrom][1]))
        
        for window_start in range(0,int(max_pos+window_size+(lam_win_size/2+1)),step_size):
            
            num_exp_hops = findinsertionslen(curChromnp, window_start, window_start+window_size - 1)
            if num_exp_hops>1:
                num_lam_win_size_hops = findinsertionslen(curChromnp, window_start-(lam_win_size/2-1), 
                                                       window_start+window_size +(lam_win_size/2) - 1)
                num_TTAAs_window = findinsertionslen(curTTAAframe, window_start, window_start+window_size - 1)
                num_TTAAs_lam_win_size = findinsertionslen(curTTAAframe, window_start-(lam_win_size/2-1), 
                                                       window_start+window_size +(lam_win_size/2) - 1)
                lambda_lam_win_size = (num_lam_win_size_hops/(max(num_TTAAs_lam_win_size,1))) #expected number of hops per TTAA

                #is this window significant?
                pvalue = 1-scistat.poisson.cdf((num_exp_hops+pseudocounts),lambda_lam_win_size*max(num_TTAAs_window,1)+pseudocounts)
            else:
                pvalue = 1
                
            if pvalue < pvalue_cutoff :
                #was last window significant?
                if sig_flag:
                    #if so, extend end of windows
                    sig_end = window_start+window_size-1
                else:
                    #otherwise, define new start and end and set flag
                    sig_start = window_start
                    sig_end = window_start+window_size-1
                    sig_flag = 1

            else:
                #current window not significant.  Was last window significant?
                if sig_flag:
                    

                    #compute peak center and add to frame
                    overlap = findinsertions(curChromnp, sig_start, sig_end)
                    peak_center = np.median(overlap)
            
                    #add number of experiment hops in peak to frame
                    num_exp_hops = len(overlap)
                    
                    num_TTAAs_peak = findinsertionslen(curTTAAframe, sig_start, sig_end)

                    #compute lambda in lam_win_size
                    num_exp_hops_lam_win_size = findinsertionslen(curChromnp, peak_center-(lam_win_size/2-1), peak_center+(lam_win_size/2))
                    num_TTAAs_lam_win_size = findinsertionslen(curTTAAframe, peak_center-(lam_win_size/2-1), peak_center+(lam_win_size/2))
                    lambda_lam_win_size = float(num_exp_hops_lam_win_size)/(max(num_TTAAs_lam_win_size,1))


                    lambda_f = lambda_lam_win_size

                    #compute pvalue and record it

                    pvalue = 1-scistat.poisson.cdf((num_exp_hops+pseudocounts),lambda_f*max(num_TTAAs_peak,1)+pseudocounts)
                    pvalue_list.append(pvalue)              
                    #number of hops that are a user-defined distance from peak center
                    
                    if record:
                        
                        #add chr, peak start, peak end
                        chr_list.append(chrom) #add chr to frame
                        start_list.append(sig_start) #add peak start to frame
                        end_list.append(sig_end) #add peak end to frame
                        center_list.append(peak_center) #add peak center to frame
                        num_exp_hops_list.append(num_exp_hops)

                        #add fraction of experiment hops in peak to frame
                        frac_exp_list.append(float(num_exp_hops)/total_experiment_hops)
                        tph_exp_list.append(float(num_exp_hops)*100000/total_experiment_hops)
                        #record type of lambda used
                        lambda_type_list.append(list_of_l_names)
                        #record lambda
                        lambda_list.append(lambda_f)
                        
                    sig_flag = 0


                    
    if record:
        peaks_frame = pd.DataFrame(columns = ["Chr","Start","End","Center","Experiment Hops",
                "Fraction Experiment","TPH Experiment","Lambda Type",
                "Lambda","Poisson pvalue"])

        peaks_frame["Center"] = center_list
        peaks_frame["Experiment Hops"] = num_exp_hops_list 
        peaks_frame["Fraction Experiment"] = frac_exp_list 
        peaks_frame["TPH Experiment"] = tph_exp_list
        peaks_frame["Lambda Type"] = lambda_type_list
        peaks_frame["Lambda"] = lambda_list

    else:
        peaks_frame = pd.DataFrame(columns = ["Chr","Start","End","Poisson pvalue"])

    peaks_frame["Chr"] = chr_list
    peaks_frame["Start"] = start_list
    peaks_frame["End"] = end_list
    peaks_frame["Poisson pvalue"] = pvalue_list
                    
    peaks_frame = peaks_frame[peaks_frame["Poisson pvalue"] <= pvalue_cutoff]
    
    if record:
        return peaks_frame
    else:                    
        return peaks_frame[['Chr','Start','End']]
    
    
def callpeaksMACS2(expdata, background, TTAAframe, 
                      window_size = 1000,lam_win_size=100000,step_size = 500,pseudocounts = 0.2,pvalue_cutoff = 0.01,
                     record = False):
    
    # function for MACS2 with background 
    
    # The chromosomes we need to consider
    chrm = list(expdata[0].unique())
    
    chr_list = []
    start_list = []
    end_list = []
    list_of_l_names = ["bg","1k","5k","10k"]
    pvalue_list = []
    
    if record:
        chr_list = []
        start_list = []
        end_list = []
        center_list = []
        num_exp_hops_list = []
        num_bg_hops_list = []
        frac_exp_list = []
        tph_exp_list = []
        frac_bg_list = []
        tph_bg_list = []
        tph_bgs_list = []
        lambda_type_list =[]
        lambda_list = []
        
        
    total_experiment_hops = len(expdata)
    total_background_hops = len(background)
    
    
    
    for chrom in tqdm.tqdm(chrm):

        # find out the insertion data of our current chromosome
        curChrom = list(expdata[expdata[0] == chrom][1]) 
        curChromnp = np.array(curChrom)
        # sort it so that we could accelarate the searching afterwards
        curChrom.sort()
        
        max_pos = curChrom[-1] + 4
        sig_start = 0
        sig_end = 0
        sig_flag = 0
        
        curTTAAframe = np.array(list(TTAAframe[TTAAframe[0]==chrom][1]))
        
        curbackgroundframe = np.array(list(background[background[0]==chrom][1]))
        
        sig_start = 0
        sig_end = 0
        sig_flag = 0
        
        for window_start in range(0,int(max_pos+window_size),int(step_size)):

            num_exp_hops = findinsertionslen(curChromnp, window_start, window_start+window_size - 1)
            if num_exp_hops > 1:
                num_bg_hops = findinsertionslen(curbackgroundframe, window_start, window_start+window_size - 1)
                p = compute_cumulative_poisson(num_exp_hops,num_bg_hops,total_experiment_hops,total_background_hops,pseudocounts)
            else:
                p = 1
                

            #is this window significant?
            if p < pvalue_cutoff:
                #was last window significant?
                if sig_flag:
                    #if so, extend end of windows
                    sig_end = window_start+window_size-1
                else:
                    #otherwise, define new start and end and set flag
                    sig_start = window_start
                    sig_end = window_start+window_size-1
                    sig_flag = 1

            else:
                #current window not significant.  Was last window significant?
                if sig_flag:
                    
                    #add full sig window to the frame of peaks 
                    #add chr, peak start, peak end
                    chr_list.append(chrom) #add chr to frame
                    start_list.append(sig_start) #add peak start to frame
                    end_list.append(sig_end) #add peak end to frame
                    
                    #compute peak center and add to frame
                    overlap = findinsertions(curChromnp, sig_start, sig_end)
                    peak_center = np.median(overlap)

                    #add number of experiment hops in peak to frame
                    num_exp_hops = len(overlap)

                    #add number of background hops in peak to frame
                    num_bg_hops = findinsertionslen(curbackgroundframe, sig_start, sig_end)
                    
                  
                    if record:
                        center_list.append(peak_center) #add peak center to frame
                        num_exp_hops_list.append(num_exp_hops)
                        #add fraction of experiment hops in peak to frame
                        frac_exp_list.append(float(num_exp_hops)/total_experiment_hops)
                        tph_exp_list.append(float(num_exp_hops)*100000/total_experiment_hops)
                        num_bg_hops_list.append(num_bg_hops)
                        frac_bg_list.append(float(num_bg_hops)/total_background_hops)
                        tph_bg_list.append(float(num_bg_hops)*100000/total_background_hops)

                    #find lambda and compute significance of peak
                    if total_background_hops >= total_experiment_hops: #scale bg hops down
                        #compute lambda bg
                        num_TTAAs = findinsertionslen(curTTAAframe, sig_start, sig_end)
                        lambda_bg = ((num_bg_hops*(float(total_experiment_hops)/total_background_hops))/max(num_TTAAs,1)) 


                        #compute lambda 1k
                        num_bg_hops_1k = findinsertionslen(curbackgroundframe, peak_center-499, peak_center+500)
                        num_TTAAs_1k = findinsertionslen(curTTAAframe, peak_center-499, peak_center+500)
                        lambda_1k = (num_bg_hops_1k*(float(total_experiment_hops)/total_background_hops))/(max(num_TTAAs_1k,1))


                        #compute lambda 5k
                        num_bg_hops_5k = findinsertionslen(curbackgroundframe, peak_center-2499, peak_center+2500)
                        num_TTAAs_5k = findinsertionslen(curTTAAframe, peak_center-2499, peak_center+2500)
                        lambda_5k = (num_bg_hops_5k*(float(total_experiment_hops)/total_background_hops))/(max(num_TTAAs_5k,1))


                        #compute lambda 10k
                        num_bg_hops_10k = findinsertionslen(curbackgroundframe, peak_center-4999, peak_center+5000)
                        num_TTAAs_10k = findinsertionslen(curTTAAframe, peak_center-4999, peak_center+5000)
                        lambda_10k = (num_bg_hops_10k*(float(total_experiment_hops)/total_background_hops))/(max(num_TTAAs_10k,1))
                        lambda_f = max([lambda_bg,lambda_1k,lambda_5k,lambda_10k])


                        #record type of lambda used
                        index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k].index(max([lambda_bg,lambda_1k,lambda_5k,lambda_10k]))
                        lambda_type_list.append(list_of_l_names[index])
                        #record lambda
                        lambda_list.append(lambda_f)
                        #compute pvalue and record it

                        pvalue = 1-scistat.poisson.cdf((num_exp_hops+pseudocounts),lambda_f*max(num_TTAAs,1)+pseudocounts)
                        pvalue_list.append(pvalue)


                        tph_bgs = float(num_exp_hops)*100000/total_experiment_hops-float(num_bg_hops)*100000/total_background_hops
                        if record:
                            tph_bgs_list.append(tph_bgs)


                        index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k].index(max([lambda_bg,lambda_1k,lambda_5k,lambda_10k]))
                        lambdatype = list_of_l_names[index]
                        #l = [pvalue,tph_bgs,lambda_f,lambdatype]

                    else: #scale experiment hops down
                        

                        #compute lambda bg
                        num_TTAAs = findinsertionslen(curTTAAframe, sig_start, sig_end)
                        lambda_bg = (float(num_bg_hops)/max(num_TTAAs,1)) 


                        #compute lambda 1k
                        num_bg_hops_1k = findinsertionslen(curbackgroundframe, peak_center-499, peak_center+500)
                        num_TTAAs_1k = findinsertionslen(curTTAAframe, peak_center-499, peak_center+500)
                        lambda_1k = (float(num_bg_hops_1k)/(max(num_TTAAs_1k,1)))


                        #compute lambda 5k
                        num_bg_hops_5k = findinsertionslen(curbackgroundframe, peak_center-2499, peak_center+2500)
                        num_TTAAs_5k = findinsertionslen(curTTAAframe, peak_center-2499, peak_center+2500)
                        lambda_5k = (float(num_bg_hops_5k)/(max(num_TTAAs_5k,1)))


                        #compute lambda 10k
                        num_bg_hops_10k = findinsertionslen(curbackgroundframe, peak_center-4999, peak_center+5000)
                        num_TTAAs_10k = findinsertionslen(curTTAAframe, peak_center-4999, peak_center+5000)
                        lambda_10k = (float(num_bg_hops_10k)/(max(num_TTAAs_10k,1)))
                        lambda_f = max([lambda_bg,lambda_1k,lambda_5k,lambda_10k])


                        #record type of lambda used
                        index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k].index(max([lambda_bg,
                                                                                      lambda_1k,lambda_5k,lambda_10k]))

                        #compute pvalue and record it
                        pvalue = 1-scistat.poisson.cdf(((float(total_background_hops)/total_experiment_hops)*num_exp_hops+ pseudocounts),lambda_f*max(num_TTAAs,1)+pseudocounts)
                        pvalue_list.append(pvalue)

                        tph_bgs = float(num_exp_hops)*100000/total_experiment_hops -float(num_bg_hops)*100000/total_background_hops

                        if record:
                            lambda_type_list.append(list_of_l_names[index])
                            lambda_list.append(lambda_f)
                            tph_bgs_list.append(tph_bgs)


                        index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k].index(max([lambda_bg,
                                                                                      lambda_1k,lambda_5k,lambda_10k]))
                        lambdatype = list_of_l_names[index]


                    #number of hops that are a user-defined distance from peak center
                    sig_flag = 0
                        
     


    if record:
        peaks_frame = pd.DataFrame(columns = ["Chr","Start","End","Center","Experiment Hops",
                "Fraction Experiment","TPH Experiment","Lambda Type",
                "Lambda","Poisson pvalue"])


        peaks_frame["Center"] = center_list
        peaks_frame["Experiment Hops"] = num_exp_hops_list 
        peaks_frame["Fraction Experiment"] = frac_exp_list 
        peaks_frame["TPH Experiment"] = tph_exp_list
        peaks_frame["Background Hops"] = num_bg_hops_list 
        peaks_frame["Fraction Background"] = frac_bg_list
        peaks_frame["TPH Background"] = tph_bg_list
        peaks_frame["TPH Background subtracted"] = tph_bgs_list
        peaks_frame["Lambda Type"] = lambda_type_list
        peaks_frame["Lambda"] = lambda_list
        
    else:
        peaks_frame = pd.DataFrame(columns = ["Chr","Start","End"])

    peaks_frame["Chr"] = chr_list
    peaks_frame["Start"] = start_list
    peaks_frame["End"] = end_list
    peaks_frame["Poisson pvalue"] = pvalue_list
        

    peaks_frame = peaks_frame[peaks_frame["Poisson pvalue"] <= pvalue_cutoff]
    
    
    if record:
        return peaks_frame
    else:                    
        return peaks_frame[['Chr','Start','End']]

def checkint(number,name):
    try:
        number = int(number)
    except ValueError:
        print('Please enter a valid positive number or 0 for' + name)
    if number <0 :
        raise ValueError('Please enter a valid positive number or 0 for' + name)
       
    return number
    
    
def callpeaks(expdata, background = None, method = "test", 
              reference = "hg38", pvalue_cutoff = 0.01, 
              mininser = 5, minlen = 0, extend = 150, maxbetween = 2800, lam_size = 10000, 
              window_size = 1000, lam_win_size=100000, step_size = 500, pseudocounts = 0.2, record = False
):
    
    # function for calling peaks
    
    # first make sure everything is correct:
    # check expdata
    if type(expdata) != pd.DataFrame :
        raise ValueError("Please input a pandas dataframe as the expression data.")
        
    # check pvalue_cutoff
    try:
        pvalue_cutoff = float(pvalue_cutoff)
    except ValueError:
        print('Please enter a valid number (0,1) for pvalue_cutoff')

    if pvalue_cutoff<= 0 or pvalue_cutoff >=1:
        raise ValueError('Please enter a valid number (0,1) for pvalue_cutoff')
        
    
    if method == "test":
        #print("For the test method, [expdata, background/reference, pvalue_cutoff, mininser, minlen, extend, maxbetween] would be utilized.")
        checkint(mininser,"mininser")
        checkint(minlen,"minlen")
        extend = checkint(extend,"extend")
        checkint(maxbetween,"maxbetween")
        lam_size = checkint(lam_size,"lam_size")
        
    elif  method == "Blockify":
        print("For the Blockify method, [expdata, background/reference, pvalue_cutoff] would be utilized.")
    
    elif  method == "MACS2":
        #print("For the MACS2 method, [expdata, background (could be ignbored), reference, pvalue_cutoff, window_size, lam_win_size, step_size, pseudocounts, record] would be utilized.")
        window_size = checkint(window_size,"window_size")
        lam_win_size = checkint(lam_win_size,"lam_win_size")
        step_size = checkint(step_size,"step_size")
        checkint(pseudocounts,"pseudocounts")
        
        if type(record) != bool:
            raise ValueError('Please enter a True/ False for record')
            
            
    else:
        raise ValueError("Not valid method.")
         
            
        
    if method == "test" or method == "Blockify":
        
        if type(background) != pd.DataFrame :
            
            if reference == "hg38":
                TTAAframe = pd.read_csv("hg38_TTAA.bed",delimiter="\t",header=None)
            elif reference == "mm10":
                TTAAframe = pd.read_csv("mm10_TTAA.bed",delimiter="\t",header=None)
            else:
                raise ValueError("Not valid reference.")
           
            scaleFactor = len(expdata) / len(TTAAframe)
            
        else:
            scaleFactor = len(expdata) / len(background)
        
    
        # The chromosomes we need to consider
        chrm = list(expdata[0].unique())

        # create a list to record every peaks
        boundnew = []

        # going through one chromosome from another
        for chrom in tqdm.tqdm(chrm):

            # find out the insertion data of our current chromosome
            curChrom = list(expdata[expdata[0] == chrom][1]) 
            curChromnp = np.array(curChrom)
            # sort it so that we could accelarate the searching afterwards
            curChrom.sort()

            if type(background) != pd.DataFrame:
                # find out the insertion data of our current chromosome
                curTTAAframe = np.array(list(TTAAframe[TTAAframe[0]==chrom][1]))

            else:
                # find out the insertion data of the backgrround chromosome
                curbackgroundframe = np.array(list(background[background[0]==chrom][1]))
                curbackgroundframe.sort()


            # make a summary of our current insertion start points
            unique, counts = np.unique(np.array(curChrom), return_counts=True)


            # create a list to find out the protintial peak region of their bounds
            bound = []

            if method == "test":

                # initial the start point, end point and the totol number of insertions
                startbound = 0
                endbound = 0
                insertionbound = 0

                # calculate the distance between each points
                dif1 = np.diff(unique, axis=0)
                # add a zero to help end the following loop at the end
                dif1 = np.concatenate((dif1,np.array([maxbetween+1])))


                # look for the uique insertion points
                for i in range(len(unique)):
                    if startbound == 0:
                        startbound = unique[i]
                        insertionbound += counts[i]
                        if dif1[i] > maxbetween :
                            endbound = unique[i] 
                            if (insertionbound >= mininser) and ((endbound- startbound) >= minlen):
                                bound.append([max(startbound-extend,0),endbound+4+extend,insertionbound, endbound-startbound])
                            startbound = 0
                            endbound = 0
                            insertionbound = 0 
                    else:
                        insertionbound += counts[i]
                        if dif1[i] > maxbetween :
                            endbound = unique[i]
                            if (insertionbound >= mininser) and ((endbound- startbound) >= minlen):
                                bound.append([max(startbound-extend,0),endbound+4+extend,insertionbound, endbound-startbound])
                            startbound = 0
                            endbound = 0
                            insertionbound = 0


                
                if type(background) != pd.DataFrame:
                    boundnew = testCompare(bound, curChromnp, curTTAAframe, boundnew, scaleFactor, pvalue_cutoff, chrom, record)
                else:
                    boundnew = testCompare(bound, curChromnp, curbackgroundframe, boundnew, scaleFactor, pvalue_cutoff, chrom, record)



            elif method == "Blockify":

                import astropy.stats as astrostats
                hist, bin_edges = astrostats.histogram(expdata[expdata[0] == chrom][1], bins="blocks")

                hist = list(hist)
                bin_edges = list(bin_edges.astype(int))
                for bins in range(len(bin_edges)-1):
                    bound.append([bin_edges[bins],bin_edges[bins+1], hist[bins], bin_edges[bins+1]-bin_edges[bins]])


                if type(background) != pd.DataFrame:
                    boundnew = blockifyCompare(bound, curChrom, curTTAAframe, boundnew, scaleFactor, pvalue_cutoff, chrom, record)
                else:
                    boundnew = blockifyCompare(bound, curChrom, curbackgroundframe, boundnew, scaleFactor, pvalue_cutoff, chrom, record)

        if record:
            return pd.DataFrame(boundnew, columns=["Chr","Start", "End", "pvalue"])
        else:
            return pd.DataFrame(boundnew, columns=["Chr","Start", "End"])
    
        
    elif method == "MACS2":
        
        if reference == "hg38":
            TTAAframe = pd.read_csv("hg38_TTAA.bed",delimiter="\t",header=None)
        elif reference == "mm10":
            TTAAframe = pd.read_csv("mm10_TTAA.bed",delimiter="\t",header=None)
        else:
            raise ValueError("Not valid reference.")

        if type(background) != pd.DataFrame :
            return callpeaksMACS2_bf(expdata, TTAAframe, window_size = window_size, 
                                     lam_win_size= lam_win_size, step_size = step_size,
                                     pseudocounts = pseudocounts, pvalue_cutoff = pvalue_cutoff, record = record)
        else:
            return callpeaksMACS2(expdata, background, TTAAframe, window_size = window_size, 
                                     lam_win_size= lam_win_size, step_size = step_size,
                                     pseudocounts = pseudocounts, pvalue_cutoff = pvalue_cutoff, record = record)
            
            
# a sorting function to sort chromosome
def myFuncsorting(e):
    try:
        return int(e.split('_')[0][3:])
    except ValueError:
        return int(ord(e.split('_')[0][3:]))
'''
# pair function 
def peakToCell_function(insertionschr, peakschr, chromo, cellpeaks, peak_name_dict, barcodes_dict):
    peak_id = 0
    
    for insertion in range(len(insertionschr)):
        
        if (insertionschr.iat[insertion,1] >= peakschr.iat[peak_id,1]-3) and insertionschr.iat[insertion,2] <= (peakschr.iat[peak_id,2]+3):
            cellpeaks[peak_name_dict[chromo+"_"+str(peakschr.iat[peak_id,1])+"_"+str(peakschr.iat[peak_id,2])],barcodes_dict[insertionschr.iat[insertion,5]]]+= 1
        elif peak_id +1 == len(peakschr):
            break
        elif insertionschr.iat[insertion,1] >= peakschr.iat[peak_id+1,1]-3 and insertionschr.iat[insertion,2] <= peakschr.iat[peak_id+1,2]+3:
            cellpeaks[peak_name_dict[chromo+"_"+str(peakschr.iat[peak_id+1,1])+"_"+str(peakschr.iat[peak_id+1,2])],barcodes_dict[insertionschr.iat[insertion,5]]]+= 1
            peak_id = peak_id +1

    return cellpeaks

def peakToCell(insertions, peaks, barcodes):
    
    
    barcodes = list(barcodes[0])
    barcodes_dict = {}
    for i in range(len(barcodes)):
        barcodes_dict[barcodes[i]] = i

    # peaks to order
    peak_name = peaks[[0,1,2]].T.astype(str).agg('_'.join)
    peak_name_unique = list()
    for number in peak_name:
        if number not in peak_name_unique:
            peak_name_unique.append(number)
    peak_name_unique.sort(key=myFuncsorting)

 
    peak_name_dict = {}
    for i in range(len(peak_name_unique)):
        peak_name_dict[peak_name_unique[i]] = i
        
    # create an empty matrix to store peaks * cell
    cellpeaks = lil_matrix((len(peak_name_unique),len(barcodes)), dtype=np.int32)
    
    #pairing
    for chromosome in list(peaks[0].unique()):

        insertionschr = insertions[insertions[0] == chromosome]
        peakschr = peaks[peaks[0] == chromosome]

        cellpeaks = peakToCell_function(insertionschr, peakschr, chromosome, cellpeaks, peak_name_dict, barcodes_dict)
        
    adata = ad.AnnData(csr_matrix(cellpeaks).T,obs= pd.DataFrame(index=barcodes),var = pd.DataFrame(index=peak_name_unique))
    
    return adata
        
'''     
