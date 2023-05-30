import GetRefGenome
import PeakSeq
import random as rand
import scipy.stats
import sys
import Jaspar

arg_len = len(sys.argv)
args = sys.argv
j_flag = False
o_flag = False

# checks for correct command 
if arg_len < 3:
    raise Exception("Incorrect number of arguments, 3 expected.")
elif arg_len > 4:
    if args[3] == '-j':
        # jasper input file 
        jasperfile = args[4]
        j_flag = True
    elif args[3] == '-o':
        if args[4].strip()[:-4] != '.txt':
            raise Exception("Missing .txt after the file name.")
        else:
            # peak txt output file 
            o_flag = True
            peakseq_output = args[4]
elif arg_len == 7:
     # jasper input file
    if args[5] == '-j':
        j_flag = True
        jasperfile = args[6]
    elif args[5] == '-o':
        if args[6].strip()[:-4] != '.txt':
            raise Exception("Missing .txt after the file name.")
        else:
            # peak txt output file 
            o_flag = True
            peakseq_output = args[6]
else:
    raise Exception("Incorrect arguments.")

peakfile = args[1]
outputfile = open(args, "w")
    

def getBackground(peaks):
    """
    Helper method to get background sequences for analysis.
    """
    bg_seqs = []

    #Method 1: Random locations in genome
    for peak in peaks:
        peak_len = peak[3] - peak[2]
        start = rand.randint(0, 57000001)
        end = start + peak_len
        chrom = peak[0]
        RG = GetRefGenome(chrom)
        bg_seq = RG[start : end]
        bg_seqs.append(bg_seq)

    return bg_seqs

# helper function to do the comparisons 
def FindMatch(matrix, sequence):
    count = 0
    return count

def ComputeEnrichment(peak_total, peak_motif, bg_total, bg_motif):
    """ Compute fisher exact test to test whether motif enriched in bound sequences
    
    Parameters
    ----------
    peak_total : int
       Number of total peaks
    peak_motif : int
       Number of peaks matching the motif
    bg_total : int
       Number of background sequences
    bg_motif : int
       Number of background sequences matching the motif
       
    Returns
    -------
    pval : float
       Fisher Exact Test p-value    
    """
    pval = -1
    table = [[peak_motif, bg_motif], [peak_total - peak_motif, bg_total - bg_motif]]
    odds, pval = scipy.stats.fisher_exact(table)
    return pval

# still need to figure out where to get the motif and stuff 
def MotifFind(motif):
    motif_list = []

    # add the motif name to the list
    motif_list.append(motif.key()) 
    motif_pwm = motif[0]

    # find the number of peak sequences that match the motif 
    peakseqs = PeakSeq(peakfile)
    peaks_count = FindMatch(motif_pwm, peakseqs)
    motif_list.append(peaks_count)

    # find the number of background sequences that match the motif
    bgseqs = getBackground(peakseqs)
    bg_count = FindMatch(motif_pwm, bgseqs)
    motif_list.append(bg_count)

    # Fisher Exact test to calculate the p-value
    pval = ComputeEnrichment(len(peakseqs), peaks_count, len(bgseqs), bg_count)
    motif_list.append(pval)

    return motif_list
