import GetRefGenome
import PeakSeq
import random as rand
import scipy.stats

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

def MotifFind(motif, peaks_file):
    motif_list = []

    # add the motif name to the list
    motif_list.append(motif.key()) 
    motif_pwm = motif[0]

    # find the number of peak sequences that match the motif 
    peakseqs = PeakSeq(peaks_file)
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
