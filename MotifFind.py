import GetRefGenome
import random as rand

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

