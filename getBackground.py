import random as rand
def getBackground(ref_genome, peaks):
    bg_seqs = []
    
    #Method 1: Random locations in genome
    for peak in peaks:
        peak_len = peak[3] - peak[2]
        start = rand.randint(0, 57000001)
        end = start + peak_len
        chrom = peak[0]
        #find the corresponding chromosome in the ref genome and extract the sequence using start and end
        bg_seq = ""
        bg_seqs.append(bg_seq)
