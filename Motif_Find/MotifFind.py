import random as rand
from matplotlib import pyplot as plt
import numpy as np
import scipy.stats
import sys
import pickle
import pandas as pd
import logomaker as lm
import seqlogo
import os

arg_len = len(sys.argv)
args = sys.argv
j_flag = False
o_flag = False

jasparfile = "jaspar.p"
# checks for correct command 
if arg_len < 3:
    raise Exception("Incorrect number of arguments, 3 expected.")
elif arg_len > 4:
    if args[3] == '-j':
        # jasper input file 
        jasperfile = args[4]
        j_flag = True
    elif args[3] == '-o':
        # if args[4].strip()[:-4] != '.txt':
        #     raise Exception("Missing .txt after the file name.")
        # else:
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
# else:
#     raise Exception("Incorrect arguments.")

peakfile = args[1]
# outputfile = open(args, "w")
    
class SequenceData:

    def __init__(self):
        self.PeakSeq = []
        self.BackgroundSeq = []
        self.RefG = {}

    def GetRefGenome(self):
        # self.RefG = pickle.load(open("/Users/lxppc/desktop/Spring 2023/CSE 185/GRCh38p14.p", "rb"))
        self.RefG = pickle.load(open("GRCh38p14test.p", "rb"))
        # print(self.RefG)

    def FindSeq(self, peaks_file):

        # open up the file and extract the lines 
        peaks = open(peaks_file)
        peakslines = peaks.readlines()

        # for each peak, extracts the specified region from the reference genome
        for pl in peakslines:
            # get basic peak sequence information
            pllist = pl.split("\t")
            chrom = pllist[0]
            chromID = chrom[3:]
            ref_genome = self.RefG["chromosome" + chromID]

            # get peak sequence
            startSite = int(pllist[1])
            endSite = int(pllist[2])
            self.PeakSeq.append(ref_genome[startSite - 1 : endSite])

            # get a random background sequence
            seqlen = endSite - startSite
            chrlen = len(ref_genome)
            bgstart = rand.randint(0, chrlen-seqlen)
            self.BackgroundSeq.append(ref_genome[bgstart: bgstart+seqlen])
        return


def ScoreSeq(matrix,seq):
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t":3}
    score = 0
    for i in range(len(seq)):
        if (seq[i]=="N"):
            return -1
        score += matrix[nucs[seq[i]]][i]
    return score

# takes params
#   motifinfo: information of a specific motif
#   peaklist: list of all peak sequence from a peak file
# output number of peaks that hits a score higher than
# the threshold of the target motif
def FindMatch(motifinfo, peaklist):
    count = 0
    matrix = motifinfo[0]
    thres = motifinfo[1]
    length = motifinfo[2]
    for peak in peaklist:
        for i in range(len(peak)-length+1):
            seq = peak[i:i+length]
            # Larger than or larger or equal than
            # confirm if each peak only count once
            if ScoreSeq(matrix,seq) > thres:
                count += 1
                break
    return count

def ComputeEnrichment(peak_total, peak_motif, bg_total, bg_motif):
    table = [[peak_motif, bg_motif], [peak_total - peak_motif, bg_total - bg_motif]]
    odds, pval = scipy.stats.fisher_exact(table)
    return pval

def sortFn(list):
    return list[3]

def Output(top5motif, motifs):
    path = './Graphs'
    if not os.path.exists(path):
        os.mkdir(path)
    for motif in top5motif:
        pwm = np.array(motifs[motif[0]][0]).transpose()
        rows = []
        for i in range(0, len(pwm)):
            rows.append(str(i + 1))
        pwm = pd.DataFrame(data=pwm, columns=["A","C","G","T"])
        lm.Logo(pwm, font_name='Arial')
        plt.savefig('Graphs/' + motif[0] + '.jpg', format='jpg')

    
    

# still need to figure out where to get the motif and stuff 
def MotifFind():
    SeqData = SequenceData()
    SeqData.GetRefGenome()
    SeqData.FindSeq(peakfile)
    motifs = pickle.load(open(jasparfile, "rb"))

    motif_list = []
    nummotif = len(motifs)
    i = 0
    # loop through all the motifs
    for motif in motifs:
        i += 1
        list = []

        # add the motif name to the list
        list.append(motif) 
        motifinfo = motifs[motif]

        # find the number of peak sequences that match the motif 
        peaks_count = FindMatch(motifinfo, SeqData.PeakSeq)
        list.append(peaks_count)

        # find the number of background sequences that match the motif
        bg_count = FindMatch(motifinfo, SeqData.BackgroundSeq)
        list.append(bg_count)

        # Fisher Exact test to calculate the p-value
        pval = ComputeEnrichment(len(SeqData.PeakSeq), peaks_count, len(SeqData.BackgroundSeq), bg_count)
        list.append(pval)

        motif_list.append(list)

        print("Processing......  " + str(i) + "/" + str(nummotif) + "\n")

    motif_list.sort(reverse=False, key=sortFn)
    
    result = motif_list[0:20]

    # Output the peak sequences in a file specified by user
    if (o_flag == True):
        output = open(peakseq_output + ".txt", "w")
        for i in SeqData.PeakSeq:
            output.write(i+"\n")
    
    print(SeqData.PeakSeq)

    # Output(motif_list[:5], motifs)

    return result

result = MotifFind()
# df = pd.DataFrame(result, columns = ['Motif Name', '# of peak sequences', '# of background sequences', 'p-value'])
# print(df)
