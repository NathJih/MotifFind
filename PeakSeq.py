import sys
import GetRefGenome

if len(sys.argv != 3):
    raise Exception("Incorrect number of arguments, 3 expected.")

peaks_file = sys.argv[1]

def PeakSeq(peaks_file):
    peakseqs = []
    # open up the file and extract the lines 
    peaks = open(peaks_file)
    peakslines = peaks.readlines()
    # for each peak, extracts the specified region from the reference genome
    for pl in peakslines:
        pl.split()
        chrom = pl[0].split(":")
        chromID = chrom[1]
        ref_genome = GetRefGenome(chromID)
        start = pl[1].split(":")
        startSite = start[1]
        end = pl[2].split(":")
        endSite = end[1]
        peakseqs.append(ref_genome[startSite - 1 : endSite])
    return peakseqs
PeakSeq(peaks_file)
