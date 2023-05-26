import sys

if len(sys.argv != 3):
    raise Exception("Incorrect number of arguments, 3 expected.")

peaks_file = sys.argv[1]
ref_genome_file = sys.argv[2]

def PeakSeq(peaks_file, ref_genome_file):
    # open up the files and extract the lines 
    peaks = open(peaks_file)
    ref_genome = open(ref_genome_file)

PeakSeq(peaks_file, ref_genome_file)