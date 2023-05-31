import sys

if len(sys.argv) != 3:
    raise Exception("Incorrect number of arguments, 3 expected.")

peaks_file = sys.argv[1]

class SequenceData:
    PeakSeq = []
    BackgroundSeq = []
<<<<<<< HEAD
    RefGen = {}
=======
>>>>>>> 33a6ba012b1be0b0930d31727f288dd403c3e2e1

def GetRefGenome(chromosome):
    RGFile = open("ReferenceGenome/chromosome" + str(chromosome) + ".fna")
    RGLines = RGFile.readlines()
    RefGen = ""
    for i in range(1, len(RGLines)):
        RefGen = RefGen + RGLines[i].strip()
    return RefGen

def PeakSeq(peaks_file):
    peakseqs = []
    # open up the file and extract the lines 
    peaks = open(peaks_file)
    peakslines = peaks.readlines()
    # for each peak, extracts the specified region from the reference genome
    for pl in peakslines:
        pllist = pl.split("\t")
        chrom = pllist[0]
        chromID = chrom[3:]
        print(chromID)
        ref_genome = GetRefGenome(chromID)
        startSite = pllist[1]
        endSite = pllist[2]
        peakseqs.append(ref_genome[int(startSite) - 1 : int(endSite)])
    return peakseqs

output = open("peakouttest1.txt", "w")
for i in PeakSeq(peaks_file):
    output.write(i + "\n")