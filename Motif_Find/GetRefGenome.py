import pickle
import sys

args = sys.argv 
instruction = "expected: python3 GetRefGenome.py <Reference Genome Fasta File> <Output File Name>"

if len(args) != 3:
    raise Exception("Incorrect number of arguments, 3 expected.\n" + instruction)

refFasta = args[1]
print(refFasta)
output = args[2]

def GetRefGenome(refFasta):
    '''
    RefG = {}
    chromosome = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"]
    for chr in chromosome:
        RGFile = open("ReferenceGenome/chromosome" + str(chr) + ".fna")
        RGLines = RGFile.readlines()
        RefGen = ""
        for i in range(1, len(RGLines)):
            RefGen = RefGen + RGLines[i].strip()
        RefG["chromosome" + str(chr)] = RefGen
    '''
     # updated GetRefGenome to combine preprocessing steps for reference genome
    RefG = {}
    chromosome = ["1"]
    # "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
    RGFile = open(refFasta)
    RGLines = RGFile.readlines()
    for chr in chromosome:
        key = "chromosome" + chr
        for i in range(0, len(RGLines)):
            if chr in RGLines[i]:
                header_lines = i
        start = header_lines + 1
        chrom = ""
        for i in range(start, len(RGLines)):
            if "chromosome" not in RGLines[i]:
                chrom += str(RGLines[i])
            else:
                break
        RefG[key] = chrom

    file = open(output + ".p", 'wb')
    pickle.dump(RefG, file)
    file.close()

    return

GetRefGenome(refFasta)