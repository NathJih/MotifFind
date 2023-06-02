import pickle
import sys

args = sys.argv 
instruction = "expected: python3 GetRefGenome.py <Reference_Genome_Fasta_File> <Output_File_Name>"

if len(args) != 3:
    raise Exception("Incorrect number of arguments, 3 expected.\n" + instruction)

refFasta = args[1]
output = args[2]

def GetRefGenome(refFasta):
     # updated GetRefGenome to combine preprocessing steps for reference genome
    RefG = {}
    chromosome = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
    RGFile = open(refFasta)
    RGLines = RGFile.readlines()
    for chr in chromosome:
        key = "chromosome" + chr
        for i in range(0, len(RGLines)):
            if "chromosome "+chr in RGLines[i]:
                header_lines = i
                break
        start = header_lines + 1
        chrom = ""
        for i in range(start, len(RGLines)):
            if "chr" not in RGLines[i]:
                chrom += str(RGLines[i].strip())
            else:
                break
        RefG[key] = chrom

    file = open(output + ".p", 'wb')
    pickle.dump(RefG, file)
    file.close()

    return

GetRefGenome(refFasta)