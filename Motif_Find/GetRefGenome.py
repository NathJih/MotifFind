import pickle

def GetRefGenome():
    RefG = {}
    chromosome = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"]
    for chr in chromosome:
        RGFile = open("ReferenceGenome/chromosome" + str(chr) + ".fna")
        RGLines = RGFile.readlines()
        RefGen = ""
        for i in range(1, len(RGLines)):
            RefGen = RefGen + RGLines[i].strip()
        RefG["chromosome" + str(chr)] = RefGen
    
    file = open('GRCh38p14.p', 'wb')
    pickle.dump(RefG, file)
    file.close()
    
    return

GetRefGenome()