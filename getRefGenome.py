def getRefGenome(chromosome):
    RGFile = open("chromosome" + str(chromosome) + ".fna")
    RGLines = RGFile.readlines()
    RefGen = ""
    for i in range(1, len(RGLines)):
        RefGen.append(RGLines[i])
    return RefGen
