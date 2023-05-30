import pickle
import sys
import random
import MotifFind
# read in the arguments
args = sys.argv 

# check for the correct number of arguments
if len(args) != 3:
    raise Exception("Incorrect number of arguments, 3 expected.")

filename = args[1]
peaks_file = args[2]

class Jaspar:
    # initialize an empty dictionary
    def __init__(self):
        self.matrixdict = {}
    
    # add JASPAR elements to the dictionary
    def add(self, name, matrix, threshold, length, TF, url):
        self.matrixdict[name] = [matrix, threshold, length, TF, url]

    # write dictionary using pickle and out put a pickle file for later usage
    def store(self):
        file = open('Jaspar.p', 'wb')
        pickle.dump(self.matrixdict, file)
        file.close()

def JasparMatrices(file):
    jasparfile = open(file, 'r')
    line = jasparfile.readline()
    jaspar = Jaspar.Jaspar()
    
    n = 1
    # read the entire file
    while (len(line) != 0):

        # get information from each motif
        if (line[0:5]) == 'MOTIF':

            # get motif name
            title = line.split(" ")
            title = [x.rstrip() for x in title]
            motifname = title[1]

            # get documented motuif related transcription factor
            tflist = title[2].split(".")
            tf = tflist[-1]

            # get motif length
            line = jasparfile.readline()
            title = line.split(" ")
            for i in range(len(title)):
                if title[i] == "w=":
                    length = int(title[i+1])
                    break

            # get motif PWM matrix
            # A - first row
            # C - second row
            # G - third row
            # T - forth row
            matrix = []
            alist = []
            clist = []
            glist = []
            tlist = []
            line = jasparfile.readline()
            while line[0] != "U":
                numbers = line.split("  ")
                numbers = [x.rstrip() for x in numbers]
                alist.append(float(numbers[0]))
                clist.append(float(numbers[1]))
                glist.append(float(numbers[2]))
                tlist.append(float(numbers[3]))
                line = jasparfile.readline()
            matrix.append(alist)
            matrix.append(clist)
            matrix.append(glist)
            matrix.append(tlist)

            # get url from the last line
            urllist = line.split(" ")
            urllist = [x.rstrip() for x in urllist]
            url = urllist[1]

            # get threshold using JasparThreshold()
            threshold = JasparThreshold(matrix, length)

            # add motif info into dictionary in Jaspar class
            jaspar.add(motifname, matrix, threshold, length, tf, url)
            print("finished " + str(n))
            n += 1

            # compare motif to peak and background sequences
            MotifFind(jaspar, peaks_file)

        line = jasparfile.readline()

    # After adding all the motif into Jaspar class, output the motif dict
    # into a pickle file
    jaspar.store()
    return


def JasparThreshold(matrix, length):
    rowinfo = {"A": 0, "C": 1, "G": 2, "T": 3}
    bases = ["A", "C", "G", "T"]
    scores = []
    for i in range(500):
        score = 0
        for j in range(length):
            r = random.randint(0,3)
            base = bases[r]
            score += matrix[rowinfo[base]][j]
        scores.append(score)
    
    scores.sort()
    index = int(len(scores)*(1 - 0.025))
    thresh = scores[index]
    return thresh

JasparMatrices(filename)