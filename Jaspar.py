import pickle

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
