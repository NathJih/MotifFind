class Jasper:
    # initialize an empty dictionary
    def __init__(self):
        self.matrixdict = {}
    
    # add JASPAR elements to the dictionary
    def add(self, name, matrix, threshold, length, TF):
        self.matrixdict[name] = [matrix, threshold, length, TF]
