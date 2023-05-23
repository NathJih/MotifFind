import Jaspar 
import sys

# read in the arguments
args = sys.argv 

# check for the correct number of arguments
if len(args) != 2:
    raise Exception("Incorrect number of arguments")

filename = args[1]

