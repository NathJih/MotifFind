# MotifFind

This tool is a motif-finding tool that will plot the top 5 best-match motif sequences with proability information on each location in order with the p-value, log p-value, percent of target, and percent of background for each motif sequence.

## Requirements 
This tool requires the following packages:
- 

## Details 

Jaspar.py 
Creates the Jaspar class and initializes the Jaspar object

JasparMatrices.py
(will add description after changes get pushed)
usage: py JasperMatrices <JASPER file>

getBackground.py
Takes in a peaks file and outputs a list of background sequences for each peak for result validation

getRefGenome.py
Open the chromosome FASTA Format DNA and Protein Sequence Alignment file (fna) and create a single string of the chromosome sequence 
usage: py getRefGenome <fna file> (helper method)

SlidingScore.py
