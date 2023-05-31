# MotifFind

This tool is a motif-finding tool that will plot the top 5 best-match motif sequences with proability information on each location in order with the p-value, log p-value, percent of target, and percent of background for each motif sequence.
<br>

## Requirements 
### Download the tool:
This tool requires to download the package named `Motif_Find`. <br>

To run this tool, make sure to have the `random`, `scipy`, `pickle`, and `pandas` libraries installed. You can install them using `pip`:

```
pip install random scipy pickle pandas 
```

### File details:
#### Reference Files:
`ReferenceGenome`: zipped folder with all the human GRCh38.p14 reference genome fasta file by chromosomes <br>
`Jaspar.p`: pickle files with a dictionary that stored all the motif informations from Jaspar

#### Python Files:
`Jaspar.py`: use to prepare known motif with PWM matrices, threshold, and basic informations. Running Jaspar.py will output `Jaspar.p` <br>
`MotifFind.py`: take a peak file specified by user, and output the top 5 most enriched motif with motif graph, p-value, and basic motif informations<br>
<br>

## Usage 

### Step 1: Prepare for reference genome
To run motif finding, you will first need to unzip the `.gz` chromosome reference sequence fasta files in the `ReferenceGenome` folder. You can unzip them using the command:
```
unzip ReferenceGenome
gunzip ReferenceGenome/chromosome*
```
<br>

### Step 2 (optional): Prepare for jaspar reference motif
***Note***: you can skip this step if you want to use the defalt human motifs dataset

To specify your own Jaspar motif dataset, you can browse and download jaspar file using the link <https://jaspar.genereg.net/collection/core/> Run the following command to output a compiled jaspar `.p` file
```
python Jasper <JASPER_reference_file> <Jasper_Output_file>
```
`<JASPER_reference_file>`: jaspar file downloaded from <https://jaspar.genereg.net/> <br>
`<Jasper_Output_file>`: name of the output jaspar pickle file. File name must have a `.p` suffix<br>

<br>

### Step 3: Run motif finding
Use your peak file to run our motif finding tool using the following command:
```
python MotifFind <peak_input_file> <output_file> [-j <Jasper_input_file>] [-o <peakseq_output_file>]
```
`<peak_input_file>`: specify input file in peak file format <br>
`<output_file>`: specify output `.html` file name for motif table <br>
`[-j <Jasper_input_file>]`: 
* optional: specify a specific jaspar pickle file you outputed from step 2
* defalt: `Jasper.p` as the jaspar pickle file for human motifs <br>

`[-o <peakseq_output_file>]`: output a `.txt` file with specific peak file sequences. Apecifty the `.txt` file name after `-o`
