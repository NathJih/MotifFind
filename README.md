# MotifFind

This tool is a motif-finding tool that will plot the top 5 best-match motif sequences with proability information on each location in order and a CSV file with the motif name, p-value, number of target matched, percent of target, number of bachground matched, percent of background, and the motif documentation for each motif sequence.
<br>

## Requirements 
### Download the tool:
To download and use the tool, please run the following command on your terminal:
```
git clone https://github.com/NathJih/MotifFind
cd Motif_Find
```

Make sure to have the `scipy`, `logomaker` and `pandas` libraries installed. You can install them using `pip`:

```
pip install scipy logomaker pandas
```

### File details:
#### Motif_Find
Motif_Find is the mean directory that you need to use to run the tool. Motif_Find including python script and datafiles including: <br>
`Jaspar.p`: pickle files with a dictionary that stored all the motif informations from Jaspar <br>
`Jaspar.py`: use to prepare known motif with PWM matrices, threshold, and basic informations. Running Jaspar.py will output `Jaspar.p` <br>
`MotifFind.py`: take a peak file specified by user, and output the top 5 most enriched motif with motif graph, p-value, and basic motif informations <br>
`GetRefGenome.py`: takes in a reference genome fasta file and returns a pickle file with only the "primary assemblies" of each chromosome to be used by other methods <br>
#### Analysis: contain sample peak files of TF FOXA1 and POLR2A, and the corresponding motif finding output result.
#### Benchmark: contain the peak file, Jaspar dataset, and output to run our benchmarking result.
#### TestFiles: contain several small test files we used to test our tool during implementation

<br>

## Usage 

### Step 1: Prepare the reference genome
To run motif finding, you will first need the reference genome pickle file, which can be obtained either by downloading it from this link: https://drive.google.com/uc?export=download&id=1TnI-1oExDtkpbal4kAWvOTWnuh7ksdmW, or by recreating it locally with your own reference genome fasta file (GRCh38p14 recommended). You can  create the file by using the command:
```
python GetRefGenome.py <reference_genome_fasta_file> <Output_path>
```
We recommand to store your reference genome output in your current (Motif_Find) directory <br>
Example usage: (note, make sure you are in the Motif_Find directory) `python GetRefGenome.py GRCh38.p14.fna GRCh38p14` 
<br>
<br>

### Step 2 (optional): Prepare the jaspar reference motif
***Note***: you can skip this step if you want to use the default human motifs dataset

To specify your own Jaspar motif dataset, you can browse and download jaspar file using the link <https://jaspar.genereg.net/collection/core/> Run the following command to output a compiled jaspar `.p` file
```
python Jasper.py <JASPER_reference_file> <Jasper_Output_file>
```
`<JASPER_reference_file>`: jaspar file downloaded from <https://jaspar.genereg.net/> <br>
`<Jasper_Output_file>`: name of the output jaspar pickle file. File name must have a `.p` suffix<br>

<br>

### Step 3: Run motif finding
Use your peak file to run our motif finding tool using the following command:
```
python MotifFind.py <peak_input_file> <output_file> -r <reference_sequence_pickle_file> [-j <Jasper_input_file>] [-o <peakseq_output_file>]
```
`<peak_input_file>`: specify input peak file in bed format <br>
`<output_file>`: specify output folder name for motif table (Note: do not specify the path. Your output folder must be in the currect directory.)<br>
`-r <reference_sequence_pickle_file>`: specify the reference genome dataset you want to use for running motif finding.<br>
`[-j <Jasper_input_file>]`: 
* optional: specify a specific jaspar pickle file you outputed from step 2
* default: `Jasper.p` as the jaspar pickle file for human motifs <br>

`[-o <peakseq_output_file>]`: output a `.txt` file with specific peak file sequences. Specify the `.txt` file name after `-o`

Example usage: (note, make sure you are in the Motif_Find directory)
`python MotifFind.py ../TestFiles/CTCF.bed CTCF_MFresults -r GRCh38p14.p` 

<br>

## Examples to run the tool
### Run our benchmarking tool
To run our benchmarking tool, please first download the preprocessed mouse mm39 reference genome from this link (https://drive.google.com/file/d/10EQK9ig_P59iKwBHlHM9BuOUDyTsXlGo/view?usp=drive_link), move the downloaded file into your local Benchmark folder, and run the following command: 
<br>
<br>
(Note: in these commands, you will use preprocessed files inside folder Benchmark. We have already processed the klf4 peak file using HOMER in Jupyter Hub by running the command `pos2bed.pl peaks.txt > klf4.bed`)
```
python Jaspar.py ../Benchmark/JASPARall.txt ../Benchmark/JASPARall.p
python MotifFind.py ../benchmark/klf4.bed klf4_MFresults -r ../Benchmark/GRCm39.p -j ../Benchmark/JASPARall.p
```

### Run our sample peak file of transcription factor POLR2A
If you downloaded the human reference genome fasta file instead of the pickle file we provided, please run the following command before processing the data:
```
python GetRefGenome.py <human_ref_gen_fasta> GRCh38p14.p
```

Then run the commands below:
```
python MotifFind.py ../Analysis/polr2a.bed polr2a_MFresults -r GRCh38p14.p
```

## Contributor
This tool is generated by Cici Bu, Jessica Yu, Nathan Jih. <br>

Please email this address [xbu@ucsd.edu] with any problems or questions you have, we will get back to you within 24 hours.
