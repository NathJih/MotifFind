ref_genome = open(input("Reference Genome File Directory: "))
lines = ref_genome.readlines()
#finds all header indices
header_lines = []
for i in range(0, len(lines)):
    if "chromosome" in lines[i]:
        header_lines.append(i)
for line in header_lines:
    print(line)

#finds all headers for a specified chromosome
header_lines = []
chrom = str(input("Chromosome: "))
chromo = "chromosome " + chrom
for i in range(0, len(lines)):
    if chromo in lines[i]:
        print(lines[i])
        header_lines.append(i)
for line in header_lines:
    print(line)
