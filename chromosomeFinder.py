ref_genome = open(input("Reference Genome File Directory: "))
lines = ref_genome.readlines()
chromfile = open(input("Empty file directory for chromosome sequence: "), 'w')
start = int(input("Start index of Primary Assembly: "))
chromfile.write(lines[start])
start = start + 1
for i in range(start, len(lines)):
    if "chromosome" not in lines[i]:
        chromfile.write(lines[i])
    else:
        break
chromfile.close()
