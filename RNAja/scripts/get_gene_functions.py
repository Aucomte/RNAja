import re
import csv
import sys
maxInt = sys.maxsize

while True:
    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

listeGene = "/home/comtea/Documents/PHIM/AUTRES/RNASEQ/DE/DownRegFDR0.01logFC1.5.csv"
gff = "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/DATA/ref/msu7.gff3"
output = "/home/comtea/Documents/PHIM/AUTRES/RNASEQ/DE/TEST.csv"

gene_functions={}
with open(gff, 'r') as f2:
    spamreader = csv.reader(f2, delimiter="\t")
    next(spamreader, None)
    for line in spamreader:
        gene = ""
        if line[2] == "gene":
            gene=re.search(r"ID=([^;,]+)[,;]",line[8]).group(1)
            gene_functions[gene] = line[8]


with open(output, 'w') as fichier, open(listeGene, 'r') as f1:
    spamreader = csv.reader(f1, delimiter=",")
    for line in spamreader:
        if line[0] == "genes":
            line = f"{','.join(line)},functions\n"
            fichier.write(line)
        else:
            gene1 = line[0]
            line = f"{','.join(line)},{gene_functions[gene1]}\n"
            fichier.write(line)

