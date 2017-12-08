import glob
import sys
import re

pdb_str = ".+/(.+?)_\d+_\d+\.pdb"
pdb_re = re.compile(pdb_str)
complete_genes = []
for path in glob.glob(sys.argv[1]+"*_*_*.pdb"):
    m = pdb_re.search(path)
    if m:
        complete_genes.append(m.group(1))

fasta_str = ".+/(.+?).fasta"
fasta_re = re.compile(fasta_str)
for path in glob.glob(sys.argv[2]+"/*"):
    for fasta in glob.glob(path+"/*"):
        m = fasta.re.search(fasta)
        if m:
            print(m.group(1))
