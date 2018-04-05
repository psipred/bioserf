import glob
import sys
import re
import os

#calculate_missing $HOME/working/genome3d/pdbaa_based_models $HOME/working/genome3d/Genome3D.2017-09-05/all_fasta

pdb_str = ".+/(.+?)_\d+_\d+\.pdb"
pdb_re = re.compile(pdb_str)
complete_genes = []
for path in glob.glob(sys.argv[1]+"*_*_*.pdb"):
    # print(path)
    m = pdb_re.search(path)
    if m:
        complete_genes.append(m.group(1))

print(complete_genes)
fasta_str = ".+/(.+?).fasta"
fasta_re = re.compile(fasta_str)
todo_count = 0
for path in glob.glob(sys.argv[2]+"/*"):
    for fasta in glob.glob(path+"/*"):
        m = fasta_re.search(fasta)
        if m:
            if m.group(1) in complete_genes:
                print("Found: "+m.group(1))
                os.remove(fasta)
            else:
                print("TODO: "+m.group(1))
                todo_count += 1

print("TODO: "+str(todo_count))
