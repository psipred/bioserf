import glob
import sys
import re

pdb_str = "/(.+?)_\d+_\d+\.pdb"
pdb_re = re.compile(pdb_str)
for path in glob.glob(sys.argv[1]+"*_*_*.pdb"):
    m = pdb_re.search(path)
    if m:
        print(m.group(1))
