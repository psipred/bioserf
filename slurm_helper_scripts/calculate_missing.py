import glob
import sys


for path in glob.glob(sys.argv[1]+"*_*_*.pdb"):
    print(path)
