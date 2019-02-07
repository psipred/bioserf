import glob
import re

for name in glob.glob('*.py'):
    print(name)
    with open(name) as fh:
        contents = fh.readlines()
        print(contents)
