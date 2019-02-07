import glob
import re

for name in glob.glob('*.py'):
    print(name)
    with open(“filename”) as file:
        line = file.readline()
        if "knowns=(" in line:
            print(line)
