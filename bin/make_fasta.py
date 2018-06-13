import sys
input_file = sys.argv[1]

print(">Query")
with open(input_file, 'r') as fafh:
    for line in fafh:
        print(line)
