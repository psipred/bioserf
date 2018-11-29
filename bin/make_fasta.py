import sys
input_file = sys.argv[1]

with open(input_file) as f:
    content = f.readlines()
    content = [x.strip() for x in content]
    content_string = "".join(content)
    header_count = content_string.count(">")
    seq_count = 0
    if header_count > 1:
        for line in content:
            if line.startswith(">"):
                print(">Query_"+str(seq_count))
                seq_count+=1
            else:
                print(line)
    else:
        print(">Query")
        for line in content:
            if not line.startswith(">"):
                print(line)
