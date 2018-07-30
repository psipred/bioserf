# now that all the domthreaders have run we need to check that all results
# have been produced

import glob
from sets import Set

dir_count = 1
seq_count = 0

results = '/home/camp/buchand/working/genome3d/domthreader_output/'
result_list = []
for dom_file in glob.glob(results+"*.presults"):
    # print(dom_file)
    dom_id = dom_file[dom_file.rfind("/")+1:-14]
    result_list.append(dom_id)

print(result_list)

targets = '/home/camp/buchand/working/genome3d/Genome3D.2017-09-05/all_fasta/'
target_list = []
count_list = 0
for target_dir in glob.glob(targets+"*"):
    for target in glob.glob(target_dir+"/*"):
        prot_id = target[target.rfind("/")+1:-6]
        if prot_id in result_list:
            print(found)
        else:
            count_list+=1

print(count_list)
