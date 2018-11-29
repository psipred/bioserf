# now that all the domthreaders have run we need to check that all results
# have been produced

import glob
import os
from sets import Set
import shutil

dir_count = 1
seq_count = 0

results = '/home/camp/buchand/working/genome3d/domthreader_output/'
result_list = []
for dom_file in glob.glob(results+"*.presults"):
    # print(dom_file)
    dom_id = dom_file[dom_file.rfind("/")+1:-14]
    result_list.append(dom_id)

# print(result_list)

targets = '/home/camp/buchand/working/genome3d/Genome3D.2017-09-05/all_fasta/'
target_list = []
count_list = 0
prot_count = 0
dir_count = 1

output_dir = '/home/camp/buchand/working/genome3d/Genome3D.2017-09-05/collapsed_fasta/'
if not os.path.isdir(output_dir+str(dir_count)):
    os.makedirs(output_dir+str(dir_count))

for target_dir in glob.glob(targets+"*"):
    for target in glob.glob(target_dir+"/*"):
        prot_id = target[target.rfind("/")+1:-6]
        if prot_id in result_list:
            print(found)
        else:
            prot_count += 1
            shutil.copy2(target, output_dir+str(dir_count)+"/"+prot_id+".fasta")
            if prot_count == 1000:
                dir_count += 1
                prot_count = 0
                if not os.path.isdir(output_dir+str(dir_count)):
                    os.makedirs(output_dir+str(dir_count))

# print(count_list)
