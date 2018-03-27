# now that all the domthreaders have run we need to check that all results
# have been produced

import glob
from sets import Set

# get a list of all the targets
targets = '/home/camp/buchand/working/genome3d/Genome3D.2017-09-05/all_fasta/'
target_list = []
for target_dir in glob.glob(targets+"*"):
    for target in glob.glob(target_dir+"/*"):
        prot_id = target[target.rfind("/")+1:-6]
        # print(target)
        # print(prot_id)
        target_list.append(prot_id)

# print(target_list)

results = '/home/camp/buchand/working/genome3d/domthreader_output/'
result_list = []
for dom_file in glob.glob(results+"*.presults"):
    print(dom_file)
    dom_id = dom_file[dom_file.rfind("/")+1:-14]
    result_list.append(dom_id)


target_set = Set(target_list)
result_set = Set(result_list)
missing_set = target.difference(result_set)
print(missing_set)
