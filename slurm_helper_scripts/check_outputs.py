# now that all the domthreaders have run we need to check that all results
# have been produced

import glob

# get a list of all the targets
targets = '/home/camp/buchand/working/genome3d/Genome3D.2017-09-05/all_fasta/'
target_list = {}
for target_dir in glob.glob(targets+"*"):
    for target in glob.glob(target_dir+"/*"):
        prot_id = target[target.rfind("/")+1:-6]
        # print(target)
        # print(prot_id)
        target_list[prot_id] = 0

# print(target_list)

results = '/home/camp/buchand/working/genome3d/domthreader_output/'
result_list = {}
for dom_file in glob.glob(results+"*.presults"):
    print(dom_file)
    dom_id = dom_file[dom_file.rfind("/"+1:-14)]
