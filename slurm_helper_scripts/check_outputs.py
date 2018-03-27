# now that all the domthreaders have run we need to check that all results
# have been produced

import glob

# get a list of all the targets
targets = '/home/camp/buchand/working/genome3d/Genome3D.2017-09-05/all_fasta/'
target_list = []
for target_dir in glob.glob(targets+"*"):
    for target in glob.glob(target_dir+"/*"):
        prot_id = target[target.rfind("/"):-5]
        print(target)
        print(prot_id)
