# now that all the domthreaders have run we need to check that all results
# have been produced

import glob

# get a list of all the targets
targets = '/home/camp/buchand/working/genome3d/Genome3D.2017-09-05/all_fasta/'
target_list = []
for target in glob.glob(targets+"*"):
    print(target)
