import subprocess
import glob
import time
import sys

#
# usage:
#
# python submitter [path to batch directories] [path to slurm bash script]
#
#
#

#
# SLURM script must makre use of something like:
#
#  #SBATCH --array=1-1000

#  FILES=($HOME/working/genome3d/Genome3D.2017-09-05/all_fasta/$1/*)
#  out_name=${FILES[$SLURM_ARRAY_TASK_ID-1]}

path = sys.argv[1]
path = path.rstrip("/")
command = sys.argv[2]

for directory in glob.glob(path+'/*'):
    set_number = directory[len(path)+1:]
    # print(set_number)
    # continue
    wait = True
    while wait:
        p = subprocess.Popen(['squeue', '--user', 'buchand'],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        queue_deets, err = p.communicate()
        # print(queue_deets.split('\n'))
        if len(queue_deets.split('\n')) < 3:
            print("sending set: "+set_number)
            batch = subprocess.Popen(['sbatch', command, set_number],
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
            batch_deets, err = batch.communicate()
            print(batch_deets)
            sys.stdout.flush()
            wait = False
    time.sleep(60)
