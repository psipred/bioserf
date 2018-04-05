#!/bin/bash
#
#SBATCH --job-name=domthreader
#SBATCH --output=messages/domthreader_%J.out
#SBATCH --error=messages/domthreader_%J.err
#
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#
#SBATCH --array=1-1000

path=$HOME/working/genome3d/Genome3D.2017-09-05/all_fasta/$1/
size=${#path}
FILES=($HOME/working/genome3d/Genome3D.2017-09-05/all_fasta/$1/*)
out_name=${FILES[$SLURM_ARRAY_TASK_ID-1]}
out_name=${out_name:size:-6}
result_path=$HOME/working/genome3d/domthreader_output/$out_name.presults

if [ -f $result_path ]; then
   echo "skipping $out_name"
else
  echo "calculating $out_name"
  $HOME/Applications/pGenTHREADER/GenThreader.sh -i ${FILES[$SLURM_ARRAY_TASK_ID-1]} -j $HOME/working/genome3d/domthreader_output/$out_name -d
fi
