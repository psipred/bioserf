#!/bin/bash
#
#SBATCH --job-name=domthreader
#SBATCH --output=domthreader.out
#SBATCH --error=domthreader.err
#
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#
#SBATCH --array=1-1000

FILES=($HOME/working/genome3d/Genome3D.2017-09-05/all_fasta/$1/*)
out_name=${FILES[$SLURM_ARRAY_TASK_ID-1]}
out_name=${out_name:70:-6}
$HOME/Applications/pGenTHREADER/GenThreader.sh -i ${FILES[$SLURM_ARRAY_TASK_ID-1]} -j $HOME/working/genome3d/domthreader_output/$out_name -d
