#!/bin/bash
#
#SBATCH --job-name=run_modeller
#SBATCH --output=run_modeller.out
#SBATCH --error=run_modeller.err
#
#SBATCH --ntasks=1
#SBATCH --array=1-30

cd $HOME/working/genome3d/parse_cath_domth_output/
export PYTHONPATH=$HOME/Applications/modeller9.19/modlib/:$HOME/Applications/modeller9.19/lib/x86_64-intel8/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Application/modeller9.19/lib/x86_64-intel8/

FILES=($HOME/modeller_targets/*)
FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}
while read -r line || [[ -n "$line" ]]; do
    $HOME/Applications/modeller9.19/bin/mod9.19 $HOME$line
done < "$FILE"

# for i in `find -iname "*.py" | sed -e 's/\.py//'`; do echo $i; $HOME/Applications/modeller9.19/bin/mod9.19 $i.py; done;
