#!/bin/bash
#
#SBATCH --job-name=run_rewrite
#SBATCH --output=run_rewrite.out
#SBATCH --error=run_rewrite.err
#
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#
cd $HOME/working/genome3d/parse_cath_domth_output/
export PYTHONPATH=$HOME/Applications/modeller9.19/modlib/:$HOME/Applications/modeller9.19/lib/x86_64-intel8/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Application/modeller9.19/lib/x86_64-intel8/

FILES=($HOME/rewrite_targets/*)
FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}
while read -r line || [[ -n "$line" ]]; do
    i=$line
    j=`echo $i | cut -d. -f1`
    fasta=`echo $i | rev | cut -d"_" -f2- | rev`
    #echo $i
    echo $HOME/Applications/bioserf/bin/rewrite_modeller.pl $HOME/working/genome3d/parse_cath_domth_output/ $HOME/$j.mod_lookups $HOME/$j.blastaligns $HOME/$j.pdomaligns $HOME/$fasta.fasta $HOME/$i $HOME/Applications/bioserf/bin/reformat.pl
done < "$FILE"
