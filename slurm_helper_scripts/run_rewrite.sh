#!/bin/bash
#
#SBATCH --job-name=run_rewrite
#SBATCH --output=run_rewrite.out
#SBATCH --error=run_rewrite.err
##SBATCH --ntasks=1
##SBATCH --array=1-27

cd $HOME/working/genome3d/parse_cath_domth_output/
export PYTHONPATH=$HOME/Applications/modeller9.19/modlib/:$HOME/Applications/modeller9.19/lib/x86_64-intel8/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Application/modeller9.19/lib/x86_64-intel8/

FILES=($HOME/rewrite_targets/*)
FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}
while read -r line || [[ -n "$line" ]]; do
    i=$line
    j=`echo $i | cut -d. -f1`
    ali=`echo $i | cut -d"/" -f4`
    fasta_full=`echo $j | rev | cut -d"_" -f2- | rev`
    fasta=`echo $fasta_full | cut -d"/" -f4`
    # echo $i
    # echo $j
    # echo $ali
    # echo $fasta_full
    # echo $fasta
    echo $HOME/Applications/bioserf/bin/rewrite_modeller.pl $HOME/working/genome3d/parse_cath_domth_output/ $HOME/$fasta_full.mod_lookups $HOME/$fasta_full.blastaligns $HOME/$fasta_full.pdomaligns $HOME/working/genome3d/Genome3D.2017-09-05/all/$fasta.fasta $ali $HOME/Applications/bioserf/bin/reformat.pl
done < "$FILE"
