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
    $HOME/Applications/bioserf/bin/rewrite_modeller.pl $HOME/working/genome3d/parse_cath_domth_output/ $HOME/working/genome3d/parse_cath_domth_output/$j.mod_lookups $HOME/working/genome3d/parse_cath_domth_output/$j.blastaligns $HOME/working/genome3d/parse_cath_domth_output/$j.pdomaligns $HOME/working/genome3d/Genome3D.2017-09-05/ecoli_fasta/all/$j.fasta $i $HOME/Applications/bioserf/bin/reformat.pl
done < "$FILE"



cd $HOME/working/genome3d/parse_cath_domth_output/
for i in `find -iname "*.ali" | sed -e 's/\.ali//'`; do echo $i; j=`echo $i | cut -d. -f1`; $HOME/Applications/bioserf/bin/rewrite_modeller.pl $HOME/working/genome3d/parse_cath_domth_output/ $HOME/working/genome3d/parse_cath_domth_output/$j.mod_lookups $HOME/working/genome3d/parse_cath_domth_output/$j.blastaligns $HOME/working/genome3d/parse_cath_domth_output/$j.pdomaligns $HOME/working/genome3d/Genome3D.2017-09-05/ecoli_fasta/all/$j.fasta $i.ali $HOME/Applications/bioserf/bin/reformat.pl; done;
#for i in `ls *.ali | sed -e 's/\.ali//'`; do echo $i; j=`echo $i | cut -d. -f1`; echo $j; done;
