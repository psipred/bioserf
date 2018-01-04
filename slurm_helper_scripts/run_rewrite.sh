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
for i in `ls *.ali | sed -e 's/\.ali//'`; do echo $i; j=`echo $i | cut -d. -f1`; echo $HOME/Applications/bioserf/bin/rewrite_modeller.pl $HOME/working/genome3d/parse_cath_domth_ouput/ $HOME/working/genome3d/parse_cath_domth_ouput/$j.mod_lookups $HOME/working/genome3d/parse_cath_domth_ouput/$j.blastaligns $HOME/working/genome3d/parse_cath_domth_ouput/$j.pdomaligns $HOME/working/genome3d/Genome3D.2017-09-05/ecoli_fasta/all/$j.fasta $i.ali $HOME/Applictions/bioserf/bin/reformat.pl; done;
#for i in `ls *.ali | sed -e 's/\.ali//'`; do echo $i; j=`echo $i | cut -d. -f1`; echo $j; done;
